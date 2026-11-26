#include "rlwe.h"
#include "prg.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* rlwe.c -> BGV over R_0[Y]/(Y^h) with a negacyclic NTT for R_0.

   Arithmetic modulo each prime q < 2^62 uses Montgomery multiplication
   (R = 2^64) for pointwise products and Shoup multiplication for the NTT
   twiddle factors.  Products in R are block convolutions truncated at
   Y^h, computed pointwise in the NTT domain. */

typedef rlwe_u128 u128;

/* Arithmetic modulo q */

static inline uint64_t add_mod(uint64_t a, uint64_t b, uint64_t q)
{
    uint64_t s = a + b;
    return s >= q ? s - q : s;
}

static inline uint64_t sub_mod(uint64_t a, uint64_t b, uint64_t q)
{
    return a >= b ? a - b : a + q - b;
}

static uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t q)
{
    return (uint64_t)(((u128)a * b) % q);
}

static uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t q)
{
    uint64_t r = 1 % q;
    a %= q;
    while (e) {
        if (e & 1) r = mul_mod(r, a, q);
        a = mul_mod(a, a, q);
        e >>= 1;
    }
    return r;
}

/* a b 2^{-64} mod q, for a, b < q < 2^62 */
static inline uint64_t mont_mul(uint64_t a, uint64_t b, uint64_t q, uint64_t qinv)
{
    u128 T = (u128)a * b;
    uint64_t m = (uint64_t)T * qinv;
    uint64_t u = (uint64_t)((T + (u128)m * q) >> 64);
    return u >= q ? u - q : u;
}

/* a w mod q with the Shoup companion wp = floor(w 2^64 / q) */
static inline uint64_t shoup_mul(uint64_t a, uint64_t w, uint64_t wp, uint64_t q)
{
    uint64_t hi = (uint64_t)(((u128)a * wp) >> 64);
    uint64_t r = a * w - hi * q;
    return r >= q ? r - q : r;
}

static int is_prime_u64(uint64_t n)
{
    if (n < 2) return 0;
    static const uint64_t small[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (size_t i = 0; i < sizeof(small) / sizeof(small[0]); i++) {
        if (n == small[i]) return 1;
        if (n % small[i] == 0) return 0;
    }
    uint64_t d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }
    for (size_t i = 0; i < sizeof(small) / sizeof(small[0]); i++) {
        uint64_t x = pow_mod(small[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int composite = 1;
        for (int r = 1; r < s; r++) {
            x = mul_mod(x, x, n);
            if (x == n - 1) { composite = 0; break; }
        }
        if (composite) return 0;
    }
    return 1;
}

static unsigned brv(unsigned x, int bits)
{
    unsigned r = 0;
    for (int i = 0; i < bits; i++) { r = (r << 1) | (x & 1); x >>= 1; }
    return r;
}

static void ntt_table_init(NTTTable *T, uint64_t q, int N0, int logN0)
{
    T->q = q;
    uint64_t inv = 1;                           /* q^{-1} mod 2^64 by Newton */
    for (int i = 0; i < 6; i++) inv *= 2 - q * inv;
    T->qinv = (uint64_t)(-inv);
    uint64_t r1 = (uint64_t)(((u128)1 << 64) % q);
    T->r2 = mul_mod(r1, r1, q);
    T->ninv = pow_mod((uint64_t)N0, q - 2, q);

    /* psi: a primitive 2N0-th root of unity, psi^{N0} = -1 */
    uint64_t psi = 0;
    for (uint64_t x = 2; x < q; x++) {
        psi = pow_mod(x, (q - 1) / (2 * (uint64_t)N0), q);
        if (pow_mod(psi, (uint64_t)N0, q) == q - 1) break;
    }
    uint64_t ipsi = pow_mod(psi, q - 2, q);
    T->psi = malloc((size_t)N0 * sizeof(uint64_t));
    T->psip = malloc((size_t)N0 * sizeof(uint64_t));
    T->ipsi = malloc((size_t)N0 * sizeof(uint64_t));
    T->ipsip = malloc((size_t)N0 * sizeof(uint64_t));
    if (!T->psi || !T->psip || !T->ipsi || !T->ipsip) abort();
    for (int k = 0; k < N0; k++) {
        unsigned e = brv((unsigned)k, logN0);
        T->psi[k] = pow_mod(psi, e, q);
        T->ipsi[k] = pow_mod(ipsi, e, q);
        T->psip[k] = (uint64_t)(((u128)T->psi[k] << 64) / q);
        T->ipsip[k] = (uint64_t)(((u128)T->ipsi[k] << 64) / q);
    }
}

static void ntt_table_free(NTTTable *T)
{
    free(T->psi); free(T->psip); free(T->ipsi); free(T->ipsip);
}

/* Forward negacyclic NTT, output in bit-reversed order (Cooley-Tukey). */
static void ntt_fwd(const NTTTable *T, uint64_t *a, int N0)
{
    const uint64_t q = T->q;
    int t = N0;
    for (int m = 1; m < N0; m <<= 1) {
        t >>= 1;
        for (int i = 0; i < m; i++) {
            int j1 = 2 * i * t;
            uint64_t S = T->psi[m + i], Sp = T->psip[m + i];
            for (int j = j1; j < j1 + t; j++) {
                uint64_t U = a[j];
                uint64_t V = shoup_mul(a[j + t], S, Sp, q);
                a[j] = add_mod(U, V, q);
                a[j + t] = sub_mod(U, V, q);
            }
        }
    }
}

/* Inverse negacyclic NTT (Gentleman-Sande), input in bit-reversed order. */
static void ntt_inv(const NTTTable *T, uint64_t *a, int N0)
{
    const uint64_t q = T->q;
    int t = 1;
    for (int m = N0; m > 1; m >>= 1) {
        int j1 = 0, hm = m >> 1;
        for (int i = 0; i < hm; i++) {
            uint64_t S = T->ipsi[hm + i], Sp = T->ipsip[hm + i];
            for (int j = j1; j < j1 + t; j++) {
                uint64_t U = a[j], V = a[j + t];
                a[j] = add_mod(U, V, q);
                a[j + t] = shoup_mul(sub_mod(U, V, q), S, Sp, q);
            }
            j1 += 2 * t;
        }
        t <<= 1;
    }
    uint64_t ninvp = (uint64_t)(((u128)T->ninv << 64) / q);
    for (int j = 0; j < N0; j++) a[j] = shoup_mul(a[j], T->ninv, ninvp, q);
}

/* Elements of R_Q */

static size_t ring_len(const RLWEParams *pp)
{
    return (size_t)pp->h * (size_t)pp->L * (size_t)pp->N0;
}

static uint64_t *ring_alloc(const RLWEParams *pp)
{
    uint64_t *x = calloc(ring_len(pp), sizeof(uint64_t));
    if (!x) abort();
    return x;
}

static inline uint64_t *blk(const RLWEParams *pp, uint64_t *x, int k, int i)
{
    return x + ((size_t)k * pp->L + (size_t)i) * pp->N0;
}

static inline const uint64_t *cblk(const RLWEParams *pp, const uint64_t *x, int k, int i)
{
    return x + ((size_t)k * pp->L + (size_t)i) * pp->N0;
}

static void ring_ntt(const RLWEParams *pp, uint64_t *x)
{
    for (int k = 0; k < pp->h; k++)
        for (int i = 0; i < pp->L; i++) ntt_fwd(&pp->ntt[i], blk(pp, x, k, i), pp->N0);
}

static void ring_intt(const RLWEParams *pp, uint64_t *x)
{
    for (int k = 0; k < pp->h; k++)
        for (int i = 0; i < pp->L; i++) ntt_inv(&pp->ntt[i], blk(pp, x, k, i), pp->N0);
}

/* out = a b in R_Q, NTT domain, with a in Montgomery form: block k of the
   product is sum_{i + j = k} a_i b_j, and blocks at Y^h and above vanish. */
static void ring_mul_ntt(const RLWEParams *pp, uint64_t *out, const uint64_t *a_mont,
                         const uint64_t *b)
{
    const int N0 = pp->N0;
    memset(out, 0, ring_len(pp) * sizeof(uint64_t));
    for (int i = 0; i < pp->h; i++)
        for (int j = 0; i + j < pp->h; j++)
            for (int p = 0; p < pp->L; p++) {
                const uint64_t q = pp->ntt[p].q, qinv = pp->ntt[p].qinv;
                const uint64_t *x = cblk(pp, a_mont, i, p), *y = cblk(pp, b, j, p);
                uint64_t *z = blk(pp, out, i + j, p);
                for (int c = 0; c < N0; c++)
                    z[c] = add_mod(z[c], mont_mul(x[c], y[c], q, qinv), q);
            }
}

static void ring_to_mont(const RLWEParams *pp, uint64_t *x)
{
    for (int k = 0; k < pp->h; k++)
        for (int p = 0; p < pp->L; p++) {
            const NTTTable *T = &pp->ntt[p];
            uint64_t *y = blk(pp, x, k, p);
            for (int c = 0; c < pp->N0; c++) y[c] = mont_mul(y[c], T->r2, T->q, T->qinv);
        }
}

static void ring_add(const RLWEParams *pp, uint64_t *out, const uint64_t *a, const uint64_t *b)
{
    for (int k = 0; k < pp->h; k++)
        for (int p = 0; p < pp->L; p++) {
            const uint64_t q = pp->ntt[p].q;
            const uint64_t *x = cblk(pp, a, k, p), *y = cblk(pp, b, k, p);
            uint64_t *z = blk(pp, out, k, p);
            for (int c = 0; c < pp->N0; c++) z[c] = add_mod(x[c], y[c], q);
        }
}

/* Load small signed coefficients (h * N0 of them) into every prime. */
static void ring_from_small(const RLWEParams *pp, uint64_t *out, const int64_t *v)
{
    for (int k = 0; k < pp->h; k++)
        for (int p = 0; p < pp->L; p++) {
            const uint64_t q = pp->ntt[p].q;
            uint64_t *z = blk(pp, out, k, p);
            const int64_t *s = v + (size_t)k * pp->N0;
            for (int c = 0; c < pp->N0; c++)
                z[c] = s[c] >= 0 ? (uint64_t)s[c] % q : q - ((uint64_t)(-s[c]) % q);
        }
}

/* Sampling */

static void sample_ternary(const RLWEParams *pp, int64_t *v)
{
    size_t n = (size_t)pp->h * pp->N0, i = 0;
    uint8_t buf[256];
    while (i < n) {
        rnd_bytes(buf, sizeof(buf));
        for (size_t j = 0; j < sizeof(buf) && i < n; j++)
            if (buf[j] < 255) v[i++] = (int64_t)(buf[j] % 3) - 1;
    }
}

/* Discrete Gaussian of standard deviation sigma, truncated at |e| <= Be,
   by a cumulative distribution table of |e| and a random sign. */
static void sample_gauss(const RLWEParams *pp, int64_t *v)
{
    size_t n = (size_t)pp->h * pp->N0;
    for (size_t i = 0; i < n; i++) {
        uint64_t u = rnd_u64();
        int64_t k = 0;
        while (k < pp->Be && u >= pp->cdt[k]) k++;
        if (k && (rnd_u64() & 1)) k = -k;
        v[i] = k;
    }
}

static void gauss_table(RLWEParams *pp)
{
    double rho[64], S = 0;
    for (int k = 0; k <= pp->Be; k++) {
        rho[k] = exp(-(double)k * k / (2.0 * pp->sigma * pp->sigma));
        S += k ? 2 * rho[k] : rho[k];
    }
    double cum = 0;
    for (int k = 0; k <= pp->Be; k++) {
        cum += (k ? 2 * rho[k] : rho[k]) / S;
        pp->cdt[k] = cum >= 1.0 ? UINT64_MAX : (uint64_t)ldexp(cum, 64);
    }
    pp->cdt[pp->Be] = UINT64_MAX;
}

/* Setup */

static uint64_t find_prime(uint64_t start, uint64_t step, uint64_t avoid)
{
    /* largest prime p <= start with p = 1 mod step, p != avoid */
    uint64_t p = start - ((start - 1) % step);
    for (; p > step; p -= step)
        if (p != avoid && is_prime_u64(p)) return p;
    abort();
}

void rlwe_setup(RLWEParams *pp, int N0, int h, int logQ, int L, uint64_t t)
{
    memset(pp, 0, sizeof(*pp));
    pp->N0 = N0;
    pp->h = h;
    pp->L = L;
    pp->n = N0;
    pp->t = t;
    pp->logN0 = 0;
    while ((1 << pp->logN0) < N0) pp->logN0++;
    if (L < 1 || L > RLWE_MAX_PRIMES || logQ > 120 || (t - 1) % (2 * (uint64_t)N0)) abort();

    /* Q = q_1 ... q_L just below 2^logQ, each q_i = 1 mod 2N0 and below 2^62. */
    double share = (double)logQ / L;
    uint64_t start = (uint64_t)floor(pow(2.0, share));
    if (share >= 62) abort();
    pp->Q = 1;
    for (int i = 0; i < L; i++) {
        pp->q[i] = find_prime(start, 2 * (uint64_t)N0, i ? pp->q[0] : 0);
        pp->Q *= pp->q[i];
        ntt_table_init(&pp->ntt[i], pp->q[i], N0, pp->logN0);
    }
    if (L == 2) pp->crt_inv = pow_mod(pp->q[0] % pp->q[1], pp->q[1] - 2, pp->q[1]);
    pp->logQ = 0;
    while (pp->logQ < 128 && ((u128)1 << pp->logQ) < pp->Q) pp->logQ++;
    ntt_table_init(&pp->ntt_t, t, N0, pp->logN0);

    pp->sigma = 3.2;
    pp->Be = (int)ceil(6 * pp->sigma);
    gauss_table(pp);

    /* KeyGen: s ternary, a uniform, e Gaussian; pk = (-(a s + t e), a). */
    size_t n = (size_t)h * N0;
    int64_t *sv = malloc(n * sizeof(int64_t)), *ev = malloc(n * sizeof(int64_t));
    if (!sv || !ev) abort();
    pp->pk0 = ring_alloc(pp);
    pp->pk1 = ring_alloc(pp);
    pp->sk = ring_alloc(pp);
    uint64_t *tmp = ring_alloc(pp), *E = ring_alloc(pp);
    for (int k = 0; k < h; k++)
        for (int p = 0; p < L; p++) {
            uint64_t *a = blk(pp, pp->pk1, k, p);
            for (int c = 0; c < N0; c++) a[c] = rnd_below(pp->q[p]);   /* uniform, NTT domain */
        }
    sample_ternary(pp, sv);
    sample_gauss(pp, ev);
    for (size_t i = 0; i < n; i++) ev[i] *= (int64_t)t;
    ring_from_small(pp, pp->sk, sv);
    ring_ntt(pp, pp->sk);
    ring_from_small(pp, E, ev);
    ring_ntt(pp, E);
    ring_to_mont(pp, pp->pk1);
    ring_mul_ntt(pp, tmp, pp->pk1, pp->sk);      /* a s */
    ring_add(pp, tmp, tmp, E);                   /* a s + t e */
    for (int k = 0; k < h; k++)
        for (int p = 0; p < L; p++) {
            uint64_t *x = blk(pp, tmp, k, p), *y = blk(pp, pp->pk0, k, p);
            for (int c = 0; c < N0; c++) y[c] = x[c] ? pp->q[p] - x[c] : 0;
        }
    ring_to_mont(pp, pp->pk0);
    ring_to_mont(pp, pp->sk);
    free(sv); free(ev); free(tmp); free(E);

    sig_keygen(&pp->sig, "ML-DSA-44");
}

void rlwe_free(RLWEParams *pp)
{
    for (int i = 0; i < pp->L; i++) ntt_table_free(&pp->ntt[i]);
    ntt_table_free(&pp->ntt_t);
    free(pp->pk0); free(pp->pk1); free(pp->sk);
    sig_free(&pp->sig);
}

void rlwe_ct_init(const RLWEParams *pp, RLWECt *c)
{
    c->c0 = ring_alloc(pp);
    c->c1 = ring_alloc(pp);
}

void rlwe_ct_free(RLWECt *c)
{
    free(c->c0); free(c->c1);
}

void rlwe_ct_copy(const RLWEParams *pp, RLWECt *dst, const RLWECt *src)
{
    memcpy(dst->c0, src->c0, ring_len(pp) * sizeof(uint64_t));
    memcpy(dst->c1, src->c1, ring_len(pp) * sizeof(uint64_t));
}

size_t rlwe_tag_bytes(const RLWEParams *pp)
{
    return ((size_t)2 * pp->N0 * pp->h * pp->logQ + 7) / 8;
}

/* Encryption: c0 = p0 u + t e0 + m, c1 = p1 u + t e1. */
void rlwe_enc(const RLWEParams *pp, RLWECt *c, const uint64_t *m)
{
    size_t n = (size_t)pp->h * pp->N0;
    int64_t *v = malloc(n * sizeof(int64_t));
    uint64_t *U = ring_alloc(pp), *E = ring_alloc(pp);
    if (!v) abort();
    sample_ternary(pp, v);
    ring_from_small(pp, U, v);
    ring_ntt(pp, U);

    for (int comp = 0; comp < 2; comp++) {
        uint64_t *out = comp ? c->c1 : c->c0;
        ring_mul_ntt(pp, out, comp ? pp->pk1 : pp->pk0, U);
        ring_intt(pp, out);
        sample_gauss(pp, v);
        for (size_t i = 0; i < n; i++) {
            v[i] *= (int64_t)pp->t;
            if (comp == 0 && m) v[i] += (int64_t)m[i];
        }
        ring_from_small(pp, E, v);
        ring_add(pp, out, out, E);
    }
    free(v); free(U); free(E);
}

/* Centered representative of the CRT combination of the residues. */
static rlwe_i128 crt_center(const RLWEParams *pp, const uint64_t *res)
{
    u128 x;
    if (pp->L == 1) {
        x = res[0];
    } else {
        uint64_t q0 = pp->q[0], q1 = pp->q[1];
        uint64_t k = mul_mod(sub_mod(res[1] % q1, res[0] % q1, q1), pp->crt_inv, q1);
        x = (u128)res[0] + (u128)q0 * k;
    }
    if (x > pp->Q / 2) return (rlwe_i128)x - (rlwe_i128)pp->Q;
    return (rlwe_i128)x;
}

double rlwe_dec(const RLWEParams *pp, uint64_t *m, const RLWECt *c)
{
    uint64_t *C1 = ring_alloc(pp), *V = ring_alloc(pp);
    memcpy(C1, c->c1, ring_len(pp) * sizeof(uint64_t));
    ring_ntt(pp, C1);
    ring_mul_ntt(pp, V, pp->sk, C1);   /* s c1 */
    ring_intt(pp, V);
    ring_add(pp, V, V, c->c0);          /* c0 + s c1 */
    u128 maxabs = 0;
    for (int k = 0; k < pp->h; k++)
        for (int c2 = 0; c2 < pp->N0; c2++) {
            uint64_t res[RLWE_MAX_PRIMES];
            for (int p = 0; p < pp->L; p++) res[p] = blk(pp, V, k, p)[c2];
            rlwe_i128 x = crt_center(pp, res);
            u128 ax = x < 0 ? (u128)(-x) : (u128)x;
            if (ax > maxabs) maxabs = ax;
            rlwe_i128 r = x % (rlwe_i128)pp->t;
            if (r < 0) r += pp->t;
            m[(size_t)k * pp->N0 + c2] = (uint64_t)r;
        }
    free(C1); free(V);
    double lq = log2((double)pp->Q) - 1.0;
    return maxabs ? lq - log2((double)maxabs) : lq;
}

void rlwe_slot_unit(const RLWEParams *pp, uint64_t *block, int i)
{
    memset(block, 0, (size_t)pp->N0 * sizeof(uint64_t));
    block[i] = 1;
    ntt_inv(&pp->ntt_t, block, pp->N0);
}

void rlwe_slot_decode(const RLWEParams *pp, uint64_t *slots, const uint64_t *block)
{
    memcpy(slots, block, (size_t)pp->N0 * sizeof(uint64_t));
    ntt_fwd(&pp->ntt_t, slots, pp->N0);
}

static size_t ct_raw_bytes(const RLWEParams *pp)
{
    return 2 * ring_len(pp) * sizeof(uint64_t);
}

static uint8_t *ct_raw(const RLWEParams *pp, const RLWECt *w)
{
    uint8_t *buf = malloc(ct_raw_bytes(pp));
    if (!buf) abort();
    memcpy(buf, w->c0, ring_len(pp) * sizeof(uint64_t));
    memcpy(buf + ring_len(pp) * sizeof(uint64_t), w->c1, ring_len(pp) * sizeof(uint64_t));
    return buf;
}

size_t rlwe_tag(const RLWEParams *pp, RLWECt *w, uint8_t *sig, int i, int b, int hprime)
{
    uint64_t *m = calloc((size_t)pp->h * pp->N0, sizeof(uint64_t));
    if (!m) abort();
    if (b) rlwe_slot_unit(pp, m + (size_t)(pp->h - hprime) * pp->N0, i - 1);
    rlwe_enc(pp, w, m);
    free(m);
    uint8_t *buf = ct_raw(pp, w);
    size_t siglen = sig_sign(&pp->sig, buf, ct_raw_bytes(pp), sig, pp->sig.sig_bytes);
    free(buf);
    return siglen;
}

int rlwe_verify_tag(const RLWEParams *pp, const RLWECt *w, const uint8_t *sig, size_t siglen)
{
    uint8_t *buf = ct_raw(pp, w);
    int ok = sig_verify(&pp->sig, buf, ct_raw_bytes(pp), sig, siglen);
    free(buf);
    return ok;
}

void rlwe_degrade(const RLWEParams *pp, RLWECt *out, const RLWECt *w)
{
    RLWECt z;
    rlwe_ct_init(pp, &z);
    rlwe_enc(pp, &z, NULL);
    size_t bl = (size_t)pp->L * pp->N0;
    for (int comp = 0; comp < 2; comp++) {
        uint64_t *o = comp ? out->c1 : out->c0;
        const uint64_t *in = comp ? w->c1 : w->c0;
        /* multiplication by Y: block k <- block k-1, block 0 <- 0 (in may alias o) */
        memmove(o + bl, in, (size_t)(pp->h - 1) * bl * sizeof(uint64_t));
        memset(o, 0, bl * sizeof(uint64_t));
        ring_add(pp, o, o, comp ? z.c1 : z.c0);
    }
    rlwe_ct_free(&z);
}

void rlwe_merge(const RLWEParams *pp, RLWECt *out, RLWECt *const *ws, int m)
{
    RLWECt acc;
    rlwe_ct_init(pp, &acc);
    rlwe_enc(pp, &acc, NULL);
    for (int j = 0; j < m; j++) {
        ring_add(pp, acc.c0, acc.c0, ws[j]->c0);
        ring_add(pp, acc.c1, acc.c1, ws[j]->c1);
    }
    rlwe_ct_copy(pp, out, &acc);
    rlwe_ct_free(&acc);
}

int rlwe_trace(const RLWEParams *pp, const RLWECt *w, uint8_t *in_S)
{
    size_t n = (size_t)pp->h * pp->N0;
    uint64_t *v = malloc(n * sizeof(uint64_t)), *slots = malloc((size_t)pp->N0 * sizeof(uint64_t));
    if (!v || !slots) abort();
    memset(in_S, 0, (size_t)pp->n);
    rlwe_dec(pp, v, w);
    int count = 0, nonzero = 0;
    const uint64_t bound = (uint64_t)1 << (pp->h - 1);
    for (int d = 0; d < pp->h && count >= 0; d++) {
        const uint64_t *vd = v + (size_t)d * pp->N0;
        int zero = 1;
        for (int c = 0; c < pp->N0; c++) if (vd[c]) { zero = 0; break; }
        if (zero) continue;
        nonzero = 1;
        rlwe_slot_decode(pp, slots, vd);
        for (int i = 0; i < pp->n; i++) {
            if (slots[i] > bound) { count = -1; break; }
            if (slots[i] && !in_S[i]) { in_S[i] = 1; count++; }
        }
    }
    if (!nonzero) count = -1;
    if (count < 0) memset(in_S, 0, (size_t)pp->n);
    free(v); free(slots);
    return count;
}
