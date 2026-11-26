/* OpenSSL first: params.h defines the macros N and H, which OpenSSL headers
   use as identifiers. */
#include <openssl/evp.h>
#include "eg_board.h"
#include "board_util.h"
#include "prg.h"
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

/* eg_board.c -> exponential ElGamal board proofs on GMP. */

/* Challenges lie in [0, 2^128) and add modulo 2^128.  Their differences are
   below q, hence invertible modulo q^h, which special soundness needs. */
#define CHAL_BITS  128
#define CHAL_BYTES (CHAL_BITS / 8)
/* Window of the per-entry tables of e1^{-1} and e2^{-1}. */
#define ENTRY_W    4
#define RESP_BYTES ((H * Q_BITS + 7) / 8)

/* ---- arithmetic ---------------------------------------------------------- */

static void from_bn(mpz_t out, const bn_t x)
{
    int n = bn_size_bin(x);
    uint8_t *buf = malloc((size_t)n + 1);
    if (!buf) abort();
    if (bn_is_zero(x)) { mpz_set_ui(out, 0); free(buf); return; }
    bn_write_bin(buf, (size_t)n, x);
    mpz_import(out, (size_t)n, 1, 1, 1, 0, buf);
    free(buf);
}

/* Uniform in [0, mod), up to a statistical distance of 2^-128. */
static void rand_mod(mpz_t x, const mpz_t mod)
{
    size_t len = (mpz_sizeinbase(mod, 2) + 7) / 8 + 16;
    uint8_t *buf = malloc(len);
    if (!buf) abort();
    rnd_bytes(buf, len);
    mpz_import(x, len, 1, 1, 1, 0, buf);
    mpz_mod(x, x, mod);
    free(buf);
}

static void fixb_init(EGFixedBase *fb, const mpz_t base, const mpz_t p, int w, size_t bits)
{
    fb->w = w;
    fb->nwin = (int)((bits + (size_t)w - 1) / (size_t)w);
    size_t per = ((size_t)1 << w) - 1;
    fb->t = malloc((size_t)fb->nwin * per * sizeof(mpz_t));
    if (!fb->t) abort();
    mpz_t b;
    mpz_init_set(b, base);
    for (int i = 0; i < fb->nwin; i++) {
        mpz_t *row = fb->t + (size_t)i * per;
        mpz_init_set(row[0], b);
        for (size_t d = 1; d < per; d++) {
            mpz_init(row[d]);
            mpz_mul(row[d], row[d - 1], b);
            mpz_mod(row[d], row[d], p);
        }
        /* b <- b^{2^w} = row[per-1] * b */
        mpz_mul(b, row[per - 1], b);
        mpz_mod(b, b, p);
    }
    mpz_clear(b);
}

static void fixb_free(EGFixedBase *fb)
{
    size_t per = ((size_t)1 << fb->w) - 1;
    for (size_t k = 0; k < (size_t)fb->nwin * per; k++) mpz_clear(fb->t[k]);
    free(fb->t);
}

/* out = out * base^e mod p for 0 <= e < 2^{w nwin}. */
static void fixb_mulpow(const EGFixedBase *fb, mpz_t out, const mpz_t e, const mpz_t p)
{
    size_t per = ((size_t)1 << fb->w) - 1;
    size_t nb = mpz_sizeinbase(e, 2);
    if ((size_t)fb->nwin * (size_t)fb->w < nb) abort();
    for (int i = 0; (size_t)i * (size_t)fb->w < nb; i++) {
        unsigned long d = 0;
        for (int k = fb->w - 1; k >= 0; k--)
            d = (d << 1) | (unsigned long)mpz_tstbit(e, (mp_bitcnt_t)(i * fb->w + k));
        if (d) {
            mpz_mul(out, out, fb->t[(size_t)i * per + d - 1]);
            mpz_mod(out, out, p);
        }
    }
}

static int below(const mpz_t x, const mpz_t bound)
{
    return mpz_sgn(x) >= 0 && mpz_cmp(x, bound) < 0;
}

/* ---- hashing ------------------------------------------------------------- */

static void put(EVP_MD_CTX *h, const mpz_t x, size_t len)
{
    uint8_t buf[P_BYTES];
    size_t n = (mpz_sizeinbase(x, 2) + 7) / 8;
    if (len > sizeof buf || n > len) abort();
    memset(buf, 0, len);
    if (mpz_sgn(x)) mpz_export(buf + len - n, NULL, 1, 1, 1, 0, x);
    EVP_DigestUpdate(h, buf, len);
}

static void put_u64(EVP_MD_CTX *h, uint64_t v)
{
    uint8_t b[8];
    for (int i = 0; i < 8; i++) b[i] = (uint8_t)(v >> (56 - 8 * i));
    EVP_DigestUpdate(h, b, 8);
}

static EVP_MD_CTX *chal_init(const EGBoardParams *bp, const char *label)
{
    EVP_MD_CTX *h = EVP_MD_CTX_new();
    if (!h || EVP_DigestInit_ex(h, EVP_sha512(), NULL) != 1) abort();
    EVP_DigestUpdate(h, label, strlen(label) + 1);
    EVP_DigestUpdate(h, bp->ctx, sizeof bp->ctx);
    return h;
}

static void chal_final(EVP_MD_CTX *h, mpz_t e, const mpz_t q)
{
    uint8_t d[64];
    unsigned int len = 0;
    EVP_DigestFinal_ex(h, d, &len);
    EVP_MD_CTX_free(h);
    mpz_import(e, 64, 1, 1, 1, 0, d);
    mpz_mod(e, e, q);
}

/* ---- parameters ---------------------------------------------------------- */

void eg_board_params_init(EGBoardParams *bp, const ECTParams *pp, int w)
{
    bp->n = pp->n;
    bp->h = pp->h;
    mpz_inits(bp->p, bp->q, bp->qh, bp->g, bp->pk, bp->cmod, NULL);
    mpz_setbit(bp->cmod, CHAL_BITS);
    mpz_init(bp->cof);
    from_bn(bp->p, pp->p);
    from_bn(bp->q, pp->q);
    from_bn(bp->qh, pp->qh);
    from_bn(bp->g, pp->g);
    from_bn(bp->pk, pp->pk);
    mpz_sub_ui(bp->cof, bp->p, 1);
    if (!mpz_divisible_p(bp->cof, bp->qh)) abort();
    mpz_divexact(bp->cof, bp->cof, bp->qh);
    bp->m = malloc((size_t)bp->n * sizeof(mpz_t));
    if (!bp->m) abort();
    for (int i = 0; i < bp->n; i++) {
        mpz_init(bp->m[i]);
        mpz_ui_pow_ui(bp->m[i], (unsigned long)pp->r, (unsigned long)i);
    }
    size_t bits = mpz_sizeinbase(bp->qh, 2);
    fixb_init(&bp->fg, bp->g, bp->p, w, bits);
    fixb_init(&bp->fpk, bp->pk, bp->p, w, bits);
    EVP_MD_CTX *h = EVP_MD_CTX_new();
    unsigned int len = 0;
    if (!h || EVP_DigestInit_ex(h, EVP_sha512(), NULL) != 1) abort();
    put(h, bp->p, P_BYTES); put(h, bp->g, P_BYTES); put(h, bp->pk, P_BYTES);
    put_u64(h, (uint64_t)bp->n); put_u64(h, (uint64_t)bp->h);
    EVP_DigestFinal_ex(h, bp->ctx, &len);
    EVP_MD_CTX_free(h);
    sig_keygen(&bp->sig, "ED25519");
}

void eg_board_params_free(EGBoardParams *bp)
{
    fixb_free(&bp->fg);
    fixb_free(&bp->fpk);
    for (int i = 0; i < bp->n; i++) mpz_clear(bp->m[i]);
    free(bp->m);
    mpz_clears(bp->p, bp->q, bp->qh, bp->g, bp->pk, bp->cmod, bp->cof, NULL);
    sig_free(&bp->sig);
}

void eg_entry_init(const EGBoardParams *bp, EGEntry *E)
{
    mpz_inits(E->e1, E->e2, E->C1, E->C2, E->c0, NULL);
    E->c = malloc((size_t)bp->n * sizeof(mpz_t));
    E->z = malloc((size_t)EG_NZ(bp->n) * sizeof(mpz_t));
    if (!E->c || !E->z) abort();
    for (int i = 0; i < bp->n; i++) mpz_init(E->c[i]);
    for (int i = 0; i < EG_NZ(bp->n); i++) mpz_init(E->z[i]);
    E->T = malloc((size_t)EG_NZ(bp->n) * sizeof(*E->T));
    if (!E->T) abort();
    for (int i = 0; i < EG_NZ(bp->n); i++) mpz_inits(E->T[i][0], E->T[i][1], NULL);
}

void eg_entry_free(const EGBoardParams *bp, EGEntry *E)
{
    mpz_clears(E->e1, E->e2, E->C1, E->C2, E->c0, NULL);
    for (int i = 0; i < bp->n; i++) mpz_clear(E->c[i]);
    for (int i = 0; i < EG_NZ(bp->n); i++) mpz_clear(E->z[i]);
    for (int i = 0; i < EG_NZ(bp->n); i++) mpz_clears(E->T[i][0], E->T[i][1], NULL);
    free(E->c); free(E->z); free(E->T);
}

void eg_board_init(EGBoard *B, const EGBoardParams *bp, size_t M)
{
    B->M = M;
    B->E = malloc(M * sizeof(EGEntry));
    B->sig = malloc(bp->sig.sig_bytes);
    if (!B->E || !B->sig) abort();
    for (size_t j = 0; j < M; j++) eg_entry_init(bp, &B->E[j]);
    B->siglen = 0;
    B->c = 0;
    mpz_inits(B->ce, B->cz, NULL);
}

void eg_board_free(EGBoard *B, const EGBoardParams *bp)
{
    for (size_t j = 0; j < B->M; j++) eg_entry_free(bp, &B->E[j]);
    free(B->E); free(B->sig);
    mpz_clears(B->ce, B->cz, NULL);
}

/* ---- the entry proof ----------------------------------------------------- */

/* Scratch of one entry proof: the inverses of the four group elements, the
   first messages T[k] = (T1, T2) of the EG_NZ(n) atoms, and temporaries. */
typedef struct {
    mpz_t e1i, e2i, C1i, C2i;
    EGFixedBase fe1, fe2;   /* tables of e1^{-1} and e2^{-1} */
    int has_tables;
    mpz_t (*T)[2];
    mpz_t t, u;
} Scratch;

static void scratch_init(const EGBoardParams *bp, Scratch *s)
{
    mpz_inits(s->e1i, s->e2i, s->C1i, s->C2i, s->t, s->u, NULL);
    s->T = malloc((size_t)EG_NZ(bp->n) * sizeof(*s->T));
    if (!s->T) abort();
    for (int k = 0; k < EG_NZ(bp->n); k++) mpz_inits(s->T[k][0], s->T[k][1], NULL);
    s->has_tables = 0;
}

static void scratch_free(const EGBoardParams *bp, Scratch *s)
{
    mpz_clears(s->e1i, s->e2i, s->C1i, s->C2i, s->t, s->u, NULL);
    for (int k = 0; k < EG_NZ(bp->n); k++) mpz_clears(s->T[k][0], s->T[k][1], NULL);
    free(s->T);
    if (s->has_tables) { fixb_free(&s->fe1); fixb_free(&s->fe2); }
}

/* out = x^c mod p, from the table of x if there is one. */
static void pow_short(const EGBoardParams *bp, mpz_t out, const EGFixedBase *fx,
                      const mpz_t x, const mpz_t c)
{
    if (fx) { mpz_set_ui(out, 1); fixb_mulpow(fx, out, c, bp->p); }
    else    mpz_powm(out, x, c, bp->p);
}

/* The atom (X, Y g^{-shift}) has the same discrete logarithm for g and pk.
   Its first message for challenge c and response z is
     T1 = g^z X^{-c},  T2 = pk^z Y^{-c} g^{shift c},
   computed from Xi = X^{-1} and Yi = Y^{-1}, or from their tables fx and fy
   when these are not NULL; shift is NULL for 0. */
static void atom_sim(const EGBoardParams *bp, mpz_t T[2],
                     const EGFixedBase *fx, const mpz_t Xi,
                     const EGFixedBase *fy, const mpz_t Yi,
                     const mpz_t shift, const mpz_t c, const mpz_t z, mpz_t t)
{
    pow_short(bp, T[0], fx, Xi, c);
    fixb_mulpow(&bp->fg, T[0], z, bp->p);
    pow_short(bp, T[1], fy, Yi, c);
    fixb_mulpow(&bp->fpk, T[1], z, bp->p);
    if (shift) {
        mpz_mul(t, shift, c);
        mpz_mod(t, t, bp->qh);
        fixb_mulpow(&bp->fg, T[1], t, bp->p);
    }
}

/* First message of a real atom: (g^k, pk^k). */
static void atom_commit(const EGBoardParams *bp, mpz_t T[2], const mpz_t k)
{
    mpz_set_ui(T[0], 1);
    fixb_mulpow(&bp->fg, T[0], k, bp->p);
    mpz_set_ui(T[1], 1);
    fixb_mulpow(&bp->fpk, T[1], k, bp->p);
}

static void respond(const EGBoardParams *bp, mpz_t z, const mpz_t k,
                    const mpz_t c, const mpz_t x)
{
    mpz_mul(z, c, x);
    mpz_add(z, z, k);
    mpz_mod(z, z, bp->qh);
}

static void entry_challenge(const EGBoardParams *bp, const EGEntry *E, uint64_t j,
                            mpz_t (*T)[2], mpz_t e)
{
    EVP_MD_CTX *h = chal_init(bp, "ECT-EG-board-entry");
    put_u64(h, j);
    put(h, E->e1, P_BYTES); put(h, E->e2, P_BYTES);
    put(h, E->C1, P_BYTES); put(h, E->C2, P_BYTES);
    for (int k = 0; k < EG_NZ(bp->n); k++) {
        put(h, T[k][0], P_BYTES);
        put(h, T[k][1], P_BYTES);
    }
    chal_final(h, e, bp->cmod);
}

static void inverses(const EGBoardParams *bp, const EGEntry *E, Scratch *s)
{
    if (!mpz_invert(s->e1i, E->e1, bp->p) || !mpz_invert(s->e2i, E->e2, bp->p) ||
        !mpz_invert(s->C1i, E->C1, bp->p) || !mpz_invert(s->C2i, E->C2, bp->p))
        abort();
    /* e1^{-1} and e2^{-1} are raised to n + 1 challenges each */
    fixb_init(&s->fe1, s->e1i, bp->p, ENTRY_W, CHAL_BITS);
    fixb_init(&s->fe2, s->e2i, bp->p, ENTRY_W, CHAL_BITS);
    s->has_tables = 1;
}

/* Atom indices: 0 = e enc 0, 1 = C enc 0, 2 = C enc 1, 3 + i = e enc m_{i+1}. */

void eg_entry_prove(const EGBoardParams *bp, EGEntry *E, uint64_t j, int slot,
                    int b, mpz_t tau_out)
{
    const int n = bp->n;
    Scratch s;
    scratch_init(bp, &s);
    mpz_t rho, tau, k1, k2, e, one, c1;
    mpz_inits(rho, tau, k1, k2, e, one, c1, NULL);
    mpz_set_ui(one, 1);
    int ia = b ? slot - 1 : -1;
    if (b && (slot < 1 || slot > n)) abort();

    /* e = (g^rho, pk^rho g^{b m}) and C = (g^tau, pk^tau g^b) */
    rand_mod(rho, bp->qh);
    rand_mod(tau, bp->qh);
    mpz_set_ui(E->e1, 1); fixb_mulpow(&bp->fg, E->e1, rho, bp->p);
    mpz_set_ui(E->e2, 1); fixb_mulpow(&bp->fpk, E->e2, rho, bp->p);
    if (b) fixb_mulpow(&bp->fg, E->e2, bp->m[ia], bp->p);
    mpz_set_ui(E->C1, 1); fixb_mulpow(&bp->fg, E->C1, tau, bp->p);
    mpz_set_ui(E->C2, 1); fixb_mulpow(&bp->fpk, E->C2, tau, bp->p);
    if (b) { mpz_mul(E->C2, E->C2, bp->g); mpz_mod(E->C2, E->C2, bp->p); }
    inverses(bp, E, &s);

    if (!b) {
        /* Right side simulated: every c[i] at random, c1 their sum. */
        mpz_set_ui(c1, 0);
        for (int i = 0; i < n; i++) {
            rand_mod(E->c[i], bp->cmod);
            mpz_add(c1, c1, E->c[i]);
            rand_mod(E->z[3 + i], bp->qh);
            atom_sim(bp, s.T[3 + i], &s.fe1, s.e1i, &s.fe2, s.e2i, bp->m[i], E->c[i], E->z[3 + i], s.t);
        }
        mpz_mod(c1, c1, bp->cmod);
        rand_mod(E->z[2], bp->qh);
        atom_sim(bp, s.T[2], NULL, s.C1i, NULL, s.C2i, one, c1, E->z[2], s.t);
        /* Left side real. */
        rand_mod(k1, bp->qh);
        rand_mod(k2, bp->qh);
        atom_commit(bp, s.T[0], k1);
        atom_commit(bp, s.T[1], k2);
        entry_challenge(bp, E, j, s.T, e);
        mpz_sub(E->c0, e, c1);
        mpz_mod(E->c0, E->c0, bp->cmod);
        respond(bp, E->z[0], k1, E->c0, rho);
        respond(bp, E->z[1], k2, E->c0, tau);
    } else {
        /* Left side simulated. */
        rand_mod(E->c0, bp->cmod);
        rand_mod(E->z[0], bp->qh);
        rand_mod(E->z[1], bp->qh);
        atom_sim(bp, s.T[0], &s.fe1, s.e1i, &s.fe2, s.e2i, NULL, E->c0, E->z[0], s.t);
        atom_sim(bp, s.T[1], NULL, s.C1i, NULL, s.C2i, NULL, E->c0, E->z[1], s.t);
        /* Right side: the inner "or" simulated except at ia. */
        mpz_set_ui(c1, 0);
        for (int i = 0; i < n; i++) {
            if (i == ia) continue;
            rand_mod(E->c[i], bp->cmod);
            mpz_add(c1, c1, E->c[i]);
            rand_mod(E->z[3 + i], bp->qh);
            atom_sim(bp, s.T[3 + i], &s.fe1, s.e1i, &s.fe2, s.e2i, bp->m[i], E->c[i], E->z[3 + i], s.t);
        }
        rand_mod(k1, bp->qh);
        rand_mod(k2, bp->qh);
        atom_commit(bp, s.T[2], k1);
        atom_commit(bp, s.T[3 + ia], k2);
        entry_challenge(bp, E, j, s.T, e);
        /* the right side takes e - c0, and c[ia] completes the inner sum */
        mpz_sub(e, e, E->c0);
        mpz_mod(e, e, bp->cmod);
        mpz_sub(E->c[ia], e, c1);
        mpz_mod(E->c[ia], E->c[ia], bp->cmod);
        respond(bp, E->z[2], k1, e, tau);
        respond(bp, E->z[3 + ia], k2, E->c[ia], rho);
    }

    for (int k = 0; k < EG_NZ(n); k++) {
        mpz_set(E->T[k][0], s.T[k][0]);
        mpz_set(E->T[k][1], s.T[k][1]);
    }
    if (tau_out) mpz_set(tau_out, tau);
    mpz_clears(rho, tau, k1, k2, e, one, c1, NULL);
    scratch_free(bp, &s);
}

/* out = prod_k x[k]^{e[k]} mod p for e[k] >= 0, by interleaved windows of
   width 3 (Straus). */
static void multiexp(mpz_t out, mpz_t *x, mpz_t *e, int k, const mpz_t p)
{
    enum { W = 3, PER = (1 << W) - 1 };
    mpz_t *t = malloc((size_t)k * PER * sizeof(mpz_t));
    if (!t) abort();
    size_t nb = 0;
    for (int i = 0; i < k; i++) {
        size_t b = mpz_sizeinbase(e[i], 2);
        if (b > nb) nb = b;
        mpz_init_set(t[i * PER], x[i]);
        for (int d = 1; d < PER; d++) {
            mpz_init(t[i * PER + d]);
            mpz_mul(t[i * PER + d], t[i * PER + d - 1], x[i]);
            mpz_mod(t[i * PER + d], t[i * PER + d], p);
        }
    }
    int nwin = (int)((nb + W - 1) / W);
    mpz_set_ui(out, 1);
    for (int w = nwin - 1; w >= 0; w--) {
        for (int q = 0; q < W && w != nwin - 1; q++) {
            mpz_mul(out, out, out);
            mpz_mod(out, out, p);
        }
        for (int i = 0; i < k; i++) {
            unsigned d = 0;
            for (int q = W - 1; q >= 0; q--)
                d = (d << 1) | (unsigned)mpz_tstbit(e[i], (mp_bitcnt_t)(w * W + q));
            if (d) {
                mpz_mul(out, out, t[i * PER + d - 1]);
                mpz_mod(out, out, p);
            }
        }
    }
    for (int i = 0; i < k * PER; i++) mpz_clear(t[i]);
    free(t);
}

static void rand_bits(mpz_t x, int bits)
{
    uint8_t buf[64];
    size_t len = (size_t)(bits + 7) / 8;
    if (len > sizeof buf) abort();
    rnd_bytes(buf, len);
    mpz_import(x, len, 1, 1, 1, 0, buf);
    mpz_fdiv_r_2exp(x, x, (mp_bitcnt_t)bits);
}

void eg_batch_init(EGBatch *Bt)
{
    mpz_inits(Bt->GE, Bt->PE, Bt->R, NULL);
    mpz_set_ui(Bt->R, 1);
    for (int r = 0; r < EG_MEM_REPS; r++) mpz_init_set_ui(Bt->mem[r], 1);
    Bt->bad = 0;
}

void eg_batch_free(EGBatch *Bt)
{
    mpz_clears(Bt->GE, Bt->PE, Bt->R, NULL);
    for (int r = 0; r < EG_MEM_REPS; r++) mpz_clear(Bt->mem[r]);
}

static int unit(const EGBoardParams *bp, const mpz_t x)
{
    return mpz_sgn(x) > 0 && mpz_cmp(x, bp->p) < 0;
}

/* The atom k with statement (X, Y g^{-shift}), challenge c, response z and
   first message (T1, T2) has the equations g^z = T1 X^c and
   pk^z g^{shift c} = T2 Y^c.  With coefficients l1 and l2, the batch gains
   GE += l1 z + l2 shift c and PE += l2 z, and R gains T1^l1 T2^l2 X^{l1 c}
   Y^{l2 c}; the powers of X and Y are collected in xe[0..3] for (e1, e2,
   C1, C2). */
static void add_atom(EGBatch *Bt, mpz_t *xe, int xi, int yi,
                     const mpz_t shift, const mpz_t c, const mpz_t z,
                     mpz_t l1, mpz_t l2, mpz_t t)
{
    rand_bits(l1, 128);
    rand_bits(l2, 128);
    mpz_addmul(Bt->GE, l1, z);
    mpz_addmul(Bt->PE, l2, z);
    if (shift) {
        mpz_mul(t, l2, shift);
        mpz_addmul(Bt->GE, t, c);
    }
    mpz_addmul(xe[xi], l1, c);
    mpz_addmul(xe[yi], l2, c);
}

void eg_batch_add(const EGBoardParams *bp, EGBatch *Bt, const EGEntry *E, uint64_t j)
{
    const int n = bp->n, K = EG_NZ(n);
    if (Bt->bad) return;
    if (!unit(bp, E->e1) || !unit(bp, E->e2) || !unit(bp, E->C1) || !unit(bp, E->C2) ||
        !below(E->c0, bp->cmod)) { Bt->bad = 1; return; }
    for (int i = 0; i < n; i++) if (!below(E->c[i], bp->cmod)) { Bt->bad = 1; return; }
    for (int k = 0; k < K; k++)
        if (!below(E->z[k], bp->qh) || !unit(bp, E->T[k][0]) || !unit(bp, E->T[k][1])) {
            Bt->bad = 1; return;
        }

    /* the challenges add up to the hash */
    mpz_t e, c1, one, l1, l2, t;
    mpz_inits(e, c1, one, l1, l2, t, NULL);
    mpz_set_ui(one, 1);
    for (int i = 0; i < n; i++) mpz_add(c1, c1, E->c[i]);
    mpz_mod(c1, c1, bp->cmod);
    entry_challenge(bp, E, j, E->T, e);
    mpz_add(t, c1, E->c0);
    mpz_mod(t, t, bp->cmod);
    if (mpz_cmp(t, e) != 0) { Bt->bad = 1; goto done; }

    /* bases: the 2K first messages, then e1, e2, C1, C2 */
    {
        mpz_t *x = malloc((size_t)(2 * K + 4) * sizeof(mpz_t));
        mpz_t *ex = malloc((size_t)(2 * K + 4) * sizeof(mpz_t));
        if (!x || !ex) abort();
        for (int i = 0; i < 2 * K + 4; i++) mpz_inits(x[i], ex[i], NULL);
        mpz_t *xe = ex + 2 * K;
        mpz_set(x[2 * K + 0], E->e1); mpz_set(x[2 * K + 1], E->e2);
        mpz_set(x[2 * K + 2], E->C1); mpz_set(x[2 * K + 3], E->C2);
        for (int k = 0; k < K; k++) {
            mpz_srcptr c = k < 2 ? E->c0 : k == 2 ? c1 : E->c[k - 3];
            int xi = (k == 1 || k == 2) ? 2 : 0;
            add_atom(Bt, xe, xi, xi + 1, k == 2 ? one : k >= 3 ? bp->m[k - 3] : NULL,
                     c, E->z[k], l1, l2, t);
            mpz_set(x[2 * k], E->T[k][0]);   mpz_set(ex[2 * k], l1);
            mpz_set(x[2 * k + 1], E->T[k][1]); mpz_set(ex[2 * k + 1], l2);
        }
        multiexp(t, x, ex, 2 * K + 4, bp->p);
        mpz_mul(Bt->R, Bt->R, t);
        mpz_mod(Bt->R, Bt->R, bp->p);
        mpz_mod(Bt->GE, Bt->GE, bp->qh);
        mpz_mod(Bt->PE, Bt->PE, bp->qh);
        for (int i = 0; i < 2 * K + 4; i++) mpz_clears(x[i], ex[i], NULL);
        free(x); free(ex);
    }

    /* membership: each of the four elements joins a random subset of the
       EG_MEM_REPS products */
    {
        const mpz_t *el[4] = { &E->e1, &E->e2, &E->C1, &E->C2 };
        uint8_t sel[EG_MEM_REPS / 8];
        for (int q = 0; q < 4; q++) {
            rnd_bytes(sel, sizeof sel);
            for (int r = 0; r < EG_MEM_REPS; r++)
                if (sel[r / 8] >> (r % 8) & 1) {
                    mpz_mul(Bt->mem[r], Bt->mem[r], *el[q]);
                    mpz_mod(Bt->mem[r], Bt->mem[r], bp->p);
                }
        }
    }
done:
    mpz_clears(e, c1, one, l1, l2, t, NULL);
}

void eg_batch_merge(const EGBoardParams *bp, EGBatch *dst, const EGBatch *src)
{
    dst->bad |= src->bad;
    mpz_add(dst->GE, dst->GE, src->GE); mpz_mod(dst->GE, dst->GE, bp->qh);
    mpz_add(dst->PE, dst->PE, src->PE); mpz_mod(dst->PE, dst->PE, bp->qh);
    mpz_mul(dst->R, dst->R, src->R);    mpz_mod(dst->R, dst->R, bp->p);
    for (int r = 0; r < EG_MEM_REPS; r++) {
        mpz_mul(dst->mem[r], dst->mem[r], src->mem[r]);
        mpz_mod(dst->mem[r], dst->mem[r], bp->p);
    }
}

int eg_batch_check(const EGBoardParams *bp, const EGBatch *Bt)
{
    if (Bt->bad) return 0;
    int ok = 1;
    mpz_t L, t;
    mpz_inits(L, t, NULL);
    mpz_set_ui(L, 1);
    fixb_mulpow(&bp->fg, L, Bt->GE, bp->p);
    fixb_mulpow(&bp->fpk, L, Bt->PE, bp->p);
    if (!mpz_invert(t, Bt->R, bp->p)) ok = 0;
    if (ok) {
        mpz_mul(L, L, t);
        mpz_mod(L, L, bp->p);
        mpz_powm(L, L, bp->cof, bp->p);
        ok = mpz_cmp_ui(L, 1) == 0;
    }
    for (int r = 0; ok && r < EG_MEM_REPS; r++) {
        mpz_powm(t, Bt->mem[r], bp->qh, bp->p);
        ok = mpz_cmp_ui(t, 1) == 0;
    }
    mpz_clears(L, t, NULL);
    return ok;
}

int eg_entry_verify(const EGBoardParams *bp, const EGEntry *E, uint64_t j)
{
    EGBatch Bt;
    eg_batch_init(&Bt);
    eg_batch_add(bp, &Bt, E, j);
    int ok = eg_batch_check(bp, &Bt);
    eg_batch_free(&Bt);
    return ok;
}

void eg_entry_leaf(const EGEntry *E, uint8_t leaf[32])
{
    EVP_MD_CTX *h = EVP_MD_CTX_new();
    unsigned int len = 0;
    uint8_t d[32];
    static const uint8_t tag_d = 0x02, tag_leaf = 0x00;
    if (!h || EVP_DigestInit_ex(h, EVP_sha256(), NULL) != 1) abort();
    EVP_DigestUpdate(h, &tag_d, 1);
    put(h, E->e1, P_BYTES);
    put(h, E->e2, P_BYTES);
    EVP_DigestFinal_ex(h, d, &len);
    if (EVP_DigestInit_ex(h, EVP_sha256(), NULL) != 1) abort();
    EVP_DigestUpdate(h, &tag_leaf, 1);
    put(h, E->C1, P_BYTES);
    put(h, E->C2, P_BYTES);
    EVP_DigestUpdate(h, d, sizeof d);
    EVP_DigestFinal_ex(h, leaf, &len);
    EVP_MD_CTX_free(h);
}

/* ---- the count proof ----------------------------------------------------- */

static void count_challenge(const EGBoardParams *bp, const mpz_t X, const mpz_t Y,
                            uint64_t M, unsigned long c, mpz_t T[2], mpz_t e)
{
    EVP_MD_CTX *h = chal_init(bp, "ECT-EG-board-count");
    put_u64(h, M);
    put_u64(h, (uint64_t)c);
    put(h, X, P_BYTES); put(h, Y, P_BYTES);
    put(h, T[0], P_BYTES); put(h, T[1], P_BYTES);
    chal_final(h, e, bp->cmod);
}

/* (X, Y g^{-c}) has the same discrete logarithm tau for g and pk. */
void eg_count_prove(const EGBoardParams *bp, const mpz_t X, const mpz_t Y,
                    uint64_t M, unsigned long c, const mpz_t tau,
                    mpz_t ce, mpz_t cz)
{
    mpz_t k, T[2];
    mpz_inits(k, T[0], T[1], NULL);
    rand_mod(k, bp->qh);
    atom_commit(bp, T, k);
    count_challenge(bp, X, Y, M, c, T, ce);
    respond(bp, cz, k, ce, tau);
    mpz_clears(k, T[0], T[1], NULL);
}

int eg_count_verify(const EGBoardParams *bp, const mpz_t X, const mpz_t Y,
                    uint64_t M, unsigned long c, const mpz_t ce, const mpz_t cz)
{
    int ok = 0;
    mpz_t Xi, Yi, sh, e, t, T[2];
    mpz_inits(Xi, Yi, sh, e, t, T[0], T[1], NULL);
    if (c > M || !below(ce, bp->cmod) || !below(cz, bp->qh)) goto done;
    if (!mpz_invert(Xi, X, bp->p) || !mpz_invert(Yi, Y, bp->p)) goto done;
    mpz_set_ui(sh, c);
    atom_sim(bp, T, NULL, Xi, NULL, Yi, sh, ce, cz, t);
    count_challenge(bp, X, Y, M, c, T, e);
    ok = mpz_cmp(e, ce) == 0;
done:
    mpz_clears(Xi, Yi, sh, e, t, T[0], T[1], NULL);
    return ok;
}

/* ---- Board and VerifyBoard ----------------------------------------------- */

typedef struct {
    const EGBoardParams *bp;
    EGBoard *B;
    const EGBoard *Bc;
    const int *slot, *b;
    mpz_t *tau;          /* per entry, for the count proof */
    uint8_t (*leaf)[32];
    EGBatch *bt;         /* one batch per thread */
    atomic_int bad;
} Job;

static void prove_range(void *ctx, int tid, size_t lo, size_t hi)
{
    (void)tid;
    Job *J = ctx;
    for (size_t j = lo; j < hi; j++) {
        eg_entry_prove(J->bp, &J->B->E[j], j, J->slot[j], J->b[j], J->tau[j]);
        eg_entry_leaf(&J->B->E[j], J->leaf[j]);
    }
}

static void verify_range(void *ctx, int tid, size_t lo, size_t hi)
{
    Job *J = ctx;
    for (size_t j = lo; j < hi; j++) {
        eg_batch_add(J->bp, &J->bt[tid], &J->Bc->E[j], j);
        eg_entry_leaf(&J->Bc->E[j], J->leaf[j]);
    }
}

static void products(const EGBoardParams *bp, const EGBoard *B, mpz_t X, mpz_t Y)
{
    mpz_set_ui(X, 1);
    mpz_set_ui(Y, 1);
    for (size_t j = 0; j < B->M; j++) {
        mpz_mul(X, X, B->E[j].C1); mpz_mod(X, X, bp->p);
        mpz_mul(Y, Y, B->E[j].C2); mpz_mod(Y, Y, bp->p);
    }
}

void eg_board(const EGBoardParams *bp, EGBoard *B, const int *slot, const int *b,
              int threads)
{
    Job J = { bp, B, B, slot, b, malloc(B->M * sizeof(mpz_t)),
              malloc(B->M * 32), NULL, 0 };
    if (!J.tau || !J.leaf) abort();
    for (size_t j = 0; j < B->M; j++) mpz_init(J.tau[j]);
    board_parallel(threads, B->M, prove_range, &J);

    merkle_root((const uint8_t (*)[32])J.leaf, B->M, B->root);
    uint8_t msg[BOARD_MSG_BYTES];
    board_message(msg, B->M, B->root);
    B->siglen = sig_sign(&bp->sig, msg, sizeof msg, B->sig, bp->sig.sig_bytes);

    mpz_t X, Y, tau;
    mpz_inits(X, Y, tau, NULL);
    B->c = 0;
    for (size_t j = 0; j < B->M; j++) {
        mpz_add(tau, tau, J.tau[j]);
        B->c += (unsigned long)(b[j] != 0);
        mpz_clear(J.tau[j]);
    }
    mpz_mod(tau, tau, bp->qh);
    products(bp, B, X, Y);
    eg_count_prove(bp, X, Y, B->M, B->c, tau, B->ce, B->cz);
    mpz_clears(X, Y, tau, NULL);
    free(J.tau); free(J.leaf);
}

int eg_verify_board(const EGBoardParams *bp, const EGBoard *B, int threads)
{
    if (B->M == 0 || B->c > B->M) return 0;
    int nt = threads < 1 ? 1 : threads;
    Job J = { bp, NULL, B, NULL, NULL, NULL, malloc(B->M * 32),
              malloc((size_t)nt * sizeof(EGBatch)), 0 };
    if (!J.leaf || !J.bt) abort();
    for (int t = 0; t < nt; t++) eg_batch_init(&J.bt[t]);
    board_parallel(threads, B->M, verify_range, &J);
    for (int t = 1; t < nt; t++) eg_batch_merge(bp, &J.bt[0], &J.bt[t]);
    int ok = eg_batch_check(bp, &J.bt[0]);
    for (int t = 0; t < nt; t++) eg_batch_free(&J.bt[t]);
    free(J.bt);
    if (ok) {
        uint8_t root[32], msg[BOARD_MSG_BYTES];
        merkle_root((const uint8_t (*)[32])J.leaf, B->M, root);
        board_message(msg, B->M, root);
        ok = memcmp(root, B->root, 32) == 0 &&
             sig_verify(&bp->sig, msg, sizeof msg, B->sig, B->siglen);
    }
    if (ok) {
        mpz_t X, Y;
        mpz_inits(X, Y, NULL);
        products(bp, B, X, Y);
        ok = eg_count_verify(bp, X, Y, B->M, B->c, B->ce, B->cz);
        mpz_clears(X, Y, NULL);
    }
    free(J.leaf);
    return ok;
}

size_t eg_entry_bytes(const EGBoardParams *bp)
{
    /* e_j and C_j (four elements), n + 1 challenges of 16 bytes, n + 3
       responses mod q^h and n + 3 first messages of two elements. */
    return 4u * P_BYTES + (size_t)(bp->n + 1) * CHAL_BYTES
         + (size_t)EG_NZ(bp->n) * (RESP_BYTES + 2u * P_BYTES);
}

size_t eg_count_bytes(const EGBoardParams *bp)
{
    (void)bp;
    return 8u + CHAL_BYTES + RESP_BYTES;
}

size_t eg_root_bytes(const EGBoardParams *bp)
{
    return 32u + bp->sig.sig_bytes;
}
