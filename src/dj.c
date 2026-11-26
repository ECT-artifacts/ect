#include "dj.h"
#include "prg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* dj.c -> Damgard-Jurik instantiation of ECT on GMP. */

static void rand_bits(mpz_t x, int bits)
{
    size_t len = (size_t)(bits + 7) / 8;
    uint8_t *buf = malloc(len);
    if (!buf) abort();
    rnd_bytes(buf, len);
    mpz_import(x, len, 1, 1, 1, 0, buf);
    free(buf);
    mpz_fdiv_r_2exp(x, x, (mp_bitcnt_t)bits);
}

/* Prime of exactly `bits` bits with its two top bits set, so that the
   product of two such primes has exactly 2*bits bits. */
static void rand_prime(mpz_t p, int bits)
{
    rand_bits(p, bits);
    mpz_setbit(p, (mp_bitcnt_t)(bits - 1));
    mpz_setbit(p, (mp_bitcnt_t)(bits - 2));
    mpz_nextprime(p, p);
    if ((int)mpz_sizeinbase(p, 2) != bits) rand_prime(p, bits);
}

void dj_setup(DJParams *pp, int h, int nbits)
{
    pp->h = h;
    pp->nbits = nbits;
    mpz_inits(pp->N, pp->r, pp->rn, pp->p1, pp->p2, pp->lambda,
              pp->lambda_inv, NULL);

    do {
        rand_prime(pp->p1, nbits / 2);
        rand_prime(pp->p2, nbits / 2);
    } while (mpz_cmp(pp->p1, pp->p2) == 0);
    mpz_mul(pp->N, pp->p1, pp->p2);
    if ((int)mpz_sizeinbase(pp->N, 2) != nbits) abort();

    mpz_t a, b, g;
    mpz_inits(a, b, g, NULL);
    mpz_sub_ui(a, pp->p1, 1);
    mpz_sub_ui(b, pp->p2, 1);
    mpz_lcm(pp->lambda, a, b);
    mpz_mul(a, a, b);                  /* phi(N) */
    mpz_gcd(g, pp->N, a);
    if (mpz_cmp_ui(g, 1) != 0) abort(); /* gcd(N, phi(N)) = 1 */

    pp->Np = malloc((size_t)(h + 2) * sizeof(mpz_t));
    pp->inv_fact = malloc((size_t)(h + 1) * sizeof(mpz_t));
    if (!pp->Np || !pp->inv_fact) abort();
    mpz_init_set_ui(pp->Np[0], 1);
    for (int k = 1; k <= h + 1; k++) {
        mpz_init(pp->Np[k]);
        mpz_mul(pp->Np[k], pp->Np[k - 1], pp->N);
    }
    pp->P1p = malloc((size_t)(h + 2) * sizeof(mpz_t));
    pp->P2p = malloc((size_t)(h + 2) * sizeof(mpz_t));
    if (!pp->P1p || !pp->P2p) abort();
    mpz_init_set_ui(pp->P1p[0], 1);
    mpz_init_set_ui(pp->P2p[0], 1);
    for (int k = 1; k <= h + 1; k++) {
        mpz_inits(pp->P1p[k], pp->P2p[k], NULL);
        mpz_mul(pp->P1p[k], pp->P1p[k - 1], pp->p1);
        mpz_mul(pp->P2p[k], pp->P2p[k - 1], pp->p2);
    }
    mpz_set_ui(a, 1);
    for (int k = 0; k <= h; k++) {
        if (k > 0) mpz_mul_ui(a, a, (unsigned long)k);
        mpz_init(pp->inv_fact[k]);
        if (!mpz_invert(pp->inv_fact[k], a, pp->Np[h + 1])) abort();
    }
    if (!mpz_invert(pp->lambda_inv, pp->lambda, pp->Np[h])) abort();

    /* r = 2^(h-1) + 1, the smallest base above the multiplicity bound. */
    mpz_set_ui(pp->r, 0);
    mpz_setbit(pp->r, (mp_bitcnt_t)(h - 1));
    mpz_add_ui(pp->r, pp->r, 1);
    /* n = max{k : r^k < N} */
    pp->n = 0;
    mpz_set_ui(pp->rn, 1);
    for (;;) {
        mpz_mul(a, pp->rn, pp->r);
        if (mpz_cmp(a, pp->N) >= 0) break;
        mpz_set(pp->rn, a);
        pp->n++;
    }
    pp->elem_bytes = (size_t)(nbits + 7) / 8;
    mpz_clears(a, b, g, NULL);

    sig_keygen(&pp->sig, "ED25519");
}

void dj_free(DJParams *pp)
{
    for (int k = 0; k <= pp->h + 1; k++) mpz_clear(pp->Np[k]);
    for (int k = 0; k <= pp->h + 1; k++) mpz_clears(pp->P1p[k], pp->P2p[k], NULL);
    free(pp->P1p);
    free(pp->P2p);
    for (int k = 0; k <= pp->h; k++) mpz_clear(pp->inv_fact[k]);
    free(pp->Np);
    free(pp->inv_fact);
    mpz_clears(pp->N, pp->r, pp->rn, pp->p1, pp->p2, pp->lambda,
               pp->lambda_inv, NULL);
    sig_free(&pp->sig);
}

void dj_rand_unit(const DJParams *pp, mpz_t x)
{
    mpz_t g;
    mpz_init(g);
    for (;;) {
        rand_bits(x, pp->nbits + 64);
        mpz_mod(x, x, pp->N);
        mpz_gcd(g, x, pp->N);
        if (mpz_cmp_ui(g, 1) == 0) break;
    }
    mpz_clear(g);
}

void dj_pow_Ns(const DJParams *pp, mpz_t out, const mpz_t x, int s)
{
    /* If y = x^{N^{k-1}} mod N^k, then y^N = x^{N^k} mod N^{k+1}. */
    mpz_mod(out, x, pp->N);
    for (int k = 1; k <= s; k++)
        mpz_powm(out, out, pp->N, pp->Np[k + 1]);
}

/* out = the root of X^{p-1} = 1 modulo p^{s+1} that is congruent to
   x^{N^s} modulo p, where Pp[k] = p^k. */
static void root_of_unity(const DJParams *pp, mpz_t out, const mpz_t x,
                          const mpz_t p, mpz_t *Pp, int s)
{
    mpz_t pm1, ex, c, inv;
    mpz_inits(pm1, ex, c, inv, NULL);
    mpz_sub_ui(pm1, p, 1);
    /* x^{N^s} mod p, with the exponent reduced modulo p-1 */
    mpz_powm_ui(ex, pp->N, (unsigned long)s, pm1);
    mpz_powm(out, x, ex, p);
    /* Newton iteration for f(X) = X^{p-1} - 1: if out^{p-1} = 1 mod p^k,
       then out - (out^{p-1} - 1) out / (p-1) is a root mod p^{2k}. */
    for (int prec = 1; prec < s + 1; ) {
        prec = 2 * prec < s + 1 ? 2 * prec : s + 1;
        mpz_powm(c, out, pm1, Pp[prec]);
        mpz_sub_ui(c, c, 1);
        if (!mpz_invert(inv, pm1, Pp[prec])) abort();
        mpz_mul(c, c, out);
        mpz_mod(c, c, Pp[prec]);
        mpz_mul(c, c, inv);
        mpz_sub(out, out, c);
        mpz_mod(out, out, Pp[prec]);
    }
    mpz_clears(pm1, ex, c, inv, NULL);
}

void dj_pow_Ns_sk(const DJParams *pp, mpz_t out, const mpz_t x, int s)
{
    mpz_t a1, a2, inv;
    mpz_inits(a1, a2, inv, NULL);
    root_of_unity(pp, a1, x, pp->p1, pp->P1p, s);
    root_of_unity(pp, a2, x, pp->p2, pp->P2p, s);
    /* out = a1 + p1^{s+1} ((a2 - a1) / p1^{s+1} mod p2^{s+1}) */
    if (!mpz_invert(inv, pp->P1p[s + 1], pp->P2p[s + 1])) abort();
    mpz_sub(a2, a2, a1);
    mpz_mul(a2, a2, inv);
    mpz_mod(a2, a2, pp->P2p[s + 1]);
    mpz_mul(a2, a2, pp->P1p[s + 1]);
    mpz_add(out, a1, a2);
    mpz_clears(a1, a2, inv, NULL);
}

void dj_one_plus_N_pow(const DJParams *pp, mpz_t out, const mpz_t m, int s)
{
    /* (1+N)^m = sum_{k=0}^{s} binom(m,k) N^k mod N^{s+1}. */
    mpz_srcptr mod = pp->Np[s + 1];
    mpz_t fall, term, mk;
    mpz_inits(fall, term, mk, NULL);
    mpz_set_ui(out, 1);
    mpz_set_ui(fall, 1);
    for (int k = 1; k <= s; k++) {
        mpz_sub_ui(mk, m, (unsigned long)(k - 1));
        mpz_mul(fall, fall, mk);
        mpz_mod(fall, fall, mod);
        if (mpz_sgn(fall) == 0) break;     /* binom(m,k) = 0 for k > m */
        mpz_mul(term, fall, pp->inv_fact[k]);
        mpz_mod(term, term, pp->Np[s + 1 - k]);
        mpz_mul(term, term, pp->Np[k]);
        mpz_add(out, out, term);
    }
    mpz_mod(out, out, mod);
    mpz_clears(fall, term, mk, NULL);
}

static void enc(const DJParams *pp, mpz_t c, const mpz_t m, int s,
                mpz_t rho_out, int sk)
{
    mpz_t rho, t;
    mpz_inits(rho, t, NULL);
    dj_rand_unit(pp, rho);
    if (sk) dj_pow_Ns_sk(pp, t, rho, s);
    else    dj_pow_Ns(pp, t, rho, s);
    dj_one_plus_N_pow(pp, c, m, s);
    mpz_mul(c, c, t);
    mpz_mod(c, c, pp->Np[s + 1]);
    if (rho_out) mpz_set(rho_out, rho);
    mpz_clears(rho, t, NULL);
}

void dj_enc(const DJParams *pp, mpz_t c, const mpz_t m, int s, mpz_t rho_out)
{
    enc(pp, c, m, s, rho_out, 0);
}

void dj_enc_sk(const DJParams *pp, mpz_t c, const mpz_t m, int s, mpz_t rho_out)
{
    enc(pp, c, m, s, rho_out, 1);
}

/* Given a = (1+N)^i mod N^{s+1}, return i mod N^s (Damgard and Jurik,
   BRICS RS-00-45, Section 3). */
static void dj_extract(const DJParams *pp, mpz_t i, const mpz_t a, int s)
{
    mpz_t t1, t2, tmp;
    mpz_inits(t1, t2, tmp, NULL);
    mpz_set_ui(i, 0);
    for (int j = 1; j <= s; j++) {
        mpz_mod(t1, a, pp->Np[j + 1]);
        mpz_sub_ui(t1, t1, 1);
        mpz_divexact(t1, t1, pp->N);          /* L(a mod N^{j+1}) */
        mpz_set(t2, i);
        for (int k = 2; k <= j; k++) {
            mpz_sub_ui(i, i, 1);
            mpz_mul(t2, t2, i);
            mpz_mod(t2, t2, pp->Np[j]);
            mpz_mul(tmp, t2, pp->Np[k - 1]);
            mpz_mul(tmp, tmp, pp->inv_fact[k]);
            mpz_sub(t1, t1, tmp);
            mpz_mod(t1, t1, pp->Np[j]);
        }
        mpz_set(i, t1);
    }
    mpz_clears(t1, t2, tmp, NULL);
}

void dj_dec(const DJParams *pp, mpz_t v, const mpz_t c, int s)
{
    mpz_t a;
    mpz_init(a);
    mpz_powm(a, c, pp->lambda, pp->Np[s + 1]);   /* (1+N)^{v lambda mod N^s} */
    dj_extract(pp, v, a, s);
    mpz_mul(v, v, pp->lambda_inv);
    mpz_mod(v, v, pp->Np[s]);
    mpz_clear(a);
}

size_t dj_elem_len(const DJParams *pp, int s)
{
    return (size_t)(s + 1) * pp->elem_bytes;
}

void dj_write(const DJParams *pp, uint8_t *buf, const mpz_t x, int s)
{
    size_t len = dj_elem_len(pp, s), count = 0;
    size_t used = (mpz_sizeinbase(x, 2) + 7) / 8;
    if (mpz_sgn(x) == 0) used = 0;
    if (used > len) abort();
    memset(buf, 0, len - used);
    if (used) mpz_export(buf + (len - used), &count, 1, 1, 1, 0, x);
}

size_t dj_tag(const DJParams *pp, mpz_t w, uint8_t *sig, int i, int b, int hprime)
{
    mpz_t m;
    mpz_init(m);
    if (b) {
        mpz_pow_ui(m, pp->r, (unsigned long)(i - 1));
        mpz_mul(m, m, pp->Np[pp->h - hprime]);
    }
    dj_enc_sk(pp, w, m, pp->h, NULL);
    size_t len = dj_elem_len(pp, pp->h);
    uint8_t *buf = malloc(len);
    if (!buf) abort();
    dj_write(pp, buf, w, pp->h);
    size_t siglen = sig_sign(&pp->sig, buf, len, sig, pp->sig.sig_bytes);
    free(buf);
    mpz_clear(m);
    return siglen;
}

int dj_verify_tag(const DJParams *pp, const mpz_t w, const uint8_t *sig, size_t siglen)
{
    size_t len = dj_elem_len(pp, pp->h);
    uint8_t *buf = malloc(len);
    if (!buf) abort();
    dj_write(pp, buf, w, pp->h);
    int ok = sig_verify(&pp->sig, buf, len, sig, siglen);
    free(buf);
    return ok;
}

void dj_rerandomizer(const DJParams *pp, mpz_t t, mpz_t xi)
{
    dj_rand_unit(pp, xi);
    dj_pow_Ns(pp, t, xi, pp->h);
}

void dj_degrade_pre(const DJParams *pp, mpz_t w_out, const mpz_t t, const mpz_t w)
{
    mpz_powm(w_out, w, pp->N, pp->Np[pp->h + 1]);
    mpz_mul(w_out, w_out, t);
    mpz_mod(w_out, w_out, pp->Np[pp->h + 1]);
}

void dj_merge_pre(const DJParams *pp, mpz_t w_out, const mpz_t t, mpz_t *ws, int m)
{
    mpz_t acc;
    mpz_init_set(acc, t);
    for (int j = 0; j < m; j++) {
        mpz_mul(acc, acc, ws[j]);
        mpz_mod(acc, acc, pp->Np[pp->h + 1]);
    }
    mpz_set(w_out, acc);
    mpz_clear(acc);
}

void dj_degrade(const DJParams *pp, mpz_t w_out, mpz_t xi, const mpz_t w)
{
    mpz_t t;
    mpz_init(t);
    dj_rerandomizer(pp, t, xi);
    dj_degrade_pre(pp, w_out, t, w);
    mpz_clear(t);
}

void dj_merge(const DJParams *pp, mpz_t w_out, mpz_t xi, mpz_t *ws, int m)
{
    mpz_t t;
    mpz_init(t);
    dj_rerandomizer(pp, t, xi);
    dj_merge_pre(pp, w_out, t, ws, m);
    mpz_clear(t);
}

int dj_trace(const DJParams *pp, const mpz_t w, uint8_t *in_S)
{
    int count = 0;
    memset(in_S, 0, (size_t)pp->n);
    mpz_t v, d, c;
    mpz_inits(v, d, c, NULL);
    dj_dec(pp, v, w, pp->h);
    if (mpz_sgn(v) == 0) { count = -1; goto done; }
    /* First pass: every base-N digit must be below r^n. */
    for (int t = 0; t < pp->h; t++) {
        mpz_tdiv_qr(v, d, v, pp->N);
        if (mpz_cmp(d, pp->rn) >= 0) { count = -1; goto done; }
        /* Second step for this depth: base-r digits are identifiers. */
        for (int i = 0; i < pp->n && mpz_sgn(d) != 0; i++) {
            mpz_tdiv_qr(d, c, d, pp->r);
            if (mpz_sgn(c) != 0 && !in_S[i]) { in_S[i] = 1; count++; }
        }
    }
done:
    mpz_clears(v, d, c, NULL);
    return count;
}
