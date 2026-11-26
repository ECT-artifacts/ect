#pragma once
#include <gmp.h>
#include <stddef.h>
#include <stdint.h>
#include "sig.h"

/* dj.h -> Damgard-Jurik instantiation of ECT (Section "Damgard-Jurik
   instantiation"): tags in Z*_{N^{h+1}}, plaintexts in Z_{N^h},
   degradation w -> w^N rho^{N^h}, merging by multiplication followed by
   rerandomization, tracing by decryption and base-N / base-r digits. */

typedef struct {
    int    h;          /* hop budget; the scheme uses s = h */
    int    nbits;      /* bit length of N */
    int    n;          /* identifier capacity floor(log_r N), r^n < N */
    mpz_t  N, r, rn;   /* modulus, base r = 2^(h-1)+1, and r^n */
    mpz_t *Np;         /* Np[k] = N^k for k = 0..h+1 */
    size_t elem_bytes; /* bytes of N; an element of Z_{N^{s+1}} takes (s+1) of these */
    /* authority secrets */
    mpz_t  p1, p2, lambda;
    mpz_t *P1p, *P2p;  /* P1p[k] = p1^k and P2p[k] = p2^k for k = 0..h+1 */
    mpz_t  lambda_inv; /* lambda^{-1} mod N^h */
    mpz_t *inv_fact;   /* inv_fact[k] = (k!)^{-1} mod N^{h+1}, k = 0..h */
    SigKey sig;        /* Ed25519 issuance key */
} DJParams;

/* Setup: distinct primes p1, p2 of nbits/2 bits each, N = p1 p2 of nbits
   bits, r = 2^(h-1)+1, n = floor(log_r N). */
void dj_setup(DJParams *pp, int h, int nbits);
void dj_free(DJParams *pp);

/* Uniform element of Z*_N. */
void dj_rand_unit(const DJParams *pp, mpz_t x);

/* out = x^{N^s} mod N^{s+1}.  The value depends only on x mod N, and it is
   computed as h successive N-th powers on growing moduli. */
void dj_pow_Ns(const DJParams *pp, mpz_t out, const mpz_t x, int s);

/* The same value as dj_pow_Ns, computed with the factorization of N
   (authority only).  Modulo p^{s+1}, for a prime factor p of N, x^{N^s} is
   the root of X^{p-1} = 1 that is congruent to x^{N^s} modulo p: write
   x = omega u with omega^{p-1} = 1 and u = 1 mod p; then u has order
   dividing p^s, which divides N^s.  The root is found by Newton iteration
   from x^{N^s mod (p-1)} mod p, and the two roots are combined by the
   Chinese remainder theorem. */
void dj_pow_Ns_sk(const DJParams *pp, mpz_t out, const mpz_t x, int s);

/* out = (1+N)^m mod N^{s+1}, by the binomial expansion (s+1 terms). */
void dj_one_plus_N_pow(const DJParams *pp, mpz_t out, const mpz_t m, int s);

/* c = (1+N)^m rho^{N^s} mod N^{s+1} with rho uniform in Z*_N; rho is
   returned when rho_out is not NULL.  Sampling rho in Z*_N gives the same
   distribution of rho^{N^s} as sampling it in Z*_{N^{s+1}}. */
void dj_enc(const DJParams *pp, mpz_t c, const mpz_t m, int s, mpz_t rho_out);

/* dj_enc with the rerandomizer computed by dj_pow_Ns_sk (authority only);
   the ciphertext has the same distribution. */
void dj_enc_sk(const DJParams *pp, mpz_t c, const mpz_t m, int s, mpz_t rho_out);

/* v = Dec(c) in Z_{N^s} (authority only). */
void dj_dec(const DJParams *pp, mpz_t v, const mpz_t c, int s);

/* Fixed-width big-endian encoding of an element of Z_{N^{s+1}}. */
size_t dj_elem_len(const DJParams *pp, int s);
void   dj_write(const DJParams *pp, uint8_t *buf, const mpz_t x, int s);

/* Tag(sk, id_i, b, h'): encrypts N^{h-h'} r^{i-1} b and signs the ciphertext.
   i is 1-based.  sig must hold pp->sig.sig_bytes bytes; returns its length. */
size_t dj_tag(const DJParams *pp, mpz_t w, uint8_t *sig, int i, int b, int hprime);
int    dj_verify_tag(const DJParams *pp, const mpz_t w, const uint8_t *sig, size_t siglen);

/* Degrade(w) -> (w', xi): w' = w^N rho^{N^h}, xi = rho. */
void dj_degrade(const DJParams *pp, mpz_t w_out, mpz_t xi, const mpz_t w);

/* Merge(w_1..w_m) -> (w', xi): w' = rho^{N^h} prod_j w_j, xi = rho. */
void dj_merge(const DJParams *pp, mpz_t w_out, mpz_t xi, mpz_t *ws, int m);

/* The rerandomizer of Degrade and Merge, which does not depend on the tag
   and can be computed in advance: t = rho^{N^h} mod N^{h+1}, xi = rho. */
void dj_rerandomizer(const DJParams *pp, mpz_t t, mpz_t xi);

/* Degrade and Merge with a precomputed rerandomizer t: w' = w^N t and
   w' = t prod_j w_j.  dj_degrade and dj_merge are dj_rerandomizer followed
   by these. */
void dj_degrade_pre(const DJParams *pp, mpz_t w_out, const mpz_t t, const mpz_t w);
void dj_merge_pre(const DJParams *pp, mpz_t w_out, const mpz_t t, mpz_t *ws, int m);

/* Trace(sk, w): returns -1 for bottom, otherwise |S|, and sets in_S[i-1]
   for every id_i in S (in_S has pp->n entries). */
int dj_trace(const DJParams *pp, const mpz_t w, uint8_t *in_S);
