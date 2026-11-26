#pragma once
#include <stddef.h>
#include <stdint.h>
#include "sig.h"

/* rlwe.h -> Ring-LWE instantiation of ECT (Section "Ring-LWE instantiation"):
   BGV ciphertexts over R = R_0[Y]/(Y^h), R_0 = Z[X]/(X^{N_0}+1), with the
   ciphertext modulus Q a product of at most two NTT primes q_i = 1 mod 2N_0
   and the plaintext modulus t = 1 mod 2N_0, so that each block of the
   plaintext splits into N_0 slots.

   An element of R_Q is stored as h blocks, each in RNS form (one residue
   polynomial per prime): index ((k * L) + i) * N_0 + c for block k, prime i
   and coefficient c. */

#define RLWE_MAX_PRIMES 2

__extension__ typedef unsigned __int128 rlwe_u128;
__extension__ typedef __int128 rlwe_i128;

typedef struct {
    uint64_t  q;        /* prime modulus, q < 2^62 */
    uint64_t  qinv;     /* -q^{-1} mod 2^64 (Montgomery) */
    uint64_t  r2;       /* 2^128 mod q */
    uint64_t  ninv;     /* N_0^{-1} mod q */
    uint64_t *psi;      /* psi^{brv(k)}, k < N_0 */
    uint64_t *psip;     /* Shoup companions of psi */
    uint64_t *ipsi;     /* psi^{-brv(k)} */
    uint64_t *ipsip;
} NTTTable;

typedef struct {
    int       N0, logN0, h, L, n;
    uint64_t  t;
    uint64_t  q[RLWE_MAX_PRIMES];
    int       logQ;             /* ceil(log2 Q) */
    rlwe_u128 Q;
    uint64_t  crt_inv;          /* q_0^{-1} mod q_1 when L = 2 */
    NTTTable  ntt[RLWE_MAX_PRIMES];
    NTTTable  ntt_t;            /* NTT modulo t, for slot encoding */
    double    sigma;            /* error standard deviation */
    int       Be;               /* error tail cut, ceil(6 sigma) */
    uint64_t  cdt[64];          /* cumulative table of |e| */
    /* keys, in NTT and Montgomery form */
    uint64_t *pk0, *pk1, *sk;
    SigKey    sig;              /* ML-DSA-44 issuance key */
} RLWEParams;

typedef struct {
    uint64_t *c0, *c1;          /* coefficient form, h * L * N0 residues each */
} RLWECt;

/* Setup: N0 a power of two, h the hop budget, logQ the target bit length of
   Q (split over L primes), t the plaintext prime. */
void rlwe_setup(RLWEParams *pp, int N0, int h, int logQ, int L, uint64_t t);
void rlwe_free(RLWEParams *pp);

void rlwe_ct_init(const RLWEParams *pp, RLWECt *c);
void rlwe_ct_free(RLWECt *c);
void rlwe_ct_copy(const RLWEParams *pp, RLWECt *dst, const RLWECt *src);

/* Size of a tag by the formula 2 N_0 h ceil(log2 Q) bits, in bytes. */
size_t rlwe_tag_bytes(const RLWEParams *pp);

/* Plaintexts are h * N0 coefficients modulo t (block-major). */
void rlwe_enc(const RLWEParams *pp, RLWECt *c, const uint64_t *m);
/* Returns the noise margin log2(Q/2) - log2 max |m + t e| in bits. */
double rlwe_dec(const RLWEParams *pp, uint64_t *m, const RLWECt *c);

/* Slot encoding: the plaintext block that is 1 in slot i (0-based) and 0
   elsewhere, and the slot values of a block. */
void rlwe_slot_unit(const RLWEParams *pp, uint64_t *block, int i);
void rlwe_slot_decode(const RLWEParams *pp, uint64_t *slots, const uint64_t *block);

/* Tag(sk, id_i, b, h'), i 1-based; sig holds pp->sig.sig_bytes bytes. */
size_t rlwe_tag(const RLWEParams *pp, RLWECt *w, uint8_t *sig, int i, int b, int hprime);
int    rlwe_verify_tag(const RLWEParams *pp, const RLWECt *w, const uint8_t *sig, size_t siglen);

/* Degrade(w): (Y c0, Y c1) + Enc(0).  Merge(w_1..w_m): sum + Enc(0). */
void rlwe_degrade(const RLWEParams *pp, RLWECt *out, const RLWECt *w);
void rlwe_merge(const RLWEParams *pp, RLWECt *out, RLWECt *const *ws, int m);

/* Trace: -1 for bottom, otherwise |S|, with in_S[i-1] set for id_i in S. */
int rlwe_trace(const RLWEParams *pp, const RLWECt *w, uint8_t *in_S);
