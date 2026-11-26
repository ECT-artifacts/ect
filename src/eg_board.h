#pragma once
#include <gmp.h>
#include <stddef.h>
#include <stdint.h>
#include "sig.h"
#include "ect.h"

/* eg_board.h -> Board and VerifyBoard for exponential ElGamal (Section
   "Board proofs with a public count").  The entry e_j = Enc(b_j m(a_j);
   rho_j) is accompanied by the bit ciphertext C_j = Enc(b_j; tau_j) and by a
   proof of the formula

     (e_j encrypts 0 and C_j encrypts 0)
       or (C_j encrypts 1 and (e_j encrypts m_1 or ... or e_j encrypts m_n)),

   where m_i = r^{i-1} is the plaintext of identifier i.  The statement does
   not contain the identifier a_j.  Each atom is a proof of equal discrete
   logarithms for the bases (g, pk), and the formula is proved by the
   composition of Cramer, Damgard and Schoenmakers: the challenges of the
   two sides of an "or" sum to the challenge of the "or", and both sides of
   an "and" take its challenge.  Challenges lie in [0, 2^128) and add
   modulo 2^128; the Fiat-Shamir challenge is SHA-512 reduced modulo 2^128.  One proof that prod_j C_j (1, g^{-c}) encrypts zero shows the
   count c.  The authority signs the Merkle root of the entries once, with
   Ed25519.

   The proof carries the first messages of its atoms.  The verifier checks
   the hash on them and all atom equations of a batch of entries at once: it
   raises each equation to a random 128-bit coefficient, multiplies them, and
   raises the product to the cofactor (p-1)/q^h, so that only the component
   in the subgroup of order q^h must vanish; an atom equation that fails
   there passes with probability at most 2^-128, and one that fails only
   outside it does not affect special soundness, since its statement lies in
   the subgroup.  The four group elements of each entry are tested for
   membership in that subgroup by 128 random-subset products, each raised to
   q^h; an element outside the subgroup passes each with probability at most
   1/2, since the cofactor has the small factors 2, 3, 7 and 19.

   The arithmetic is on GMP, so that the entries can be proved and verified
   on several threads; fixed-base tables speed up the powers of g and pk,
   and per-entry tables the powers of e1^{-1} and e2^{-1}. */

typedef struct {
    int    w, nwin;   /* window width and number of windows */
    mpz_t *t;         /* t[i (2^w - 1) + d - 1] = base^{d 2^{w i}} mod p */
} EGFixedBase;

typedef struct {
    int    n, h;
    mpz_t  p, q, qh, g, pk;
    mpz_t  cmod;       /* 2^128: challenges lie in [0, cmod) */
    mpz_t  cof;        /* (p-1)/q^h */
    mpz_t *m;          /* m[i] = r^i, the plaintext of identifier i + 1 */
    EGFixedBase fg, fpk;
    uint8_t ctx[64];   /* SHA-512 of (p, g, pk), hashed into every challenge */
    SigKey sig;        /* Ed25519 key that signs the root of each board */
} EGBoardParams;

/* Number of responses of an entry proof: e enc 0, C enc 0, C enc 1, and
   e enc m_i for i = 1..n. */
#define EG_NZ(n) ((n) + 3)

typedef struct {
    mpz_t  e1, e2;     /* the entry */
    mpz_t  C1, C2;     /* the bit ciphertext */
    mpz_t  c0;         /* challenge of the left side of the outer "or" */
    mpz_t *c;          /* c[i]: challenge of "e enc m_{i+1}", i = 0..n-1 */
    mpz_t *z;          /* EG_NZ(n) responses mod q^h, in the order above */
    mpz_t (*T)[2];     /* EG_NZ(n) first messages (T1, T2) */
} EGEntry;

/* Number of random-subset products of the membership test. */
#define EG_MEM_REPS 128

/* Accumulated checks of a batch of entries: the check of the atom equations
   is (g^GE pk^PE R^{-1})^{(p-1)/q^h} = 1, that of membership mem[r]^{q^h} = 1
   for every r. */
typedef struct {
    mpz_t GE, PE, R;
    mpz_t mem[EG_MEM_REPS];
    int   bad;
} EGBatch;

typedef struct {
    size_t   M;
    EGEntry *E;
    uint8_t  root[32];
    uint8_t *sig;
    size_t   siglen;
    unsigned long c;   /* public count */
    mpz_t    ce, cz;   /* count proof */
} EGBoard;

/* Converts the public key of ect_setup and builds the fixed-base tables of
   window width w (8 takes about 64 MB). */
void eg_board_params_init(EGBoardParams *bp, const ECTParams *pp, int w);
void eg_board_params_free(EGBoardParams *bp);

void eg_entry_init(const EGBoardParams *bp, EGEntry *E);
void eg_entry_free(const EGBoardParams *bp, EGEntry *E);
void eg_board_init(EGBoard *B, const EGBoardParams *bp, size_t M);
void eg_board_free(EGBoard *B, const EGBoardParams *bp);

/* The entry at position j for identifier slot in [1, n] and flag b; slot is
   ignored when b = 0.  tau_out, if not NULL, receives tau_j. */
void eg_entry_prove(const EGBoardParams *bp, EGEntry *E, uint64_t j, int slot,
                    int b, mpz_t tau_out);
int  eg_entry_verify(const EGBoardParams *bp, const EGEntry *E, uint64_t j);

/* Batch verification: eg_batch_add checks the ranges and the hash of the
   entry at position j and adds its equations to the batch; eg_batch_merge
   adds the batch src to dst; eg_batch_check returns 1 if every entry added
   verifies.  eg_entry_verify is a batch of one entry. */
void eg_batch_init(EGBatch *Bt);
void eg_batch_free(EGBatch *Bt);
void eg_batch_add(const EGBoardParams *bp, EGBatch *Bt, const EGEntry *E, uint64_t j);
void eg_batch_merge(const EGBoardParams *bp, EGBatch *dst, const EGBatch *src);
int  eg_batch_check(const EGBoardParams *bp, const EGBatch *Bt);
/* Leaf of the Merkle tree: SHA-256(0x00 || C1 || C2 || SHA-256(0x02 || e1 || e2)),
   so that the root fixes the bit ciphertext of each entry, and the bit
   ciphertexts can be checked against the root from the digests of the
   entries. */
void eg_entry_leaf(const EGEntry *E, uint8_t leaf[32]);

/* Count proof for X = prod_j C1_j and Y = prod_j C2_j, with tau the sum of
   the tau_j mod q^h. */
void eg_count_prove(const EGBoardParams *bp, const mpz_t X, const mpz_t Y,
                    uint64_t M, unsigned long c, const mpz_t tau,
                    mpz_t ce, mpz_t cz);
int  eg_count_verify(const EGBoardParams *bp, const mpz_t X, const mpz_t Y,
                     uint64_t M, unsigned long c, const mpz_t ce, const mpz_t cz);

/* Board and VerifyBoard on the given number of threads. */
void eg_board(const EGBoardParams *bp, EGBoard *B, const int *slot, const int *b,
              int threads);
int  eg_verify_board(const EGBoardParams *bp, const EGBoard *B, int threads);

size_t eg_entry_bytes(const EGBoardParams *bp);   /* e_j, C_j and its proof, with its first messages */
size_t eg_count_bytes(const EGBoardParams *bp);   /* c and its proof */
size_t eg_root_bytes(const EGBoardParams *bp);    /* root and its signature */
