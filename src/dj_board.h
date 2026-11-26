#pragma once
#include "dj.h"

/* dj_board.h -> Board and VerifyBoard for Damgard-Jurik (Section "Board
   proofs with a public count").  The entry e_j = (1+N)^{b_j m(a_j)}
   rho_j^{N^h} mod N^{h+1} is accompanied by the bit ciphertext
   C_j = (1+N)^{b_j} tau_j^N mod N^2 and by a proof of the formula

     (e_j encrypts 0 and C_j encrypts 0)
       or (C_j encrypts 1 and (e_j encrypts m_1 or ... or e_j encrypts m_n)),

   where m_i = r^{i-1} is the plaintext of identifier i.  The statement does
   not contain the identifier a_j.  Each atom is a proof that a unit is an
   N^s-th power modulo N^{s+1}, and the formula is proved by the composition
   of Cramer, Damgard and Schoenmakers, with challenges in [0, 2^128) that
   add modulo 2^128, and the first 128 bits of SHA-256 as the random oracle.
   The differences of challenges are below the prime factors of N, hence
   units, which special soundness needs.  One proof that
   prod_j C_j (1+N)^{-c} is an N-th power shows the count c.  The authority
   signs the Merkle root of the entries once, with Ed25519.

   A proof is sent as its challenges and responses: the verifier recomputes
   each first message a = z^{N^s} u^{-e} and checks the hash.  Since
   z^{N^s} mod N^{s+1} depends only on z mod N, responses are sent mod N. */

#define DJ_CHAL_BITS 128

/* Number of responses of an entry proof: e enc 0, C enc 0, C enc 1, and
   e enc m_i for i = 1..n. */
#define DJ_NZ(n) ((n) + 3)

typedef struct {
    mpz_t  e;          /* the entry, Z*_{N^{h+1}} */
    mpz_t  C;          /* the bit ciphertext, Z*_{N^2} */
    mpz_t  c0;         /* challenge of the left side of the outer "or" */
    mpz_t *c;          /* c[i]: challenge of "e enc m_{i+1}", i = 0..n-1 */
    mpz_t *z;          /* DJ_NZ(n) responses mod N, in the order above */
} DJEntry;

typedef struct {
    size_t   M;
    DJEntry *E;
    uint8_t  root[32];
    uint8_t *sig;
    size_t   siglen;
    unsigned long c;   /* public count */
    mpz_t    ce, cz;   /* count proof: challenge and response mod N */
} DJBoard;

void dj_entry_init(const DJParams *pp, DJEntry *E);
void dj_entry_free(const DJParams *pp, DJEntry *E);
void dj_board_init(DJBoard *B, const DJParams *pp, size_t M);
void dj_board_free(DJBoard *B, const DJParams *pp);

/* The entry at position j for identifier slot in [1, n] and flag b; slot is
   ignored when b = 0.  tau_out, if not NULL, receives tau_j. */
void dj_entry_prove(const DJParams *pp, DJEntry *E, uint64_t j, int slot, int b,
                    mpz_t tau_out);
int  dj_entry_verify(const DJParams *pp, const DJEntry *E, uint64_t j);
/* Leaf of the Merkle tree: SHA-256(0x00 || C || SHA-256(0x02 || e)), so that
   the root fixes the bit ciphertext of each entry, and the bit ciphertexts
   can be checked against the root from the digests of the entries. */
void dj_entry_leaf(const DJParams *pp, const DJEntry *E, uint8_t leaf[32]);

/* Count proof for u = prod_j C_j mod N^2, with tau the product of the tau_j
   mod N. */
void dj_count_prove(const DJParams *pp, const mpz_t u, uint64_t M, unsigned long c,
                    const mpz_t tau, mpz_t ce, mpz_t cz);
int  dj_count_verify(const DJParams *pp, const mpz_t u, uint64_t M, unsigned long c,
                     const mpz_t ce, const mpz_t cz);

/* Board and VerifyBoard on the given number of threads. */
void dj_board(const DJParams *pp, DJBoard *B, const int *slot, const int *b,
              int threads);
int  dj_verify_board(const DJParams *pp, const DJBoard *B, int threads);

/* Serialized sizes in bytes. */
size_t dj_entry_bytes(const DJParams *pp);   /* e_j, C_j and its proof */
size_t dj_count_bytes(const DJParams *pp);   /* c and its proof */
size_t dj_root_bytes(const DJParams *pp);    /* root and its signature */
