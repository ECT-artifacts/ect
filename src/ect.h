#pragma once
#include <relic/relic.h>
#include <stdint.h>
#include "params.h"
#include "nizk.h"
#include "schnorr.h"

/* ect.h -> Ephemeral Coin Tracing scheme  */

/* Maximum number of Degrade steps per tag (= hop budget H) */
#define MAX_PROOFS H

/* ECTProofStep -> per-hop ciphertext snapshot for NIZK chain replay.  */
typedef struct {
    bn_t c1_in,  c2_in;    /* ciphertext BEFORE this hop (input to Degrade) */
    bn_t c1_out, c2_out;   /* ciphertext AFTER this hop  (output of Degrade) */
} ECTProofStep;

/* ECTTag -> the watermark attached to a coin.  */
typedef struct {
    bn_t c1, c2;        /* current ciphertext -> changes with each Degrade */
    bn_t c1_0, c2_0;    /* original ciphertext at issuance -> never changes */
    int  hops;
    NIZKProof    proofs[MAX_PROOFS];
    ECTProofStep steps [MAX_PROOFS]; /* per-hop (c_before, c_after) snapshots */
    int          n_proofs;
    SchnorrSig sig;
    int        has_sig;
} ECTTag;

/*
  ECTParams -> public parameters and authority key material.
 
    p, g        : group modulus and generator (from params.h)
    q, qh       : base prime and q^H (exponent space modulus)
    pk, sk      : ElGamal encryption key pair (authority only: sk)
    pk_sig, sk_sig : Schnorr signing key pair (authority only: sk_sig)
    h, r, n     : scheme parameters H, R_BASE, N
 
    trace_table : precomputed lookup table for Trace decryption.
                  trace_table[i][t] = g^(q^t * r^i) mod p
                  for i=0..n-1 (user i+1), t=0..h-1 (depth t).
                  Computed once at ect_setup.
 
    tagged_users  : list of user IDs (1-indexed) the authority has tagged
                    with b=1 (traceable). Maintained by ect_tag calls.
    n_tagged      : number of entries in tagged_users.
 */
typedef struct {
    bn_t p, g, q, qh;
    bn_t pk, sk;
    bn_t pk_sig, sk_sig;
    int  h, r, n;
    bn_t trace_table[N][H];  /* [user_idx 0-based][depth] = g^(q^t * r^i) */
    int  tagged_users[N];    /* user IDs (1-based) the authority tagged b=1 */
    int  n_tagged;           /* number of tagged users currently tracked */
} ECTParams;

 
/* ect_mult_max -> largest per-depth multiplicity Trace will decode.
   Defaults to MULT_MAX; raising it widens the meet-in-the-middle table by a
   factor of (ect_mult_max+1)^(K/2), lowering it narrows the recoverable set.
   Clamped to [1, R_BASE-1] on use. */
extern int ect_mult_max;

void ect_tag_new (ECTTag *w);
void ect_tag_free(ECTTag *w);
void ect_tag_copy(ECTTag *dst, const ECTTag *src);  /* deep copy */

void ect_params_new (ECTParams *pp);
void ect_params_free(ECTParams *pp);

/* ect_clear_tagged -> reset the authority's tagged users list. */
void ect_clear_tagged(ECTParams *pp);

/* ect_setup -> initialise ECTParams from the constants in params.h. */
void ect_setup(ECTParams *pp);

 
/* Scheme algorithms*/
 

/* ect_tag -> Tag(sktrace, id, b) ->  */
void ect_tag(ECTTag *w, ECTParams *pp, int id, int b);

/* ect_degrade -> Degrade(w) -> (w', proof) */
void ect_degrade(ECTTag *w, const ECTParams *pp,
                 bn_t *C1s, bn_t *C2s, int m, int real_idx);

/* ect_merge -> Merge(w_1, ..., w_k) -> w_out.
   Component-wise ciphertext multiplication, then rerandomization with a
   fresh encryption of zero.  The host ledger, not Merge, authenticates
   issued inputs and validates transitions. */
void ect_merge(ECTTag *w_out, ECTTag *tags[], int k, const ECTParams *pp);

/*  * ect_trace -> Trace(sktrace, w) -> S ∪ {bot} */
int ect_trace(const ECTParams *pp, const ECTTag *w);

/* Silent version for benchmarking */
int ect_trace_silent(const ECTParams *pp, const ECTTag *w);

/* ect_trace_ciphertext -> run Trace on a bare (c1, c2) ciphertext.*/
int ect_trace_ciphertext(const ECTParams *pp,
                          const bn_t c1, const bn_t c2,
                          int verbose);
 
 
/* Verification     */
 

/* ect_verify_proofs -> verify the full NIZK proof chain on a tag. */
int ect_verify_proofs(const ECTTag *w, const ECTParams *pp);

/* ect_verify_tag_sig -> verify the authority's Schnorr signature on a tag */
int ect_verify_tag_sig(const ECTTag *w, const ECTParams *pp);
