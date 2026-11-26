#pragma once
#include <relic/relic.h>

/* nizk.h -> CDS one-out-of-many NIZK proof for the Degrade relation. */

#define MAX_INPUTS 1     /* input ciphertexts in ECTTag proof chain -> always 1 */
#define MAX_INPUTS_BENCH 128  /* upper bound for standalone NIZK benchmark only */

/* NIZKProof -> CDS proof for MAX_INPUTS clauses.
    m       : actual number of clauses (1 <= m <= MAX_INPUTS)
    A1[j]   : commitment g^k_j mod p           (Z*_p element)
    A2[j]   : commitment pk^k_j mod p          (Z*_p element)
    z[j]    : response z_j mod q^H             (Z_{q^H} scalar)
    e[j]    : per-clause challenge e_j mod q^H (Z_{q^H} scalar)
              constraint: sum e[j] = global Fiat-Shamir challenge
 */
typedef struct {
    int  m;
    bn_t A1[MAX_INPUTS];
    bn_t A2[MAX_INPUTS];
    bn_t z [MAX_INPUTS];
    bn_t e [MAX_INPUTS];
} NIZKProof;

void nizk_proof_new (NIZKProof *pi);
void nizk_proof_free(NIZKProof *pi);

/*
  nizk_prove -> generate a CDS one-out-of-many proof. 
    pi         : output proof s
    pk         : ElGamal public key
    C1s, C2s   : arrays of m input ciphertext components (C1s[j]=c1_j)
    m          : number of input ciphertexts (1 <= m <= MAX_INPUTS)
    real_idx   : index of the ciphertext actually degraded (0-based)
    c1p, c2p   : output ciphertext c' after degradation + rerandomization
    r_prime    : rerandomization scalar used for the real clause
    g,q,qh,p   : group parameters
 */
void nizk_prove(NIZKProof *pi,
                const bn_t pk,
                bn_t *C1s, bn_t *C2s, int m, int real_idx,
                const bn_t c1p, const bn_t c2p,
                const bn_t r_prime,
                const bn_t g, const bn_t q, const bn_t qh, const bn_t p);

/* nizk_verify -> verify a CDS one-out-of-many proof. */
int nizk_verify(const NIZKProof *pi,
                const bn_t pk,
                bn_t *C1s, bn_t *C2s, int m,
                const bn_t c1p, const bn_t c2p,
                const bn_t g, const bn_t q, const bn_t qh, const bn_t p);

/*nizk_bench_prove / nizk_bench_verify -> benchmark variants for m > MAX_INPUTS. */
void nizk_bench_prove(const bn_t pk,
                      bn_t *C1s, bn_t *C2s, int m,
                      const bn_t c1p, const bn_t c2p,
                      const bn_t r_prime,
                      const bn_t g, const bn_t q, const bn_t qh, const bn_t p,
                      /* outputs -> caller must free */
                      bn_t **A1_out, bn_t **A2_out,
                      bn_t **z_out,  bn_t **e_out);

int nizk_bench_verify(const bn_t pk,
                      bn_t *C1s, bn_t *C2s, int m,
                      const bn_t c1p, const bn_t c2p,
                      bn_t *A1, bn_t *A2, bn_t *z, bn_t *e,
                      const bn_t g, const bn_t q, const bn_t qh, const bn_t p);

