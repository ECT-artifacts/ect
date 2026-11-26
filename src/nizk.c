#include "nizk.h"
#include "block.h"
#include <string.h>
#include <stdlib.h>

/* nizk.c -> CDS one-out-of-many proof for R_exp. */

void nizk_proof_new(NIZKProof *pi)
{
    pi->m = 0;
    for (int j = 0; j < MAX_INPUTS; j++) {
        bn_null(pi->A1[j]); bn_new(pi->A1[j]);
        bn_null(pi->A2[j]); bn_new(pi->A2[j]);
        bn_null(pi->z[j]);  bn_new(pi->z[j]);
        bn_null(pi->e[j]);  bn_new(pi->e[j]);
    }
}

void nizk_proof_free(NIZKProof *pi)
{
    for (int j = 0; j < MAX_INPUTS; j++) {
        bn_free(pi->A1[j]); bn_free(pi->A2[j]);
        bn_free(pi->z[j]);  bn_free(pi->e[j]);
    }
    (void)pi;
}

 
/* Fiat-Shamir global challenge  */
 

/*
  fiat_shamir_cds -> compute the global Fiat-Shamir challenge.
 
  e = SHA-256( pk || c1_0||c2_0 || ... || c1_{m-1}||c2_{m-1}
                 || c1'||c2'
                 || A1_0||A2_0 || ... || A1_{m-1}||A2_{m-1} )  mod q^H
 */
static void fiat_shamir_cds(bn_t e,
                             const bn_t pk,
                             bn_t *C1s, bn_t *C2s, int m,
                             const bn_t c1p, const bn_t c2p,
                             bn_t *A1s, bn_t *A2s,
                             const bn_t qh)
{
    int     n_elems = 1 + 2*m + 2 + 2*m;   /* pk, C, c', commitments */
    size_t  total   = (size_t)n_elems * WM_ELEM_BYTES;
    uint8_t *buf    = (uint8_t *)malloc(total);
    uint8_t  digest[RLC_MD_LEN];

    memset(buf, 0, total);

    #define WRITE(idx, x) do {                                          \
        int _n = bn_size_bin(x);                                        \
        if (_n > WM_ELEM_BYTES) _n = WM_ELEM_BYTES;                    \
        if (_n > 0)                                                     \
            bn_write_bin(buf + (idx)*WM_ELEM_BYTES + WM_ELEM_BYTES-_n, \
                         _n, (x));                                      \
    } while(0)

    int slot = 0;
    WRITE(slot++, pk);
    for (int j = 0; j < m; j++) { WRITE(slot++, C1s[j]); WRITE(slot++, C2s[j]); }
    WRITE(slot++, c1p);
    WRITE(slot++, c2p);
    for (int j = 0; j < m; j++) { WRITE(slot++, A1s[j]); WRITE(slot++, A2s[j]); }

    #undef WRITE

    md_map(digest, buf, total);
    free(buf);

    bn_read_bin(e, digest, RLC_MD_LEN);
    bn_mod(e, e, qh);
}

 
/* Prove */
 

void nizk_prove(NIZKProof *pi,
                const bn_t pk,
                bn_t *C1s, bn_t *C2s, int m, int real_idx,
                const bn_t c1p, const bn_t c2p,
                const bn_t r_prime,
                const bn_t g, const bn_t q, const bn_t qh, const bn_t p)
{
    pi->m = m;

    bn_t k, global_e, sum_sim_e, tmp;
    bn_null(k); bn_null(global_e); bn_null(sum_sim_e); bn_null(tmp);

    RLC_TRY {
        bn_new(k); bn_new(global_e); bn_new(sum_sim_e); bn_new(tmp);

        /* Simulate all non-real clauses.
          For j != real_idx: pick e_j, z_j at random and compute
          commitments backwards from the verification equation:
            g^z_j * c1'^{e_j} == A1_j * (c1_j^q)^{e_j}
          =>  A1_j = g^z_j * c1'^{e_j} * (c1_j^q)^{-e_j}  mod p
          and similarly for A2_j. */
        bn_zero(sum_sim_e);
        for (int j = 0; j < m; j++) {
            if (j == real_idx) continue;   /* handle real index later */

            bn_rand_mod(pi->e[j], qh);
            bn_rand_mod(pi->z[j], qh);

            bn_t c1q, c2q, inv1, inv2;
            bn_null(c1q); bn_null(c2q);
            bn_null(inv1); bn_null(inv2);
            bn_new(c1q);  bn_new(c2q);
            bn_new(inv1); bn_new(inv2);
            bn_mxp(c1q, C1s[j], q, p);
            bn_mxp(c2q, C2s[j], q, p);

            /* A1_j = g^z_j * c1'^{e_j} * (c1_j^q)^{-e_j} mod p */
            bn_mxp(pi->A1[j], g,   pi->z[j], p);
            bn_mxp(tmp,       c1p, pi->e[j], p);
            bn_mul(pi->A1[j], pi->A1[j], tmp);
            bn_mod(pi->A1[j], pi->A1[j], p);
            bn_mxp(inv1, c1q, pi->e[j], p);
            bn_mod_inv(inv1, inv1, p);
            bn_mul(pi->A1[j], pi->A1[j], inv1);
            bn_mod(pi->A1[j], pi->A1[j], p);

            /* A2_j = pk^z_j * c2'^{e_j} * (c2_j^q)^{-e_j} mod p */
            bn_mxp(pi->A2[j], pk,  pi->z[j], p);
            bn_mxp(tmp,       c2p, pi->e[j], p);
            bn_mul(pi->A2[j], pi->A2[j], tmp);
            bn_mod(pi->A2[j], pi->A2[j], p);
            bn_mxp(inv2, c2q, pi->e[j], p);
            bn_mod_inv(inv2, inv2, p);
            bn_mul(pi->A2[j], pi->A2[j], inv2);
            bn_mod(pi->A2[j], pi->A2[j], p);

            bn_add(sum_sim_e, sum_sim_e, pi->e[j]);
            bn_mod(sum_sim_e, sum_sim_e, qh);

            bn_free(c1q); bn_free(c2q);
            bn_free(inv1); bn_free(inv2);
        }

        /* Commit for real clause.
         * A1_i = g^k,  A2_i = pk^k */
        bn_rand_mod(k, qh);
        bn_mxp(pi->A1[real_idx], g,  k, p);
        bn_mxp(pi->A2[real_idx], pk, k, p);

        /* Compute global Fiat-Shamir challenge */
        fiat_shamir_cds(global_e, pk,
                        C1s, C2s, m, c1p, c2p,
                        pi->A1, pi->A2, qh);

        /* Real clause challenge = global_e - sum of simulated */
        bn_sub(pi->e[real_idx], global_e, sum_sim_e);
        if (bn_sign(pi->e[real_idx]) == RLC_NEG)
            bn_add(pi->e[real_idx], pi->e[real_idx], qh);
        bn_mod(pi->e[real_idx], pi->e[real_idx], qh);

        /* Real clause response: z_i = k - e_i * r' mod q^H */
        bn_mul(tmp, pi->e[real_idx], r_prime);
        bn_mod(tmp, tmp, qh);
        bn_sub(pi->z[real_idx], k, tmp);
        if (bn_sign(pi->z[real_idx]) == RLC_NEG)
            bn_add(pi->z[real_idx], pi->z[real_idx], qh);
        bn_mod(pi->z[real_idx], pi->z[real_idx], qh);

    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(k); bn_free(global_e); bn_free(sum_sim_e); bn_free(tmp);
    }
}

 
/* Verify */
 

int nizk_verify(const NIZKProof *pi,
                const bn_t pk,
                bn_t *C1s, bn_t *C2s, int m,
                const bn_t c1p, const bn_t c2p,
                const bn_t g, const bn_t q, const bn_t qh, const bn_t p)
{
    volatile int ok = 0;
    bn_t global_e, sum_e, lhs, rhs, tmp, c1q, c2q;
    bn_null(global_e); bn_null(sum_e);
    bn_null(lhs); bn_null(rhs); bn_null(tmp);
    bn_null(c1q); bn_null(c2q);

    RLC_TRY {
        bn_new(global_e); bn_new(sum_e);
        bn_new(lhs); bn_new(rhs); bn_new(tmp);
        bn_new(c1q); bn_new(c2q);

        /* Recompute global Fiat-Shamir challenge */
        fiat_shamir_cds(global_e, pk,
                        C1s, C2s, m, c1p, c2p,
                        (bn_t *)pi->A1, (bn_t *)pi->A2, qh);

        /* Check sum constraint: sum e[j] == global_e mod q^H */
        bn_zero(sum_e);
        for (int j = 0; j < m; j++) {
            bn_add(sum_e, sum_e, pi->e[j]);
            bn_mod(sum_e, sum_e, qh);
        }
        if (bn_cmp(sum_e, global_e) != RLC_EQ) goto done;

        /* Check each clause j */
        for (int j = 0; j < m; j++) {
            bn_mxp(c1q, C1s[j], q, p);
            bn_mxp(c2q, C2s[j], q, p);

            /* Check 1: g^z[j] * c1'^{e[j]} == A1[j] * (c1_j^q)^{e[j]} */
            bn_mxp(lhs, g,   pi->z[j], p);
            bn_mxp(tmp, c1p, pi->e[j], p);
            bn_mul(lhs, lhs, tmp);
            bn_mod(lhs, lhs, p);

            bn_mxp(rhs, c1q, pi->e[j], p);
            bn_mul(rhs, rhs, pi->A1[j]);
            bn_mod(rhs, rhs, p);

            if (bn_cmp(lhs, rhs) != RLC_EQ) goto done;
 
            /* Check 2: pk^z[j] * c2'^{e[j]} == A2[j] * (c2_j^q)^{e[j]} */
            bn_mxp(lhs, pk,  pi->z[j], p);
            bn_mxp(tmp, c2p, pi->e[j], p);
            bn_mul(lhs, lhs, tmp);
            bn_mod(lhs, lhs, p);

            bn_mxp(rhs, c2q, pi->e[j], p);
            bn_mul(rhs, rhs, pi->A2[j]);
            bn_mod(rhs, rhs, p);

            if (bn_cmp(lhs, rhs) != RLC_EQ) goto done;
        }

        ok = 1;
done:;
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(global_e); bn_free(sum_e);
        bn_free(lhs); bn_free(rhs); bn_free(tmp);
        bn_free(c1q); bn_free(c2q);
    }
    return (int)ok;
}

 
/* Benchmark variants: dynamic m*/
 

void nizk_bench_prove(const bn_t pk,
                      bn_t *C1s, bn_t *C2s, int m,
                      const bn_t c1p, const bn_t c2p,
                      const bn_t r_prime,
                      const bn_t g, const bn_t q, const bn_t qh, const bn_t p,
                      bn_t **A1_out, bn_t **A2_out,
                      bn_t **z_out,  bn_t **e_out)
{
    bn_t *A1 = (bn_t*)malloc((size_t)m * sizeof(bn_t));
    bn_t *A2 = (bn_t*)malloc((size_t)m * sizeof(bn_t));
    bn_t *z  = (bn_t*)malloc((size_t)m * sizeof(bn_t));
    bn_t *e  = (bn_t*)malloc((size_t)m * sizeof(bn_t));
    for (int j = 0; j < m; j++) {
        bn_null(A1[j]); bn_new(A1[j]);
        bn_null(A2[j]); bn_new(A2[j]);
        bn_null(z[j]);  bn_new(z[j]);
        bn_null(e[j]);  bn_new(e[j]);
    }

    bn_t k, global_e, sum_sim_e, tmp;
    bn_null(k); bn_null(global_e); bn_null(sum_sim_e); bn_null(tmp);

    RLC_TRY {
        bn_new(k); bn_new(global_e); bn_new(sum_sim_e); bn_new(tmp);
        bn_zero(sum_sim_e);

        /* Simulate all non-real clauses (j != 0). */
        for (int j = 1; j < m; j++) {
            bn_rand_mod(e[j], qh);
            bn_rand_mod(z[j], qh);
            bn_t c1q, c2q, inv1, inv2;
            bn_null(c1q); bn_null(c2q);
            bn_null(inv1); bn_null(inv2);
            bn_new(c1q); bn_new(c2q);
            bn_new(inv1); bn_new(inv2);
            bn_mxp(c1q, C1s[j], q, p);
            bn_mxp(c2q, C2s[j], q, p);

            /* A1[j] = g^z[j] * c1'^{e[j]} * (c1_j^q)^{-e[j]} */
            bn_mxp(A1[j], g,   z[j], p);
            bn_mxp(tmp,   c1p, e[j], p);
            bn_mul(A1[j], A1[j], tmp); bn_mod(A1[j], A1[j], p);
            bn_mxp(inv1, c1q, e[j], p); bn_mod_inv(inv1, inv1, p);
            bn_mul(A1[j], A1[j], inv1); bn_mod(A1[j], A1[j], p);

            /* A2[j] = pk^z[j] * c2'^{e[j]} * (c2_j^q)^{-e[j]} */
            bn_mxp(A2[j], pk,  z[j], p);
            bn_mxp(tmp,   c2p, e[j], p);
            bn_mul(A2[j], A2[j], tmp); bn_mod(A2[j], A2[j], p);
            bn_mxp(inv2, c2q, e[j], p); bn_mod_inv(inv2, inv2, p);
            bn_mul(A2[j], A2[j], inv2); bn_mod(A2[j], A2[j], p);

            bn_add(sum_sim_e, sum_sim_e, e[j]);
            bn_mod(sum_sim_e, sum_sim_e, qh);
            bn_free(c1q); bn_free(c2q);
            bn_free(inv1); bn_free(inv2);
        }

        /* Real clause (j=0): honest commitment */
        bn_rand_mod(k, qh);
        bn_mxp(A1[0], g,  k, p);
        bn_mxp(A2[0], pk, k, p);

        /* Global Fiat-Shamir challenge */
        fiat_shamir_cds(global_e, pk, C1s, C2s, m,
                        c1p, c2p, A1, A2, qh);

        /* Real clause: e[0] = global_e - sum(e[j], j>0) */
        bn_sub(e[0], global_e, sum_sim_e);
        if (bn_sign(e[0]) == RLC_NEG) bn_add(e[0], e[0], qh);
        bn_mod(e[0], e[0], qh);

        /* Real clause response: z[0] = k - e[0]*r' */
        bn_mul(tmp, e[0], r_prime);
        bn_mod(tmp, tmp, qh);
        bn_sub(z[0], k, tmp);
        if (bn_sign(z[0]) == RLC_NEG) bn_add(z[0], z[0], qh);
        bn_mod(z[0], z[0], qh);

    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(k); bn_free(global_e); bn_free(sum_sim_e); bn_free(tmp);
    }

    *A1_out = A1; *A2_out = A2; *z_out = z; *e_out = e;
}

int nizk_bench_verify(const bn_t pk,
                      bn_t *C1s, bn_t *C2s, int m,
                      const bn_t c1p, const bn_t c2p,
                      bn_t *A1, bn_t *A2, bn_t *z, bn_t *e,
                      const bn_t g, const bn_t q, const bn_t qh, const bn_t p)
{
    volatile int ok = 0;
    bn_t global_e, sum_e, lhs, rhs, tmp, c1q, c2q;
    bn_null(global_e); bn_null(sum_e);
    bn_null(lhs); bn_null(rhs); bn_null(tmp);
    bn_null(c1q); bn_null(c2q);

    RLC_TRY {
        bn_new(global_e); bn_new(sum_e);
        bn_new(lhs); bn_new(rhs); bn_new(tmp);
        bn_new(c1q); bn_new(c2q);

        fiat_shamir_cds(global_e, pk, C1s, C2s, m,
                        c1p, c2p, A1, A2, qh);

        /* Sum constraint */
        bn_zero(sum_e);
        for (int j = 0; j < m; j++) {
            bn_add(sum_e, sum_e, e[j]);
            bn_mod(sum_e, sum_e, qh);
        }
        if (bn_cmp(sum_e, global_e) != RLC_EQ) goto done;

        /* Check each clause */
        for (int j = 0; j < m; j++) {
            bn_mxp(c1q, C1s[j], q, p);
            bn_mxp(c2q, C2s[j], q, p);

            /* g^z[j] * c1'^{e[j]} == A1[j] * (c1_j^q)^{e[j]} */
            bn_mxp(lhs, g,   z[j], p);
            bn_mxp(tmp, c1p, e[j], p);
            bn_mul(lhs, lhs, tmp); bn_mod(lhs, lhs, p);
            bn_mxp(rhs, c1q, e[j], p);
            bn_mul(rhs, rhs, A1[j]); bn_mod(rhs, rhs, p);
            if (bn_cmp(lhs, rhs) != RLC_EQ) goto done;

            /* pk^z[j] * c2'^{e[j]} == A2[j] * (c2_j^q)^{e[j]} */
            bn_mxp(lhs, pk,  z[j], p);
            bn_mxp(tmp, c2p, e[j], p);
            bn_mul(lhs, lhs, tmp); bn_mod(lhs, lhs, p);
            bn_mxp(rhs, c2q, e[j], p);
            bn_mul(rhs, rhs, A2[j]); bn_mod(rhs, rhs, p);
            if (bn_cmp(lhs, rhs) != RLC_EQ) goto done;
        }
        ok = 1;
done:;
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(global_e); bn_free(sum_e);
        bn_free(lhs); bn_free(rhs); bn_free(tmp);
        bn_free(c1q); bn_free(c2q);
    }
    return (int)ok;
}
