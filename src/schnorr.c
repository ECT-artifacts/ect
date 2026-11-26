#include "schnorr.h"
#include <stdlib.h>
#include <string.h>

/* schnorr.c -> Schnorr signature scheme (nizk fiat shamir) */

 
/* Internal: Fiat-Shamir challenge  */
 

/* SCH_ELEM_BYTES: fixed slot width for each group element in the hash.
    Hash input: pk || R || c1 || c2  (4 elements of P_BYTES each)
    Hash output: SHA-256 -> 32 bytes -> reduced mod q^H
 */
#define SCH_ELEM_BYTES P_BYTES

static void schnorr_challenge(bn_t e,
                               const bn_t pk, const bn_t R_commit,
                               const bn_t c1, const bn_t c2,
                               const bn_t qh)
{
    size_t   total = 4 * SCH_ELEM_BYTES;
    uint8_t *buf   = (uint8_t *)malloc(total);
    uint8_t  digest[RLC_MD_LEN];

    memset(buf, 0, total);

    /* Write each element right-aligned (big-endian, zero-padded left) */
    #define WRITE_SLOT(i, x) do {                                       \
        int _n = bn_size_bin(x);                                        \
        if (_n > SCH_ELEM_BYTES) _n = SCH_ELEM_BYTES;                  \
        if (_n > 0)                                                     \
            bn_write_bin(buf + (i)*SCH_ELEM_BYTES + SCH_ELEM_BYTES - _n,\
                         _n, (x));                                      \
    } while(0)

    WRITE_SLOT(0, pk);
    WRITE_SLOT(1, R_commit);
    WRITE_SLOT(2, c1);
    WRITE_SLOT(3, c2);

    #undef WRITE_SLOT

    md_map(digest, buf, total);
    free(buf);

    /* e = SHA-256(...) mod q^H */
    bn_read_bin(e, digest, RLC_MD_LEN);
    bn_mod(e, e, qh);
}


void schnorr_sig_new(SchnorrSig *sig)
{
    bn_null(sig->R); bn_null(sig->s);
    bn_new(sig->R);  bn_new(sig->s);
}

void schnorr_sig_free(SchnorrSig *sig)
{
    bn_free(sig->R); bn_free(sig->s);
    (void)sig;
}

 
/* Key generation   */
 

void schnorr_keygen(bn_t sk, bn_t pk,
                    const bn_t g, const bn_t qh, const bn_t p)
{
    /* sk <- Z_{q^H},  pk = g^sk mod p */
    bn_rand_mod(sk, qh);
    bn_mxp(pk, g, sk, p);
}

 
/* Sign  */
 

void schnorr_sign(SchnorrSig *sig,
                  const bn_t sk, const bn_t pk,
                  const bn_t c1, const bn_t c2,
                  const bn_t g, const bn_t qh, const bn_t p)
{
    bn_t k, e, esk;
    bn_null(k); bn_null(e); bn_null(esk);
    RLC_TRY {
        bn_new(k); bn_new(e); bn_new(esk);

        /* Step 1: commit -> k -< Z_{q^H}, R = g^k mod p */
        bn_rand_mod(k, qh);
        bn_mxp(sig->R, g, k, p);

        /* Step 2: Fiat-Shamir challenge
         * e = H(pk || R || c1 || c2) mod q^H */
        schnorr_challenge(e, pk, sig->R, c1, c2, qh);

        /* Step 3: response -> s = k - e*sk mod q^H */
        bn_mul(esk, e, sk);
        bn_mod(esk, esk, qh);
        bn_sub(sig->s, k, esk);
        if (bn_sign(sig->s) == RLC_NEG)
            bn_add(sig->s, sig->s, qh);
        bn_mod(sig->s, sig->s, qh);

    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(k); bn_free(e); bn_free(esk);
    }
}

 
/* Verify */
 

int schnorr_verify(const SchnorrSig *sig,
                   const bn_t pk,
                   const bn_t c1, const bn_t c2,
                   const bn_t g, const bn_t qh, const bn_t p)
{
    volatile int ok = 0;
    bn_t e, lhs, rhs;
    bn_null(e); bn_null(lhs); bn_null(rhs);
    RLC_TRY {
        bn_new(e); bn_new(lhs); bn_new(rhs);

        /* Recompute challenge from the statement */
        schnorr_challenge(e, pk, sig->R, c1, c2, qh);

        /* Check: g^s * pk^e == R mod p */
        bn_mxp(lhs, g,  sig->s, p);   /* g^s     */
        bn_mxp(rhs, pk, e,      p);   /* pk^e    */
        bn_mul(lhs, lhs, rhs);
        bn_mod(lhs, lhs, p);           /* g^s * pk^e */

        ok = (bn_cmp(lhs, sig->R) == RLC_EQ);

    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(e); bn_free(lhs); bn_free(rhs);
    }
    return (int)ok;
}
