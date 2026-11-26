#include "elgamal.h"
#include "params.h"

/*  elgamal.c -> Exponential ElGamal encryption over Z*_p. */

/* Key generation   */
 

void elgamal_keygen(bn_t sk, bn_t pk,
                    const bn_t g, const bn_t qh, const bn_t p)
{
    /* sk ← Z_{q^H} uniformly at random */
    bn_rand_mod(sk, qh);
    /* pk = g^sk mod p */
    bn_mxp(pk, g, sk, p);
}

 
/* Encryption   */
 

void elgamal_enc(bn_t c1, bn_t c2, bn_t r_out,
                 const bn_t m, const bn_t pk,
                 const bn_t g, const bn_t qh, const bn_t p)
{
    bn_t gm;
    bn_null(gm);
    RLC_TRY {
        bn_new(gm);

        /* r ← Z_{q^H} (fresh randomness, also returned for NIZK use) */
        bn_rand_mod(r_out, qh);

        /* c1 = g^r mod p  (ephemeral key -> hides the randomness) */
        bn_mxp(c1, g, r_out, p);

        /* c2 = pk^r * g^m mod p  (masked plaintext) */
        bn_mxp(c2,  pk, r_out, p);   /* pk^r */
        bn_mxp(gm,  g,  m,     p);   /* g^m  */
        bn_mul(c2, c2, gm);
        bn_mod(c2, c2, p);

    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(gm);
    }
}

 
/* Decryption -> group element only */
 

void elgamal_dec_gv(bn_t gv,
                    const bn_t c1, const bn_t c2, const bn_t sk,
                    const bn_t p)
{
    bn_t s;
    bn_null(s);
    RLC_TRY {
        bn_new(s);
        /* gv = c2 * (c1^sk)^{-1} mod p */
        bn_mxp(s,  c1, sk, p);   /* s  = c1^sk = pk^r */
        bn_mod_inv(s, s, p);      /* s  = (pk^r)^{-1} */
        bn_mul(gv, c2, s);
        bn_mod(gv, gv, p);        /* gv = g^v */
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(s);
    }
}

 
/* Decryption -> brute-force DLOG for small stuff */
 

int elgamal_dec(bn_t m,
                const bn_t c1, const bn_t c2, const bn_t sk,
                const bn_t g, const bn_t p)
{
    volatile int ret = -1;
    bn_t s, gm, t;
    bn_null(s); bn_null(gm); bn_null(t);
    RLC_TRY {
        bn_new(s); bn_new(gm); bn_new(t);

        /* Recover g^m = c2 / c1^sk = c2 * (c1^sk)^{-1} mod p */
        bn_mxp(s,  c1, sk, p);   /* s  = c1^sk = g^(r*sk) = pk^r */
        bn_mod_inv(s, s, p);      /* s  = (pk^r)^{-1} */
        bn_mul(gm, c2, s);
        bn_mod(gm, gm, p);        /* gm = c2 * s = g^m */

        /* Special case: g^0 = 1 -> m = 0 */
        bn_set_dig(t, 1);
        if (bn_cmp(gm, t) == RLC_EQ) {
            bn_zero(m);
            ret = 0;
        } else {
            /* Brute-force DLOG: find i such that g^i ≡ gm (mod p).
             * t is maintained as g^i via incremental multiplication,
             * avoiding a full modular exponentiation per iteration. */
            bn_set_dig(t, 1);
            for (int i = 1; i < DLOG_MAX; i++) {
                bn_mul(t, t, g);
                bn_mod(t, t, p);
                if (bn_cmp(t, gm) == RLC_EQ) {
                    bn_set_dig(m, (dig_t)i);
                    ret = 0;
                    break;
                }
            }
        }
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(s); bn_free(gm); bn_free(t);
    }
    return ret;
}

 
/* Rerandomization  */
 

void elgamal_rerand(bn_t c1, bn_t c2, bn_t r_out,
                    const bn_t pk,
                    const bn_t g, const bn_t qh, const bn_t p)
{
    bn_t tmp;
    bn_null(tmp);
    RLC_TRY {
        bn_new(tmp);

        /* r' ← Z_{q^H} */
        bn_rand_mod(r_out, qh);

        /* c1 ← c1 * g^r' mod p */
        bn_mxp(tmp, g, r_out, p);
        bn_mul(c1, c1, tmp);
        bn_mod(c1, c1, p);

        /* c2 ← c2 * pk^r' mod p */
        bn_mxp(tmp, pk, r_out, p);
        bn_mul(c2, c2, tmp);
        bn_mod(c2, c2, p);

    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(tmp);
    }
}

 
/* Homomorphic multiplication  */
 

void elgamal_mul(bn_t c1, bn_t c2,
                 const bn_t d1, const bn_t d2,
                 const bn_t p)
{
    /* (c1,c2) <- (c1*d1, c2*d2) mod p
    Decrypts to the sum of the two plaintexts: m1 + m2.*/
    bn_mul(c1, c1, d1); bn_mod(c1, c1, p);
    bn_mul(c2, c2, d2); bn_mod(c2, c2, p);
}
