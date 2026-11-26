#include "utils.h"

/* utils.c -> Arithmetic utilities for the ECT scheme.  */

 
/* Group arithmetic */
 

void compute_qh(bn_t qh, const bn_t q, int h)
{
    /* qh = q^H.  Start with q and multiply h-1 more times.*/
    bn_copy(qh, q);
    for (int i = 1; i < h; i++)
        bn_mul(qh, qh, q);
}

void ciphertext_pow_q(bn_t c1, bn_t c2, const bn_t q, const bn_t p)
{
    /* Raise both components of an ElGamal ciphertext to the q-th power. */
    bn_mxp(c1, c1, q, p);
    bn_mxp(c2, c2, q, p);
}

 
/* Identifier encoding */
 

void user_identifier(bn_t m, int id, int r)
{
    /* m_id = r^(id-1). */
    bn_set_dig(m, 1);
    for (int i = 1; i < id; i++)
        bn_mul_dig(m, m, (dig_t)r);
}

 
/* Plaintext decomposition  */
 

void base_q_decompose(bn_t *d, const bn_t v, const bn_t q, int h)
{
    /* Extract base-q digits: v = d[0] + d[1]*q + d[2]*q^2 + ... */
    bn_t tmp, quot;
    bn_null(tmp); bn_null(quot);
    RLC_TRY {
        bn_new(tmp); bn_new(quot);
        bn_copy(tmp, v);
        for (int t = 0; t < h; t++) {
            bn_div_rem(quot, d[t], tmp, q);
            bn_copy(tmp, quot);
        }
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(tmp); bn_free(quot);
    }
}

void base_r_decompose(uint32_t *c, const bn_t v, int r, int n)
{
    /* Extract base-r digits: v = c[0] + c[1]*r + c[2]*r^2 + ...
      For an honestly generated set tag, c[i] is 0 or 1.  Larger values
      represent multiplicities created by merges.
     */
    bn_t tmp, base, quot, rem;
    bn_null(tmp); bn_null(base); bn_null(quot); bn_null(rem);
    RLC_TRY {
        bn_new(tmp); bn_new(base); bn_new(quot); bn_new(rem);
        bn_copy(tmp, v);
        bn_set_dig(base, (dig_t)r);
        for (int i = 0; i < n; i++) {
            dig_t digit;
            bn_div_rem(quot, rem, tmp, base);
            bn_get_dig(&digit, rem);
            c[i] = (uint32_t)digit;
            bn_copy(tmp, quot);
        }
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(tmp); bn_free(base); bn_free(quot); bn_free(rem);
    }
}
