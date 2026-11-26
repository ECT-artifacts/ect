#pragma once
#include <relic/relic.h>
#include <stdint.h>

/*
  utils.h -> Arithmetic utilities 
 */
void compute_qh(bn_t qh, const bn_t q, int h);

/*  ciphertext_pow_q -> raise an ElGamal ciphertext to the q-th power  */
void ciphertext_pow_q(bn_t c1, bn_t c2, const bn_t q, const bn_t p);

/*  user_identifier -> compute the base-r identifier for user id  */
void user_identifier(bn_t m, int id, int r);

/*  base_q_decompose -> extract the base-q digit representation of v.  */
void base_q_decompose(bn_t *d, const bn_t v, const bn_t q, int h);

/*  base_r_decompose -> extract the base-r digit representation of v. */
void base_r_decompose(uint32_t *c, const bn_t v, int r, int n);
