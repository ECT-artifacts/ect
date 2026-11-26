#pragma once
#include <relic/relic.h>

/* elgamal.h -> Exponential ElGamal encryption over Z*_p. */

/* elgamal_keygen -> generate a fresh ElGamal key pair. */
void elgamal_keygen(bn_t sk, bn_t pk,
                    const bn_t g, const bn_t qh, const bn_t p);

/*  elgamal_enc -> encrypt plaintext m under public key pk. */
void elgamal_enc(bn_t c1, bn_t c2, bn_t r_out,
                 const bn_t m, const bn_t pk,
                 const bn_t g, const bn_t qh, const bn_t p);

/* elgamal_dec_gv -> decrypt (c1, c2) to recover the group element g^v. */
void elgamal_dec_gv(bn_t gv,
                    const bn_t c1, const bn_t c2, const bn_t sk,
                    const bn_t p);

/* elgamal_dec -> decrypt (c1, c2) to recover m via brute-force DLOG. */
int elgamal_dec(bn_t m,
                const bn_t c1, const bn_t c2, const bn_t sk,
                const bn_t g, const bn_t p);

/* elgamal_rerand -> rerandomize ciphertext (c1, c2) in-place. */
void elgamal_rerand(bn_t c1, bn_t c2, bn_t r_out,
                    const bn_t pk,
                    const bn_t g, const bn_t qh, const bn_t p);

/* elgamal_mul -> homomorphic ciphertext multiplication  */
void elgamal_mul(bn_t c1, bn_t c2,
                 const bn_t d1, const bn_t d2,
                 const bn_t p);
