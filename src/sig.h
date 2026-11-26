#pragma once
#include <stddef.h>
#include <stdint.h>
#include <openssl/evp.h>

/* sig.h -> issuance signatures for the Damgard-Jurik and Ring-LWE prototypes,
   through OpenSSL's EVP interface.  The Damgard-Jurik prototype uses Ed25519
   and the Ring-LWE prototype uses ML-DSA-44 (FIPS 204), so that its issuance
   signature is post-quantum as well. */

typedef struct {
    EVP_PKEY *key;
    size_t    sig_bytes;   /* maximal signature length */
} SigKey;

/* alg is an OpenSSL key type name, "ED25519" or "ML-DSA-44". */
void   sig_keygen(SigKey *k, const char *alg);
void   sig_free(SigKey *k);
size_t sig_sign(const SigKey *k, const uint8_t *msg, size_t len,
                uint8_t *sig, size_t cap);
int    sig_verify(const SigKey *k, const uint8_t *msg, size_t len,
                  const uint8_t *sig, size_t siglen);
