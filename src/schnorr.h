#pragma once
#include <relic/relic.h>
#include "params.h"

/*
 * schnorr.h -> Schnorr signature scheme for authority tag authentication */

typedef struct {
    bn_t R;   /* commitment g^k mod p  (Z*_p element, <= P_BYTES) */
    bn_t s;   /* response   k - e*sk  mod q^H  (scalar, <= 160 bytes) */
} SchnorrSig;

void schnorr_sig_new (SchnorrSig *sig);
void schnorr_sig_free(SchnorrSig *sig);

/* schnorr_keygen -> generate a Schnorr signing key pair. */
void schnorr_keygen(bn_t sk, bn_t pk,
                    const bn_t g, const bn_t qh, const bn_t p);

/* schnorr_sign -> sign the ciphertext (c1, c2) under key (sk, pk). */
void schnorr_sign(SchnorrSig *sig,
                  const bn_t sk, const bn_t pk,
                  const bn_t c1, const bn_t c2,
                  const bn_t g, const bn_t qh, const bn_t p);

/* schnorr_verify -> verify a Schnorr signature on (c1, c2). */
int schnorr_verify(const SchnorrSig *sig,
                   const bn_t pk,
                   const bn_t c1, const bn_t c2,
                   const bn_t g, const bn_t qh, const bn_t p);
