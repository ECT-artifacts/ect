#include "sig.h"
#include <openssl/err.h>
#include <stdio.h>
#include <stdlib.h>

/* sig.c -> one-shot signatures over the whole message (EVP_DigestSign with
   no external digest, which is the interface of both Ed25519 and ML-DSA). */

static void die(const char *what)
{
    fprintf(stderr, "%s failed\n", what);
    ERR_print_errors_fp(stderr);
    abort();
}

void sig_keygen(SigKey *k, const char *alg)
{
    k->key = EVP_PKEY_Q_keygen(NULL, NULL, alg);
    if (!k->key) die("EVP_PKEY_Q_keygen");
    k->sig_bytes = (size_t)EVP_PKEY_get_size(k->key);
}

void sig_free(SigKey *k)
{
    EVP_PKEY_free(k->key);
    k->key = NULL;
}

size_t sig_sign(const SigKey *k, const uint8_t *msg, size_t len,
                uint8_t *sig, size_t cap)
{
    EVP_MD_CTX *ctx = EVP_MD_CTX_new();
    size_t siglen = cap;
    if (!ctx) die("EVP_MD_CTX_new");
    if (EVP_DigestSignInit_ex(ctx, NULL, NULL, NULL, NULL, k->key, NULL) != 1)
        die("EVP_DigestSignInit_ex");
    if (EVP_DigestSign(ctx, sig, &siglen, msg, len) != 1)
        die("EVP_DigestSign");
    EVP_MD_CTX_free(ctx);
    return siglen;
}

int sig_verify(const SigKey *k, const uint8_t *msg, size_t len,
               const uint8_t *sig, size_t siglen)
{
    EVP_MD_CTX *ctx = EVP_MD_CTX_new();
    if (!ctx) die("EVP_MD_CTX_new");
    if (EVP_DigestVerifyInit_ex(ctx, NULL, NULL, NULL, NULL, k->key, NULL) != 1)
        die("EVP_DigestVerifyInit_ex");
    int ok = EVP_DigestVerify(ctx, sig, siglen, msg, len) == 1;
    EVP_MD_CTX_free(ctx);
    ERR_clear_error();
    return ok;
}
