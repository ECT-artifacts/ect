#define _POSIX_C_SOURCE 200809L
#include "board_util.h"
#include <openssl/evp.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

/* board_util.c -> Merkle root, signed message and parallel loop. */

static void sha256(uint8_t out[32], const uint8_t *a, size_t alen,
                   const uint8_t *b, size_t blen, const uint8_t *c, size_t clen)
{
    EVP_MD_CTX *ctx = EVP_MD_CTX_new();
    unsigned int len = 0;
    if (!ctx || EVP_DigestInit_ex(ctx, EVP_sha256(), NULL) != 1) abort();
    EVP_DigestUpdate(ctx, a, alen);
    if (blen) EVP_DigestUpdate(ctx, b, blen);
    if (clen) EVP_DigestUpdate(ctx, c, clen);
    EVP_DigestFinal_ex(ctx, out, &len);
    EVP_MD_CTX_free(ctx);
}

void merkle_root(const uint8_t (*leaves)[32], size_t M, uint8_t root[32])
{
    uint8_t (*lv)[32] = malloc(M * 32);
    if (!lv) abort();
    memcpy(lv, leaves, M * 32);
    static const uint8_t inner = 0x01;
    while (M > 1) {
        size_t k = 0;
        for (size_t i = 0; i + 1 < M; i += 2)
            sha256(lv[k++], &inner, 1, lv[i], 32, lv[i + 1], 32);
        if (M % 2) memmove(lv[k++], lv[M - 1], 32);
        M = k;
    }
    memcpy(root, lv[0], 32);
    free(lv);
}

void board_message(uint8_t *buf, uint64_t M, const uint8_t root[32])
{
    memset(buf, 0, 16);
    memcpy(buf, "ECT-board-root", 14);
    for (int i = 0; i < 8; i++) buf[16 + i] = (uint8_t)(M >> (56 - 8 * i));
    memcpy(buf + 24, root, 32);
}

typedef struct {
    board_range_fn fn;
    void *ctx;
    size_t M;
    atomic_size_t next;
} Pool;

typedef struct { Pool *pool; int tid; } Job;

/* Each thread takes the next entry until none is left, so that the threads
   on the faster cores do more entries. */
static void *run_job(void *arg)
{
    Job *j = arg;
    Pool *P = j->pool;
    for (;;) {
        size_t i = atomic_fetch_add(&P->next, 1);
        if (i >= P->M) break;
        P->fn(P->ctx, j->tid, i, i + 1);
    }
    return NULL;
}

void board_parallel(int threads, size_t M, board_range_fn fn, void *ctx)
{
    if (threads <= 1 || M <= 1) { fn(ctx, 0, 0, M); return; }
    if ((size_t)threads > M) threads = (int)M;
    Pool P = { fn, ctx, M, 0 };
    Job *jobs = malloc((size_t)threads * sizeof(Job));
    pthread_t *tid = malloc((size_t)threads * sizeof(pthread_t));
    if (!jobs || !tid) abort();
    for (int t = 0; t < threads; t++) {
        jobs[t] = (Job){ &P, t };
        if (pthread_create(&tid[t], NULL, run_job, &jobs[t]) != 0) abort();
    }
    for (int t = 0; t < threads; t++) pthread_join(tid[t], NULL);
    free(jobs); free(tid);
}
