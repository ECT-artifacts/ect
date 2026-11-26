#include "prg.h"
#include <openssl/rand.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* prg.c -> batched RAND_bytes, with one pool per thread. */

#define PRG_BATCH 65536

static _Thread_local uint8_t pool[PRG_BATCH];
static _Thread_local size_t  pool_pos = PRG_BATCH;

static void refill(void)
{
    if (RAND_bytes(pool, PRG_BATCH) != 1) {
        fprintf(stderr, "RAND_bytes failed\n");
        abort();
    }
    pool_pos = 0;
}

void rnd_bytes(uint8_t *buf, size_t len)
{
    if (len >= PRG_BATCH) {
        if (RAND_bytes(buf, (int)len) != 1) {
            fprintf(stderr, "RAND_bytes failed\n");
            abort();
        }
        return;
    }
    while (len > 0) {
        if (pool_pos == PRG_BATCH) refill();
        size_t take = PRG_BATCH - pool_pos;
        if (take > len) take = len;
        memcpy(buf, pool + pool_pos, take);
        pool_pos += take;
        buf += take;
        len -= take;
    }
}

uint64_t rnd_u64(void)
{
    uint64_t x;
    rnd_bytes((uint8_t *)&x, sizeof(x));
    return x;
}

uint64_t rnd_below(uint64_t bound)
{
    if (bound <= 1) return 0;
    /* Reject the top partial interval so the result is exactly uniform. */
    uint64_t limit = UINT64_MAX - (UINT64_MAX % bound);
    for (;;) {
        uint64_t x = rnd_u64();
        if (x < limit) return x % bound;
    }
}
