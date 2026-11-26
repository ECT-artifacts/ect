#pragma once
#include <stddef.h>
#include <stdint.h>

/* prg.h -> randomness for the Damgard-Jurik and Ring-LWE prototypes.
   All bytes come from OpenSSL's RAND_bytes (a CTR-DRBG seeded by the
   operating system), drawn in 64 KiB batches. */

void     rnd_bytes(uint8_t *buf, size_t len);
uint64_t rnd_u64(void);

/* Uniform integer in [0, bound) by rejection sampling, bound >= 1. */
uint64_t rnd_below(uint64_t bound);
