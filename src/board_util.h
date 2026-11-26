#pragma once
#include <stddef.h>
#include <stdint.h>

/* board_util.h -> parts shared by the exponential ElGamal and Damgard-Jurik
   board proofs: the Merkle root of the entries, the message that the
   authority signs, and a parallel loop over the entries. */

/* Merkle root under SHA-256 of M >= 1 leaves of 32 bytes.  An inner node is
   SHA-256(0x01 || left || right); a node without a sibling moves up
   unchanged.  The leaves are not modified. */
void merkle_root(const uint8_t (*leaves)[32], size_t M, uint8_t root[32]);

/* The message signed for a board: a label, the number of entries and the
   root.  buf must hold BOARD_MSG_BYTES bytes. */
#define BOARD_MSG_BYTES (16 + 8 + 32)
void board_message(uint8_t *buf, uint64_t M, const uint8_t root[32]);

/* Runs fn(ctx, tid, i, i+1) for every i in [0, M) on threads 0..threads-1,
   each thread taking the next i, and waits for all of them.  threads <= 1
   runs fn(ctx, 0, 0, M) in the calling thread. */
typedef void (*board_range_fn)(void *ctx, int tid, size_t lo, size_t hi);
void board_parallel(int threads, size_t M, board_range_fn fn, void *ctx);
