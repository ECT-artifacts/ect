#pragma once
#include <relic/relic.h>
#include <stdint.h>
#include "params.h"
#include "nizk.h"

/*
 * block.h -> Fixed-size block and append-only ledger for the ECT.
 */
/* Group elements live in Z*_p, so P_BYTES each (726 at the bundled 5807-bit
   p).  Derived from params.h so new parameters cannot silently overflow. */
#define WM_ELEM_BYTES   P_BYTES

/* The NIZK scalars z and e are elements of Z_{q^H}, so q^H needs
   ceil(H*Q_BITS/8) bytes; the +1 keeps a byte of slack.  Derived from H and
   Q_BITS so that changing either cannot silently overflow the buffer. */
#define WM_Z_BYTES      ((((H) * (Q_BITS) + 7) / 8) + 1)

#define WM_OFF_C1       0
#define WM_OFF_C2       (WM_OFF_C1      + WM_ELEM_BYTES)
#define WM_OFF_A1       (WM_OFF_C2      + WM_ELEM_BYTES)
#define WM_OFF_A2       (WM_OFF_A1      + WM_ELEM_BYTES)
#define WM_OFF_Z        (WM_OFF_A2      + WM_ELEM_BYTES)
#define WM_OFF_ECLAUSE  (WM_OFF_Z       + WM_Z_BYTES)
#define WM_OFF_HOPS     (WM_OFF_ECLAUSE + WM_Z_BYTES)
#define WM_OFF_HAS      (WM_OFF_HOPS    + 1)
#define WM_OFF_C1_0     (WM_OFF_HAS     + 1)
#define WM_OFF_C2_0     (WM_OFF_C1_0    + WM_ELEM_BYTES)
#define WATERMARK_BYTES (WM_OFF_C2_0    + WM_ELEM_BYTES)

typedef struct {
    uint8_t data[WATERMARK_BYTES];
} Watermark;


/*
  Block byte layout.  BLOCK_BYTES is derived from the watermark rather than
  fixed, so a larger group does not silently overflow the block.
 */
#define HASH_BYTES     32   /* SHA-256 output length */

#define BLK_OFF_PREV_HASH    0
#define BLK_OFF_TX_HASH     32
#define BLK_OFF_HEIGHT      64
#define BLK_OFF_SENDER      68
#define BLK_OFF_RECEIVER    72
#define BLK_OFF_AMOUNT      76
#define BLK_OFF_WATERMARK   84
#define BLK_PADDING_START  (BLK_OFF_WATERMARK + WATERMARK_BYTES)

/* Padded up to the next 4 KiB boundary. */
#define BLOCK_BYTES  ((((BLK_PADDING_START) + 4095) / 4096) * 4096)

/* Compile-time check: watermark + header must fit within BLOCK_BYTES */
_Static_assert(BLK_OFF_WATERMARK + WATERMARK_BYTES <= BLOCK_BYTES,
               "Watermark overflows block -> increase BLOCK_BYTES");

typedef struct {
    uint8_t raw[BLOCK_BYTES];
} Block;

 
/* Ledger */
 

/*
    The ledger is a flat array of Block structs linked by their prev_hash
    fields.  It is append-only: once committed, a block is never modified.

 */
#define MAX_BLOCKS 1024

typedef struct {
    Block    blocks[MAX_BLOCKS];
    uint32_t height;   /* number of blocks appended so far (next index) */
} Ledger;

 
/* Watermark    */
 

/*
    wm_pack -> serialise ciphertext, optional NIZK proof, and original
    ciphertext (c1_0, c2_0) into a Watermark.
 */
void wm_pack(Watermark *wm,
             const bn_t c1, const bn_t c2,
             const NIZKProof *pi,
             uint8_t hop_count,
             const bn_t c1_0, const bn_t c2_0);

/*
    wm_unpack -> deserialize a Watermark.
    
    c1_0, c2_0 : output for original ciphertext fields; pass NULL to skip.
 */
void wm_unpack(bn_t c1, bn_t c2,
               NIZKProof *pi,
               uint8_t *hop_count,
               bn_t c1_0, bn_t c2_0,
               const Watermark *wm);

 
/* Block                                             */
 

/*
    block_set_header -> write scalar fields into block->raw.

    prev_hash : SHA-256 of the previous block, or all-zero for init.
    height    : block index (0-based)
    sender_id, receiver_id : 1-indexed user IDs
    amount    : coin value

 */
void block_set_header(Block *b,
                      const uint8_t prev_hash[HASH_BYTES],
                      uint32_t height,
                      uint32_t sender_id,
                      uint32_t receiver_id,
                      uint64_t amount);

/* block_set_watermark -> copy a serialized Watermark into the block. */
void block_set_watermark(Block *b, const Watermark *wm);

/*
    block_compute_tx_hash -> compute and store the transaction hash.

    tx_hash = SHA-256(sender_id || receiver_id || amount || watermark)

 */
void block_compute_tx_hash(Block *b);

 
/* Block  */
 

uint32_t block_get_height   (const Block *b);
uint32_t block_get_sender   (const Block *b);
uint32_t block_get_receiver (const Block *b);
uint64_t block_get_amount   (const Block *b);
void     block_get_watermark(Watermark *wm, const Block *b);

/*
    block_hash -> compute SHA-256 of the entire 4096-byte block.
 */
void block_hash(uint8_t out[HASH_BYTES], const Block *b);

/* Print a one-line human-readable summary to stdout. */
void block_print(const Block *b);

 
/* Ledger */
 

void    ledger_init(Ledger *l);
Ledger *ledger_new (void);   /* heap-allocate; always use this, not stack */
void    ledger_free(Ledger *l);

/*
    ledger_append -> finalize and append a block to the ledger.
 
    returns 1 on success, 0 if the ledger is full (height >= MAX_BLOCKS).
 */
int  ledger_append(Ledger *l, Block *b);

/* Print a summary of every block in the ledger. */
void ledger_print(const Ledger *l);
