#include "block.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* block.c -> Block and ledger implementation (All serialization uses explicit big-endian byte)*/



/*  write_u32 / write_u64 -> write an integer into a byte buffer, 
    most-significant byte first (big-endian / network byte order).
*/

static void write_u32(uint8_t *dst, uint32_t v)
{
    dst[0] = (uint8_t)(v >> 24);
    dst[1] = (uint8_t)(v >> 16);
    dst[2] = (uint8_t)(v >>  8);
    dst[3] = (uint8_t)(v);
}

static void write_u64(uint8_t *dst, uint64_t v)
{
    dst[0] = (uint8_t)(v >> 56);
    dst[1] = (uint8_t)(v >> 48);
    dst[2] = (uint8_t)(v >> 40);
    dst[3] = (uint8_t)(v >> 32);
    dst[4] = (uint8_t)(v >> 24);
    dst[5] = (uint8_t)(v >> 16);
    dst[6] = (uint8_t)(v >>  8);
    dst[7] = (uint8_t)(v);
}

static uint32_t read_u32(const uint8_t *src)
{
    return ((uint32_t)src[0] << 24) | ((uint32_t)src[1] << 16)
         | ((uint32_t)src[2] <<  8) |  (uint32_t)src[3];
}

static uint64_t read_u64(const uint8_t *src)
{
    return ((uint64_t)src[0] << 56) | ((uint64_t)src[1] << 48)
         | ((uint64_t)src[2] << 40) | ((uint64_t)src[3] << 32)
         | ((uint64_t)src[4] << 24) | ((uint64_t)src[5] << 16)
         | ((uint64_t)src[6] <<  8) |  (uint64_t)src[7];
}

/* bn_to_slot -> serialize a RELIC big integer into a fixed-width byte slot. */
static void bn_to_slot(uint8_t *dst, size_t len, const bn_t x)
{
    memset(dst, 0, len);
    int actual = bn_size_bin(x);
    if (actual <= 0) return;
    if ((size_t)actual > len) actual = (int)len;  /* truncation guard */
    bn_write_bin(dst + len - (size_t)actual, actual, x);
}

/*
    slot_to_bn -> deserialize a fixed-width byte slot into a RELIC big integer.
    bn_read_bin interprets the bytes as a big-endian unsigned integer.
 */
static void slot_to_bn(bn_t x, const uint8_t *src, size_t len)
{
    bn_read_bin(x, src, len);
}

 
/* Watermark */
 

void wm_pack(Watermark *wm,
             const bn_t c1, const bn_t c2,
             const NIZKProof *pi,
             uint8_t hop_count,
             const bn_t c1_0, const bn_t c2_0)
{
    memset(wm->data, 0, WATERMARK_BYTES);

    /* Current ciphertext */
    bn_to_slot(wm->data + WM_OFF_C1, WM_ELEM_BYTES, c1);
    bn_to_slot(wm->data + WM_OFF_C2, WM_ELEM_BYTES, c2);

    /* CDS proof (clause 0) */
    if (pi != NULL) {
        bn_to_slot(wm->data + WM_OFF_A1,      WM_ELEM_BYTES, pi->A1[0]);
        bn_to_slot(wm->data + WM_OFF_A2,      WM_ELEM_BYTES, pi->A2[0]);
        bn_to_slot(wm->data + WM_OFF_Z,       WM_Z_BYTES,    pi->z[0]);
        bn_to_slot(wm->data + WM_OFF_ECLAUSE, WM_Z_BYTES,    pi->e[0]);
        wm->data[WM_OFF_HAS] = 1;
    }

    wm->data[WM_OFF_HOPS] = hop_count;

    /*  Original ciphertext at issuance -> carried in every block.
        Enables validators to verify the authority sig and first-hop NIZK
        without accessing the ECTTag struct. */
    if (c1_0 != NULL) bn_to_slot(wm->data + WM_OFF_C1_0, WM_ELEM_BYTES, c1_0);
    if (c2_0 != NULL) bn_to_slot(wm->data + WM_OFF_C2_0, WM_ELEM_BYTES, c2_0);
}

void wm_unpack(bn_t c1, bn_t c2,
               NIZKProof *pi,
               uint8_t *hop_count,
               bn_t c1_0, bn_t c2_0,
               const Watermark *wm)
{
    if (c1) slot_to_bn(c1, wm->data + WM_OFF_C1, WM_ELEM_BYTES);
    if (c2) slot_to_bn(c2, wm->data + WM_OFF_C2, WM_ELEM_BYTES);

    if (hop_count) *hop_count = wm->data[WM_OFF_HOPS];

    if (pi) {
        slot_to_bn(pi->A1[0], wm->data + WM_OFF_A1,      WM_ELEM_BYTES);
        slot_to_bn(pi->A2[0], wm->data + WM_OFF_A2,      WM_ELEM_BYTES);
        slot_to_bn(pi->z[0],  wm->data + WM_OFF_Z,       WM_Z_BYTES);
        slot_to_bn(pi->e[0],  wm->data + WM_OFF_ECLAUSE, WM_Z_BYTES);
    }

    if (c1_0) slot_to_bn(c1_0, wm->data + WM_OFF_C1_0, WM_ELEM_BYTES);
    if (c2_0) slot_to_bn(c2_0, wm->data + WM_OFF_C2_0, WM_ELEM_BYTES);
}

 
/* Block construction */
 

void block_set_header(Block *b,
                      const uint8_t prev_hash[HASH_BYTES],
                      uint32_t height,
                      uint32_t sender_id,
                      uint32_t receiver_id,
                      uint64_t amount)
{
    memset(b->raw, 0, BLOCK_BYTES);
    memcpy(b->raw + BLK_OFF_PREV_HASH, prev_hash, HASH_BYTES);
    write_u32(b->raw + BLK_OFF_HEIGHT,   height);
    write_u32(b->raw + BLK_OFF_SENDER,   sender_id);
    write_u32(b->raw + BLK_OFF_RECEIVER, receiver_id);
    write_u64(b->raw + BLK_OFF_AMOUNT,   amount);
}

void block_set_watermark(Block *b, const Watermark *wm)
{
    memcpy(b->raw + BLK_OFF_WATERMARK, wm->data, WATERMARK_BYTES);
}

void block_compute_tx_hash(Block *b)
{
    /*
        tx_hash commits to the transaction content: SHA-256(sender_id[4] || receiver_id[4] || amount[8] || watermark[2048])
     */
    uint8_t digest[HASH_BYTES];
    size_t  input_len = 4 + 4 + 8 + WATERMARK_BYTES;
    md_map(digest, b->raw + BLK_OFF_SENDER, input_len);
    memcpy(b->raw + BLK_OFF_TX_HASH, digest, HASH_BYTES);
}

 
/* Block read  */
 

uint32_t block_get_height   (const Block *b) { return read_u32(b->raw + BLK_OFF_HEIGHT);   }
uint32_t block_get_sender   (const Block *b) { return read_u32(b->raw + BLK_OFF_SENDER);   }
uint32_t block_get_receiver (const Block *b) { return read_u32(b->raw + BLK_OFF_RECEIVER); }
uint64_t block_get_amount   (const Block *b) { return read_u64(b->raw + BLK_OFF_AMOUNT);   }

void block_get_watermark(Watermark *wm, const Block *b)
{
    memcpy(wm->data, b->raw + BLK_OFF_WATERMARK, WATERMARK_BYTES);
}

void block_hash(uint8_t out[HASH_BYTES], const Block *b)
{
    /* Hash the entire raw block -> used to build the prev_hash chain */
    md_map(out, b->raw, BLOCK_BYTES);
}

void block_print(const Block *b)
{
    uint8_t h[HASH_BYTES];
    block_hash(h, b);
    printf("block #%u  sender=%u -> receiver=%u  amount=%llu  "
           "hops=%u  hash=%02x%02x..%02x%02x\n",
           block_get_height(b),
           block_get_sender(b),
           block_get_receiver(b),
           (unsigned long long)block_get_amount(b),
           b->raw[BLK_OFF_WATERMARK + WM_OFF_HOPS],
           h[0], h[1], h[HASH_BYTES-2], h[HASH_BYTES-1]);
}

 
/* Ledger */
 

void ledger_init(Ledger *l)
{
    memset(l, 0, sizeof(Ledger));
    l->height = 0;
}

Ledger *ledger_new(void)
{
    Ledger *l = (Ledger *)malloc(sizeof(Ledger));
    if (l) ledger_init(l);
    return l;
}

void ledger_free(Ledger *l)
{
    free(l);
}

int ledger_append(Ledger *l, Block *b)
{
    if (l->height >= MAX_BLOCKS) return 0;

    /* Set prev_hash: zero for init, otherwise SHA-256 of previous block */
    if (l->height == 0) {
        memset(b->raw + BLK_OFF_PREV_HASH, 0, HASH_BYTES);
    } else {
        block_hash(b->raw + BLK_OFF_PREV_HASH, &l->blocks[l->height - 1]);
    }

    /* Write the block's position in the chain */
    write_u32(b->raw + BLK_OFF_HEIGHT, l->height);

    /* Compute and store the transaction content hash */
    block_compute_tx_hash(b);

    /* Append */
    l->blocks[l->height] = *b;
    l->height++;
    return 1;
}

void ledger_print(const Ledger *l)
{
    printf("=== Ledger (%u blocks) ===\n", l->height);
    for (uint32_t i = 0; i < l->height; i++)
        block_print(&l->blocks[i]);
}
