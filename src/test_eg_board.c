#include "eg_board.h"
#include "board_util.h"
#include "elgamal.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* test_eg_board.c -> accept/reject tests for the exponential ElGamal board
   proof and for the rerandomized Merge.  Usage: ./test_eg_board */

static int failures = 0;

#define EXPECT(cond, ...) do {                                   \
    if (cond) { printf("ok    "); }                             \
    else { printf("FAIL  "); failures++; }                      \
    printf(__VA_ARGS__); printf("\n"); fflush(stdout);          \
} while (0)

int main(void)
{
    if (core_init() != RLC_OK) return 1;
    ECTParams pp;
    ect_params_new(&pp);
    ect_setup(&pp);
    printf("exponential ElGamal board: h=%d n=%d |p|=%d\n", H, N, P_BITS);

    /* Rerandomized merge keeps the plaintext sum and changes the ciphertext. */
    {
        ECTTag a, b, m1, m2;
        ect_tag_new(&a); ect_tag_new(&b); ect_tag_new(&m1); ect_tag_new(&m2);
        ect_tag(&a, &pp, 4, 1);
        ect_tag(&b, &pp, 9, 1);
        ECTTag *in[2] = { &a, &b };
        ect_merge(&m1, in, 2, &pp);
        ect_merge(&m2, in, 2, &pp);
        EXPECT(bn_cmp(m1.c1, m2.c1) != RLC_EQ, "two merges of the same inputs differ");
        EXPECT(ect_trace_silent(&pp, &m1) == 2 && ect_trace_silent(&pp, &m2) == 2,
              "both merges trace two identifiers");
        ect_tag_free(&a); ect_tag_free(&b); ect_tag_free(&m1); ect_tag_free(&m2);
    }

    EGBoardParams bp;
    eg_board_params_init(&bp, &pp, 6);
    mpz_t t;
    mpz_init(t);

    /* Every identifier, and the zero entry, gives an entry that verifies. */
    {
        EGEntry E;
        eg_entry_init(&bp, &E);
        int all = 1;
        for (int s = 1; s <= bp.n; s++) {
            eg_entry_prove(&bp, &E, 7, s, 1, NULL);
            all &= eg_entry_verify(&bp, &E, 7);
        }
        EXPECT(all, "a nonzero entry verifies for each of the %d identifiers", bp.n);
        eg_entry_prove(&bp, &E, 7, 0, 0, NULL);
        EXPECT(eg_entry_verify(&bp, &E, 7), "a zero entry verifies");
        EXPECT(!eg_entry_verify(&bp, &E, 8), "an entry verified at another position is rejected");

        eg_entry_prove(&bp, &E, 7, 3, 1, NULL);
        /* e2 * g^{m_3} encrypts 2 m_3, which no identifier encodes */
        mpz_powm(t, bp.g, bp.m[2], bp.p);
        mpz_mul(E.e2, E.e2, t); mpz_mod(E.e2, E.e2, bp.p);
        EXPECT(!eg_entry_verify(&bp, &E, 7), "an entry changed to encrypt 2 m_3 is rejected");
        eg_entry_prove(&bp, &E, 7, 3, 1, NULL);
        mpz_add_ui(E.z[3 + 2], E.z[3 + 2], 1);
        EXPECT(!eg_entry_verify(&bp, &E, 7), "a modified response is rejected");
        eg_entry_prove(&bp, &E, 7, 3, 1, NULL);
        mpz_mul(E.T[5][1], E.T[5][1], bp.g); mpz_mod(E.T[5][1], E.T[5][1], bp.p);
        EXPECT(!eg_entry_verify(&bp, &E, 7), "a modified first message is rejected");
        eg_entry_prove(&bp, &E, 7, 3, 1, NULL);
        /* e1 * (p-1): the subgroup component is unchanged, so only the
           membership test can reject it */
        mpz_sub_ui(t, bp.p, 1);
        mpz_mul(E.e1, E.e1, t); mpz_mod(E.e1, E.e1, bp.p);
        EXPECT(!eg_entry_verify(&bp, &E, 7), "an entry outside the subgroup is rejected");
        eg_entry_prove(&bp, &E, 7, 3, 0, NULL);
        /* C2 * g encrypts 1 while e encrypts 0 */
        mpz_mul(E.C2, E.C2, bp.g); mpz_mod(E.C2, E.C2, bp.p);
        EXPECT(!eg_entry_verify(&bp, &E, 7), "a zero entry with a bit ciphertext of 1 is rejected");
        eg_entry_prove(&bp, &E, 7, 3, 1, NULL);
        /* multiply C2 by p-1, an element of order 2 outside the subgroup */
        mpz_sub_ui(t, bp.p, 1);
        mpz_mul(E.C2, E.C2, t); mpz_mod(E.C2, E.C2, bp.p);
        EXPECT(!eg_entry_verify(&bp, &E, 7), "a bit ciphertext outside the subgroup is rejected");
        eg_entry_free(&bp, &E);
    }

    /* A board of 8 entries, 3 of them nonzero. */
    size_t M = 8;
    int slot[8], bits[8];
    for (size_t j = 0; j < M; j++) { slot[j] = 1 + (int)(j * 5) % bp.n; bits[j] = (j % 3 == 0); }
    EGBoard B;
    eg_board_init(&B, &bp, M);
    eg_board(&bp, &B, slot, bits, 4);
    EXPECT(B.c == 3, "board count c = 3");
    EXPECT(eg_verify_board(&bp, &B, 1), "honest board verifies on one thread");
    EXPECT(eg_verify_board(&bp, &B, 4), "honest board verifies on four threads");
    B.c++;
    EXPECT(!eg_verify_board(&bp, &B, 4), "board with count c+1 is rejected");
    B.c--;
    B.root[0] ^= 1;
    EXPECT(!eg_verify_board(&bp, &B, 4), "board with another root is rejected");
    B.root[0] ^= 1;
    B.sig[0] ^= 1;
    EXPECT(!eg_verify_board(&bp, &B, 4), "board with a changed signature is rejected");
    B.sig[0] ^= 1;
    /* entry 2 replaced by a valid entry that is not under the signed root */
    eg_entry_prove(&bp, &B.E[2], 2, slot[2], bits[2], NULL);
    EXPECT(!eg_verify_board(&bp, &B, 4), "a valid entry that is not under the signed root is rejected");
    /* the same, with the root and signature recomputed: the count proof fails */
    eg_board(&bp, &B, slot, bits, 4);
    {
        uint8_t (*leaf)[32] = malloc(M * 32);
        uint8_t msg[BOARD_MSG_BYTES];
        eg_entry_prove(&bp, &B.E[2], 2, slot[2], bits[2], NULL);
        for (size_t j = 0; j < M; j++) eg_entry_leaf(&B.E[j], leaf[j]);
        merkle_root((const uint8_t (*)[32])leaf, M, B.root);
        board_message(msg, M, B.root);
        B.siglen = sig_sign(&bp.sig, msg, sizeof msg, B.sig, bp.sig.sig_bytes);
        EXPECT(!eg_verify_board(&bp, &B, 4), "count proof over a changed bit ciphertext is rejected");
        free(leaf);
    }
    eg_board(&bp, &B, slot, bits, 4);
    EXPECT(eg_verify_board(&bp, &B, 4), "recomputed board verifies");

    mpz_clear(t);
    eg_board_free(&B, &bp);
    eg_board_params_free(&bp);
    ect_params_free(&pp);
    core_clean();
    printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
