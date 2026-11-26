#include "dj.h"
#include "dj_board.h"
#include "board_util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* test_dj_board.c -> accept/reject tests for the Damgard-Jurik board proof.
   Usage: ./test_dj_board [h] [bits]   (defaults: 10 384, a modulus too
   small for security that keeps the n + 3 atoms of each entry fast, with
   prime factors above the challenge bound 2^128) */

static int failures = 0;

#define EXPECT(cond, ...) do {                                   \
    if (cond) { printf("ok    "); }                             \
    else { printf("FAIL  "); failures++; }                      \
    printf(__VA_ARGS__); printf("\n"); fflush(stdout);          \
} while (0)

int main(int argc, char **argv)
{
    int h = argc > 1 ? atoi(argv[1]) : 10;
    int bits = argc > 2 ? atoi(argv[2]) : 384;
    DJParams pp;
    dj_setup(&pp, h, bits);
    printf("Damgard-Jurik board: h=%d |N|=%d n=%d\n", h, bits, pp.n);
    mpz_t t, keep;
    mpz_inits(t, keep, NULL);

    {
        DJEntry E;
        dj_entry_init(&pp, &E);
        int all = 1;
        for (int s = 1; s <= pp.n; s++) {
            dj_entry_prove(&pp, &E, 7, s, 1, NULL);
            all &= dj_entry_verify(&pp, &E, 7);
        }
        EXPECT(all, "a nonzero entry verifies for each of the %d identifiers", pp.n);
        dj_entry_prove(&pp, &E, 7, 0, 0, NULL);
        EXPECT(dj_entry_verify(&pp, &E, 7), "a zero entry verifies");
        EXPECT(!dj_entry_verify(&pp, &E, 8), "an entry verified at another position is rejected");

        dj_entry_prove(&pp, &E, 7, 3, 1, NULL);
        /* e (1+N)^{m_3} encrypts 2 m_3, which no identifier encodes */
        mpz_pow_ui(t, pp.r, 2);
        dj_one_plus_N_pow(&pp, t, t, h);
        mpz_mul(E.e, E.e, t); mpz_mod(E.e, E.e, pp.Np[h + 1]);
        EXPECT(!dj_entry_verify(&pp, &E, 7), "an entry changed to encrypt 2 m_3 is rejected");
        dj_entry_prove(&pp, &E, 7, 3, 1, NULL);
        mpz_add_ui(E.z[3 + 2], E.z[3 + 2], 1);
        EXPECT(!dj_entry_verify(&pp, &E, 7), "a modified response is rejected");
        dj_entry_prove(&pp, &E, 7, 3, 0, NULL);
        /* C (1+N) encrypts 1 while e encrypts 0 */
        mpz_add_ui(t, pp.N, 1);
        mpz_mul(E.C, E.C, t); mpz_mod(E.C, E.C, pp.Np[2]);
        EXPECT(!dj_entry_verify(&pp, &E, 7), "a zero entry with a bit ciphertext of 1 is rejected");
        dj_entry_prove(&pp, &E, 7, 3, 1, NULL);
        mpz_mul(E.C, E.C, pp.p1); mpz_mod(E.C, E.C, pp.Np[2]);
        EXPECT(!dj_entry_verify(&pp, &E, 7), "a bit ciphertext that is not a unit is rejected");
        dj_entry_free(&pp, &E);
    }

    size_t M = 8;
    int slot[8], b[8];
    for (size_t j = 0; j < M; j++) { slot[j] = 1 + (int)(j * 37) % pp.n; b[j] = (j % 3 == 0); }
    DJBoard B;
    dj_board_init(&B, &pp, M);
    dj_board(&pp, &B, slot, b, 4);
    EXPECT(B.c == 3, "board count c = 3");
    EXPECT(dj_verify_board(&pp, &B, 1), "honest board verifies on one thread");
    EXPECT(dj_verify_board(&pp, &B, 4), "honest board verifies on four threads");
    B.c++;
    EXPECT(!dj_verify_board(&pp, &B, 4), "board with count c+1 is rejected");
    B.c--;
    B.root[0] ^= 1;
    EXPECT(!dj_verify_board(&pp, &B, 4), "board with another root is rejected");
    B.root[0] ^= 1;
    B.sig[0] ^= 1;
    EXPECT(!dj_verify_board(&pp, &B, 4), "board with a changed signature is rejected");
    B.sig[0] ^= 1;
    dj_entry_prove(&pp, &B.E[2], 2, slot[2], b[2], NULL);
    EXPECT(!dj_verify_board(&pp, &B, 4), "a valid entry that is not under the signed root is rejected");
    dj_board(&pp, &B, slot, b, 4);
    {
        uint8_t (*leaf)[32] = malloc(M * 32);
        uint8_t msg[BOARD_MSG_BYTES];
        dj_entry_prove(&pp, &B.E[2], 2, slot[2], b[2], NULL);
        for (size_t j = 0; j < M; j++) dj_entry_leaf(&pp, &B.E[j], leaf[j]);
        merkle_root((const uint8_t (*)[32])leaf, M, B.root);
        board_message(msg, M, B.root);
        B.siglen = sig_sign(&pp.sig, msg, sizeof msg, B.sig, pp.sig.sig_bytes);
        EXPECT(!dj_verify_board(&pp, &B, 4), "count proof over a changed bit ciphertext is rejected");
        free(leaf);
    }
    dj_board(&pp, &B, slot, b, 4);
    EXPECT(dj_verify_board(&pp, &B, 4), "recomputed board verifies");

    mpz_clears(t, keep, NULL);
    dj_board_free(&B, &pp);
    dj_free(&pp);
    printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
