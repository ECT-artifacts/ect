#include <stdio.h>
#include <stdlib.h>
#include "ect.h"
#include "nizk.h"
#include "schnorr.h"
#include "block.h"
#include "sim.h"
#include "scenario.h"

/* main.c -> Test suite and simulation entry point for the ECT prototype. */

static void sep(const char *label) { printf("\n=== %s ===\n", label); }

/* Convenience wrapper for single-input Degrade (account-based model).
    Builds C = {current ciphertext} and calls ect_degrade with m=1, idx=0. */
static void degrade_single(ECTTag *w, const ECTParams *pp)
{
    bn_t C1, C2;
    bn_null(C1); bn_null(C2);
    RLC_TRY {
        bn_new(C1); bn_new(C2);
        bn_copy(C1, w->c1);
        bn_copy(C2, w->c2);
        ect_degrade(w, pp, &C1, &C2, 1, 0);
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(C1); bn_free(C2);
    }
}



int main(void)
{
    core_init();

    ECTParams pp;
    ect_params_new(&pp);
    ect_setup(&pp);
    printf("Setup done. H=%d R_BASE=%d N=%d\n", pp.h, pp.r, pp.n);

    /*  Tests 1-3: Trace correctness (Theorem 6.3) and expiry (Theorem 6.5).*/
       

    sep("Test 1: fresh tag, user 1 -> recovered at depth 0");
    ect_clear_tagged(&pp);
    /* Verifies: Tag(id=1, b=1) followed by Trace returns {user 1, depth 0} */
    { ECTTag w; ect_tag_new(&w); ect_tag(&w, &pp, 1, 1);
      ect_trace(&pp, &w); ect_tag_free(&w); }

    sep("Test 2: tag degraded h-1 times -> still recoverable at depth h-1");
    ect_clear_tagged(&pp);
    /* Verifies: after h-1 Degrade steps, user is recovered at depth h-1. Confirms the tag has not expired prematurely. */
    { ECTTag w; ect_tag_new(&w); ect_tag(&w, &pp, 2, 1);
      for (int i = 0; i < pp.h-1; i++) degrade_single(&w, &pp);
      ect_trace(&pp, &w); ect_tag_free(&w); }

    sep("Test 3: tag degraded h times -> expires exactly at hop budget");
    ect_clear_tagged(&pp);
    /* Verifies that after exactly h Degrade steps, the plaintext collapses to 0 mod q^H and Trace returns ⊥. Unconditional expiry. */
    { ECTTag w; ect_tag_new(&w); ect_tag(&w, &pp, 2, 1);
      for (int i = 0; i < pp.h; i++) degrade_single(&w, &pp);
      ect_trace(&pp, &w); ect_tag_free(&w); }

    /* Test 4: Indistinguishability (Definition 4.3). */
       

    sep("Test 4: dummy tag (b=0) -> trace returns expired");
    ect_clear_tagged(&pp);
    /* Verifies: Tag(id=3, b=0) encrypts 0.  Trace sees the same result  as an expired tag */
    { ECTTag w; ect_tag_new(&w); ect_tag(&w, &pp, 3, 0);
      ect_trace(&pp, &w); ect_tag_free(&w); }

    /* Tests 5-6: Merge correctness (Theorem 6.4).*/
       

    sep("Test 5: merge user 1 + user 3 at depth 0 -> both recovered");
    ect_clear_tagged(&pp);
    /* Verifies: Merge of two depth-0 tags recovers both users at depth 0.
        Tests the homomorphic multiplication: Enc(m1) * Enc(m2) = Enc(m1+m2). */
    { ECTTag w1, w3, merged;
      ect_tag_new(&w1); ect_tag_new(&w3); ect_tag_new(&merged);
      ect_tag(&w1, &pp, 1, 1); ect_tag(&w3, &pp, 3, 1);
      ECTTag *arr[2] = { &w1, &w3 };
      ect_merge(&merged, arr, 2, &pp);
      ect_trace(&pp, &merged);
      ect_tag_free(&w1); ect_tag_free(&w3); ect_tag_free(&merged); }

    sep("Test 6: merge user 1 (depth 2) + user 3 (depth 0) -> depth-aware");
    ect_clear_tagged(&pp);
    /* Verifies: Merge of tags at different depths recovers each user at their correct depth.  Tests the base-q digit separation: user 1
        appears in digit d[2], user 3 in digit d[0], with no collision. */
    { ECTTag w1, w3, merged;
      ect_tag_new(&w1); ect_tag_new(&w3); ect_tag_new(&merged);
      ect_tag(&w1, &pp, 1, 1);
      degrade_single(&w1, &pp); degrade_single(&w1, &pp);
      ect_tag(&w3, &pp, 3, 1);
      ECTTag *arr[2] = { &w1, &w3 };
      ect_merge(&merged, arr, 2, &pp);
      ect_trace(&pp, &merged);
      ect_tag_free(&w1); ect_tag_free(&w3); ect_tag_free(&merged); }

    /* Test 7: NIZK proof for R_exp */
       

    sep("Test 7: NIZK proof on one Degrade step");
    ect_clear_tagged(&pp);
    /* Verifies: nizk_prove generates a valid proof for (c_before -> c_after), and nizk_verify correctly detects tampering with A1. */
    {
        ECTTag w;
        ect_tag_new(&w);
        ect_tag(&w, &pp, 2, 1);

        bn_t c1_before, c2_before;
        bn_null(c1_before); bn_null(c2_before);
        bn_new(c1_before);  bn_new(c2_before);
        bn_copy(c1_before, w.c1);
        bn_copy(c2_before, w.c2);
        degrade_single(&w, &pp);

        if (w.n_proofs == 1) {
            /* C = {c_before} -> single input set used during degrade_single */
            bn_t C1[1], C2[1];
            bn_null(C1[0]); bn_null(C2[0]);
            bn_new(C1[0]);  bn_new(C2[0]);
            bn_copy(C1[0], c1_before);
            bn_copy(C2[0], c2_before);

            int ok = nizk_verify(&w.proofs[0], pp.pk,
                                 C1, C2, 1,
                                 w.c1, w.c2,
                                 pp.g, pp.q, pp.qh, pp.p);
            printf("  NIZK verify:              %s\n", ok ? "PASS" : "FAIL");

            /* Tamper with A1[0] -> verification must fail */
            bn_add_dig(w.proofs[0].A1[0], w.proofs[0].A1[0], 1);
            int ok2 = nizk_verify(&w.proofs[0], pp.pk,
                                  C1, C2, 1,
                                  w.c1, w.c2,
                                  pp.g, pp.q, pp.qh, pp.p);
            printf("  Tamper test (expect FAIL): %s\n", ok2 ? "PASS" : "FAIL");
            bn_free(C1[0]); bn_free(C2[0]);
        }
        bn_free(c1_before); bn_free(c2_before);
        ect_tag_free(&w);
    }

    /* Test 8: Block serialization */
       

    sep("Test 8: block serialize / deserialize round-trip");
    ect_clear_tagged(&pp);
    /* Verifies: wm_pack -> block -> wm_unpack preserves (c1, c2) and the NIZK proof fields.  
        After unpacking, nizk_verify still passes and the ledger hash chain is consistent across two blocks. */
    {
        ECTTag w;
        ect_tag_new(&w);
        ect_tag(&w, &pp, 3, 1);

        bn_t c1_before, c2_before;
        bn_null(c1_before); bn_null(c2_before);
        bn_new(c1_before);  bn_new(c2_before);
        bn_copy(c1_before, w.c1);
        bn_copy(c2_before, w.c2);
        degrade_single(&w, &pp);

        Watermark wm;
        wm_pack(&wm, w.c1, w.c2,
                w.n_proofs > 0 ? &w.proofs[0] : NULL,
                (uint8_t)w.hops,
                w.c1_0, w.c2_0);

        Block b;
        uint8_t zero_hash[HASH_BYTES] = {0};
        block_set_header(&b, zero_hash, 0, 1, 3, 100);
        block_set_watermark(&b, &wm);
        block_print(&b);

        Watermark wm2;
        block_get_watermark(&wm2, &b);
        bn_t rc1, rc2;
        /* Heap-allocate: NIZKProof is large with MAX_INPUTS=128 */
        NIZKProof *rpi = (NIZKProof *)malloc(sizeof(NIZKProof));
        uint8_t hops;
        bn_null(rc1); bn_null(rc2);
        bn_new(rc1);  bn_new(rc2);
        nizk_proof_new(rpi);
        wm_unpack(rc1, rc2, rpi, &hops, NULL, NULL, &wm2);

        int ok;
        {
            bn_t C1[1], C2[1];
            bn_null(C1[0]); bn_null(C2[0]);
            bn_new(C1[0]);  bn_new(C2[0]);
            bn_copy(C1[0], c1_before);
            bn_copy(C2[0], c2_before);
            ok = nizk_verify(rpi, pp.pk,
                             C1, C2, 1,
                             rc1, rc2,
                             pp.g, pp.q, pp.qh, pp.p);
            bn_free(C1[0]); bn_free(C2[0]);
        }
        printf("  Round-trip NIZK verify: %s  hops=%u\n",
               ok ? "PASS" : "FAIL", hops);

        Ledger *ledger = ledger_new();
        ledger_append(ledger, &b);
        Block b2;
        Watermark wm3;
        wm_pack(&wm3, rc1, rc2, NULL, hops, NULL, NULL);
        block_set_header(&b2, zero_hash, 0, 3, 7, 100);
        block_set_watermark(&b2, &wm3);
        ledger_append(ledger, &b2);
        ledger_print(ledger);
        ledger_free(ledger);

        bn_free(c1_before); bn_free(c2_before);
        bn_free(rc1); bn_free(rc2);
        nizk_proof_free(rpi);
        free(rpi);
        ect_tag_free(&w);
    }

    /* Test 9: Authority signature (Schnorr).*/
       

    sep("Test 9: authority signature on (c1_0, c2_0)");
    ect_clear_tagged(&pp);
    /* Verifies:
     - Fresh tag: ect_verify_tag_sig returns PASS (sig covers c1_0, c2_0)
     - After Degrade: sig is CLEARED or Degrade output is (c', null, {proof}) -> no signature forwarded). FAIL is correct.
     - Tamper c1_0 on fresh tag: FAIL (sig no longer matches) */
    {
        ECTTag w;
        ect_tag_new(&w);
        ect_tag(&w, &pp, 1, 1);

        int ok = ect_verify_tag_sig(&w, &pp);
        printf("  Fresh tag sig verify:                    %s\n",
               ok ? "PASS" : "FAIL");

        /* Degrade clears has_sig per paper -> sig is not forwarded */
        degrade_single(&w, &pp);
        int ok2 = ect_verify_tag_sig(&w, &pp);
        printf("  After Degrade, sig cleared (expect FAIL): %s\n",
               ok2 ? "PASS" : "FAIL");

        /* Verify on a fresh tag with tampered c1_0 -> must fail */
        ECTTag w2;
        ect_tag_new(&w2);
        ect_tag(&w2, &pp, 1, 1);
        bn_add_dig(w2.c1_0, w2.c1_0, 1);
        int ok3 = ect_verify_tag_sig(&w2, &pp);
        printf("  Tamper c1_0 on fresh tag (expect FAIL):  %s\n",
               ok3 ? "PASS" : "FAIL");

        ect_tag_free(&w);
        ect_tag_free(&w2);
    }

    /* Tests 10-13: full-capacity (N) scale tests. */
       
    sep("Test 10: N=28, tag 6 users, all recovered at depth 0");
    ect_clear_tagged(&pp);
    /* Verifies the encoding capacity claim: 6 simultaneous traces out of
      27 users, all at depth 0.  Users spread across the range to confirm
      no identifier collisions (users 1, 6, 11, 16, 22, 27). */
    {
        int traced[] = {1, 6, 11, 16, 22, 27};
        int n_traced = 6;
        ECTTag *tags[6];
        for (int i = 0; i < n_traced; i++) {
            tags[i] = (ECTTag *)malloc(sizeof(ECTTag));
            ect_tag_new(tags[i]);
            ect_tag(tags[i], &pp, traced[i], 1);
        }
        ECTTag merged;
        ect_tag_new(&merged);
        ect_merge(&merged, tags, n_traced, &pp);
        printf("  Traced users: ");
        for (int i = 0; i < n_traced; i++)
            printf("%d%s", traced[i], i < n_traced-1 ? ", " : "\n");
        ect_trace(&pp, &merged);

        /* Verify all authority signatures survive the merge */
        int all_sigs_ok = 1;
        for (int i = 0; i < n_traced; i++)
            if (!ect_verify_tag_sig(tags[i], &pp)) all_sigs_ok = 0;
        printf("  All authority sigs valid: %s\n", all_sigs_ok ? "YES" : "NO");

        ect_tag_free(&merged);
        for (int i = 0; i < n_traced; i++) {
            ect_tag_free(tags[i]);
            free(tags[i]);
        }
    }

    sep("Test 11: N=28, 6 users at mixed depths -> depth-aware recovery");
    ect_clear_tagged(&pp);
    /* Verifies merge correctness at scale: each of the 6
      users is degraded a different number of times before merging.
      After merge, each user must appear at their exact depth. */
    {
        int users[]  = {  2,  7, 12, 17, 22, 26 };
        int depths[] = {  0,  2,  4,  1,  3,  0 };
        int n = 6;
        ECTTag *tags[6];
        for (int i = 0; i < n; i++) {
            tags[i] = (ECTTag *)malloc(sizeof(ECTTag));
            ect_tag_new(tags[i]);
            ect_tag(tags[i], &pp, users[i], 1);
            for (int d = 0; d < depths[i]; d++)
                degrade_single(tags[i], &pp);
        }
        ECTTag merged;
        ect_tag_new(&merged);
        ect_merge(&merged, tags, n, &pp);
        printf("  Expected: ");
        for (int i = 0; i < n; i++)
            printf("user %d@depth %d%s",
                   users[i], depths[i], i < n-1 ? ", " : "\n");
        ect_trace(&pp, &merged);
        ect_tag_free(&merged);
        for (int i = 0; i < n; i++) { ect_tag_free(tags[i]); free(tags[i]); }
    }

    sep("Test 12: N=28, tag degraded h times -> unconditional expiry");
    ect_clear_tagged(&pp);
    /* Verifies at full capacity: expiry is independent of N. */
    { ECTTag w; ect_tag_new(&w); ect_tag(&w, &pp, 7, 1);
      for (int i = 0; i < pp.h; i++) degrade_single(&w, &pp);
      ect_trace(&pp, &w); ect_tag_free(&w); }

    sep("Test 13: N=28, dummy tags in merge -> only traceable user recovered");
    ect_clear_tagged(&pp);
    /* Verifies: when a traceable tag is merged with dummy tags (b=0),  only the traceable user appears in the Trace output. */
    {
        ECTTag w8, w9, w10, w11, merged;
        ect_tag_new(&w8);  ect_tag_new(&w9);
        ect_tag_new(&w10); ect_tag_new(&w11);
        ect_tag_new(&merged);
        ect_tag(&w8,  &pp,  8, 1);
        ect_tag(&w9,  &pp,  9, 0);
        ect_tag(&w10, &pp, 10, 0);
        ect_tag(&w11, &pp, 11, 0);
        ECTTag *arr[4] = { &w8, &w9, &w10, &w11 };
        ect_merge(&merged, arr, 4, &pp);
        printf("  Expected: only user 8 recovered\n");
        ect_trace(&pp, &merged);
        ect_tag_free(&w8);  ect_tag_free(&w9);
        ect_tag_free(&w10); ect_tag_free(&w11);
        ect_tag_free(&merged);
    }

    /* Test 14: complete per-epoch tag board. */

    sep("Test 14: full tag board, one ciphertext per enrolled user");
    ect_clear_tagged(&pp);
    {
        const int traced[] = {5, 11, 16, 21, 24, 27};
        const int n_traced = (int)(sizeof(traced) / sizeof(traced[0]));
        ECTTag *board = (ECTTag *)malloc((size_t)pp.n * sizeof(ECTTag));
        int signatures_ok = 1;
        int recovered_real = 0;
        int recovered_dummy = 0;

        if (board == NULL) {
            fprintf(stderr, "  Failed to allocate full tag board\n");
        } else {
            /* Publish one indistinguishable ciphertext for every identity. */
            for (int id = 1; id <= pp.n; id++) {
                int is_traced = 0;
                for (int j = 0; j < n_traced; j++) {
                    if (traced[j] == id) {
                        is_traced = 1;
                        break;
                    }
                }
                ect_tag_new(&board[id - 1]);
                ect_tag(&board[id - 1], &pp, id, is_traced);
                if (!ect_verify_tag_sig(&board[id - 1], &pp))
                    signatures_ok = 0;
            }

            /* The authority must recover exactly the nonzero board entries. */
            for (int id = 1; id <= pp.n; id++) {
                int found = ect_trace_silent(&pp, &board[id - 1]);
                int expected = 0;
                for (int j = 0; j < n_traced; j++)
                    if (traced[j] == id) expected = 1;
                if (found > 0 && expected) recovered_real++;
                if (found > 0 && !expected) recovered_dummy++;
            }

            printf("  Board entries:              %d\n", pp.n);
            printf("  Real / dummy entries:       %d / %d\n",
                   n_traced, pp.n - n_traced);
            printf("  All authority signatures:   %s\n",
                   signatures_ok ? "PASS" : "FAIL");
            printf("  Exact traced-set recovery:  %s (%d/%d real, %d dummy)\n",
                   recovered_real == n_traced && recovered_dummy == 0
                       ? "PASS" : "FAIL",
                   recovered_real, n_traced, recovered_dummy);

            for (int i = 0; i < pp.n; i++) ect_tag_free(&board[i]);
            free(board);
        }
    }

    /* Full simulation: 10 validators, N identifier slots, multi-hops.*/
       
    run_ect_scenario(&pp); 

    /* Benchmarks: timing measurements */
       
    // run_benchmarks(&pp);

    ect_params_free(&pp);
    core_clean();
    return 0;
}
