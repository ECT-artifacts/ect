/* test_mitm.c -> correctness tests for meet-in-the-middle Trace recovery.
 *
 *   make test_mitm && ./test_mitm
 *
 * Covers the cases the benchmark harness never exercises: repeated identifiers
 * at one depth (multiplicity > 1), identifiers spread across several depths,
 * and expiry at the hop budget.
 */
#include "ect.h"
#include "elgamal.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static ECTParams pp;
static int failures = 0;

static void degrade_once(ECTTag *w)
{
    bn_t C1, C2;
    bn_null(C1); bn_null(C2);
    RLC_TRY {
        bn_new(C1); bn_new(C2);
        bn_copy(C1, w->c1); bn_copy(C2, w->c2);
        ect_degrade(w, &pp, &C1, &C2, 1, 0);
    } RLC_CATCH_ANY {
    } RLC_FINALLY {
        bn_free(C1); bn_free(C2);
    }
}

static void check(const char *name, int got, int want)
{
    int ok = (got == want);
    if (!ok) failures++;
    printf("  [%s] %-52s got %d, want %d\n",
           ok ? "PASS" : "FAIL", name, got, want);
}

/* Tag `ids[i]` (1-based), degrade it `depth[i]` times, merge all, then trace.
   Returns the number of (depth, identifier) pairs Trace recovered. */
static int tag_merge_trace(const int *ids, const int *depth, int k)
{
    ECTTag *tags = (ECTTag *)malloc((size_t)k * sizeof(ECTTag));
    ECTTag **ptr = (ECTTag **)malloc((size_t)k * sizeof(ECTTag *));
    ECTTag merged;

    for (int i = 0; i < k; i++) {
        ect_tag_new(&tags[i]);
        ect_tag(&tags[i], &pp, ids[i], 1);
        for (int d = 0; d < depth[i]; d++) degrade_once(&tags[i]);
        ptr[i] = &tags[i];
    }
    ect_tag_new(&merged);
    ect_merge(&merged, ptr, k, &pp);

    int found = ect_trace(&pp, &merged);

    ect_tag_free(&merged);
    for (int i = 0; i < k; i++) ect_tag_free(&tags[i]);
    free(tags); free(ptr);
    return found;
}

int main(void)
{
    core_init();
    ect_params_new(&pp);
    ect_setup(&pp);
    printf("Setup: H=%d R_BASE=%d N=%d MULT_MAX=%d\n\n",
           pp.h, pp.r, pp.n, MULT_MAX);

    /* T1: distinct identifiers, all fresh.  One (depth 0, id) pair each. */
    {
        int k = 6;
        int ids[6]   = {1, 2, 3, 5, 8, 13};
        int depth[6] = {0, 0, 0, 0, 0, 0};
        ect_clear_tagged(&pp);
        printf("T1 distinct identifiers, depth 0:\n");
        check("6 distinct ids recovered", tag_merge_trace(ids, depth, k), 6);
    }

    /* T2: the same identifier arriving twice at one depth.  This is what
       R_BASE > 2^(H-1) exists for, and what subset enumeration could not
       decode: the digit is 2*r^4, which is not a sum of distinct r^i. */
    {
        int k = 2;
        int ids[2]   = {5, 5};
        int depth[2] = {0, 0};
        ect_clear_tagged(&pp);
        printf("\nT2 multiplicity 2 at one depth (expect 'user 5 (x2)'):\n");
        check("repeated id recovered as one pair", tag_merge_trace(ids, depth, k), 1);
    }

    /* T3: multiplicity 3, still within MULT_MAX. */
    {
        int k = 3;
        int ids[3]   = {7, 7, 7};
        int depth[3] = {0, 0, 0};
        ect_clear_tagged(&pp);
        printf("\nT3 multiplicity 3 at one depth (expect 'user 7 (x3)'):\n");
        check("triple id recovered as one pair", tag_merge_trace(ids, depth, k), 1);
    }

    /* T4: identifiers spread over distinct depths. */
    {
        int k = 4;
        int ids[4]   = {2, 4, 6, 9};
        int depth[4] = {0, 1, 2, 3};
        ect_clear_tagged(&pp);
        printf("\nT4 four ids at depths 0..3:\n");
        check("4 pairs across 4 depths", tag_merge_trace(ids, depth, k), 4);
    }

    /* T5: same identifier at two different depths -> two separate pairs. */
    {
        int k = 2;
        int ids[2]   = {11, 11};
        int depth[2] = {0, 2};
        ect_clear_tagged(&pp);
        printf("\nT5 one id at two depths:\n");
        check("same id at depths 0 and 2", tag_merge_trace(ids, depth, k), 2);
    }

    /* T6: expiry.  After H degradations the plaintext is 0. */
    {
        int k = 1;
        int ids[1]   = {3};
        int depth[1] = {H};
        ect_clear_tagged(&pp);
        printf("\nT6 expiry at the hop budget:\n");
        check("tag expired, nothing recovered", tag_merge_trace(ids, depth, k), 0);
    }

    /* T7: the highest identifier the capacity allows. */
    {
        int k = 2;
        int ids[2]   = {1, N};
        int depth[2] = {0, 0};
        ect_clear_tagged(&pp);
        printf("\nT7 boundary identifiers 1 and N=%d:\n", N);
        check("first and last id recovered", tag_merge_trace(ids, depth, k), 2);
    }

    printf("\n%s (%d failure%s)\n", failures ? "FAILURES" : "ALL PASS",
           failures, failures == 1 ? "" : "s");
    ect_params_free(&pp);
    core_clean();
    return failures ? 1 : 0;
}
