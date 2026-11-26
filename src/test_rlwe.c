#include "rlwe.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* test_rlwe.c -> correctness tests for the Ring-LWE instantiation.
   Usage: ./test_rlwe N0 h log2Q primes d [L]
   (defaults: the two-state profile 4096 16 65 2 2 32).  d is the fan-in of
   the correctness-horizon test and L its number of transitions. */

static int failures = 0;

#define EXPECT(cond, ...) do {                                  \
    if (cond) { printf("ok    "); }                             \
    else { printf("FAIL  "); failures++; }                      \
    printf(__VA_ARGS__); printf("\n"); fflush(stdout);          \
} while (0)

static int trace_is(const RLWEParams *pp, const RLWECt *w, const uint8_t *want, int expect_bot)
{
    uint8_t *got = calloc((size_t)pp->n, 1);
    int k = rlwe_trace(pp, w, got);
    int ok = expect_bot ? (k == -1) : (k >= 0 && memcmp(got, want, (size_t)pp->n) == 0);
    free(got);
    return ok;
}

/* Noise bound of the parameter section: F = Be (1 + 4 sqrt(N0 h)),
   B_{j+1} = d (B_j + F), B_0 = F. */
static double log2_bound(int d, int L, double F)
{
    if (d == 1) return log2((L + 1) * F);
    return log2(F) + log2((pow(d, L) * (2 * d - 1) - d) / (d - 1));
}

int main(int argc, char **argv)
{
    int N0 = argc > 1 ? atoi(argv[1]) : 4096;
    int h = argc > 2 ? atoi(argv[2]) : 16;
    int logQ = argc > 3 ? atoi(argv[3]) : 65;
    int L = argc > 4 ? atoi(argv[4]) : 2;
    int d = argc > 5 ? atoi(argv[5]) : 2;
    int steps = argc > 6 ? atoi(argv[6]) : 32;
    const uint64_t t = 65537;

    RLWEParams pp;
    rlwe_setup(&pp, N0, h, logQ, L, t);
    double F = pp.Be * (1 + 4 * sqrt((double)N0 * h));
    double lb = log2_bound(d, steps, F);
    printf("Ring-LWE: N0=%d h=%d t=%llu log2Q=%d (%d primes) d=%d L=%d tag=%zu bytes\n",
           N0, h, (unsigned long long)t, pp.logQ, L, d, steps, rlwe_tag_bytes(&pp));
    EXPECT(log2((double)pp.Q) > 1 + log2((double)t) + log2(pow(2, lb) + 1),
           "Q > 2t(B_L + 1): log2 Q = %.2f, log2 B_L = %.2f", log2((double)pp.Q), lb);

    size_t n = (size_t)h * N0;
    uint64_t *m = calloc(n, sizeof(uint64_t)), *v = calloc(n, sizeof(uint64_t));
    uint64_t *expect = calloc(n, sizeof(uint64_t)), *next = calloc(n, sizeof(uint64_t));
    uint8_t *want = calloc((size_t)pp.n, 1), *sig = malloc(pp.sig.sig_bytes);
    RLWECt w, x, ws[4];
    rlwe_ct_init(&pp, &w);
    rlwe_ct_init(&pp, &x);
    for (int i = 0; i < 4; i++) rlwe_ct_init(&pp, &ws[i]);

    for (size_t i = 0; i < n; i++) m[i] = (uint64_t)rand() % t;
    rlwe_enc(&pp, &w, m);
    double margin = rlwe_dec(&pp, v, &w);
    EXPECT(memcmp(m, v, n * sizeof(uint64_t)) == 0, "Dec(Enc(m)) = m, noise margin %.1f bits", margin);

    size_t sl = rlwe_tag(&pp, &w, sig, 3, 1, h);
    EXPECT(rlwe_verify_tag(&pp, &w, sig, sl), "issuance signature verifies (ML-DSA-44, %zu bytes)", sl);
    w.c0[0] ^= 1;
    EXPECT(!rlwe_verify_tag(&pp, &w, sig, sl), "signature on a modified tag is rejected");

    /* Trace after k hops, below and at the budget. */
    int budgets[2] = { h, h / 2 };
    for (int bi = 0; bi < 2; bi++) {
        int hp = budgets[bi];
        rlwe_tag(&pp, &w, sig, 5, 1, hp);
        memset(want, 0, (size_t)pp.n);
        want[4] = 1;
        int ok_live = 1, ok_exp = 1;
        for (int k = 0; k <= hp + 1; k++) {
            if (k > 0) rlwe_degrade(&pp, &w, &w);
            if (k < hp) ok_live &= trace_is(&pp, &w, want, 0);
            else        ok_exp  &= trace_is(&pp, &w, want, 1);
        }
        EXPECT(ok_live, "h'=%d: id_5 traced after k = 0..%d hops", hp, hp - 1);
        EXPECT(ok_exp, "h'=%d: Trace returns bottom after k = %d and %d hops", hp, hp, hp + 1);
    }

    /* Merge with live, expired and dummy inputs. */
    rlwe_tag(&pp, &ws[0], sig, 7, 1, h);
    for (int k = 0; k < 2; k++) rlwe_degrade(&pp, &ws[0], &ws[0]);
    rlwe_tag(&pp, &ws[1], sig, 9, 1, 2);
    for (int k = 0; k < 2; k++) rlwe_degrade(&pp, &ws[1], &ws[1]);
    rlwe_tag(&pp, &ws[2], sig, 11, 0, h);
    for (int k = 0; k < 3; k++) rlwe_degrade(&pp, &ws[2], &ws[2]);
    rlwe_tag(&pp, &ws[3], sig, pp.n, 1, h);
    RLWECt *in[4] = { &ws[0], &ws[1], &ws[2], &ws[3] };
    rlwe_merge(&pp, &w, in, 4);
    memset(want, 0, (size_t)pp.n);
    want[6] = 1;
    want[pp.n - 1] = 1;
    EXPECT(trace_is(&pp, &w, want, 0), "Merge(live id_7, expired id_9, dummy, live id_%d) traces {id_7, id_%d}", pp.n, pp.n);
    rlwe_merge(&pp, &w, in + 1, 2);
    EXPECT(trace_is(&pp, &w, want, 1), "Merge(expired, dummy) traces to bottom");

    /* Slot counts at and above the multiplicity bound 2^(h-1): block 3 holds
       count 2^(h-1) in slot 1 and count 1 in slot 2. */
    {
        uint64_t *unit = calloc((size_t)N0, sizeof(uint64_t));
        uint64_t *b3 = m + (size_t)3 * N0;
        uint64_t count[3] = { 0, (uint64_t)1 << (h - 1), 1 };
        memset(m, 0, n * sizeof(uint64_t));
        for (int s = 1; s <= 2; s++) {
            rlwe_slot_unit(&pp, unit, s);
            for (int i = 0; i < N0; i++) b3[i] = (b3[i] + count[s] * unit[i]) % t;
        }
        rlwe_enc(&pp, &w, m);
        memset(want, 0, (size_t)pp.n);
        want[1] = 1; want[2] = 1;
        EXPECT(trace_is(&pp, &w, want, 0), "slot count 2^(h-1) next to count 1 traces {id_2, id_3}");
        rlwe_slot_unit(&pp, unit, 1);
        for (int i = 0; i < N0; i++) b3[i] = (b3[i] + unit[i]) % t;
        rlwe_enc(&pp, &w, m);
        EXPECT(trace_is(&pp, &w, want, 1), "slot count 2^(h-1)+1 traces to bottom");
        free(unit);
    }

    /* Correctness horizon: L transitions of fan-in d.  For d = 1 a
       transition is one degradation; for d >= 2 it merges the tag with d-1
       copies of its own degradation, so that the noise terms add up
       aligned, the worst case of the recurrence B_{j+1} = d (B_j + F). */
    rlwe_tag(&pp, &w, sig, 13, 1, h);
    memset(expect, 0, n * sizeof(uint64_t));
    rlwe_slot_unit(&pp, expect, 12);
    double minmargin = 1e9;
    int horizon_ok = 1;
    for (int j = 1; j <= steps; j++) {
        rlwe_degrade(&pp, &x, &w);
        memset(next, 0, n * sizeof(uint64_t));
        if (d == 1) {
            rlwe_ct_copy(&pp, &w, &x);
            for (int k = 1; k < h; k++)
                memcpy(next + (size_t)k * N0, expect + (size_t)(k - 1) * N0, (size_t)N0 * sizeof(uint64_t));
        } else {
            RLWECt *mi[5] = { &w, &x, &x, &x, &x };
            rlwe_merge(&pp, &w, mi, d);
            for (int k = 0; k < h; k++)
                for (int i = 0; i < N0; i++) {
                    uint64_t a = expect[(size_t)k * N0 + i];
                    uint64_t b = k ? expect[(size_t)(k - 1) * N0 + i] : 0;
                    next[(size_t)k * N0 + i] = (a + (uint64_t)(d - 1) * b) % t;
                }
        }
        memcpy(expect, next, n * sizeof(uint64_t));
        double mg = rlwe_dec(&pp, v, &w);
        if (mg < minmargin) minmargin = mg;
        if (memcmp(v, expect, n * sizeof(uint64_t)) != 0) horizon_ok = 0;
    }
    EXPECT(horizon_ok, "%d transitions of fan-in %d decrypt correctly, smallest noise margin %.1f bits", steps, d, minmargin);
    printf("horizon,%d,%d,%d,%d,%d,%d,%.2f,%.2f,%.2f\n", N0, h, pp.logQ, d, steps,
           horizon_ok, minmargin, lb, log2((double)pp.Q));

    rlwe_ct_free(&w); rlwe_ct_free(&x);
    for (int i = 0; i < 4; i++) rlwe_ct_free(&ws[i]);
    free(m); free(v); free(expect); free(next); free(want); free(sig);
    rlwe_free(&pp);
    printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
