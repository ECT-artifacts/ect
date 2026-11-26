#include "dj.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* test_dj.c -> correctness tests for the Damgard-Jurik instantiation; the
   board proof is tested in test_dj_board.c.  Usage: ./test_dj [h] [bits]
   (defaults: 10 3072) */

static int failures = 0;

#define EXPECT(cond, ...) do {                                   \
    if (cond) { printf("ok    "); }                             \
    else { printf("FAIL  "); failures++; }                      \
    printf(__VA_ARGS__); printf("\n"); fflush(stdout);          \
} while (0)

/* Does Trace return exactly the identifiers in want (want[i] = 1 for id_{i+1})? */
static int trace_is(const DJParams *pp, const mpz_t w, const uint8_t *want, int expect_bot)
{
    uint8_t *got = calloc((size_t)pp->n, 1);
    int k = dj_trace(pp, w, got);
    int ok;
    if (expect_bot) ok = (k == -1);
    else ok = (k >= 0) && memcmp(got, want, (size_t)pp->n) == 0;
    free(got);
    return ok;
}

int main(int argc, char **argv)
{
    int h = argc > 1 ? atoi(argv[1]) : 10;
    int bits = argc > 2 ? atoi(argv[2]) : 3072;
    DJParams pp;
    dj_setup(&pp, h, bits);
    printf("Damgard-Jurik: h=%d |N|=%d n=%d tag=%zu bytes\n",
           h, bits, pp.n, dj_elem_len(&pp, h));
    uint8_t *want = calloc((size_t)pp.n, 1);
    uint8_t *sig = malloc(pp.sig.sig_bytes);
    mpz_t w, w2, xi, ws[4];
    mpz_inits(w, w2, xi, NULL);
    for (int i = 0; i < 4; i++) mpz_init(ws[i]);

    /* Decryption of a random plaintext. */
    {
        mpz_t m, v;
        mpz_inits(m, v, NULL);
        mpz_set_ui(m, 123456789);
        mpz_mul(m, m, pp.Np[h - 1]);
        mpz_add_ui(m, m, 987654321);
        dj_enc(&pp, w, m, h, NULL);
        dj_dec(&pp, v, w, h);
        EXPECT(mpz_cmp(v, m) == 0, "Dec(Enc(m)) = m");
        mpz_clears(m, v, NULL);
    }

    /* The authority's computation of x^{N^s} equals the public one. */
    {
        mpz_t x, y1, y2;
        mpz_inits(x, y1, y2, NULL);
        int same = 1;
        for (int k = 0; k < 3; k++) {
            dj_rand_unit(&pp, x);
            dj_pow_Ns(&pp, y1, x, h);
            dj_pow_Ns_sk(&pp, y2, x, h);
            same &= mpz_cmp(y1, y2) == 0;
            dj_pow_Ns(&pp, y1, x, 1);
            dj_pow_Ns_sk(&pp, y2, x, 1);
            same &= mpz_cmp(y1, y2) == 0;
        }
        EXPECT(same, "x^{N^s} with the factorization equals x^{N^s} without it, s = h and s = 1");
        mpz_clears(x, y1, y2, NULL);
    }

    /* Signature on an issued tag. */
    size_t sl = dj_tag(&pp, w, sig, 3, 1, h);
    EXPECT(dj_verify_tag(&pp, w, sig, sl), "issuance signature verifies");
    sig[0] ^= 1;
    EXPECT(!dj_verify_tag(&pp, w, sig, sl), "modified signature is rejected");
    sig[0] ^= 1;

    /* Trace after k hops, below and at the budget, for h' = h and h' = h/2. */
    int budgets[2] = { h, h / 2 > 0 ? h / 2 : 1 };
    for (int bi = 0; bi < 2; bi++) {
        int hp = budgets[bi];
        dj_tag(&pp, w, sig, 5, 1, hp);
        memset(want, 0, (size_t)pp.n);
        want[4] = 1;
        int ok_live = 1, ok_exp = 1;
        for (int k = 0; k <= hp + 1; k++) {
            if (k > 0) { dj_degrade(&pp, w2, xi, w); mpz_set(w, w2); }
            if (k < hp) ok_live &= trace_is(&pp, w, want, 0);
            else        ok_exp  &= trace_is(&pp, w, want, 1);
        }
        EXPECT(ok_live, "h'=%d: id_5 traced after k = 0..%d hops", hp, hp - 1);
        EXPECT(ok_exp, "h'=%d: Trace returns bottom after k = %d and %d hops", hp, hp, hp + 1);
    }

    /* Merge with live, expired and dummy inputs. */
    {
        int hp_exp = 2;
        dj_tag(&pp, ws[0], sig, 7, 1, h);          /* live, degraded twice */
        for (int k = 0; k < 2; k++) { dj_degrade(&pp, w2, xi, ws[0]); mpz_set(ws[0], w2); }
        dj_tag(&pp, ws[1], sig, 9, 1, hp_exp);     /* expired after 2 hops */
        for (int k = 0; k < hp_exp; k++) { dj_degrade(&pp, w2, xi, ws[1]); mpz_set(ws[1], w2); }
        dj_tag(&pp, ws[2], sig, 11, 0, h);         /* dummy, degraded 3 times */
        for (int k = 0; k < 3; k++) { dj_degrade(&pp, w2, xi, ws[2]); mpz_set(ws[2], w2); }
        dj_tag(&pp, ws[3], sig, pp.n, 1, h);       /* live, last slot, fresh */
        dj_merge(&pp, w, xi, ws, 4);
        memset(want, 0, (size_t)pp.n);
        want[6] = 1;
        want[pp.n - 1] = 1;
        EXPECT(trace_is(&pp, w, want, 0), "Merge(live id_7, expired id_9, dummy, live id_%d) traces {id_7, id_%d}", pp.n, pp.n);
        dj_merge(&pp, w, xi, ws + 1, 2);
        EXPECT(trace_is(&pp, w, want, 1), "Merge(expired, dummy) traces to bottom");
    }

    /* Multiplicity: 2^(h-1) copies of one contribution at one depth. */
    if (h <= 12) {
        int copies = 1 << (h - 1);
        dj_tag(&pp, w, sig, 2, 1, h);
        mpz_t acc;
        mpz_init_set(acc, w);
        for (int c = 1; c < copies; c++) {
            mpz_mul(acc, acc, w);
            mpz_mod(acc, acc, pp.Np[h + 1]);
        }
        dj_tag(&pp, ws[0], sig, 3, 1, h);
        mpz_set(ws[1], acc);
        dj_merge(&pp, w, xi, ws, 2);
        memset(want, 0, (size_t)pp.n);
        want[1] = 1; want[2] = 1;
        EXPECT(trace_is(&pp, w, want, 0), "%d copies of id_2 next to id_3 trace {id_2, id_3}", copies);
        mpz_clear(acc);
    }

    /* A depth digit at or above r^n makes Trace return bottom. */
    dj_enc(&pp, w, pp.rn, h, NULL);
    EXPECT(trace_is(&pp, w, want, 1), "digit r^n at depth 0 traces to bottom");

    mpz_clears(w, w2, xi, NULL);
    for (int i = 0; i < 4; i++) mpz_clear(ws[i]);
    free(want); free(sig);
    dj_free(&pp);
    printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
