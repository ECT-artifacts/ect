#include "bench_util.h"
#include "dj.h"
#include "dj_board.h"
#include <string.h>

/* bench_dj.c -> Damgard-Jurik benchmarks.
   Usage: ./bench_dj h runs trace_runs board_M threads [bits]
   board_M = 0 skips the board proof.  One entry is proved and verified
   board_runs = min(runs, 3) times on one thread, and the full board of
   board_M entries once on the given number of threads.  For h > 20 the untimed warm-up is
   skipped, since one operation takes minutes.  Output: CSV rows of
   bench_util.h. */

typedef struct {
    DJParams *pp;
    mpz_t w, w2, xi, ws[2], wn;
    uint8_t *sig;
    size_t siglen;
    uint8_t *in_S;
    DJBoard B;
    DJEntry E;
    int M, *slot, *b, threads;
    mpz_t tau, prod, t;
} Ctx;

static void op_tag(void *o)    { Ctx *c = o; c->siglen = dj_tag(c->pp, c->w, c->sig, 1, 1, c->pp->h); }
static void op_verify(void *o) { Ctx *c = o; if (!dj_verify_tag(c->pp, c->w, c->sig, c->siglen)) abort(); }
static void op_degrade(void *o){ Ctx *c = o; dj_degrade(c->pp, c->w2, c->xi, c->w); }
static void op_merge(void *o)  { Ctx *c = o; dj_merge(c->pp, c->w2, c->xi, c->ws, 2); }
static void op_rerand(void *o) { Ctx *c = o; dj_rerandomizer(c->pp, c->t, c->xi); }
static void op_degrade_pre(void *o) { Ctx *c = o; dj_degrade_pre(c->pp, c->w2, c->t, c->w); }
static void op_merge_pre(void *o)   { Ctx *c = o; dj_merge_pre(c->pp, c->w2, c->t, c->ws, 2); }
static void op_trace1(void *o) { Ctx *c = o; if (dj_trace(c->pp, c->w, c->in_S) != 1) abort(); }
static void op_tracen(void *o) { Ctx *c = o; if (dj_trace(c->pp, c->wn, c->in_S) != c->pp->n) abort(); }
static void op_entry(void *o)  { Ctx *c = o; dj_entry_prove(c->pp, &c->E, 0, 5, 1, c->tau); }
static void op_ventry(void *o) { Ctx *c = o; if (!dj_entry_verify(c->pp, &c->E, 0)) abort(); }
static void op_board(void *o)  { Ctx *c = o; dj_board(c->pp, &c->B, c->slot, c->b, c->threads); }
static void op_vboard(void *o) { Ctx *c = o; if (!dj_verify_board(c->pp, &c->B, c->threads)) abort(); }

int main(int argc, char **argv)
{
    int h = argc > 1 ? atoi(argv[1]) : 10;
    int runs = argc > 2 ? atoi(argv[2]) : 10;
    int trace_runs = argc > 3 ? atoi(argv[3]) : 3;
    int M = argc > 4 ? atoi(argv[4]) : 0;
    int threads = argc > 5 ? atoi(argv[5]) : 1;
    int bits = argc > 6 ? atoi(argv[6]) : 3072;

    if (h > 20) bench_warmup = 0;
    DJParams pp;
    dj_setup(&pp, h, bits);
    char label[64];
    snprintf(label, sizeof label, "dj_h%d_N%d_n%d", h, bits, pp.n);
    fprintf(stderr, "%s: tag %zu bytes\n", label, dj_elem_len(&pp, h));

    Ctx c;
    memset(&c, 0, sizeof c);
    c.pp = &pp;
    mpz_inits(c.w, c.w2, c.xi, c.ws[0], c.ws[1], c.wn, c.tau, c.prod, c.t, NULL);
    c.sig = malloc(pp.sig.sig_bytes);
    c.in_S = malloc((size_t)pp.n);
    size_t tagb = dj_elem_len(&pp, h);

    bench_measure(label, "tag", runs, tagb + pp.sig.sig_bytes, op_tag, NULL, &c);
    bench_measure(label, "verify_signature", runs, tagb + pp.sig.sig_bytes, op_verify, NULL, &c);
    bench_measure(label, "degrade", runs, tagb, op_degrade, NULL, &c);
    dj_tag(&pp, c.ws[0], c.sig, 2, 1, h);
    c.siglen = dj_tag(&pp, c.ws[1], c.sig, 3, 1, h);
    bench_measure(label, "merge_k2", runs, tagb, op_merge, NULL, &c);
    /* The split of degrade and merge: the rerandomizer, computed in advance,
       and the work that depends on the tags, with that rerandomizer. */
    bench_measure(label, "rerandomizer", runs, tagb, op_rerand, NULL, &c);
    bench_measure(label, "degrade_online", runs, tagb, op_degrade_pre, NULL, &c);
    bench_measure(label, "merge_k2_online", runs, tagb, op_merge_pre, NULL, &c);
    c.siglen = dj_tag(&pp, c.w, c.sig, 1, 1, h);
    bench_measure(label, "trace_K1", trace_runs, tagb, op_trace1, NULL, &c);
    {
        /* A tag with all n identifiers live, identifier i at depth (i-1) mod h.
           Merge rerandomizes with a fresh N^h-th residue, so a merged tag with
           this plaintext has the distribution of a fresh encryption of it. */
        mpz_t m, term;
        mpz_inits(m, term, NULL);
        for (int i = 1; i <= pp.n; i++) {
            mpz_pow_ui(term, pp.r, (unsigned long)(i - 1));
            mpz_mul(term, term, pp.Np[(i - 1) % h]);
            mpz_add(m, m, term);
        }
        dj_enc(&pp, c.wn, m, h, NULL);
        mpz_clears(m, term, NULL);
    }
    char name[64];
    snprintf(name, sizeof name, "trace_K%d", pp.n);
    bench_measure(label, name, trace_runs, tagb, op_tracen, NULL, &c);

    if (M > 0) {
        int er = runs < 3 ? runs : 3;
        int warm = bench_warmup;
        bench_warmup = 0;
        c.M = M;
        c.threads = threads;
        c.slot = malloc((size_t)M * sizeof(int));
        c.b = malloc((size_t)M * sizeof(int));
        for (int j = 0; j < M; j++) { c.slot[j] = 1 + j % pp.n; c.b[j] = (j % 10 == 0); }
        size_t eb = dj_entry_bytes(&pp);
        dj_entry_init(&pp, &c.E);
        op_entry(&c);
        bench_measure(label, "board_entry", er, eb, op_entry, NULL, &c);
        bench_measure(label, "verify_board_entry", er, eb, op_ventry, NULL, &c);
        dj_entry_free(&pp, &c.E);
        dj_board_init(&c.B, &pp, (size_t)M);
        size_t bb = (size_t)M * eb + dj_count_bytes(&pp) + dj_root_bytes(&pp);
        snprintf(name, sizeof name, "board_M%d_t%d", M, threads);
        bench_measure(label, name, 1, bb, op_board, NULL, &c);
        snprintf(name, sizeof name, "verify_board_M%d_t%d", M, threads);
        bench_measure(label, name, 1, bb, op_vboard, NULL, &c);
        dj_board_free(&c.B, &pp);
        free(c.slot); free(c.b);
        bench_warmup = warm;
    }

    mpz_clears(c.w, c.w2, c.xi, c.ws[0], c.ws[1], c.wn, c.tau, c.prod, c.t, NULL);
    free(c.sig); free(c.in_S);
    dj_free(&pp);
    return 0;
}
