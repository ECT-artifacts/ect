#include "bench_util.h"
#include "rlwe.h"
#include <string.h>

/* bench_rlwe.c -> Ring-LWE benchmarks at one parameter set.
   Usage: ./bench_rlwe profile N0 h log2Q primes runs
   Output: CSV rows of bench_util.h. */

typedef struct {
    const char *profile;
    int N0, h, logQ, L;
    RLWEParams *pp;
    RLWECt w, w2, a, b, wn;
    RLWECt *in[2];
    uint8_t *sig, *in_S;
    size_t siglen;
} Ctx;

static void op_setup(void *o)
{
    Ctx *c = o;
    RLWEParams q;
    rlwe_setup(&q, c->N0, c->h, c->logQ, c->L, 65537);
    rlwe_free(&q);
}
static void op_tag(void *o)    { Ctx *c = o; c->siglen = rlwe_tag(c->pp, &c->w, c->sig, 1, 1, c->pp->h); }
static void op_verify(void *o) { Ctx *c = o; if (!rlwe_verify_tag(c->pp, &c->w, c->sig, c->siglen)) abort(); }
static void op_degrade(void *o){ Ctx *c = o; rlwe_degrade(c->pp, &c->w2, &c->w); }
static void op_merge(void *o)  { Ctx *c = o; rlwe_merge(c->pp, &c->w2, c->in, 2); }
static void op_trace1(void *o) { Ctx *c = o; if (rlwe_trace(c->pp, &c->w, c->in_S) != 1) abort(); }
static void op_tracen(void *o) { Ctx *c = o; if (rlwe_trace(c->pp, &c->wn, c->in_S) != c->pp->n) abort(); }

int main(int argc, char **argv)
{
    if (argc < 7) {
        fprintf(stderr, "usage: %s profile N0 h log2Q primes runs\n", argv[0]);
        return 2;
    }
    Ctx c;
    memset(&c, 0, sizeof c);
    c.profile = argv[1];
    c.N0 = atoi(argv[2]);
    c.h = atoi(argv[3]);
    c.logQ = atoi(argv[4]);
    c.L = atoi(argv[5]);
    int runs = atoi(argv[6]);

    RLWEParams pp;
    rlwe_setup(&pp, c.N0, c.h, c.logQ, c.L, 65537);
    c.pp = &pp;
    char label[96];
    snprintf(label, sizeof label, "rlwe_%s_N0%d_h%d_logQ%d", c.profile, c.N0, c.h, pp.logQ);
    rlwe_ct_init(&pp, &c.w); rlwe_ct_init(&pp, &c.w2);
    rlwe_ct_init(&pp, &c.a); rlwe_ct_init(&pp, &c.b); rlwe_ct_init(&pp, &c.wn);
    c.sig = malloc(pp.sig.sig_bytes);
    c.in_S = malloc((size_t)pp.n);
    size_t tagb = rlwe_tag_bytes(&pp);

    bench_measure(label, "setup", runs < 5 ? runs : 5, 0, op_setup, NULL, &c);
    bench_measure(label, "tag", runs, tagb + pp.sig.sig_bytes, op_tag, NULL, &c);
    bench_measure(label, "verify_signature", runs, tagb + pp.sig.sig_bytes, op_verify, NULL, &c);
    bench_measure(label, "degrade", runs, tagb, op_degrade, NULL, &c);
    rlwe_tag(&pp, &c.a, c.sig, 2, 1, c.h);
    rlwe_tag(&pp, &c.b, c.sig, 3, 1, c.h);
    c.in[0] = &c.a; c.in[1] = &c.b;
    bench_measure(label, "merge_k2", runs, tagb, op_merge, NULL, &c);
    c.siglen = rlwe_tag(&pp, &c.w, c.sig, 1, 1, c.h);
    bench_measure(label, "trace_K1", runs, tagb, op_trace1, NULL, &c);
    {
        /* All N0 identifiers live, identifier i at depth (i-1) mod h; Merge
           rerandomizes, and this tag is a fresh encryption of that plaintext. */
        size_t n = (size_t)c.h * c.N0;
        uint64_t *m = calloc(n, sizeof(uint64_t));
        uint64_t *slots = calloc((size_t)c.h * c.N0, sizeof(uint64_t));
        for (int i = 0; i < c.N0; i++) slots[(size_t)(i % c.h) * c.N0 + i] = 1;
        for (int d = 0; d < c.h; d++) {
            /* inverse slot transform of the depth-d slot vector */
            uint64_t *blkm = m + (size_t)d * c.N0;
            uint64_t *unit = calloc((size_t)c.N0, sizeof(uint64_t));
            for (int i = 0; i < c.N0; i++) {
                if (!slots[(size_t)d * c.N0 + i]) continue;
                rlwe_slot_unit(&pp, unit, i);
                for (int k = 0; k < c.N0; k++) blkm[k] = (blkm[k] + unit[k]) % pp.t;
            }
            free(unit);
        }
        rlwe_enc(&pp, &c.wn, m);
        free(m); free(slots);
    }
    char name[64];
    snprintf(name, sizeof name, "trace_K%d", pp.n);
    bench_measure(label, name, runs, tagb, op_tracen, NULL, &c);

    rlwe_ct_free(&c.w); rlwe_ct_free(&c.w2); rlwe_ct_free(&c.a);
    rlwe_ct_free(&c.b); rlwe_ct_free(&c.wn);
    free(c.sig); free(c.in_S);
    rlwe_free(&pp);
    return 0;
}
