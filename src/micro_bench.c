#define _POSIX_C_SOURCE 200809L

#include "ect.h"
#include "block.h"
#include "elgamal.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef void (*bench_fn)(void *);

static double now_ms(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec * 1000.0 + (double)t.tv_nsec / 1e6;
}

static void measure(const char *name, int runs, size_t bytes,
                    bench_fn fn, void *ctx)
{
    double *samples = malloc((size_t)runs * sizeof(*samples));
    if (!samples) exit(1);

    fn(ctx); /* warm-up */
    double sum = 0.0, min = INFINITY;
    for (int i = 0; i < runs; i++) {
        double begin = now_ms();
        fn(ctx);
        double elapsed = now_ms() - begin;
        samples[i] = elapsed;
        sum += elapsed;
        if (elapsed < min) min = elapsed;
    }
    double mean = sum / runs, variance = 0.0;
    for (int i = 0; i < runs; i++) {
        double delta = samples[i] - mean;
        variance += delta * delta;
    }
    double sd = runs > 1 ? sqrt(variance / (runs - 1)) : 0.0;
    printf("%s,%d,%.3f,%.3f,%.3f,%zu\n",
           name, runs, mean, sd, min, bytes);
    fflush(stdout);
    free(samples);
}

typedef struct { ECTParams *pp; ECTTag tag; int id; } TagCtx;

static void run_tag(void *opaque)
{
    TagCtx *ctx = opaque;
    ect_tag(&ctx->tag, ctx->pp, ctx->id, 1);
}

static void run_sig_verify(void *opaque)
{
    TagCtx *ctx = opaque;
    if (!ect_verify_tag_sig(&ctx->tag, ctx->pp)) abort();
}

typedef struct { ECTParams *pp; ECTTag tag; bn_t c1, c2, rho; } DegradeCtx;

static void reset_degrade(DegradeCtx *ctx)
{
    ctx->tag.n_proofs = 0;
    ctx->tag.hops = 0;
    ect_tag(&ctx->tag, ctx->pp, 1, 1);
    bn_copy(ctx->c1, ctx->tag.c1);
    bn_copy(ctx->c2, ctx->tag.c2);
}

static void run_degrade_prove(void *opaque)
{
    DegradeCtx *ctx = opaque;
    ect_degrade(&ctx->tag, ctx->pp, &ctx->c1, &ctx->c2, 1, 0);
}

static void run_degrade(void *opaque)
{
    DegradeCtx *ctx = opaque;
    ciphertext_pow_q(ctx->tag.c1, ctx->tag.c2, ctx->pp->q, ctx->pp->p);
    elgamal_rerand(ctx->tag.c1, ctx->tag.c2, ctx->rho,
                   ctx->pp->pk, ctx->pp->g, ctx->pp->qh, ctx->pp->p);
}

static void measure_degrade(const char *name, int runs, size_t bytes,
                            bench_fn fn, DegradeCtx *ctx)
{
    double *samples = malloc((size_t)runs * sizeof(*samples));
    if (!samples) exit(1);

    reset_degrade(ctx);
    fn(ctx); /* warm-up */
    double sum = 0.0, min = INFINITY;
    for (int i = 0; i < runs; i++) {
        reset_degrade(ctx); /* deliberately outside the timed interval */
        double begin = now_ms();
        fn(ctx);
        double elapsed = now_ms() - begin;
        samples[i] = elapsed;
        sum += elapsed;
        if (elapsed < min) min = elapsed;
    }
    double mean = sum / runs, variance = 0.0;
    for (int i = 0; i < runs; i++) {
        double delta = samples[i] - mean;
        variance += delta * delta;
    }
    double sd = runs > 1 ? sqrt(variance / (runs - 1)) : 0.0;
    printf("%s,%d,%.3f,%.3f,%.3f,%zu\n",
           name, runs, mean, sd, min, bytes);
    fflush(stdout);
    free(samples);
}

static void run_proof_verify(void *opaque)
{
    DegradeCtx *ctx = opaque;
    if (!ect_verify_proofs(&ctx->tag, ctx->pp)) abort();
}

typedef struct { ECTParams *pp; ECTTag *tags; ECTTag **inputs; ECTTag out; int k; } MergeCtx;

static void run_merge(void *opaque)
{
    MergeCtx *ctx = opaque;
    ect_merge(&ctx->out, ctx->inputs, ctx->k, ctx->pp);
}

typedef struct { ECTParams *pp; ECTTag merged; int expected; } TraceCtx;

static void run_trace(void *opaque)
{
    TraceCtx *ctx = opaque;
    if (ect_trace_silent(ctx->pp, &ctx->merged) != ctx->expected) abort();
}

int main(int argc, char **argv)
{
    int runs = argc > 1 ? atoi(argv[1]) : 10;
    int trace_runs = argc > 2 ? atoi(argv[2]) : 3;
    int trace_k = argc > 3 ? atoi(argv[3]) : N;
    if (runs < 1 || trace_runs < 1 || trace_k < 1 || trace_k > N) {
        fprintf(stderr, "usage: %s [runs>=1] [trace-runs>=1] [1<=K<=N]\n", argv[0]);
        return 2;
    }

    if (core_init() != RLC_OK) return 1;
    ECTParams pp;
    ect_params_new(&pp);
    ect_setup(&pp);

    fprintf(stderr, "parameters: H=%d N=%d r=%d |p|=%d |q|=%d\n",
            H, N, R_BASE, P_BITS, Q_BITS);
    printf("operation,runs,mean_ms,sd_ms,min_ms,bytes\n");

    TagCtx tag = {.pp = &pp, .id = 1};
    ect_tag_new(&tag.tag);
    measure("tag", runs,
            3u * P_BYTES + (size_t)((H * Q_BITS + 7) / 8),
            run_tag, &tag);
    measure("verify_signature", runs,
            2u * P_BYTES + P_BYTES + (size_t)((H * Q_BITS + 7) / 8),
            run_sig_verify, &tag);

    DegradeCtx degrade = {.pp = &pp};
    ect_tag_new(&degrade.tag);
    bn_null(degrade.c1); bn_new(degrade.c1);
    bn_null(degrade.c2); bn_new(degrade.c2);
    bn_null(degrade.rho); bn_new(degrade.rho);
    reset_degrade(&degrade);
    measure_degrade("degrade", runs, 2u * P_BYTES,
                    run_degrade, &degrade);
    measure_degrade("degrade_and_prove", runs,
                    4u * P_BYTES + 2u * (size_t)((H * Q_BITS + 7) / 8),
                    run_degrade_prove, &degrade);
    measure("verify_degradation_proof", runs,
            2u * P_BYTES + 2u * (size_t)((H * Q_BITS + 7) / 8),
            run_proof_verify, &degrade);

    const int max_k = 10;
    const int fixture_k = trace_k > max_k ? trace_k : max_k;
    MergeCtx merge = {.pp = &pp, .k = 2};
    merge.tags = calloc((size_t)fixture_k, sizeof(*merge.tags));
    merge.inputs = calloc((size_t)fixture_k, sizeof(*merge.inputs));
    if (!merge.tags || !merge.inputs) return 1;
    ect_tag_new(&merge.out);
    for (int i = 0; i < fixture_k; i++) {
        ect_tag_new(&merge.tags[i]);
        ect_tag(&merge.tags[i], &pp, i + 1, 1);
        merge.inputs[i] = &merge.tags[i];
    }
    measure("merge_k2", runs, 2u * P_BYTES, run_merge, &merge);
    merge.k = max_k;
    measure("merge_k10", runs, 2u * P_BYTES, run_merge, &merge);

    TraceCtx trace = {.pp = &pp, .expected = trace_k};
    ect_tag_new(&trace.merged);
    ect_clear_tagged(&pp);
    for (int i = 0; i < trace_k; i++)
        ect_tag(&merge.tags[i], &pp, i + 1, 1);
    ect_merge(&trace.merged, merge.inputs, trace_k, &pp);
    ect_mult_max = 1; /* this workload contains distinct identifiers once */
    char trace_name[32];
    snprintf(trace_name, sizeof(trace_name), "trace_k%d", trace_k);
    measure(trace_name, trace_runs, 2u * P_BYTES, run_trace, &trace);

    ect_tag_free(&trace.merged);
    for (int i = 0; i < fixture_k; i++) ect_tag_free(&merge.tags[i]);
    ect_tag_free(&merge.out);
    free(merge.tags); free(merge.inputs);
    bn_free(degrade.c1); bn_free(degrade.c2); bn_free(degrade.rho);
    ect_tag_free(&degrade.tag);
    ect_tag_free(&tag.tag);
    ect_params_free(&pp);
    core_clean();
    return 0;
}
