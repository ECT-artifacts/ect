#pragma once
#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* bench_util.h -> timing helpers shared by the benchmarks.  Each measurement
   runs the operation once untimed (unless bench_warmup is 0), then `runs`
   timed times, and prints one
   CSV row: label,operation,runs,mean_ms,sd_ms,min_ms,bytes. */

typedef void (*bench_op)(void *ctx);

/* Set to 0 to skip the untimed execution, for operations of minutes. */
static int bench_warmup = 1;

static inline double bench_now_ms(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec * 1000.0 + (double)t.tv_nsec / 1e6;
}

/* reset (may be NULL) runs before each execution, outside the timed interval. */
static inline void bench_measure(const char *label, const char *name, int runs,
                                 size_t bytes, bench_op op, bench_op reset, void *ctx)
{
    double *s = malloc((size_t)runs * sizeof(double));
    if (!s) exit(1);
    if (bench_warmup) {
        if (reset) reset(ctx);
        op(ctx);
    }
    double sum = 0, min = INFINITY;
    for (int i = 0; i < runs; i++) {
        if (reset) reset(ctx);
        double t0 = bench_now_ms();
        op(ctx);
        s[i] = bench_now_ms() - t0;
        sum += s[i];
        if (s[i] < min) min = s[i];
    }
    double mean = sum / runs, var = 0;
    for (int i = 0; i < runs; i++) var += (s[i] - mean) * (s[i] - mean);
    double sd = runs > 1 ? sqrt(var / (runs - 1)) : 0;
    printf("%s,%s,%d,%.3f,%.3f,%.3f,%zu\n", label, name, runs, mean, sd, min, bytes);
    fflush(stdout);
    free(s);
}
