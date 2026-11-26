#define _POSIX_C_SOURCE 200809L
/* bench_trace.c -> Trace cost vs number of simultaneously traced users K,
 * swept over the multiplicity bound M.
 *
 *   make bench_trace && ./bench_trace
 *
 * Emits CSV: mult_max,K,table_entries,probes,trace_ms
 */
#include "ect.h"
#include "elgamal.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

static ECTParams pp;

static double now_ms(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec * 1000.0 + t.tv_nsec / 1e6;
}

static double trace_once(int K, int *found)
{
    ECTTag *tags = (ECTTag *)malloc((size_t)K * sizeof(ECTTag));
    ECTTag **ptr = (ECTTag **)malloc((size_t)K * sizeof(ECTTag *));
    ECTTag merged;

    ect_clear_tagged(&pp);
    for (int i = 0; i < K; i++) {
        ect_tag_new(&tags[i]);
        ect_tag(&tags[i], &pp, i + 1, 1);
        ptr[i] = &tags[i];
    }
    ect_tag_new(&merged);
    ect_merge(&merged, ptr, K, &pp);

    double t0 = now_ms();
    *found = ect_trace_silent(&pp, &merged);
    double dt = now_ms() - t0;

    ect_tag_free(&merged);
    for (int i = 0; i < K; i++) ect_tag_free(&tags[i]);
    free(tags); free(ptr);
    return dt;
}

int main(int argc, char **argv)
{
    int runs   = (argc > 1) ? atoi(argv[1]) : 3;
    int max_k  = (argc > 2) ? atoi(argv[2]) : 1 << 30;
    core_init();
    ect_params_new(&pp);
    ect_setup(&pp);

    fprintf(stderr, "H=%d R_BASE=%d N=%d  |p|=%d bits  (%d run%s each)\n\n",
            pp.h, pp.r, pp.n, P_BITS, runs, runs == 1 ? "" : "s");
    printf("mult_max,K,table_entries,probes_per_depth,trace_ms,recovered\n");

    int m_vals[] = {1, 2, 4};
    int k_vals[] = {1, 2, 4, 8, 12, 16, 20, 24, 28, 32};

    for (size_t mi = 0; mi < sizeof(m_vals) / sizeof(*m_vals); mi++) {
        ect_mult_max = m_vals[mi];
        double lg = log2((double)(m_vals[mi] + 1));
        for (size_t ki = 0; ki < sizeof(k_vals) / sizeof(*k_vals); ki++) {
            int K = k_vals[ki];
            if (K > pp.n || K > max_k) continue;

            /* mirror mitm_split so the CSV records the actual shape */
            int a = (K + 1) / 2;
            while (a > 0 && (double)a * lg > (double)MITM_CAP_LOG2) a--;
            int b = K - a;
            double tab = pow(m_vals[mi] + 1.0, a);
            double prb = pow(m_vals[mi] + 1.0, b);
            /* skip cells that would run for hours */
            if (tab * 1e-6 * 5 + prb * 1e-6 * 5 * pp.h > 120.0) continue;

            double best = 1e30;
            int found = 0;
            for (int r = 0; r < runs; r++) {
                int f;
                double dt = trace_once(K, &f);
                if (dt < best) best = dt;
                found = f;
            }
            printf("%d,%d,%.0f,%.0f,%.1f,%d\n",
                   m_vals[mi], K, tab, prb, best, found);
            fflush(stdout);
        }
    }

    ect_params_free(&pp);
    core_clean();
    return 0;
}
