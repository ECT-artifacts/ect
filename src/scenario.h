#pragma once
#include "sim.h"
#include <time.h>

double timespec_ms(struct timespec start, struct timespec end);

typedef struct {
    int    n_blocks_committed;
    int    n_traced_expected;
    int    n_traced_recovered;
    int    n_expired;
    double ms_issue;
    double ms_spend;
    double ms_trace;
    double ms_total;
    double ms_validator;
} ScenarioReport;

void report_print(const ScenarioReport *r);
void run_ect_scenario(ECTParams *pp);

/* BENCH_RUNS: repetitions for the per-operation benchmarks (B1-B3, B8, B9).
   Was 100 for the 3071-bit parameters; every group operation is roughly 7x
   more expensive at the 5807-bit modulus, so this is dialled down. */
#define BENCH_RUNS 10

/* TRACE_RUNS: repetitions for the sweeps that call Trace (B5-B7).  A single
   Trace now costs ~1 s, dominated by the h Pohlig-Hellman projections, so
   these need far fewer repetitions to stay tractable. */
#define TRACE_RUNS 3

/* BENCH_POP_ALLOC_MAX: cap on how many background users B5/B6 materialise.
   No timed section of either benchmark reads the background set -- it is
   allocated, filled from a prototype tag, and freed -- so capping it leaves
   every measurement unchanged while avoiding 4-8 GB of bn_t allocation per
   run (each bn_t reserves ~2 KB under BN_PRECI 8192 / BN_MAGNI DOUBLE).
   The nominal population is still reported, since it is what sets K. */
#define BENCH_POP_ALLOC_MAX 16384

/* B1: Authority -> Tag (enc + sign) */
void bench_authority_tag(ECTParams *pp);

/* B2: User -> Degrade (pow_q+rerand+nizk), Merge,  Degrade after merge */
void bench_user_degrade_merge(ECTParams *pp);

/* B3: Validator -> schnorr_verify, nizk_verify, hash_check, full block verify */
void bench_validator_verify(ECTParams *pp);

/* B4: Sizes -> 6 quantities in bytes/KB */
void bench_sizes(ECTParams *pp);

/* B5: Population 2^4..2^13, K=log2(pop) */
void bench_population_scale(ECTParams *pp);

/* B6: Fixed pop=2^10, vary K 1..64 */
void bench_ktraced_scale(ECTParams *pp);

/* B7: Fixed pop=2^10, K=6, vary depth 1..H */
void bench_degradation_depth(ECTParams *pp);

/* B8: Vary merge size 1..10, then degrade */
void bench_merge_then_degrade(ECTParams *pp);

/* B9: Full with overhead*/
void bench_full_fidelity(ECTParams *pp);

/* B_nizk: NIZK prove/verify scaling with m=1,2,4,8,16,32,64,128 */
void bench_nizk_scale(ECTParams *pp);

/* run_benchmarks: B1..B8 */
void run_benchmarks(ECTParams *pp);
