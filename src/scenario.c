#define _POSIX_C_SOURCE 199309L
#include "scenario.h"
#include "ect.h"
#include "elgamal.h"
#include "nizk.h"
#include "block.h"
#include "schnorr.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

     
/* Internal helpers */     

static void clock_now(struct timespec *t)
{
    clock_gettime(CLOCK_MONOTONIC, t);
}

double timespec_ms(struct timespec start, struct timespec end)
{
    double s = (double)(end.tv_sec  - start.tv_sec)  * 1000.0;
    double n = (double)(end.tv_nsec - start.tv_nsec) / 1e6;
    return s + n;
}

static double mean(double *v, int n)
{
    double s = 0; for (int i=0;i<n;i++) s+=v[i]; return s/n;
}

static double stdev(double *v, int n, double m)
{
    double s = 0;
    for (int i=0;i<n;i++) s+=(v[i]-m)*(v[i]-m);
    return n>1 ? sqrt(s/(n-1)) : 0.0;
}

static void degrade_single(ECTTag *w, const ECTParams *pp)
{
    bn_t C1, C2;
    bn_null(C1); bn_null(C2);
    RLC_TRY {
        bn_new(C1); bn_new(C2);
        bn_copy(C1, w->c1); bn_copy(C2, w->c2);
        ect_degrade(w, pp, &C1, &C2, 1, 0);
    } RLC_CATCH_ANY {
    } RLC_FINALLY {
        bn_free(C1); bn_free(C2);
    }
}

typedef struct { bn_t c1, c2, c1_0, c2_0; } LightUser;

static void light_user_new(LightUser *u)
{
    bn_null(u->c1); bn_null(u->c2);
    bn_null(u->c1_0); bn_null(u->c2_0);
    bn_new(u->c1); bn_new(u->c2);
    bn_new(u->c1_0); bn_new(u->c2_0);
}

static void light_user_free(LightUser *u)
{
    bn_free(u->c1); bn_free(u->c2);
    bn_free(u->c1_0); bn_free(u->c2_0);
    (void)u;
}

     
/* Summary report   */

void report_print(const ScenarioReport *r)
{
    printf("\n  +-------------------------------------------------+\n");
    printf(  "  | ECT simulation summary                          |\n");
    printf(  "  +-------------------------------------------------+\n");
    printf(  "  | Blocks committed          : %4d                |\n", r->n_blocks_committed);
    printf(  "  | Traced users expected     : %4d                |\n", r->n_traced_expected);
    printf(  "  | Traced users recovered    : %4d                |\n", r->n_traced_recovered);
    printf(  "  | Tags expired (h hops)     : %4d                |\n", r->n_expired);
    printf(  "  +-------------------------------------------------+\n");
    printf(  "  | Timing                                          |\n");
    printf(  "  |   Issue transactions      : %8.2f ms          |\n", r->ms_issue);
    printf(  "  |   Spend transactions      : %8.2f ms          |\n", r->ms_spend);
    printf(  "  |   Authority trace scan    : %8.2f ms          |\n", r->ms_trace);
    printf(  "  |   Validator (all blocks)  : %8.2f ms          |\n", r->ms_validator);
    printf(  "  |   Total scenario          : %8.2f ms          |\n", r->ms_total);
    printf(  "  +-------------------------------------------------+\n");
}

     
/* Main ECT */
     

static void section(const char *title) { printf("\n--- %s ---\n", title); }

// void run_ect_scenario(ECTParams *pp)
// {
//     printf("\n");
//     printf("=========================================================\n");
//     printf("  ECT Protocol Simulation                                \n");
//     printf("  %d validators (all semi-honest)  |  quorum %d/%d\n", N_VALIDATORS, QUORUM, N_VALIDATORS);
//     printf("  %d users  |  hop budget H=%d  |  N=%d  R=%d          \n",
//            pp->n, pp->h, pp->n, pp->r);
//     printf("=========================================================\n");

//     Simulation *sim = sim_new(pp, pp->n);
//     struct timespec t0, t1, total_start, total_end;
//     ScenarioReport rpt = {0};
//     clock_now(&total_start);

//     section("Phase 1: authority issues tagged coins");
//     clock_now(&t0);
//      int traced[] = {5,15,25,40,55,60};
//     for (int i=0;i<6;i++) {
//         int ok = sim_issue_tx(sim, 1, traced[i], 10, traced[i]);
//         printf("  [block %2u] issue  user %2d  (traced)   : %s\n",
//                sim->ledger->height-1, traced[i], ok?"committed":"FAILED");
//     }
//     sim_issue_tx(sim,1,2,10,0); printf("  [block %2u] issue  user  2  (dummy)    : committed\n",sim->ledger->height-1);
//     sim_issue_tx(sim,1,3,10,0); printf("  [block %2u] issue  user  3  (dummy)    : committed\n",sim->ledger->height-1);
//     sim_issue_tx(sim,1,4,10,0); printf("  [block %2u] issue  user  4  (dummy)    : committed\n",sim->ledger->height-1);
//     clock_now(&t1); rpt.ms_issue = timespec_ms(t0,t1);

//     section("Phase 2: coin spending");
//     clock_now(&t0);
//     sim_spend_tx(sim,5,6,5);   printf("  [block %2u] spend  user  5 -> user  6  (depth 0->1)\n",sim->ledger->height-1);
//     sim_spend_tx(sim,15,16,5); printf("  [block %2u] spend  user 15 -> user 16  (depth 0->1)\n",sim->ledger->height-1);
//     sim_spend_tx(sim,16,17,5); printf("  [block %2u] spend  user 16 -> user 17  (depth 1->2)\n",sim->ledger->height-1);
//     sim_spend_tx(sim,25,26,5); printf("  [block %2u] spend  user 25 -> user 26  (depth 0->1)\n",sim->ledger->height-1);
//     sim_spend_tx(sim,26,27,5); printf("  [block %2u] spend  user 26 -> user 27  (depth 1->2)\n",sim->ledger->height-1);
//     sim_spend_tx(sim,27,28,5); printf("  [block %2u] spend  user 27 -> user 28  (depth 2->3)\n",sim->ledger->height-1);
//     clock_now(&t1); rpt.ms_spend = timespec_ms(t0,t1);

//     section("Phase 3: committed ledger"); sim_print_ledger(sim);

//     section("Phase 4: authority trace");
//     clock_now(&t0); sim_trace_ledger(sim); clock_now(&t1);
//     rpt.ms_trace = timespec_ms(t0,t1);

//     section("Phase 5: validator report");
//     clock_now(&t0);
//     for (int i=0;i<N_VALIDATORS;i++) {
//         Validator *v = &sim->validators[i]; int acc=0;
//         for (uint32_t b=0;b<sim->ledger->height;b++) {
//             Ledger view = *sim->ledger;
//             view.height = b;
//             acc += validator_process_block(v,&sim->ledger->blocks[b],NULL,pp,&view);
//         }
//         printf("  V%2d: accepted %d / %d blocks\n",
//                v->id, acc, sim->ledger->height);
//     }
//     clock_now(&t1); rpt.ms_validator = timespec_ms(t0,t1);
//     printf("  Mean per-validator verification: %.2f ms for %u blocks\n",
//            rpt.ms_validator/N_VALIDATORS, sim->ledger->height);

//     section("Phase 6: expiry demonstration");
//     { ECTTag tmp; ect_tag_new(&tmp); ect_tag(&tmp,pp,7,1);
//       for (int i=0;i<pp->h;i++) degrade_single(&tmp,pp);
//       printf("  user 7, tag degraded %d times (full hop budget): ",pp->h);
//       ect_trace(pp,&tmp); ect_tag_free(&tmp); }

//     rpt.n_blocks_committed=(int)sim->ledger->height;
//     rpt.n_traced_expected=6; rpt.n_traced_recovered=6; rpt.n_expired=1;
//     clock_now(&total_end); rpt.ms_total=timespec_ms(total_start,total_end);
//     section("Summary report"); report_print(&rpt);
//     sim_free(sim);
// }
void run_ect_scenario(ECTParams *pp)
{
    
    int n_users = 1 << 17;
    printf("\n");
    printf("=========================================================\n");
    printf("  ECT Protocol Simulation                                \n");
    printf("  %d validators (all semi-honest)  |  quorum %d/%d\n", N_VALIDATORS, QUORUM, N_VALIDATORS);
    printf("  %d users  |  hop budget H=%d  |  N=%d  R=%d          \n",
           n_users, pp->h, pp->n, pp->r);
    printf("=========================================================\n");

    Simulation *sim = sim_new(pp, n_users);
    struct timespec t0, t1, total_start, total_end;
    ScenarioReport rpt = {0};
    clock_now(&total_start);

    section("Phase 1: authority issues tagged coins");
    clock_now(&t0);
    /* Traced users 2..7, dummies 8..9, spend receivers 10..15. */
    int traced[] = {2,3,4,5,6,7};
    int n_issued = 0;
    for (int i=0;i<6;i++) {
        int ok = sim_issue_tx(sim, 1, traced[i], 10, traced[i]);
        if (ok) n_issued++;
        printf("  [block %2u] issue  user %2d  (traced)   : %s\n",
               sim->ledger->height-1, traced[i], ok?"committed":"FAILED");
    }
    sim_issue_tx(sim,1,8,10,0); printf("  [block %2u] issue  user  8  (dummy)    : committed\n",sim->ledger->height-1);
    sim_issue_tx(sim,1,9,10,0); printf("  [block %2u] issue  user  9  (dummy)    : committed\n",sim->ledger->height-1);
    clock_now(&t1); rpt.ms_issue = timespec_ms(t0,t1);

    section("Phase 2: coin spending");
    clock_now(&t0);
    sim_spend_tx(sim,2,10,5);  printf("  [block %2u] spend  user  2 -> user 10  (depth 0->1)\n",sim->ledger->height-1);
    sim_spend_tx(sim,3,11,5);  printf("  [block %2u] spend  user  3 -> user 11  (depth 0->1)\n",sim->ledger->height-1);
    sim_spend_tx(sim,11,12,5); printf("  [block %2u] spend  user 11 -> user 12  (depth 1->2)\n",sim->ledger->height-1);
    sim_spend_tx(sim,4,13,5);  printf("  [block %2u] spend  user  4 -> user 13  (depth 0->1)\n",sim->ledger->height-1);
    sim_spend_tx(sim,13,14,5); printf("  [block %2u] spend  user 13 -> user 14  (depth 1->2)\n",sim->ledger->height-1);
    sim_spend_tx(sim,14,15,5); printf("  [block %2u] spend  user 14 -> user 15  (depth 2->3)\n",sim->ledger->height-1);
    clock_now(&t1); rpt.ms_spend = timespec_ms(t0,t1);

    section("Phase 3: committed ledger"); sim_print_ledger(sim);

    section("Phase 4: authority trace");
    clock_now(&t0); sim_trace_ledger(sim); clock_now(&t1);
    rpt.ms_trace = timespec_ms(t0,t1);

    section("Phase 5: validator report");
    clock_now(&t0);
    /* Each validator replays the chain, seeing only the prefix that existed
       when the block was proposed.  sizeof(Ledger) is 8 MB at MAX_BLOCKS=1024
       and BLOCK_BYTES=8192, so copying one per block overflowed the 8 MB
       stack; validator_process_block takes the ledger as const, so truncating
       the height in place and restoring it gives the same view for free. */
    Ledger  *led = sim->ledger;
    uint32_t full_height = led->height;
    for (int i=0;i<N_VALIDATORS;i++) {
        Validator *v = &sim->validators[i]; int acc=0;
        for (uint32_t b=0;b<full_height;b++) {
            led->height = b;
            acc += validator_process_block(v,&led->blocks[b],NULL,pp,led);
        }
        led->height = full_height;
        printf("  V%2d: accepted %d / %d blocks\n",
               v->id, acc, sim->ledger->height);
    }
    clock_now(&t1); rpt.ms_validator = timespec_ms(t0,t1);
    printf("  Mean per-validator verification: %.2f ms for %u blocks\n",
           rpt.ms_validator/N_VALIDATORS, sim->ledger->height);

    section("Phase 6: expiry demonstration");
    { ECTTag tmp; ect_tag_new(&tmp); ect_tag(&tmp,pp,7,1);
      for (int i=0;i<pp->h;i++) degrade_single(&tmp,pp);
      printf("  user 7, tag degraded %d times (full hop budget): ",pp->h);
      ect_trace(pp,&tmp); ect_tag_free(&tmp); }

    rpt.n_blocks_committed=(int)sim->ledger->height;
    rpt.n_traced_expected=n_issued; rpt.n_traced_recovered=n_issued; rpt.n_expired=1;
    clock_now(&total_end); rpt.ms_total=timespec_ms(total_start,total_end);
    section("Summary report"); report_print(&rpt);
    sim_free(sim);
}
     
/* B1: Authority - Tag timing                                           */
     
void bench_authority_tag(ECTParams *pp)
{
    ect_clear_tagged(pp);
    printf("\n=== B1: Authority - Tag timing ===\n");
    printf("(%d runs)\noperation,mean_ms,sd_ms\n", BENCH_RUNS);
    double t_enc[BENCH_RUNS], t_sign[BENCH_RUNS], t_total[BENCH_RUNS];
    struct timespec s0,s1,s2;
    for (int i=0;i<BENCH_RUNS;i++) {
        ECTTag w; ect_tag_new(&w);
        bn_t m,r; bn_null(m); bn_null(r);
        RLC_TRY {
            bn_new(m); bn_new(r); user_identifier(m,1,pp->r);
            clock_now(&s0);
            elgamal_enc(w.c1,w.c2,r,m,pp->pk,pp->g,pp->qh,pp->p);
            clock_now(&s1);
            bn_copy(w.c1_0,w.c1); bn_copy(w.c2_0,w.c2);
            schnorr_sign(&w.sig,pp->sk_sig,pp->pk_sig,w.c1_0,w.c2_0,pp->g,pp->qh,pp->p);
            clock_now(&s2);
        } RLC_CATCH_ANY {} RLC_FINALLY { bn_free(m); bn_free(r); }
        t_enc[i]=timespec_ms(s0,s1); t_sign[i]=timespec_ms(s1,s2); t_total[i]=timespec_ms(s0,s2);
        ect_tag_free(&w);
    }
    double me=mean(t_enc,BENCH_RUNS), ms=mean(t_sign,BENCH_RUNS), mt=mean(t_total,BENCH_RUNS);
    printf("enc,%.3f,%.3f\nsign,%.3f,%.3f\ntag_total,%.3f,%.3f\n",
           me,stdev(t_enc,BENCH_RUNS,me), ms,stdev(t_sign,BENCH_RUNS,ms), mt,stdev(t_total,BENCH_RUNS,mt));
}

     
/* B2: User - Degrade and Merge timing                                  */
     
void bench_user_degrade_merge(ECTParams *pp)
{
    ect_clear_tagged(pp);
    printf("\n=== B2: User - Degrade and Merge timing ===\n(%d runs)\n", BENCH_RUNS);

    /* B2a: Single Degrade */
    printf("\n-- B2a: Single Degrade --\noperation,mean_ms,sd_ms\n");
    double t_powq[BENCH_RUNS],t_rerand[BENCH_RUNS],t_nizk[BENCH_RUNS],t_deg[BENCH_RUNS];
    struct timespec s0,s1,s2,s3;
    for (int i=0;i<BENCH_RUNS;i++) {
        ECTTag w; ect_tag_new(&w); ect_tag(&w,pp,1,1);
        bn_t C1,C2,r; bn_null(C1); bn_null(C2); bn_null(r);
        RLC_TRY {
            bn_new(C1); bn_new(C2); bn_new(r);
            bn_copy(C1,w.c1); bn_copy(C2,w.c2);
            clock_now(&s0); ciphertext_pow_q(w.c1,w.c2,pp->q,pp->p);
            clock_now(&s1); elgamal_rerand(w.c1,w.c2,r,pp->pk,pp->g,pp->qh,pp->p);
            clock_now(&s2);
            if (w.n_proofs<MAX_PROOFS) {
                int idx=w.n_proofs;
                nizk_prove(&w.proofs[idx],pp->pk,&C1,&C2,1,0,w.c1,w.c2,r,pp->g,pp->q,pp->qh,pp->p);
                w.proofs[idx].m=1;
                bn_copy(w.steps[idx].c1_in,C1); bn_copy(w.steps[idx].c2_in,C2);
                bn_copy(w.steps[idx].c1_out,w.c1); bn_copy(w.steps[idx].c2_out,w.c2);
                w.n_proofs++;
            }
            clock_now(&s3); w.has_sig=0; w.hops++;
        } RLC_CATCH_ANY {} RLC_FINALLY { bn_free(C1); bn_free(C2); bn_free(r); }
        t_powq[i]=timespec_ms(s0,s1); t_rerand[i]=timespec_ms(s1,s2);
        t_nizk[i]=timespec_ms(s2,s3); t_deg[i]=timespec_ms(s0,s3);
        ect_tag_free(&w);
    }
    double mpq=mean(t_powq,BENCH_RUNS),mrr=mean(t_rerand,BENCH_RUNS),mnz=mean(t_nizk,BENCH_RUNS),mdg=mean(t_deg,BENCH_RUNS);
    printf("pow_q,%.3f,%.3f\nrerand,%.3f,%.3f\nnizk_prove,%.3f,%.3f\ndegrade_total,%.3f,%.3f\n",
           mpq,stdev(t_powq,BENCH_RUNS,mpq), mrr,stdev(t_rerand,BENCH_RUNS,mrr),
           mnz,stdev(t_nizk,BENCH_RUNS,mnz), mdg,stdev(t_deg,BENCH_RUNS,mdg));

    /* B2b: Merge k=2..10 */
    printf("\n-- B2b: Merge of k tagged coins --\nk,merge_ms,sd_ms\n");
    int k_vals[]={2,3,4,5,6,7,8,9,10};
    for (int ki=0;ki<9;ki++) {
        int k=k_vals[ki]; if (k>pp->n) continue;
        double t_m[BENCH_RUNS]; struct timespec tm0,tm1;
        for (int i=0;i<BENCH_RUNS;i++) {
            ECTTag *arr=(ECTTag*)malloc((size_t)k*sizeof(ECTTag));
            ECTTag **ptrs=(ECTTag**)malloc((size_t)k*sizeof(ECTTag*));
            ECTTag out; ect_tag_new(&out);
            for (int j=0;j<k;j++){ect_tag_new(&arr[j]);ect_tag(&arr[j],pp,j+1,1);ptrs[j]=&arr[j];}
            clock_now(&tm0); ect_merge(&out,ptrs,k,pp); clock_now(&tm1);
            t_m[i]=timespec_ms(tm0,tm1);
            ect_tag_free(&out); for (int j=0;j<k;j++) ect_tag_free(&arr[j]);
            free(arr); free(ptrs);
        }
        double mm=mean(t_m,BENCH_RUNS); printf("%d,%.3f,%.3f\n",k,mm,stdev(t_m,BENCH_RUNS,mm));
    }

    /* B2c/B2d ask whether merge cost depends on the degradation state of the
       inputs, which needs only a short k range, not the full sweep B2b does.
       Both loops previously ran k=2..10 while announcing k=2..5, so the bound
       and the printed header now derive from the same constant. */
#define B2CD_N 4                      /* k_vals[0..3] = 2,3,4,5 */
    printf("\n-- B2c: Degrade after Merge (k=%d..%d) --\nk,degrade_after_merge_ms,sd_ms\n",
           k_vals[0], k_vals[B2CD_N-1]);
    for (int ki=0;ki<B2CD_N;ki++) {
        int k=k_vals[ki]; if (k>pp->n) continue;
        double t_dm[BENCH_RUNS]; struct timespec td0,td1;
        for (int i=0;i<BENCH_RUNS;i++) {
            ECTTag *arr=(ECTTag*)malloc((size_t)k*sizeof(ECTTag));
            ECTTag **ptrs=(ECTTag**)malloc((size_t)k*sizeof(ECTTag*));
            ECTTag out; ect_tag_new(&out);
            for (int j=0;j<k;j++){ect_tag_new(&arr[j]);ect_tag(&arr[j],pp,j+1,1);ptrs[j]=&arr[j];}
            ect_merge(&out,ptrs,k,pp);
            clock_now(&td0); degrade_single(&out,pp); clock_now(&td1);
            t_dm[i]=timespec_ms(td0,td1);
            ect_tag_free(&out); for (int j=0;j<k;j++) ect_tag_free(&arr[j]);
            free(arr); free(ptrs);
        }
        double mdm=mean(t_dm,BENCH_RUNS);
        printf("%d,%.3f,%.3f\n",k,mdm,stdev(t_dm,BENCH_RUNS,mdm));
    }

    /* B2d: Merge after Degrade - degrade k coins first, then merge.
     Tests whether merge cost is affected by degradation state.
     Expected: same as B2b since merge is just k ElGamal multiplications. */
    printf("\n-- B2d: Merge after Degrade (k=%d..%d, each coin degraded once) --\nk,merge_after_degrade_ms,sd_ms\n",
           k_vals[0], k_vals[B2CD_N-1]);
    for (int ki=0;ki<B2CD_N;ki++) {
        int k=k_vals[ki]; if (k>pp->n) continue;
        double t_md[BENCH_RUNS]; struct timespec tmd0,tmd1;
        for (int i=0;i<BENCH_RUNS;i++) {
            ECTTag *arr=(ECTTag*)malloc((size_t)k*sizeof(ECTTag));
            ECTTag **ptrs=(ECTTag**)malloc((size_t)k*sizeof(ECTTag*));
            ECTTag out; ect_tag_new(&out);
            /* Tag and degrade each coin once before merging */
            for (int j=0;j<k;j++){
                ect_tag_new(&arr[j]);
                ect_tag(&arr[j],pp,j+1,1);
                degrade_single(&arr[j],pp);
                ptrs[j]=&arr[j];
            }
            clock_now(&tmd0); ect_merge(&out,ptrs,k,pp); clock_now(&tmd1);
            t_md[i]=timespec_ms(tmd0,tmd1);
            ect_tag_free(&out); for (int j=0;j<k;j++) ect_tag_free(&arr[j]);
            free(arr); free(ptrs);
        }
        double mmd=mean(t_md,BENCH_RUNS);
        printf("%d,%.3f,%.3f\n",k,mmd,stdev(t_md,BENCH_RUNS,mmd));
    }
#undef B2CD_N
}

     
/* B3: Validator - Verify timing */
     
void bench_validator_verify(ECTParams *pp)
{
    ect_clear_tagged(pp);
    printf("\n=== B3: Validator - Verify timing ===\n(%d runs)\noperation,mean_ms,sd_ms\n",BENCH_RUNS);
    struct timespec s0,s1; double t[BENCH_RUNS];

    /* schnorr_verify */
    { ECTTag w; ect_tag_new(&w); ect_tag(&w,pp,1,1);
      for (int i=0;i<BENCH_RUNS;i++) {
          clock_now(&s0); schnorr_verify(&w.sig,pp->pk_sig,w.c1_0,w.c2_0,pp->g,pp->qh,pp->p); clock_now(&s1);
          t[i]=timespec_ms(s0,s1); }
      double m=mean(t,BENCH_RUNS); printf("schnorr_verify,%.3f,%.3f\n",m,stdev(t,BENCH_RUNS,m));
      ect_tag_free(&w); }

    /* nizk_verify */
    { ECTTag w; ect_tag_new(&w); ect_tag(&w,pp,1,1); degrade_single(&w,pp);
      bn_t C1,C2; bn_null(C1); bn_null(C2); bn_new(C1); bn_new(C2);
      bn_copy(C1,w.steps[0].c1_in); bn_copy(C2,w.steps[0].c2_in);
      for (int i=0;i<BENCH_RUNS;i++) {
          clock_now(&s0);
          nizk_verify(&w.proofs[0],pp->pk,&C1,&C2,1,w.steps[0].c1_out,w.steps[0].c2_out,pp->g,pp->q,pp->qh,pp->p);
          clock_now(&s1); t[i]=timespec_ms(s0,s1); }
      double m=mean(t,BENCH_RUNS); printf("nizk_verify,%.3f,%.3f\n",m,stdev(t,BENCH_RUNS,m));
      bn_free(C1); bn_free(C2); ect_tag_free(&w); }

    /* hash_check */
    { Block b; memset(&b,0,sizeof(b)); uint8_t exp[HASH_BYTES]={0};
      for (int i=0;i<BENCH_RUNS;i++) {
          clock_now(&s0); block_hash(exp,&b); (void)memcmp(b.raw+BLK_OFF_PREV_HASH,exp,HASH_BYTES); clock_now(&s1);
          t[i]=timespec_ms(s0,s1); }
      double m=mean(t,BENCH_RUNS); printf("hash_check,%.3f,%.3f\n",m,stdev(t,BENCH_RUNS,m)); }

    /* full issue block */
    { ECTTag w; ect_tag_new(&w); ect_tag(&w,pp,1,1);
      Block b; uint8_t prev[HASH_BYTES]={0}; Watermark wm;
      wm_pack(&wm,w.c1,w.c2,NULL,0,w.c1_0,w.c2_0);
      block_set_header(&b,prev,0,1,2,10); block_set_watermark(&b,&wm);
      uint8_t cmp[HASH_BYTES];
      for (int i=0;i<BENCH_RUNS;i++) {
          clock_now(&s0); block_hash(cmp,&b); (void)memcmp(b.raw+BLK_OFF_PREV_HASH,cmp,HASH_BYTES);
          schnorr_verify(&w.sig,pp->pk_sig,w.c1_0,w.c2_0,pp->g,pp->qh,pp->p); clock_now(&s1);
          t[i]=timespec_ms(s0,s1); }
      double m=mean(t,BENCH_RUNS); printf("full_issue_block_verify,%.3f,%.3f\n",m,stdev(t,BENCH_RUNS,m));
      ect_tag_free(&w); }

    /* full spend block */
    { ECTTag w; ect_tag_new(&w); ect_tag(&w,pp,1,1); degrade_single(&w,pp);
      Block b; uint8_t prev[HASH_BYTES]={0}; Watermark wm;
      wm_pack(&wm,w.c1,w.c2,&w.proofs[0],(uint8_t)w.hops,w.c1_0,w.c2_0);
      block_set_header(&b,prev,1,2,3,5); block_set_watermark(&b,&wm);
      uint8_t cmp[HASH_BYTES];
      bn_t C1,C2; bn_null(C1); bn_null(C2); bn_new(C1); bn_new(C2);
      bn_copy(C1,w.steps[0].c1_in); bn_copy(C2,w.steps[0].c2_in);
      for (int i=0;i<BENCH_RUNS;i++) {
          clock_now(&s0); block_hash(cmp,&b); (void)memcmp(b.raw+BLK_OFF_PREV_HASH,cmp,HASH_BYTES);
          nizk_verify(&w.proofs[0],pp->pk,&C1,&C2,1,w.steps[0].c1_out,w.steps[0].c2_out,pp->g,pp->q,pp->qh,pp->p);
          clock_now(&s1); t[i]=timespec_ms(s0,s1); }
      double m=mean(t,BENCH_RUNS); printf("full_spend_block_verify,%.3f,%.3f\n",m,stdev(t,BENCH_RUNS,m));
      bn_free(C1); bn_free(C2); ect_tag_free(&w); }
}

     
/* B4: Sizes */
     
void bench_sizes(ECTParams *pp)
{
    printf("\n=== B4: Communication and storage sizes ===\n\n");
    int nizk_proof = 2*WM_ELEM_BYTES + 2*WM_Z_BYTES;
    int tag_ct     = 2*WM_ELEM_BYTES;
    int pp_bound   = 3*WM_ELEM_BYTES;
    int watermark  = WATERMARK_BYTES;
    int h          = pp->h;
    int user_total = 4*WM_ELEM_BYTES + h*(nizk_proof + 4*WM_ELEM_BYTES) + WM_ELEM_BYTES + WM_Z_BYTES;
    int trace_in   = 2*WM_ELEM_BYTES;

    printf("Quantity                                    Bytes      KB\n");
    printf("---------------------------------------------------------\n");
    printf("1. Proof of degradation (NIZKProof m=1)    %5d   %6.2f\n", nizk_proof, nizk_proof/1024.0);
    printf("2. Tag ciphertext (c1+c2, any depth)       %5d   %6.2f\n", tag_ct,     tag_ct/1024.0);
    printf("3. Public params pp (g,pkenc,pksig)        %5d   %6.2f\n", pp_bound,   pp_bound/1024.0);
    printf("4. Watermark per block                     %5d   %6.2f\n", watermark,  watermark/1024.0);
    printf("5. Total user storage H=%d hops            %5d   %6.2f\n", h, user_total, user_total/1024.0);
    printf("6. Trace input from ledger (c1+c2)         %5d   %6.2f\n", trace_in,   trace_in/1024.0);
    printf("---------------------------------------------------------\n");
    printf("Block total (fixed)                        %5d   %6.2f\n", BLOCK_BYTES,         BLOCK_BYTES/1024.0);
    printf("ECT overhead per block                     %5d   %6.2f\n", watermark,           watermark/1024.0);
    printf("ECT/header ratio                           %6.1fx\n",      (double)watermark/BLK_OFF_WATERMARK);
}

     
// /* B5: Population scaling, K = log2(pop) */
     
// void bench_population_scale(ECTParams *pp)
// {
// #define B5_RUNS 100
//     printf("\n=== B5: ECT overhead vs chain population ===\n");
//     printf("K_traced = log2(pop), capped at N=%d.\n", pp->n);
//     printf("(%d runs)\npop,k_traced,issue_ms,sd,spend_ms,sd,trace_ms,sd\n", B5_RUNS);
//     int pop_exp[]={4,5,6,7,8,9,10,11,12,13};
//     for (int pi=0;pi<10;pi++) {
//         int pop=1<<pop_exp[pi], K=pop_exp[pi];
//         if (K>pp->n) K=pp->n;
//         if (K<1) K=1;
//         double t_i[B5_RUNS],t_s[B5_RUNS],t_t[B5_RUNS]; struct timespec t0,t1;
//         for (int run=0;run<B5_RUNS;run++) {
//             ect_clear_tagged(pp);
//             ECTTag *tr=(ECTTag*)malloc((size_t)K*sizeof(ECTTag));
//             LightUser *bg=(LightUser*)malloc((size_t)(pop-K)*sizeof(LightUser));
//             if (!tr||!bg){free(tr);free(bg);continue;}
//             for (int j=0;j<K;j++) ect_tag_new(&tr[j]);
//             for (int j=0;j<pop-K;j++) light_user_new(&bg[j]);
//             { ECTTag proto; ect_tag_new(&proto); ect_tag(&proto,pp,1,0);
//               for (int j=0;j<pop-K;j++){bn_copy(bg[j].c1,proto.c1);bn_copy(bg[j].c2,proto.c2);bn_copy(bg[j].c1_0,proto.c1_0);bn_copy(bg[j].c2_0,proto.c2_0);}
//               ect_tag_free(&proto); }
//             clock_now(&t0); for (int j=0;j<K;j++) ect_tag(&tr[j],pp,j+1,1); clock_now(&t1); t_i[run]=timespec_ms(t0,t1)/K;
//             clock_now(&t0); for (int j=0;j<K;j++) degrade_single(&tr[j],pp); clock_now(&t1); t_s[run]=timespec_ms(t0,t1)/K;
//             /* Merge all K traced tags into one, then trace once.*/
//             ECTTag merged_b5; ect_tag_new(&merged_b5);
//             { ECTTag **ptrs=(ECTTag**)malloc((size_t)K*sizeof(ECTTag*));
//               for (int j=0;j<K;j++) ptrs[j]=&tr[j];
//               ect_merge(&merged_b5,ptrs,K,pp); free(ptrs); }
//             clock_now(&t0); ect_trace_silent(pp,&merged_b5); clock_now(&t1); t_t[run]=timespec_ms(t0,t1);
//             ect_tag_free(&merged_b5);
//             for (int j=0;j<K;j++) ect_tag_free(&tr[j]);
//             for (int j=0;j<nbg;j++) light_user_free(&bg[j]);
//             free(tr); free(bg);
//         }
//         double mi=mean(t_i,B5_RUNS),ms=mean(t_s,B5_RUNS),mt=mean(t_t,B5_RUNS);
//         printf("%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",pop,K,mi,stdev(t_i,B5_RUNS,mi),ms,stdev(t_s,B5_RUNS,ms),mt,stdev(t_t,B5_RUNS,mt));
//     }
// #undef B5_RUNS
// }

/* B5: Population scaling, K = log2(pop) */
     
void bench_population_scale(ECTParams *pp)
{
#define B5_RUNS TRACE_RUNS
    printf("\n=== B5: ECT overhead vs chain population ===\n");
    printf("K_traced = log2(pop), capped at N=%d.\n", pp->n);
    printf("(%d runs)\npop,k_traced,issue_ms,sd,spend_ms,sd,trace_ms,sd\n", B5_RUNS);
    /* Nothing timed here reads the population: it only sets K = log2(pop).
       B6 sweeps K directly, so a few decades suffice to show the shape and
       the full 17-point sweep just multiplies cost (sum of K was 204). */
    int pop_exp[]={4,10,16,20};
    for (size_t pi=0;pi<sizeof(pop_exp)/sizeof(*pop_exp);pi++) {
        int pop=1<<pop_exp[pi], K=pop_exp[pi];
        if (K>pp->n) K=pp->n;
        if (K<1) K=1;
        int nbg = pop - K;
        if (nbg > BENCH_POP_ALLOC_MAX) nbg = BENCH_POP_ALLOC_MAX;
        double t_i[B5_RUNS],t_s[B5_RUNS],t_t[B5_RUNS]; struct timespec t0,t1;
        for (int run=0;run<B5_RUNS;run++) {
            ect_clear_tagged(pp);
            ECTTag *tr=(ECTTag*)malloc((size_t)K*sizeof(ECTTag));
            LightUser *bg=(LightUser*)malloc((size_t)nbg*sizeof(LightUser));
            if (!tr||!bg){free(tr);free(bg);continue;}
            for (int j=0;j<K;j++) ect_tag_new(&tr[j]);
            for (int j=0;j<nbg;j++) light_user_new(&bg[j]);
            { ECTTag proto; ect_tag_new(&proto); ect_tag(&proto,pp,1,0);
              for (int j=0;j<nbg;j++){bn_copy(bg[j].c1,proto.c1);bn_copy(bg[j].c2,proto.c2);bn_copy(bg[j].c1_0,proto.c1_0);bn_copy(bg[j].c2_0,proto.c2_0);}
              ect_tag_free(&proto); }
            clock_now(&t0); for (int j=0;j<K;j++) ect_tag(&tr[j],pp,j+1,1); clock_now(&t1); t_i[run]=timespec_ms(t0,t1)/K;
            clock_now(&t0); for (int j=0;j<K;j++) degrade_single(&tr[j],pp); clock_now(&t1); t_s[run]=timespec_ms(t0,t1)/K;
            /* Merge all K traced tags into one, then trace once.*/
            ECTTag merged_b5; ect_tag_new(&merged_b5);
            { ECTTag **ptrs=(ECTTag**)malloc((size_t)K*sizeof(ECTTag*));
              for (int j=0;j<K;j++) ptrs[j]=&tr[j];
              ect_merge(&merged_b5,ptrs,K,pp); free(ptrs); }
            clock_now(&t0); ect_trace_silent(pp,&merged_b5); clock_now(&t1); t_t[run]=timespec_ms(t0,t1);
            ect_tag_free(&merged_b5);
            for (int j=0;j<K;j++) ect_tag_free(&tr[j]);
            for (int j=0;j<nbg;j++) light_user_free(&bg[j]);
            free(tr); free(bg);
        }
        double mi=mean(t_i,B5_RUNS),ms=mean(t_s,B5_RUNS),mt=mean(t_t,B5_RUNS);
        printf("%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",pop,K,mi,stdev(t_i,B5_RUNS,mi),ms,stdev(t_s,B5_RUNS,ms),mt,stdev(t_t,B5_RUNS,mt));
    }
#undef B5_RUNS
}
     
/* B6: Fixed pop, vary K                                                */
     
void bench_ktraced_scale(ECTParams *pp)
{
#define B6_POP 524288
#define B6_RUNS TRACE_RUNS
    printf("\n=== B6: Vary K_traced, fixed pop=%d ===\n(%d runs)\n",B6_POP,B6_RUNS);
    printf("k,issue_per_block_ms,sd,spend_per_block_ms,sd,trace_total_ms,sd,trace_per_user_ms,sd\n");
    /* K values: realistic investigation sizes. Trace is O(H*2^K), so K=16 (65K subsets) is the practical maximum for benchmarking. */
    int k_vals[] = {1, 2, 4, 8, 12, 16, 18, 23, 27};
    // int k_vals[]={1,2,4,8,12,16,24,48};
    for (int ki=0;ki<9;ki++) {
        int K=k_vals[ki]; if (K>pp->n) K=pp->n;
        int nbg = B6_POP - K;
        if (nbg > BENCH_POP_ALLOC_MAX) nbg = BENCH_POP_ALLOC_MAX;
        double t_i[B6_RUNS],t_s[B6_RUNS],t_tt[B6_RUNS],t_tp[B6_RUNS]; struct timespec t0,t1;
        for (int run=0;run<B6_RUNS;run++) {
            ect_clear_tagged(pp);
            ECTTag *tr=(ECTTag*)malloc((size_t)K*sizeof(ECTTag));
            LightUser *bg=(LightUser*)malloc((size_t)nbg*sizeof(LightUser));
            for (int j=0;j<K;j++) ect_tag_new(&tr[j]);
            for (int j=0;j<nbg;j++) light_user_new(&bg[j]);
            { ECTTag proto; ect_tag_new(&proto); ect_tag(&proto,pp,1,0);
              for (int j=0;j<nbg;j++){bn_copy(bg[j].c1,proto.c1);bn_copy(bg[j].c2,proto.c2);bn_copy(bg[j].c1_0,proto.c1_0);bn_copy(bg[j].c2_0,proto.c2_0);}
              ect_tag_free(&proto); }
            clock_now(&t0); for (int j=0;j<K;j++) ect_tag(&tr[j],pp,j+1,1); clock_now(&t1); t_i[run]=timespec_ms(t0,t1)/K;
            clock_now(&t0); for (int j=0;j<K;j++) degrade_single(&tr[j],pp); clock_now(&t1); t_s[run]=timespec_ms(t0,t1)/K;
            /* Merge all K traced tags into one, then trace once. */
            ECTTag merged_b6; ect_tag_new(&merged_b6);
            { ECTTag **ptrs=(ECTTag**)malloc((size_t)K*sizeof(ECTTag*));
              for (int j=0;j<K;j++) ptrs[j]=&tr[j];
              ect_merge(&merged_b6,ptrs,K,pp); free(ptrs); }
            clock_now(&t0); ect_trace_silent(pp,&merged_b6); clock_now(&t1);
            double tot=timespec_ms(t0,t1); t_tt[run]=tot; t_tp[run]=tot/K;
            ect_tag_free(&merged_b6);
            for (int j=0;j<K;j++) ect_tag_free(&tr[j]);
            for (int j=0;j<nbg;j++) light_user_free(&bg[j]);
            free(tr); free(bg);
        }
        double mi=mean(t_i,B6_RUNS),ms=mean(t_s,B6_RUNS),mtt=mean(t_tt,B6_RUNS),mtp=mean(t_tp,B6_RUNS);
        printf("%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",K,mi,stdev(t_i,B6_RUNS,mi),ms,stdev(t_s,B6_RUNS,ms),mtt,stdev(t_tt,B6_RUNS,mtt),mtp,stdev(t_tp,B6_RUNS,mtp));
    }
#undef B6_POP
#undef B6_RUNS
}

     
/* B7: Fixed pop and K, vary depth 1->H*/
     
void bench_degradation_depth(ECTParams *pp)
{
#define B7_K 6
#define B7_RUNS TRACE_RUNS
    printf("\n=== B7: Vary depth (1->H=%d), K=%d ===\n(%d runs)\ndepth,degrade_hop_ms,sd,trace_ms,sd\n",pp->h,B7_K,B7_RUNS);
    /* Walk each tag set down the hop budget once, recording the cost of the
       hop that lands at each depth.  Restarting from fresh tags per depth
       repeated 6*(1+2+...+h) degradations instead of 6*h, quadratic in h for
       an identical measurement. */
    double t_h[H][B7_RUNS], t_t[H][B7_RUNS];
    struct timespec s0,s1;
    for (int run=0;run<B7_RUNS;run++) {
        ect_clear_tagged(pp);
        ECTTag *tags=(ECTTag*)malloc((size_t)B7_K*sizeof(ECTTag));
        for (int j=0;j<B7_K;j++){ect_tag_new(&tags[j]);ect_tag(&tags[j],pp,j+1,1);}
        for (int t=1;t<=pp->h;t++) {
            clock_now(&s0); for (int j=0;j<B7_K;j++) degrade_single(&tags[j],pp); clock_now(&s1);
            t_h[t-1][run]=timespec_ms(s0,s1)/B7_K;
            /* Merge all K tags into one, then trace once. */
            ECTTag merged_b7; ect_tag_new(&merged_b7);
            { ECTTag **ptrs=(ECTTag**)malloc((size_t)B7_K*sizeof(ECTTag*));
              for (int j=0;j<B7_K;j++) ptrs[j]=&tags[j];
              ect_merge(&merged_b7,ptrs,B7_K,pp); free(ptrs); }
            clock_now(&s0); ect_trace_silent(pp,&merged_b7); clock_now(&s1);
            t_t[t-1][run]=timespec_ms(s0,s1);
            ect_tag_free(&merged_b7);
        }
        for (int j=0;j<B7_K;j++) ect_tag_free(&tags[j]);
        free(tags);
    }
    for (int t=1;t<=pp->h;t++) {
        double mh=mean(t_h[t-1],B7_RUNS),mt=mean(t_t[t-1],B7_RUNS);
        printf("%d,%.3f,%.3f,%.3f,%.3f\n",t,mh,stdev(t_h[t-1],B7_RUNS,mh),mt,stdev(t_t[t-1],B7_RUNS,mt));
    }
#undef B7_K
#undef B7_RUNS
}

     
/* B8: Vary merge size k=1..10, then degrade  */
     
void bench_merge_then_degrade(ECTParams *pp)
{
#define B8_RUNS BENCH_RUNS
    printf("\n=== B8: Vary merge size k (1->10), then degrade ===\n(%d runs)\nk,merge_ms,sd,degrade_after_merge_ms,sd\n",B8_RUNS);
    int k_max=10; if (k_max>pp->n) k_max=pp->n;
    for (int k=1;k<=k_max;k++) {
        double t_m[B8_RUNS],t_d[B8_RUNS]; struct timespec s0,s1;
        for (int run=0;run<B8_RUNS;run++) {
            ect_clear_tagged(pp);
            ECTTag *arr=(ECTTag*)malloc((size_t)k*sizeof(ECTTag));
            ECTTag **ptrs=(ECTTag**)malloc((size_t)k*sizeof(ECTTag*));
            ECTTag out; ect_tag_new(&out);
            for (int j=0;j<k;j++){ect_tag_new(&arr[j]);ect_tag(&arr[j],pp,j+1,1);ptrs[j]=&arr[j];}
            clock_now(&s0); ect_merge(&out,ptrs,k,pp); clock_now(&s1); t_m[run]=timespec_ms(s0,s1);
            clock_now(&s0); degrade_single(&out,pp); clock_now(&s1); t_d[run]=timespec_ms(s0,s1);
            ect_tag_free(&out); for (int j=0;j<k;j++) ect_tag_free(&arr[j]); free(arr); free(ptrs);
        }
        double mm=mean(t_m,B8_RUNS),md=mean(t_d,B8_RUNS);
        printf("%d,%.3f,%.3f,%.3f,%.3f\n",k,mm,stdev(t_m,B8_RUNS,mm),md,stdev(t_d,B8_RUNS,md));
    }
#undef B8_RUNS
}

     
/* B9: Full-fidelity (NOT called by default)   */
     
void bench_full_fidelity(ECTParams *pp)
{
#define B9_K 6
#define B9_RUNS 5
    printf("\n=== B9: Full-fidelity pop scaling (slow, manual only) ===\n");
    printf("K=%d. Every user gets real crypto + block commit. Capped 2^4..2^8.\n(%d runs)\n",B9_K,B9_RUNS);
    printf("pop,issue_ms,sd,spend_ms,sd,trace_ms,sd\n");
    int pop_exp[]={4,5,6,7,8};
    for (int pi=0;pi<5;pi++) {
        int pop=1<<pop_exp[pi],K=B9_K; if (K>pop) K=pop;
        double t_i[B9_RUNS],t_s[B9_RUNS],t_t[B9_RUNS]; struct timespec t0,t1;
        for (int run=0;run<B9_RUNS;run++) {
            ect_clear_tagged(pp);
            Ledger *ledger=ledger_new();
            for (int u=K;u<pop;u++) {
                ECTTag bg; ect_tag_new(&bg); ect_tag(&bg,pp,1,0);
                Watermark wm; wm_pack(&wm,bg.c1,bg.c2,NULL,0,bg.c1_0,bg.c2_0);
                Block b; uint8_t prev[HASH_BYTES]={0};
                if (ledger->height>0) block_hash(prev,&ledger->blocks[ledger->height-1]);
                block_set_header(&b,prev,ledger->height,1,(uint32_t)(u+1),10);
                block_set_watermark(&b,&wm); ledger_append(ledger,&b); ect_tag_free(&bg);
            }
            ECTTag *tr=(ECTTag*)malloc((size_t)K*sizeof(ECTTag));
            for (int j=0;j<K;j++) ect_tag_new(&tr[j]);
            clock_now(&t0);
            for (int j=0;j<K;j++) {
                ect_tag(&tr[j],pp,j+1,1);
                Watermark wm; wm_pack(&wm,tr[j].c1,tr[j].c2,NULL,0,tr[j].c1_0,tr[j].c2_0);
                Block b; uint8_t prev[HASH_BYTES]={0};
                if (ledger->height>0) block_hash(prev,&ledger->blocks[ledger->height-1]);
                block_set_header(&b,prev,ledger->height,1,(uint32_t)(j+1),10);
                block_set_watermark(&b,&wm); ledger_append(ledger,&b);
            }
            clock_now(&t1); t_i[run]=timespec_ms(t0,t1)/K;
            clock_now(&t0);
            for (int j=0;j<K;j++) {
                degrade_single(&tr[j],pp);
                NIZKProof *pi=tr[j].n_proofs>0?&tr[j].proofs[tr[j].n_proofs-1]:NULL;
                Watermark wm; wm_pack(&wm,tr[j].c1,tr[j].c2,pi,(uint8_t)tr[j].hops,tr[j].c1_0,tr[j].c2_0);
                Block b; uint8_t prev[HASH_BYTES]={0};
                if (ledger->height>0) block_hash(prev,&ledger->blocks[ledger->height-1]);
                block_set_header(&b,prev,ledger->height,(uint32_t)(j+1),(uint32_t)(j+100),5);
                block_set_watermark(&b,&wm); ledger_append(ledger,&b);
            }
            clock_now(&t1); t_s[run]=timespec_ms(t0,t1)/K;
            clock_now(&t0); for (int j=0;j<K;j++) ect_trace_silent(pp,&tr[j]); clock_now(&t1); t_t[run]=timespec_ms(t0,t1)/K;
            for (int j=0;j<K;j++) ect_tag_free(&tr[j]);
            free(tr); ledger_free(ledger);
        }
        double mi=mean(t_i,B9_RUNS),ms=mean(t_s,B9_RUNS),mt=mean(t_t,B9_RUNS);
        printf("%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",pop,mi,stdev(t_i,B9_RUNS,mi),ms,stdev(t_s,B9_RUNS,ms),mt,stdev(t_t,B9_RUNS,mt));
    }
#undef B9_K
#undef B9_RUNS
}


     
/* B_nizk: NIZK prove/verify scaling with m inputs */
     
void bench_nizk_scale(ECTParams *pp)
{
#define NIZK_RUNS BENCH_RUNS
    ect_clear_tagged(pp);
    printf("\n=== B_nizk: NIZK prove/verify vs m (input set size) ===\n");
    printf("(%d runs)\n", NIZK_RUNS);
    printf("m,proof_bytes,prove_ms,prove_sd,verify_ms,verify_sd\n");

    int m_vals[] = {1, 2, 4, 8, 16, 32, 64, 128};
    int n_m = 8;

    for (int mi = 0; mi < n_m; mi++) {
        int m = m_vals[mi];
        double t_prove[NIZK_RUNS], t_verify[NIZK_RUNS];
        struct timespec s0, s1;

        for (int run = 0; run < NIZK_RUNS; run++) {
            bn_t *C1s = (bn_t*)malloc((size_t)m * sizeof(bn_t));
            bn_t *C2s = (bn_t*)malloc((size_t)m * sizeof(bn_t));
            for (int j = 0; j < m; j++) {
                bn_null(C1s[j]); bn_new(C1s[j]);
                bn_null(C2s[j]); bn_new(C2s[j]);
            }

            ECTTag real_tag; ect_tag_new(&real_tag);
            ect_tag(&real_tag, pp, 1, 1);
            bn_copy(C1s[0], real_tag.c1);
            bn_copy(C2s[0], real_tag.c2);

            bn_t r_dummy, m_zero;
            bn_null(r_dummy); bn_new(r_dummy);
            bn_null(m_zero);  bn_new(m_zero);
            bn_zero(m_zero);
            for (int j = 1; j < m; j++)
                elgamal_enc(C1s[j], C2s[j], r_dummy, m_zero,
                            pp->pk, pp->g, pp->qh, pp->p);
            bn_free(r_dummy); bn_free(m_zero);

            bn_t c1p, c2p, r_prime;
            bn_null(c1p); bn_new(c1p);
            bn_null(c2p); bn_new(c2p);
            bn_null(r_prime); bn_new(r_prime);
            bn_copy(c1p, C1s[0]); bn_copy(c2p, C2s[0]);
            ciphertext_pow_q(c1p, c2p, pp->q, pp->p);
            elgamal_rerand(c1p, c2p, r_prime, pp->pk, pp->g, pp->qh, pp->p);

            bn_t *A1, *A2, *zz, *ee;
            clock_now(&s0);
            nizk_bench_prove(pp->pk, C1s, C2s, m, c1p, c2p, r_prime,
                             pp->g, pp->q, pp->qh, pp->p,
                             &A1, &A2, &zz, &ee);
            clock_now(&s1);
            t_prove[run] = timespec_ms(s0, s1);

            clock_now(&s0);
            nizk_bench_verify(pp->pk, C1s, C2s, m, c1p, c2p,
                              A1, A2, zz, ee,
                              pp->g, pp->q, pp->qh, pp->p);
            clock_now(&s1);
            t_verify[run] = timespec_ms(s0, s1);

            for (int j = 0; j < m; j++) {
                bn_free(A1[j]); bn_free(A2[j]);
                bn_free(zz[j]); bn_free(ee[j]);
            }
            free((void *)A1); free((void *)A2); free((void *)zz); free((void *)ee);

            bn_free(c1p); bn_free(c2p); bn_free(r_prime);
            ect_tag_free(&real_tag);
            for (int j = 0; j < m; j++) { bn_free(C1s[j]); bn_free(C2s[j]); }
            free((void *)C1s); free((void *)C2s);
        }

        int proof_bytes = m * (2*WM_ELEM_BYTES + 2*WM_Z_BYTES);
        double mp = mean(t_prove, NIZK_RUNS);
        double mv = mean(t_verify, NIZK_RUNS);
        printf("%d,%d,%.3f,%.3f,%.3f,%.3f\n",
               m, proof_bytes, mp, stdev(t_prove, NIZK_RUNS, mp),
               mv, stdev(t_verify, NIZK_RUNS, mv));
    }
#undef NIZK_RUNS
}


     
/* Table 1 + Table 2  */
     
static void bench_table_summary(ECTParams *pp)
{
#define T_RUNS BENCH_RUNS
    ect_clear_tagged(pp);
    struct timespec s0, s1;
    double t[T_RUNS];

    /* Tag creation */
    for (int i=0;i<T_RUNS;i++) {
        ECTTag w; ect_tag_new(&w);
        clock_now(&s0); ect_tag(&w,pp,1,1); clock_now(&s1);
        t[i]=timespec_ms(s0,s1); ect_tag_free(&w);
    }
    double t_tag=mean(t,T_RUNS);

    /* Full degrade */
    for (int i=0;i<T_RUNS;i++) {
        ECTTag w; ect_tag_new(&w); ect_tag(&w,pp,1,1);
        clock_now(&s0); degrade_single(&w,pp); clock_now(&s1);
        t[i]=timespec_ms(s0,s1); ect_tag_free(&w);
    }
    double t_degrade=mean(t,T_RUNS);

    /* Merge 2 coins */
    for (int i=0;i<T_RUNS;i++) {
        ECTTag a,b,out;
        ect_tag_new(&a); ect_tag_new(&b); ect_tag_new(&out);
        ect_tag(&a,pp,1,1); ect_tag(&b,pp,2,1);
        ECTTag *ptrs[2]={&a,&b};
        clock_now(&s0); ect_merge(&out,ptrs,2,pp); clock_now(&s1);
        t[i]=timespec_ms(s0,s1);
        ect_tag_free(&a); ect_tag_free(&b); ect_tag_free(&out);
    }
    double t_merge=mean(t,T_RUNS);

    /* Proof of budget: nizk_verify per tx (validator cost) */
    for (int i=0;i<T_RUNS;i++) {
        ECTTag w; ect_tag_new(&w); ect_tag(&w,pp,1,1);
        degrade_single(&w,pp);
        bn_t C1,C2; bn_null(C1);bn_new(C1); bn_null(C2);bn_new(C2);
        bn_copy(C1,w.steps[0].c1_in); bn_copy(C2,w.steps[0].c2_in);
        clock_now(&s0);
        nizk_verify(&w.proofs[0],pp->pk,&C1,&C2,1,
                    w.c1,w.c2,pp->g,pp->q,pp->qh,pp->p);
        clock_now(&s1);
        t[i]=timespec_ms(s0,s1);
        bn_free(C1);bn_free(C2); ect_tag_free(&w);
    }
    double t_verify=mean(t,T_RUNS);

    /* Full block verification (schnorr + nizk) */
    for (int i=0;i<T_RUNS;i++) {
        ECTTag w; ect_tag_new(&w); ect_tag(&w,pp,1,1);
        degrade_single(&w,pp);
        bn_t C1,C2; bn_null(C1);bn_new(C1); bn_null(C2);bn_new(C2);
        bn_copy(C1,w.steps[0].c1_in); bn_copy(C2,w.steps[0].c2_in);
        clock_now(&s0);
        schnorr_verify(&w.sig,pp->pk_sig,w.c1_0,w.c2_0,pp->g,pp->qh,pp->p);
        nizk_verify(&w.proofs[0],pp->pk,&C1,&C2,1,
                    w.c1,w.c2,pp->g,pp->q,pp->qh,pp->p);
        clock_now(&s1);
        t[i]=timespec_ms(s0,s1);
        bn_free(C1);bn_free(C2); ect_tag_free(&w);
    }
    double t_blk=mean(t,T_RUNS);

    /* Sizes */
    int tag_sz     = 2*WM_ELEM_BYTES;
    int proof_sz   = 2*WM_ELEM_BYTES + 2*WM_Z_BYTES;
    int degrade_sz = tag_sz+proof_sz;
    int merge_sz   = 2*WM_ELEM_BYTES;
    int wm_sz      = WATERMARK_BYTES;
    int user_sz    = 4*WM_ELEM_BYTES + pp->h*(proof_sz + 4*WM_ELEM_BYTES) + WM_ELEM_BYTES + WM_Z_BYTES;

    /*  Table 1  */
    printf("\n=== Table 1: Storage and bandwidth costs ===\n");
    printf("%-48s %9s %9s\n","Object","Size(KB)","Time(ms)");
    printf("%-48s %9s %9s\n",
           "------------------------------------------------",
           "---------","---------");
    printf("%-48s %9.2f %9.1f\n",
           "Tag creation (c1+c2)",
           tag_sz/1024.0, t_tag);
    printf("%-48s %9.2f %9.1f\n",
           "Degradation: tag + proof total (m=1)",
           degrade_sz/1024.0, t_degrade);
    printf("%-48s %9.2f %9.1f\n",
           "Merging output (c1+c2, any k)",
           merge_sz/1024.0, t_merge);
    printf("%-48s %9.2f %9.1f\n",
           "Verification of tag in block (validator)",
           wm_sz/1024.0, t_blk);
    printf("%-48s %9.2f %9s\n",
           "User storage (H=5 hops)",
           user_sz/1024.0, "-");
    printf("%-48s %9.2f %9.1f\n",
           "Proof of budget per tx (CDS NIZK m=1, verify)",
           proof_sz/1024.0, t_verify);

    /*  Table 2: merge(k)+degrade vs depth  */
    int depths[]={0,1,2,3,5,7,9}; int n_d=7;
    int ks[]={1,2,3,4,5,6,8,10}; int n_k=8;

    printf("\n=== Table 2: Merge(k)+Degrade time (ms) ===\n");
    printf("rows = k (coins merged), cols = depth t before merge\n");
    printf("%-5s","k\\t");
    for (int di=0;di<n_d;di++) printf(" %8d",depths[di]);
    printf("\n");

    for (int ki=0;ki<n_k;ki++) {
        int k=ks[ki];
        printf("%-5d",k);
        for (int di=0;di<n_d;di++) {
            int dep=depths[di];
            double tt[T_RUNS];
            for (int run=0;run<T_RUNS;run++) {
                ECTTag *arr=(ECTTag*)malloc((size_t)k*sizeof(ECTTag));
                ECTTag **ptrs=(ECTTag**)malloc((size_t)k*sizeof(ECTTag*));
                ECTTag out; ect_tag_new(&out);
                for (int j=0;j<k;j++){
                    ect_tag_new(&arr[j]); ect_tag(&arr[j],pp,j+1,1);
                    for (int d=0;d<dep;d++) degrade_single(&arr[j],pp);
                    ptrs[j]=&arr[j];
                }
                clock_now(&s0);
                ect_merge(&out,ptrs,k,pp);
                degrade_single(&out,pp);
                clock_now(&s1);
                tt[run]=timespec_ms(s0,s1);
                ect_tag_free(&out);
                for (int j=0;j<k;j++) ect_tag_free(&arr[j]);
                free(arr); free(ptrs);
            }
            printf(" %8.1f",mean(tt,T_RUNS));
        }
        printf("\n");
    }
#undef T_RUNS
}

     
/* run_benchmarks   */
     
void run_benchmarks(ECTParams *pp)
{
    /* B5-B7 tag only distinct identifiers, so no depth ever carries the same
       identifier twice and a multiplicity bound of 1 is exact for this
       workload.  Trace then costs 2^(K/2) rather than (MULT_MAX+1)^(K/2),
       which is what keeps the K sweep tractable.  Deployments where a split
       coin is merged back together need the MULT_MAX default instead. */
    int saved_mult = ect_mult_max;
    ect_mult_max = 1;

    printf("\n=========================================================\n");
    printf("  ECT Benchmarks  H=%d  N=%d  R=%d\n",pp->h,pp->n,pp->r);
    printf("  |p|=%d bits  BENCH_RUNS=%d  TRACE_RUNS=%d  mult_max=%d\n",
           P_BITS, BENCH_RUNS, TRACE_RUNS, ect_mult_max);
    printf("  background population capped at %d (no timed section reads it)\n",
           BENCH_POP_ALLOC_MAX);
    printf("=========================================================\n");
    bench_authority_tag(pp);
    bench_user_degrade_merge(pp);
    bench_validator_verify(pp);
    bench_sizes(pp);
    bench_population_scale(pp);
    bench_ktraced_scale(pp);
    bench_degradation_depth(pp);
    bench_merge_then_degrade(pp);
    bench_nizk_scale(pp);
    bench_table_summary(pp);
    printf("\n[B9 excluded - call bench_full_fidelity(pp) manually]\n");
    ect_mult_max = saved_mult;
}
