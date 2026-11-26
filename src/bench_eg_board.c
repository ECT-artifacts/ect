#include "bench_util.h"
#include "eg_board.h"
#include "board_util.h"
#include "prg.h"
#include <pthread.h>
#include <string.h>

/* bench_eg_board.c -> exponential ElGamal board-proof benchmarks.
   Usage: ./bench_eg_board runs M threads [chunk] [w]
   runs: timed repetitions of one entry on one thread; M: entries of the
   full board; threads: threads of the full board; chunk: entries held in
   memory at a time (default 4096); w: window of the fixed-base tables
   (default 8).  The full board is proved and verified once, chunk by chunk:
   each chunk is proved, its leaves and products are recorded, and it is
   added to the verification batch of its thread, so that a board of 10^6
   entries does not have to fit in memory; the batches are checked at the
   end.
   The times of the chunks are summed into the time of Board and of
   VerifyBoard, together with the Merkle root, the signature and the count
   proof.  Output: CSV rows of bench_util.h; the full board is labelled with
   the number of threads. */

typedef struct {
    const EGBoardParams *bp;
    EGEntry *E;
    size_t base;          /* position of E[0] in the board */
    int *slot, *b;
    mpz_t *tau;
    uint8_t (*leaf)[32];
    EGBatch *bt;          /* one verification batch per thread, across chunks */
    int bad;
    pthread_mutex_t mu;
} Chunk;

typedef struct { const EGBoardParams *bp; EGEntry E; mpz_t tau; } One;

static void op_prove1(void *o)  { One *c = o; eg_entry_prove(c->bp, &c->E, 0, 5, 1, c->tau); }
static void op_verify1(void *o) { One *c = o; if (!eg_entry_verify(c->bp, &c->E, 0)) abort(); }

static void prove_range(void *ctx, int tid, size_t lo, size_t hi)
{
    (void)tid;
    Chunk *C = ctx;
    for (size_t i = lo; i < hi; i++) {
        size_t j = C->base + i;
        eg_entry_prove(C->bp, &C->E[i], j, C->slot[j], C->b[j], C->tau[i]);
        eg_entry_leaf(&C->E[i], C->leaf[j]);
    }
}

static void verify_range(void *ctx, int tid, size_t lo, size_t hi)
{
    Chunk *C = ctx;
    uint8_t leaf[32];
    for (size_t i = lo; i < hi; i++) {
        size_t j = C->base + i;
        eg_batch_add(C->bp, &C->bt[tid], &C->E[i], j);
        eg_entry_leaf(&C->E[i], leaf);
        if (memcmp(leaf, C->leaf[j], 32) != 0) {
            pthread_mutex_lock(&C->mu); C->bad = 1; pthread_mutex_unlock(&C->mu);
        }
    }
}

static void accumulate(const EGBoardParams *bp, const EGEntry *E, size_t k, mpz_t X, mpz_t Y)
{
    for (size_t i = 0; i < k; i++) {
        mpz_mul(X, X, E[i].C1); mpz_mod(X, X, bp->p);
        mpz_mul(Y, Y, E[i].C2); mpz_mod(Y, Y, bp->p);
    }
}

int main(int argc, char **argv)
{
    int runs = argc > 1 ? atoi(argv[1]) : 10;
    size_t M = argc > 2 ? (size_t)atol(argv[2]) : 100;
    int threads = argc > 3 ? atoi(argv[3]) : 1;
    size_t chunk = argc > 4 ? (size_t)atol(argv[4]) : 4096;
    int w = argc > 5 ? atoi(argv[5]) : 8;
    if (core_init() != RLC_OK) return 1;
    ECTParams pp;
    ect_params_new(&pp);
    ect_setup(&pp);
    EGBoardParams bp;
    eg_board_params_init(&bp, &pp, w);
    char label[64];
    snprintf(label, sizeof label, "eg_h%d_p%d_n%d", H, P_BITS, N);
    size_t eb = eg_entry_bytes(&bp);

    One one;
    one.bp = &bp;
    eg_entry_init(&bp, &one.E);
    mpz_init(one.tau);
    bench_measure(label, "board_entry", runs, eb, op_prove1, NULL, &one);
    bench_measure(label, "verify_board_entry", runs, eb, op_verify1, NULL, &one);
    eg_entry_free(&bp, &one.E);
    mpz_clear(one.tau);

    /* The full board: one entry in ten nonzero, identifiers in turn. */
    if (chunk > M) chunk = M;
    Chunk C;
    memset(&C, 0, sizeof C);
    C.bp = &bp;
    pthread_mutex_init(&C.mu, NULL);
    C.slot = malloc(M * sizeof(int));
    C.b = malloc(M * sizeof(int));
    C.leaf = malloc(M * 32);
    C.E = malloc(chunk * sizeof(EGEntry));
    C.tau = malloc(chunk * sizeof(mpz_t));
    int nt = threads < 1 ? 1 : threads;
    C.bt = malloc((size_t)nt * sizeof(EGBatch));
    if (!C.slot || !C.b || !C.leaf || !C.E || !C.tau || !C.bt) abort();
    for (int t = 0; t < nt; t++) eg_batch_init(&C.bt[t]);
    unsigned long c = 0;
    for (size_t j = 0; j < M; j++) {
        C.slot[j] = 1 + (int)(j % (size_t)bp.n);
        C.b[j] = (j % 10 == 0);
        c += (unsigned long)C.b[j];
    }
    for (size_t i = 0; i < chunk; i++) { eg_entry_init(&bp, &C.E[i]); mpz_init(C.tau[i]); }

    mpz_t Xp, Yp, Xv, Yv, tau, ce, cz;
    mpz_inits(Xp, Yp, Xv, Yv, tau, ce, cz, NULL);
    mpz_set_ui(Xp, 1); mpz_set_ui(Yp, 1); mpz_set_ui(Xv, 1); mpz_set_ui(Yv, 1);
    double tp = 0, tv = 0, t0;
    for (C.base = 0; C.base < M; C.base += chunk) {
        size_t k = M - C.base < chunk ? M - C.base : chunk;
        t0 = bench_now_ms();
        board_parallel(threads, k, prove_range, &C);
        accumulate(&bp, C.E, k, Xp, Yp);
        for (size_t i = 0; i < k; i++) mpz_add(tau, tau, C.tau[i]);
        mpz_mod(tau, tau, bp.qh);
        tp += bench_now_ms() - t0;
        t0 = bench_now_ms();
        board_parallel(threads, k, verify_range, &C);
        accumulate(&bp, C.E, k, Xv, Yv);
        tv += bench_now_ms() - t0;
        if (C.bad) { fprintf(stderr, "entry rejected\n"); abort(); }
        fprintf(stderr, "%zu/%zu entries, prove %.0f s, verify %.0f s\n",
                C.base + k, M, tp / 1000, tv / 1000);
    }

    uint8_t root[32], root2[32], msg[BOARD_MSG_BYTES], *sig = malloc(bp.sig.sig_bytes);
    t0 = bench_now_ms();
    merkle_root((const uint8_t (*)[32])C.leaf, M, root);
    board_message(msg, M, root);
    size_t siglen = sig_sign(&bp.sig, msg, sizeof msg, sig, bp.sig.sig_bytes);
    eg_count_prove(&bp, Xp, Yp, M, c, tau, ce, cz);
    tp += bench_now_ms() - t0;
    t0 = bench_now_ms();
    merkle_root((const uint8_t (*)[32])C.leaf, M, root2);
    board_message(msg, M, root2);
    for (int t = 1; t < nt; t++) eg_batch_merge(&bp, &C.bt[0], &C.bt[t]);
    int ok = eg_batch_check(&bp, &C.bt[0])
          && memcmp(root, root2, 32) == 0 && sig_verify(&bp.sig, msg, sizeof msg, sig, siglen)
          && eg_count_verify(&bp, Xv, Yv, M, c, ce, cz);
    tv += bench_now_ms() - t0;
    if (!ok) { fprintf(stderr, "board rejected\n"); abort(); }

    size_t bb = M * eb + eg_count_bytes(&bp) + eg_root_bytes(&bp);
    char name[96];
    snprintf(name, sizeof name, "board_M%zu_t%d", M, threads);
    printf("%s,%s,1,%.3f,0.000,%.3f,%zu\n", label, name, tp, tp, bb);
    snprintf(name, sizeof name, "verify_board_M%zu_t%d", M, threads);
    printf("%s,%s,1,%.3f,0.000,%.3f,%zu\n", label, name, tv, tv, bb);

    for (size_t i = 0; i < chunk; i++) { eg_entry_free(&bp, &C.E[i]); mpz_clear(C.tau[i]); }
    for (int t = 0; t < nt; t++) eg_batch_free(&C.bt[t]);
    free(C.slot); free(C.b); free(C.leaf); free(C.E); free(C.tau); free(C.bt); free(sig);
    mpz_clears(Xp, Yp, Xv, Yv, tau, ce, cz, NULL);
    eg_board_params_free(&bp);
    ect_params_free(&pp);
    core_clean();
    return 0;
}
