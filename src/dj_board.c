#include "dj_board.h"
#include "board_util.h"
#include "prg.h"
#include <openssl/evp.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

/* dj_board.c -> Damgard-Jurik board proofs. */

/* ---- hashing ------------------------------------------------------------- */

static void put_elem(EVP_MD_CTX *h, const DJParams *pp, const mpz_t x, int s)
{
    size_t len = dj_elem_len(pp, s);
    uint8_t *buf = malloc(len);
    if (!buf) abort();
    dj_write(pp, buf, x, s);
    EVP_DigestUpdate(h, buf, len);
    free(buf);
}

static void put_u64(EVP_MD_CTX *h, uint64_t v)
{
    uint8_t b[8];
    for (int i = 0; i < 8; i++) b[i] = (uint8_t)(v >> (56 - 8 * i));
    EVP_DigestUpdate(h, b, 8);
}

static EVP_MD_CTX *chal_init(const DJParams *pp, const char *label)
{
    EVP_MD_CTX *h = EVP_MD_CTX_new();
    if (!h || EVP_DigestInit_ex(h, EVP_sha256(), NULL) != 1) abort();
    EVP_DigestUpdate(h, label, strlen(label) + 1);
    put_elem(h, pp, pp->N, 0);
    put_u64(h, (uint64_t)pp->h);
    put_u64(h, (uint64_t)pp->n);
    return h;
}

static void chal_final(EVP_MD_CTX *h, mpz_t e)
{
    uint8_t d[32];
    unsigned int len = 0;
    EVP_DigestFinal_ex(h, d, &len);
    EVP_MD_CTX_free(h);
    mpz_import(e, DJ_CHAL_BITS / 8, 1, 1, 1, 0, d);
}

static void rand_chal(mpz_t e)
{
    uint8_t d[DJ_CHAL_BITS / 8];
    rnd_bytes(d, sizeof d);
    mpz_import(e, sizeof d, 1, 1, 1, 0, d);
}

/* ---- atoms --------------------------------------------------------------- */

/* The atom "x (1+N)^{-shift} is an N^s-th power mod N^{s+1}".  Its first
   message for challenge e and response z is
     a = z^{N^s} x^{-e} (1+N)^{shift e} mod N^{s+1},
   computed from xi = x^{-1}; shift is NULL for 0.  sk selects the
   computation of z^{N^s} with the factorization of N.  Returns 0 if z is
   not a unit. */
static int atom_first(const DJParams *pp, mpz_t a, const mpz_t xi, const mpz_t shift,
                      const mpz_t e, const mpz_t z, int s, int sk, mpz_t t)
{
    mpz_gcd(t, z, pp->N);
    if (mpz_sgn(z) <= 0 || mpz_cmp_ui(t, 1) != 0) return 0;
    if (sk) dj_pow_Ns_sk(pp, a, z, s);
    else    dj_pow_Ns(pp, a, z, s);
    mpz_powm(t, xi, e, pp->Np[s + 1]);
    mpz_mul(a, a, t);
    mpz_mod(a, a, pp->Np[s + 1]);
    if (shift) {
        /* dj_one_plus_N_pow overwrites its output before it reads m */
        mpz_t x;
        mpz_init(x);
        mpz_mul(x, shift, e);
        mpz_mod(x, x, pp->Np[s]);
        dj_one_plus_N_pow(pp, t, x, s);
        mpz_mul(a, a, t);
        mpz_mod(a, a, pp->Np[s + 1]);
        mpz_clear(x);
    }
    return 1;
}

/* z = k x^e mod N */
static void respond(const DJParams *pp, mpz_t z, const mpz_t k, const mpz_t e,
                    const mpz_t x, mpz_t t)
{
    mpz_powm(t, x, e, pp->N);
    mpz_mul(z, k, t);
    mpz_mod(z, z, pp->N);
}

static int is_unit(const DJParams *pp, const mpz_t x, int s, mpz_t t)
{
    if (mpz_sgn(x) <= 0 || mpz_cmp(x, pp->Np[s + 1]) >= 0) return 0;
    mpz_gcd(t, x, pp->N);
    return mpz_cmp_ui(t, 1) == 0;
}

static int in_range(const mpz_t x, const mpz_t bound)
{
    return mpz_sgn(x) >= 0 && mpz_cmp(x, bound) < 0;
}

/* ---- allocation ---------------------------------------------------------- */

void dj_entry_init(const DJParams *pp, DJEntry *E)
{
    mpz_inits(E->e, E->C, E->c0, NULL);
    E->c = malloc((size_t)pp->n * sizeof(mpz_t));
    E->z = malloc((size_t)DJ_NZ(pp->n) * sizeof(mpz_t));
    if (!E->c || !E->z) abort();
    for (int i = 0; i < pp->n; i++) mpz_init(E->c[i]);
    for (int i = 0; i < DJ_NZ(pp->n); i++) mpz_init(E->z[i]);
}

void dj_entry_free(const DJParams *pp, DJEntry *E)
{
    mpz_clears(E->e, E->C, E->c0, NULL);
    for (int i = 0; i < pp->n; i++) mpz_clear(E->c[i]);
    for (int i = 0; i < DJ_NZ(pp->n); i++) mpz_clear(E->z[i]);
    free(E->c); free(E->z);
}

void dj_board_init(DJBoard *B, const DJParams *pp, size_t M)
{
    B->M = M;
    B->E = malloc(M * sizeof(DJEntry));
    B->sig = malloc(pp->sig.sig_bytes);
    if (!B->E || !B->sig) abort();
    for (size_t j = 0; j < M; j++) dj_entry_init(pp, &B->E[j]);
    B->siglen = 0;
    B->c = 0;
    mpz_inits(B->ce, B->cz, NULL);
}

void dj_board_free(DJBoard *B, const DJParams *pp)
{
    for (size_t j = 0; j < B->M; j++) dj_entry_free(pp, &B->E[j]);
    free(B->E); free(B->sig);
    mpz_clears(B->ce, B->cz, NULL);
}

/* ---- the entry proof ----------------------------------------------------- */

/* Atom indices: 0 = e enc 0, 1 = C enc 0, 2 = C enc 1, 3 + i = e enc m_{i+1}. */

static void entry_challenge(const DJParams *pp, const DJEntry *E, uint64_t j,
                            mpz_t *a, mpz_t e)
{
    const int h = pp->h;
    EVP_MD_CTX *H = chal_init(pp, "ECT-DJ-board-entry");
    put_u64(H, j);
    put_elem(H, pp, E->e, h);
    put_elem(H, pp, E->C, 1);
    for (int k = 0; k < DJ_NZ(pp->n); k++)
        put_elem(H, pp, a[k], (k == 1 || k == 2) ? 1 : h);
    chal_final(H, e);
}

static mpz_t *new_vec(int k)
{
    mpz_t *v = malloc((size_t)k * sizeof(mpz_t));
    if (!v) abort();
    for (int i = 0; i < k; i++) mpz_init(v[i]);
    return v;
}

static void free_vec(mpz_t *v, int k)
{
    for (int i = 0; i < k; i++) mpz_clear(v[i]);
    free(v);
}

void dj_entry_prove(const DJParams *pp, DJEntry *E, uint64_t j, int slot, int b,
                    mpz_t tau_out)
{
    const int h = pp->h, n = pp->n;
    if (b && (slot < 1 || slot > n)) abort();
    int ia = b ? slot - 1 : -1;
    mpz_t *a = new_vec(DJ_NZ(n)), *m = new_vec(n);
    mpz_t rho, tau, k1, k2, e, c1, one, bm, ei, Ci, t, mask;
    mpz_inits(rho, tau, k1, k2, e, c1, one, bm, ei, Ci, t, mask, NULL);
    mpz_set_ui(one, 1);
    mpz_set_ui(mask, 0);
    mpz_setbit(mask, DJ_CHAL_BITS);
    for (int i = 0; i < n; i++) mpz_pow_ui(m[i], pp->r, (unsigned long)i);

    if (b) mpz_set(bm, m[ia]); else mpz_set_ui(bm, 0);
    dj_enc_sk(pp, E->e, bm, h, rho);
    mpz_set_ui(t, (unsigned long)b);
    dj_enc_sk(pp, E->C, t, 1, tau);
    if (!mpz_invert(ei, E->e, pp->Np[h + 1]) || !mpz_invert(Ci, E->C, pp->Np[2])) abort();

    if (!b) {
        mpz_set_ui(c1, 0);
        for (int i = 0; i < n; i++) {
            rand_chal(E->c[i]);
            mpz_add(c1, c1, E->c[i]);
            dj_rand_unit(pp, E->z[3 + i]);
            atom_first(pp, a[3 + i], ei, m[i], E->c[i], E->z[3 + i], h, 1, t);
        }
        mpz_mod(c1, c1, mask);
        dj_rand_unit(pp, E->z[2]);
        atom_first(pp, a[2], Ci, one, c1, E->z[2], 1, 1, t);
        dj_rand_unit(pp, k1);
        dj_rand_unit(pp, k2);
        dj_pow_Ns_sk(pp, a[0], k1, h);
        dj_pow_Ns_sk(pp, a[1], k2, 1);
        entry_challenge(pp, E, j, a, e);
        mpz_sub(E->c0, e, c1);
        mpz_mod(E->c0, E->c0, mask);
        respond(pp, E->z[0], k1, E->c0, rho, t);
        respond(pp, E->z[1], k2, E->c0, tau, t);
    } else {
        rand_chal(E->c0);
        dj_rand_unit(pp, E->z[0]);
        dj_rand_unit(pp, E->z[1]);
        atom_first(pp, a[0], ei, NULL, E->c0, E->z[0], h, 1, t);
        atom_first(pp, a[1], Ci, NULL, E->c0, E->z[1], 1, 1, t);
        mpz_set_ui(c1, 0);
        for (int i = 0; i < n; i++) {
            if (i == ia) continue;
            rand_chal(E->c[i]);
            mpz_add(c1, c1, E->c[i]);
            dj_rand_unit(pp, E->z[3 + i]);
            atom_first(pp, a[3 + i], ei, m[i], E->c[i], E->z[3 + i], h, 1, t);
        }
        dj_rand_unit(pp, k1);
        dj_rand_unit(pp, k2);
        dj_pow_Ns_sk(pp, a[2], k1, 1);
        dj_pow_Ns_sk(pp, a[3 + ia], k2, h);
        entry_challenge(pp, E, j, a, e);
        mpz_sub(e, e, E->c0);
        mpz_mod(e, e, mask);
        mpz_sub(E->c[ia], e, c1);
        mpz_mod(E->c[ia], E->c[ia], mask);
        respond(pp, E->z[2], k1, e, tau, t);
        respond(pp, E->z[3 + ia], k2, E->c[ia], rho, t);
    }

    if (tau_out) mpz_set(tau_out, tau);
    mpz_clears(rho, tau, k1, k2, e, c1, one, bm, ei, Ci, t, mask, NULL);
    free_vec(a, DJ_NZ(n));
    free_vec(m, n);
}

int dj_entry_verify(const DJParams *pp, const DJEntry *E, uint64_t j)
{
    const int h = pp->h, n = pp->n;
    int ok = 0;
    mpz_t *a = new_vec(DJ_NZ(n));
    mpz_t m, e, c1, one, ei, Ci, t, mask;
    mpz_inits(m, e, c1, one, ei, Ci, t, mask, NULL);
    mpz_set_ui(one, 1);
    mpz_set_ui(mask, 0);
    mpz_setbit(mask, DJ_CHAL_BITS);

    if (!is_unit(pp, E->e, h, t) || !is_unit(pp, E->C, 1, t)) goto done;
    if (!in_range(E->c0, mask)) goto done;
    for (int i = 0; i < n; i++) if (!in_range(E->c[i], mask)) goto done;
    for (int k = 0; k < DJ_NZ(n); k++) if (!in_range(E->z[k], pp->N)) goto done;
    if (!mpz_invert(ei, E->e, pp->Np[h + 1]) || !mpz_invert(Ci, E->C, pp->Np[2])) goto done;

    mpz_set_ui(c1, 0);
    mpz_set_ui(m, 1);
    for (int i = 0; i < n; i++) {
        mpz_add(c1, c1, E->c[i]);
        if (!atom_first(pp, a[3 + i], ei, m, E->c[i], E->z[3 + i], h, 0, t)) goto done;
        mpz_mul(m, m, pp->r);
    }
    mpz_mod(c1, c1, mask);
    if (!atom_first(pp, a[0], ei, NULL, E->c0, E->z[0], h, 0, t)) goto done;
    if (!atom_first(pp, a[1], Ci, NULL, E->c0, E->z[1], 1, 0, t)) goto done;
    if (!atom_first(pp, a[2], Ci, one, c1, E->z[2], 1, 0, t)) goto done;
    entry_challenge(pp, E, j, a, e);
    mpz_add(c1, c1, E->c0);
    mpz_mod(c1, c1, mask);
    ok = mpz_cmp(c1, e) == 0;
done:
    mpz_clears(m, e, c1, one, ei, Ci, t, mask, NULL);
    free_vec(a, DJ_NZ(n));
    return ok;
}

void dj_entry_leaf(const DJParams *pp, const DJEntry *E, uint8_t leaf[32])
{
    EVP_MD_CTX *h = EVP_MD_CTX_new();
    unsigned int len = 0;
    uint8_t d[32];
    static const uint8_t tag_d = 0x02, tag_leaf = 0x00;
    if (!h || EVP_DigestInit_ex(h, EVP_sha256(), NULL) != 1) abort();
    EVP_DigestUpdate(h, &tag_d, 1);
    put_elem(h, pp, E->e, pp->h);
    EVP_DigestFinal_ex(h, d, &len);
    if (EVP_DigestInit_ex(h, EVP_sha256(), NULL) != 1) abort();
    EVP_DigestUpdate(h, &tag_leaf, 1);
    put_elem(h, pp, E->C, 1);
    EVP_DigestUpdate(h, d, sizeof d);
    EVP_DigestFinal_ex(h, leaf, &len);
    EVP_MD_CTX_free(h);
}

/* ---- the count proof ----------------------------------------------------- */

static void count_challenge(const DJParams *pp, const mpz_t u, uint64_t M,
                            unsigned long c, const mpz_t a, mpz_t e)
{
    EVP_MD_CTX *h = chal_init(pp, "ECT-DJ-board-count");
    put_u64(h, M);
    put_u64(h, (uint64_t)c);
    put_elem(h, pp, u, 1);
    put_elem(h, pp, a, 1);
    chal_final(h, e);
}

/* u (1+N)^{-c} is an N-th power mod N^2 with root tau. */
void dj_count_prove(const DJParams *pp, const mpz_t u, uint64_t M, unsigned long c,
                    const mpz_t tau, mpz_t ce, mpz_t cz)
{
    mpz_t k, a, t;
    mpz_inits(k, a, t, NULL);
    dj_rand_unit(pp, k);
    dj_pow_Ns_sk(pp, a, k, 1);
    count_challenge(pp, u, M, c, a, ce);
    respond(pp, cz, k, ce, tau, t);
    mpz_clears(k, a, t, NULL);
}

int dj_count_verify(const DJParams *pp, const mpz_t u, uint64_t M, unsigned long c,
                    const mpz_t ce, const mpz_t cz)
{
    int ok = 0;
    mpz_t ui, sh, a, e, t, mask;
    mpz_inits(ui, sh, a, e, t, mask, NULL);
    mpz_setbit(mask, DJ_CHAL_BITS);
    if (c > M || !in_range(ce, mask) || !in_range(cz, pp->N)) goto done;
    if (!mpz_invert(ui, u, pp->Np[2])) goto done;
    mpz_set_ui(sh, c);
    if (!atom_first(pp, a, ui, sh, ce, cz, 1, 0, t)) goto done;
    count_challenge(pp, u, M, c, a, e);
    ok = mpz_cmp(e, ce) == 0;
done:
    mpz_clears(ui, sh, a, e, t, mask, NULL);
    return ok;
}

/* ---- Board and VerifyBoard ----------------------------------------------- */

typedef struct {
    const DJParams *pp;
    DJBoard *B;
    const DJBoard *Bc;
    const int *slot, *b;
    mpz_t *tau;
    uint8_t (*leaf)[32];
    atomic_int bad;
} Job;

static void prove_range(void *ctx, int tid, size_t lo, size_t hi)
{
    (void)tid;
    Job *J = ctx;
    for (size_t j = lo; j < hi; j++) {
        dj_entry_prove(J->pp, &J->B->E[j], j, J->slot[j], J->b[j], J->tau[j]);
        dj_entry_leaf(J->pp, &J->B->E[j], J->leaf[j]);
    }
}

static void verify_range(void *ctx, int tid, size_t lo, size_t hi)
{
    (void)tid;
    Job *J = ctx;
    for (size_t j = lo; j < hi; j++) {
        if (atomic_load(&J->bad)) return;
        if (!dj_entry_verify(J->pp, &J->Bc->E[j], j)) atomic_store(&J->bad, 1);
        dj_entry_leaf(J->pp, &J->Bc->E[j], J->leaf[j]);
    }
}

static void product(const DJParams *pp, const DJBoard *B, mpz_t u)
{
    mpz_set_ui(u, 1);
    for (size_t j = 0; j < B->M; j++) {
        mpz_mul(u, u, B->E[j].C);
        mpz_mod(u, u, pp->Np[2]);
    }
}

void dj_board(const DJParams *pp, DJBoard *B, const int *slot, const int *b,
              int threads)
{
    Job J = { pp, B, B, slot, b, malloc(B->M * sizeof(mpz_t)),
              malloc(B->M * 32), 0 };
    if (!J.tau || !J.leaf) abort();
    for (size_t j = 0; j < B->M; j++) mpz_init(J.tau[j]);
    board_parallel(threads, B->M, prove_range, &J);

    merkle_root((const uint8_t (*)[32])J.leaf, B->M, B->root);
    uint8_t msg[BOARD_MSG_BYTES];
    board_message(msg, B->M, B->root);
    B->siglen = sig_sign(&pp->sig, msg, sizeof msg, B->sig, pp->sig.sig_bytes);

    mpz_t u, tau;
    mpz_inits(u, tau, NULL);
    mpz_set_ui(tau, 1);
    B->c = 0;
    for (size_t j = 0; j < B->M; j++) {
        mpz_mul(tau, tau, J.tau[j]);
        mpz_mod(tau, tau, pp->N);
        B->c += (unsigned long)(b[j] != 0);
        mpz_clear(J.tau[j]);
    }
    product(pp, B, u);
    dj_count_prove(pp, u, B->M, B->c, tau, B->ce, B->cz);
    mpz_clears(u, tau, NULL);
    free(J.tau); free(J.leaf);
}

int dj_verify_board(const DJParams *pp, const DJBoard *B, int threads)
{
    if (B->M == 0 || B->c > B->M) return 0;
    Job J = { pp, NULL, B, NULL, NULL, NULL, malloc(B->M * 32), 0 };
    if (!J.leaf) abort();
    board_parallel(threads, B->M, verify_range, &J);
    int ok = !atomic_load(&J.bad);
    if (ok) {
        uint8_t root[32], msg[BOARD_MSG_BYTES];
        merkle_root((const uint8_t (*)[32])J.leaf, B->M, root);
        board_message(msg, B->M, root);
        ok = memcmp(root, B->root, 32) == 0 &&
             sig_verify(&pp->sig, msg, sizeof msg, B->sig, B->siglen);
    }
    if (ok) {
        mpz_t u;
        mpz_init(u);
        product(pp, B, u);
        ok = dj_count_verify(pp, u, B->M, B->c, B->ce, B->cz);
        mpz_clear(u);
    }
    free(J.leaf);
    return ok;
}

size_t dj_entry_bytes(const DJParams *pp)
{
    /* e_j in Z_{N^{h+1}}, C_j in Z_{N^2}, n + 1 challenges of 16 bytes and
       n + 3 responses mod N. */
    return dj_elem_len(pp, pp->h) + dj_elem_len(pp, 1)
         + (size_t)(pp->n + 1) * (DJ_CHAL_BITS / 8)
         + (size_t)DJ_NZ(pp->n) * pp->elem_bytes;
}

size_t dj_count_bytes(const DJParams *pp)
{
    return 8 + DJ_CHAL_BITS / 8 + pp->elem_bytes;
}

size_t dj_root_bytes(const DJParams *pp)
{
    return 32 + pp->sig.sig_bytes;
}
