#include "ect.h"
#include "elgamal.h"
#include "utils.h"
#include "nizk.h"
#include "schnorr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

int ect_mult_max = MULT_MAX;

/*    ECTTag    */

void ect_tag_new(ECTTag *w)
{
    bn_null(w->c1);   bn_null(w->c2);
    bn_null(w->c1_0); bn_null(w->c2_0);
    bn_new(w->c1);   bn_new(w->c2);
    bn_new(w->c1_0); bn_new(w->c2_0);
    w->hops = 0;
    w->n_proofs = 0;
    w->has_sig = 0;
    for (int i = 0; i < MAX_PROOFS; i++) {
        nizk_proof_new(&w->proofs[i]);
        bn_null(w->steps[i].c1_in);  bn_new(w->steps[i].c1_in);
        bn_null(w->steps[i].c2_in);  bn_new(w->steps[i].c2_in);
        bn_null(w->steps[i].c1_out); bn_new(w->steps[i].c1_out);
        bn_null(w->steps[i].c2_out); bn_new(w->steps[i].c2_out);
    }
    schnorr_sig_new(&w->sig);
}

void ect_tag_free(ECTTag *w)
{
    bn_free(w->c1);   bn_free(w->c2);
    bn_free(w->c1_0); bn_free(w->c2_0);
    for (int i = 0; i < MAX_PROOFS; i++) {
        nizk_proof_free(&w->proofs[i]);
        bn_free(w->steps[i].c1_in);  bn_free(w->steps[i].c2_in);
        bn_free(w->steps[i].c1_out); bn_free(w->steps[i].c2_out);
    }
    schnorr_sig_free(&w->sig);
    (void)w;
}

void ect_tag_copy(ECTTag *dst, const ECTTag *src)
{
    bn_copy(dst->c1,   src->c1);
    bn_copy(dst->c2,   src->c2);
    bn_copy(dst->c1_0, src->c1_0);
    bn_copy(dst->c2_0, src->c2_0);
    dst->hops     = src->hops;
    dst->n_proofs = src->n_proofs;
    dst->has_sig  = src->has_sig;
    for (int i = 0; i < src->n_proofs; i++) {
        dst->proofs[i].m = src->proofs[i].m;
        for (int j = 0; j < src->proofs[i].m; j++) {
            bn_copy(dst->proofs[i].A1[j], src->proofs[i].A1[j]);
            bn_copy(dst->proofs[i].A2[j], src->proofs[i].A2[j]);
            bn_copy(dst->proofs[i].z[j],  src->proofs[i].z[j]);
            bn_copy(dst->proofs[i].e[j],  src->proofs[i].e[j]);
        }
        bn_copy(dst->steps[i].c1_in,  src->steps[i].c1_in);
        bn_copy(dst->steps[i].c2_in,  src->steps[i].c2_in);
        bn_copy(dst->steps[i].c1_out, src->steps[i].c1_out);
        bn_copy(dst->steps[i].c2_out, src->steps[i].c2_out);
    }
    bn_copy(dst->sig.R, src->sig.R);
    bn_copy(dst->sig.s, src->sig.s);
}

/*    ECTParams    */

void ect_params_new(ECTParams *pp)
{
    bn_null(pp->p);      bn_new(pp->p);
    bn_null(pp->g);      bn_new(pp->g);
    bn_null(pp->q);      bn_new(pp->q);
    bn_null(pp->qh);     bn_new(pp->qh);
    bn_null(pp->pk);     bn_new(pp->pk);
    bn_null(pp->sk);     bn_new(pp->sk);
    bn_null(pp->pk_sig); bn_new(pp->pk_sig);
    bn_null(pp->sk_sig); bn_new(pp->sk_sig);
    for (int i = 0; i < N; i++)
        for (int t = 0; t < H; t++) {
            bn_null(pp->trace_table[i][t]);
            bn_new(pp->trace_table[i][t]);
        }
    pp->n_tagged = 0;
    memset(pp->tagged_users, 0, sizeof(pp->tagged_users));
}

void ect_params_free(ECTParams *pp)
{
    bn_free(pp->p); bn_free(pp->g); bn_free(pp->q);
    bn_free(pp->qh); bn_free(pp->pk); bn_free(pp->sk);
    bn_free(pp->pk_sig); bn_free(pp->sk_sig);
    for (int i = 0; i < N; i++)
        for (int t = 0; t < H; t++)
            bn_free(pp->trace_table[i][t]);
    (void)pp;
}

void ect_clear_tagged(ECTParams *pp)
{
    pp->n_tagged = 0;
    memset(pp->tagged_users, 0, sizeof(pp->tagged_users));
}

void ect_setup(ECTParams *pp)
{
    bn_read_str(pp->p, P_STR, strlen(P_STR), 10);
    bn_read_str(pp->g, G_STR, strlen(G_STR), 10);
    bn_read_str(pp->q, Q_STR, strlen(Q_STR), 10);
    pp->h = H; pp->r = R_BASE; pp->n = N;
    compute_qh(pp->qh, pp->q, pp->h);

    /* The serialisation widths in params.h/block.h are compile-time constants
       derived from P_BITS and Q_BITS.  If new literals are dropped in without
       updating those, every fixed-width buffer silently truncates, so check
       once here rather than failing deep inside a Fiat-Shamir hash. */
    if (bn_bits(pp->p) != P_BITS || bn_bits(pp->q) != Q_BITS) {
        fprintf(stderr,
                "ect_setup: P_BITS/Q_BITS are %d/%d but P_STR/Q_STR decode to "
                "%d/%d; update params.h\n",
                P_BITS, Q_BITS, (int)bn_bits(pp->p), (int)bn_bits(pp->q));
        exit(1);
    }

    /* G_STR generates the order-q^GENERATED_H subgroup.  Project it into
       the order-q^H subgroup when benchmarking a smaller hop budget. */
    if (pp->h < GENERATED_H) {
        bn_t projection;
        bn_null(projection);
        bn_new(projection);
        compute_qh(projection, pp->q, GENERATED_H - pp->h);
        bn_mxp(pp->g, pp->g, projection, pp->p);
        bn_free(projection);
    }
    elgamal_keygen(pp->sk, pp->pk, pp->g, pp->qh, pp->p);
    schnorr_keygen(pp->sk_sig, pp->pk_sig, pp->g, pp->qh, pp->p);

    bn_t m;
    bn_null(m);
    RLC_TRY {
        bn_new(m);
        for (int i = 0; i < N; i++) {
            /* m = r^i = identifier for user i+1 */
            user_identifier(m, i + 1, pp->r);

            /* depth 0: g^(m_i) mod p */
            bn_mxp(pp->trace_table[i][0], pp->g, m, pp->p);

            /* depth t: raise previous entry to q-th power */
            for (int t = 1; t < H; t++)
                bn_mxp(pp->trace_table[i][t],
                       pp->trace_table[i][t-1],
                       pp->q, pp->p);
        }
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(m);
    }
}

/*    Tag    */

void ect_tag(ECTTag *w, ECTParams *pp, int id, int b)
{
    bn_t m, r_out;
    bn_null(m); bn_null(r_out);
    RLC_TRY {
        bn_new(m); bn_new(r_out);
        if (b == 1)
            user_identifier(m, id, pp->r);
        else
            bn_zero(m);
        elgamal_enc(w->c1, w->c2, r_out,
                    m, pp->pk, pp->g, pp->qh, pp->p);
        bn_copy(w->c1_0, w->c1);
        bn_copy(w->c2_0, w->c2);
        w->hops = 0;
        schnorr_sign(&w->sig, pp->sk_sig, pp->pk_sig,
                     w->c1_0, w->c2_0, pp->g, pp->qh, pp->p);
        w->has_sig = 1;

        /* Authority records which real users it tagged for Trace lookup */
        if (b == 1 && pp->n_tagged < N) {
            /* Only add if not already recorded */
            int already = 0;
            for (int k = 0; k < pp->n_tagged; k++)
                if (pp->tagged_users[k] == id) { already = 1; break; }
            if (!already)
                pp->tagged_users[pp->n_tagged++] = id;
        }
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(m); bn_free(r_out);
    }
}

/*    Degrade    */

void ect_degrade(ECTTag *w, const ECTParams *pp,
                 bn_t *C1s, bn_t *C2s, int m, int real_idx)
{
    bn_t r_out, c1_old, c2_old;
    bn_null(r_out); bn_null(c1_old); bn_null(c2_old);
    RLC_TRY {
        bn_new(r_out); bn_new(c1_old); bn_new(c2_old);

        /* Save input ciphertext before modification */
        bn_copy(c1_old, w->c1);
        bn_copy(c2_old, w->c2);

        /* Raise to q-th power */
        ciphertext_pow_q(w->c1, w->c2, pp->q, pp->p);

        /* Rerandomise: c' = (c1^q * g^r', c2^q * pk^r') */
        elgamal_rerand(w->c1, w->c2, r_out,
                       pp->pk, pp->g, pp->qh, pp->p);

        /* Generate CDS one-out-of-many NIZK proof over input set C */
        if (w->n_proofs < MAX_PROOFS) {
            int idx = w->n_proofs;

            nizk_prove(&w->proofs[idx],
                       pp->pk,
                       C1s, C2s, m, real_idx,
                       w->c1, w->c2,
                       r_out,
                       pp->g, pp->q, pp->qh, pp->p);
            w->proofs[idx].m = m;

            /*  c1_in/c2_in is the ciphertext before this hop 
                c1_out/c2_out is the ciphertext after this hop. */
            bn_copy(w->steps[idx].c1_in,  C1s[real_idx]);
            bn_copy(w->steps[idx].c2_in,  C2s[real_idx]);
            bn_copy(w->steps[idx].c1_out, w->c1);
            bn_copy(w->steps[idx].c2_out, w->c2);

            w->n_proofs++;
        }

        w->has_sig = 0;
        w->hops++;

    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(r_out); bn_free(c1_old); bn_free(c2_old);
    }
}

/*    Merge    */

void ect_merge(ECTTag *out, ECTTag *tags[], int k, const ECTParams *pp)
{
    /* Merge is the public ECT operation (Algorithm eg-merge): component-wise
       ciphertext multiplication followed by rerandomization with a fresh
       encryption of zero.  Authentication of fresh board entries and
       enforcement of transitions belong to the host ledger's validity
       predicate. */
    bn_t xi;
    bn_null(xi);
    bn_new(xi);
    bn_set_dig(out->c1, 1);
    bn_set_dig(out->c2, 1);
    out->hops = 0;
    out->n_proofs = 0;
    out->has_sig = 0;

    for (int j = 0; j < k; j++) {
        elgamal_mul(out->c1, out->c2,
                    tags[j]->c1, tags[j]->c2, pp->p);
    }
    elgamal_rerand(out->c1, out->c2, xi, pp->pk, pp->g, pp->qh, pp->p);
    bn_free(xi);
}

/*    ect_verify_proofs -> full NIZK chain replay    */

int ect_verify_proofs(const ECTTag *w, const ECTParams *pp)
{
    if (w->n_proofs == 0) return 1;

    volatile int all_ok = 1;
    bn_t C1, C2;
    bn_null(C1); bn_null(C2);

    RLC_TRY {
        bn_new(C1); bn_new(C2);

        for (int i = 0; i < w->n_proofs; i++) {
            /* Build single-element input set C = {steps[i].c1_in} */
            bn_copy(C1, w->steps[i].c1_in);
            bn_copy(C2, w->steps[i].c2_in);

            int ok = nizk_verify(&w->proofs[i],
                                 pp->pk,
                                 &C1, &C2, 1,
                                 w->steps[i].c1_out,
                                 w->steps[i].c2_out,
                                 pp->g, pp->q, pp->qh, pp->p);
            if (!ok) {
                printf("  [verify_proofs] hop %d: NIZK proof INVALID\n", i);
                all_ok = 0;
                break;
            }
        }
    } RLC_CATCH_ANY {
        all_ok = 0;
    } RLC_FINALLY {
        bn_free(C1); bn_free(C2);
    }
    return (int)all_ok;
}

/*    Trace    */

/*    Meet-in-the-middle digit recovery    */

/*
  At depth t the projection X_t = gv_rem^{q^{h-1-t}} equals gamma^{d_t}, where
  gamma = g^{q^{h-1}} generates the order-q subgroup and

      d_t = sum_i c_i r^{live[i]},    c_i in [0, MULT_MAX].

  The authority knows the live positions (it issued the tags) but not the
  coefficients: c_i counts how many merged inputs carried identifier live[i]
  into depth t, and is 1 in the common case.

  Splitting the live set into halves A and B and writing d_t = e_A + e_B lets
  us tabulate gamma^{e_A} over every coefficient vector on A, then scan B for a
  match against X_t * gamma^{-e_B}.  That costs (M+1)^|A| memory and
  h*(M+1)^|B| probes instead of the (M+1)^K of direct enumeration.

  The base elements gamma^{r^i} = trace_table[i][h-1] do not depend on t, so the
  table is built once and reused across all h depths.

  Recovery is unambiguous: coefficients are bounded by MULT_MAX <= R_BASE-1, so
  sum_i c_i r^i < r^N < q and the base-r encoding is injective mod q.
 */

typedef struct {
    uint64_t *fp;      /* fingerprint of gamma^{e_A}; 0 marks an empty slot */
    uint64_t *code;    /* mixed-radix code of the coefficient vector on A */
    size_t    slots;   /* power of two */
    size_t    mask;
} MitmTable;

/* bn_fingerprint -> 64-bit digest of a group element, for hash lookup only.
   Equal bn_t values are normalised identically, so equal values agree here;
   collisions are resolved by a full bn_cmp against the recomputed product. */
static uint64_t bn_fingerprint(bn_t x)
{
    uint64_t hv = 1469598103934665603ULL;
    int limbs = x->used < 4 ? x->used : 4;
    for (int i = 0; i < limbs; i++)
        hv = (hv ^ (uint64_t)x->dp[i]) * 1099511628211ULL;
    hv = (hv ^ (uint64_t)x->used) * 1099511628211ULL;
    hv ^= hv >> 29;
    return hv ? hv : 1;              /* 0 is reserved for the empty slot */
}

static void mitm_free(MitmTable *T)
{
    free(T->fp);
    free(T->code);
    T->fp = NULL; T->code = NULL; T->slots = 0; T->mask = 0;
}

static int mitm_alloc(MitmTable *T, uint64_t entries)
{
    size_t slots = 1;
    while (slots < (size_t)(entries * 2)) slots <<= 1;
    T->fp   = (uint64_t *)calloc(slots, sizeof(uint64_t));
    T->code = (uint64_t *)malloc(slots * sizeof(uint64_t));
    if (!T->fp || !T->code) { mitm_free(T); return 0; }
    T->slots = slots;
    T->mask  = slots - 1;
    return 1;
}

static void mitm_insert(MitmTable *T, uint64_t fp, uint64_t code)
{
    size_t i = (size_t)fp & T->mask;
    while (T->fp[i]) {
        if (T->fp[i] == fp && T->code[i] == code) return;
        i = (i + 1) & T->mask;
    }
    T->fp[i] = fp;
    T->code[i] = code;
}

/* mitm_split -> how many live positions go in the tabulated half.  Balanced at
   K/2, then shrunk while (M+1)^a would exceed the 2^MITM_CAP_LOG2 entry
   ceiling, which trades probe time for bounded memory rather than failing. */
static int mitm_split(int K, int M)
{
    int a = (K + 1) / 2;
    double lg = log2((double)(M + 1));
    while (a > 0 && (double)a * lg > (double)MITM_CAP_LOG2) a--;
    return a;
}

/* ipow_u64 -> (M+1)^len, saturating so the caller can reject a hopeless split
   instead of overflowing. */
static uint64_t ipow_u64(int radix, int len, int *overflow)
{
    uint64_t v = 1;
    *overflow = 0;
    for (int i = 0; i < len; i++) {
        if (v > UINT64_MAX / (uint64_t)radix) { *overflow = 1; return UINT64_MAX; }
        v *= (uint64_t)radix;
    }
    return v;
}

/* code_digits -> expand a mixed-radix code into its coefficient vector. */
static void code_digits(uint64_t code, int len, int radix, uint8_t *out)
{
    for (int i = 0; i < len; i++) {
        out[i] = (uint8_t)(code % (uint64_t)radix);
        code /= (uint64_t)radix;
    }
}

static int ect_trace_internal(const ECTParams *pp, const ECTTag *w,
                               int verbose)
{
    volatile int found = 0;
    int M = ect_mult_max;
    if (M < 1) M = 1;
    if (M > R_BASE - 1) M = R_BASE - 1;
    const int radix = M + 1;

    bn_t gv, one, gv_rem, v_int, exp, Xt, accA, accB, Y, prodA, tmp;
    bn_t * volatile pw   = NULL;   /* [i*radix + c] = gamma_i^{ c} */
    bn_t * volatile ipw  = NULL;   /* [i*radix + c] = gamma_i^{-c} */
    bn_t * volatile digits = NULL;
    uint32_t * volatile coeffs = NULL;
    int      * volatile live   = NULL;
    uint8_t  * volatile digA   = NULL;
    uint8_t  * volatile digB   = NULL;
    uint8_t  * volatile cvA    = NULL;
    uint8_t  * volatile cvB    = NULL;
    MitmTable T = {NULL, NULL, 0, 0};
    int K = 0, a = 0, b = 0;
    volatile int npw = 0;

    bn_null(gv); bn_null(one); bn_null(gv_rem); bn_null(v_int);
    bn_null(exp); bn_null(Xt); bn_null(accA); bn_null(accB);
    bn_null(Y); bn_null(prodA); bn_null(tmp);

    RLC_TRY {
        bn_new(gv); bn_new(one); bn_new(gv_rem); bn_new(v_int);
        bn_new(exp); bn_new(Xt); bn_new(accA); bn_new(accB);
        bn_new(Y); bn_new(prodA); bn_new(tmp);
        bn_set_dig(one, 1);

        /* Live positions: the identifiers the authority holds traceable.  Any
           tagged id outside the capacity is a caller error, so say so rather
           than silently dropping it from the search. */
        live = (int *)malloc((size_t)pp->n * sizeof(int));
        if (!live) RLC_THROW(ERR_NO_MEMORY);
        for (int j = 0; j < pp->n_tagged; j++) {
            int idx = pp->tagged_users[j] - 1;
            if (idx < 0 || idx >= pp->n) {
                if (verbose)
                    printf("Trace: user %d outside capacity N=%d, not traceable\n",
                           pp->tagged_users[j], pp->n);
                continue;
            }
            live[K++] = idx;
        }

        /* g^v = c2 * (c1^sk)^{-1} mod p */
        elgamal_dec_gv(gv, w->c1, w->c2, pp->sk, pp->p);

        if (bn_cmp(gv, one) == RLC_EQ || K == 0) {
            if (verbose) printf("Trace: tag expired\n");
        } else {
            a = mitm_split(K, M);
            b = K - a;

            int ovA = 0, ovB = 0;
            uint64_t sizeA = ipow_u64(radix, a, &ovA);
            uint64_t sizeB = ipow_u64(radix, b, &ovB);
            if (ovA || ovB) {
                if (verbose)
                    printf("Trace: K=%d with MULT_MAX=%d exceeds the search "
                           "budget\n", K, M);
                RLC_THROW(ERR_NO_VALID);
            }

            /* Per-position powers gamma_i^{+-c}, c in [0, M]. */
            npw = K * radix;
            pw  = (bn_t *)malloc((size_t)npw * sizeof(bn_t));
            ipw = (bn_t *)malloc((size_t)npw * sizeof(bn_t));
            if (!pw || !ipw) RLC_THROW(ERR_NO_MEMORY);
            for (int i = 0; i < npw; i++) {
                bn_null(pw[i]);  bn_new(pw[i]);
                bn_null(ipw[i]); bn_new(ipw[i]);
            }
            for (int i = 0; i < K; i++) {
                bn_set_dig(pw[i * radix], 1);
                bn_set_dig(ipw[i * radix], 1);
                if (M >= 1) {
                    bn_copy(pw[i * radix + 1],
                            pp->trace_table[live[i]][pp->h - 1]);
                    bn_mod_inv(ipw[i * radix + 1],
                               pp->trace_table[live[i]][pp->h - 1], pp->p);
                }
                for (int c = 2; c <= M; c++) {
                    bn_mul(pw[i * radix + c], pw[i * radix + c - 1],
                           pp->trace_table[live[i]][pp->h - 1]);
                    bn_mod(pw[i * radix + c], pw[i * radix + c], pp->p);
                    bn_mul(ipw[i * radix + c], ipw[i * radix + c - 1],
                           ipw[i * radix + 1]);
                    bn_mod(ipw[i * radix + c], ipw[i * radix + c], pp->p);
                }
            }

            digA = (uint8_t *)calloc((size_t)(a > 0 ? a : 1), 1);
            digB = (uint8_t *)calloc((size_t)(b > 0 ? b : 1), 1);
            cvA  = (uint8_t *)calloc((size_t)(a > 0 ? a : 1), 1);
            cvB  = (uint8_t *)calloc((size_t)(b > 0 ? b : 1), 1);
            if (!digA || !digB || !cvA || !cvB) RLC_THROW(ERR_NO_MEMORY);

            /* Build the table over A once; it is independent of the depth. */
            if (!mitm_alloc(&T, sizeA)) RLC_THROW(ERR_NO_MEMORY);
            bn_set_dig(accA, 1);
            memset((void *)digA, 0, (size_t)(a > 0 ? a : 1));
            for (uint64_t code = 0; code < sizeA; code++) {
                mitm_insert(&T, bn_fingerprint(accA), code);
                /* odometer step: +1 on a digit, or wrap M -> 0 and carry */
                for (int j = 0; j < a; j++) {
                    if (digA[j] < M) {
                        digA[j]++;
                        bn_mul(accA, accA, pw[j * radix + 1]);
                        bn_mod(accA, accA, pp->p);
                        break;
                    }
                    digA[j] = 0;
                    bn_mul(accA, accA, ipw[j * radix + M]);
                    bn_mod(accA, accA, pp->p);
                }
            }

            bn_copy(gv_rem, gv);
            bn_zero(v_int);

            for (int t = 0; t < pp->h; t++) {
                /* X_t = gv_rem^{q^{h-1-t}} = gamma^{d_t} */
                bn_set_dig(exp, 1);
                for (int k = 0; k < (pp->h - 1 - t); k++)
                    bn_mul(exp, exp, pp->q);
                bn_mxp(Xt, gv_rem, exp, pp->p);
                if (bn_cmp(Xt, one) == RLC_EQ) continue;

                /* Scan B, maintaining accB = gamma^{-e_B}. */
                int hit = 0;
                uint64_t hitA = 0, hitB = 0;
                bn_set_dig(accB, 1);
                memset((void *)digB, 0, (size_t)(b > 0 ? b : 1));

                for (uint64_t codeB = 0; codeB < sizeB && !hit; codeB++) {
                    bn_mul(Y, Xt, accB);
                    bn_mod(Y, Y, pp->p);

                    uint64_t fp = bn_fingerprint(Y);
                    size_t   i  = (size_t)fp & T.mask;
                    while (T.fp[i]) {
                        if (T.fp[i] == fp) {
                            /* Confirm the fingerprint hit exactly. */
                            code_digits(T.code[i], a, radix, (uint8_t *)cvA);
                            bn_set_dig(prodA, 1);
                            for (int j = 0; j < a; j++) {
                                if (!cvA[j]) continue;
                                bn_mul(prodA, prodA, pw[j * radix + cvA[j]]);
                                bn_mod(prodA, prodA, pp->p);
                            }
                            if (bn_cmp(prodA, Y) == RLC_EQ) {
                                hit = 1; hitA = T.code[i]; hitB = codeB;
                                break;
                            }
                        }
                        i = (i + 1) & T.mask;
                    }
                    if (hit) break;

                    for (int j = 0; j < b; j++) {
                        if (digB[j] < M) {
                            digB[j]++;
                            bn_mul(accB, accB, ipw[(a + j) * radix + 1]);
                            bn_mod(accB, accB, pp->p);
                            break;
                        }
                        digB[j] = 0;
                        bn_mul(accB, accB, pw[(a + j) * radix + M]);
                        bn_mod(accB, accB, pp->p);
                    }
                }

                if (!hit) {
                    /* No coefficient vector within [0, MULT_MAX] explains this
                       depth.  Leaving gv_rem unchanged would corrupt every
                       later projection, so stop and report what is certain. */
                    if (verbose)
                        printf("Trace: depth %d not decodable with "
                               "MULT_MAX=%d\n", t, M);
                    break;
                }

                code_digits(hitA, a, radix, (uint8_t *)cvA);
                code_digits(hitB, b, radix, (uint8_t *)cvB);

                /* v_int += q^t * d_t and remove the depth-t term from gv_rem */
                bn_t qt, ri, contrib, prod, inv;
                bn_null(qt); bn_null(ri); bn_null(contrib);
                bn_null(prod); bn_null(inv);
                bn_new(qt); bn_new(ri); bn_new(contrib);
                bn_new(prod); bn_new(inv);

                bn_set_dig(qt, 1);
                for (int k = 0; k < t; k++) bn_mul(qt, qt, pp->q);
                bn_set_dig(prod, 1);

                for (int i2 = 0; i2 < K; i2++) {
                    int c = (i2 < a) ? cvA[i2] : cvB[i2 - a];
                    if (!c) continue;
                    user_identifier(ri, live[i2] + 1, pp->r);
                    bn_mul(contrib, qt, ri);
                    bn_mul_dig(contrib, contrib, (dig_t)c);
                    bn_add(v_int, v_int, contrib);
                    for (int e = 0; e < c; e++) {
                        bn_mul(prod, prod, pp->trace_table[live[i2]][t]);
                        bn_mod(prod, prod, pp->p);
                    }
                }
                bn_mod_inv(inv, prod, pp->p);
                bn_mul(gv_rem, gv_rem, inv);
                bn_mod(gv_rem, gv_rem, pp->p);

                bn_free(qt); bn_free(ri); bn_free(contrib);
                bn_free(prod); bn_free(inv);
            }

            if (bn_is_zero(v_int)) {
                if (verbose) printf("Trace: tag expired\n");
            } else {
                digits = (bn_t *)malloc((size_t)pp->h * sizeof(bn_t));
                coeffs = (uint32_t *)malloc((size_t)pp->n * sizeof(uint32_t));
                if (!digits || !coeffs) RLC_THROW(ERR_NO_MEMORY);
                for (int t = 0; t < pp->h; t++) {
                    bn_null(digits[t]); bn_new(digits[t]);
                }
                base_q_decompose(digits, v_int, pp->q, pp->h);
                for (int t = 0; t < pp->h; t++) {
                    if (bn_is_zero(digits[t])) continue;
                    memset((void *)coeffs, 0,
                           (size_t)pp->n * sizeof(uint32_t));
                    base_r_decompose((uint32_t *)coeffs, digits[t],
                                     pp->r, pp->n);
                    for (int i = 0; i < pp->n; i++) {
                        if (coeffs[i] > 0) {
                            if (verbose)
                                printf("  depth %d -> user %d (x%u)\n",
                                       t, i + 1, coeffs[i]);
                            found++;
                        }
                    }
                }
                if (found == 0 && verbose)
                    printf("Trace: no users found\n");
            }
        }

    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        bn_free(gv); bn_free(one); bn_free(gv_rem); bn_free(v_int);
        bn_free(exp); bn_free(Xt); bn_free(accA); bn_free(accB);
        bn_free(Y); bn_free(prodA); bn_free(tmp);
        if (pw)  { for (int i = 0; i < npw; i++) bn_free(pw[i]);  free((void *)pw);  }
        if (ipw) { for (int i = 0; i < npw; i++) bn_free(ipw[i]); free((void *)ipw); }
        if (digits) {
            for (int t = 0; t < pp->h; t++) bn_free(digits[t]);
            free((void *)digits);
        }
        free((void *)coeffs);
        free((void *)live);
        free((void *)digA); free((void *)digB);
        free((void *)cvA);  free((void *)cvB);
        mitm_free(&T);
    }
    return (int)found;
}

int ect_trace(const ECTParams *pp, const ECTTag *w)
{
    return ect_trace_internal(pp, w, 1);
}

int ect_trace_silent(const ECTParams *pp, const ECTTag *w)
{
    return ect_trace_internal(pp, w, 0);
}

int ect_trace_ciphertext(const ECTParams *pp,
                          const bn_t c1, const bn_t c2,
                          int verbose)
{
    ECTTag tmp;
    ect_tag_new(&tmp);
    bn_copy(tmp.c1, c1);
    bn_copy(tmp.c2, c2);
    int r = ect_trace_internal(pp, &tmp, verbose);
    ect_tag_free(&tmp);
    return r;
}

/*    ect_verify_tag_sig    */

int ect_verify_tag_sig(const ECTTag *w, const ECTParams *pp)
{
    if (!w->has_sig) return 0;
    /* Verify against the ORIGINAL ciphertext (c1_0, c2_0)*/
    return schnorr_verify(&w->sig, pp->pk_sig,
                          w->c1_0, w->c2_0,
                          pp->g, pp->qh, pp->p);
}
