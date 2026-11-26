#include "sim.h"
#include "nizk.h"
#include "schnorr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


Simulation *sim_new(ECTParams *pp, int n_users)
{
    Simulation *sim = (Simulation *)malloc(sizeof(Simulation));
    sim->pp      = pp;
    sim->n_users = n_users;
    sim->ledger  = ledger_new();

    for (int i = 0; i < N_VALIDATORS; i++)
        sim->validators[i].id = i + 1;

    /* All users start with balance 1000 and no active tag */
    sim->users = (User *)malloc((size_t)n_users * sizeof(User));
    for (int i = 0; i < n_users; i++) {
        sim->users[i].id         = i + 1;
        sim->users[i].balance    = 1000;
        sim->users[i].tag_active = 0;
        ect_tag_new(&sim->users[i].tag);
    }

    return sim;
}

void sim_free(Simulation *sim)
{
    for (int i = 0; i < sim->n_users; i++)
        ect_tag_free(&sim->users[i].tag);
    free(sim->users);
    ledger_free(sim->ledger);
    free(sim);
}

 
/* Validator */
 

int validator_process_block(Validator *v, const Block *b,
                             const ECTTag *tag,
                             const ECTParams *pp,
                             const Ledger *ledger)
{
    /* Check 1: hash chain integrity */
    if (ledger->height > 0) {
        uint8_t expected[HASH_BYTES];
        block_hash(expected, &ledger->blocks[ledger->height - 1]);
        if (memcmp(b->raw + BLK_OFF_PREV_HASH, expected, HASH_BYTES) != 0) {
            printf("  [V%2d] REJECT block #%u: bad prev_hash\n",
                   v->id, block_get_height(b));
            return 0;
        }
    }


    /* Check 2: authority signature on (c1_0, c2_0). */
    {
        Watermark wm2;
        block_get_watermark(&wm2, b);
        bn_t wc1_0, wc2_0;
        bn_null(wc1_0); bn_null(wc2_0);
        bn_new(wc1_0);  bn_new(wc2_0);
        wm_unpack(NULL, NULL, NULL, NULL, wc1_0, wc2_0, &wm2);

        if (!bn_is_zero(wc1_0) && tag != NULL && tag->has_sig) {
            if (!ect_verify_tag_sig(tag, pp)) {
                printf("  [V%2d] REJECT block #%u: invalid authority sig\n",
                       v->id, block_get_height(b));
                bn_free(wc1_0); bn_free(wc2_0);
                return 0;
            }
        }
        bn_free(wc1_0); bn_free(wc2_0);
    }

    /* Check 3: NIZK proof content verification. */
    {
        Watermark wm;
        block_get_watermark(&wm, b);
        uint8_t hops = wm.data[WM_OFF_HOPS];

        if (hops > 0 && wm.data[WM_OFF_HAS] == 0) {
            printf("  [V%2d] REJECT block #%u: missing NIZK proof at hop %u\n",
                   v->id, block_get_height(b), hops);
            return 0;
        }

        if (hops > 0 && wm.data[WM_OFF_HAS] == 1) {
            volatile int proof_ok = 0;
            bn_t c1_cur, c2_cur, c1_before, c2_before, C1, C2;
            /* Heap-allocate: NIZKProof is ~12KB with MAX_INPUTS=128 */
            NIZKProof *pi = (NIZKProof *)malloc(sizeof(NIZKProof));

            bn_null(c1_cur);    bn_null(c2_cur);
            bn_null(c1_before); bn_null(c2_before);
            bn_null(C1);        bn_null(C2);

            RLC_TRY {
                bn_new(c1_cur);    bn_new(c2_cur);
                bn_new(c1_before); bn_new(c2_before);
                bn_new(C1);        bn_new(C2);
                nizk_proof_new(pi);

                /* Unpack current ciphertext and proof from this block */
                uint8_t h;
                wm_unpack(c1_cur, c2_cur, pi, &h, NULL, NULL, &wm);

                /* Determine c_before for this hop */
                if (hops == 1) {
                    /* hop 1: c_before = (c1_0, c2_0) from this block */
                    wm_unpack(NULL, NULL, NULL, NULL,
                              c1_before, c2_before, &wm);
                } else {
                    /* hop t>1: c_before = c1/c2 from the previous block */
                    uint32_t cur_height = block_get_height(b);
                    if (cur_height > 0 && ledger->height >= cur_height) {
                        Watermark prev_wm;
                        block_get_watermark(&prev_wm,
                                            &ledger->blocks[cur_height - 1]);
                        wm_unpack(c1_before, c2_before, NULL, NULL,
                                  NULL, NULL, &prev_wm);
                    }
                }

                /* Build single-element input set C = {c_before} */
                bn_copy(C1, c1_before);
                bn_copy(C2, c2_before);

                proof_ok = nizk_verify(pi, pp->pk,
                                       &C1, &C2, 1,
                                       c1_cur, c2_cur,
                                       pp->g, pp->q, pp->qh, pp->p);

            } RLC_CATCH_ANY {
                proof_ok = 0;
            } RLC_FINALLY {
                bn_free(c1_cur);    bn_free(c2_cur);
                bn_free(c1_before); bn_free(c2_before);
                bn_free(C1);        bn_free(C2);
                nizk_proof_free(pi);
                free(pi);
            }

            if (!proof_ok) {
                printf("  [V%2d] REJECT block #%u: invalid NIZK proof\n",
                       v->id, block_get_height(b));
                return 0;
            }
        }
    }

    return 1;
}

 
/* Broadcast and commit */
 

int sim_broadcast(Simulation *sim, Block *b, const ECTTag *tag)
{
    int accepts = 0;
    for (int i = 0; i < N_VALIDATORS; i++)
        accepts += validator_process_block(&sim->validators[i], b,
                                           tag, sim->pp, sim->ledger);
    return accepts;
}

int sim_commit(Simulation *sim, Block *b, const ECTTag *tag)
{
    int accepts = sim_broadcast(sim, b, tag);
    if (accepts >= QUORUM)
        return ledger_append(sim->ledger, b);
    printf("  [commit] FAILED: only %d/%d accepts (need %d)\n",
           accepts, N_VALIDATORS, QUORUM);
    return 0;
}

 
/* Transactions     */
 

int sim_issue_tx(Simulation *sim,
                 int sender_id, int receiver_id,
                 uint64_t amount, int trace_id)
{
    if (sender_id < 1 || sender_id > sim->n_users) return 0;
    if (receiver_id < 1 || receiver_id > sim->n_users) return 0;

    User *sender   = &sim->users[sender_id   - 1];
    User *receiver = &sim->users[receiver_id - 1];

    if (sender->balance < amount) {
        printf("  [issue_tx] sender %d: insufficient balance\n", sender_id);
        return 0;
    }

    /* Authority creates a fresh tag (signed, with c1_0/c2_0 stored) */
    ECTTag tag;
    ect_tag_new(&tag);
    ect_tag(&tag, sim->pp,
            trace_id > 0 ? trace_id : receiver_id,
            trace_id > 0 ? 1 : 0);

    /* Pack watermark -> fresh tag has no proof yet. */
    Watermark wm;
    wm_pack(&wm, tag.c1, tag.c2, NULL, (uint8_t)tag.hops,
            tag.c1_0, tag.c2_0);

    Block b;
    uint8_t prev[HASH_BYTES] = {0};
    if (sim->ledger->height > 0)
        block_hash(prev, &sim->ledger->blocks[sim->ledger->height - 1]);
    block_set_header(&b, prev,
                     sim->ledger->height,
                     (uint32_t)sender_id,
                     (uint32_t)receiver_id,
                     amount);
    block_set_watermark(&b, &wm);

    /* Pass the tag to validators so they can check the authority sig */
    int committed = sim_commit(sim, &b, &tag);
    if (committed) {
        sender->balance   -= amount;
        receiver->balance += amount;
        ect_tag_copy(&receiver->tag, &tag);
        receiver->tag_active = 1;
    }

    ect_tag_free(&tag);
    return committed;
}

int sim_spend_tx(Simulation *sim,
                 int sender_id, int receiver_id,
                 uint64_t amount)
{
    if (sender_id < 1 || sender_id > sim->n_users) return 0;
    if (receiver_id < 1 || receiver_id > sim->n_users) return 0;

    User *sender   = &sim->users[sender_id   - 1];
    User *receiver = &sim->users[receiver_id - 1];

    if (sender->balance < amount) {
        printf("  [spend_tx] sender %d: insufficient balance\n", sender_id);
        return 0;
    }

    /* Degrade the tag one hop with C = {current ciphertext}, m=1, idx=0 */
    if (sender->tag_active) {
        bn_t C1, C2;
        bn_null(C1); bn_null(C2);
        RLC_TRY {
            bn_new(C1); bn_new(C2);
            bn_copy(C1, sender->tag.c1);
            bn_copy(C2, sender->tag.c2);
            ect_degrade(&sender->tag, sim->pp, &C1, &C2, 1, 0);
        } RLC_CATCH_ANY {
            RLC_THROW(ERR_CAUGHT);
        } RLC_FINALLY {
            bn_free(C1); bn_free(C2);
        }
    }

    /* Pack watermark with degraded ciphertext, NIZK proof, and original  ciphertext */
    Watermark wm;
    NIZKProof *proof_ptr = (sender->tag_active && sender->tag.n_proofs > 0)
                         ? &sender->tag.proofs[sender->tag.n_proofs - 1]
                         : NULL;
    wm_pack(&wm,
            sender->tag.c1, sender->tag.c2,
            proof_ptr,
            (uint8_t)sender->tag.hops,
            sender->tag.c1_0, sender->tag.c2_0);

    Block b;
    uint8_t prev[HASH_BYTES] = {0};
    if (sim->ledger->height > 0)
        block_hash(prev, &sim->ledger->blocks[sim->ledger->height - 1]);
    block_set_header(&b, prev,
                     sim->ledger->height,
                     (uint32_t)sender_id,
                     (uint32_t)receiver_id,
                     amount);
    block_set_watermark(&b, &wm);

    /* Pass the sender's tag so validators can check sig and proof */
    int committed = sim_commit(sim, &b, sender->tag_active ? &sender->tag : NULL);
    if (committed) {
        sender->balance   -= amount;
        receiver->balance += amount;
        if (sender->tag_active) {
            ect_tag_copy(&receiver->tag, &sender->tag);
            receiver->tag_active = 1;
            sender->tag_active   = 0;
        }
    }

    return committed;
}

 
/* Tracing */
 

// void sim_trace_ledger(Simulation *sim)
// {
//     /* Authority scans active user tags and calls ect_trace on each. */
//     printf("  Scanning %u committed blocks...\n", sim->ledger->height);
//     printf("  Active tagged users:\n");

//     int found_any = 0;
//     for (int i = 0; i < sim->n_users; i++) {
//         User *u = &sim->users[i];
//         if (!u->tag_active) continue;
//         int n = ect_trace(sim->pp, &u->tag);
//         if (n > 0) found_any = 1;
//     }
//     if (!found_any)
//         printf("  No traced users found.\n");
// }
void sim_trace_ledger(Simulation *sim)
{
    /* Authority scans the committed ledger and traces each coin at its  most recent hop by decrypting the watermark ciphertext of each 
    terminal block using sk_enc. */
    printf("  Scanning %u committed blocks...\n", sim->ledger->height);
    printf("  Active tagged users:\n");
 
    int found_any = 0;
 
    for (uint32_t b = 0; b < sim->ledger->height; b++) {
        Watermark wm;
        block_get_watermark(&wm, &sim->ledger->blocks[b]);
        uint8_t hops_b = wm.data[WM_OFF_HOPS];
 
        /* Skip if this block is superseded by a later spend of the same coin */
        int is_terminal = 1;
        uint32_t recv_b = block_get_receiver(&sim->ledger->blocks[b]);
        for (uint32_t bp = b + 1; bp < sim->ledger->height; bp++) {
            uint32_t send_bp = block_get_sender(&sim->ledger->blocks[bp]);
            Watermark wm2;
            block_get_watermark(&wm2, &sim->ledger->blocks[bp]);
            uint8_t hops_bp = wm2.data[WM_OFF_HOPS];
            if (send_bp == recv_b && hops_bp == hops_b + 1) {
                is_terminal = 0;
                break;
            }
        }
        if (!is_terminal) continue;
 
        /* Decrypt and trace this terminal block */
        bn_t c1, c2;
        bn_null(c1); bn_null(c2);
        bn_new(c1);  bn_new(c2);
        wm_unpack(c1, c2, NULL, NULL, NULL, NULL, &wm);
        int n = ect_trace_ciphertext(sim->pp, c1, c2, 1);
        if (n > 0) found_any = 1;
        bn_free(c1); bn_free(c2);
    }
 
    if (!found_any)
        printf("  No traced users found.\n");
}
 
/* Prints */
 

void sim_print_ledger(const Simulation *sim)
{
    ledger_print(sim->ledger);
}

void sim_print_balances(const Simulation *sim)
{
    printf("  User balances:\n");
    for (int i = 0; i < sim->n_users; i++)
        printf("    user %2d: %4llu%s\n",
               sim->users[i].id,
               (unsigned long long)sim->users[i].balance,
               sim->users[i].tag_active ? "  [tagged]" : "");
}
