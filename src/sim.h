#pragma once
#include <stdint.h>
#include "ect.h"
#include "block.h"

/*  * sim.h -> Simulation of the ECT protocol. */

#define N_VALIDATORS   10
#define QUORUM         N_VALIDATORS

 
/* Validator        */
 

typedef struct {
    int   id;
} Validator;

 
/* User             */
 

typedef struct {
    int      id;          /* 1-indexed*/
    uint64_t balance;
    ECTTag   tag;         /* watermark currently attached to this user's coins */
    int      tag_active;  /* 1 if a live (non-expired, non-dummy) tag is present */
} User;

 
/* Simulation state */
 

typedef struct {
    ECTParams   *pp;               /* authority parameters (owns all keys) */
    Validator    validators[N_VALIDATORS];
    User        *users;            /* array of n_users User structs */
    int          n_users;
    Ledger      *ledger;           /* shared append-only ledger */
} Simulation;



Simulation *sim_new (ECTParams *pp, int n_users);
void        sim_free(Simulation *sim);

 
/* Validator operations */
 

/*
  validator_process_block -> single validator's accept/reject decision.
 
  Each validator performs three checks:
    1. Hash chain: prev_hash matches SHA-256 of the last block.
    2. Authority signature: schnorr_verify on (c1_0, c2_0).
    3. NIZK proof: nizk_verify on the degradation proof (if present).
 
  Returns 1 (accept) or 0 (reject).
 */
int validator_process_block(Validator *v, const Block *b,
                            const ECTTag *tag,
                            const ECTParams *pp,
                            const Ledger *ledger);

/*
    sim_broadcast -> send a block to all validators and count accepts.
    Returns the number of validators that accepted (0..N_VALIDATORS).
 */
int sim_broadcast(Simulation *sim, Block *b, const ECTTag *tag);

/*
    sim_commit -> broadcast a block and append it if quorum is reached.
    Returns 1 on success (committed), 0 if quorum not met.
 */
int sim_commit(Simulation *sim, Block *b, const ECTTag *tag);

 
/* Transactions     */
 

/*
    sim_issue_tx -> authority issues a coin to receiver.

    trace_id > 0 : tagged coin for user trace_id (b=1)
    trace_id = 0 : dummy watermark (b=0, untraced)

    Creates a fresh ECTTag, packs it into a block, commits. On success, attaches the tag to the receiver's User struct.
    Returns 1 on success, 0 on failure.
 */
int sim_issue_tx(Simulation *sim,
                 int sender_id, int receiver_id,
                 uint64_t amount, int trace_id);

/*
  sim_spend_tx -> user-to-user coin transfer.
 
  Degrades the sender's tag by one hop (generating a NIZK proof), packs the degraded tag into a block, commits.
  On success, transfers the tag to the receiver's User struct.
  Returns 1 on success, 0 on failure.
 */
int sim_spend_tx(Simulation *sim,
                 int sender_id, int receiver_id,
                 uint64_t amount);

 
/* Tracing  */
 

/*
  sim_trace_ledger -> authority scans active user tags and reports.
*/
void sim_trace_ledger(Simulation *sim);

 
/* Reporting        */
 

void sim_print_ledger  (const Simulation *sim);
void sim_print_balances(const Simulation *sim);
