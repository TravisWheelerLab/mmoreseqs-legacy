/*******************************************************************************
 *  - FILE:      work_scoring.h
 *  - DESC:    Workflow Subroutines for Final Scoring.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../parsers/_parsers.h"
#include "../algs_linear/_algs_linear.h"
#include "../algs_quad/_algs_quad.h"
#include "../algs_naive/_algs_naive.h"
#include "../algs_sparse/_algs_sparse.h"
#include "../reporting/_reporting.h"

/* easel library */
#include "easel.h"
#include "esl_gumbel.h"
#include "esl_exponential.h"

/* header */
#include "_work.h"
#include "work_posterior.h"

/*! FUNCTION:  WORK_construct_scores()
 *  SYNOPSIS:  Construct final scores and evalues from results.
 */
void WORK_construct_scores(WORKER* worker) {
  FILE* fp = NULL;
  TASKS* tasks = worker->tasks;
  ARGS* args = worker->args;
  CLOCK* timer = worker->timer;
  /* input data */
  HMM_PROFILE* t_prof = worker->t_prof;
  int db_size = worker->stats->n_query_db;
  /* output data */
  TIMES* times = worker->times;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  DOMAIN_DEF* dom_def = worker->dom_def;
  STATS* stats = worker->stats;
  /* working vars */
  float natsc;
  float null_omega;
  float null1_hmm_bias;
  float null2_seq_bias;
  float presc;
  float seqsc;
  float pval;
  float eval;

  /* old assignments */
  natsc = scores->sparse_bound_fwd;
  null_omega = scores->null_omega;
  null1_hmm_bias = scores->null1_hmm_bias;
  null2_seq_bias = scores->null2_seq_bias;
  /* TODO: new assignments */
  // natsc             = finalsc->bound_fwdback_natsc;
  // null_omega        = finalsc->null_omega_natsc;
  // null1_hmm_bias    = finalsc->null1_hmm_bias_natsc;
  // null2_seq_bias    = finalsc->null2_seq_bias_natsc;

  /* TODO: Fix issue with null1 bias computation */
  // null1_hmm_bias = 0.0f;

  /* add prior probability to bias correction */
  null2_seq_bias = MATH_Sum(0.0f, null2_seq_bias + null_omega);

  /* compute bitscore, pval, and eval from natscore */
  STATS_Fwdback_Nats_to_Eval(natsc, &presc, &seqsc, &pval, &eval,
                             t_prof->forward_dist, db_size, null1_hmm_bias, null2_seq_bias);

#if DEBUG
  {
    printf("presc: %f, null2: %f, seqsc: %f, eval: %e\n",
    presc, null2_seq_bias, seqsc, eval);
  }
#endif
  

  /* capture results */
  finalsc->null2_seq_bias_natsc = null2_seq_bias;
  finalsc->null2_seq_bias_bitsc = STATS_Nats_to_Bits(null2_seq_bias);
  finalsc->nat_sc = natsc;
  finalsc->pre_sc = presc;
  finalsc->seq_sc = seqsc;
  finalsc->pval = pval;
  finalsc->eval = eval;
}

/*! FUNCTION:  WORK_construct_scores_bydom()
 *  SYNOPSIS:  Construct final scores and evalues from results.
 */
void WORK_construct_scores_bydom(WORKER* worker) {
}
