/*******************************************************************************
 *  - FILE:      work_fwdback.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Tests if scores pass threshold filters.
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

/* header */
#include "_work.h"
#include "work_threshold.h"

/*! FUNCTION:  	WORK_thresholds_pval_to_eval()
 *  SYNOPSIS:  	Converts threshold scores from P-values to E-values.
 */
void WORK_thresholds_pval_to_eval(WORKER* worker) {
  ARGS* args = worker->args;
  STATS* stats = worker->stats;
  int db_size = stats->n_query_db;

  /* if thresholds given as p-values (default) */
  if (args->use_pvals == true) {
    // fprintf(stdout, "THRESHOLDS -> PVALS: %.3e => %.3e => %.3e => %.3e\n",
    //         args->threshold_vit, args->threshold_cloud, args->threshold_boundfwd, args->threshold_fwd);

    args->threshold_vit = STATS_Pval_to_Eval(args->threshold_vit, db_size);
    args->threshold_cloud = STATS_Pval_to_Eval(args->threshold_cloud, db_size);
    args->threshold_boundfwd = STATS_Pval_to_Eval(args->threshold_boundfwd, db_size);
    args->threshold_fwd = STATS_Pval_to_Eval(args->threshold_fwd, db_size);
  }
  /* if threshold given as e-values */
  else {
    fprintf(stdout, "WARNING: THRESHOLDS ALREADY GIVEN AS E-VALUES.\n");
  }

  // fprintf(stdout, "THRESHOLDS -> EVALS: %.3e => %.3e => %.3e => %.3e\n",
  //         args->threshold_vit, args->threshold_cloud, args->threshold_boundfwd, args->threshold_fwd);
}

/*! FUNCTION:  	WORK_viterbi_natsc_to_eval()
 *  SYNOPSIS:  	Converts Viterbi natscore to e-value.
 */
void WORK_viterbi_mmore_natsc_to_eval(WORKER* worker) {
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  HMM_PROFILE* t_prof = worker->t_prof;
  SEQUENCE* q_seq = worker->q_seq;
  STATS* stats = worker->stats;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  int db_size = worker->stats->n_query_db;
  float natsc, bitsc, pval, eval;

  natsc = finalsc->viterbi_mmore_natsc;

  STATS_Viterbi_Nats_to_Eval(
      natsc, &bitsc, NULL, &pval, &eval,
      t_prof->viterbi_dist, db_size, 0.0f, 0.0f);

  // printf_vhi("VITERBI_MMORE SCORES: %.3f %.3f %.3e %.3e\n",
  //            natsc, bitsc, pval, eval);

  finalsc->viterbi_mmore_bitsc = bitsc;
  finalsc->viterbi_mmore_eval = eval;
}

/*! FUNCTION:  	WORK_viterbi_natsc_to_eval()
 *  SYNOPSIS:  	Converts Viterbi natscore to e-value.
 */
void WORK_viterbi_natsc_to_eval(WORKER* worker) {
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  HMM_PROFILE* t_prof = worker->t_prof;
  SEQUENCE* q_seq = worker->q_seq;
  STATS* stats = worker->stats;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  int db_size = worker->stats->n_query_db;
  float natsc, bitsc, pval, eval;

  natsc = finalsc->viterbi_natsc;

  STATS_Viterbi_Nats_to_Eval(
      finalsc->viterbi_natsc, &bitsc, NULL, &pval, &eval,
      t_prof->viterbi_dist, db_size, 0.0f, 0.0f);

  // printf_vhi("VITERBI SCORES: %.3f %.3f %.3e %.3e\n",
  //            natsc, bitsc, pval, eval);

  finalsc->viterbi_bitsc = bitsc;
  finalsc->viterbi_eval = eval;
}

/*! FUNCTION:  	WORK_viterbi_test_threshold()
 *  SYNOPSIS:  	Tests viterbi score against threshold.
 *  RETURN:       <true/false> if passed threshold.
 */
bool WORK_viterbi_test_threshold(WORKER* worker) {
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  STATS* stats = worker->stats;

  bool is_passed;
  float vit_threshold = args->threshold_vit;
  float vit_eval = result->score_viterbi;

  /* check if passes */
  is_passed = (vit_eval < vit_threshold);
  // printf_vhi("VITERBI THRESHOLD (%s):  score = %.4e vs threshold = %.4e => %s\n",
  //            (args->is_run_viterbi_filter ? "ON" : "OFF"), vit_eval, vit_threshold, (is_passed ? "PASS" : "FAIL"));

  /* if filter is off, immediately passes */
  if (args->is_run_viterbi_filter == false) {
    result->is_passed_viterbi = true;
  }
  /* otherwise passes if eval is better (smaller) */
  else {
    result->is_passed_viterbi = is_passed;
  }
  /* stats tracks passing scores even when filter is off */
  if (is_passed) {
    stats->n_passed_viterbi += 1;
  }

  // return result->is_passed_viterbi;

  /*! TODO: Temporary bypass of Viterbi filter */
  return true;
}

/*! FUNCTION:  	WORK_cloud_natsc_to_eval()
 *  SYNOPSIS:  	Converts Cloud natscore to e-value.
 */
void WORK_cloud_natsc_to_eval(WORKER* worker) {
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  HMM_PROFILE* t_prof = worker->t_prof;
  SEQUENCE* q_seq = worker->q_seq;
  STATS* stats = worker->stats;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  int db_size = worker->stats->n_query_db;
  float natsc, bitsc, pval, eval;

  natsc = finalsc->cloud_natsc;

  STATS_Fwdback_Nats_to_Eval(
      natsc, &bitsc, NULL, &pval, &eval,
      t_prof->forward_dist, db_size, 0.0f, 0.0f);

  // printf_vhi("CLOUD SCORES: %.3f %.3f %.3e %.3e\n",
  //            finalsc->cloud_natsc, bitsc, pval, eval);

  finalsc->cloud_bitsc = bitsc;
  finalsc->cloud_eval = eval;
}

/*! FUNCTION:  	WORK_cloud_test_threshold()
 *  SYNOPSIS:  	Tests cloud eval against threshold.
 *  RETURN:       <true/false> if passed threshold.
 */
bool WORK_cloud_test_threshold(WORKER* worker) {
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  STATS* stats = worker->stats;

  bool is_passed;
  float cloud_threshold = args->threshold_cloud;
  float cloud_eval = finalsc->cloud_eval;

  is_passed = (cloud_eval < cloud_threshold);
  // printf_vhi("CLOUD THRESHOLD (%s):  score = %.4e vs threshold = %.4e => %s\n",
  //            (args->is_run_cloud_filter ? "ON" : "OFF"), cloud_eval, cloud_threshold, (is_passed ? "PASS" : "FAIL"));

  /* if filter is off, immediately passes */
  if (args->is_run_cloud_filter == false) {
    result->is_passed_cloud = true;
  }
  /* otherwise passes if eval is better (smaller) */
  else {
    result->is_passed_cloud = is_passed;
  }
  /* stats tracks passing scores even when filter is off */
  if (is_passed) {
    stats->n_passed_cloud += 1;
  }

  return result->is_passed_cloud;
}

/*! FUNCTION:  	WORK_bound_fwdback_natsc_to_eval()
 *  SYNOPSIS:  	Converts Forward natscore to E-value.
 */
void WORK_bound_fwdback_natsc_to_eval(WORKER* worker) {
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  HMM_PROFILE* t_prof = worker->t_prof;
  SEQUENCE* q_seq = worker->q_seq;
  STATS* stats = worker->stats;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  int db_size = worker->stats->n_query_db;
  float natsc, bitsc, pval, eval;

  natsc = finalsc->fwdback_natsc;

  STATS_Fwdback_Nats_to_Eval(
      natsc, &bitsc, NULL, &pval, &eval,
      t_prof->forward_dist, db_size, 0.0f, 0.0f);

  // printf_vhi("FWDBACK SCORES: %.3f %.3f %.3e %.3e\n",
  //            finalsc->fwdback_natsc, bitsc, pval, eval);

  finalsc->fwdback_bitsc = bitsc;
  finalsc->fwdback_eval = eval;
}

/*! FUNCTION:  	WORK_bound_fwdback_test_threshold()
 *  SYNOPSIS:  	Tests forward eval against threshold.
 *  RETURN:       <true/false> if passed threshold.
 */
bool WORK_bound_fwdback_test_threshold(WORKER* worker) {
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  STATS* stats = worker->stats;

  bool is_passed;
  float fwdback_threshold = args->threshold_boundfwd;
  float fwdback_eval = finalsc->fwdback_eval;

  is_passed = (fwdback_eval < fwdback_threshold);
  // printf_vhi("FWDBACK THRESHOLD (%s):  score = %.4e vs threshold = %.4e => %s\n",
  //            (args->is_run_boundfwd_filter ? "ON" : "OFF"), fwdback_eval, fwdback_threshold, (is_passed ? "PASS" : "FAIL"));

  /* if filter is off, immediately passes */
  if (args->is_run_boundfwd_filter == false) {
    result->is_passed_fwdback = true;
  }
  /* otherwise passes if eval is better (smaller) */
  else {
    result->is_passed_fwdback = is_passed;
  }
  /* stats tracks passing scores even when filter is off */
  if (is_passed) {
    stats->n_passed_fwdback += 1;
  }

  return result->is_passed_fwdback;
}

/*! FUNCTION:  	WORK_posterior_test_threshold()
 *  SYNOPSIS:  	Tests report eval against threshold (report is the bias-corrected, domain-corrected score).
 *  RETURN:       <true/false> if passed threshold.
 */
bool WORK_report_test_threshold(WORKER* worker) {
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  STATS* stats = worker->stats;

  bool is_passed;
  float report_threshold = args->threshold_report_eval;
  float report_eval = finalsc->eval;

  is_passed = (report_eval < report_threshold);
  // printf_vhi("REPORT THRESHOLD (%s):  score = %.4e vs threshold = %.4e => %s\n",
  //            (args->is_run_report_filter ? "ON" : "OFF"), report_eval, report_threshold, (is_passed ? "PASS" : "FAIL"));

  /* if filter is off, immediately passes */
  if (args->is_run_report_filter == false) {
    result->is_passed_report = true;
  }
  /* otherwise passes if eval is better (smaller) */
  else {
    result->is_passed_report = is_passed;
  }
  /* stats tracks passing scores even when filter is off */
  if (is_passed) {
    stats->n_passed_report += 1;
  }

  return result->is_passed_report;
}
