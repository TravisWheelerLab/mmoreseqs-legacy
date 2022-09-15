/*******************************************************************************
 *  - FILE:      work_maintenance.h
 *  - DESC:    Pipelines Workflow Subroutines.
 *             Computes optimal accuracy and traceback alignment.
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
#include "work_optacc.h"

/*! FUNCTION:  	WORK_optimal_accuracy()
 *  SYNOPSIS:  	Run Optimal Accuracy algorithm of Posterior.
 *                Depends on <task> settings in <worker>.
 */
void WORK_optimal_accuracy(WORKER* worker) {
  FILE* fp = NULL;
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  CLOCK* timer = worker->timer;
  /* input data */
  SEQUENCE* q_seq = worker->q_seq;
  int Q = q_seq->N;
  HMM_PROFILE* t_prof = worker->t_prof;
  int T = t_prof->N;
  EDGEBOUNDS* edg = worker->edg_row;
  /* working data */
  MATRIX_3D_SPARSE* st_SMX = worker->st_SMX;
  MATRIX_3D_SPARSE* st_SMX_post = worker->st_SMX_post;
  MATRIX_3D_SPARSE* st_SMX_opt = worker->st_SMX_optacc;
  MATRIX_2D* sp_MX = worker->sp_MX;
  MATRIX_2D* sp_MX_post = worker->sp_MX_post;
  /** TODO: this should use the optimal accuracy special mx, but is throwing errors */
  MATRIX_2D* sp_MX_opt = worker->sp_MX_fwd;
  /* output data */
  TIMES* times = worker->times;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* final_scores = &result->final_scores;
  float opt_sc;

  CLOCK_Start(timer);

  /* we need posterior in log space for computing the optimal accuraccy */
  MATRIX_3D_SPARSE_Log(st_SMX_post);
  MATRIX_2D_Log(sp_MX_post);
  MATRIX_3D_SPARSE_Fill(st_SMX_opt, 0.0f);
  MATRIX_2D_Fill(sp_MX_opt, 0.0f);

#if DEBUG
{
  printf("debug_folder: " DEBUG_FOLDER);
  fp = fopen(DEBUG_FOLDER "/my.optacc_postout.sp.000.mx", "w+");
  MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_post, debugger->test_MX);
  DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_post, fp);
  fclose(fp);

  fp = fopen(DEBUG_FOLDER "/my.optacc_optout.sp.000.mx", "w+");
  MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_opt, debugger->test_MX);
  DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_opt, fp);
  fclose(fp);
}
#endif

  run_OptimalAccuracy_Sparse(
      q_seq, t_prof, Q, T, edg, NULL,
      st_SMX_post, sp_MX_post, st_SMX_opt, sp_MX_opt, &opt_sc);

  CLOCK_Stop(timer);
  times->sp_optacc = CLOCK_Duration(timer);

  fprintf(stdout, "# ==> Optimal Accuracy (full cloud): %11.4f\n", opt_sc);
#if DEBUG
  {
    fp = fopen(DEBUG_FOLDER "/my.optacc.sp.000.mx", "w+");
    MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_opt, debugger->test_MX);
    DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_opt, fp);
    fclose(fp);
  }
#endif
}

/*! FUNCTION:  	WORK_optacc_traceback()
 *  SYNOPSIS:  	Run Traceback algorithm of Optimal Accuracy / Posterior.
 *                Depends on <task> settings in <worker>.
 */
void WORK_optacc_traceback(WORKER* worker) {
  FILE* fp = NULL;
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  CLOCK* timer = worker->timer;
  /* input data */
  SEQUENCE* q_seq = worker->q_seq;
  int Q = q_seq->N;
  HMM_PROFILE* t_prof = worker->t_prof;
  int T = t_prof->N;
  EDGEBOUNDS* edg = worker->edg_row;
  ALIGNMENT* aln = worker->trace_post;
  /* working data */
  MATRIX_3D_SPARSE* st_SMX = worker->st_SMX;
  MATRIX_3D_SPARSE* st_SMX_post = worker->st_SMX_post;
  MATRIX_3D_SPARSE* st_SMX_opt = worker->st_SMX_optacc;
  MATRIX_2D* sp_MX = worker->sp_MX;
  MATRIX_2D* sp_MX_post = worker->sp_MX_post;
  /** TODO: this should use the optimal accuracy special mx, but is throwing errors */
  MATRIX_2D* sp_MX_opt = worker->sp_MX_fwd;
  /* output data */
  TIMES* times = worker->times;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* final_scores = &result->final_scores;
  float opt_sc;

  /* Optimal Alignment Traceback */
  CLOCK_Start(timer);

  run_OptimalAccuracy_Traceback_Sparse(
      q_seq, t_prof, Q, T, edg, NULL,
      st_SMX_post, sp_MX_post, st_SMX_opt, sp_MX_opt, aln);
  ALIGNMENT_Build_MMSEQS_Style(aln, q_seq, t_prof);
  ALIGNMENT_Build_HMMER_Style(aln, q_seq, t_prof);

  CLOCK_Stop(timer);
  times->sp_optacc += CLOCK_Duration(timer);

#if DEBUG
  {
    fp = fopen(DEBUG_FOLDER "/my.optacc_traceback.sp.000.mx", "w+");
    ALIGNMENT_Dump(aln, fp);
    fclose(fp);
  }
#endif
}
