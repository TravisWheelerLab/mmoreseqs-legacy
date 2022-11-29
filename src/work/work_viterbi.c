/*******************************************************************************
 *  - FILE:      work_fwdback.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for forward-backward.
 *******************************************************************************/

/* imports */
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
#include "work_viterbi.h"

/*! FUNCTION:  	WORK_viterbi_mmore()
 *  SYNOPSIS:  	Run Viterbi algorithm.
 *                Depends on <task> settings in <worker>.
 */
void WORK_viterbi_mmore(WORKER* worker) {
  CLOCK* timer = worker->timer;
  /* input data */
  SEQUENCE* q_seq = worker->q_seq;
  int Q = q_seq->N;
  HMM_PROFILE* t_prof = worker->t_prof;
  int T = t_prof->N;
  /* working data */
  MATRIX_3D* st_MX3 = worker->st_MX3;
  MATRIX_2D* sp_MX = worker->sp_MX;
  /* output data */
  TIMES* times = worker->times;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  float sc;

  /* Viterbi */
  CLOCK_Start(timer);
  run_Viterbi_Linear(q_seq, t_prof, Q, T, st_MX3, sp_MX, &sc);
  CLOCK_Stop(timer);

  times->lin_vit = CLOCK_Duration(timer);
  scores->lin_vit = sc;
  finalsc->viterbi_mmore_natsc = sc;

#if DEBUG
  {
    // printf("# printing viterbi linear...\n");
    DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, DEBUG_FOLDER "/my.vit.lin.000.mx");
  }
#endif
}

/*! FUNCTION:  	WORK_viterbi_sparse()
 *  SYNOPSIS:  	Run Viterbi algorithm.
 *                Depends on <task> settings in <worker>.
 */
void WORK_viterbi_sparse(WORKER* worker) {
  CLOCK* timer = worker->timer;
  /* input data */
  SEQUENCE* q_seq = worker->q_seq;
  int Q = q_seq->N;
  HMM_PROFILE* t_prof = worker->t_prof;
  int T = t_prof->N;
  EDGEBOUNDS* edg = worker->edg_row;
  /* working data */
  MATRIX_3D_SPARSE* st_SMX_vit = worker->st_SMX_optacc;
  MATRIX_2D* sp_MX_vit = worker->sp_MX_optacc;
  /* output data */
  TIMES* times = worker->times;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;
  float sc;

  CLOCK_Start(timer);
  run_Bound_Viterbi_Sparse(
      q_seq, t_prof, Q, T, edg, NULL,
      st_SMX_vit, sp_MX_vit, &sc);
  CLOCK_Stop(timer);

  times->sp_vit = CLOCK_Duration(timer);
  scores->sparse_vit = sc;
  finalsc->viterbi_natsc = sc;

#if DEBUG
  {
    fp = fopen(DEBUG_FOLDER "/my.vit.sparse.000.mx", "w+");
    MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_vit, debugger->test_MX);
    DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_vit, fp);
    fclose(fp);
  }
#endif
}

/*! FUNCTION:  WORK_viterbi_traceback_sparse()
 *  SYNOPSIS:  Run Viterbi Traceback algorithm.
 *             Depends on <task> settings in <worker>.
 */
void WORK_viterbi_traceback_sparse(WORKER* worker) {
  CLOCK* timer = worker->timer;
  /* input data */
  SEQUENCE* q_seq = worker->q_seq;
  int Q = q_seq->N;
  HMM_PROFILE* t_prof = worker->t_prof;
  int T = t_prof->N;
  EDGEBOUNDS* edg = worker->edg_row;
  ALIGNMENT* aln = worker->trace_vit;
  /* working data */
  MATRIX_3D_SPARSE* st_SMX_vit = worker->st_SMX_optacc;
  MATRIX_2D* sp_MX_vit = worker->sp_MX_optacc;
  /* output data */
  TIMES* times = worker->times;

  /* Optimal Alignment Traceback */
  CLOCK_Start(timer);

  run_Viterbi_Traceback_Sparse(
      q_seq, t_prof, Q, T, edg, NULL,
      st_SMX_vit, sp_MX_vit, aln);
  ALIGNMENT_Build_MMSEQS_Style(aln, q_seq, t_prof);
  ALIGNMENT_Build_HMMER_Style(aln, q_seq, t_prof);

  CLOCK_Stop(timer);
  times->sp_vit_trace += CLOCK_Duration(timer);

#if DEBUG
  {
    fp = fopen(DEBUG_FOLDER "/my.vit.traceback.000.mx", "w+");
    ALIGNMENT_Dump(aln, fp);
    fclose(fp);
  }
#endif
}
