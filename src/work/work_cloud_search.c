/*******************************************************************************
 *  - FILE:      work_cloud_search.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for cloud search / pruned / adaptive-banding forward-backward.
 *******************************************************************************/

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
#include "work_cloud_search.h"

/*! FUNCTION:  	WORK_cloud_search_linear()
 *  SYNOPSIS:  	Run linear-space "cloud search" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void WORK_cloud_search_linear(WORKER* worker) {
  TASKS* tasks = worker->tasks;
  /* input data */
  SEQUENCE* q_seq = worker->q_seq;
  int Q = q_seq->N;
  HMM_PROFILE* t_prof = worker->t_prof;
  int T = t_prof->N;
  ALIGNMENT* tr = worker->trace_vit;
  CLOUD_PARAMS* cloud_params = &(worker->cloud_params);
  /* working data */
  MATRIX_3D* st_MX3 = worker->st_MX3;
  MATRIX_2D* sp_MX = worker->sp_MX;
  EDGEBOUND_ROWS* edg_rows_tmp = worker->edg_rows_tmp;
  EDGEBOUNDS* edg_fwd = worker->edg_fwd;
  EDGEBOUNDS* edg_bck = worker->edg_bck;
  /* output data */
  TIMES* times = worker->times;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;

  float inner_fwdsc, inner_bcksc, outer_fwdsc, outer_bcksc;
  float max_fwdsc, max_bcksc, inner_maxsc;
  float max_sc, sum_sc, compo_sc2;

  /* if running linear cloud search  */
  if (tasks->lin_cloud_fwd || tasks->lin_cloud_bck) {
    /* cloud forward */
    CLOCK_Start(worker->timer);
    run_Cloud_Forward_Linear(
        q_seq, t_prof, Q, T, st_MX3, sp_MX, tr, edg_rows_tmp, edg_fwd, cloud_params, &inner_fwdsc, &max_fwdsc);
    CLOCK_Stop(worker->timer);
    times->lin_cloud_fwd = CLOCK_Duration(worker->timer);
    scores->lin_cloud_fwd = max_fwdsc;
#if DEBUG
    {
      printf("CLD FWD SCORES:: inner_max = %f, outer_max = %f\n", inner_fwdsc, max_fwdsc);
      DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, DEBUG_FOLDER "/my.cloud_fwd.lin.000.mx");
      EDGEBOUNDS_Save(edg_fwd, DEBUG_FOLDER "/my.fwd.000.edg");
    }
#endif

    /* cloud backward */
    CLOCK_Start(worker->timer);
    run_Cloud_Backward_Linear(
        q_seq, t_prof, Q, T, st_MX3, sp_MX, tr, edg_rows_tmp, edg_bck, cloud_params, &inner_bcksc, &max_bcksc);
    CLOCK_Stop(worker->timer);
    times->lin_cloud_bck = CLOCK_Duration(worker->timer);
    scores->lin_cloud_bck = max_bcksc;
#if DEBUG
    {
      printf("CLD BCK SCORES:: inner_max = %f, outer_max = %f\n", inner_bcksc, max_bcksc);
      DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, DEBUG_FOLDER "/my.cloud_bck.lin.000.mx");
      EDGEBOUNDS_Save(edg_bck, DEBUG_FOLDER "/my.bck.000.edg");
    }
#endif

    /* Maximum score for threshold test */
    max_sc = MATH_Max(max_fwdsc, max_bcksc);

    /* Sum of scores for threshold test */
    sum_sc = MATH_Logsum_explicit(max_fwdsc, max_bcksc);

    /* Composite score for threshold test: */
    /* max score inside viterbi alignment range plus each search's score outside viterbi alignment range */
    inner_maxsc = MATH_Max(inner_fwdsc, inner_bcksc);

    /* Composite via Multiplication - This is a decent and simple approximater, but not really justified */

    /* Composite via Sum - (if outer is negative, it would create an error, but that shouldn't happen) */
    outer_fwdsc = MATH_Logdiff_explicit(max_fwdsc, inner_fwdsc);
    outer_bcksc = MATH_Logdiff_explicit(max_bcksc, inner_bcksc);
    compo_sc2 = MATH_Logsum_explicit(inner_maxsc, MATH_Logsum_explicit(outer_fwdsc, outer_bcksc));

#if DEBUG
    {
      printf("CLOUD SCORES: (compo1) %.3f (compo2) %.3f (max): %.3f (sum): %.3f\n", compo_sc1, compo_sc2, max_sc, sum_sc);
    }
#endif

    /* save thresholds */
    scores->threshold_cloud_max = max_sc;
    scores->threshold_cloud_compo = compo_sc2;
    finalsc->cloud_natsc = sum_sc;
  }
}
