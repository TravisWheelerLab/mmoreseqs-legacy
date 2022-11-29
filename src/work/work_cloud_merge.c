/*******************************************************************************
 *  - FILE:      work_cloud_merge.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for cloud search / pruned / adaptive-banding forward-backward.
 *******************************************************************************/

/* imports */
#include <stdio.h>


/* local imports */
#include "../objects/structs.h"
#include "../objects/_objects.h"

/* header */
#include "work_cloud_merge.h"

/*! FUNCTION:  	WORK_cloud_merge()
 *  SYNOPSIS:  	Run "cloud merge" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 *                Caller must have run WORK_cloud_search().
 */
void WORK_cloud_merge_and_reorient(WORKER* worker) {
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  CLOCK* timer = worker->timer;
  /* input data */
  SEQUENCE* q_seq = worker->q_seq;
  int Q = q_seq->N;
  HMM_PROFILE* t_prof = worker->t_prof;
  int T = t_prof->N;
  /* working data */
  EDGEBOUNDS* edg_fwd = worker->edg_fwd;
  EDGEBOUNDS* edg_bck = worker->edg_bck;
  EDGEBOUNDS* edg_diag = worker->edg_diag;
  EDGEBOUNDS* edg_row = worker->edg_row;
  EDGEBOUND_ROWS* edg_builder = worker->edg_rows_tmp;
  /* output data */
  TIMES* times = worker->times;
  RESULT* result = worker->result;

  /* if performing linear fb-pruner, run cloud search  */
  if (tasks->lin_cloud_fwd || tasks->lin_cloud_bck) {
    /* merge edgebounds */
    printf_vall("# ==> merge...\n");
    CLOCK_Start(worker->timer);
    EDGEBOUNDS_Union(Q, T, edg_fwd, edg_bck, edg_diag);
    CLOCK_Stop(worker->timer);
    times->lin_merge = CLOCK_Duration(worker->timer);
#if DEBUG
    {
      EDGEBOUNDS_Save(edg_diag, DEBUG_FOLDER "/my.cloud.diags.000.edg");
      int pre_diag = EDGEBOUNDS_Count(edg_diag);
      int pre_row = EDGEBOUNDS_Count(edg_row);
      int pre_diag_brt = EDGEBOUNDS_BruteCount(edg_diag);
      int pre_row_brt = EDGEBOUNDS_BruteCount(edg_row);
      /* correctness check */
      printf("PRE-REORIENT CELL_COUNTS: %d, %d => %d, %d\n", pre_diag, pre_row, pre_diag_brt, pre_row_brt);
    }
#endif

    /* reorient edgebounds */
    printf_vall("# ==> reorient...\n");
    CLOCK_Start(timer);
    EDGEBOUNDS_ReorientToRow(Q, T, edg_diag, edg_builder, edg_row);
    EDGEBOUNDS_Index(edg_row);
    CLOCK_Stop(timer);
    times->lin_reorient = CLOCK_Duration(timer);

    /* compute the number of cells in matrix computed */
    result->cloud_cells = EDGEBOUNDS_Count(edg_row);
    result->total_cells = (Q + 1) * (T + 1);
    result->perc_cells = (float)result->cloud_cells / (float)result->total_cells;

#if DEBUG
    {
      int post_diag = EDGEBOUNDS_Count(edg_diag);
      int post_row = EDGEBOUNDS_Count(edg_row);
      int post_diag_brt = EDGEBOUNDS_BruteCount(edg_diag);
      int post_row_brt = EDGEBOUNDS_BruteCount(edg_row);
      /* correctness check */
      printf("POST-REORIENT CELL_COUNTS: %d, %d => %d, %d\n", post_diag, post_row, post_diag_brt, post_row_brt);
      EDGEBOUNDS_Save(edg_row, DEBUG_FOLDER "/my.cloud.rows.000.edg");
    }
#endif
  }
}
