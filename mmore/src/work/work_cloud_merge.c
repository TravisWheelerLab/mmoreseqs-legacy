/*******************************************************************************
 *  FILE:      work_cloud_merge.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for cloud search / pruned / adaptive-banding forward-backward.
 *
 *  AUTHOR:    Dave Rich
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
#include "work_cloud_merge.h"

/*! FUNCTION:  	WORK_cloud_merge()
 *  SYNOPSIS:  	Run "cloud merge" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 *                Caller must have run WORK_cloud_search().
 */
void 
WORK_cloud_merge_and_reorient( WORKER*  worker )
{
   ARGS*             args           = worker->args;
   TASKS*            tasks          = worker->tasks;
   CLOCK*            timer          = worker->timer;
   /* input data */
   SEQUENCE*         q_seq          = worker->q_seq;
   int               Q              = q_seq->N;
   HMM_PROFILE*      t_prof         = worker->t_prof;
   int               T              = t_prof->N;
   /* working data */
   EDGEBOUNDS*       edg_fwd        = worker->edg_fwd;
   EDGEBOUNDS*       edg_bck        = worker->edg_bck;
   EDGEBOUNDS*       edg_diag       = worker->edg_diag;
   EDGEBOUNDS*       edg_row        = worker->edg_row;
   EDGEBOUND_ROWS*   edg_builder    = worker->edg_rows_tmp;
   /* output data */
   TIMES*            times          = worker->times;
   RESULT*           result         = worker->result;
   ALL_SCORES*       scores         = &result->scores;
   SCORES*           final          = &result->final_scores;

   /* if performing linear fb-pruner, run cloud search  */
   if ( tasks->lin_cloud_fwd || tasks->lin_cloud_bck ) 
   {
      /* merge edgebounds */
      printf_vall("# ==> merge...\n");
      CLOCK_Start( worker->timer );
      EDGEBOUNDS_Union( Q, T, edg_fwd, edg_bck, edg_diag );
      CLOCK_Stop( worker->timer );
      times->lin_merge = CLOCK_Duration( worker->timer );
      #if DEBUG
      {
         EDGEBOUNDS_Save( edg_diag, "test_output/my.cloud.diags.edg");
      }
      #endif

      int pre_diag   = EDGEBOUNDS_Count( edg_diag );
      int pre_row    = EDGEBOUNDS_Count( edg_row );

      /* reorient edgebounds */
      printf_vall("# ==> reorient...\n");
      CLOCK_Start(timer);
      EDGEBOUNDS_ReorientToRow( Q, T, edg_diag, edg_builder, edg_row );
      CLOCK_Stop(timer);
      times->lin_reorient = CLOCK_Duration(timer);

      int post_diag   = EDGEBOUNDS_Count( edg_diag );
      int post_row    = EDGEBOUNDS_Count( edg_row );

      /* correctness check */
      printf("CELL COUNTS: %d, %d => %d, %d\n", pre_diag, pre_row, post_diag, post_row );

      /* compute the number of cells in matrix computed */
      result->cloud_cells  = EDGEBOUNDS_Count( edg_row );
      result->total_cells  = (Q+1) * (T+1);
      result->perc_cells   = (float) result->cloud_cells / (float) result->total_cells;

      #if DEBUG
      {
         // MATRIX_2D_Fill( debugger->cloud_MX, 0 );
         // MATRIX_2D_Cloud_Fill( debugger->cloud_MX, edg_fwd, 1 );
         // MATRIX_2D_Cloud_Fill( debugger->cloud_MX, edg_bck, 2 );
         // DP_MATRIX_VIZ_Trace( debugger->cloud_MX, tr );
         // DP_MATRIX_VIZ_Color_Dump( debugger->cloud_MX, stdout );
         // DP_MATRIX_VIZ_Dump( debugger->cloud_MX, stdout );
         EDGEBOUNDS_Save( edg_row, "test_output/my.cloud.rows.edg");
      }
      #endif
   }
}
