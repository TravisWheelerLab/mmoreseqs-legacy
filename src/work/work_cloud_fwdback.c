/*******************************************************************************
 *  FILE:      work_cloud.c
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
#include "work_cloud_fwdback.h"

/*! FUNCTION:  	WORK_bound_forward_backward()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 *                Caller must have run WORK_cloud_merge_and_reorient().
 */
void 
WORK_bound_fwdback( WORKER*    worker )
{
   /* linear bound forward/backward */
   WORK_bound_fwdback_linear( worker );
   /* quadratic bound forward/backward */
   WORK_bound_fwdback_quadratic( worker );
   /* sparse bound forward/backward */
   WORK_bound_fwdback_sparse( worker );
}

/*! FUNCTION:  	WORK_bound_forward_backward_linear()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Uses the linear space matrix implementation.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_bound_fwdback_linear( WORKER* worker )
{
   ARGS*             args           = worker->args;
   TASKS*            tasks          = worker->tasks;
   CLOCK*            timer          = worker->timer;
   /* input data */
   SEQUENCE*         q_seq          = worker->q_seq;
   int               Q              = q_seq->N;
   HMM_PROFILE*      t_prof         = worker->t_prof;
   int               T              = t_prof->N;
   EDGEBOUNDS*       edg_row        = worker->edg_row;
   /* working data */
   MATRIX_3D*        st_MX3         = worker->st_MX3;
   MATRIX_3D*        st_MX3_fwd     = worker->st_MX3_fwd;
   MATRIX_3D*        st_MX3_bck     = worker->st_MX3_bck;
   MATRIX_2D*        sp_MX          = worker->sp_MX;
   MATRIX_2D*        sp_MX_fwd      = worker->sp_MX_fwd;
   MATRIX_2D*        sp_MX_bck      = worker->sp_MX_bck;
   /* output data */
   TIMES*            times          = worker->times;
   RESULT*           result         = worker->result;
   ALL_SCORES*       scores         = &result->scores;
   SCORES*           finalsc        = &result->final_scores; 
   float             sc;

   /* linear bounded forward */
   if ( tasks->lin_bound_fwd ) 
   {
      printf_vall("# ==> bound forward (linear)...\n");
      CLOCK_Start( worker->timer );
      run_Bound_Forward_Linear( q_seq, t_prof, Q, T, st_MX3_fwd, sp_MX_fwd, edg_row, &sc );
      CLOCK_Stop( worker->timer );
      times->lin_bound_fwd    = CLOCK_Duration( worker->timer );
      scores->lin_bound_fwd   = sc;

      #if DEBUG
      {
         printf("# printing linear bound forward...\n");
         printf("# lin bound forward: %f\n", scores->lin_bound_fwd );
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, DEBUG_FOLDER "/my.bound_fwd.lin.mx");
         // DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, sp_MX_fwd, tr, DEBUG_FOLDER "/my.bound_fwd.lin.viz.mx");
      }
      #endif
   }
   /* linear bounded backward */
   if ( tasks->lin_bound_bck ) 
   {
      printf_vall("# ==> bound backward (linear)...\n");
      CLOCK_Start( worker->timer );
      run_Bound_Backward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX_bck, edg_row, &sc );
      CLOCK_Stop( worker->timer );
      times->lin_bound_bck    = CLOCK_Duration( worker->timer );
      scores->lin_bound_bck   = sc;

      #if DEBUG
      {
         printf_vall("# printing linear bound backward...\n");
         printf("# lin bound backward: %f\n", scores->lin_bound_bck );
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_bck, DEBUG_FOLDER "/my.bound_bck.lin.mx");
         // DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, sp_MX, tr, DEBUG_FOLDER "/my.bound_bck.lin.viz.mx");
      }
      #endif
   }
   /* get threshold score */
   if ( tasks->lin_bound_fwd || tasks->lin_bound_bck ) 
   {
      float fwdsc = ( tasks->lin_bound_fwd == true ? scores->lin_bound_fwd : -INF );
      float bcksc = ( tasks->lin_bound_fwd == true ? scores->lin_bound_bck : -INF );
      float maxsc = MAX( fwdsc, bcksc );
      finalsc->fwdback_natsc        = fwdsc;
      finalsc->bound_fwdback_natsc  = fwdsc;
      scores->threshold_bound_max   = maxsc;
      result->score_fwdback         = maxsc;
   }
}

/*! FUNCTION:  	WORK_bound_forward_backward_linear()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Uses the quadratic space matrix implementation.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_bound_fwdback_quadratic( WORKER* worker )
{
   ARGS*             args           = worker->args;
   TASKS*            tasks          = worker->tasks;
   CLOCK*            timer           = worker->timer;
   /* imput data */
   SEQUENCE*         q_seq          = worker->q_seq;
   int               Q              = q_seq->N;
   HMM_PROFILE*      t_prof         = worker->t_prof;
   int               T              = t_prof->N;
   EDGEBOUNDS*       edg_row        = worker->edg_row;
   /* working data */
   MATRIX_3D*        st_MX          = worker->st_MX;
   MATRIX_3D*        st_MX_fwd      = worker->st_MX_fwd;
   MATRIX_3D*        st_MX_bck      = worker->st_MX_bck;
   MATRIX_2D*        sp_MX          = worker->sp_MX;
   MATRIX_2D*        sp_MX_fwd      = worker->sp_MX_fwd;
   MATRIX_2D*        sp_MX_bck      = worker->sp_MX_bck;
   /* output data */
   TIMES*            times          = worker->times;
   RESULT*           result         = worker->result;
   ALL_SCORES*       scores         = &result->scores;
   SCORES*           final_scores   = &result->final_scores; 
   float             sc;

   /* quadratic bound forward */
   if ( tasks->quad_bound_fwd ) 
   {
      printf_vall("# ==> bound forward (quad)...\n");
      CLOCK_Start(timer);
      run_Bound_Forward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, edg_row, &sc );
      CLOCK_Stop(timer);
      times->quad_bound_fwd   = CLOCK_Duration(timer);
      scores->quad_bound_fwd  = sc;
   }
   /* quadratic bounded backward */
   if ( tasks->quad_bound_bck ) 
   {
      printf_vall("# ==> bound backward (quad)...\n");
      CLOCK_Start(timer);
      run_Bound_Backward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, edg_row, &sc );
      CLOCK_Stop(timer);
      times->quad_bound_bck   = CLOCK_Duration(timer);
      scores->quad_bound_bck  = sc;
   }
}

/*! FUNCTION:  	WORK_bound_forward_backward_sparse()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Uses the sparse matrix implementation.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_bound_fwdback_sparse( WORKER* worker )
{
   ARGS*                args           = worker->args;
   TASKS*               tasks          = worker->tasks;
   CLOCK*               timer          = worker->timer;
   /* input data */
   SEQUENCE*            q_seq          = worker->q_seq;
   int                  Q              = q_seq->N;
   HMM_PROFILE*         t_prof         = worker->t_prof;
   int                  T              = t_prof->N;
   EDGEBOUNDS*          edg_row        = worker->edg_row;
   /* working data */
   MATRIX_3D_SPARSE*    st_SMX         = worker->st_SMX;
   MATRIX_3D_SPARSE*    st_SMX_fwd     = worker->st_SMX_fwd;
   MATRIX_3D_SPARSE*    st_SMX_bck     = worker->st_SMX_bck;
   MATRIX_2D*           sp_MX          = worker->sp_MX;
   MATRIX_2D*           sp_MX_fwd      = worker->sp_MX_fwd;
   MATRIX_2D*           sp_MX_bck      = worker->sp_MX_bck;
   /* output data */
   TIMES*               times          = worker->times;
   RESULT*              result         = worker->result;
   ALL_SCORES*          scores         = &result->scores;
   SCORES*              finalsc        = &result->final_scores; 
   float                sc;

   /* sparse bounded forward */
   if ( tasks->sparse_bound_fwd ) 
   {
      printf_vall("# ==> bound forward (sparse)...\n");
      CLOCK_Start( timer );
      run_Bound_Forward_Sparse( 
         q_seq, t_prof, Q, T, st_SMX_fwd, sp_MX_fwd, edg_row, NULL, &sc );
      scores->sparse_bound_fwd = sc;
      CLOCK_Stop( timer );
      times->sp_bound_fwd        = CLOCK_Duration( timer );
      scores->sparse_bound_fwd   = sc;
      #if DEBUG
      {
         printf("# printing sparse bound forward...\n");
         printf("# sparse bound forward: %f\n", scores->sparse_bound_fwd );
         
         MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_fwd, debugger->test_MX );
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_fwd, DEBUG_FOLDER "/my.bound_fwd.sp.mx");
      }
      #endif
   }
   /* sparse bounded backward */
   if ( tasks->sparse_bound_bck ) 
   {
      printf_vall("# ==> bound backward (sparse)...\n");
      CLOCK_Start(timer);
      run_Bound_Backward_Sparse( 
         q_seq, t_prof, Q, T, st_SMX_bck, sp_MX_bck, edg_row, NULL, &sc );
      CLOCK_Stop(timer);
      times->sp_bound_bck        = CLOCK_Duration( timer );
      scores->sparse_bound_bck   = sc;
      #if DEBUG
      {
         printf_vall("# printing sparse bound backward...\n");
         printf("# sparse bound backward: %f\n", scores->sparse_bound_bck );
         MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_bck, debugger->test_MX );
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_bck, DEBUG_FOLDER "/my.bound_bck.sp.mx");
      }
      #endif
   }

   finalsc->fwdback_natsc = MAX( scores->sparse_bound_fwd, scores->sparse_bound_bck );
}
