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
#include "work_cloud.h"

/*! FUNCTION:  	WORK_cloud_search()
 *  SYNOPSIS:  	Run "cloud search" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_cloud_search( WORKER*  worker )
{
   /* workflow */
   ARGS*             args           = worker->args;
   TASKS*            tasks          = worker->tasks;
   TIMES*            times          = worker->times;
   CLOCK*            clok           = worker->clok;
   RESULT*           result         = worker->result;
   /* search parameters */
   SEQUENCE*         q_seq          = worker->q_seq;
   HMM_PROFILE*      t_prof         = worker->t_prof;
   int               Q              = q_seq->N;
   int               T              = t_prof->N;
   CLOUD_PARAMS*     cloud_params   = &(worker->cloud_params);
   /* working data */
   MATRIX_3D*        st_cloud_MX    = worker->st_cloud_MX;
   MATRIX_3D*        st_MX          = worker->st_MX;
   MATRIX_3D_SPARSE* st_SMX_fwd     = worker->st_SMX_fwd;
   MATRIX_3D_SPARSE* st_SMX_bck     = worker->st_SMX_bck;

   MATRIX_3D*        st_MX3         = worker->st_MX3;
   MATRIX_3D*        st_MX3_fwd     = worker->st_MX3_fwd;
   MATRIX_3D*        st_MX3_bck     = worker->st_MX3_bck;

   MATRIX_2D*        sp_MX          = worker->sp_MX;
   MATRIX_2D*        sp_MX_fwd      = worker->sp_MX_fwd;
   MATRIX_2D*        sp_MX_bck      = worker->sp_MX_bck;
    
   ALIGNMENT*        tr             = worker->traceback;

   EDGEBOUNDS*       edg_fwd        = worker->edg_fwd;
   EDGEBOUNDS*       edg_bck        = worker->edg_bck;
   EDGEBOUNDS*       edg_diag       = worker->edg_diag;
   EDGEBOUNDS*       edg_row        = worker->edg_row;
   EDGEBOUND_ROWS*   edg_rows_tmp   = worker->edg_rows_tmp;
   VECTOR_INT**      lb_vec         = worker->lb_vec;
   VECTOR_INT**      rb_vec         = worker->rb_vec;

   float             max_fwdsc, max_bcksc;
   float             inner_fwdsc, inner_bcksc;
   float             outer_fwdsc, outer_bcksc;
   float             max_sc, inner_maxsc, compo_sc;
   

   /* if performing linear fb-pruner, run cloud search  */
   if ( tasks->lin_cloud_fwd || tasks->lin_cloud_bck ) 
   {
      /* cloud forward */
      printf_vall("# ==> cloud forward (linear)...\n");
      CLOCK_Start(clok);
      run_Cloud_Forward_Linear( 
         q_seq, t_prof, Q, T, st_MX3, sp_MX, tr, edg_rows_tmp, edg_fwd, cloud_params, &inner_fwdsc, &max_fwdsc );
      CLOCK_Stop(clok);
      times->lin_cloud_fwd = CLOCK_Duration(clok);
      worker->scores->lin_cloud_fwd = max_fwdsc;
      #if DEBUG
      {
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, "test_output/my.cloud_fwd.lin.mx");
         EDGEBOUNDS_Save( edg_fwd, "test_output/my.fwd.edg" );
      }
      #endif

      /* cloud backward */
      printf_vall("# ==> cloud backward (linear)...\n");
      CLOCK_Start(clok);
      run_Cloud_Backward_Linear( 
         q_seq, t_prof, Q, T, st_MX3, sp_MX, tr, edg_rows_tmp, edg_bck, cloud_params, &inner_bcksc, &max_bcksc );
      CLOCK_Stop(clok);
      times->lin_cloud_bck = CLOCK_Duration(clok);
      worker->scores->lin_cloud_bck = max_bcksc;
      #if DEBUG
      {
         DP_MATRIX_Save( Q, T, debugger->test_MX, sp_MX, "test_output/my.cloud_bck.lin.mx" );
         EDGEBOUNDS_Save( edg_bck, "test_output/my.bck.edg" );
      }
      #endif

      /* take maximum score for threshold test */
      max_sc      = MAX( max_fwdsc, max_bcksc );
      /* take composite score for threshold test: */
      /* max score inside viterbi alignment range plus each search's score outside viterbi alignment range */
      inner_maxsc = MAX( inner_fwdsc, inner_bcksc );
      outer_fwdsc = max_fwdsc - inner_fwdsc;
      outer_bcksc = max_bcksc - inner_bcksc;
      compo_sc    = inner_maxsc + outer_fwdsc + outer_bcksc;

      worker->scores->threshold_cloud_max    = max_sc;
      worker->scores->threshold_cloud_compo  = compo_sc;
   }

   /* if performing quadratic bounded forward or backward  */
   if ( tasks->quad_bound_fwd || tasks->quad_bound_bck ) 
   {
      /* cloud forward */
      CLOCK_Start(clok);
      run_Cloud_Forward_Quad( 
         q_seq, t_prof, Q, T, st_MX, sp_MX, tr, edg_rows_tmp, lb_vec, rb_vec, edg_fwd, cloud_params, &inner_fwdsc, &max_fwdsc );
      CLOCK_Stop(clok);
      times->quad_cloud_fwd = CLOCK_Duration(clok);

      /* cloud backward */
      CLOCK_Start(clok);
      run_Cloud_Backward_Quad( 
         q_seq, t_prof, Q, T, st_MX, sp_MX, tr, edg_rows_tmp, lb_vec, rb_vec, edg_bck, cloud_params, &inner_bcksc, &max_bcksc );
      CLOCK_Stop(clok);
      times->quad_cloud_bck = CLOCK_Duration(clok);
   }
}


/*! FUNCTION:  	WORK_cloud_merge()
 *  SYNOPSIS:  	Run "cloud merge" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_cloud_merge( WORKER*  worker )
{

}

/*! FUNCTION:  	WORK_bound_forward_backward()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_bound_forward_backward( WORKER*    worker )
{
   ARGS*             args           = worker->args;
   TASKS*            tasks          = worker->tasks;

   SEQUENCE*         q_seq          = worker->q_seq;
   HMM_PROFILE*      t_prof         = worker->t_prof;
   int               Q              = q_seq->N;
   int               T              = t_prof->N;

   TIMES*            times          = worker->times;
   CLOCK*            clok           = worker->clok;
   RESULT*           result         = worker->result;
   NAT_SCORES*       scores         = worker->scores;

   float             sc             = 0.0f;

   MATRIX_3D*        st_cloud_MX    = worker->st_cloud_MX; 

   MATRIX_3D*        st_MX          = worker->st_MX;
   MATRIX_3D_SPARSE* st_SMX_fwd     = worker->st_SMX_fwd;
   MATRIX_3D_SPARSE* st_SMX_bck     = worker->st_SMX_bck;

   MATRIX_3D*        st_MX3         = worker->st_MX3;
   MATRIX_3D*        st_MX3_fwd     = worker->st_MX3_fwd;
   MATRIX_3D*        st_MX3_bck     = worker->st_MX3_bck;

   MATRIX_2D*        sp_MX          = worker->sp_MX;
   MATRIX_2D*        sp_MX_fwd      = worker->sp_MX_fwd;
   MATRIX_2D*        sp_MX_bck      = worker->sp_MX_bck;

   EDGEBOUNDS*       edg_fwd        = worker->edg_fwd;
   EDGEBOUNDS*       edg_bck        = worker->edg_bck;
   EDGEBOUNDS*       edg_diag       = worker->edg_diag;
   EDGEBOUNDS*       edg_row        = worker->edg_row;

   int status;

   /* if performing linear fb-pruner, run cloud search  */
   if ( tasks->lin_cloud_fwd || tasks->lin_cloud_bck ) 
   {
      /* merge edgebounds */
      printf_vall("# ==> merge (linear)...\n");
      CLOCK_Start(clok);
      EDGEBOUNDS_Union( Q, T, edg_fwd, edg_bck, edg_diag );
      CLOCK_Stop(clok);
      times->lin_merge = CLOCK_Duration(clok);
      #if DEBUG
      {
         EDGEBOUNDS_Save( edg_diag, "test_output/my.cloud.diags.edg");
      }
      #endif


      /* reorient edgebounds */
      printf_vall("# ==> reorient (linear)...\n");
      CLOCK_Start(clok);
      int precount  = EDGEBOUNDS_Count( edg_row );
      EDGEBOUNDS_Reorient_to_Row( Q, T, edg_diag, edg_row );
      CLOCK_Stop(clok);
      times->lin_reorient = CLOCK_Duration(clok);

      /* compute the number of cells in matrix computed */
      result->cloud_cells = EDGEBOUNDS_Count( edg_row );
      result->total_cells  = (Q+1) * (T+1);

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
    
   /* linear bounded forward */
   if ( tasks->lin_bound_fwd ) 
   {
      printf_vall("# ==> bound forward (linear)...\n");
      CLOCK_Start(clok);
      run_Bound_Forward_Linear( q_seq, t_prof, Q, T, st_MX3_fwd, sp_MX_fwd, edg_row, &sc );
      CLOCK_Stop(clok);
      times->lin_bound_fwd    = CLOCK_Duration(clok);
      scores->lin_bound_fwd   = sc;

      #if DEBUG
      {
         printf("# printing linear bound forward...\n");
         printf("# lin bound forward: %f\n", scores->lin_bound_fwd );
         // DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, "test_output/my.bound_fwd.lin.mx");
         // DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, sp_MX_fwd, tr, "test_output/my.bound_fwd.lin.viz.mx");
      }
      #endif
   }
   /* linear bounded backward */
   if ( tasks->lin_bound_bck ) 
   {
      printf_vall("# ==> bound backward (linear)...\n");
      CLOCK_Start(clok);
      run_Bound_Backward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX_bck, edg_row, &sc );
      CLOCK_Stop(clok);
      times->lin_bound_bck    = CLOCK_Duration(clok);
      scores->lin_bound_bck   = sc;
      #if DEBUG
      {
         printf_vall("# printing linear bound backward...\n");
         printf("# lin bound backward: %f\n", scores->lin_bound_bck );
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_bck, "test_output/my.bound_bck.lin.mx");
         // DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, sp_MX, tr, "test_output/my.bound_bck.lin.viz.mx");
      }
      #endif
   }
   /* get threshold score */
   if ( tasks->lin_bound_fwd || tasks->lin_bound_bck ) 
   {
      float fwdsc = ( tasks->lin_bound_fwd == true ? scores->lin_bound_fwd : -INF );
      float bcksc = ( tasks->lin_bound_fwd == true ? scores->lin_bound_bck : -INF );
      scores->threshold_bound_max = MAX( fwdsc, bcksc );
   }

   /* if performing sparse fb-pruner, initialize sparse matrix */
   if (tasks->sparse_bound_fwd || tasks->sparse_bound_bck)
   {
      #if (DEBUG)
      {
         /* for testing: cloud fills entire dp matrix */
         // EDGEBOUNDS_Cover_Matrix(edg_row, Q, T);
      }
      #endif

      MATRIX_3D_SPARSE_Shape_Like_Edgebounds( st_SMX_fwd, edg_row );
      MATRIX_3D_SPARSE_Copy( st_SMX_bck, st_SMX_fwd );
   }

   /* sparse bounded forward */
   if ( tasks->sparse_bound_fwd ) 
   {
      printf_vall("# ==> bound forward (sparse)...\n");
      CLOCK_Start(clok);
      run_Bound_Forward_Sparse( 
         q_seq, t_prof, Q, T, st_SMX_fwd, sp_MX_fwd, edg_row, NULL, &sc );
      scores->sparse_bound_fwd = sc;
      CLOCK_Stop(clok);

      /* compute the number of cells in matrix computed */
      result->cloud_cells = EDGEBOUNDS_Count( edg_row );
      result->total_cells  = (Q+1) * (T+1);

      #if DEBUG
      {
         printf("# printing sparse bound forward...\n");
         printf("# sparse bound forward: %f\n", scores->sparse_bound_fwd );
         
         MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_fwd, debugger->test_MX );
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_fwd, "test_output/my.bound_fwd.sp.mx");
      }
      #endif
      times->sp_bound_fwd        = CLOCK_Duration(clok);
      scores->sparse_bound_fwd   = sc;
   }
   /* sparse bounded backward */
   if ( tasks->sparse_bound_bck ) 
   {
      printf_vall("# ==> bound backward (sparse)...\n");
      CLOCK_Start(clok);
      run_Bound_Backward_Sparse( 
         q_seq, t_prof, Q, T, st_SMX_bck, sp_MX_bck, edg_row, NULL, &sc );
      CLOCK_Stop(clok);
      times->sp_bound_bck        = CLOCK_Duration(clok);
      scores->sparse_bound_bck   = sc;
      #if DEBUG
      {
         printf_vall("# printing sparse bound backward...\n");
         printf("# sparse bound backward: %f\n", scores->sparse_bound_bck );
         MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_bck, debugger->test_MX );
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_bck, "test_output/my.bound_bck.sp.mx");
      }
      #endif
   }

   /* if performing quadratic bounded forward or backward  */
   if ( tasks->quad_bound_fwd || tasks->quad_bound_bck ) 
   {
      /* merge edgebounds */
      CLOCK_Start(clok);
      EDGEBOUNDS_Union( Q, T, edg_fwd, edg_bck, edg_diag );
      CLOCK_Stop(clok);
      times->quad_merge = CLOCK_Duration(clok);

      /* reorient edgebounds */
      CLOCK_Start(clok);
      EDGEBOUNDS_Reorient_to_Row( Q, T, edg_diag, edg_row );
      CLOCK_Stop(clok);
      times->quad_reorient = CLOCK_Duration(clok);
   }
   /* quadratic bounded forward */
   if ( tasks->quad_bound_fwd ) 
   {
      printf_vall("# ==> bound forward (quad)...\n");
      CLOCK_Start(clok);
      run_Bound_Forward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, edg_row, &sc );
      CLOCK_Stop(clok);
      times->quad_bound_fwd   = CLOCK_Duration(clok);
      scores->quad_bound_fwd  = sc;
   }
   /* quadratic bounded backward */
   if ( tasks->quad_bound_bck ) 
   {
      printf_vall("# ==> bound backward (quad)...\n");
      CLOCK_Start(clok);
      run_Bound_Backward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, edg_row, &sc );
      CLOCK_Stop(clok);
      times->quad_bound_bck   = CLOCK_Duration(clok);
      scores->quad_bound_bck  = sc;
   }
}

/*! FUNCTION:  	WORK_bound_forward_backward_linear()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Uses the linear space matrix implementation.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_bound_forward_backward_linear( WORKER* worker )
{

}

/*! FUNCTION:  	WORK_bound_forward_backward_sparse()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Uses the sparse matrix implementation.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_bound_forward_backward_sparse( WORKER* worker )
{

}