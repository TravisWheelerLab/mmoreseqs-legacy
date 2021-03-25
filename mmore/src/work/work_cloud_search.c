/*******************************************************************************
 *  FILE:      work_cloud_search.c
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
#include "work_cloud_search.h"

/*! FUNCTION:  	WORK_cloud_search()
 *  SYNOPSIS:  	Run "cloud search" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_cloud_search( WORKER*  worker )
{
   /* linear cloud search */
   WORK_cloud_search_linear( worker );
   /* quadratic cloud search */      
   WORK_cloud_search_quadratic( worker );
}

/*! FUNCTION:  	WORK_cloud_search_linear()
 *  SYNOPSIS:  	Run linear-space "cloud search" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_cloud_search_linear( WORKER*  worker )
{
   ARGS*             args           = worker->args;
   TASKS*            tasks          = worker->tasks;
   CLOCK*            timer          = worker->timer;
   /* input data */
   SEQUENCE*         q_seq          = worker->q_seq;
   int               Q              = q_seq->N;
   HMM_PROFILE*      t_prof         = worker->t_prof;
   int               T              = t_prof->N;
   ALIGNMENT*        tr             = worker->trace_vit;
   CLOUD_PARAMS*     cloud_params   = &(worker->cloud_params);
   /* working data */
   MATRIX_3D*        st_MX3         = worker->st_MX3;
   MATRIX_2D*        sp_MX          = worker->sp_MX;
   EDGEBOUND_ROWS*   edg_rows_tmp   = worker->edg_rows_tmp;
   EDGEBOUNDS*       edg_fwd        = worker->edg_fwd;
   EDGEBOUNDS*       edg_bck        = worker->edg_bck;
   /* output data */
   TIMES*            times          = worker->times;
   RESULT*           result         = worker->result;
   ALL_SCORES*       scores         = &result->scores;
   SCORES*           finalsc        = &result->final_scores;
   float             inner_fwdsc, inner_bcksc, outer_fwdsc, outer_bcksc;
   float             max_fwdsc, max_bcksc, inner_maxsc, outer_maxsc;
   float             max_sc, compo_sc;

   /* if running linear cloud search  */
   if ( tasks->lin_cloud_fwd || tasks->lin_cloud_bck ) 
   {
      /* cloud forward */
      printf_vall("# ==> cloud forward (linear)...\n");
      CLOCK_Start( worker->timer );
      run_Cloud_Forward_Linear( 
         q_seq, t_prof, Q, T, st_MX3, sp_MX, tr, edg_rows_tmp, edg_fwd, cloud_params, &inner_fwdsc, &max_fwdsc );
      CLOCK_Stop( worker->timer );
      times->lin_cloud_fwd = CLOCK_Duration( worker->timer );
      scores->lin_cloud_fwd = max_fwdsc;
      printf("CLD FWD SCORES:: inner_max = %f, outer_max = %f\n", inner_fwdsc, max_fwdsc );
      #if DEBUG
      {
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, DEBUG_FOLDER "/my.cloud_fwd.lin.mx");
         EDGEBOUNDS_Save( edg_fwd, DEBUG_FOLDER "/my.fwd.edg" );
      }
      #endif

      /* cloud backward */
      printf_vall("# ==> cloud backward (linear)...\n");
      CLOCK_Start( worker->timer );
      run_Cloud_Backward_Linear( 
         q_seq, t_prof, Q, T, st_MX3, sp_MX, tr, edg_rows_tmp, edg_bck, cloud_params, &inner_bcksc, &max_bcksc );
      CLOCK_Stop( worker->timer );
      times->lin_cloud_bck = CLOCK_Duration( worker->timer );
      scores->lin_cloud_bck = max_bcksc;
      printf("CLD BCK SCORES:: inner_max = %f, outer_max = %f\n", inner_bcksc, max_bcksc );
      #if DEBUG
      {
         DP_MATRIX_Save( Q, T, debugger->test_MX, sp_MX, DEBUG_FOLDER "/my.cloud_bck.lin.mx" );
         EDGEBOUNDS_Save( edg_bck, DEBUG_FOLDER "/my.bck.edg" );
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
      printf("CLOUD SCORE: %f %f\n", compo_sc, max_sc);
      /* save thresholds */
      scores->threshold_cloud_max      = max_sc;
      scores->threshold_cloud_compo    = compo_sc;
      finalsc->cloud_natsc             = compo_sc;
   }
}

/*! FUNCTION:  	WORK_cloud_search_quadratic()
 *  SYNOPSIS:  	Run quadratic-space "cloud search" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_cloud_search_quadratic( WORKER*  worker )
{
   ARGS*             args           = worker->args;
   TASKS*            tasks          = worker->tasks;
   CLOCK*            timer          = worker->timer;
   /* input data */
   SEQUENCE*         q_seq          = worker->q_seq;
   int               Q              = q_seq->N;
   HMM_PROFILE*      t_prof         = worker->t_prof;
   int               T              = t_prof->N;
   ALIGNMENT*        tr             = worker->trace_vit;
   CLOUD_PARAMS*     cloud_params   = &(worker->cloud_params);
   /* working data */
   MATRIX_3D*        st_MX          = worker->st_MX;
   MATRIX_2D*        sp_MX          = worker->sp_MX;
   EDGEBOUND_ROWS*   edg_rows_tmp   = worker->edg_rows_tmp;
   EDGEBOUNDS*       edg_fwd        = worker->edg_fwd;
   EDGEBOUNDS*       edg_bck        = worker->edg_bck;
   VECTOR_INT**      lb_vec         = worker->lb_vec;
   VECTOR_INT**      rb_vec         = worker->rb_vec;
   /* output data */
   TIMES*            times          = worker->times;
   RESULT*           result         = worker->result;
   ALL_SCORES*       scores         = &result->scores;
   SCORES*           finalsc        = &result->final_scores;
   float             inner_fwdsc, inner_bcksc, outer_fwdsc, outer_bcksc;
   float             max_fwdsc, max_bcksc, inner_maxsc, outer_maxsc;
   float             max_sc, compo_sc;

   /* if running quadratic cloud search  */
   if ( tasks->quad_bound_fwd || tasks->quad_bound_bck ) 
   {
      /* cloud forward */
      CLOCK_Start(timer);
      run_Cloud_Forward_Quad( 
         q_seq, t_prof, Q, T, st_MX, sp_MX, tr, edg_rows_tmp, lb_vec, rb_vec, edg_fwd, cloud_params, &inner_fwdsc, &max_fwdsc );
      CLOCK_Stop(timer);
      times->quad_cloud_fwd = CLOCK_Duration(timer);

      /* cloud backward */
      CLOCK_Start(timer);
      run_Cloud_Backward_Quad( 
         q_seq, t_prof, Q, T, st_MX, sp_MX, tr, edg_rows_tmp, lb_vec, rb_vec, edg_bck, cloud_params, &inner_bcksc, &max_bcksc );
      CLOCK_Stop(timer);
      times->quad_cloud_bck = CLOCK_Duration(timer);

      /* take maximum score for threshold test */
      max_sc      = MAX( max_fwdsc, max_bcksc );
      /* take composite score for threshold test: */
      /* max score inside viterbi alignment range plus each search's score outside viterbi alignment range */
      inner_maxsc = MAX( inner_fwdsc, inner_bcksc );
      outer_fwdsc = max_fwdsc - inner_fwdsc;
      outer_bcksc = max_bcksc - inner_bcksc;
      compo_sc    = inner_maxsc + outer_fwdsc + outer_bcksc;
      printf("CLOUD SCORE: %f %f\n", compo_sc, max_sc);
      /* save thresholds */
      scores->threshold_cloud_max      = max_sc;
      scores->threshold_cloud_compo    = compo_sc;
      finalsc->cloud_natsc             = compo_sc;
   }
}
