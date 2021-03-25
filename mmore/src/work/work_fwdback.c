/*******************************************************************************
 *  FILE:      work_fwdback.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for forward-backward.
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
#include "work_fwdback.h"

/*! FUNCTION:  	WORK_forward_backward()
 *  SYNOPSIS:  	Run "cloud search" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_forward_backward( WORKER*  worker )
{
   ARGS*          args     = worker->args;
   TASKS*         tasks    = worker->tasks;
   CLOCK*         timer    = worker->timer;
   
   /* input data */
   SEQUENCE*      q_seq    = worker->q_seq;
   HMM_PROFILE*   t_prof   = worker->t_prof;
   int            Q        = q_seq->N;
   int            T        = t_prof->N;
   /* working data */
   MATRIX_3D*     st_MX    = worker->st_MX;
   MATRIX_3D*     st_MX3   = worker->st_MX3;
   MATRIX_2D*     sp_MX    = worker->sp_MX;
   ALIGNMENT*     tr       = worker->trace_vit;
   /* output data */
   TIMES*         times    = worker->times;
   RESULT*        result   = worker->result;
   ALL_SCORES*    scores   = &result->scores;
   SCORES*        finalsc  = &result->final_scores;
   float          sc;

   /* forward */
   if ( tasks->lin_fwd ) {
      printf_vall("# ==> forward (lin)...\n");
      CLOCK_Start(timer);
      run_Forward_Linear( 
         q_seq, t_prof, Q, T, worker->st_MX3_fwd, worker->sp_MX_fwd, &sc );
      CLOCK_Stop(timer);
      times->lin_fwd    = CLOCK_Duration(timer);
      scores->lin_fwd   = sc;
      #if DEBUG 
      {
         printf("# printing forward linear...\n");
         printf("# lin forward score: %f\n", scores->lin_fwd);
         DP_MATRIX_Save(Q, T, debugger->test_MX, worker->sp_MX_fwd, DEBUG_FOLDER "/my.fwd.lin.mx");
         DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, worker->sp_MX_fwd, tr, DEBUG_FOLDER "/my.fwd.lin.viz.mx");
      }
      #endif
   } 
   if ( tasks->quad_fwd ) {
      printf_vall("# ==> forward (quad)...\n");
      CLOCK_Start(timer);
      run_Forward_Quad( 
         q_seq, t_prof, Q, T, worker->st_MX_fwd, worker->sp_MX_fwd, &sc );
      CLOCK_Stop(timer);
      times->quad_fwd   = CLOCK_Duration(timer);
      scores->quad_fwd  = sc;
      #if DEBUG 
      {
         printf("# printing forward quadratic...\n");
         printf("# quad forward score: %f\n", scores->quad_fwd);
         DP_MATRIX_Save(Q, T, debugger->test_MX, worker->sp_MX_fwd, DEBUG_FOLDER "/my.fwd.quad.mx");
      }
      #endif
   }

   /* backward */
   if ( tasks->lin_bck ) {
      printf_vall("# ==> backward (lin)...\n");
      CLOCK_Start(timer);
      run_Backward_Linear( 
         q_seq, t_prof, Q, T, worker->st_MX3_bck, worker->sp_MX_bck, &sc );
      CLOCK_Stop(timer);
      times->lin_bck    = CLOCK_Duration(timer);
      scores->lin_bck   = sc;
      #if DEBUG 
      {
         printf("# printing backward linear...\n");
         printf("# lin backward score: %f\n", scores->lin_bck);
         DP_MATRIX_Save(Q, T, debugger->test_MX, worker->sp_MX_bck, DEBUG_FOLDER "/my.bck.lin.mx");
         DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, worker->sp_MX_bck, tr, DEBUG_FOLDER "/my.bck.lin.viz.mx");
      }
      #endif
   }
   if ( tasks->quad_bck ) {
      printf_vall("# ==> backward (quad)...\n");
      CLOCK_Start(timer);
      run_Backward_Quad( 
         q_seq, t_prof, Q, T, worker->st_MX_bck, worker->sp_MX_bck, &sc );
      CLOCK_Stop(timer);
      times->quad_bck   = CLOCK_Duration(timer);
      scores->quad_bck  = sc;
      #if DEBUG 
      {
         printf("# printing backward quadratic...\n");
         printf("# quad backward score: %f\n", scores->quad_bck);
         DP_MATRIX_Save(Q, T, debugger->test_MX, worker->sp_MX_bck, DEBUG_FOLDER "/my.bck.quad.mx");
      }
      #endif
   }
}