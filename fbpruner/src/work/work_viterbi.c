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
#include "work_viterbi.h"

/*! FUNCTION:  	WORK_viterbi_and_traceback()
 *  SYNOPSIS:  	Run Viterbi algorithm.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_viterbi_and_traceback( WORKER*  worker )
{
   ARGS*          args     = worker->args;
   TASKS*         tasks    = worker->tasks;
   NAT_SCORES*    scores   = worker->scores;
   TIMES*         times    = worker->times;
   CLOCK*         clok     = worker->clok;

   SEQUENCE*      q_seq    = worker->q_seq;
   HMM_PROFILE*   t_prof   = worker->t_prof;

   int   Q  = q_seq->N;
   int   T  = t_prof->N;

   MATRIX_3D*     st_MX    = worker->st_MX;
   MATRIX_3D*     st_MX3   = worker->st_MX3;
   MATRIX_2D*     sp_MX    = worker->sp_MX;
   ALIGNMENT*     tr       = worker->traceback;
   float          sc       = -1;

   /* Viterbi */
   if ( tasks->lin_vit ) {
      printf_vall("# ==> viterbi (lin)...\n");
      CLOCK_Start(clok);
      run_Viterbi_Linear( q_seq, t_prof, Q, T, st_MX, sp_MX, &sc );
      CLOCK_Stop(clok);
      times->lin_vit = CLOCK_Duration(clok);
      scores->lin_vit = sc;
      #if DEBUG 
      {
         printf("# printing viterbi linear...\n");
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, "test_output/my.vit.lin.mx");
      }
      #endif
   }
   if ( tasks->quad_vit ) {
      printf_vall("# ==> viterbi (quad)...\n");
      CLOCK_Start(clok);
      run_Viterbi_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, &sc );
      CLOCK_Stop(clok);
      times->quad_vit = CLOCK_Duration(clok);
      scores->quad_vit = sc;
      #if DEBUG 
      {
         printf("# printing viterbi quad...\n");
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, "test_output/my.vit.quad.mx");
      }
      #endif
   }

   /* Traceback */
   if ( tasks->lin_trace ) {
      fprintf(stderr, "ERROR: Operation not supported.\n");
      exit(EXIT_FAILURE);
      // TODO: linear traceback goes here
      // printf_vall("# ==> traceback (lin)...\n");
      // CLOCK_Start(clok);
      // run_Traceback_Quad(q_seq, t_prof, Q, T, st_MX, sp_MX, tr);
      // CLOCK_Stop(clok);
      // times->lin_trace = CLOCK_Duration(clok);
   }
   if ( tasks->quad_trace ) {
      printf_vall("# ==> traceback (quad)...\n");
      CLOCK_Start(clok);
      run_Traceback_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, tr );
      CLOCK_Stop(clok);
      times->quad_trace = CLOCK_Duration(clok);
   }
}