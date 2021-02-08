/*******************************************************************************
 *  FILE:      work_sparse_mx.h
 *  PURPOSE:   Workflow Subroutines for Posterior Algorithms.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
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
#include "work_sparse_mx.h"

/* build sparse matrix */
/*! FUNCTION:  	
 *  SYNOPSIS:  	
 */
void 
WORK_build_sparse_matrix(  WORKER* worker )
{
   FILE*          fp             = NULL;
   TASKS*         tasks          = worker->tasks;
   ARGS*          args           = worker->args;
   CLOCK*         timer          = worker->timer;
   /* input data */
   SEQUENCE*      q_seq          = worker->q_seq;
   int            Q              = q_seq->N;
   HMM_PROFILE*   t_prof         = worker->t_prof;
   int            T              = t_prof->N;
   EDGEBOUNDS*    edg            = worker->edg_row;
   /* output data */
   TIMES*         times          = worker->times;
   RESULT*        result         = worker->result;
   ALL_SCORES*    scores         = &result->scores;
   SCORES*        final_scores   = &result->final_scores;

   /* if running full fwdbackward or pruned (for comparison testing) */
   if ( worker->args->is_run_full == true ) {
      /* cloud fills entire dp matrix */
      EDGEBOUNDS_Cover_Matrix( edg, Q, T );
      result->cloud_cells  = EDGEBOUNDS_Count( edg );
   }

   CLOCK_Start( worker->timer );

   /* to build first matrix using <edg_row> as template */
   MATRIX_3D_SPARSE_Shape_Like_Edgebounds( worker->st_SMX_fwd, edg );
   /* TODO: this fill step may be optimized out later */
   // MATRIX_3D_SPARSE_Fill_Outer( worker->st_SMX_fwd, -INF );
   /* for other sparse matrices, simply copy first sparse matrix */
   MATRIX_3D_SPARSE_Copy( worker->st_SMX_bck, worker->st_SMX_fwd );
   // if ( worker->st_SMX_post != worker->st_SMX_bck ) {
   //    MATRIX_3D_SPARSE_Copy( worker->st_SMX_post, worker->st_SMX_fwd );
   // }
   // if ( worker->st_SMX_optacc != worker->st_SMX_fwd ) {
   //    MATRIX_3D_SPARSE_Copy( worker->st_SMX_post, worker->st_SMX_fwd );
   // }

   CLOCK_Stop( worker->timer );
   worker->times->sp_build_mx = CLOCK_Duration( worker->timer );
}

