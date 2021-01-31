/*******************************************************************************
 *  FILE:      work.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *
 *  AUTHOR:    Dave Rich
 *  NOTES:
 *    - This is being phased out.       
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

/* easel library */
#include "easel.h"
#include "esl_gumbel.h"
#include "esl_exponential.h"

/* header */
#include "_pipelines.h"

/* generic workflow */
void WORK_main_workflow( WORKER*  worker )
{

}

/* debug workflow */
void WORK_debug_workflow( WORKER*  worker )
{

}

/* search workflow */
void WORK_search_all_v_all( WORKER* worker )
{

}

/* search via list */
void WORK_search_list( WORKER* worker )
{

}

/* use sparse matrix to cover alignment */
void 
WORK_capture_alignment( WORKER* worker )
{
   FILE*          fp       = NULL;
   SEQUENCE*      q_seq    = worker->q_seq;
   HMM_PROFILE*   t_prof   = worker->t_prof;
   NAT_SCORES*    scores   = worker->scores;
   ALIGNMENT*     tr       = worker->traceback;
   float          sc;

   int   Q  = q_seq->N;
   int   T  = t_prof->N;

   /* create sparse matrix to fill edgebounds determined by cloud search */
   MATRIX_3D_SPARSE_Shape_Like_Edgebounds( worker->st_SMX_fwd, worker->edg_row );
   MATRIX_3D_SPARSE_Copy( worker->st_SMX_bck, worker->st_SMX_fwd );

   #if DEBUG 
   {
      MATRIX_2D_Fill( debugger->cloud_MX, 0.0 );
      MATRIX_2D_Cloud_Fill( debugger->cloud_MX, worker->st_SMX_fwd->edg_outer, 1.0 );
      MATRIX_2D_Cloud_Fill( debugger->cloud_MX, worker->st_SMX_fwd->edg_inner, 0.1 );
      fp = fopen("test_output/sparse_mx.mx", "w");
      MATRIX_2D_Dump( debugger->cloud_MX, fp );
      fclose(fp);
   }
   #endif

   // /* run forward */
   // run_Bound_Forward_Sparse( 
   //    worker->q_seq, worker->t_prof, q_seq->N, t_prof->N, worker->st_SMX_fwd, worker->sp_MX_fwd, worker->edg_row, NULL, &sc);
   // scores->sparse_bound_fwd = sc;

   // /* run backward */
   // run_Bound_Backward_Sparse( 
   //    worker->q_seq, worker->t_prof, q_seq->N, t_prof->N, worker->st_SMX_bck, worker->sp_MX_bck, worker->edg_row, NULL, &sc );
   // scores->sparse_bound_bck = sc;

   // /* Recover alignment */
   // run_Posterior_Sparse( 
   //    q_seq, t_prof, q_seq->N, t_prof->N, worker->st_SMX_fwd, worker->sp_MX, worker->edg_row, worker->trace_post );
   // ALIGNMENT_Dump(worker->trace_post, stdout);

   // /* generate alignments */
   // ALIGNMENT_Build_MMSEQS_Style( 
   //    worker->traceback, worker->q_seq, worker->t_prof );
   // ALIGNMENT_Build_HMMER_Style( 
   //    worker->traceback, worker->q_seq, worker->t_prof );
}



