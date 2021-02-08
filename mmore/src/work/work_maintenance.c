/*******************************************************************************
 *  FILE:      work_maintenance.h
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Initialization and Cleanup.
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
#include "work_maintenance.h"

/*! FUNCTION:  	WORK_init()
 *  SYNOPSIS:  	Initialize <worker> parts that are not handled by WORKER_Create().
 *                Allocate data structs according to settings in its <args> and <tasks>.
 */
void 
WORK_init( WORKER* worker )
{
   TASKS*   tasks = worker->tasks;
   ARGS*    args  = worker->args;

   /* initialize logrithmic sum table */
   logsum_Init();

   STR read_mode     = "r";
   STR write_mode    = "w";
   STR append_mode   = "a";
   /* input files */
   worker->q_seq_file         = FILER_Create( args->q_filepath, read_mode );
   worker->t_prof_file        = FILER_Create( args->t_filepath, read_mode );
   worker->q_index_file       = FILER_Create( args->q_indexpath, read_mode );
   worker->t_index_file       = FILER_Create( args->t_indexpath, read_mode );
   worker->mmseqs_file        = FILER_Create( args->mmseqs_res_filepath, read_mode );
   worker->hitlist_file       = FILER_Create( args->hitlist_filepath, read_mode );
   /* output files */
   worker->output_file        = FILER_Create( args->output_filepath, write_mode );
   if ( args->is_redirect_stdout == false ) {
      worker->output_file->fp = stdout;
   }
   worker->tblout_file        = FILER_Create( args->tblout_filepath, write_mode );
   worker->m8out_file         = FILER_Create( args->m8out_filepath, write_mode );
   worker->myout_file         = FILER_Create( args->myout_filepath, write_mode );
   worker->mydomout_file      = FILER_Create( args->mydomout_filepath, write_mode );
   worker->mytimeout_file     = FILER_Create( args->mytimeout_filepath, write_mode );
   worker->mythreshout_file   = FILER_Create( args->mythreshout_filepath, write_mode );

   /* target and profile structures */
   worker->q_seq           = SEQUENCE_Create();
   worker->t_seq           = SEQUENCE_Create();
   worker->t_prof          = HMM_PROFILE_Create();
   worker->hmm_bg          = HMM_BG_Create();
   /* target and profile indexes */
   worker->q_index         = F_INDEX_Create();
   worker->t_index         = F_INDEX_Create();
   /* results in from mmseqs and out for general searches */
   worker->mmseqs_data     = M8_RESULTS_Create();
   worker->results         = RESULTS_Create();
   worker->result          = ERROR_malloc( sizeof(RESULT) );
   /* data structs for viterbi alignment search */
   worker->trace_vit       = ALIGNMENT_Create();
   worker->trace_post      = ALIGNMENT_Create();
   /* data structs for cloud edgebounds */
   worker->edg_fwd         = EDGEBOUNDS_Create();
   worker->edg_bck         = EDGEBOUNDS_Create();
   worker->edg_diag        = EDGEBOUNDS_Create();
   worker->edg_row         = EDGEBOUNDS_Create();
   /* row-wise edgebounds */
   worker->edg_rows_tmp    = EDGEBOUND_ROWS_Create();
   for ( int i=0; i<3; i++ ) {
      worker->lb_vec[i]    = VECTOR_INT_Create();
      worker->rb_vec[i]    = VECTOR_INT_Create();
   }
   /* cloud search parameters */
   worker->cloud_params.alpha    = worker->args->alpha;
   worker->cloud_params.beta     = worker->args->beta;
   worker->cloud_params.gamma    = worker->args->gamma;
   /* create necessary dp matrices */
   /* quadratic */
   worker->st_MX_fwd       = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   worker->st_MX_bck       = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   // worker->st_MX_post      = MATRIX_3D_Create( NUM_NORMAL_STATES, 1, 1 );
   // worker->st_MX_optacc    = MATRIX_3D_Create( NUM_NORMAL_STATES, 1, 1 );
   worker->st_MX_post      = worker->st_MX_bck;
   worker->st_MX_optacc    = worker->st_MX_fwd;
   worker->st_MX           = worker->st_MX_bck;
   /* linear */
   worker->st_MX3_fwd      = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   worker->st_MX3_bck      = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   worker->st_MX3          = worker->st_MX_fwd;
   /* sparse */
   worker->st_SMX_fwd      = MATRIX_3D_SPARSE_Create();
   worker->st_SMX_bck      = MATRIX_3D_SPARSE_Create();
   /* posterior and optacc matrices can be recycled */
   // worker->st_SMX_post     = MATRIX_3D_SPARSE_Create();
   // worker->st_SMX_optacc   = MATRIX_3D_SPARSE_Create();
   worker->st_SMX_post     = worker->st_SMX_bck;
   worker->st_SMX_optacc   = worker->st_SMX_fwd;
   worker->st_SMX          = worker->st_SMX_fwd;
   /* special state */
   worker->sp_MX_fwd       = MATRIX_2D_Create( NUM_SPECIAL_STATES, 1 );
   worker->sp_MX_bck       = MATRIX_2D_Create( NUM_SPECIAL_STATES, 1 );
   // worker->sp_MX_post      = MATRIX_2D_Create( NUM_SPECIAL_STATES, 1 );
   // worker->sp_MX_optacc    = MATRIX_2D_Create( NUM_SPECIAL_STATES, 1 );
   worker->sp_MX_post      = worker->sp_MX_bck;
   worker->sp_MX_optacc    = worker->sp_MX_fwd;
   worker->sp_MX           = worker->sp_MX_fwd;
   /* domain definition */
   worker->dom_def         = DOMAIN_DEF_Create();

   /* init totals */
   WORK_times_init( worker, worker->times_totals );
}

/*! FUNCTION:  	WORK_reuse()
 *  SYNOPSIS:  	Resize and reallocate data structs in <worker> for problem size.
 */
void 
WORK_reuse( WORKER* worker )
{
   TASKS*   tasks    = worker->tasks; 
   int      T        = worker->t_prof->N;
   int      Q        = worker->q_seq->N;

   /* clear trace_vit and resize */
   ALIGNMENT_Reuse( worker->trace_vit, Q, T );
   ALIGNMENT_Reuse( worker->trace_post, Q, T );
   /* clear edgebounds and resize */
   EDGEBOUNDS_Reuse( worker->edg_fwd, Q, T );
   EDGEBOUNDS_Reuse( worker->edg_bck, Q, T );
   EDGEBOUNDS_Reuse( worker->edg_diag, Q, T );
   EDGEBOUNDS_Reuse( worker->edg_row, Q, T );
   /* clear row-wise edgebounds and resize */
   EDGEBOUND_ROWS_Reuse( worker->edg_rows_tmp, Q, T );
   /* domain definitions */
   DOMAIN_DEF_Reuse( worker->dom_def );

   /* matrix for quadratic algs */
   if ( tasks->quadratic ) {
      MATRIX_3D_Reuse_Clean( worker->st_MX_fwd,    NUM_NORMAL_STATES,  Q+1, T+1 );
      MATRIX_3D_Reuse_Clean( worker->st_MX_bck,    NUM_NORMAL_STATES,  Q+1, T+1 );
      // MATRIX_3D_Reuse_Clean( worker->st_MX_post,   NUM_NORMAL_STATES,  Q+1, T+1 );
      // MATRIX_3D_Reuse_Clean( worker->st_MX_optacc, NUM_NORMAL_STATES,  Q+1, T+1 );
   }
   /* matrix for linear algs */
   if ( tasks->linear ) {
      MATRIX_3D_Reuse_Clean( worker->st_MX3_fwd, NUM_NORMAL_STATES,  3, (Q+1)+(T+1) );
      MATRIX_3D_Reuse_Clean( worker->st_MX3_bck, NUM_NORMAL_STATES,  3, (Q+1)+(T+1) );
      worker->st_MX3 = worker->st_MX3_fwd;
   }
   /* matrices for special states */
   if ( tasks->quadratic || tasks->linear ) {
      MATRIX_2D_Reuse_Clean( worker->sp_MX_fwd,    NUM_SPECIAL_STATES, Q+1);
      MATRIX_2D_Reuse_Clean( worker->sp_MX_bck,    NUM_SPECIAL_STATES, Q+1);
      // MATRIX_2D_Reuse_Clean( worker->sp_MX_post,   NUM_SPECIAL_STATES, Q+1);
      // MATRIX_2D_Reuse_Clean( worker->sp_MX_optacc, NUM_SPECIAL_STATES, Q+1);
      worker->sp_MX = worker->sp_MX_fwd;
   }
   /* sparse matrices */
   MATRIX_3D_SPARSE_Reuse( worker->st_SMX_fwd );
   MATRIX_3D_SPARSE_Reuse( worker->st_SMX_bck );
   // MATRIX_3D_SPARSE_Reuse( worker->st_SMX_post )
   // MATRIX_3D_SPARSE_Reuse( worker->st_SMX_optacc );
   worker->st_SMX = worker->st_SMX_fwd;

   #if DEBUG 
   {
      DEBUGGER_Reuse( debugger, Q, T );
   }
   #endif

   /* TODO: remove this? */
   // MATRIX_3D_Clean( worker->st_MX );
   // MATRIX_3D_Clean( worker->st_MX3 );
}

/*! FUNCTION:  	WORK_cleanup()
 *  SYNOPSIS:  	Clean up and free data allocated by WORK_init().
 */
void 
WORK_cleanup( WORKER* worker )
{
   ARGS* args; 

   /* input files */
   worker->q_seq_file         = FILER_Destroy( worker->q_seq_file );
   worker->t_prof_file        = FILER_Destroy( worker->t_prof_file );
   worker->q_index_file       = FILER_Destroy( worker->q_index_file );
   worker->t_index_file       = FILER_Destroy( worker->t_index_file );
   worker->mmseqs_file        = FILER_Destroy( worker->mmseqs_file );
   worker->hitlist_file       = FILER_Destroy( worker->hitlist_file );
   /* output files */
   worker->output_file        = FILER_Destroy( worker->output_file );
   worker->tblout_file        = FILER_Destroy( worker->tblout_file );
   worker->m8out_file         = FILER_Destroy( worker->m8out_file );
   worker->myout_file         = FILER_Destroy( worker->myout_file );
   worker->mydomout_file      = FILER_Destroy( worker->mydomout_file );
   worker->mytimeout_file     = FILER_Destroy( worker->mytimeout_file );
   worker->mythreshout_file   = FILER_Destroy( worker->mythreshout_file );

   /* target and profile structures */
   worker->q_seq        = SEQUENCE_Destroy( worker->q_seq );
   worker->t_seq        = SEQUENCE_Destroy( worker->t_seq );
   worker->t_prof       = HMM_PROFILE_Destroy( worker->t_prof );
   worker->hmm_bg       = HMM_BG_Destroy( worker->hmm_bg );
   /* target and profile indexes */
   worker->q_index      = F_INDEX_Destroy( worker->q_index );
   worker->t_index      = F_INDEX_Destroy( worker->t_index );
   /* results in from mmseqs and out for general searches */
   worker->mmseqs_data  = M8_RESULTS_Destroy( worker->mmseqs_data );
   worker->results      = RESULTS_Destroy( worker->results );
   /* free single result */
   ERROR_free( worker->result );
   worker->result       = NULL;
   /* data structs for viterbi alignment */
   worker->trace_vit    = ALIGNMENT_Destroy( worker->trace_vit );
   worker->trace_post   = ALIGNMENT_Destroy( worker->trace_post );
   /* data structs for cloud edgebounds */
   worker->edg_fwd      = EDGEBOUNDS_Destroy( worker->edg_fwd );
   worker->edg_bck      = EDGEBOUNDS_Destroy( worker->edg_bck );
   worker->edg_diag     = EDGEBOUNDS_Destroy( worker->edg_diag );
   worker->edg_row      = EDGEBOUNDS_Destroy( worker->edg_row );
   /* row-wise edgebounds */
   worker->edg_rows_tmp  = EDGEBOUND_ROWS_Destroy( worker->edg_rows_tmp );
   for ( int i=0; i<3; i++ ) {
      worker->lb_vec[i]    = VECTOR_INT_Destroy( worker->lb_vec[i] );
      worker->rb_vec[i]    = VECTOR_INT_Destroy( worker->rb_vec[i] );
   }
   /* necessary dp matrices */
   /* quadratic space */
   // worker->st_MX           = MATRIX_3D_Destroy( worker->st_MX );
   worker->st_MX_fwd       = MATRIX_3D_Destroy( worker->st_MX_fwd );
   worker->st_MX_bck       = MATRIX_3D_Destroy( worker->st_MX_bck );
   // worker->st_MX_post      = MATRIX_3D_Destroy( worker->st_MX_post );
   // worker->st_MX_optacc    = MATRIX_3D_Destroy( worker->st_MX_optacc );
   /* linear space */
   // worker->st_MX3          = MATRIX_3D_Destroy( worker->st_MX3 );
   worker->st_MX3_fwd      = MATRIX_3D_Destroy( worker->st_MX3_fwd );
   worker->st_MX3_bck      = MATRIX_3D_Destroy( worker->st_MX3_bck );
   /* sparse */
   // worker->st_SMX          = MATRIX_3D_SPARSE_Destroy( worker->st_SMX );
   worker->st_SMX_fwd      = MATRIX_3D_SPARSE_Destroy( worker->st_SMX_fwd );
   worker->st_SMX_bck      = MATRIX_3D_SPARSE_Destroy( worker->st_SMX_bck );
   // worker->st_SMX_post     = MATRIX_3D_SPARSE_Destroy( worker->st_SMX_post );
   // worker->st_SMX_optacc   = MATRIX_3D_SPARSE_Destroy( worker->st_SMX_optacc );
   /* special states */
   worker->sp_MX_fwd       = MATRIX_2D_Destroy( worker->sp_MX_fwd );
   worker->sp_MX_bck       = MATRIX_2D_Destroy( worker->sp_MX_bck );
   // worker->sp_MX_post      = MATRIX_2D_Destroy( worker->sp_MX_post );
   // worker->sp_MX_optacc    = MATRIX_2D_Destroy( worker->sp_MX_optacc );
   /* domain definition */
   worker->dom_def         = DOMAIN_DEF_Destroy( worker->dom_def );
}