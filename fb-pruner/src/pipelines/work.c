/*******************************************************************************
 *  FILE:      work.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
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
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../parsers/parsers.h"
#include "../algs_linear/algs_linear.h"
#include "../algs_quad/algs_quad.h"
#include "../algs_naive/algs_naive.h"
#include "../algs_sparse/algs_sparse.h"

/* easel library */
#include "easel.h"
#include "esl_exponential.h"

/* header */
#include "pipelines.h"

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

/* initialize data structs: dynamic programming matrices, edgebounds, etc */
void WORK_init( WORKER* worker )
{
   TASKS*   tasks = worker->tasks;
   ARGS*    args  = worker->args;

   /* initialize logrithmic sum table */
   logsum_Init();

   /* target and profile structures */
   worker->q_seq           = SEQUENCE_Create();
   worker->t_seq           = SEQUENCE_Create();
   worker->t_prof          = HMM_PROFILE_Create();
   worker->hmm_bg          = HMM_BG_Create();
   /* target and profile indexes */
   worker->q_index         = F_INDEX_Create();
   worker->t_index         = F_INDEX_Create();
   /* results in from mmseqs and out for general searches */
   worker->results_in      = RESULTS_Create();
   worker->results         = RESULTS_Create();
   worker->result          = (RESULT*) ERROR_malloc( sizeof(RESULT) );
   /* data structs for viterbi alignment search */
   worker->traceback       = ALIGNMENT_Create();
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
   worker->st_MX_fwd    = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   worker->st_MX_bck    = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   worker->st_MX_post   = MATRIX_3D_Create( NUM_NORMAL_STATES, 1, 1 );
   worker->st_MX        = worker->st_MX_fwd;
   worker->st_MX3_fwd   = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   worker->st_MX3_bck   = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   worker->st_MX3       = worker->st_MX_fwd;
   worker->st_SMX_fwd   = MATRIX_3D_SPARSE_Create();
   worker->st_SMX_bck   = MATRIX_3D_SPARSE_Create();
   // worker->st_SMX_post  = MATRIX_3D_SPARSE_Create();
   worker->st_SMX       = worker->st_SMX_fwd;
   worker->sp_MX_fwd    = MATRIX_2D_Create( NUM_SPECIAL_STATES, 1 );
   worker->sp_MX_bck    = MATRIX_2D_Create( NUM_SPECIAL_STATES, 1 );
   worker->sp_MX_post   = MATRIX_2D_Create( NUM_SPECIAL_STATES, 1 );
   worker->sp_MX        = worker->sp_MX_fwd;
   /* domain definition */
   worker->dom_def      = DOMAIN_DEF_Create();
}

/* clean up data structs */
void WORK_cleanup( WORKER* worker )
{
   TASKS*   tasks       = worker->tasks;
   ARGS*    args        = worker->args;

   /* target and profile structures */
   worker->q_seq        = SEQUENCE_Destroy( worker->q_seq );
   worker->t_seq        = SEQUENCE_Destroy( worker->t_seq );
   worker->t_prof       = HMM_PROFILE_Destroy( worker->t_prof );
   worker->hmm_bg       = HMM_BG_Destroy( worker->hmm_bg );
   /* target and profile indexes */
   worker->q_index      = F_INDEX_Destroy( worker->q_index );
   worker->t_index      = F_INDEX_Destroy( worker->t_index );
   /* results in from mmseqs and out for general searches */
   worker->results_in   = RESULTS_Destroy( worker->results_in );
   worker->results      = RESULTS_Destroy( worker->results );
   /* free single result */
   ERROR_free( worker->result );
   worker->result       = NULL;
   /* data structs for viterbi alignment */
   worker->traceback    = ALIGNMENT_Destroy( worker->traceback );
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
   printf("test\n");
   /* necessary dp matrices */
   worker->st_MX_fwd       = MATRIX_3D_Destroy( worker->st_MX_fwd );
   printf("test\n");
   worker->st_MX_bck       = MATRIX_3D_Destroy( worker->st_MX_bck );
   printf("test\n");
   worker->st_MX_post      = MATRIX_3D_Destroy( worker->st_MX_post );
   printf("test\n");
   worker->st_MX3_fwd      = MATRIX_3D_Destroy( worker->st_MX3_fwd );
   worker->st_MX3_bck      = MATRIX_3D_Destroy( worker->st_MX3_bck );
   printf("test\n");
   worker->st_SMX_fwd      = MATRIX_3D_SPARSE_Destroy( worker->st_SMX_fwd );
   worker->st_SMX_bck      = MATRIX_3D_SPARSE_Destroy( worker->st_SMX_bck );
   worker->sp_MX_fwd       = MATRIX_2D_Destroy( worker->sp_MX_fwd );
   worker->sp_MX_bck       = MATRIX_2D_Destroy( worker->sp_MX_bck );
   worker->sp_MX_bck       = MATRIX_2D_Destroy( worker->sp_MX_post );
   /* domain definition */
   worker->dom_def         = DOMAIN_DEF_Destroy( worker->dom_def );
}

/* open worker files */
void WORK_open( WORKER* worker )
{
   TASKS*   tasks = worker->tasks;
   ARGS*    args  = worker->args;

   /* TODO: need to handle results in bulk. currently report one-at-a-time. */
   /* open file pointers */
   if ( args->is_redirect_stdout ) {
      worker->output_fp    = ERROR_fopen( args->output_filepath, "w" );
   } else {
      worker->output_fp    = stdout;
   }
   if ( args->is_tblout ) {
      worker->tblout_fp    = ERROR_fopen( args->tblout_filepath, "w" );
   }
   if ( args->is_m8out ) {
      worker->m8out_fp     = ERROR_fopen( args->m8out_filepath, "w" );
   }
   if ( args->is_myout ) {
      worker->myout_fp     = ERROR_fopen( args->myout_filepath, "w" );
   }
}

/* close worker files */
void WORK_close( WORKER* worker )
{
   TASKS*   tasks =  worker->tasks;
   ARGS*    args  =  worker->args;  

   /* open file pointers */
   if ( args->is_redirect_stdout ) {
      if ( worker->output_fp != stdout ) {
         worker->output_fp    = ERROR_fclose( worker->output_fp );
      } 
   }
   if ( args->is_tblout ) {
      worker->tblout_fp    = ERROR_fclose( worker->tblout_fp );
   }
   if ( args->is_m8out ) {
      worker->m8out_fp     = ERROR_fclose( worker->m8out_fp );
   }
   if ( args->is_myout ) {
      worker->myout_fp    = ERROR_fclose( worker->myout_fp );
   }
}

/* initialize dynamic programming matrices */
void WORK_reuse( WORKER* worker )
{
   TASKS*   tasks    = worker->tasks; 
   int      T        = worker->t_prof->N;
   int      Q        = worker->q_seq->N;

   /* clear traceback and resize */
   ALIGNMENT_Reuse( worker->traceback, Q, T );
   ALIGNMENT_Reuse( worker->trace_post, Q, T );
   /* clear edgebounds and resize */
   EDGEBOUNDS_Reuse( worker->edg_fwd, Q, T );
   EDGEBOUNDS_Reuse( worker->edg_bck, Q, T );
   EDGEBOUNDS_Reuse( worker->edg_diag, Q, T );
   EDGEBOUNDS_Reuse( worker->edg_row, Q, T );
   /* clear row-wise edgebounds and resize */
   EDGEBOUND_ROWS_Reuse( worker->edg_rows_tmp, Q, T );

   /* matrix for quadratic algs */
   if ( tasks->quadratic ) {
      MATRIX_3D_Reuse_Clean( worker->st_MX_fwd, NUM_NORMAL_STATES,  Q+1, T+1 );
      MATRIX_3D_Reuse_Clean( worker->st_MX_bck, NUM_NORMAL_STATES,  Q+1, T+1 );
      MATRIX_3D_Reuse_Clean( worker->st_MX_post, NUM_NORMAL_STATES,  Q+1, T+1 );
   }
   /* matrix for linear algs */
   if ( tasks->linear ) {
      MATRIX_3D_Reuse_Clean( worker->st_MX3_fwd, NUM_NORMAL_STATES,  3, (Q+1)+(T+1) );
      MATRIX_3D_Reuse_Clean( worker->st_MX3_bck, NUM_NORMAL_STATES,  3, (Q+1)+(T+1) );
      worker->st_MX3 = worker->st_MX3_fwd;
   }
   /* matrices for special states */
   if ( tasks->quadratic || tasks->linear ) {
      MATRIX_2D_Reuse_Clean( worker->sp_MX_fwd, NUM_SPECIAL_STATES, Q+1);
      MATRIX_2D_Reuse_Clean( worker->sp_MX_bck, NUM_SPECIAL_STATES, Q+1);
      MATRIX_2D_Reuse_Clean( worker->sp_MX_post, NUM_SPECIAL_STATES, Q+1);
      worker->sp_MX = worker->sp_MX_fwd;
   }

   /* sparse matrices */
   MATRIX_3D_SPARSE_Reuse( worker->st_SMX_fwd );
   MATRIX_3D_SPARSE_Reuse( worker->st_SMX_bck );
   // MATRIX_3D_SPARSE_Reuse( worker->st_SMX_post );
   worker->st_SMX = worker->st_SMX_fwd;

   #if DEBUG 
   {
      DEBUGGER_Reuse( debugger, Q, T );
   }
   #endif

   /* TODO: remove this? */
   MATRIX_3D_Clean( worker->st_MX );
   MATRIX_3D_Clean( worker->st_MX3 );
}

/* load or build target and query index files */
void WORK_index( WORKER* worker )
{
   /* build file indexes */
   WORK_load_target_index( worker );
   WORK_load_query_index( worker );
   /* sort file indexes */
   F_INDEX_Sort_by_Id( worker->t_index );
   F_INDEX_Sort_by_Id( worker->q_index );
}

/* load or build target and query index files */
void WORK_load_index_by_id( WORKER* worker )
{
   /* build file indexes */
   WORK_load_target_index( worker );
   WORK_load_query_index( worker );
   /* sort file indexes */
   F_INDEX_Sort_by_Id( worker->t_index );
   F_INDEX_Sort_by_Id( worker->q_index );
}

/* load or build target and query index files */
void WORK_load_index_by_name( WORKER* worker )
{
   /* load target index */
   WORK_load_target_index( worker );
   F_INDEX_Sort_by_Name( worker->t_index );
   /* load query index */
   WORK_load_query_index( worker );
   F_INDEX_Sort_by_Name( worker->q_index );
}

/* load target index (or build them if argument missing) */
void WORK_load_target_index(  WORKER*     worker ) 
{
   FILE*    fp       = NULL;

   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;
   TIMES*   times    = worker->times;
   CLOCK*   clok     = worker->clok;

   /* begin time */
   CLOCK_Start(clok);

   /* default index location  (same as main file but with .idx extension) */
   char* t_indexpath_tmp = NULL;
   if (args->t_indexpath == NULL) {
      char* ext = ".idx";
      t_indexpath_tmp = (char*) malloc( sizeof(char) * (strlen(args->t_filepath) + strlen(ext) + 1) );
      if (t_indexpath_tmp == NULL) {
           fprintf(stderr, "ERROR: malloc failed.\n");
       }
      strcpy( t_indexpath_tmp, args->t_filepath );
      strcat( t_indexpath_tmp, ext );
   }

   /* load or build target file index */
   if (args->t_indexpath != NULL) {
      /* load file passed by commandline */
      printf_vhi("# loading indexpath from commandline: '%s'...\n", args->t_indexpath );
      worker->t_index = F_INDEX_Load( worker->t_index, args->t_indexpath );
   }
   else if ( access( t_indexpath_tmp, F_OK ) == 0 ) {
      /* check if standard extension index file exists */
      printf_vhi("# found index at database location: '%s'...\n", t_indexpath_tmp );
      args->t_indexpath = strdup( t_indexpath_tmp );
      worker->t_index = F_INDEX_Load( worker->t_index, t_indexpath_tmp );
   }
   else {
      /* build index on the fly */
      printf("# building index of file...\n");
      if (args->t_filetype == FILE_HMM) {
         worker->t_index = F_INDEX_Hmm_Build( worker->t_index, args->t_filepath );
      }
      else if (args->t_filetype == FILE_FASTA) {
         worker->t_index = F_INDEX_Fasta_Build( worker->t_index, args->t_filepath );
      }
      else {
         fprintf(stderr, "ERROR: target filetype is not supported.\n" );
         exit(EXIT_FAILURE);
      }
      /* identify the query file being indexed */
      worker->t_index->source_path = strdup(args->t_filepath);
      args->t_indexpath = strdup(t_indexpath_tmp);

      /* save index file */
      fp = fopen( t_indexpath_tmp, "w+" );
      F_INDEX_Dump( worker->t_index, fp );
      fclose(fp);
   }
   ERROR_free(t_indexpath_tmp);

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_target_index = CLOCK_Secs(clok);
}

/* load query index (or build them) */
void WORK_load_query_index(   WORKER*     worker ) 
{
   FILE*    fp       = NULL;

   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;
   TIMES*   times    = worker->times;

   CLOCK*   clok    = worker->clok;

   /* begin time */
   CLOCK_Start(clok);

   /* default index location (same as main file but with .idx extenstion) */
   char* q_indexpath_tmp = NULL;
   if (args->q_indexpath == NULL) {
      char* ext = ".idx";
      q_indexpath_tmp = (char*) malloc( sizeof(char) * (strlen(args->q_filepath) + strlen(ext) + 1) );
      if (q_indexpath_tmp == NULL) {
           fprintf(stderr, "ERROR: malloc failed.\n");
       }
      strcpy( q_indexpath_tmp, args->q_filepath );
      strcat( q_indexpath_tmp, ext );
   }

   /* load or build target file index */
   if (args->q_indexpath != NULL) {
      /* load file passed by commandline */
      printf_vhi("# loading indexpath from commandline: '%s'...\n", args->q_indexpath );
      worker->q_index = F_INDEX_Load( worker->q_index, args->q_indexpath );
   } 
   else if ( access( q_indexpath_tmp, F_OK ) == 0 ) {
      /* check if standard extension index file exists */
      printf_vhi("# found index at database location: '%s'...\n", q_indexpath_tmp );
      args->q_indexpath = strdup( q_indexpath_tmp );
      worker->q_index = F_INDEX_Load( worker->q_index, q_indexpath_tmp );
   }  
   else {
      /* build index on the fly */
      printf_vhi("# building index of file...\n");
      if (args->q_filetype == FILE_HMM) {
         worker->q_index = F_INDEX_Hmm_Build( worker->q_index, args->q_filepath );
      }
      else if (args->q_filetype == FILE_FASTA) {
         worker->q_index = F_INDEX_Fasta_Build( worker->q_index, args->q_filepath );
      }
      else {
         fprintf(stderr, "ERROR: query filetype is not supported.\n" );
         exit(EXIT_FAILURE);
      }
      /* identify the query file being indexed */
      worker->q_index->source_path = strdup(args->q_filepath);
      args->q_indexpath = strdup(q_indexpath_tmp);

      /* save index file */
      fp = fopen( q_indexpath_tmp, "w+" );
      F_INDEX_Dump( worker->q_index, fp );
      fclose(fp);
   }
   ERROR_free(q_indexpath_tmp);

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_query_index = CLOCK_Secs(clok);
}

/* build target index */
void WORK_build_target_index(    WORKER*     worker ) 
{
   FILE*    fp       = NULL;

   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;
   TIMES*   times    = worker->times;
   CLOCK*   clok     = worker->clok;

   /* begin time */
   CLOCK_Start(clok);

   /* build index on the fly */
   if ( args->t_filetype == FILE_HMM ) {
      worker->t_index = F_INDEX_Hmm_Build( worker->t_index, args->t_filepath );
   }
   else if (args->t_filetype == FILE_FASTA) {
      worker->t_index = F_INDEX_Fasta_Build( worker->t_index, args->t_filepath );
   }
   else {
      fprintf(stderr, "ERROR: target filetype is not supported.\n" );
      exit(EXIT_FAILURE);
   }
   /* identify the query file being indexed */
   worker->t_index->source_path = strdup(args->t_filepath);

   /* if output location not specified, then use default naming scheme */
   /* default index location  (same as main file but with .idx extension) */
   if ( args->t_indexpath == NULL ) {
      char* t_indexpath_tmp = NULL;
      char* ext = ".idx";
      t_indexpath_tmp = (char*) malloc( sizeof(char) * (strlen(args->t_filepath) + strlen(ext) + 1) );
      if (t_indexpath_tmp == NULL) {
           fprintf(stderr, "ERROR: malloc failed.\n");
      }
      strcpy( t_indexpath_tmp, args->t_filepath );
      strcat( t_indexpath_tmp, ext );
      args->t_indexpath = strdup(t_indexpath_tmp);
      ERROR_free(t_indexpath_tmp);
   }

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_target_index = CLOCK_Secs(clok);
}

/* load query index (or build them) */
void WORK_build_query_index(  WORKER*   worker ) 
{
   FILE*    fp       = NULL;

   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;
   TIMES*   times    = worker->times;

   CLOCK*   clok    = worker->clok;

   /* begin time */
   CLOCK_Start(clok);

   /* build index on the fly */
   if ( args->q_filetype == FILE_HMM ) {
      worker->q_index = F_INDEX_Hmm_Build( worker->q_index, args->q_filepath );
   }
   else if ( args->q_filetype == FILE_FASTA ) {
      worker->q_index = F_INDEX_Fasta_Build( worker->q_index, args->q_filepath );
   }
   else {
      fprintf(stderr, "ERROR: query filetype is not supported.\n" );
      exit(EXIT_FAILURE);
   }
   /* identify the query file being indexed */
   worker->q_index->source_path = strdup(args->q_filepath);

   /* if output location not specified, then use default naming scheme */
   /* default index location  (same as main file but with .idx extension) */
   if ( args->q_indexpath == NULL ) {
      char* q_indexpath_tmp = NULL;
      char* ext = ".idx";
      q_indexpath_tmp = (char*) malloc( sizeof(char) * (strlen(args->q_filepath) + strlen(ext) + 1) );
      if (q_indexpath_tmp == NULL) {
           fprintf(stderr, "ERROR: malloc failed.\n");
      }
      strcpy( q_indexpath_tmp, args->q_filepath );
      strcat( q_indexpath_tmp, ext );
      args->q_indexpath = strdup(q_indexpath_tmp);
      ERROR_free(q_indexpath_tmp);
   }

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_query_index = CLOCK_Secs(clok);
}

/* output target index to file */
void WORK_output_target_index(   WORKER*    worker )
{
   FILE*    fp    = NULL;
   ARGS*    args  = worker->args;

   /* determine the output file to save to */
   if ( args->t_indexpath == NULL ) {
      /* if no output name is given, append ".idx" and save in same directory */
      const char* ext = ".idx";
      args->t_indexpath = (char*) malloc( sizeof(char) * (strlen(args->t_filepath) + strlen(ext) + 1) );
      if ( args->t_indexpath == NULL ) {
         fprintf(stderr, "ERROR: malloc failed.\n");
      }
      strcpy( args->t_indexpath, args->t_filepath );
      strcat( args->t_indexpath, ext );
   }
   /* output target index */
   fp = fopen( args->t_indexpath, "w" );
   F_INDEX_Dump( worker->t_index, fp );
   fclose(fp);
}

/* output query index to file */
void WORK_output_query_index(    WORKER*  worker )
{
   FILE*    fp    = NULL;
   ARGS*    args  = worker->args;

   /* determine the output file to save to */
   if ( args->q_indexpath == NULL ) {
      /* if no output name is given, append ".idx" and save in same directory */
      const char* ext = ".idx";
      args->q_indexpath = (char*) malloc( sizeof(char) * (strlen(args->q_filepath) + strlen(ext) + 1) );
      if ( args->q_indexpath == NULL ) {
         fprintf(stderr, "ERROR: malloc failed.\n");
      }
      strcpy( args->q_indexpath, args->q_filepath );
      strcat( args->q_indexpath, ext );
   }
   /* output query index */
   fp = fopen( args->q_indexpath, "w" );
   F_INDEX_Dump( worker->q_index, fp );
   fclose(fp);
}

/* set and verify ranges */
void WORK_set_ranges(   WORKER*    worker )
{
   ARGS* args = worker->args;

   /* verify range of targets to search */
   /* if beginning is less than zero, default to entire range */
   if ( args->t_range.beg < 0 ) {
      args->t_range = (RANGE) { 0, worker->t_index->N };
   } 
   /* check for invalid ranges ( assert min > max, and min/max are less than the total in file ) */
   else if ( args->t_range.beg > args->t_range.end  || 
           args->t_range.beg > worker->t_index->N || 
           args->t_range.end > worker->t_index->N ) {
      fprintf(stderr, "ERROR: Invalid target id range (%d,%d).\n", args->t_range.beg, args->t_range.end );
      exit(EXIT_FAILURE);
   }

   /* verify range of queries to search */
   /* if beginning is less than zero, default to entire range */
   if ( args->q_range.beg < 0 ) {
      args->q_range = (RANGE) { 0, worker->q_index->N };
   } 
   /* check for invalid ranges ( assert min > max, and min/max are less than the total in file ) */
   else if ( args->q_range.beg > args->q_range.end  || 
           args->q_range.beg > worker->q_index->N || 
           args->q_range.end > worker->q_index->N ) {
      fprintf(stderr, "ERROR: Invalid target id range (%d,%d).\n", args->q_range.beg, args->q_range.end );
      exit(EXIT_FAILURE);
   }
}

/* Establish and verify that result range is valid */
void WORK_load_mmseqs_list( WORKER* worker )
{
   ARGS* args = worker->args;

   /* If range values are negative, then range is set to full list set */
   if ( args->list_range.beg < 0 && args->list_range.end < 0 ) {
      args->list_range.beg = 0;
      args->list_range.end = INT_MAX;
   }

   /* m8+ file contains target_id, query_id, and result_id fields */
   RESULTS_M8_Parse( 
      worker->results_in, args->mmseqs_res_filepath, args->list_range.beg, args->list_range.end );

   /* Truncate or extract valid result range */
   args->list_range.beg = MAX(args->list_range.beg, 0);
   args->list_range.end = MIN(args->list_range.end, worker->results_in->N + args->list_range.beg);
}

/* load target by file index id */
void WORK_load_target_by_id(  WORKER*     worker,
                              int         id )
{
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   CLOCK*         clok           = worker->clok;
   int            index_id       = 0;
   int            index_offset   = 0;

   /* begin time */
   CLOCK_Start(clok);

   /* search and load target by offset */
   index_id = F_INDEX_Search_Id( worker->t_index, id );
   if ( index_id == -1 ) {
      fprintf(stderr, "ERROR: Target id '%d' not found in F_INDEX.\n", id );
      exit(EXIT_FAILURE);
   }
   WORK_load_target_by_index_id( worker, index_id );

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_target = CLOCK_Secs(clok);
}

/* load query by file index id */
void WORK_load_query_by_id(   WORKER*     worker,
                              int         id )
{
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   CLOCK*         clok           = worker->clok;
   int            index_id       = 0;
   int            index_offset   = 0;

   /* begin time */
   CLOCK_Start(clok);

   /* search and load target by offset */
   index_id       = F_INDEX_Search_Id( worker->q_index, id );
   if ( index_id == -1 ) {
      fprintf(stderr, "ERROR: Query id '%d' not found in F_INDEX.\n", id );
      exit(EXIT_FAILURE);
   }
   WORK_load_query_by_index_id( worker, index_id );

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_query = CLOCK_Secs(clok);
}

/* load target by file index name */
void WORK_load_target_by_name(   WORKER*    worker,
                                 char*      name )
{
   // printf_vall("==> WORK_load_target_by_name()\n");
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   TIMES*         t_times        = worker->times_totals;
   CLOCK*         clok           = worker->clok;
   int            index_id       = 0;
   int            index_offset   = 0;

   /* begin time */
   CLOCK_Start(clok);

   /* search and load target by offset */
   index_id = F_INDEX_Search_Name( worker->t_index, name );
   if ( index_id == -1 ) {
      fprintf(stderr, "ERROR: Target name '%s' not found in F_INDEX.\n", name );
      exit(EXIT_FAILURE);
   }
   WORK_load_target_by_index_id( worker, index_id );

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_target = CLOCK_Secs(clok);
   t_times->load_target += times->load_target;
}

/* load target by file index name */
void WORK_load_query_by_name( WORKER*     worker,
                              char*       name )
{
   // printf_vall("==> WORK_load_query_by_name()\n");
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   TIMES*         t_times        = worker->times_totals;
   CLOCK*         clok           = worker->clok;
   int            index_id       = 0;
   int            index_offset   = 0;

   /* begin time */
   CLOCK_Start(clok);

   /* search and load target by offset */
   index_id = F_INDEX_Search_Name( worker->q_index, name );
   if ( index_id == -1 ) {
      fprintf(stderr, "ERROR: Query name '%s' not found in F_INDEX.\n", name );
      exit(EXIT_FAILURE);
   }
   WORK_load_query_by_index_id( worker, index_id);

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_query = CLOCK_Secs(clok);
   t_times->load_query += times->load_query;
}

/* load target by file index id */
void WORK_load_target_by_index_id(  WORKER*     worker,
                                    int         index_id )
{
   ARGS*          args     = worker->args;
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      t_seq    = worker->t_seq;
   SEQUENCE*      q_seq    = worker->q_seq;
   F_INDEX_NODE*  my_idx   = &worker->t_index->nodes[index_id];

   worker->t_id = index_id;

   /* load target profile by file type */
   switch ( args->t_filetype )
   {
      case FILE_HMM: {
         HMM_PROFILE_Parse( worker->t_prof, args->t_filepath, my_idx->offset ); 
         HMM_PROFILE_Convert_NegLog_To_Real( worker->t_prof );
         HMM_PROFILE_Config( worker->t_prof, args->search_mode );
         // HMM_PROFILE_Dump( worker->t_prof, stdout );
      } break;
      case FILE_FASTA: {
         SEQUENCE_Fasta_Parse( worker->t_seq, args->t_filepath, my_idx->offset );
         SEQUENCE_to_HMM_PROFILE( worker->t_seq, worker->t_prof );
         HMM_PROFILE_Dump( worker->t_prof, stdout );
         exit(EXIT_SUCCESS);
      } break;
      default: {
         fprintf(stderr, "ERROR: Only HMM and FASTA filetypes are supported for targets.\n");
         exit(EXIT_FAILURE);
      }
   }
}

/* load target by file index id */
void WORK_load_query_by_index_id(   WORKER*     worker,
                                    int         index_id )
{
   ARGS*          args     = worker->args;
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      t_seq    = worker->t_seq;
   SEQUENCE*      q_seq    = worker->q_seq;
   F_INDEX_NODE*  my_idx   = &worker->q_index->nodes[index_id];

   worker->q_id = index_id;

   /* load query by file type */
   switch ( args->q_filetype )
   {
      /* fasta only supported file type */
      case FILE_FASTA: {
         SEQUENCE_Fasta_Parse( worker->q_seq, args->q_filepath, my_idx->offset );
         // SEQUENCE_Dump( worker->q_seq, stdout );
      } break;
      case FILE_HMM: {

      } 
      default: {
         fprintf(stderr, "ERROR: Only FASTA filetypes are supported for queries.\n");
         exit(EXIT_FAILURE);
      }
   }

   /* set special state transitions based on query sequence length */
   if ( t_prof != NULL ) 
   {
      HMM_PROFILE_ReconfigLength( worker->t_prof, worker->q_seq->N );
   } else {
      fprintf(stderr, "ERROR: Target profile must be loaded before Query Sequence. Currently NULL.\n");
      exit(EXIT_FAILURE);
   }
}

/* viterbi and traceback */
void WORK_viterbi_and_traceback( WORKER*  worker )
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
      times->lin_vit = CLOCK_Secs(clok);
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
      times->quad_vit = CLOCK_Secs(clok);
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
      // times->lin_trace = CLOCK_Secs(clok);
   }
   if ( tasks->quad_trace ) {
      printf_vall("# ==> traceback (quad)...\n");
      CLOCK_Start(clok);
      run_Traceback_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, tr );
      CLOCK_Stop(clok);
      times->quad_trace = CLOCK_Secs(clok);
   }
}

/* forward/backward */
void WORK_forward_backward( WORKER*  worker )
{
   ARGS*          args     = worker->args;
   TASKS*         tasks    = worker->tasks;
   NAT_SCORES*    scores   = worker->scores;
   TIMES*         times    = worker->times;
   CLOCK*         clok     = worker->clok;
   RESULT*        result   = worker->result;

   SEQUENCE*      q_seq    = worker->q_seq;
   HMM_PROFILE*   t_prof   = worker->t_prof;

   int            Q        = q_seq->N;
   int            T        = t_prof->N;

   MATRIX_3D*     st_MX    = worker->st_MX;
   MATRIX_3D*     st_MX3   = worker->st_MX3;
   MATRIX_2D*     sp_MX    = worker->sp_MX;
   ALIGNMENT*     tr       = worker->traceback;
   float          sc       = -1;

   /* forward */
   if ( tasks->lin_fwd ) {
      printf_vall("# ==> forward (lin)...\n");
      CLOCK_Start(clok);
      run_Forward_Linear( q_seq, t_prof, Q, T, worker->st_MX3_fwd, worker->sp_MX_fwd, &sc );
      CLOCK_Stop(clok);
      times->lin_fwd    = CLOCK_Secs(clok);
      scores->lin_fwd   = sc;
      #if DEBUG 
      {
         printf("# printing forward linear...\n");
         printf("# lin forward score: %f\n", scores->lin_fwd);
         DP_MATRIX_Save(Q, T, debugger->test_MX, worker->sp_MX_fwd, "test_output/my.fwd.lin.mx");
         DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, worker->sp_MX_fwd, tr, "test_output/my.fwd.lin.viz.mx");
      }
      #endif
   } 
   if ( tasks->quad_fwd ) {
      printf_vall("# ==> forward (quad)...\n");
      CLOCK_Start(clok);
      run_Forward_Quad( q_seq, t_prof, Q, T, worker->st_MX_fwd, worker->sp_MX_fwd, &sc );
      CLOCK_Stop(clok);
      times->quad_fwd   = CLOCK_Secs(clok);
      scores->quad_fwd  = sc;
      #if DEBUG 
      {
         printf("# printing forward quadratic...\n");
         printf("# quad forward score: %f\n", scores->quad_fwd);
         DP_MATRIX_Save(Q, T, debugger->test_MX, worker->sp_MX_fwd, "test_output/my.fwd.quad.mx");
      }
      #endif
   }

   /* backward */
   if ( tasks->lin_bck ) {
      printf_vall("# ==> backward (lin)...\n");
      CLOCK_Start(clok);
      run_Backward_Linear( q_seq, t_prof, Q, T, worker->st_MX3_bck, worker->sp_MX_bck, &sc );
      CLOCK_Stop(clok);
      times->lin_bck    = CLOCK_Secs(clok);
      scores->lin_bck   = sc;
      #if DEBUG 
      {
         printf("# printing backward linear...\n");
         printf("# lin backward score: %f\n", scores->lin_bck);
         DP_MATRIX_Save(Q, T, debugger->test_MX, worker->sp_MX_bck, "test_output/my.bck.lin.mx");
         DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, worker->sp_MX_bck, tr, "test_output/my.bck.lin.viz.mx");
      }
      #endif
   }
   if ( tasks->quad_bck ) {
      printf_vall("# ==> backward (quad)...\n");
      CLOCK_Start(clok);
      run_Backward_Quad( q_seq, t_prof, Q, T, worker->st_MX_bck, worker->sp_MX_bck, &sc );
      CLOCK_Stop(clok);
      times->quad_bck   = CLOCK_Secs(clok);
      scores->quad_bck  = sc;
      #if DEBUG 
      {
         printf("# printing backward quadratic...\n");
         printf("# quad backward score: %f\n", scores->quad_bck);
         DP_MATRIX_Save(Q, T, debugger->test_MX, worker->sp_MX_bck, "test_output/my.bck.quad.mx");
      }
      #endif
   }
}

/* cloud search */
void WORK_cloud_search( WORKER* worker )
{
   ARGS*             args           = worker->args;
   TASKS*            tasks          = worker->tasks;

   TIMES*            times          = worker->times;
   TIMES*            t_times        = worker->times_totals;
   CLOCK*            clok           = worker->clok;
   RESULT*           result         = worker->result;

   NAT_SCORES*       scores         = worker->scores;

   SEQUENCE*         q_seq          = worker->q_seq;
   HMM_PROFILE*      t_prof         = worker->t_prof;

   F_INDEX*          q_index        = worker->q_index;
   F_INDEX*          t_index        = worker->t_index;

   int               Q              = q_seq->N;
   int               T              = t_prof->N;

   float             alpha          = args->alpha;
   float             beta           = args->beta;
   int               gamma          = args->gamma;
   float             sc             = 0.0;

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
   CLOUD_PARAMS*     cloud_params   = &(worker->cloud_params);

   int status;

   /* if performing linear fb-pruner, run cloud search  */
   if ( tasks->lin_bound_fwd || tasks->lin_bound_bck ) 
   {
      /* cloud forward */
      printf_vall("# ==> cloud forward (linear)...\n");
      CLOCK_Start(clok);
      run_Cloud_Forward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX, tr, edg_rows_tmp, edg_fwd, cloud_params );
      CLOCK_Stop(clok);
      times->lin_cloud_fwd = CLOCK_Secs(clok);
      t_times->lin_cloud_fwd += times->lin_cloud_fwd;
      #if DEBUG
      {
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, "test_output/my.cloud_fwd.lin.mx");
      }
      #endif

      /* cloud backward */
      printf_vall("# ==> cloud backward (linear)...\n");
      CLOCK_Start(clok);
      run_Cloud_Backward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX, tr, edg_rows_tmp, edg_bck, cloud_params );
      CLOCK_Stop(clok);
      times->lin_cloud_bck = CLOCK_Secs(clok);
      t_times->lin_cloud_bck += times->lin_cloud_bck ;
      #if DEBUG
      {
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, "test_output/my.cloud_bck.lin.mx");
      }
      #endif

      /* merge edgebounds */
      printf_vall("# ==> merge (linear)...\n");
      CLOCK_Start(clok);
      EDGEBOUNDS_Merge_Together( Q, T, edg_fwd, edg_bck, edg_diag );
      CLOCK_Stop(clok);
      times->lin_merge = CLOCK_Secs(clok);
      t_times->lin_merge += times->lin_merge;

      /* reorient edgebounds */
      printf_vall("# ==> reorient (linear)...\n");
      CLOCK_Start(clok);
      int precount  = EDGEBOUNDS_Count( edg_row );
      EDGEBOUNDS_Reorient_to_Row( Q, T, edg_diag, edg_row );
      CLOCK_Stop(clok);
      times->lin_reorient = CLOCK_Secs(clok);
      t_times->lin_reorient += times->lin_reorient;

      /* compute the number of cells in matrix computed */
      result->cloud_cells = EDGEBOUNDS_Count( edg_row );
      result->total_cells  = (Q+1) * (T+1);

      /* if debugging, print the cloud */
      #if DEBUG
      {
         // MATRIX_2D_Fill( debugger->cloud_MX, 0 );
         // MATRIX_2D_Cloud_Fill( debugger->cloud_MX, edg_fwd, 1 );
         // MATRIX_2D_Cloud_Fill( debugger->cloud_MX, edg_bck, 2 );
         // DP_MATRIX_VIZ_Trace( debugger->cloud_MX, tr );
         // DP_MATRIX_VIZ_Color_Dump( debugger->cloud_MX, stdout );
         // DP_MATRIX_VIZ_Dump( debugger->cloud_MX, stdout );
         printf("# printing cloud rows...\n");
         /* TEMPORARY OVERRIDE */
         EDGEBOUNDS_Cover_Matrix(edg_row, Q, T);
         EDGEBOUNDS_Save( edg_row, "test_output/my.cloud.rows.edg");
      }
      #endif

      /* compute the number of cells in matrix computed */
      result->cloud_cells  = EDGEBOUNDS_Count( edg_row );
      result->total_cells  = (Q+1) * (T+1);
   }
    
   /* linear bounded forward */
   if ( tasks->lin_bound_fwd ) 
   {
      printf_vall("# ==> bound forward (linear)...\n");
      CLOCK_Start(clok);
      run_Bound_Forward_Linear( q_seq, t_prof, Q, T, st_MX3_fwd, sp_MX_fwd, edg_row, &sc );
      scores->lin_bound_fwd = sc;
      CLOCK_Stop(clok);

      #if DEBUG
      {
         printf("# printing linear bound forward...\n");
         printf("# lin bound forward: %f\n", scores->lin_bound_fwd );
         printf("# cells => total: %d, cloud: %d, perc: %f\n", result->total_cells, result->cloud_cells, (float)result->cloud_cells/(float)result->total_cells );
         // DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX, "test_output/my.bound_fwd.lin.mx");
         // DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, sp_MX_fwd, tr, "test_output/my.bound_fwd.lin.viz.mx");
      }
      #endif
      times->lin_bound_fwd    = CLOCK_Secs(clok);
      t_times->lin_bound_fwd  += times->lin_bound_fwd;
      scores->lin_bound_fwd   = sc;

      times->lin_fbpruner_total  =  times->lin_cloud_fwd + times->lin_cloud_bck  + 
                                    times->lin_merge     + times->lin_reorient   +
                                    times->lin_bound_fwd + times->lin_bound_bck;
   }
   /* linear bounded backward */
   if ( tasks->lin_bound_bck ) 
   {
      printf_vall("# ==> bound backward (linear)...\n");
      CLOCK_Start(clok);
      run_Bound_Backward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX_bck, edg_row, &sc );
      scores->lin_bound_bck = sc;
      CLOCK_Stop(clok);
      #if DEBUG
      {
         printf_vall("# printing linear bound backward...\n");
         printf("# lin bound backward: %f\n", scores->lin_bound_bck );
         // DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_bck, "test_output/my.bound_bck.lin.mx");
         // DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, sp_MX, tr, "test_output/my.bound_bck.lin.viz.mx");
      }
      #endif
      times->lin_bound_bck    = CLOCK_Secs(clok);
      scores->lin_bound_bck   = sc;

      /* TODO: Remove this!!! */
      DP_MATRIX_Clean( Q, T, st_MX3, sp_MX );
   }

   /* if performing sparse fb-pruner, initialize sparse matrix */
   if (tasks->sparse_bound_fwd || tasks->sparse_bound_bck)
   {
      MATRIX_3D_SPARSE_Shape_Like_Edgebounds( st_SMX_fwd, edg_row );
      MATRIX_3D_SPARSE_Copy( st_SMX_bck, st_SMX_fwd );
   }

   /* sparse bounded forward */
   if ( tasks->sparse_bound_fwd ) 
   {
      printf_vall("# ==> bound forward (sparse)...\n");
      CLOCK_Start(clok);
      run_Bound_Forward_Sparse( q_seq, t_prof, Q, T, st_SMX_fwd, sp_MX_fwd, edg_row, &sc );
      scores->sparse_bound_fwd = sc;
      CLOCK_Stop(clok);

      /* compute the number of cells in matrix computed */
      result->cloud_cells = EDGEBOUNDS_Count( edg_row );
      result->total_cells  = (Q+1) * (T+1);

      #if DEBUG
      {
         printf("# printing sparse bound forward...\n");
         printf("# sparse bound forward: %f\n", scores->sparse_bound_fwd );
         printf("# cells => total: %d, cloud: %d, perc: %f\n", result->total_cells, result->cloud_cells, (float)result->cloud_cells/(float)result->total_cells );
         // DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_fwd, "test_output/my.bound_fwd.sp.mx");
         // DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, sp_MX_fwd, tr, "test_output/my.bound_fwd.sp.viz.mx");
      }
      #endif
      times->sp_bound_fwd        = CLOCK_Secs(clok);
      t_times->sp_bound_fwd      += times->lin_bound_fwd;
      scores->sparse_bound_fwd   = sc;

      times->sp_fbpruner_total   =  times->lin_cloud_fwd + times->lin_cloud_bck  + 
                                    times->lin_merge     + times->lin_reorient   +
                                    times->sp_bound_fwd + times->sp_bound_bck;
   }
   /* sparse bounded backward */
   if ( tasks->sparse_bound_bck ) 
   {
      printf_vall("# ==> bound backward (sparse)...\n");
      CLOCK_Start(clok);
      run_Bound_Backward_Sparse( q_seq, t_prof, Q, T, st_SMX_bck, sp_MX_bck, edg_row, &sc );
      scores->sparse_bound_bck = sc;
      CLOCK_Stop(clok);
      #if DEBUG
      {
         printf_vall("# printing sparse bound backward...\n");
         printf("# sparse bound backward: %f\n", scores->sparse_bound_bck );
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_bck, "test_output/my.bound_bck.sp.mx");
         DP_MATRIX_Trace_Save(Q, T, debugger->test_MX, sp_MX_bck, tr, "test_output/my.bound_bck.sp.viz.mx");
      }
      #endif
      times->sp_bound_bck        = CLOCK_Secs(clok);
      scores->sparse_bound_bck   = sc;
   }

   /* if performing quadratic bounded forward or backward  */
   if ( tasks->quad_bound_fwd || tasks->quad_bound_bck ) 
   {
      /* cloud forward */
      CLOCK_Start(clok);
      run_Cloud_Forward_Quad( 
         q_seq, t_prof, Q, T, st_MX, sp_MX, tr, edg_rows_tmp, lb_vec, rb_vec, edg_fwd, cloud_params );
      CLOCK_Stop(clok);
      times->quad_cloud_fwd = CLOCK_Secs(clok);

      /* cloud backward */
      CLOCK_Start(clok);
      run_Cloud_Backward_Quad( 
         q_seq, t_prof, Q, T, st_MX, sp_MX, tr, edg_rows_tmp, lb_vec, rb_vec, edg_bck, cloud_params );
      CLOCK_Stop(clok);
      times->quad_cloud_bck = CLOCK_Secs(clok);

      /* merge edgebounds */
      CLOCK_Start(clok);
      EDGEBOUNDS_Merge_Together( Q, T, edg_fwd, edg_bck, edg_diag );
      CLOCK_Stop(clok);
      times->quad_merge = CLOCK_Secs(clok);

      /* reorient edgebounds */
      CLOCK_Start(clok);
      EDGEBOUNDS_Reorient_to_Row( Q, T, edg_diag, edg_row );
      CLOCK_Stop(clok);
      times->quad_merge = CLOCK_Secs(clok);
   }
   /* quadratic bounded forward */
   if ( tasks->quad_bound_fwd ) 
   {
      printf_vall("# ==> bound forward (quad)...\n");
      CLOCK_Start(clok);
      run_Bound_Forward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, edg_row, &sc );
      CLOCK_Stop(clok);
      times->quad_bound_fwd   = CLOCK_Secs(clok);
      scores->quad_bound_fwd  = sc;

      times->quad_fbpruner_total =     times->quad_cloud_fwd   + times->quad_cloud_bck + 
                                       times->quad_merge       +
                                       times->quad_bound_fwd   + times->quad_bound_bck;
   }
   /* quadratic bounded backward */
   if ( tasks->quad_bound_bck ) 
   {
      printf_vall("# ==> bound backward (quad)...\n");
      CLOCK_Start(clok);
      run_Bound_Backward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, edg_row, &sc );
      CLOCK_Stop(clok);
      times->quad_bound_bck   = CLOCK_Secs(clok);
      scores->quad_bound_bck  = sc;
   }
}

/* use sparse matrix to cover alignment */
void WORK_capture_alignment( WORKER* worker )
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
      // fp = fopen("test_output/sparse_mx.tbl", "w");
      // MATRIX_3D_SPARSE_Dump( worker->st_SMX_fwd, fp );
      // fclose(fp);
      // fp = fopen("test_output/sparse_mx.inner.edg", "w");
      // EDGEBOUNDS_Dump( worker->st_SMX_fwd->edg_inner, fp );
      // fclose(fp);
      // fp = fopen("test_output/sparse_mx.outer.edg", "w");
      // EDGEBOUNDS_Dump( worker->st_SMX_fwd->edg_outer, fp );
      // fclose(fp);
      // fp = fopen("test_output/sparse_mx.fullmap.tbl", "w");
      // MATRIX_3D_SPARSE_Map_to_Inner_Dump( worker->st_SMX_fwd, NULL, fp );
      // fclose(fp);

      MATRIX_2D_Fill( debugger->cloud_MX, 0.0 );
      MATRIX_2D_Cloud_Fill( debugger->cloud_MX, worker->st_SMX_fwd->edg_outer, 1.0 );
      MATRIX_2D_Cloud_Fill( debugger->cloud_MX, worker->st_SMX_fwd->edg_inner, 0.1 );
      fp = fopen("test_output/sparse_mx.mx", "w");
      MATRIX_2D_Dump( debugger->cloud_MX, fp );
      fclose(fp);
   }
   #endif

   /* run forward */
   run_Bound_Forward_Sparse( 
      worker->q_seq, worker->t_prof, q_seq->N, t_prof->N, worker->st_SMX_fwd, worker->sp_MX_fwd, worker->edg_row, &sc);
   scores->sparse_bound_fwd = sc;
   #if DEBUG
   {
      printf("# printing sparse bound forward...\n");
      // DP_MATRIX_Copy( Q, T, debugger->test_MX, worker->sp_MX, st_MX_dst, sp_MX_dst );
      // DP_MATRIX_Save( Q, T, debugger->test_MX, worker->sp_MX, "test_output/my.bound_fwd.sparse.mx");
      DP_MATRIX_Trace_Save( Q, T, debugger->test_MX, worker->sp_MX, tr, "test_output/my.bound_fwd.sparse.viz.mx");
   }
   #endif

   /* run backward */
   run_Bound_Backward_Sparse( 
      worker->q_seq, worker->t_prof, q_seq->N, t_prof->N, worker->st_SMX_bck, worker->sp_MX_bck, worker->edg_row, &sc );
   scores->sparse_bound_bck = sc;
   #if DEBUG
   {
      printf("# printing sparse bound backward...\n");
      // DP_MATRIX_Save( Q, T, debugger->test_MX, worker->sp_MX, "test_output/my.bound_bck.sparse.mx" );
      DP_MATRIX_Save( Q, T, debugger->test_MX, worker->sp_MX, "test_output/my.bound_bck.sparse.mx" );
      DP_MATRIX_Trace_Save( Q, T, debugger->test_MX, worker->sp_MX, tr, "test_output/my.bound_bck.sparse.viz.mx" );
   }
   #endif

   /* Recover alignment */
   // run_Posterior_Sparse( 
   //    q_seq, t_prof, q_seq->N, t_prof->N, worker->st_SMX_fwd, worker->sp_MX, worker->edg_row, worker->trace_post );
   // ALIGNMENT_Dump(worker->trace_post, stdout);

   /* generate alignments */
   // ALIGNMENT_Build_MMSEQS_Style( 
   //    worker->traceback, worker->q_seq, worker->t_prof );
   // ALIGNMENT_Build_HMMER_Style( 
   //    worker->traceback, worker->q_seq, worker->t_prof );
}

/* compute correction bias and convert natscore -> bitscore -> pval -> eval */
void WORK_convert_scores( WORKER* worker )
{
   HMM_BG*        bg       = worker->hmm_bg;
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      q_seq    = worker->q_seq;
   RESULT*        result   = worker->result;
   TASKS*         tasks    = worker->tasks;

   int   T           = worker->t_prof->N;
   int   Q           = worker->q_seq->N;
   /* alignment scores */
   float nat_sc      = 0.0f;  /* score in NATS */
   float pre_sc      = 0.0f;  /* adjusted for compo bias / score in BITS */
   float seq_sc      = 0.0f;  /* adjusted for sequence bias / score in BITS */
   float ln_pval     = 0.0f;  /* natural log of p-value */
   float pval        = 0.0f;  /* p-value */
   float eval        = 0.0f;  /* e-value */
   /* bias correction */
   float null_sc     = 0.0f;  /* in NATS */
   float filter_sc   = 0.0f;  /* in NATS */
   float seq_bias    = 0.0f;  /* in NATS */
   /* parameters for exponential distribution, for converting bitscore -> p-value */
   float tau         = worker->t_prof->forward_dist.param1;
   float lambda      = worker->t_prof->forward_dist.param2;
   /* number of sequences in database, for computing e-value */
   int   n_seqs      = worker->q_index->N;

   /* initialize hmm_bg */
   HMM_BG_SetSequence( bg, q_seq );
   HMM_BG_SetFilter( bg, t_prof->N, t_prof->bg_model->compo );
   HMM_BG_SetLength( bg, q_seq->N );
   /* compute null one */
   HMM_BG_NullOne( bg, q_seq->N, &null_sc );
   /* compute nullscore for bias */
   HMM_BG_FilterScore( bg, q_seq, &filter_sc );
   /* free digitized sequence TODO: move to sequence */
   // HMM_BG_UnsetSequence( bg, seq );

   /* composition bias */
   seq_bias = 0.0f;
   if ( worker->args->is_compo_bias )
   {
      /* Find domains and assesses domain-specific correction bias */
      /* see p7_domaindef_ByPosteriorHeuristics() */

      if ( tasks->quad_bias_corr == true )
      {
         run_Posterior_Quad(
            worker->q_seq, worker->t_prof, worker->q_seq->N, worker->t_prof->N, worker->hmm_bg, worker->edg_row,
            worker->st_MX_fwd, worker->sp_MX_fwd, worker->st_MX_bck, worker->sp_MX_bck, worker->st_MX_bck, worker->sp_MX_bck, 
            worker->dom_def );
      }
      if ( tasks->sparse_bias_corr == true )
      {
         run_Posterior_Sparse(
            worker->q_seq, worker->t_prof, worker->q_seq->N, worker->t_prof->N, worker->hmm_bg, worker->edg_row,
            worker->st_SMX_fwd, worker->sp_MX_fwd, worker->st_SMX_bck, worker->sp_MX_bck, worker->st_SMX_bck, worker->sp_MX_bck, 
            worker->dom_def );
      }

      
      /* compute sequence bias */
      seq_bias = logsum(0.0, worker->hmm_bg->omega + worker->dom_def->seq_bias);
      printf("nat_sc: %7.3f, null_sc: %7.3f, (final) seq_bias: %9.5f,\n", nat_sc, null_sc, seq_bias);
      seq_bias = worker->dom_def->seq_bias;
   }

   /* get cloud forward score */
   nat_sc = worker->scores->lin_bound_fwd;
   /* compute pre_score and sequence_score by accounting for bias and convert from nats -> bits */
   pre_sc = (nat_sc - null_sc) / eslCONST_LOG2;
   seq_sc = (nat_sc - (null_sc + seq_bias)) / eslCONST_LOG2;
   printf("# nat_sc = %7.4f, null_sc = %7.4f, seq_bias = %7.4f, pre_sc = %7.4f, seq_sc = %7.4f\n",
      nat_sc, null_sc, seq_bias, pre_sc, seq_sc );

   /* compute log of P-value */ 
   ln_pval  = esl_exp_logsurv( seq_sc, tau, lambda );
   pval     = exp(ln_pval);
   eval     = pval * n_seqs;
   // printf("# ln_pval = %7.4f, pval = %9.2e, eval = %9.2e\n",
   //    ln_pval, pval, eval );

   /* save scores */
   result->final_scores.nat_sc      = nat_sc;
   result->final_scores.null_sc     = null_sc;
   result->final_scores.filter_sc   = filter_sc;
   result->final_scores.seq_bias    = seq_bias;
   result->final_scores.pre_sc      = pre_sc;
   result->final_scores.seq_sc      = seq_sc;
   result->final_scores.pre_sc      = pre_sc;
   result->final_scores.ln_pval     = ln_pval;
   result->final_scores.pval        = pval;
   result->final_scores.eval        = eval;
}

/* print header for results file (default) */
void WORK_report_header( WORKER* worker )
{
   ARGS* args = worker->args;

   /* open file pointers */
   REPORT_stdout_header( worker, worker->output_fp );

   if ( args->is_tblout ) {
      REPORT_tblout_header( worker, worker->tblout_fp );
   }
   if ( args->is_m8out ) {
      REPORT_m8out_header( worker, worker->m8out_fp );
   }
   if ( args->is_myout ) {
      REPORT_myout_header( worker, worker->myout_fp );
   }
}

/* print current result (default) */
void WORK_report_result_current( WORKER* worker )
{
   ARGS* args = worker->args;

   /* open file pointers */
   REPORT_stdout_entry( worker, worker->result, worker->output_fp );
   if ( args->is_tblout ) {
      REPORT_tblout_entry( worker, worker->result, worker->tblout_fp );
   }
   if ( args->is_m8out ) {
      REPORT_m8out_entry( worker, worker->result, worker->m8out_fp );
   }
   if ( args->is_myout ) { 
      REPORT_myout_entry( worker, worker->result, worker->myout_fp );
   }
}

/* print all results */
void WORK_report_result_all( WORKER* worker )
{
   ARGS* args = worker->args;

}

/* initial output printed before search */
/* for verbose output */
void WORK_report_footer( WORKER*  worker ) 
{
   ARGS* args = worker->args;

   /* open file pointers */
   REPORT_stdout_footer( worker, worker->output_fp );

   if ( args->is_tblout ) {
      REPORT_tblout_footer( worker, worker->tblout_fp );
   }
   if ( args->is_m8out ) {
      REPORT_m8out_footer( worker, worker->m8out_fp );
   }
   if ( args->is_myout ) { 
      REPORT_myout_footer( worker, worker->myout_fp );
   }
}

