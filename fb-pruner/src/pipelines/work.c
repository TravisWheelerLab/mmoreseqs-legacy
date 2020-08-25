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
#include "structs.h"
#include "utilities.h"
#include "objects.h"
#include "parsers.h"
#include "algs_linear.h"
#include "algs_quad.h"
#include "algs_naive.h"
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
   TASKS* tasks = worker->tasks;

   /* initialize logrithmic sum table */
   logsum_Init();

   /* target and profile structures */
   worker->q_seq        = SEQUENCE_Create();
   worker->t_seq        = SEQUENCE_Create();
   worker->t_prof       = HMM_PROFILE_Create();
   worker->hmm_bg       = HMM_BG_Create();
   /* target and profile indexes */
   worker->q_index      = F_INDEX_Create();
   worker->t_index      = F_INDEX_Create();
   /* results in from mmseqs and out for general searches */
   worker->results_in   = RESULTS_Create();
   worker->results      = RESULTS_Create();
   worker->result       = (RESULT*) malloc( sizeof(RESULT) );
   /* data structs for viterbi alignment search */
   worker->traceback    = ALIGNMENT_Create();
   worker->trace_post   = ALIGNMENT_Create();
   /* data structs for cloud edgebounds */
   worker->edg_fwd      = EDGEBOUNDS_Create();
   worker->edg_bck      = EDGEBOUNDS_Create();
   worker->edg_diag     = EDGEBOUNDS_Create();
   worker->edg_row      = EDGEBOUNDS_Create();
   /* row-wise edgebounds */
   worker->edg_rows_tmp  = EDGEBOUND_ROWS_Create();
   for ( int i=0; i<3; i++ ) {
      worker->lb_vec[i] = VECTOR_INT_Create();
      worker->rb_vec[i] = VECTOR_INT_Create();
   }
   /* cloud search parameters */
   worker->cloud_params.alpha    = worker->args->alpha;
   worker->cloud_params.beta     = worker->args->beta;
   worker->cloud_params.gamma    = worker->args->gamma;

   /* create necessary dp matrices */
   worker->st_MX  = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   worker->st_MX3 = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   worker->st_SMX = MATRIX_3D_SPARSE_Create();
   worker->sp_MX  = MATRIX_2D_Create( NUM_SPECIAL_STATES, 1 );

   #if DEBUG
   {
      debugger->cloud_MX   = MATRIX_2D_Create( 1, 1 ); 
      debugger->cloud_MX3  = MATRIX_2D_Create( 1, 1 ); 
      debugger->test_MX    = MATRIX_3D_Create( NUM_NORMAL_STATES, 1, 1 );
      debugger->test_MX3   = MATRIX_3D_Create( NUM_NORMAL_STATES, 1, 1 );

      debugger->test_edg = EDGEBOUNDS_Create();
   }
   #endif
}

/* clean up data structs */
void WORK_cleanup( WORKER* worker )
{
   TASKS* tasks = worker->tasks;

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
   ERRORCHECK_free( worker->result );
   worker->result = NULL;
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
      worker->lb_vec[i] = VECTOR_INT_Destroy( worker->lb_vec[i] );
      worker->rb_vec[i] = VECTOR_INT_Destroy( worker->rb_vec[i] );
   }
   /* create necessary dp matrices */
   worker->st_MX = MATRIX_3D_Destroy( worker->st_MX );
   worker->st_MX3 = MATRIX_3D_Destroy( worker->st_MX3 );
   worker->st_SMX = MATRIX_3D_SPARSE_Destroy( worker->st_SMX );
   worker->sp_MX  = MATRIX_2D_Destroy( worker->sp_MX );

   #if DEBUG
   {
      debugger->cloud_MX   = MATRIX_2D_Destroy( debugger->cloud_MX ); 
      debugger->cloud_MX3  = MATRIX_2D_Destroy( debugger->cloud_MX3 ); 
      debugger->test_MX    = MATRIX_3D_Destroy( debugger->test_MX );
      debugger->test_MX3   = MATRIX_3D_Destroy( debugger->test_MX3 );

      debugger->test_edg   = EDGEBOUNDS_Destroy( debugger->test_edg );
   }
   #endif
}

/* initialize dynamic programming matrices */
void WORK_reuse( WORKER* worker )
{
   TASKS*   tasks    = worker->tasks; 

   int   T  = worker->t_prof->N;
   int   Q  = worker->q_seq->N;

   /* clear traceback and resize */
   ALIGNMENT_Reuse( worker->traceback, Q, T );
   /* clear edgebounds and resize */
   EDGEBOUNDS_Reuse( worker->edg_fwd, Q, T );
   EDGEBOUNDS_Reuse( worker->edg_bck, Q, T );
   EDGEBOUNDS_Reuse( worker->edg_diag, Q, T );
   EDGEBOUNDS_Reuse( worker->edg_row, Q, T );
   /* clear row-wise edgebounds and resize */
   EDGEBOUND_ROWS_Reuse( worker->edg_rows_tmp, Q, T );

   /* reuse necessary dp matrices (only reallocs if new size is larger) */
   if ( tasks->quadratic ) {
      MATRIX_3D_Reuse_Clean( worker->st_MX, NUM_NORMAL_STATES,  Q+1, T+1 );
   }
   if ( tasks->linear ) {
      MATRIX_3D_Reuse_Clean( worker->st_MX3, NUM_NORMAL_STATES,  3, (Q+1)+(T+1) );
   }
   if ( tasks->quadratic || tasks->linear ) {
      MATRIX_2D_Reuse_Clean( worker->sp_MX, NUM_SPECIAL_STATES, Q+1);
   }

   #if DEBUG 
   {
      MATRIX_2D_Reuse_Clean( debugger->cloud_MX, Q+1, T+1 );
      MATRIX_2D_Reuse_Clean( debugger->cloud_MX3, 3, (Q+1)+(T+1) );
      MATRIX_3D_Reuse_Clean( debugger->test_MX, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
      MATRIX_3D_Reuse_Clean( debugger->test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
   }
   #endif

   /* TODO: remove this */
   MATRIX_3D_Clean( worker->st_MX );
   MATRIX_3D_Clean( worker->st_MX3 );
   MATRIX_2D_Clean( worker->sp_MX );
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

/* load target index (or build them if argument missing) */
void WORK_load_target_index( WORKER* worker ) 
{
   FILE*    fp       = NULL;

   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;
   TIMES*   times    = worker->times;
   CLOCK*   clok     = worker->clok;

   /* begin time */
   CLOCK_Start(clok);

   /* defuault index location  (same as main file but with .idx extension) */
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
      printf_vhi("loading indexpath from commandline...\n");
      worker->t_index = F_INDEX_Load( worker->t_index, args->t_indexpath );
   }
   else if ( access( t_indexpath_tmp, F_OK ) == 0 ) {
      /* check if standard extension index file exists */
      printf_vhi("found index at database location...\n");
      worker->t_index = F_INDEX_Load( worker->t_index, t_indexpath_tmp );
   }
   else {
      /* build index on the fly */
      printf("building index of file...\n");
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
   ERRORCHECK_free(t_indexpath_tmp);

   /* if we have a mmseqs tmp file location, and index is not already using mmseqs names */
   if ( tasks->mmseqs_lookup && worker->t_index->mmseqs_names ) {
      printf_vhi("updating query with mmseqs lookup...\n");
      if ( args->t_lookup_filepath != NULL ) {
         /* if filepath given directly as argument, do nothing */
      }
      else if ( args->mmseqs_tmp_filepath != NULL ) {
         /* else, construct default path to lookup from mmseqs tmp folder */
         char* t_lookup_path     = "/latest/target.lookup"; 
         args->t_lookup_filepath = (char*) malloc( sizeof(char) * ( strlen(args->mmseqs_tmp_filepath) + strlen(t_lookup_path) + 1 ) );
         strcpy( args->t_lookup_filepath, args->mmseqs_tmp_filepath );
         strcat( args->t_lookup_filepath, t_lookup_path );
      }
      /*  */
      F_INDEX_Lookup_Update( worker->t_index, args->t_lookup_filepath );
   }

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_target_index = CLOCK_Secs(clok);
}

/* load query index (or build them) */
void WORK_load_query_index( WORKER* worker ) 
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
      printf_vhi("loading indexpath from commandline...\n");
      worker->q_index = F_INDEX_Load( worker->q_index, args->q_indexpath );
   } 
   else if ( access( q_indexpath_tmp, F_OK ) == 0 ) {
      /* check if standard extension index file exists */
      printf_vhi("found index at database location...\n");
      args->q_indexpath = strdup( q_indexpath_tmp );
      worker->q_index = F_INDEX_Load( worker->q_index, q_indexpath_tmp );
   }  
   else {
      /* build index on the fly */
      printf_vhi("building index of file...\n");
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
   ERRORCHECK_free(q_indexpath_tmp);

   /* if we have a mmseqs tmp file location, and index is not already using mmseqs names */
   if ( tasks->mmseqs_lookup && worker->q_index->mmseqs_names ) {
      printf_vhi("updating query with mmseqs lookup...\n");
      if ( args->q_lookup_filepath != NULL ) {
         /* if filepath given directly as argument, do nothing */
      }
      else if ( args->mmseqs_tmp_filepath != NULL ) {
         /* else, construct default path to lookup from mmseqs tmp folder */
         char* q_lookup_path     = "/latest/query.lookup"; 
         args->q_lookup_filepath = (char*) malloc( sizeof(char) * ( strlen(args->mmseqs_tmp_filepath) + strlen(q_lookup_path) + 1 ) );
         strcpy( args->q_lookup_filepath, args->mmseqs_tmp_filepath );
         strcat( args->q_lookup_filepath, q_lookup_path );
      }
      /*  */
      F_INDEX_Lookup_Update( worker->q_index, args->q_lookup_filepath );
   }

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_query_index = CLOCK_Secs(clok);
}

/* output target index to file */
void WORK_output_target_index( WORKER* worker )
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
void WORK_output_query_index( WORKER* worker )
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
void WORK_set_ranges( WORKER* worker )
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

/* load target by file index id */
void WORK_load_target_by_id( WORKER* worker,
                             int     id )
{
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   CLOCK*         clok           = worker->clok;
   RESULTS*       results        = worker->results;
   RESULT*        result         = worker->result;

   HMM_PROFILE*   t_prof         = worker->t_prof;
   SEQUENCE*      t_seq          = worker->t_seq;

   F_INDEX*       t_index        = worker->t_index;
   long           t_offset       = -1;

   char*          t_filepath     = args->t_filepath;
   int            t_filetype     = args->t_filetype;
   int            mode           = args->search_mode;

   /* set current target id */
   worker->t_id = id;
   printf("# loading target by id: %d..\n", id);

   /* begin time */
   CLOCK_Start(clok);

   /* get offset into by checking index */
   int            term;
   F_INDEX_NODE*  node;
   term  = F_INDEX_Search_Id(t_index, id);
   node = &(t_index->nodes[ term ]);
   t_offset = node->offset;

   /* print results */
   printf("# target_id: %d, result_id: %d, node_id: %d, node_name: %s, offset: %ld\n", 
      id, term, node->id, node->name, t_offset );
   F_INDEX_Node_Dump( t_index, term, stdout );

   /* load target profile by file type */
   switch ( t_filetype )
   {
      case FILE_HMM:
         HMM_PROFILE_Parse( t_prof, t_filepath, t_offset ); 
         HMM_PROFILE_Convert_NegLog_To_Real( t_prof );
         HMM_PROFILE_Config( t_prof, mode );
         // printf("=== HMM PROFILE ===\n");
         // HMM_PROFILE_Dump( t_prof, stdout );
         break;
      case FILE_FASTA:
         SEQUENCE_Fasta_Parse( t_seq, t_filepath, t_offset );
         SEQUENCE_to_HMM_PROFILE( t_seq, t_prof );
         // HMM_PROFILE_Dump( t_prof, stdout );
         break;
      default:
         fprintf(stderr, "ERROR: Only HMM and FASTA filetypes are supported for t_profs.\n");
         exit(EXIT_FAILURE);
   }

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_target = CLOCK_Secs(clok);
}

/* load query by file index id */
void WORK_load_query_by_id( WORKER* worker,
                            int     id )
{
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   CLOCK*         clok           = worker->clok;
   RESULTS*       results        = worker->results;
   RESULT*        result         = worker->result;

   HMM_PROFILE*   t_prof         = worker->t_prof;
   SEQUENCE*      q_seq          = worker->q_seq;

   F_INDEX*       q_index        = worker->q_index;
   long           q_offset       = -1;

   char*          q_filepath     = args->q_filepath;
   int            q_filetype     = args->q_filetype;

   /* set current query id */
   worker->q_id = id;
   printf_vlo("# loading query by id: %d..\n", id);

   /* begin time */
   CLOCK_Start(clok);

   /* get offset into by checking index */
   int            term;
   F_INDEX_NODE*  node;
   term  = F_INDEX_Search_Id( q_index, id );
   node = &(q_index->nodes[term]);
   q_offset = node->offset;

   /* print results */
   printf_vlo("# index_size: %d, target_id: %d, result_id: %d, node_id: %d, node_name: %s, offset: %ld\n", 
      q_index->N, id, term, node->id, node->name, q_offset );
   F_INDEX_Node_Dump( q_index, term, stdout );
   printf("test test, %d ? %d\n", q_filetype, FILE_FASTA);

   /* load query by file type */
   switch ( q_filetype )
   {
      /* fasta only supported file type */
      case FILE_FASTA: {
         printf("=== SEQUENCE ===\n");
         SEQUENCE_Fasta_Parse( q_seq, q_filepath, q_offset );
         // SEQUENCE_Dump( q_seq, stdout );
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
      HMM_PROFILE_ReconfigLength( t_prof, q_seq->N );
   }
   else
   {
      fprintf(stderr, "ERROR: Target profile must be loaded before Query Sequence. Currently NULL.\n");
      exit(EXIT_FAILURE);
   }

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_query = CLOCK_Secs(clok);
}

/* load target by file index name */
void WORK_load_target_by_name( WORKER* worker,
                               char*   name )
{
   int t_id = F_INDEX_Search_Name( worker->t_index, name );
   if ( t_id != -1 ) {
      fprintf(stderr, "ERROR: Target name '%s' not found in F_INDEX.\n", name );
      exit(EXIT_FAILURE);
   }
   WORK_load_target_by_id( worker, t_id );
}

/* load target by file index name */
void WORK_load_query_by_name( WORKER* worker,
                              char*   name )
{
   int q_id = F_INDEX_Search_Name( worker->q_index, name );
   if ( q_id != -1 ) {
      fprintf(stderr, "ERROR: Query name '%s' not found in F_INDEX.\n", name );
      exit(EXIT_FAILURE);
   }
   WORK_load_query_by_id( worker, q_id );
}

/* viterbi and traceback */
/* TODO complete this */
void WORK_viterbi_and_traceback( WORKER*  worker )
{
   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;
   SCORES*  scores   = worker->scores;
   TIMES*   times    = worker->times;
   CLOCK*   clok    = worker->clok;

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
      fprintf(stderr, "ERROR: Operation not supported.\n");
      exit(EXIT_FAILURE);
      // TODO: linear viterbi goes here
      // CLOCK_Start(clok);
      // viterbi_Lin( q_seq, t_prof, Q, T, st_MX, sp_MX, quad_sc );
      // CLOCK_Stop(clok);
      // times->lin_vit = CLOCK_Secs(clok);
   }
   if ( tasks->quad_vit ) {
      printf_vall("=> viterbi (quad)...\n");
      CLOCK_Start(clok);
      run_Viterbi_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, &sc );
      CLOCK_Stop(clok);
      times->quad_vit = CLOCK_Secs(clok);
      scores->quad_vit = sc;
   }

   /* Traceback */
   if ( tasks->lin_trace ) {
      fprintf(stderr, "ERROR: Operation not supported.\n");
      exit(EXIT_FAILURE);
      // TODO: linear traceback goes here
      // CLOCK_Start(clok);
      // run_Traceback_Quad(q_seq, t_prof, Q, T, st_MX, sp_MX, tr);
      // CLOCK_Stop(clok);
      // times->lin_trace = CLOCK_Secs(clok);
   }
   if ( tasks->quad_trace ) {
      printf_vall("=> traceback (quad)...\n");
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
   SCORES*        scores   = worker->scores;
   TIMES*         times    = worker->times;
   CLOCK*         clok    = worker->clok;
   RESULT*        result   = worker->result;

   SEQUENCE*      q_seq    = worker->q_seq;
   HMM_PROFILE*   t_prof   = worker->t_prof;

   int   Q  = q_seq->N;
   int   T  = t_prof->N;

   MATRIX_3D*     st_MX    = worker->st_MX;
   MATRIX_3D*     st_MX3   = worker->st_MX3;
   MATRIX_2D*     sp_MX    = worker->sp_MX;
   ALIGNMENT*     tr       = worker->traceback;
   float          sc       = -1;

   /* forward */
   if ( tasks->lin_fwd ) {
      printf_vall("=> forward (lin)...\n");
      CLOCK_Start(clok);
      run_Backward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX, &sc );
      CLOCK_Stop(clok);
      times->lin_fwd    = CLOCK_Secs(clok);
      scores->lin_fwd   = sc;
   } 
   if ( tasks->quad_fwd ) {
      printf_vall("=> forward (quad)...\n");
      CLOCK_Start(clok);
      run_Forward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, &sc );
      CLOCK_Stop(clok);
      times->quad_fwd   = CLOCK_Secs(clok);
      scores->quad_fwd  = sc;
      #if DEBUG 
      {
         printf("# printing forward quadratic...\n");
         DP_MATRIX_Trace_Save(Q, T, st_MX, sp_MX, tr, "test_output/my.fwd.quad.mx");
      }
      #endif
   }

   /* backward */
   if ( tasks->lin_bck ) {
      printf_vall("=> backward (lin)...\n");
      CLOCK_Start(clok);
      run_Backward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX, &sc );
      CLOCK_Stop(clok);
      times->lin_bck    = CLOCK_Secs(clok);
      scores->lin_bck   = sc;
   }
   if ( tasks->quad_bck ) {
      printf_vall("=> backward (quad)...\n");
      CLOCK_Start(clok);
      run_Forward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, &sc );
      CLOCK_Stop(clok);
      times->quad_bck   = CLOCK_Secs(clok);
      scores->quad_bck  = sc;
      #if DEBUG 
      {
         DP_MATRIX_Trace_Save(Q, T, st_MX, sp_MX, tr, "test_output/my.bck.quad.mx");
      }
      #endif
   }
}

/* cloud search */
void WORK_cloud_search( WORKER* worker )
{
   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;

   SCORES*  scores   = worker->scores;
   SCORES*  pvals    = worker->pvals;
   SCORES*  evals    = worker->evals;

   TIMES*   times    = worker->times;
   CLOCK*   clok     = worker->clok;
   RESULT*  result   = worker->result;

   SEQUENCE*      q_seq    = worker->q_seq;
   HMM_PROFILE*   t_prof   = worker->t_prof;

   F_INDEX*       q_index  = worker->q_index;
   F_INDEX*       t_index  = worker->t_index;

   int   Q  = q_seq->N;
   int   T  = t_prof->N;

   float    alpha       = args->alpha;
   float    beta        = args->beta;
   int      gamma       = args->gamma;
   float    sc;

   MATRIX_3D*     st_MX       = worker->st_MX;
   MATRIX_3D*     st_MX3      = worker->st_MX3;
   MATRIX_2D*     sp_MX       = worker->sp_MX;
   MATRIX_3D*     st_cloud_MX = worker->st_cloud_MX; 

   ALIGNMENT*     tr          = worker->traceback;

   EDGEBOUNDS*    edg_fwd     = worker->edg_fwd;
   EDGEBOUNDS*    edg_bck     = worker->edg_bck;
   EDGEBOUNDS*    edg_diag    = worker->edg_diag;
   EDGEBOUNDS*    edg_row     = worker->edg_row;

   EDGEBOUND_ROWS*   edg_rows_tmp   = worker->edg_rows_tmp;
   VECTOR_INT**      lb_vec         = worker->lb_vec;
   VECTOR_INT**      rb_vec         = worker->rb_vec;
   CLOUD_PARAMS*     cloud_params   = &(worker->cloud_params);

   int status;

   /* if performing linear bounded forward or backward  */
   if ( tasks->lin_bound_fwd || tasks->lin_bound_bck ) {
      /* cloud forward */
      printf_vall("=> cloud forward (lin)...\n");
      CLOCK_Start(clok);
      run_Cloud_Forward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX, tr, edg_rows_tmp, edg_fwd, cloud_params );
      CLOCK_Stop(clok);
      times->lin_cloud_fwd = CLOCK_Secs(clok);
      DP_MATRIX_Clean( Q, T, st_MX3, sp_MX );

      /* cloud backward */
      printf_vall("=> cloud backward (lin)...\n");
      CLOCK_Start(clok);
      run_Cloud_Backward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX, tr, edg_rows_tmp, edg_bck, cloud_params );
      CLOCK_Stop(clok);
      times->lin_cloud_bck = CLOCK_Secs(clok);
      DP_MATRIX_Clean( Q, T, st_MX3, sp_MX );

      /* merge edgebounds */
      printf_vall("=> merge (lin)...\n");
      CLOCK_Start(clok);
      EDGEBOUNDS_Merge_Together( Q, T, edg_fwd, edg_bck, edg_diag );
      CLOCK_Stop(clok);
      times->lin_merge = CLOCK_Secs(clok);

      /* reorient edgebounds */
      printf_vall("=> reorient (lin)...\n");
      CLOCK_Start(clok);
      int precount  = EDGEBOUNDS_Count( edg_row );
      EDGEBOUNDS_Reorient_to_Row( Q, T, edg_diag, edg_row );
      CLOCK_Stop(clok);
      times->lin_reorient = CLOCK_Secs(clok);

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
         EDGEBOUNDS_Save( edg_row, "test_output/my.cloud.rows.edg");
      }
      #endif

      /* compute the number of cells in matrix computed */
      result->cloud_cells  = EDGEBOUNDS_Count( edg_row );
      result->total_cells  = (Q+1) * (T+1);
   }
   /* bounded forward */
   if ( tasks->lin_bound_fwd ) {
      printf_vall("=> bound forward (lin)...\n");
      CLOCK_Start(clok);
      run_Bound_Forward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX, edg_row, &sc );
      CLOCK_Stop(clok);
      #if DEBUG
      {
         MATRIX_3D* st_MX_lin = debugger->test_MX;
         printf("# printing linear bound forward...\n");
         DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX, tr, "test_output/my.bound_fwd.lin.mx");
      }
      #endif
      times->lin_bound_fwd = CLOCK_Secs(clok);
      scores->lin_cloud_fwd = sc;

      times->lin_total_cloud =   times->lin_cloud_fwd + times->lin_cloud_bck  + 
                                 times->lin_merge     + times->lin_reorient   +
                                 times->lin_bound_fwd + times->lin_bound_bck;
   }

   /* bounded backward */
   if ( tasks->lin_bound_bck ) {
      printf_vall("=> bound backward (lin)...\n");
      CLOCK_Start(clok);
      run_Bound_Backward_Linear( q_seq, t_prof, Q, T, st_MX3, sp_MX, edg_row, &sc );
      CLOCK_Stop(clok);
      times->lin_bound_bck = CLOCK_Secs(clok);
      scores->lin_cloud_bck = sc;

      DP_MATRIX_Clean( Q, T, st_MX3, sp_MX );
   }

   /* if performing quadratic bounded forward or backward  */
   if ( tasks->quad_bound_fwd || tasks->quad_bound_bck ) {
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
      
      /* compute the number of cells in matrix computed */
      result->cloud_cells = EDGEBOUNDS_Count( edg_row );
      result->total_cells  = (Q+1) * (T+1);
   }

   /* bounded forward */
   if ( tasks->quad_bound_fwd ) {
      printf_vall("=> bound forward (quad)...\n");
      CLOCK_Start(clok);
      run_Bound_Forward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, edg_row, &sc );
      CLOCK_Stop(clok);
      times->quad_bound_fwd = CLOCK_Secs(clok);
      scores->quad_cloud_fwd = sc;

      times->quad_total_cloud =  times->quad_cloud_fwd   + times->quad_cloud_bck + 
                                 times->quad_merge +
                                 times->quad_bound_fwd   + times->quad_bound_bck;
   }
   /* bounded backward */
   if ( tasks->quad_bound_bck ) {
      printf_vall("=> bound backward (quad)...\n");
      CLOCK_Start(clok);
      run_Bound_Backward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, edg_row, &sc );
      CLOCK_Stop(clok);
      times->quad_bound_bck = CLOCK_Secs(clok);
      scores->quad_cloud_bck = sc;
   }
}

/* compute correction bias and convert natscore -> bitscore -> pval -> eval */
void WORK_convert_scores( WORKER* worker )
{
   HMM_BG*        bg    = worker->hmm_bg;
   HMM_PROFILE*   prof  = worker->t_prof;
   SEQUENCE*      seq   = worker->q_seq;

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
   HMM_BG_SetSequence( bg, seq );
   HMM_BG_SetFilter( bg, prof->N, prof->bg_model->compo );
   HMM_BG_SetLength( bg, seq->N );
   /* compute null one */
   HMM_BG_NullOne( bg, seq->N, &null_sc );
   /* compute nullscore for bias */
   HMM_BG_FilterScore( bg, seq, &filter_sc );
   /* TODO: compute sequence bias? (requires null2 for seq/domain specific correction) */
   /* see p7_domaindef_ByPosteriorHeuristics() */
   seq_bias  = 0.0f;

   /* get cloud forward score */
   nat_sc = worker->scores->lin_cloud_fwd;
   /* compute pre_score and sequence_score by accounting for bias and convert from nats -> bits */
   pre_sc = (nat_sc - null_sc) / eslCONST_LOG2;
   seq_sc = (nat_sc - (null_sc + seq_bias)) / eslCONST_LOG2;
   // printf("# nat_sc = %7.4f, null_one = %7.4f, null_sc = %7.4f, pre_sc = %7.4f, seq_sc = %7.4f\n",
   //    nat_sc, null_sc, filter_sc, pre_sc, seq_sc );

   /* compute log of P-value */ 
   ln_pval  = esl_exp_logsurv( seq_sc, tau, lambda );
   pval     = exp(ln_pval);
   eval     = pval * n_seqs;
   // printf("# ln_pval = %7.4f, pval = %9.2e, eval = %9.2e\n",
   //    ln_pval, pval, eval );

   /* save e-value */
   worker->pvals->lin_cloud_fwd = pval;
   worker->evals->lin_cloud_fwd = eval;
}

/* TODO: fill result header based on task flags */
void WORK_result_header_fill( WORKER* worker )
{
   FILE*    fp       = NULL;

   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;
   SCORES*  scores   = worker->scores;
   TIMES*   times    = worker->times;
   RESULTS* results  = worker->results;

   /* open file */
   fp = fopen(args->output_filepath, "w");
   if ( fp == NULL ) {
      fprintf( stderr, "ERROR: Bad File Pointer => '%s'\n", args->output_filepath );
      exit(EXIT_FAILURE);
   }

   /* header */
   fprintf(fp, ">");
   /* query/target */
   fprintf(fp, "{%s}\t", "Q_ID");
   fprintf(fp, "{%s}\t", "Q_NAME");

   fprintf(fp, "{%s}\t", "T_ID");
   fprintf(fp, "{%s}\t", "T_NAME");

   if ( tasks->scores ) {
      if ( tasks->quadratic ) {
         if ( tasks->quad_fwd ) fprintf(fp, "{%s}\t", "QUAD_FWD_SCORE");
         if ( tasks->quad_bck ) fprintf(fp, "{%s}\t", "QUAD_BCK_SCORE");
         if ( tasks->quad_vit ) fprintf(fp, "{%s}\t", "QUAD_VIT_SCORE");
         if ( tasks->quad_bound_fwd ) fprintf(fp, "{%s}\t", "QUAD_BND_FWD_SCORE");
         if ( tasks->quad_bound_bck ) fprintf(fp, "{%s}\t", "QUAD_BND_BCK_SCORE");
      } 
      if ( tasks->linear ) {
         if ( tasks->lin_fwd ) fprintf(fp, "{%s}\t", "LIN_FWD_SCORE");
         if ( tasks->lin_bck ) fprintf(fp, "{%s}\t", "LIN_BCK_SCORE");
         if ( tasks->lin_vit ) fprintf(fp, "{%s}\t", "LIN_VIT_SCORE");
         if ( tasks->lin_bound_fwd ) fprintf(fp, "{%s}\t", "LIN_BND_FWD_SCORE");
         if ( tasks->lin_bound_bck ) fprintf(fp, "{%s}\t", "LIN_BND_BCK_SCORE");
      }
   }
   if ( tasks->time ) {
      if ( tasks->quadratic ) {
         if ( tasks->quad_fwd ) fprintf(fp, "{%s}\t", "QUAD_FWD_TIME");
         if ( tasks->quad_bck ) fprintf(fp, "{%s}\t", "QUAD_BCK_TIME");
         if ( tasks->quad_vit ) fprintf(fp, "{%s}\t", "QUAD_VIT_TIME");
         if ( tasks->quad_trace ) fprintf(fp, "{%s}\t", "QUAD_TRACE_TIME");
         if ( tasks->quad_bound_fwd ) fprintf(fp, "{%s}\t", "QUAD_BND_FWD_TIME");
         if ( tasks->quad_bound_bck ) fprintf(fp, "{%s}\t", "QUAD_BND_BCK_TIME");
      } 
      if ( tasks->linear ) {
         if ( tasks->lin_fwd ) fprintf(fp, "{%s}\t", "LIN_FWD_TIME");
         if ( tasks->lin_bck ) fprintf(fp, "{%s}\t", "LIN_BCK_TIME");
         if ( tasks->lin_vit ) fprintf(fp, "{%s}\t", "LIN_VIT_TIME");
         if ( tasks->lin_trace ) fprintf(fp, "{%s}\t", "LIN_TRACE_TIME");
         if ( tasks->lin_bound_fwd ) fprintf(fp, "{%s}\t", "LIN_BND_FWD_TIME");
         if ( tasks->lin_bound_bck ) fprintf(fp, "{%s}\t", "LIN_BND_BCK_TIME");
      }
   }
   fprintf(fp, "\n");

   fclose(fp);
}

/* TODO: fill result with data based on task flags */
void WORK_result_fill( WORKER* worker )
{
   SCORES*  scores   = worker->scores;
   TIMES*   times    = worker->times;
   RESULT*  result   = worker->result;
}

/* print header for results file (default) */
void WORK_print_result_header( WORKER* worker )
{
   FILE*    fp    = worker->out_file;
   ARGS*    args  = worker->args;

   /* open file */
   fprintf(fp, "# TARGET_SOURCE:\t%s\n", args->t_filepath);
   fprintf(fp, "# QUERY_SOURCE:\t%s\n",  args->q_filepath);
   fprintf(fp, "# TARGET_INDEX:\t%s\n",  args->t_indexpath);
   fprintf(fp, "# QUERY_INDEX:\t%s\n",   args->q_indexpath);

   /* print header */
   int pad = 0;
   fprintf(fp, ">");

   fprintf(fp, "{%*s}\t", pad, "T_ID" );
   fprintf(fp, "{%*s}\t", pad, "Q_ID" );
   // fprintf(fp, "{%*s}\t", pad, "T_NAME" );
   // fprintf(fp, "{%*s}\t", pad, "Q_NAME" );
   fprintf(fp, "{%*s}\t", pad, "T_LEN" );
   fprintf(fp, "{%*s}\t", pad, "Q_LEN" );

   fprintf(fp, "{%*s}\t", pad, "ALPHA" );
   fprintf(fp, "{%*s}\t", pad, "BETA" );

   fprintf(fp, "{%*s}\t", pad, "TOTAL_CNT" );
   fprintf(fp, "{%*s}\t", pad, "CLOUD_CNT" );

   fprintf(fp, "{%*s}\t", pad, "VIT_SC" );
   fprintf(fp, "{%*s}\t", pad, "CLOUD_SC" );

   fprintf(fp, "{%*s}\t", pad, "VIT_T" );
   fprintf(fp, "{%*s}\t", pad, "TRACE_T" );
   fprintf(fp, "{%*s}\t", pad, "CL_FWD_T" );
   fprintf(fp, "{%*s}\t", pad, "CL_BCK_T" );
   fprintf(fp, "{%*s}\t", pad, "MERGE_T" );
   fprintf(fp, "{%*s}\t", pad, "REORIENT_T" );
   fprintf(fp, "{%*s}\t", pad, "BND_FWD_T" );

   fprintf(fp, "\n");
}

/* print current result (default) */
void WORK_print_result_current( WORKER* worker )
{
   FILE*    fp       = worker->out_file;
   ARGS*    args     = worker->args;
   TIMES*   times    = worker->times;
   SCORES*  scores   = worker->scores;

   int      name_pad = 50;

   fprintf(fp, "%d\t",     worker->t_id );
   fprintf(fp, "%d\t",     worker->q_id );
   // fprintf(fp, "%.*s\t",   name_pad, worker->t_prof->name );
   // fprintf(fp, "%.*s\t",   name_pad, worker->q_seq->name );
   fprintf(fp, "%d\t",     worker->t_prof->N );
   fprintf(fp, "%d\t",     worker->q_seq->N );

   fprintf(fp, "%.1f\t",    args->alpha );
   fprintf(fp, "%.1f\t",    args->beta );
   fprintf(fp, "%d\t",      args->gamma );

   fprintf(fp, "%d\t",     worker->result->total_cells );
   fprintf(fp, "%d\t",     worker->result->cloud_cells );

   fprintf(fp, "%f\t",     scores->quad_vit );
   fprintf(fp, "%f\t",     scores->lin_cloud_fwd );

   fprintf(fp, "%f\t",     times->quad_vit );
   fprintf(fp, "%f\t",     times->quad_trace );
   fprintf(fp, "%f\t",     times->lin_cloud_fwd );
   fprintf(fp, "%f\t",     times->lin_cloud_bck );
   fprintf(fp, "%f\t",     times->lin_merge );
   fprintf(fp, "%f\t",     times->lin_reorient );
   fprintf(fp, "%f\t",     times->lin_bound_fwd );

   fprintf(fp, "\n");
}

/* initial output printed before search */
/* for verbose output */
void WORK_print_header( WORKER*  worker ) 
{

}