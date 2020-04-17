/*******************************************************************************
 *  FILE:      work.c
 *  PURPOSE:   Pipelines Workflow Subroutines
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
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

/* data structures */
#include "objects/structs.h"
#include "utilities/utility.h"
#include "error_handler.h"

/* file parsers */
#include "parsers/arg_parser.h"
#include "parsers/hmm_parser.h"
#include "parsers/seq_parser.h"
#include "parsers/m8_parser.h"
#include "parsers/index_parser.h"
#include "parsers/seq_to_profile.h"

/* objects */
#include "objects/f_index.h"
#include "objects/results.h"
#include "objects/alignment.h"
#include "objects/sequence.h"
#include "objects/hmm_profile.h"
#include "objects/edgebound.h"
#include "objects/clock.h"
#include "objects/matrix/matrix_2d.h"
#include "objects/matrix/matrix_3d.h"
#include "objects/mystring.h"
#include "objects/vectors/vector_range_2d.h"

/* viterbi & fwdbck (quadratic) */
#include "algs_quad/viterbi_quad.h"
#include "algs_quad/traceback_quad.h"
#include "algs_quad/fwdback_quad.h"
/* viterbi & fwdbck (linear) */
#include "algs_linear/fwdback_linear.h"

/* cloud search (naive) */
#include "algs_naive/bounded_fwdbck_naive.h"
/* cloud search (quadratic space) */
#include "algs_quad/cloud_search_quad.h"
#include "algs_quad/merge_reorient_quad.h"
#include "algs_quad/bounded_fwdbck_quad.h"
/* cloud search (linear space) */
#include "algs_linear/cloud_search_linear.h"
#include "algs_linear/merge_reorient_linear.h"
#include "algs_linear/bounded_fwdbck_linear.h"
/* temp test */
#include "algs_linear/cloud_search_linear_rows.h"

/* debugging methods */
#include "testing.h"

/* header */
#include "pipeline.h"

/* generic workflow */
void WORK_workflow( WORKER*  work )
{

}

/* initial output printed before search */
/* for verbose output */
void WORK_print_header( WORKER*  worker ) 
{

}

/* initialize dynamic programming matrices, edgebounds,  */
void WORK_init( WORKER* worker )
{
   TASKS*   tasks    = worker->tasks;

   /* initialize logrithmic sum table */
   init_Logsum();

   /* target and profile structures */
   worker->q_seq        = SEQUENCE_Create();
   worker->t_seq        = SEQUENCE_Create();
   worker->t_prof       = HMM_PROFILE_Create();

   /* results in from mmseqs and out for general searches */
   worker->results_in   = RESULTS_Create();
   worker->results      = RESULTS_Create();

   /* data structs for cloud search */
   worker->traceback    = ALIGNMENT_Create();

   worker->edg_fwd      = EDGEBOUNDS_Create();
   worker->edg_bck      = EDGEBOUNDS_Create();
   worker->edg_diag     = EDGEBOUNDS_Create();
   worker->edg_row      = EDGEBOUNDS_Create();

   /* create necessary dp matrices */
   if ( tasks->quadratic ) {
      worker->st_MX = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   }
   if ( tasks->linear ) {
      worker->st_MX3 = MATRIX_3D_Create( NUM_NORMAL_STATES,  1, 1 );
   }
   if ( tasks->quadratic || tasks->linear ) {
      worker->sp_MX  = MATRIX_2D_Create( NUM_SPECIAL_STATES, 1 );
   }
}

/* initialize dynamic programming matrices */
void WORK_reuse( WORKER* worker )
{
   TASKS*   tasks    = worker->tasks; 

   int   T  = worker->t_prof->N;
   int   Q  = worker->q_seq->N;

   /* clear traceback and edgebound data */
   ALIGNMENT_Reuse( worker->traceback );
   
   EDGEBOUNDS_Reuse( worker->edg_fwd );
   EDGEBOUNDS_Reuse( worker->edg_bck );
   EDGEBOUNDS_Reuse( worker->edg_diag );
   EDGEBOUNDS_Reuse( worker->edg_row );

   /* reuse necessary dp matrices (only reallocs if new size is larger) */
   if ( tasks->quadratic ) {
      MATRIX_3D_Reuse( worker->st_MX, NUM_NORMAL_STATES,  Q+1, T+1 );
   }
   if ( tasks->linear ) {
      MATRIX_3D_Reuse( worker->st_MX3, NUM_NORMAL_STATES,  3, (Q+1)+(T+1) );
   }
   if ( tasks->quadratic || tasks->linear ) {
      MATRIX_2D_Reuse( worker->sp_MX, NUM_SPECIAL_STATES, Q+1);
   }
}

/* load target index (or build them) */
void WORK_load_target_index( WORKER* worker ) 
{
   FILE*    fp       = NULL;

   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;
   TIMES*   times    = worker->times;
   CLOCK*   clock    = worker->clock;

   /* begin time */
   CLOCK_Start(clock);

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
      printf("loading indexpath from commandline...\n");
      worker->t_index = F_INDEX_Load(args->t_indexpath);
   }
   else if ( access( t_indexpath_tmp, F_OK ) == 0 ) {
      /* check if standard extension index file exists */
      printf("found index at database location...\n");
      worker->t_index = F_INDEX_Load( t_indexpath_tmp );
   }
   else {
      /* build index on the fly */
      printf("building index of file...\n");
      if (args->t_filetype == FILE_HMM) {
         worker->t_index = F_INDEX_Hmm_Build(args->t_filepath);
      }
      else if (args->t_filetype == FILE_FASTA) {
         worker->t_index = F_INDEX_Fasta_Build(args->t_filepath);
      }
      else {
         fprintf(stderr, "ERROR: target filetype is not supported.\n" );
         exit(EXIT_FAILURE);
      }
      /* identify the query file being indexed */
      worker->t_index->source_path = strdup(args->t_filepath);

      /* save index file */
      args->t_indexpath = t_indexpath_tmp;
      fp = fopen( t_indexpath_tmp, "w+" );
      F_INDEX_Dump( worker->t_index, fp );
      fclose(fp);
   }

   /* if we have a mmseqs tmp file location, and index is not already using mmseqs names */
   if ( tasks->mmseqs_lookup && worker->q_index->mmseqs_names ) {
      printf("updating target with mmseqs lookup...\n");
      char* t_lookup_path  = "/latest/target.lookup"; 
      char* t_lookup       = (char*)malloc( sizeof(char) * ( strlen(args->mmseqs_tmp_filepath) + strlen(t_lookup_path) + 1 ) );
      strcpy( t_lookup, args->mmseqs_tmp_filepath );
      strcat( t_lookup, t_lookup_path );
      F_INDEX_Lookup_Update( worker->t_index, t_lookup );
      free( t_lookup );
   }

   /* end and save time */
   CLOCK_Stop(clock);
   times->load_target_index = CLOCK_Secs(clock);
}

/* load query index (or build them) */
void WORK_load_query_index( WORKER* worker ) 
{
   FILE*    fp       = NULL;

   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;
   TIMES*   times    = worker->times;

   CLOCK*   clock    = worker->clock;

   /* begin time */
   CLOCK_Start(clock);

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
      printf("loading indexpath from commandline...\n");
      worker->q_index = F_INDEX_Load(args->q_indexpath);
   } 
   else if ( access( q_indexpath_tmp, F_OK ) == 0 ) {
      /* check if standard extension index file exists */
      printf("found index at database location...\n");
      args->q_indexpath = strdup( q_indexpath_tmp );
      worker->q_index = F_INDEX_Load( q_indexpath_tmp );
   }  
   else {
      /* build index on the fly */
      printf("building index of file...\n");
      if (args->q_filetype == FILE_HMM) {
         worker->q_index = F_INDEX_Hmm_Build(args->q_filepath);
      }
      else if (args->q_filetype == FILE_FASTA) {
         worker->q_index = F_INDEX_Fasta_Build(args->q_filepath);
      }
      else {
         fprintf(stderr, "ERROR: query filetype is not supported.\n" );
         exit(EXIT_FAILURE);
      }
      /* identify the query file being indexed */
      worker->q_index->source_path = strdup(args->q_filepath);
      args->q_indexpath = strdup( q_indexpath_tmp );

      /* save index file */
      fp = fopen( q_indexpath_tmp, "w+" );
      F_INDEX_Dump( worker->q_index, fp );
      fclose(fp);
   }
   free(q_indexpath_tmp);

   /* if we have a mmseqs tmp file location, and index is not already using mmseqs names */
   if ( tasks->mmseqs_lookup && worker->q_index->mmseqs_names ) {
      printf("updating query with mmseqs lookup...\n");
      /* construct default path to lookup from mmseqs tmp folder */
      char* q_lookup_path  = "/latest/query.lookup"; 
      char* q_lookup       = (char*)malloc( sizeof(char) * ( strlen(args->mmseqs_tmp_filepath) + strlen(q_lookup_path) + 1 ) );
      strcpy( q_lookup, args->mmseqs_tmp_filepath );
      strcat( q_lookup, q_lookup_path );
      F_INDEX_Lookup_Update( worker->q_index, q_lookup );
      free( q_lookup );
   }

   /* end and save time */
   CLOCK_Stop(clock);
   times->load_query_index = CLOCK_Secs(clock);
}

/* output target index to file */
void WORK_output_target_index( WORKER* worker )
{
   FILE*    fp    = NULL;
   ARGS*    args  = worker->args;

   /* determine the method of output for target index */
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
   // F_INDEX_Sort_by_Name( worker->t_index );
   fclose(fp);
}

/* output query index to file */
void WORK_output_query_index( WORKER* worker )
{
   FILE*    fp    = NULL;
   ARGS*    args  = worker->args;

   /* determine the method of output for query index */
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
   // F_INDEX_Sort_by_Name( worker->q_index );
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
   CLOCK*         clock          = worker->clock;
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

   /* begin time */
   CLOCK_Start(clock);

   /* add id to result report */
   result->target_id = id;

   /* get offset into by checking index */
   t_offset = t_index->nodes[id].offset;

   /* load target profile */
   if ( t_filetype == FILE_HMM ) 
   {
      HMM_PROFILE_Parse( t_prof, t_filepath, t_offset ); 
      HMM_PROFILE_Convert_NegLog_To_Real( t_prof );
      HMM_PROFILE_Config( t_prof, mode );
   }
   else if ( t_filetype == FILE_FASTA )
   {
      fprintf(stderr, "ERROR: Fasta file not currently supported.\n");
      exit(EXIT_FAILURE);

      SEQUENCE_Fasta_Parse( t_seq, t_filepath, t_offset );
      SEQUENCE_to_HMM_PROFILE( t_seq, t_prof );
   }
   else
   {
      fprintf(stderr, "ERROR: Only HMM and FASTA filetypes are supported for t_profs.\n");
      exit(EXIT_FAILURE);
   }
   // HMM_PROFILE_ReconfigLength( t_prof, q_seq->N );

   /* end and save time */
   CLOCK_Start(clock);
   times->load_target = CLOCK_Secs(clock);
}

/* load query by file index id */
void WORK_load_query_by_id( WORKER* worker,
                            int     id )
{
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   CLOCK*         clock          = worker->clock;
   RESULTS*       results        = worker->results;
   RESULT*        result         = worker->result;

   SEQUENCE*      q_seq          = worker->q_seq;

   F_INDEX*       q_index        = worker->q_index;
   long           q_offset       = -1;

   char*          q_filepath     = args->q_filepath;
   int            q_filetype     = args->q_filetype;

   /* set current query id */
   worker->q_id = id;

   /* begin time */
   CLOCK_Start(clock);

   /* add id to result report */
   result->query_id = id;

   /* get offset into by checking index */
   q_offset = q_index->nodes[id].offset;

   /* load query */
   if ( q_filetype == FILE_FASTA ) 
   {
      SEQUENCE_Fasta_Parse( q_seq, q_filepath, q_offset );
   }
   else 
   {
      fprintf(stderr, "ERROR: Only FASTA filetypes are supported for queries.\n");
      exit(EXIT_FAILURE);
   }

   /* end and save time */
   CLOCK_Start(clock);
   times->load_query = CLOCK_Secs(clock);
}

/* load target by file index name */
void WORK_load_target_by_name( WORKER* worker,
                               char*   name )
{
   int t_id = F_INDEX_Search( worker->t_index, name );
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
   int q_id = F_INDEX_Search( worker->q_index, name );
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
   CLOCK*   clock    = worker->clock;

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
      // CLOCK_Start(clock);
      // viterbi_Lin( q_seq, t_prof, Q, T, st_MX, sp_MX, quad_sc );
      // CLOCK_Stop(clock);
      // times->lin_vit = CLOCK_Secs(clock);
   }
   if ( tasks->quad_vit ) {
      // printf("quad viterbi...\n");
      CLOCK_Start(clock);
      viterbi_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, &sc );
      CLOCK_Stop(clock);
      times->quad_vit = CLOCK_Secs(clock);
      scores->quad_vit = sc;
   }

   /* Traceback */
   if ( tasks->lin_trace ) {
      fprintf(stderr, "ERROR: Operation not supported.\n");
      exit(EXIT_FAILURE);
      // TODO: linear traceback goes here
      // CLOCK_Start(clock);
      // traceback_Build(q_seq, t_prof, Q, T, st_MX, sp_MX, tr);
      // CLOCK_Stop(clock);
      // times->lin_trace = CLOCK_Secs(clock);
   }
   if ( tasks->quad_trace ) {
      // printf("quad traceback...\n");
      CLOCK_Start(clock);
      traceback_Build( q_seq, t_prof, Q, T, st_MX, sp_MX, tr );
      CLOCK_Stop(clock);
      times->quad_trace = CLOCK_Secs(clock);
   }
}

/* forward/backward */
void WORK_forward_backward( WORKER*  worker )
{
   ARGS*          args     = worker->args;
   TASKS*         tasks    = worker->tasks;
   SCORES*        scores   = worker->scores;
   TIMES*         times    = worker->times;
   CLOCK*         clock    = worker->clock;
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
      CLOCK_Start(clock);
      backward_Linear( q_seq, t_prof, Q, T, st_MX3, st_MX, sp_MX, &sc );
      CLOCK_Stop(clock);
      times->lin_fwd    = CLOCK_Secs(clock);
      scores->lin_fwd   = sc;
   } 
   if ( tasks->quad_fwd ) {
      CLOCK_Start(clock);
      forward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, &sc );
      CLOCK_Stop(clock);
      times->quad_fwd   = CLOCK_Secs(clock);
      scores->quad_fwd  = sc;
   }

   /* backward */
   if ( tasks->lin_bck ) {
      CLOCK_Start(clock);
      backward_Linear( q_seq, t_prof, Q, T, st_MX3, st_MX, sp_MX, &sc );
      CLOCK_Stop(clock);
      times->lin_bck    = CLOCK_Secs(clock);
      scores->lin_bck   = sc;
   }
   if ( tasks->quad_bck ) {
      CLOCK_Start(clock);
      forward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, &sc );
      CLOCK_Stop(clock);
      times->lin_fwd    = CLOCK_Secs(clock);
      scores->lin_fwd   = sc;
   }
}

/* cloud search */
void WORK_cloud_search( WORKER* worker )
{
   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;
   SCORES*  scores   = worker->scores;
   TIMES*   times    = worker->times;
   CLOCK*   clock    = worker->clock;
   RESULT*  result   = worker->result;

   SEQUENCE*      q_seq    = worker->q_seq;
   HMM_PROFILE*   t_prof   = worker->t_prof;

   int   Q           = q_seq->N;
   int   T           = t_prof->N;

   float    alpha       = args->alpha;
   int      beta        = args->beta;
   bool     is_testing  = args->is_testing;
   float    sc;

   MATRIX_3D*     st_MX       = worker->st_MX;
   MATRIX_3D*     st_MX3      = worker->st_MX3;
   MATRIX_2D*     sp_MX       = worker->sp_MX;
   MATRIX_3D*     st_cloud_MX = worker->st_cloud_MX; 

   ALIGNMENT*     tr       = worker->traceback;

   EDGEBOUNDS*    edg_fwd  = worker->edg_fwd;
   EDGEBOUNDS*    edg_bck  = worker->edg_bck;
   EDGEBOUNDS*    edg_diag = worker->edg_diag;
   EDGEBOUNDS*    edg_row  = worker->edg_row;

   /* if performing linear bounded forward or backward  */
   if ( tasks->lin_bound_fwd || tasks->lin_bound_bck ) {
      /* cloud forward */
      // printf("=> cloud forward...\n");
      CLOCK_Start(clock);
      cloud_Forward_Linear( q_seq, t_prof, Q, T, st_MX3, st_MX, sp_MX, tr, edg_fwd, alpha, beta, is_testing );
      CLOCK_Stop(clock);
      times->lin_cloud_fwd = CLOCK_Secs(clock);

      /* cloud backward */
      // printf("=> cloud backward...\n");
      CLOCK_Start(clock);
      cloud_Backward_Linear( q_seq, t_prof, Q, T, st_MX3, st_MX, sp_MX, tr, edg_bck, alpha, beta, is_testing );
      CLOCK_Stop(clock);
      times->lin_cloud_bck = CLOCK_Secs(clock);

      /* merge edgebounds */
      // printf("=> merge...\n");
      CLOCK_Start(clock);
      EDGEBOUNDS_Merge( Q, T, edg_fwd, edg_bck, edg_diag );
      CLOCK_Stop(clock);
      times->lin_merge = CLOCK_Secs(clock);

      /* reorient edgebounds */
      // printf("=> reorient...\n");
      CLOCK_Start(clock);
      EDGEBOUNDS_Reorient( Q, T, edg_diag, edg_row );
      CLOCK_Stop(clock);
      times->lin_reorient = CLOCK_Secs(clock);

      /* compute the number of cells in matrix computed */
      result->cloud_cells  = EDGEBOUNDS_Count( edg_row );
      result->total_cells  = (Q+1) * (T+1);
   }
   /* bounded forward */
   if ( tasks->lin_bound_fwd ) {
      // printf("=> linear bound forward...\n");
      CLOCK_Start(clock);
      bound_Forward_Linear( q_seq, t_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row, is_testing, &sc );
      CLOCK_Stop(clock);
      times->lin_bound_fwd = CLOCK_Secs(clock);
      scores->lin_cloud_fwd = sc;
   }
   /* bounded backward */
   if ( tasks->lin_bound_bck ) {
      // printf("=> linear bound backward...\n");
      CLOCK_Start(clock);
      bound_Backward_Linear( q_seq, t_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row, is_testing, &sc );
      CLOCK_Stop(clock);
      times->lin_bound_bck = CLOCK_Secs(clock);
      scores->lin_cloud_bck = sc;
   }

   /* if performing quadratic bounded forward or backward  */
   if ( tasks->quad_bound_fwd || tasks->quad_bound_bck ) {
      /* cloud forward */
      CLOCK_Start(clock);
      cloud_Forward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, tr, edg_fwd, alpha, beta );
      CLOCK_Stop(clock);
      times->quad_cloud_fwd = CLOCK_Secs(clock);

      /* cloud backward */
      CLOCK_Start(clock);
      cloud_Backward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, tr, edg_bck, alpha, beta );
      CLOCK_Stop(clock);
      times->quad_cloud_bck = CLOCK_Secs(clock);

      /* merge edgebounds */
      CLOCK_Start(clock);
      EDGEBOUNDS_Merge( Q, T, edg_fwd, edg_bck, edg_diag );
      CLOCK_Stop(clock);
      times->quad_merge = CLOCK_Secs(clock);

      /* reorient edgebounds */
      CLOCK_Start(clock);
      EDGEBOUNDS_Reorient( Q, T, edg_diag, edg_row );
      CLOCK_Stop(clock);
      times->quad_merge = CLOCK_Secs(clock);
      
      /* compute the number of cells in matrix computed */
      result->cloud_cells = EDGEBOUNDS_Count( edg_row );
      result->total_cells  = (Q+1) * (T+1);
   }
   /* bounded forward */
   if ( tasks->quad_bound_fwd ) {
      CLOCK_Start(clock);
      bound_Forward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, edg_row, is_testing, &sc );
      CLOCK_Stop(clock);
      times->quad_bound_fwd = CLOCK_Secs(clock);
      scores->quad_cloud_fwd = sc;
   }
   /* bounded backward */
   if ( tasks->quad_bound_bck ) {
      CLOCK_Start(clock);
      bound_Backward_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, edg_row, is_testing, &sc );
      CLOCK_Stop(clock);
      times->quad_bound_bck = CLOCK_Secs(clock);
      scores->quad_cloud_bck = sc;
   }
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
   FILE* fp = worker->out_file;

   /* open file */
   fprintf(fp, "TARGET_SOURCE:\t%s\n", args->t_filepath);
   fprintf(fp, "QUERY_SOURCE:\t%s\n", args->q_filepath);
   fprintf(fp, "TARGET_INDEX:\t%s\n", args->t_indexpath);
   fprintf(fp, "QUERY_INDEX:\t%s\n", args->q_indexpath);

   /* print header */
   int pad = 0;
   fprintf(fp, ">");

   fprintf(fp, "{%*s}\t", pad, "t_id" );
   fprintf(fp, "{%*s}\t", pad, "q_id" );
   // fprintf(fp, "{%*s}\t", pad, "t_name" );
   // fprintf(fp, "{%*s}\t", pad, "q_name" );
   fprintf(fp, "{%*s}\t", pad, "t_len" );
   fprintf(fp, "{%*s}\t", pad, "q_len" );

   fprintf(fp, "{%*s}\t", pad, "alpha" );
   fprintf(fp, "{%*s}\t", pad, "beta" );

   fprintf(fp, "{%*s}\t", pad, "tot_c" );
   fprintf(fp, "{%*s}\t", pad, "cld_c" );

   fprintf(fp, "{%*s}\t", pad, "vit_sc" );
   fprintf(fp, "{%*s}\t", pad, "cl_sc" );

   fprintf(fp, "{%*s}\t", pad, "vit_t" );
   fprintf(fp, "{%*s}\t", pad, "trace_t" );
   fprintf(fp, "{%*s}\t", pad, "cl_fwd_t" );
   fprintf(fp, "{%*s}\t", pad, "cl_bck_t" );
   fprintf(fp, "{%*s}\t", pad, "merge_t" );
   fprintf(fp, "{%*s}\t", pad, "reorient_t" );
   fprintf(fp, "{%*s}\t", pad, "bnd_fwd_t" );

   fprintf(fp, "\n");
}

/* print current result (default) */
void WORK_print_result_current( WORKER* worker )
{
   FILE*    fp       = worker->out_file;
   ARGS*    args     = worker->args;
   TIMES*   times    = worker->times;
   SCORES*  scores   = worker->scores;

   fprintf(fp, "%d\t", worker->t_id );
   fprintf(fp, "%d\t", worker->q_id );
   // fprintf(fp, "%s\t", worker->t_index->nodes[i].name );
   // fprintf(fp, "%s\t", worker->q_index->nodes[i].name );
   fprintf(fp, "%d\t", worker->t_prof->N );
   fprintf(fp, "%d\t", worker->q_seq->N );

   fprintf(fp, "%9.4f\t",  args->alpha );
   fprintf(fp, "%d\t",     args->beta );

   fprintf(fp, "%d\t", worker->result->total_cells );
   fprintf(fp, "%d\t", worker->result->cloud_cells );

   fprintf(fp, "%9.4f\t", scores->quad_vit );
   fprintf(fp, "%9.4f\t", scores->lin_cloud_fwd );

   fprintf(fp, "%9.4f\t", times->quad_vit );
   fprintf(fp, "%9.4f\t", times->quad_trace );
   fprintf(fp, "%9.4f\t", times->lin_cloud_fwd );
   fprintf(fp, "%9.4f\t", times->lin_cloud_bck );
   fprintf(fp, "%9.4f\t", times->lin_merge );
   fprintf(fp, "%9.4f\t", times->lin_reorient );
   fprintf(fp, "%9.4f\t", times->lin_bound_fwd );

   fprintf(fp, "\n");
}