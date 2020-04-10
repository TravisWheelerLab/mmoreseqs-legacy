/*******************************************************************************
 *  FILE:      pipeline_main.c
 *  PURPOSE:   Pipelines and Subroutines
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

/* load query and target indexes (or build them) */
void WORK_load_indexes( WORKER* worker ) 
{
   FILE*    fp       = NULL;

   ARGS*    args  = worker->args;
   TASKS*   tasks    = worker->tasks;

   char*    t_indexpath = args->t_indexpath;
   char*    q_indexpath = args->q_indexpath;

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
      /* save index file */
      args->t_indexpath = t_indexpath_tmp;
      fp = fopen( t_indexpath_tmp, "w" );
      F_INDEX_Dump( worker->t_index, fp );
      fclose(fp);
   }

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
   else if ( access( q_indexpath_tmp, F_OK ) != -1 ) {
      /* check if standard extension index file exists */
      printf("found index at database location...\n");
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
         fprintf(stderr, "ERROR: target filetype is not supported.\n" );
         exit(EXIT_FAILURE);
      }
      /* save index file */
      args->q_indexpath = q_indexpath_tmp;
      fp = fopen( q_indexpath_tmp, "w" );
      F_INDEX_Dump( worker->q_index, fp );
      fclose(fp);
   }
}

/* load target and query from file */
void WORK_load_target_query( WORKER*  worker )
{
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;

   HMM_PROFILE*   t_prof         = worker->t_prof;
   SEQUENCE*      t_seq          = worker->t_seq;
   SEQUENCE*      q_seq          = worker->q_seq;
   char*          t_filepath     = args->t_filepath;
   int            t_filetype     = args->t_filetype;
   long           t_offset       = args->t_offset;
   char*          q_filepath     = args->q_filepath;
   int            q_filetype     = args->q_filetype;
   long           q_offset       = args->q_offset;
   int            mode           = args->search_mode;

   /* load query sequence */
   printf("building q_seq sequence...\n");

   if ( q_filetype == FILE_FASTA ) 
   {
      SEQUENCE_Fasta_Parse( q_seq, q_filepath, q_offset );
   }
   else 
   {
      fprintf(stderr, "ERROR: Only FASTA filetypes are supported for queries.\n");
      exit(EXIT_FAILURE);
   }
   
   /* load target profile */
   printf("building hmm profile...\n");
   if ( t_filetype == FILE_HMM ) 
   {
      HMM_PROFILE_Parse( t_prof, t_filepath, t_offset ); 
      HMM_PROFILE_Convert_NegLog_To_Real( t_prof );
      HMM_PROFILE_Config( t_prof, mode );
   }
   else if ( t_filetype == FILE_FASTA )
   {
      SEQUENCE_Fasta_Parse( t_seq, t_filepath, t_offset );
      SEQUENCE_to_HMM_PROFILE( t_seq, t_prof );
   }
   else
   {
      fprintf(stderr, "ERROR: Only HMM and FASTA filetypes are supported for t_profs.\n");
      exit(EXIT_FAILURE);
   }
   HMM_PROFILE_ReconfigLength( t_prof, q_seq->N );
}

/* initialize dynamic programming matrices, edgebounds,  */
void WORK_init( WORKER* worker )
{
   int   T        = worker->t_prof->N;
   int   Q        = worker->q_seq->N;

   /* are we performing linear or quadratic algorithms? */
   int   isQuad   = worker->tasks->quadratic;
   int   isLin    = worker->tasks->linear;

   /* create necessary dp matrices */
   if ( isQuad ) {
      worker->st_MX = MATRIX_3D_Create(NUM_NORMAL_STATES,  1, 1);
      worker->sp_MX = MATRIX_2D_Create(NUM_SPECIAL_STATES, 1);
   }
   if ( isLin ) {
      worker->st_MX3 = MATRIX_3D_Create(NUM_NORMAL_STATES,  1, 1);
      worker->sp_MX  = MATRIX_2D_Create(NUM_SPECIAL_STATES, 1);
   }
   
}

/* initialize dynamic programming matrices */
void WORK_reuse( WORKER* worker )
{
   int   T  = worker->t_prof->N;
   int   Q  = worker->q_seq->N;

   /* are we performing linear or quadratic algorithms? */
   int   isQuad   = worker->tasks->quadratic;
   int   isLin    = worker->tasks->linear;

   /* dynamic programming matrices */
   MATRIX_3D*  st_MX    = worker->st_MX;
   MATRIX_3D*  st_MX3   = worker->st_MX3;
   MATRIX_2D*  sp_MX    = worker->sp_MX;

   /* edgebounds for cloud search */
   EDGEBOUNDS*    edg_fwd     = worker->edg_fwd;
   EDGEBOUNDS*    edg_bck     = worker->edg_bck;
   EDGEBOUNDS*    edg_diag    = worker->edg_diag;
   EDGEBOUNDS*    edg_row     = worker->edg_row;

   if ( isQuad ) {
      MATRIX_3D_Reuse( st_MX, NUM_NORMAL_STATES,  Q+1, T+1 );
      MATRIX_2D_Reuse( sp_MX, NUM_SPECIAL_STATES, Q+1 );
   }
   if ( isLin ) {
      MATRIX_3D_Reuse( st_MX3, NUM_NORMAL_STATES,  3, (Q+T+1) );
      MATRIX_2D_Reuse( sp_MX, NUM_SPECIAL_STATES, Q+1 );
   }

   EDGEBOUNDS_Reuse( edg_fwd );
   EDGEBOUNDS_Reuse( edg_bck );
   EDGEBOUNDS_Reuse( edg_diag );
   EDGEBOUNDS_Reuse( edg_row );
}

/* initial output printed before search */
void WORK_print_header( WORKER*  worker ) 
{

}

/* viterbi and traceback */
void WORK_viterbi_and_traceback( WORKER*  worker )
{
   ARGS*          args     = worker->args;
   TASKS*         tasks    = worker->tasks;
   SCORES*        scores   = worker->scores;
   TIMES*         times    = worker->times;

   SEQUENCE*      q_seq    = worker->q_seq;
   HMM_PROFILE*   t_prof   = worker->t_prof;

   int            Q        = q_seq->N;
   int            T        = t_prof->N;

   MATRIX_3D*     st_MX    = worker->st_MX;
   MATRIX_3D*     st_MX3   = worker->st_MX3;
   MATRIX_2D*     sp_MX    = worker->sp_MX;
   ALIGNMENT*     tr       = worker->traceback;

   bool           isLinVit    = tasks->lin_vit;
   bool           isQuadVit   = tasks->quad_vit;
   bool           isLinTr     = tasks->lin_trace;
   bool           isQuadTr    = tasks->quad_trace;

   float*         quad_sc  = &(scores->viterbi_quad_sc);
   float*         lin_sc   = &(scores->viterbi_sc);

   /* Viterbi */
   if ( isQuadVit ) {
      viterbi_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, quad_sc );
      if ( isQuadTr ) {
         traceback_Build(q_seq, t_prof, Q, T, st_MX->data, sp_MX->data, tr);
      }
   }
   if ( isLinVit ) {
      // TODO: linear viterbi goes here
      // viterbi_Lin( q_seq, t_prof, Q, T, st_MX, sp_MX, quad_sc );
      if ( isLinTr ) {
         // TODO: linear traceback goes here
         // traceback_Lin_Build(q_seq, t_prof, Q, T, st_MX->data, sp_MX->data, tr);
      }
   }
}

/* forward/backward */
void WORK_forward_backward( WORKER*  worker )
{
   ARGS*          args     = worker->args;
   TASKS*         tasks    = worker->tasks;
   SCORES*        scores   = worker->scores;
   TIMES*         times    = worker->times;

   SEQUENCE*      q_seq    = worker->q_seq;
   HMM_PROFILE*   t_prof   = worker->t_prof;

   int   Q     = q_seq->N;
   int   T     = t_prof->N;

   MATRIX_3D*     st_MX    = worker->st_MX;
   MATRIX_3D*     st_MX3   = worker->st_MX3;
   MATRIX_2D*     sp_MX    = worker->sp_MX;
   ALIGNMENT*     tr       = worker->traceback;

   bool     isLinFwd    = tasks->lin_fwd;
   bool     isLinBck    = tasks->lin_bck;
   bool     isQuadFwd   = tasks->quad_fwd;
   bool     isQuadBck   = tasks->quad_bck;

   float*   quad_fwd_sc  = &(scores->fwd_quad_sc);
   float*   lin_fwd_sc   = &(scores->fwd_sc);
   float*   quad_bck_sc  = &(scores->bck_quad_sc);
   float*   lin_bck_sc   = &(scores->bck_sc);

   if ( isLinFwd ) {

   } 
}