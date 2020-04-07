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

   /* dynamic programming matrices */
   MATRIX_3D*     st_MX       = worker->st_MX;
   MATRIX_3D*     st_MX3      = worker->st_MX3;
   MATRIX_2D*     sp_MX       = worker->sp_MX;

   /* edgebounds for cloud search */
   EDGEBOUNDS*    edg_fwd     = worker->edg_fwd;
   EDGEBOUNDS*    edg_bck     = worker->edg_bck;
   EDGEBOUNDS*    edg_diag    = worker->edg_diag;
   EDGEBOUNDS*    edg_row     = worker->edg_row;

   if ( isQuad ) {
      st_MX = MATRIX_3D_Create(NUM_NORMAL_STATES,  1, 1);
      sp_MX = MATRIX_2D_Create(NUM_SPECIAL_STATES, 1);
   }
   if ( isLin ) {
      st_MX3 = MATRIX_3D_Create(NUM_NORMAL_STATES,  1, 1);
      sp_MX  = MATRIX_2D_Create(NUM_SPECIAL_STATES, 1);
   }

   edg_fwd  = EDGEBOUNDS_Create();
   edg_bck  = EDGEBOUNDS_Create();
   edg_diag = EDGEBOUNDS_Create();
   edg_row  = EDGEBOUNDS_Create();
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

   bool           isLinVit    = tasks->linear_vit;
   bool           isQuadVit   = tasks->quad_vit;
   bool           isLinTr     = tasks->linear_trace;
   bool           isQuadTr    = tasks->quad_trace;

   float*         quad_sc  = &(scores->viterbi_quad_sc);
   float*         lin_sc   = &(scores->viterbi_sc);

   /* Viterbi */
   if ( isLinVit ) {
      viterbi_Quad( q_seq, t_prof, Q, T, st_MX, sp_MX, quad_sc );
      if ( isLinTr ) {
         traceback_Build(q_seq, t_prof, Q, T, st_MX->data, sp_MX->data, tr);
      }
   }
   if ( isQuadVit ) {
      // TODO: linear viterbi goes here
      if ( isQuadTr ) {
         // TODO: linear traceback goes here
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

   bool     isLinVit    = tasks->linear_vit;
   bool     isQuadVit   = tasks->quad_vit;
   bool     isLinTr     = tasks->linear_trace;
   bool     isQuadTr    = tasks->quad_trace;

   float*   quad_fwd_sc  = &(scores->fwd_quad_sc);
   float*   lin_fwd_sc   = &(scores->fwd_sc);
   float*   quad_bck_sc  = &(scores->bck_quad_sc);
   float*   lin_bck_sc   = &(scores->bck_sc);

   
}