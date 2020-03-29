/*******************************************************************************
 *  FILE:      pipeline_main.c
 *  PURPOSE:   Main Cloud Search Pipeline.
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

/* data structures and file parsers */
#include "objects/structs.h"
#include "utility.h"
#include "hmm_parser.h"
#include "seq_parser.h"

/* objects */
#include "objects/f_index.h"
#include "objects/alignment.h"
#include "objects/sequence.h"
#include "objects/hmm_profile.h"
#include "objects/edgebound.h"
#include "objects/clock.h"
#include "objects/matrix/matrix_2d.h"
#include "objects/matrix/matrix_3d.h"
#include "objects/vectors/vector_range_2d.h"

/* viterbi & fwdbck (quadratic) */
#include "viterbi_quad.h"
#include "traceback_quad.h"
#include "fwdback_quad.h"

/* cloud search (naive) */
#include "bounded_fwdbck_naive.h"
/* cloud search (quadratic space) */
#include "cloud_search_quad.h"
#include "merge_reorient_quad.h"
#include "bounded_fwdbck_quad.h"
/* cloud search (linear space) */
#include "cloud_search_linear.h"
#include "merge_reorient_linear.h"
#include "bounded_fwdbck_linear.h"
/* temp test */
#include "cloud_search_linear_rows.h"

/* debugging methods */
#include "testing.h"

/* header */
#include "pipeline.h"

/* standard pipeline */
void main_pipeline(ARGS* args) 
{
   /* Get Arguments */
   float    alpha          = args->alpha;
   int      beta           = args->beta;

   char*    t_filepath     = args->target_filepath;
   char*    q_filepath     = args->query_filepath;

   // int      mode        = MODE_MULTILOCAL;    /* HMMER standard mode (allows jumps) */
   int      mode           = MODE_UNILOCAL;      /* Cloud Search mode (prohibiits jumps) */

   /* Target & Query */
   HMM_PROFILE* target     = NULL;
   SEQUENCE*    query      = NULL;

   /* Timing & Scoring */
   SCORES*  scores      = (SCORES*) malloc( sizeof(SCORES) );
   TIMES*   times       = (TIMES*) malloc( sizeof(TIMES) );
   CLOCK*   cl          = CLOCK_Create();

   bool     is_testing  = false;
   bool     test        = false;
   time_t   runtime     = 0; 
   time_t   runtime_tot = 0;

   /* records percentage of space computed */
   int      T           = 0;
   int      Q           = 0;
   float    sc          = 0.f;

   // printf("=== LOAD HMM_PROFILE / QUERY -> START ===\n");
   target = HMM_PROFILE_Parse(t_filepath, 0);
   HMM_PROFILE_Config(target, mode);
   T = target->N;

   query = SEQUENCE_Fasta_Parse(q_filepath, 0);
   Q = query->N;

   ALIGNMENT*  tr             = ALIGNMENT_Create();
   EDGEBOUNDS* edg_fwd_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS* edg_bck_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS* edg_row_lin    = NULL;
   EDGEBOUNDS* edg_diag_lin   = NULL;

   /* allocate memory for quadratic algs (for DEBUGGING) */
   MATRIX_3D*  st_MATRIX      = MATRIX_3D_Create(NUM_NORMAL_STATES,  Q+1, T+1);
   float*      st_MX          = st_MATRIX->data;
   MATRIX_2D*  sp_MATRIX      = MATRIX_2D_Create(NUM_SPECIAL_STATES, Q+1);
   float*      sp_MX          = sp_MATRIX->data;
   /* allocate memory for linear algs */
   MATRIX_3D*  st_MATRIX3     = MATRIX_3D_Create(NUM_NORMAL_STATES, 3, (Q+T+1) );
   float*      st_MX3         = st_MATRIX3->data;

   /* === INDEX FILES === */
   F_INDEX*    t_index        = F_INDEX_Create( t_filepath );
   F_INDEX*    q_index        = F_INDEX_Create( q_filepath );


   /* ==== CLOUD SEARCH === */

   /* run viterbi algorithm */
   viterbi_Quad(query, target, Q, T, st_MX, sp_MX, &sc);

   /* build traceback */
   traceback_Build(query, target, Q, T, st_MX, sp_MX, tr);

   /* run forward/backward algorithms */
   init_Logsum();

   cloud_Forward_Linear(query, target, Q, T, st_MX, st_MX3, sp_MX, tr, edg_fwd_lin, alpha, beta, is_testing);
   cloud_Backward_Linear(query, target, Q, T, st_MX, st_MX3, sp_MX, tr, edg_bck_lin, alpha, beta, is_testing);

   edg_diag_lin = EDGEBOUNDS_Merge(Q, T, edg_fwd_lin, edg_bck_lin);
   edg_row_lin  = EDGEBOUNDS_Reorient(Q, T, edg_diag_lin);

   bound_Forward_Linear(query, target, Q, T, st_MX3, st_MX, sp_MX, edg_row_lin, is_testing, &sc);
   scores->cloud_fwd_sc = sc;
   bound_Backward_Linear(query, target, Q, T, st_MX3, st_MX, sp_MX, edg_row_lin, is_testing, &sc);
   scores->cloud_bck_sc = sc;

   printf("=== SCORES ===\n");
   printf("BOUND_FWD: %f\n", scores->cloud_fwd_sc );
   printf("BOUND_BCK: %f\n", scores->cloud_bck_sc );
   printf("==============\n\n");

   /* free matrices */
   MATRIX_3D_Destroy(st_MATRIX);
   MATRIX_2D_Destroy(sp_MATRIX);
   MATRIX_3D_Destroy(st_MATRIX3);

   /* free hmm_profile and sequence*/
   HMM_PROFILE_Destroy(target);
   SEQUENCE_Destroy(query);

   /* free alignment */
   ALIGNMENT_Destroy(tr);

   /* free edgebounds */
   EDGEBOUNDS_Destroy(edg_fwd_lin);
   EDGEBOUNDS_Destroy(edg_bck_lin);
   EDGEBOUNDS_Destroy(edg_diag_lin);
   EDGEBOUNDS_Destroy(edg_row_lin);
}