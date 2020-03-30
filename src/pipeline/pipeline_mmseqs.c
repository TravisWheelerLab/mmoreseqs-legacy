/*******************************************************************************
 *  FILE:      pipeline_mmseqs.c
 *  PURPOSE:   Cloud Search Pipeline for MMSEQS.
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
#include "../parsers/hmm_parser.h"
#include "../parsers/seq_parser.h"
#include "../parsers/m8_parser.h"
#include "../parsers/index_parser.h"

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
void mmseqs_pipeline(ARGS* args) 
{
   /* Get Arguments */
   float          alpha             = args->alpha;
   int            beta              = args->beta;

   char*          t_filepath        = args->target_filepath;
   char*          q_filepath        = args->query_filepath;

   char*          t_indexpath       = args->target_indexpath;
   char*          q_indexpath       = args->query_indexpath;

   int            t_filetype        = args->target_filetype;
   int            q_filetype        = args->query_filetype;

   char*          m8_filepath       = args->hits_filepath;

   /* results inputted from mmseqs pipeline */
   RESULTS*       results_in        = NULL;
   RESULT*        result            = NULL;
   RESULT*        result_prv        = NULL;

   /* results outputted after cloud fwd/bck */
   RESULTS*       results_out       = RESULTS_Create();
   
   /* target & query objects */
   SEQUENCE*      query             = NULL;
   HMM_PROFILE*   target            = NULL;

   /* file indexes for target & query */
   F_INDEX*       q_index           = NULL;
   F_INDEX*       t_index           = NULL;

   /* dynamic programming matrices for computing cloud fwd/bck */
   MATRIX_3D*     st_matrix         = MATRIX_3D_Create(1, 3, NUM_NORMAL_STATES);
   float*         st_MX3            = NULL;
   MATRIX_2D*     sp_matrix         = MATRIX_2D_Create(1, NUM_SPECIAL_STATES);
   float*         sp_MX             = NULL;

   ALIGNMENT*     aln               = ALIGNMENT_Create();

   EDGEBOUNDS*    edg_fwd           = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_bck           = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_diag          = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_row           = EDGEBOUNDS_Create();

   /* clock for timing cloud forward-backward */
   CLOCK*         cl                = CLOCK_Create();
   /* file pointer for writing out to file */
   FILE*          fp                = NULL;

   /* problem size */
   int   T  = 1;
   int   Q  = 1;

   /* threshold scores */
   float threshold_sc   = 0;
   float sc             = 0;

   /* input results file from MMSEQS pipeline */
   results_in = RESULTS_M8_Parse(m8_filepath);
   RESULTS_M8_Dump(results_in, stdout);
   
   exit(0);

   /* initialize logrithmic sum table */
   init_Logsum();

   /* === INDEX FILES === */

   /* if target file doesn't have an index file, create one */
   if ( t_indexpath == NULL ) 
   {
      switch( args->target_filetype )
      {
         case FILE_HMM:
            t_index = F_INDEX_Hmm_Build( args->target_filepath );
            break;

         case FILE_FASTA:
            t_index = F_INDEX_Fasta_Build( args->target_filepath );
            break;

         default:
            fprintf(stderr, "ERROR: File type of '%s' is not supported.\n", args->target_filepath);
            exit(EXIT_FAILURE);
      }
   } 
   else 
   {
      t_index = F_INDEX_Load( args->target_indexpath );
   }
   F_INDEX_Sort( t_index );

   /* if query file doesn't have an index file, create one */
   if ( args->query_indexpath == NULL ) 
   {
      switch( args->query_filetype )
      {
         case FILE_FASTA:
            q_index = F_INDEX_Fasta_Build( args->target_filepath );
            break;

         default:
            fprintf(stderr, "ERROR: File type of '%s' is not supported.\n", args->target_filepath);
            exit(EXIT_FAILURE);
      }
   } 
   else 
   {
      q_index = F_INDEX_Load( args->target_indexpath );
   }
   F_INDEX_Sort( q_index );

   /* === ITERATE OVER EACH RESULT === */

   /* Look through each input result */
   for (int i = 0; i < results_in->N; i++) 
   {
      /* get next result from list */
      result = &(results_in->data[i]);

      /* if current result is not the same query as previous result, load it in. */
      if ( result_prv != NULL && strcmp( result->target_name, result_prv->target_name ) ) 
      {
         /* find query in index file */
         long offset = F_INDEX_Search( t_index, result->target_name );

         /* jump to position in query file and load in query */
         target = HMM_PROFILE_Parse( args->target_filepath, offset );
      }

      /* if current result is not the same query as previous result, load it in. */
      if ( result_prv != NULL && strcmp( result->query_name, result_prv->query_name ) ) 
      {
         /* find query in index file */
         long offset = F_INDEX_Search( q_index, result->query_name );

         /* jump to position in query file and load in query */
         query = SEQUENCE_Fasta_Parse( args->query_filepath, offset );
      }

      /* resize data structures for next search */
      T = target->N;
      Q = query->N;
      MATRIX_3D_Reuse( st_matrix, NUM_NORMAL_STATES, 3, (Q+T+1) );
      st_MX3 = st_matrix->data;
      MATRIX_2D_Reuse( sp_matrix, NUM_SPECIAL_STATES, Q+1 );
      sp_MX = sp_matrix->data;

      /* get initial search window */
      ALIGNMENT_Clear(aln);
      ALIGNMENT_Pushback(aln, &((TRACE){ result->query_start, result->target_start, M_ST }) );
      ALIGNMENT_Pushback(aln, &((TRACE){ result->query_end, result->target_end, M_ST }) );
      aln->beg = 0;
      aln->end = 1;

      /* empty edgebound data */
      EDGEBOUNDS_Clear(edg_fwd);
      EDGEBOUNDS_Clear(edg_bck);

      /* === CLOUD SEARCH === */

      /* perform cloud search */
      cloud_Forward_Linear(query, target, Q, T, NULL, st_MX3, sp_MX, aln, edg_fwd, alpha, beta, false);
      cloud_Backward_Linear(query, target, Q, T, NULL, st_MX3, sp_MX, aln, edg_bck, alpha, beta, false);

      edg_diag = EDGEBOUNDS_Merge(Q, T, edg_fwd, edg_bck);
      edg_row  = EDGEBOUNDS_Reorient(Q, T, edg_diag);

      bound_Forward_Linear(query, target, Q, T, st_MX3, NULL, sp_MX, edg_row, false, &sc);
      result->cloud_fwd_sc = sc;

      /* if it clears scoring threshold, add to results */
      if ( sc > threshold_sc ) {
         RESULTS_PushBack(results_out, result);
      }

      free(edg_diag);
      free(edg_row);
   }

   /* final output of results */
   // fp = fopen()
   fp = stdout;
   RESULTS_My_Dump( results_out, fp );
   if (fp != stdout) fclose(fp);
}