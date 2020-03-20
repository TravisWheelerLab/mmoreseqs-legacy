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
   char*          target_filepath   = args->target_filepath;
   char*          query_filepath    = args->query_filepath;
   char*          m8_filepath       = args->hits_filepath;

   /* results inputted from mmseqs pipeline */
   RESULTS*       results_in        = NULL;
   RESULT*        result            = NULL;
   RESULT*        result_prv        = NULL;

   /* results outputted after cloud fwd/bck */
   RESULTS*       results_out       = NULL;
   
   /* target & query objects */
   SEQUENCE*      query_seq         = NULL;
   HMM_PROFILE*   target_prof       = NULL;

   /* file indexes for target & query */
   F_INDEX*       query_idx         = NULL;
   F_INDEX*       target_idx        = NULL;

   /* dynamic programming matrices for computing cloud fwd/bck */
   MATRIX_3D*     st_matrix         = MATRIX_3D_Create(1, 3, NUM_NORMAL_STATES);
   float*         st_MX             = NULL;
   MATRIX_2D*     sp_matrix         = MATRIX_2D_Create(1, NUM_SPECIAL_STATES);
   float*         sp_MX             = NULL;

   /* input results file from MMSEQS pipeline */
   results_in = RESULTS_M8_Parse(m8_filepath);
   RESULTS_M8_Dump(results_in, stdout);

   /* determine file type of target */
   for ( int i = 0; i < NUM_FILE_TYPES; i++ ) {
      if ( STRING_EndsWith( args->target_filepath, FILE_TYPE_NAMES[i], strlen(FILE_TYPE_NAMES[i]) ) == 0 ) {
         args->target_filetype = FILE_TYPE_MAP[i];
      }
   }

   /* if target file doesn't have an index file, create one */
   if ( args->target_indexpath == NULL ) 
   {
      switch( args->target_filetype )
      {
         case FILE_HMM:
            target_idx = F_INDEX_Hmm_Build( args->target_filepath );
            break;

         case FILE_FASTA:
            target_idx = F_INDEX_Fasta_Build( args->target_filepath );
            break;

         default:
            fprintf(stderr, "ERROR: File type of '%s' is not supported.\n", args->target_filepath);
            exit(EXIT_FAILURE);
      }
   } 
   else 
   {
      target_idx = F_INDEX_Load( args->target_indexpath );
   }

   /* determine file type of query */
   for ( int i = 0; i < NUM_FILE_TYPES; i++ ) {
      if ( STRING_EndsWith( args->query_filepath, FILE_TYPE_NAMES[i], strlen(FILE_TYPE_NAMES[i]) ) == 0 ) {
         args->query_filetype = FILE_TYPE_MAP[i];
      }
   }

   /* if query file doesn't have an index file, create one */
   if ( args->query_indexpath == NULL ) 
   {
      switch( args->query_filetype )
      {
         case FILE_FASTA:
            target_idx = F_INDEX_Fasta_Build( args->target_filepath );
            break;

         default:
            fprintf(stderr, "ERROR: File type of '%s' is not supported.\n", args->target_filepath);
            exit(EXIT_FAILURE);
      }
   } 
   else 
   {
      target_idx = F_INDEX_Load( args->target_indexpath );
   }

   /* Look through each input result */
   for (int i = 0; i < results_in->N; i++) 
   {
      /* get next result from list */
      result = &(results_in->data[i]);

      /* if current result is not the same query as previous result, load it in. */
      if ( result_prv != NULL && strcmp( result->target_name, result_prv->target_name ) ) 
      {
         /* find query in index file */
         long offset = F_INDEX_Search( target_idx, result->target_name );

         /* jump to position in query file and load in query */
         target_prof = HMM_PROFILE_Parse( args->target_filepath, offset );
      }

      /* if current result is not the same query as previous result, load it in. */
      if ( result_prv != NULL && strcmp( result->query_name, result_prv->query_name ) ) 
      {
         /* find query in index file */
         long offset = F_INDEX_Search( query_idx, result->query_name );

         /* jump to position in query file and load in query */
         query_seq = SEQUENCE_Fasta_Parse( args->query_filepath, offset );
      }
     
      
   }
}