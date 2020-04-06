/*******************************************************************************
 *  FILE:      pipeline_index.c
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

/* data structures */
#include "objects/structs.h"
#include "utilities/utility.h"

/* file parsers */
#include "arg_parser.h"
#include "seq_parser.h"
#include "hmm_parser.h"
#include "index_parser.h"

/* objects */
#include "objects/sequence.h"
#include "objects/f_index.h"
#include "objects/hmm_profile.h"
#include "objects/alignment.h"
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

/* ****************************************************************************************** *
 *  
 *  FUNCTION:  index_pipeline()
 *  SYNOPSIS:  Indexing workflow pipeline. 
 *             Indexes a target and query file, then saves those indexes to file.
 *
 *  ARGS:      <args>     parsed commandline arguments
 *
 *  RETURN:    No Return.
 *
/* ****************************************************************************************** */
void index_pipeline( WORKER* worker ) 
{
   ARGS*       args           = worker->args;
   FILE*       fp             = NULL;

   F_INDEX*    t_index        = NULL; 
   F_INDEX*    q_index        = NULL;

   char*       t_filepath     = args->t_filepath;
   char*       q_filepath     = args->q_filepath;

   int         t_filetype     = args->t_filetype;
   int         q_filetype     = args->q_filetype; 

   char*       out_filepath   = args->output_filepath;

   /* index target file */
   if ( t_filetype == FILE_HMM ) {
      t_index = F_INDEX_Hmm_Build( t_filepath );
   }
   else if ( t_filetype == FILE_FASTA ) {
      t_index = F_INDEX_Fasta_Build( t_filepath ); 
   }

   /* index query file */
   if ( q_filetype == FILE_HMM ) {
      q_index = F_INDEX_Hmm_Build( q_filepath );
   }
   else if ( q_filetype == FILE_FASTA ) {
      q_index = F_INDEX_Fasta_Build( q_filepath ); 
   }

   if ( strcmp( out_filepath, "#stdout" ) == 0 ) {
      F_INDEX_Dump( t_index, stdout );
      F_INDEX_Sort( t_index );

      F_INDEX_Dump( q_index, stdout );
      F_INDEX_Sort( q_index );
   }
   else
   {
      char* t_indexpath = strdup( t_filepath );
      strcat( t_indexpath, ".idx" );
      printf("saving target_index to:\t%s...\n", t_indexpath);

      fp = fopen(t_indexpath, "w");
      F_INDEX_Dump( t_index, fp );
      F_INDEX_Sort( t_index );
      fclose(fp);

      char* q_indexpath = strdup( q_filepath );
      strcat( q_indexpath, ".idx" );
      printf("saving query_index to:\t%s...\n", q_indexpath);

      fp = fopen(q_indexpath, "w");
      F_INDEX_Dump( q_index, fp );
      F_INDEX_Sort( q_index );
      fclose(fp);
   }
   
}
