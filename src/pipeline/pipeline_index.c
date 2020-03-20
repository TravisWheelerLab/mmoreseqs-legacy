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

/* data structures and file parsers */
#include "objects/structs.h"
#include "utility.h"
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
 *  SYNOPSIS:  Runs a workflow pipeline. 
 *             Indexes a target and query file, then saves the indexes to file.
 *
 *  ARGS:      <args>     parsed commandline arguments
 *
 *  RETURN:    No Return.
 *
/* ****************************************************************************************** */
void index_pipeline(ARGS* args) 
{
   printf("test test\n");
   F_INDEX*    t_index        = NULL; 
   F_INDEX*    q_index        = NULL;

   int         t_filepath     = args->target_filepath;
   int         q_filepath     = args->query_filepath;

   int         t_filetype     = args->target_filetype;
   int         q_filetype     = args->query_filetype; 

   char*       out_filepath   = args->output_filepath;

   /* index target file */
   if ( t_filetype == FILE_HMM ) 
   {
      t_index = F_INDEX_Hmm_Build( t_filepath );
   }
   else if ( t_filetype == FILE_FASTA )
   {
      t_index = F_INDEX_Fasta_Build( t_filepath ); 
   }

   printf("test\n");

   /* index query file */
   if ( q_filetype == FILE_HMM ) 
   {
      q_index = F_INDEX_Hmm_Build( q_filepath );
   }
   else if ( q_filetype == FILE_FASTA )
   {
      q_index = F_INDEX_Fasta_Build( q_filepath ); 
   }

   printf("test\n");

   F_INDEX_Dump( t_index, stdout );
   F_INDEX_Dump( q_index, stdout );
}
