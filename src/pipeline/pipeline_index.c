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

   char*       t_indexpath    = args->t_indexpath;
   char*       q_indexpath    = args->q_indexpath;

   int         t_filetype     = args->t_filetype;
   int         q_filetype     = args->q_filetype; 

   char*       t_lookup_filepath = args->t_lookup_filepath;
   char*       q_lookup_filepath = args->q_lookup_filepath;   

   char*       out_filepath   = args->output_filepath;

   /* index target file */
   if ( t_filetype == FILE_HMM ) {
      t_index = F_INDEX_Hmm_Build( t_filepath );
   }
   else if ( t_filetype == FILE_FASTA ) {
      t_index = F_INDEX_Fasta_Build( t_filepath ); 
   }
   t_index->index_path = args->t_indexpath;

   /* index query file */
   if ( q_filetype == FILE_HMM ) {
      q_index = F_INDEX_Hmm_Build( q_filepath );
   }
   else if ( q_filetype == FILE_FASTA ) {
      q_index = F_INDEX_Fasta_Build( q_filepath ); 
   }
   q_index->index_path = args->q_indexpath;

   /* if we have a mmseqs lookup, update names to match */
   if ( t_lookup_filepath != NULL ) {
      F_INDEX_Lookup_Update( t_index, t_lookup_filepath );
   }
   if ( q_lookup_filepath != NULL ) {
      F_INDEX_Lookup_Update( q_index, q_lookup_filepath );
   }

   /* determine the method of output for target index */
   if ( t_indexpath == NULL ) {
      /* if no output name is given, append ".idx" and save in same directory */
      const char* ext = ".idx";
      t_indexpath = (char*) malloc( sizeof(char) * (strlen(t_filepath) + strlen(ext) + 1) );
      if (t_indexpath == NULL) {
         fprintf(stderr, "ERROR: malloc failed.\n");
      }
      strcpy( t_indexpath, t_filepath );
      strcat( t_indexpath, ext );
   }
   /* output target index */
   fp = fopen( t_indexpath, "w" );
   F_INDEX_Dump( t_index, fp );
   // F_INDEX_Sort_by_Name( t_index );
   fclose(fp);
   
   /* determine the method of output for query index */
   if ( q_indexpath == NULL ) {
      /* if no output name is given, append ".idx" and save in same directory */
      const char* ext = ".idx";
      q_indexpath = (char*) malloc( sizeof(char) * (strlen(q_filepath) + strlen(ext) + 1) );
      if (t_indexpath == NULL) {
         fprintf(stderr, "ERROR: malloc failed.\n");
      }
      strcpy( q_indexpath, q_filepath );
      strcat( q_indexpath, ".idx" );
   }
   /* output query index */
   fp = fopen( q_indexpath, "w" );
   F_INDEX_Dump( q_index, fp );
   // F_INDEX_Sort_by_Name( q_index );
   fclose(fp);
}
