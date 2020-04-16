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
   FILE*       fp             = NULL;

   ARGS*       args           = worker->args;
   TASKS*      tasks          = worker->tasks;
   CLOCK*      clock          = worker->clock;

   printf("load target index...\n");
   WORK_load_target_index( worker );
   printf("load query index...\n");
   WORK_load_query_index( worker );

   printf("output target index...\n");
   WORK_output_target_index( worker );
   printf("output query index...\n");
   WORK_output_query_index( worker );

   printf("\nTARGET:\n");
   F_INDEX_Dump( worker->t_index, stdout );
   printf("QUERY:\n");
   F_INDEX_Dump( worker->q_index, stdout );
}
