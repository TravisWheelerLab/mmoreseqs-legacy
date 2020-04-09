/*******************************************************************************
 *  FILE:      application.c
 *  PURPOSE:   Entry Point to Application, Argument Parsing
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

/* data structures, general utilities, and testing */
#include "objects/structs.h"
#include "utilities/utility.h"
#include "utilities/error_handler.h"
#include "utilities/testing.h"

/* file parsers */
#include "parsers/arg_parser.h"
#include "parsers/hmm_parser.h"
#include "parsers/seq_parser.h"
#include "parsers/m8_parser.h"
#include "parsers/index_parser.h"
#include "parsers/seq_to_profile.h"

/* objects */
#include "objects/worker.h"
#include "objects/f_index.h"
#include "objects/results.h"
#include "objects/alignment.h"
#include "objects/sequence.h"
#include "objects/hmm_profile.h"
#include "objects/edgebound.h"
#include "objects/clock.h"
#include "objects/mystring.h"
#include "objects/matrix/matrix_2d.h"
#include "objects/matrix/matrix_3d.h"
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

/* set debug macros */
#ifndef DEBUG
#define DEBUG false
#endif

/* === HEADER === */

/* === MAIN ENTRY-POINT TO PROGRAM === */
int main ( int argc, char *argv[] )
{
   #if DEBUG
      printf("In DEBUG MODE...\n");
   #endif

   /* parse command line arguments */
   ARGS* args  = NULL;
   args = ARGS_Parse( argc, argv );

   /* output arguments */
   ARGS_Dump( args, stdout );

   /* build worker */
   WORKER* worker = WORKER_Create();
   worker->args = args;
   worker->tasks = (TASKS*) malloc( sizeof(TASKS) );

   /* jumps to pipeline based on -p flag */
   printf("> Running %s...\n\n", PIPELINE_NAMES[args->pipeline_mode] );
   PIPELINES[ args->pipeline_mode ]( worker );

   exit(EXIT_SUCCESS);
}