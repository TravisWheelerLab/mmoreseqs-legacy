/*******************************************************************************
 *  FILE:      main.c
 *  PURPOSE:   Main Method, Parses Command Line Arguments, then
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

/* data structures and utility functions */
#include "objects/structs.h"
#include "utility.h"

/* file parsers */
#include "arg_parser.h"
#include "hmm_parser.h"
#include "seq_parser.h"
#include "index_parser.h"

/* objects */
#include "objects/alignment.h"
#include "objects/sequence.h"
#include "objects/hmm_profile.h"
#include "objects/edgebound.h"
#include "objects/clock.h"
#include "objects/f_index.h"
#include "objects/mystring.h"
#include "objects/score_matrix.h"

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

/* pipelines */
#include "pipeline.h"

/* set debug macros */
#ifndef DEBUG
   #define DEBUG false
#endif

/* === HEADER === */

/* === MAIN ENTRY-POINT TO PROGRAM === */
int main ( int argc, char *argv[] )
{
   /* parse command line arguments */
   ARGS* args  = NULL;
   args = ARGS_Parse(argc, argv);

   /* output arguments */
   ARGS_Dump( args, stdout );

   /* jumps to pipeline based on -p flag */
   printf("> Running %s...\n\n", PIPELINE_NAMES[args->pipeline_mode] );
   PIPELINES[ args->pipeline_mode ]( args );

   exit(EXIT_SUCCESS);
}
