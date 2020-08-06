/*******************************************************************************
 *  FILE:      report.c
 *  PURPOSE:   Reporting Subroutines
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

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "objects.h"
#include "parsers.h"
#include "algs_linear.h"
#include "algs_quad.h"
#include "algs_naive.h"

/* header */
#include "report.h"

/* print header to console */
void REPORTER_header(   WORKER*  worker,
                        FILE*    fp )
{
   /* TODO: Name is work in progress */
   fprintf( fp, "# fb-pruner :: %s :: heuristic pruning of Forward-Backward for faster profile/sequence search\n", PIPELINE_NAMES[args->pipeline_mode] );
   fprintf( fp, "# FB-PRUNER 0.1 (August 2020): http://github.com/TravisWheelerLab/fb-pruner/\n" );
   fprintf( fp, "# Copyright (C) 2020 Travis Wheeler Lab, University of Montana.\n" );
   REPORTER_hr( fp );
   fprintf( fp, "#%25s %s", "HMM file:", worker->args->t_filepath );
   fprintf( fp, "#%25s %s", "sequence database:", worker->args->q_filepath );
   REPORTER_hr( fp );
   printf("\n\n");
}

/* print summary statistics */
void REPORTER_summary_stats(  WORKER*  worker,
                              FILE*    fp )
{
   /* TODO */
}

/* print alignment */
void REPORTER_alignment(   WORKER*  worker,
                           FILE*    fp )
{
   HMM_PROFILE_Set_Consensus( worker->t_prof );

   ALIGNMENT*     aln      = worker->traceback;
   char*          t_in     = worker->t_prof->consensus;
   char*          q_in     = worker->q_seq->seq;

   VECTOR_CHAR*   t_out    = VECTOR_CHAR_Create();
   VECTOR_CHAR*   q_out    = VECTOR_CHAR_Create();
   TRACE*         tr;

   /* report target */
   for (int i = aln->beg; i < aln->end; i++)
   {
      if ( tr->st == M_ST ) {

      }
   }
}

/* print horizontal rule */
void REPORTER_hr( FILE* fp )
{
   fprintf( fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n" );
}