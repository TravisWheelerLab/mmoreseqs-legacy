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
   fprintf( fp, "# fb-pruner :: heuristic pruning of Forward-Backward for faster profile/sequence search\n" );
   fprintf( fp, "# FB-PRUNER 1.0 (August 2020): http://github.com/TravisWheelerLab/fb-pruner/\n" );
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

}

/* print alignment */
void REPORTER_alignment(   WORKER*  worker,
                           FILE*    fp )
{
   
}

/* print horizontal rule */
void REPORTER_hr( FILE* fp )
{
   fprintf( fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n" );
}