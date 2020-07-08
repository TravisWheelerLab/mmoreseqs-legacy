/*******************************************************************************
 *  FILE:      pipeline_hmmbuild.c
 *  PURPOSE:   Pipeline for building 
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
#include "pipelines.h"

/* mmseqs pipeline */
void hmmbuild_pipeline( WORKER* worker )
{
   /* worker objects */
   ARGS*       args        = worker->args;
   TASKS*      tasks       = worker->tasks;
   TIMES*      times       = worker->times;

   /* set file pointer */
   worker->out_file = fopen( args->output_filepath, "w+" );
   if ( worker->out_file == NULL ) {
      fprintf(stderr, "ERROR: Bad File Pointer.\n");
      exit(EXIT_FAILURE);
   }

}