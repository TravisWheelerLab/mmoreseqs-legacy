/*******************************************************************************
 *  FILE:      pipeline_hmmbuild.c
 *  PURPOSE:   Pipeline for building hmm's from single fasta sequence.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:      
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
   worker->output_fp = ERROR_fopen( args->output_filepath, "w+" );

}