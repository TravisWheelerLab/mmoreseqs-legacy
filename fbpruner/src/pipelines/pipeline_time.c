/*******************************************************************************
 *  FILE:      pipeline_time.c
 *  PURPOSE:   Time Trial Cloud Search Pipeline.
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
#include "../objects/structs.h"
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../parsers/parsers.h"
#include "../algs_linear/algs_linear.h"
#include "../algs_quad/algs_quad.h"
#include "../algs_naive/algs_naive.h"

/* header */
#include "pipelines.h"

/*  FUNCTION:  time_pipeline()
 *  SYNOPSIS:  Runs a workflow pipeline. 
 *             Takes in a single target/query pair.  
 *             Runs generic Forward-Backward algorithm
 *
 *  ARGS:      <args>     parsed commandline arguments
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
void time_pipeline( WORKER* worker ) 
{
   
}