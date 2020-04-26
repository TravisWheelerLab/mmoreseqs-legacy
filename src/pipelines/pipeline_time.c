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
#include "structs.h"
#include "utilities.h"
#include "objects.h"
#include "parsers.h"
#include "algs_linear.h"
#include "algs_quad.h"
#include "algs_naive.h"

/* header */
#include "pipelines.h"

/* ****************************************************************************************** *
 *  
 *  FUNCTION:  time_pipeline()
 *  SYNOPSIS:  Runs a workflow pipeline. 
 *             Takes in a single target/query pair.  
 *             Runs generic Forward-Backward algorithm
 *
 *  ARGS:      <args>     parsed commandline arguments
 *
 *  RETURN:    No Return.
 *
/* ****************************************************************************************** */
void time_pipeline( WORKER* worker ) 
{
   
}