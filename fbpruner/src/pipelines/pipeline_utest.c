/*******************************************************************************
 *  FILE:      pipeline_utest.c
 *  PURPOSE:   Unit Test Cloud Search Pipeline.
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
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../parsers/_parsers.h"
#include "../algs_linear/_algs_linear.h"
#include "../algs_quad/_algs_quad.h"
#include "../algs_naive/_algs_naive.h"
#include "../work/_work.h"

/* header */
#include "_pipelines.h"

/*! FUNCTION:  	utest_pipeline()
 *  SYNOPSIS:  	Pipeline runs unit tests.
 */
STATUS_FLAG 
utest_pipeline( WORKER* worker ) 
{
   printf("=== UNIT TEST PIPELINE ===\n");
   
   return STATUS_SUCCESS;
}