/*******************************************************************************
 *  FILE:      pipeline_mmseqs_plus.c
 *  PURPOSE:   Cloud Search Pipeline for MMSEQS: R.
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
void mmseqs_pipeline_plus( WORKER* worker )
{
   printf("# begin mmseqs...\n");
   
   /* set environmental variables */
}
