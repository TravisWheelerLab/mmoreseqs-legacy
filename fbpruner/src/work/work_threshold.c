/*******************************************************************************
 *  FILE:      work_fwdback.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for forward-backward.
 *
 *  AUTHOR:    Dave Rich
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
#include "../algs_sparse/_algs_sparse.h"
#include "../reporting/_reporting.h"

/* header */
#include "_work.h"
#include "work_threshold.h"

/*! FUNCTION:  	WORK_threshold_viterbi()
 *  SYNOPSIS:  	Converts Viterbi natscore to e-value, then compares to threshold.
 */
void 
WORK_threshold_viterbi( WORKER* worker )
{

}

/*! FUNCTION:  	WORK_threshold_cloud()
 *  SYNOPSIS:  	Converts Cloud natscore to e-value, then compares to threshold.
 */
void 
WORK_threshold_cloud( WORKER* worker )
{

}

/*! FUNCTION:  	WORK_threshold_bound_fwdback()
 *  SYNOPSIS:  	Converts Bound Forward natscore to e-value, then compares to threshold.
 */
void 
WORK_threshold_bound_fwdback( WORKER* worker )
{
   
}