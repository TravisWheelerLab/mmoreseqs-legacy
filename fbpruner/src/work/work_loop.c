/*******************************************************************************
 *  FILE:      work_loop.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for inner search loop.
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

/*! FUNCTION:  	WORK_preloop()
 *  SYNOPSIS:  	Prep <worker> for main loop.
 */
void 
WORK_preloop( WORKER* worker )
{

}

/*! FUNCTION:  	WORK_preiter()
 *  SYNOPSIS:  	Prep <worker> for next main loop iteration.
 */
void
WORK_preiter( WORKER* worker )
{

}

/*! FUNCTION:  	WORK_postiter()
 *  SYNOPSIS:  	Clean up <worker> at end of main loop iteration.
 */
void 
WORK_postiter( WORKER* worker )
{

}

/*! FUNCTION:  	WORK_postloop()
 *  SYNOPSIS:  	Clean up <worker> at end of main loop.
 */
void
WORK_postloop( WORKER* worker )
{
   
}