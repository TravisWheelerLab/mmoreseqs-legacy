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

/* header */
#include "_work.h"
#include "work_loop.h"

/*! FUNCTION:  	WORK_preloop()
 *  SYNOPSIS:  	Prep <worker> for main loop.
 */
void 
WORK_preloop( WORKER* worker )
{
   FILE*          fp             = NULL;
   CLOCK*         timer          = worker->timer;
   TIMES*         times          = worker->times;
   TIMES*         times_totals   = worker->times_totals;
   RESULT*        result         = worker->result;
   ALL_SCORES*    scores         = &result->scores;
   SCORES*        final          = &result->final_scores;

   
}

/*! FUNCTION:  	WORK_preiter()
 *  SYNOPSIS:  	Prep <worker> for next main loop iteration.
 */
void
WORK_preiter( WORKER* worker )
{
   FILE*          fp             = NULL;
   CLOCK*         timer          = worker->timer;
   TIMES*         times          = worker->times;
   TIMES*         times_totals   = worker->times_totals;
   RESULT*        result         = worker->result;
   ALL_SCORES*    scores         = &result->scores;
   SCORES*        final          = &result->final_scores;

   /* start timer for iteration */
   times->loop_start = CLOCK_GetTime( worker->timer );
   /* update loop value */
   worker->search_id += 1;

   /* reset timers */
   WORK_times_init( worker, times );
   /* reset scores */
   WORK_scores_init( worker, scores );
}

/*! FUNCTION:  	WORK_postiter()
 *  SYNOPSIS:  	Clean up <worker> at end of main loop iteration.
 */
void 
WORK_postiter( WORKER* worker )
{
   FILE*          fp             = NULL;
   CLOCK*         timer          = worker->timer;
   TIMES*         times          = worker->times;
   TIMES*         times_totals   = worker->times_totals;
   RESULT*        result         = worker->result;
   ALL_SCORES*    scores         = &result->scores;
   SCORES*        final          = &result->final_scores;

   /* capture total runtime for current iteration */
   times->loop_end   = CLOCK_GetTime( worker->timer );
   times->loop       = CLOCK_GetDiff( worker->timer, times->loop_start, times->loop_end );

   /* add current iteration times to totals */
   WORK_times_add( worker );
}

/*! FUNCTION:  	WORK_postloop()
 *  SYNOPSIS:  	Clean up <worker> at end of main loop.
 */
void
WORK_postloop( WORKER* worker )
{
   FILE*          fp             = NULL;
   CLOCK*         timer          = worker->timer;
   TIMES*         times          = worker->times;
   TIMES*         times_totals   = worker->times_totals;
   RESULT*        result         = worker->result;
   ALL_SCORES*    scores         = &result->scores;
   SCORES*        final          = &result->final_scores;
   
   /* capture total-ish program runtime (from clock init to end of main loop) */
   times_totals->program_end  = CLOCK_GetTime( worker->timer );
   times_totals->program      = CLOCK_GetDiff( worker->timer, times_totals->program_start, times_totals->program_end );

   /* report total times */
   REPORT_mytimeout_totals( worker, result, stdout );
}