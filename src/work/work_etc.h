/*******************************************************************************
 *  FILE:      work_etc.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Various uncategorized functions.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORK_ETC
#define _WORK_ETC

/*! FUNCTION:  	WORK_set_ranges()
 *  SYNOPSIS:  	Set and Verify <worker> ranges based on <args>.
 *                If <args> ranges are out-of-bounds, constrains range to valid.
 */
void 
WORK_set_ranges( WORKER*    worker );


/* initialize all scores to -inf */
void 
WORK_scores_init(    WORKER*        worker,
                     ALL_SCORES*    scores );

/* initialize all times to zero */
void 
WORK_times_init(  WORKER*     worker,
                  TIMES*      times );

/* initialize all times to zero */
void 
WORK_times_add(   WORKER*     worker );

#endif /* _WORK_ETC */