/*******************************************************************************
 *  FILE:      work_threshold.h
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Computes and tests threshold cutoffs for inner loop.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORK_THRESHOLD
#define _WORK_THRESHOLD

/*! FUNCTION:  	WORK_threshold_viterbi()
 *  SYNOPSIS:  	Converts Viterbi natscore to e-value, then compares to threshold.
 */
void 
WORK_threshold_viterbi( WORKER* worker );

/*! FUNCTION:  	WORK_threshold_cloud()
 *  SYNOPSIS:  	Converts Cloud natscore to e-value, then compares to threshold.
 */
void 
WORK_threshold_cloud( WORKER* worker );

/*! FUNCTION:  	WORK_threshold_bound_fwdback()
 *  SYNOPSIS:  	Converts Bound Forward natscore to e-value, then compares to threshold.
 */
void 
WORK_threshold_bound_fwdback( WORKER* worker );


#endif /* _WORK_THRESHOLD */