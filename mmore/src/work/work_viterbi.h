/*******************************************************************************
 *  FILE:      work_viterbi.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for Viterbi.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORK_VITERBI
#define _WORK_VITERBI

/*! FUNCTION:  	WORK_viterbi_and_traceback()
 *  SYNOPSIS:  	Run Viterbi algorithm.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_viterbi_and_traceback( WORKER*  worker );

/*! FUNCTION:  	WORK_viterbi_sparse()
 *  SYNOPSIS:  	Run Viterbi algorithm.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_viterbi_sparse( WORKER*  worker );

/*! FUNCTION:  	WORK_viterbi_traceback_sparse()
 *  SYNOPSIS:  	Run Viterbi Traceback algorithm.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_viterbi_traceback_sparse( WORKER*  worker );


#endif /* _WORK_VITERBI */