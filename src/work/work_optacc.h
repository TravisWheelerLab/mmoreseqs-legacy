/*******************************************************************************
 *  - FILE:   work_optacc.c
 *  - DESC:   Pipelines Workflow Subroutines.
 *            WORK interfaces between pipeline WORKER object and various functions. Subroutines for optimal accuracy.
 *******************************************************************************/

#ifndef _WORK_OPTACC
#define _WORK_OPTACC

/*! FUNCTION:  	WORK_optimal_accuracy()
 *  SYNOPSIS:  	Run optimal accuracy of posterior.
 *                Depends on <task> settings in <worker>.
 */
void WORK_optimal_accuracy(WORKER* worker);

/*! FUNCTION:  	WORK_optacc_traceback()
 *  SYNOPSIS:  	Run traceback of optimal accuracy / posterior.
 *                Depends on <task> settings in <worker>.
 */
void WORK_optacc_traceback(WORKER* worker);

#endif /* _WORK_OPTACC */
