/*******************************************************************************
 *  - FILE:      work_fwdback.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various
 *functions. Subroutines for forward-backward.
 *******************************************************************************/

#ifndef _WORK_FWDBACK
#define _WORK_FWDBACK

/*! FUNCTION:  	WORK_forward_backward()
 *  SYNOPSIS:  	Run "cloud search" step of pruned forward/backward (aka
 * adaptive-band forward/backward). Depends on <task> settings in <worker>.
 */
void WORK_forward_backward(WORKER* worker);

#endif /* _WORK_FWDBACK */
