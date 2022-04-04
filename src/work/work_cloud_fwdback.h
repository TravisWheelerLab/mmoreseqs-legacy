/*******************************************************************************
 *  - FILE:      work_cloud_fwdback.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various
 *functions. Subroutines for cloud search / pruned / adaptive-banding
 *forward-backward.
 *******************************************************************************/

#ifndef _WORK_CLOUD_FWDBACK
#define _WORK_CLOUD_FWDBACK

/*! FUNCTION:  	WORK_bound_forward_backward()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward
 * (aka adaptive-band forward/backward). Depends on <task> settings in <worker>.
 *                Caller must have run WORK_cloud_merge_and_reorient().
 */
void WORK_bound_fwdback(WORKER* worker);

/*! FUNCTION:  	WORK_bound_forward_backward_linear()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward
 * (aka adaptive-band forward/backward). Uses the linear space matrix
 * implementation. Depends on <task> settings in <worker>.
 */
void WORK_bound_fwdback_linear(WORKER* worker);

/*! FUNCTION:  	WORK_bound_forward_backward_linear()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward
 * (aka adaptive-band forward/backward). Uses the quadratic space matrix
 * implementation. Depends on <task> settings in <worker>.
 */
void WORK_bound_fwdback_quadratic(WORKER* worker);

/*! FUNCTION:  	WORK_bound_forward_backward_sparse()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward
 * (aka adaptive-band forward/backward). Uses the sparse matrix implementation.
 *                Depends on <task> settings in <worker>.
 */
void WORK_bound_fwdback_sparse(WORKER* worker);

#endif /* _WORK_CLOUD_FWDBACK */
