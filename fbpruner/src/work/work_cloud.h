/*******************************************************************************
 *  FILE:      work_cloud.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for cloud search / pruned / adaptive-banding forward-backward.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORK_CLOUD
#define _WORK_CLOUD

/*! FUNCTION:  	WORK_cloud_search()
 *  SYNOPSIS:  	Run "cloud search" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_cloud_search( WORKER*  worker );

/*! FUNCTION:  	WORK_cloud_merge()
 *  SYNOPSIS:  	Run "cloud merge" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_cloud_merge( WORKER*  worker );

/*! FUNCTION:  	WORK_bound()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_bound_forward_backward( WORKER* worker );

#endif /* _WORK_CLOUD */