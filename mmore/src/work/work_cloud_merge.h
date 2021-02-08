/*******************************************************************************
 *  FILE:      work_cloud_merge.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for cloud search / pruned / adaptive-banding forward-backward.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORK_CLOUD_MERGE
#define _WORK_CLOUD_MERGE

/*! FUNCTION:  	WORK_cloud_merge()
 *  SYNOPSIS:  	Run "cloud merge" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 *                Caller must have run WORK_cloud_search().
 */
void 
WORK_cloud_merge_and_reorient( WORKER*  worker );

#endif /* _WORK_CLOUD_MERGE */