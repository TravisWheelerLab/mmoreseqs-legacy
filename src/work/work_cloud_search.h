/*******************************************************************************
 *  FILE:      work_cloud.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for cloud search / pruned / adaptive-banding forward-backward.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORK_CLOUD_SEARCH
#define _WORK_CLOUD_SEARCH

/*! FUNCTION:  	WORK_cloud_search()
 *  SYNOPSIS:  	Run "cloud search" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_cloud_search( WORKER*  worker );

/*! FUNCTION:  	WORK_cloud_search_linear()
 *  SYNOPSIS:  	Run linear-space "cloud search" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_cloud_search_linear( WORKER*  worker );

/*! FUNCTION:  	WORK_cloud_search_quadratic()
 *  SYNOPSIS:  	Run quadratic-space "cloud search" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_cloud_search_quadratic( WORKER*  worker );

#endif /* _WORK_CLOUD_SEARCH */