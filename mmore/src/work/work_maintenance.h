/*******************************************************************************
 *  FILE:      work_maintenance.h
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Initialization and Cleanup.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORK_MAINTENANCE
#define _WORK_MAINTENANCE

/*! FUNCTION:  	WORK_init()
 *  SYNOPSIS:  	Initialize <worker> parts that are not handled by WORKER_Create().
 *                Allocate data structs according to settings in its <args> and <tasks>.
 */
void 
WORK_init( WORKER* worker );

/*! FUNCTION:  	WORK_reuse()
 *  SYNOPSIS:  	Resize and reallocate data structs in <worker> for problem size.
 */
void 
WORK_reuse( WORKER* worker );

/*! FUNCTION:  	WORK_cleanup()
 *  SYNOPSIS:  	Clean up and free data allocated by WORK_init().
 */
void 
WORK_cleanup( WORKER* worker );


#endif /* _WORK_MAINTENANCE */