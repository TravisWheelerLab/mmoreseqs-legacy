/*******************************************************************************
 *  - FILE:      work_loop.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various
 *functions. Subroutines for inner search loop.
 *******************************************************************************/

#ifndef _WORK_LOOP
#define _WORK_LOOP

/*! FUNCTION:  	WORK_preloop()
 *  SYNOPSIS:  	Prep <worker> for main loop.
 */
void WORK_preloop(WORKER* worker);

/*! FUNCTION:  	WORK_preiter()
 *  SYNOPSIS:  	Prep <worker> for next main loop iteration.
 */
void WORK_preiter(WORKER* worker);

/*! FUNCTION:  	WORK_postiter()
 *  SYNOPSIS:  	Clean up <worker> at end of main loop iteration.
 */
void WORK_postiter(WORKER* worker);

/*! FUNCTION:  	WORK_postloop()
 *  SYNOPSIS:  	Clean up <worker> at end of main loop.
 */
void WORK_postloop(WORKER* worker);

#endif /* _WORK_LOOP */
