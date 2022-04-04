/*******************************************************************************
 *  - FILE:      work.h
 *  - DESC:    Pipelines Workflow Subroutines
 *  NOTES:    - This is being phased out.
 *******************************************************************************/

#ifndef _WORKFLOW_H
#define _WORKFLOW_H

/* generic workflow */
void WORK_main_workflow(WORKER* work);

/* debugger workflow */
void WORK_debug_workflow(WORKER* work);

/* get the alignment from sparse forward/backward */
void WORK_capture_alignment(WORKER* worker);

#endif /* _WORKFLOW_H */
