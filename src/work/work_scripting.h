/*******************************************************************************
 *  - FILE:      work_scripting.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various
 *functions. Subroutines for executing scripts.
 *******************************************************************************/

#ifndef _WORK_SCRIPTING
#define _WORK_SCRIPTING

/*! FUNCTION:  	WORK_preloop()
 *  SYNOPSIS:  	Sets all environmental variables for values from <args>.
 */
void WORK_load_script_env_args(WORKER* worker);

#endif /* _WORK_SCRIPTING */
