/*******************************************************************************
 *  FILE:      work_loader.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Loads query sequences and target hmm profiles.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORK_LOADER
#define _WORK_LOADER

/*! FUNCTION:  	WORK_load_target_by_id()
 *  SYNOPSIS:  	Loads mmseqs input m8 file into <results_in>, located at <mmseqs_res_filepath>.
 *                Verifies that search index range does not exceed bounds of list.
 */
void 
WORK_load_mmseqs_list( WORKER* worker );

/*! FUNCTION:  	WORK_load_target_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_target_by_id(    WORKER*     worker,
                           int         id );

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_query_by_id(     WORKER*     worker,
                           int         id );

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <name> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_target_by_name(  WORKER*    worker,
                           char*      name );

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <name> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_query_by_name(   WORKER*     worker,
                           char*       name );

/*! FUNCTION:  	WORK_load_target_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_target_by_index_id(    WORKER*     worker,
                                 int         index_id );

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_query_by_index_id(     WORKER*     worker,
                                 int         index_id );

#endif /* _WORK_LOADER */