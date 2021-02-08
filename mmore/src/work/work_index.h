/*******************************************************************************
 *  FILE:      work_index.h
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Build, load and retrieve file indexes.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORK_INDEX
#define _WORK_INDEX

/*! FUNCTION:  	WORK_load_indexes()
 *  SYNOPSIS:  	Load or build target and query index files <t_index> for <q_index>.
 *                Stored in <worker>.
 */
void 
WORK_load_indexes( WORKER* worker );

/*! FUNCTION:  	WORK_load_indexes_by_id()
 *  SYNOPSIS:  	Load or build target and query index files <t_index> for <q_index>.
 */
void 
WORK_load_indexes_by_id( WORKER* worker );

/*! FUNCTION:  	WORK_load_indexes_by_name()
 *  SYNOPSIS:  	Load or build target index <t_index> and query index <q_index> for <worker>.
 *                Then sorts index by name.
 */
void 
WORK_load_indexes_by_name( WORKER* worker );

/*! FUNCTION:  	WORK_build_target_index()
 *  SYNOPSIS:  	Build target index <t_index> for <worker>.
 */
void 
WORK_build_target_index( WORKER* worker );

/*! FUNCTION:  	WORK_build_query_index()
 *  SYNOPSIS:  	Build query index <q_index> for <worker>.
 */
void 
WORK_build_query_index( WORKER* worker );

/*! FUNCTION:  	WORK_load_target_index()
 *  SYNOPSIS:  	Load target index <t_index> for <worker>.
 */
void 
WORK_load_target_index( WORKER* worker );

/*! FUNCTION:  	WORK_load_query_index()
 *  SYNOPSIS:  	Load query index <q_index> for <worker>.
 */
void 
WORK_load_query_index( WORKER* worker );

/*! FUNCTION:  	WORK_output_target_index()
 *  SYNOPSIS:  	Write target index out to <t_indexpath>.
 */
void 
WORK_output_target_index( WORKER* worker );

/*! FUNCTION:  	WORK_output_target_index()
 *  SYNOPSIS:  	Write query index out to <q_indexpath>.
 */
void 
WORK_output_query_index( WORKER* worker );

#endif /* _WORK_INDEX */