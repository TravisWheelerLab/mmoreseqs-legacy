/*******************************************************************************
 *  FILE:      work.h
 *  PURPOSE:   Pipelines Workflow Subroutines
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _WORK_H
#define _WORK_H

/* generic workflow */
void WORK_workflow( WORKER*  work );

/* initialize dp matrices */
void WORK_init( WORKER* worker );

/* resize dp matrices if necessary */
void WORK_reuse( WORKER* worker );

/* load query and target indexes (or build them) */
void WORK_load_target_index( WORKER* worker );

/* load query index (or build them) */
void WORK_load_query_index( WORKER* worker );

/* output target index to file */
void WORK_output_target_index( WORKER* worker );

/* output query index to file */
void WORK_output_query_index( WORKER* worker );

/* set and verify ranges */
void WORK_set_ranges( WORKER* worker );

/* load target by file index */
void WORK_load_target_by_id( WORKER* worker,
                             int     id );

/* load query by file index id */
void WORK_load_query_by_id( WORKER* worker,
                            int     id );

/* load target by file index name */
void WORK_load_target_by_name( WORKER* worker,
                               char*   name );

/* load target by file index name */
void WORK_load_query_by_name( WORKER* worker,
                              char*   name );

/* viterbi and traceback */
void WORK_viterbi_and_traceback( WORKER*  worker );

/* forward/backward */
void WORK_forward_backward( WORKER*  worker );

/* cloud search AKA pruned forward/backward */
void WORK_cloud_search( WORKER* worker );

/* print header for results file (default) */
void WORK_print_result_header( WORKER* worker );

/* print current result (default) */
void WORK_print_result_current( WORKER* worker );

#endif /* _WORK_H */