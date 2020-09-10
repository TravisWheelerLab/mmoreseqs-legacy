/*******************************************************************************
 *  FILE:      work.h
 *  PURPOSE:   Pipelines Workflow Subroutines
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORK_H
#define _WORK_H

/* generic workflow */
void WORK_main_workflow( WORKER*  work );

/* debugger workflow */
void WORK_debug_workflow( WORKER* work );

/* initialize dp matrices */
void WORK_init( WORKER* worker );

/* clean up data structs */
void WORK_cleanup( WORKER* worker );

/* resize dp matrices if necessary */
void WORK_reuse( WORKER* worker );

/* open worker files */
void WORK_open( WORKER* worker );

/* close worker files */
void WORK_close( WORKER* worker );

/* load or build target and query index files */
void WORK_index( WORKER* worker );

/* load query and target indexes (or build them) */
void WORK_build_target_index( WORKER* worker );

/* load query index (or build them) */
void WORK_build_query_index( WORKER* worker );

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

/* load target by file index id */
void WORK_load_target_by_index_id(  WORKER*     worker,
                                    int         index_id );

/* load target by file index id */
void WORK_load_query_by_index_id(   WORKER*     worker,
                                    int         index_id );

/* viterbi and traceback */
void WORK_viterbi_and_traceback( WORKER*  worker );

/* forward/backward */
void WORK_forward_backward( WORKER*  worker );

/* cloud search AKA pruned forward/backward */
void WORK_cloud_search( WORKER* worker );

/* compute correction bias and convert natscore -> bitscore -> pval -> eval */
void WORK_convert_scores( WORKER* worker );

/* print header for results file (default) */
void WORK_report_header( WORKER* worker );

/* print current result (default) */
void WORK_report_result_current( WORKER* worker );

/* print current result (default) */
void WORK_report_result_all( WORKER* worker );

/* print header for results file (default) */
void WORK_report_footer( WORKER* worker );

#endif /* _WORK_H */