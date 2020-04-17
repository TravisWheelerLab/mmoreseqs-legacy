/*******************************************************************************
 *  FILE:      pipeline_main.c
 *  PURPOSE:   Pipelines and Subroutines
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _PIPELINE_MAIN_H
#define _PIPELINE_MAIN_H

/* === MAIN PIPELINES === */

/* standard pipeline */
void main_pipeline( WORKER* worker );

/* testing pipeline */
void test_pipeline( WORKER* worker );

/* MMSEQS pipeline */
void mmseqs_pipeline( WORKER* worker );

/* time trial pipeline */
void time_pipeline( WORKER* worker );

/* indexing pipeline */
void index_pipeline( WORKER* worker );

/* unit and integration testing pipeline */
void utest_pipeline( WORKER* worker );

/* === PIPELINE SUBROUTINES === */

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

#endif /* _PIPELINE_MAIN_H */