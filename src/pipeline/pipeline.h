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

/* load query and target indexes (or build them) */
void WORK_load_indexes( WORKER* worker );

/* load target and query from file */
void WORK_load_target_query( WORKER*  worker );

/* initialize dp matrices */
void WORK_init( WORKER* worker );

/* resize dp matrices if necessary */
void WORK_reuse( WORKER* worker );

/* viterbi and traceback */
void WORK_viterbi_and_traceback( WORKER*  worker );

/* forward/backward */
void WORK_forward_backward( WORKER*  worker );

#endif /* _PIPELINE_MAIN_H */