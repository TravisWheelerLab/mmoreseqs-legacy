/*******************************************************************************
 *  FILE:      pipeline_main.c
 *  PURPOSE:   Main Cloud Search Pipeline
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _PIPELINE_MAIN_H
#define _PIPELINE_MAIN_H

/* standard pipeline */
void main_pipeline( ARGS* args );

/* testing pipeline */
void test_pipeline( ARGS* args );

/* MMSEQS pipeline */
void mmseqs_pipeline( ARGS* args );

/* time trial pipeline */
void time_pipeline( ARGS* args );

/* indexing pipeline */
void index_pipeline( ARGS* args );

/* unit and integration testing pipeline */
void utest_pipeline( ARGS* args );

#endif /* _PIPELINE_MAIN_H */