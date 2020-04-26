/*******************************************************************************
 *  FILE:      pipeline_main.c
 *  PURPOSE:   Pipelines and Subroutines
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _PIPELINE_MAIN_H
#define _PIPELINE_MAIN_H

/* === PIPELINE SUBROUTINES === */

#include "work.h"

/* === MAIN PIPELINES === */

/* standard pipeline */
void main_pipeline( WORKER* worker );

/*
 *  FUNCTION:  test_pipeline()
 *  SYNOPSIS:  Pipeline runs integration tests. 
 *             Runs optimized and unoptimized versions of search algs and compares results.
 *             For full functionality, must be compiled in DEBUG mode.
 */
void test_pipeline( WORKER* worker );

/* MMSEQS pipeline */
void mmseqs_pipeline( WORKER* worker );

/* time trial pipeline */
void time_pipeline( WORKER* worker );

/* indexing pipeline */
void index_pipeline( WORKER* worker );

/* unit and integration testing pipeline */
void utest_pipeline( WORKER* worker );

#endif /* _PIPELINE_MAIN_H */