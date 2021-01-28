/*******************************************************************************
 *  FILE:      pipeline.h
 *  PURPOSE:   Pipelines and Subroutines
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _PIPELINE_MAIN_H
#define _PIPELINE_MAIN_H

/* === PIPELINE SUBROUTINES === */
#include "work.h"

/* === MAIN PIPELINES === */

/*! FUNCTION:  	main_pipeline()
 *  SYNOPSIS:  	Main pipeline:
 * 				      Can optionally run viterbi, forward-backward, and/or pruned forward-backward.
 */
void main_pipeline( WORKER* worker );


/*! FUNCTION:  	mmseqs_pipeline()
 *  SYNOPSIS:  	Pipeline for MMSeqs-Plus: 
 *                   Runs Cloud Search using results from main MMSeqs Search..	
 */
void mmseqs_pipeline( WORKER* worker );


/*! FUNCTION:  	mmseqs_pipeline_plus()
 *  SYNOPSIS:  	Pipeline for MMSeqs-Plus: 
 *                   Runs Cloud Search using results from main MMSeqs Search..	
 */
void mmseqs_plus_pipeline( WORKER* worker );


/*
 *  FUNCTION:  	time_pipeline()
 *  SYNOPSIS:  	Time Trial Pipeline: 
 *                   Runs all algorithms and reports times.
 */
void time_pipeline( WORKER* worker );

/*
 *  FUNCTION:  	index_pipeline()
 *  SYNOPSIS:   Index Pipeline: Indexes FASTA or HMM files.
 */
void index_pipeline( WORKER* worker );


/*
 *  FUNCTION:  hmmbuild_pipeline()
 *  SYNOPSIS:  Pipeline that builds HMM model from FASTA file.
 */
void hmmbuild_pipeline( WORKER* worker );

/*
 *  FUNCTION:  itest_pipeline()
 *  SYNOPSIS:  Pipeline runs integration tests. 
 *             Runs optimized and unoptimized versions of search algs and compares results.
 *             For full functionality, must be compiled in DEBUG mode.
 */
void itest_pipeline( WORKER* worker );

/*
 *  FUNCTION:  	utest_pipeline()
 *  SYNOPSIS:  	Unit Test Pipeline: runs all unit tests.
 */
void utest_pipeline( WORKER* worker );

/*
 *  FUNCTION:  vizualization_pipeline()
 *  SYNOPSIS:  Pipeline runs viterbi, forward/backward, and cloud search.  
 *             Output visualizations for python scripts.
 */
void vizualization_pipeline( WORKER* worker );

/*
 *  FUNCTION:  null_pipeline()
 *  SYNOPSIS:  Pipeline does nothing.
 */
void null_pipeline( WORKER* worker );

#endif /* _PIPELINE_MAIN_H */