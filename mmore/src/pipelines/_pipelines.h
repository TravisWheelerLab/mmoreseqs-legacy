/*******************************************************************************
 *  FILE:      _pipeline.h
 *  PURPOSE:   Pipelines and Subroutines
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _PIPELINE_MAIN_H
#define _PIPELINE_MAIN_H

/* === PIPELINE SUBROUTINES === */
#include "workflow.h"

/* === MAIN PIPELINES === */

/*! FUNCTION:  	mmore_pipeline()
 *  SYNOPSIS:  	Pipeline for MMore. 
 *                Runs MMore pipeline. Calls MMseqs pipeline, and results are piped to mmore_tail_pipeline().
 */
STATUS_FLAG 
mmore_pipeline( WORKER* worker );

/*! FUNCTION:  	generic_pipeline()
 *  SYNOPSIS:  	Generic pipeline.
 * 				   Can optionally run viterbi, forward-backward, and/or adaptive forward-backward.
 */
STATUS_FLAG 
generic_pipeline( WORKER* worker );

/*! FUNCTION:  	interactive_pipeline()
 *  SYNOPSIS:  	Interactive, tutorialized pipeline.
 *                An interactive pipeline which shows users how to work with tool and purpose of options.
 */
STATUS_FLAG 
interactive_pipeline( WORKER* worker );

/* === HELPER / INTERNAL PIPELINES === */

/*! FUNCTION:  	mmore_main_pipeline()
 *  SYNOPSIS:  	Internal pipeline for tail of MMore, post-MMseqs.
 *                Runs Adaptive Banding / Cloud Search step of MMORE pipeline.
 */
STATUS_FLAG 
mmore_main_pipeline( WORKER* worker );

/*! FUNCTION:  	index_pipeline()
 *  SYNOPSIS:     Index Pipeline: Indexes FASTA or HMM files.
 */
STATUS_FLAG 
index_pipeline( WORKER* worker );

/*! FUNCTION:     hmmbuild_pipeline()
 *  SYNOPSIS:     Pipeline that builds HMM model from FASTA file.
 *                Calls to HMMER suit to run hmmbuild.
 */
STATUS_FLAG 
hmmbuild_pipeline( WORKER* worker );

/* === DEBUGGING PIPELINES === */

/*! FUNCTION:     itest_pipeline()
 *  SYNOPSIS:     Pipeline runs integration tests. 
 *                Runs optimized and unoptimized versions of search algs and compares results.
 *                For full functionality, must be compiled in BUILD=DEBUG mode.
 */
STATUS_FLAG 
itest_pipeline( WORKER* worker );

/*! FUNCTION:  	utest_pipeline()
 *  SYNOPSIS:  	Pipeline runs unit test.
 *                For full functionality, must be compiled in DEBUG mode.
 */
STATUS_FLAG 
utest_pipeline( WORKER* worker );

/*! FUNCTION:  null_pipeline()
 *  SYNOPSIS:  Pipeline does nothing. 
 *             For testing.
 */
STATUS_FLAG 
null_pipeline( WORKER* worker );

#endif /* _PIPELINE_MAIN_H */