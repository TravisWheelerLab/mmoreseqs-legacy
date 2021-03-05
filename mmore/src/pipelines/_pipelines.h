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
 *  SYNOPSIS:  	Full Pipeline for MMORE-SEQS Search. 
 *                Runs MMSEQS and MMORE stages. Uses Target and Query for MMORE and Target Profile, Target Sequence and Query db for MMSEQS as input.
 */
STATUS_FLAG 
mmore_pipeline( WORKER* worker );

/*! FUNCTION:     prep_pipeline()
 *  SYNOPSIS:     Helper for overarching MMORE-SEQS pipeline. 
 *                Prepares files for use in the MMORE-SEQS Pipeline
 */
STATUS_FLAG 
prep_pipeline( WORKER* worker );

/*! FUNCTION:     mmore_prep_pipeline()
 *  SYNOPSIS:     Full Pipeline for MMORE-SEQS Search. 
 *                Runs MMSEQS and MMORE stages. Uses prep folder as input.
 */
STATUS_FLAG 
mmore_prep_pipeline( WORKER* worker );

/*! FUNCTION:     mmore_easysearch_pipeline()
 *  SYNOPSIS:     Full Pipeline for MMORE-SEQS Search.  
 *                Uses only Target (MSA), Query (FASTA), and Prep Folder location.
 *                Performs all prepatory steps and runs search.
 */
STATUS_FLAG 
mmore_easysearch_pipeline( WORKER* worker );

/* === ALTERNATIVE SEARCH PIPELINES === */

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
 *                Runs Adaptive Banding / Cloud Search step of MMORE-SEQS pipeline.
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