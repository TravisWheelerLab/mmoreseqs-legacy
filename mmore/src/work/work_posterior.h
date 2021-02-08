/*******************************************************************************
 *  FILE:      work_posterior.h
 *  PURPOSE:   Workflow Subroutines for Posterior Algorithms.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _WORK_POSTERIOR_H
#define _WORK_POSTERIOR_H

/*! FUNCTION:  WORK_posterior()
 *  SYNOPSIS:  Run full posterior, for full sequence and per-domain. 
 */
void 
WORK_posterior( WORKER* worker );

/*! FUNCTION:  WORK_posterior_sparse()
 *  SYNOPSIS:  Run posterior for full sequence using sparse matrices.
 */
void 
WORK_posterior_sparse( WORKER* worker );

/*! FUNCTION:  WORK_null1_hmm_bias()
 *  SYNOPSIS:  Compute the correction bias for the hmm model.
 */
void 
WORK_null1_hmm_bias( WORKER* worker );

/*! FUNCTION:  WORK_decode_domains()
 *  SYNOPSIS:  Find domain ranges in posterior.
 */
void
WORK_decode_domains( WORKER* worker );

/*! FUNCTION:  WORK_decode_posterior()
 *  SYNOPSIS:  Compute posterior from forward and backward matrices.
 */
void 
WORK_decode_posterior( WORKER* worker );

/*! FUNCTION:  WORK_decode_posterior()
 *  SYNOPSIS:  Compute the correction bias for the sequence.
 */
void 
WORK_null2_seq_bias( WORKER* worker );

/*! FUNCTION:  WORK_posterior_construct_scores()
 *  SYNOPSIS:  Construct final scores and evalues from results.
 */
void 
WORK_posterior_construct_scores( WORKER* worker );

#endif /* _WORK_POSTERIOR_H */