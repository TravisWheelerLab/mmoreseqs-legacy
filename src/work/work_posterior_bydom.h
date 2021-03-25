/*******************************************************************************
 *  FILE:      work_posterior.h
 *  PURPOSE:   Workflow Subroutines for Posterior Algorithms.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _WORK_POSTERIOR_BYDOM_H
#define _WORK_POSTERIOR_BYDOM_H

/*! FUNCTION:  WORK_posterior()
 *  SYNOPSIS:  Run full posterior, for full sequence and per-domain. 
 */
void 
WORK_posterior_bydom( WORKER* worker );

/*! FUNCTION:  WORK_posterior_sparse()
 *  SYNOPSIS:  Run posterior for full sequence using sparse matrices.
 */
void 
WORK_posterior_sparse_bydom( WORKER* worker );

/*! FUNCTION:  WORK_null1_hmm_bias_bydom()
 *  SYNOPSIS:  Compute the correction bias for the hmm model.
 *             Only for domain region.
 */
void 
WORK_null1_hmm_bias_bydom( WORKER* worker );

/*! FUNCTION:  WORK_decode_posterior()
 *  SYNOPSIS:  Compute posterior from forward and backward matrices.
 */
void 
WORK_decode_posterior_bydom( WORKER* worker );

/*! FUNCTION:  WORK_decode_posterior()
 *  SYNOPSIS:  Compute the correction bias for the sequence.
 */
void 
WORK_null2_seq_bias_bydom( WORKER* worker );

#endif /* _WORK_POSTERIOR_BYDOM_H */