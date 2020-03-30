/*******************************************************************************
 *  FILE:      seq_to_hmm_hmmer.h
 *  PURPOSE:   Converts a SEQUENCE to HMM_PROFILE.
 *             Uses the Easel and HMMER libraries.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _SEQ_TO_HMM_HMMER_H
#define _SEQ_TO_HMM_HMMER_H

/* converts SEQUENCE to HMM_PROFILE model */
void SEQUENCE_to_HMM_PROFILE( SEQUENCE*      seq, 
                              HMM_PROFILE*   prof );

#endif /* _SEQ_TO_HMM_HMMER_H */
