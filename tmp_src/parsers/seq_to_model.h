/*******************************************************************************
 *  FILE:      seq_to_model.c
 *  PURPOSE:   Converts a SEQUENCE to HMM_PROFILE.
 *             Uses the Easel and HMMER libraries.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _SEQ_TO_MODEL_H
#define _SEQ_TO_MODEL_H

/* === FUNCTIONS === */

void SEQ_to_PROF( SEQUENCE*      seq, 
                  HMM_PROFILE*   prof );


#endif /* _SEQ_TO_MODEL_H */