/*******************************************************************************
 *  FILE:      seq_to_model.c
 *  PURPOSE:   Converts a SEQUENCE to HMM_PROFILE.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _SEQ_TO_MODEL_H
#define _SEQ_TO_MODEL_H

/* === FUNCTIONS === */

void SEQ_to_PROF( SEQUENCE*      seq, 
                  HMM_PROFILE*   prof );

HMM_BG* HMM_BG_Create();

#endif /* _SEQ_TO_MODEL_H */