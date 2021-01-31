/*******************************************************************************
 *  FILE:      seq_to_model.c
 *  PURPOSE:   Converts a SEQUENCE to HMM_PROFILE model.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _SEQ_TO_MODEL_H
#define _SEQ_TO_MODEL_H

/* === FUNCTIONS === */

/* converts single SEQUENCE to HMM_PROFILE model */
void SEQUENCE_to_HMM_PROFILE( SEQUENCE*      seq, 
                              HMM_PROFILE*   prof );

/* create HMM from SEQUENCE */
void HMM_PROFILE_From_Seq( HMM_PROFILE*   prof,
                           SEQUENCE*      seq );

/* set background composition of sequence */
void HMM_PROFILE_Set_Composition( HMM_PROFILE* prof );

/* run simulation to calibrate e-value parameters of model */
void HMM_PROFILE_Calibrate( HMM_PROFILE* prof );

#endif /* _SEQ_TO_MODEL_H */