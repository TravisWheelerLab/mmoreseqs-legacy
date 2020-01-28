/*******************************************************************************
 *  @file forward_backward.h
 *  @brief Function prototypes for the Forward-Backward algorithm.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

#ifndef _FORWARD_BACKWARD_H
#define _FORWARD_BACKWARD_H

void fwdbck_Run (const SEQ* query, 
                           const HMM_PROFILE* target,
                           int Q, int T, 
                           float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                           float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                           RESULTS* res);

float forward_Run (const SEQ* query, 
                   const HMM_PROFILE* target, 
                   int Q, int T, 
                   float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                   float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                   RESULTS* res);

float backward_Run (const SEQ* query, 
                    const HMM_PROFILE* target, 
                    int Q, int T, 
                    float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                    float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                    RESULTS* res);

#endif /* _FORWARD_BACKWARD_H */
