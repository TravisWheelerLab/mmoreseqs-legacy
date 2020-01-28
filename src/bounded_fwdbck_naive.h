/*******************************************************************************
 *  @file forward_backward.h
 *  @brief Function prototypes for the Forward-Backward algorithm.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

#ifndef _BOUNDED_FWDBCK_NAIVE_H
#define _BOUNDED_FWDBCK_NAIVE_H

float forward_Bounded_Naive_Run (const SEQ* query, 
                               const HMM_PROFILE* target, 
                               int Q, int T, 
                               float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                               float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                               float st_MX_cloud[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                               float* sc);

float backward_Bounded_Naive_Run (const SEQ* query, 
                                const HMM_PROFILE* target, 
                                int Q, int T, 
                                float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                                float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                                float st_MX_cloud[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                                float* sc);

#endif /* _BOUNDED_FWDBCK_NAIVE_H */
