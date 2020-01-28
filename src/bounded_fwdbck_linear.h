/*******************************************************************************
 *  @file bounded_fwdbck_linear.c
 *  @brief Cloud Forward-Backward Algorithm (Linear Space Implementation)
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/


#ifndef _BOUNDED_FWDBCK_LINEAR_H
#define _BOUNDED_FWDBCK_LINEAR_H

float forward_bounded_Run3(const SEQ* query, 
                           const HMM_PROFILE* target,
                           int Q, int T, 
                           float* st_MX3, 
                           float* st_MX,
                           float* sp_MX, 
                           EDGEBOUNDS* edg,
                           float *sc_final );

float backward_bounded_Run3(const SEQ* query, 
                           const HMM_PROFILE* target,
                           int Q, int T, 
                           float* st_MX3, 
                           float* st_MX,
                           float* sp_MX, 
                           EDGEBOUNDS* edg,
                           float *sc_final );


#endif /* _BOUNDED_FWDBCK_LINEAR_H */