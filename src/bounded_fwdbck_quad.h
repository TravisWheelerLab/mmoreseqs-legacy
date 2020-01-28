/*******************************************************************************
 *  @file cloud_search.h
 *  @brief Bounded Forward-Backward Algorithm (QUADRATIC SPACE)
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/


#ifndef _BOUNDED_FWDBCK_QUAD_H
#define _BOUNDED_FWDBCK_QUAD_H

float forward_bounded_Run(const SEQ* query, 
                        const HMM_PROFILE* target,
                        int Q, int T, 
                        float* st_MX, 
                        float* sp_MX, 
                        EDGEBOUNDS* edg, 
                        float *sc_final );

float backward_bounded_Run(const SEQ* query, 
                        const HMM_PROFILE* target,
                        int Q, int T, 
                        float* st_MX, 
                        float* sp_MX,
                        EDGEBOUNDS* edg,
                        float *sc_final );

#endif /* BOUNDED_FWDBCK_QUAD_H */