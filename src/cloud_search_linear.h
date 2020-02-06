/*******************************************************************************
 *  @file cloud_search_linear.h
 *  @brief Cloud Search for Forward-Backward Pruning Alg. (LINEAR SPACE)
 *
 *  @synopsis
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_LINEAR_H
#define _CLOUD_SEARCH_LINEAR_H

void cloud_forward_Run3(const SEQ* query, 
                     const HMM_PROFILE* target,
                     int Q, int T, 
                     float* st_MX, 
                     float* st_MX3,
                     float* sp_MX, 
                     TRACEBACK* tr,
                     EDGEBOUNDS* edg,
                     float alpha, int beta,
                     bool test );

void cloud_backward_Run3(const SEQ* query, 
                        const HMM_PROFILE* target,
                        int Q, int T, 
                        float* st_MX, 
                        float* st_MX3,
                        float* sp_MX, 
                        TRACEBACK* tr,
                        EDGEBOUNDS* edg,
                        float alpha, int beta,
                        bool test );


#endif /* _CLOUD_SEARCH_LINEAR_H */