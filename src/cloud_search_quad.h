/*******************************************************************************
 *  @file cloud_search_quad.c
 *  @brief Cloud Search for Forward-Backward Pruning Alg. (QUADRATIC SPACE)
 *
 *  @synopsis
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_QUAD_H
#define _CLOUD_SEARCH_QUAD_H

void cloud_forward_Run(const SEQ* query, 
                     const HMM_PROFILE* target,
                     int Q, int T, 
                     float* st_MX, 
                     float* sp_MX, 
                     TRACEBACK* tr,
                     EDGEBOUNDS* edg,
                     float alpha, int beta );

void cloud_backward_Run(const SEQ* query, 
                     const HMM_PROFILE* target,
                     int Q, int T, 
                     float* st_MX, 
                     float* sp_MX, 
                     TRACEBACK* tr,
                     EDGEBOUNDS* edg,
                     float alpha, int beta );

#endif /* _CLOUD_SEARCH_QUAD_H */