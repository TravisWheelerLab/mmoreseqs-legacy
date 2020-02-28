/*******************************************************************************
 *  @file cloud_search_linear.h
 *  @brief Cloud Search for Forward-Backward Pruning Alg. (LINEAR SPACE)
 *
 *  @synopsis
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_LINEAR_ROWS_H
#define _CLOUD_SEARCH_LINEAR_ROWS_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/edgebound.h"
// #include "objects/hmm_profile.h"
// #include "objects/sequence.h"
// #include "objects/alignment.h"

/* === FUNCTIONS === */
void cloud_forward_Run3_rows( const SEQUENCE*    query, 
                              const HMM_PROFILE* target,
                              const int          Q, 
                              const int          T, 
                              float*             st_MX, 
                              float*             st_MX3,
                              float*             sp_MX, 
                              const ALIGNMENT*   tr,
                              EDGEBOUNDS*        edg,
                              const float        alpha, 
                              const int          beta,
                              const bool         test );

void cloud_backward_Run3_rows(const SEQUENCE*    query, 
                              const HMM_PROFILE* target,
                              const int          Q, 
                              const int          T, 
                              float*             st_MX, 
                              float*             st_MX3,
                              float*             sp_MX, 
                              const ALIGNMENT*   tr,
                              EDGEBOUNDS*        edg,
                              const float        alpha, 
                              const int          beta,
                              const bool         test );


#endif /* _CLOUD_SEARCH_LINEAR_ROWS_H */