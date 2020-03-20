/*******************************************************************************
 *  @file cloud_search_quad.c
 *  @brief Cloud Search for Forward-Backward Pruning Alg. (QUADRATIC SPACE)
 *
 *  @synopsis
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_QUAD_H
#define _CLOUD_SEARCH_QUAD_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/edgebound.h"
// #include "objects/hmm_profile.h"
// #include "objects/sequence.h"
// #include "objects/alignment.h"

/* === FUNCTIONS === */
void cloud_Forward_Quad(const SEQUENCE*    query, 
                        const HMM_PROFILE* target,
                        const int          Q, 
                        const int          T, 
                        float*             st_MX, 
                        float*             sp_MX, 
                        const ALIGNMENT*   tr,
                        EDGEBOUNDS*        edg,
                        const float        alpha, 
                        const int          beta );

void cloud_Backward_Quad(const SEQUENCE*    query, 
                        const HMM_PROFILE* target,
                        const int          Q, 
                        const int          T, 
                        float*             st_MX, 
                        float*             sp_MX, 
                        const ALIGNMENT*   tr,
                        EDGEBOUNDS*        edg,
                        const float        alpha, 
                        const int          beta );

#endif /* _CLOUD_SEARCH_QUAD_H */