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

void cloud_Forward_Quad(const SEQUENCE*    query, 
                        const HMM_PROFILE* target,
                        const int          Q, 
                        const int          T, 
                        MATRIX_3D*         st_MX, 
                        MATRIX_2D*         sp_MX, 
                        const ALIGNMENT*   tr,
                        EDGEBOUNDS*        edg,
                        const float        alpha, 
                        const int          beta );

void cloud_Backward_Quad(const SEQUENCE*    query, 
                        const HMM_PROFILE* target,
                        const int          Q, 
                        const int          T, 
                        MATRIX_3D*         st_MX, 
                        MATRIX_2D*         sp_MX, 
                        const ALIGNMENT*   tr,
                        EDGEBOUNDS*        edg,
                        const float        alpha, 
                        const int          beta );

#endif /* _CLOUD_SEARCH_QUAD_H */