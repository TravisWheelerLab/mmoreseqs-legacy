/*******************************************************************************
 *  @file forward_backward.h
 *  @brief Function prototypes for the Forward-Backward algorithm.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _BOUNDED_FWDBCK_NAIVE_H
#define _BOUNDED_FWDBCK_NAIVE_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/edgebound.h"
// #include "objects/hmm_profile.h"
// #include "objects/sequence.h"

/* === FUNCTIONS === */
float bound_Forward_Naive(const SEQUENCE*    query, 
                          const HMM_PROFILE* target, 
                          const int          Q, 
                          const int          T, 
                          MATRIX_3D*         st_MX, 
                          MATRIX_2D*         sp_MX,
                          MATRIX_3D*         st_MX_cloud, 
                          float*             sc_final);

float bound_Backward_Naive(const SEQUENCE*    query, 
                          const HMM_PROFILE* target, 
                          const int          Q, 
                          const int          T, 
                          MATRIX_3D*         st_MX, 
                          MATRIX_2D*         sp_MX,
                          MATRIX_3D*         st_MX_cloud, 
                          float*             sc_final);

#endif /* _BOUNDED_FWDBCK_NAIVE_H */
