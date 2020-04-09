/*******************************************************************************
 *  FILE:      fwdback_linear.c
 *  SYNOPSIS:  The Forward-Backward Algorithm for Sequence Alignment Search.
 *             (Linear Space Alg)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _FWDBACK_LIN_H
#define _FWDBACK_LIN_H

float forward_Linear(const SEQUENCE*   query, 
                     const HMM_PROFILE* target, 
                     const int          Q, 
                     const int          T, 
                     MATRIX_3D*         st_MX3,
                     MATRIX_3D*         st_MX, 
                     MATRIX_2D*         sp_MX,
                     float*             sc_final);

float backward_Linear(const SEQUENCE*    query, 
                     const HMM_PROFILE* target, 
                     const int          Q, 
                     const int          T, 
                     MATRIX_3D*         st_MX3, 
                     MATRIX_3D*         st_MX, 
                     MATRIX_2D*         sp_MX,
                     float*             sc_final);

#endif /* _FWDBACK_LIN_H */
