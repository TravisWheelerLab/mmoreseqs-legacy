/*******************************************************************************
 *  FILE:      fwdback_quad.h
 *  PURPOSE:   The Forward-Backward Algorithm for Sequence Alignment Search.
 *             (Quadratic Space Alg)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _FWDBACK_QUAD_H
#define _FWDBACK_QUAD_H


int run_Forward_Quad(   const SEQUENCE*   query, 
                        const HMM_PROFILE* target, 
                        const int          Q, 
                        const int          T, 
                        MATRIX_3D*         st_MX, 
                        MATRIX_2D*         sp_MX,
                        float*             sc_final );


int run_Backward_Quad(  const SEQUENCE*    query, 
                        const HMM_PROFILE* target, 
                        const int          Q, 
                        const int          T, 
                        MATRIX_3D*         st_MX, 
                        MATRIX_2D*         sp_MX,
                        float*             sc_final );

#endif /* _FWDBACK_QUAD_H */
