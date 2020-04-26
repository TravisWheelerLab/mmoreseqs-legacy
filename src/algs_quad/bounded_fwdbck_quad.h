/*******************************************************************************
 *  FILE:      bounded_fwdbck_quad.c
 *  PURPOSE:   Bounded Forward-Backward Algorithm (QUADRATIC SPACE)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _BOUNDED_FWDBCK_QUAD_H
#define _BOUNDED_FWDBCK_QUAD_H


float bound_Forward_Quad(  const SEQUENCE*      query, 
                           const HMM_PROFILE*   target,
                           const int            Q, 
                           const int            T, 
                           MATRIX_3D*           st_MX,
                           MATRIX_2D*           sp_MX, 
                           EDGEBOUNDS*          edg,
                           float*               sc_final );

float bound_Backward_Quad( const SEQUENCE*      query, 
                           const HMM_PROFILE*   target,
                           const int            Q, 
                           const int            T, 
                           MATRIX_3D*           st_MX,
                           MATRIX_2D*           sp_MX, 
                           EDGEBOUNDS*          edg,
                           float*               sc_final );

#endif /* BOUNDED_FWDBCK_QUAD_H */