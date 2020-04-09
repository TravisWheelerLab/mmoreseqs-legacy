/*******************************************************************************
 *  @file forward_backward.h
 *  @brief Function prototypes for the Forward-Backward algorithm. (QUADRATIC ALGORITHM)
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _FWDBACK_QUAD_H
#define _FWDBACK_QUAD_H


float forward_Quad(const SEQUENCE*   query, 
                  const HMM_PROFILE* target, 
                  const int          Q, 
                  const int          T, 
                  MATRIX_3D*         st_MX, 
                  MATRIX_2D*         sp_MX,
                  float*             sc_final);

float backward_Quad( const SEQUENCE*    query, 
                     const HMM_PROFILE* target, 
                     const int          Q, 
                     const int          T, 
                     MATRIX_3D*         st_MX, 
                     MATRIX_2D*         sp_MX,
                     float*             sc_final);

#endif /* _FWDBACK_QUAD_H */
