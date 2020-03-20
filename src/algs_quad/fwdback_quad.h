/*******************************************************************************
 *  @file forward_backward.h
 *  @brief Function prototypes for the Forward-Backward algorithm. (QUADRATIC ALGORITHM)
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _FWDBACK_QUAD_H
#define _FWDBACK_QUAD_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/hmm_profile.h"
// #include "objects/sequence.h"

/* === FUNCTIONS === */
float forward_Quad(const SEQUENCE*   query, 
                  const HMM_PROFILE* target, 
                  const int          Q, 
                  const int          T, 
                  float*             st_MX, 
                  float*             sp_MX,
                  float*             sc_final);

float backward_Quad( const SEQUENCE*    query, 
                     const HMM_PROFILE* target, 
                     const int          Q, 
                     const int          T, 
                     float*             st_MX, 
                     float*             sp_MX,
                     float*             sc_final);

#endif /* _FWDBACK_QUAD_H */
