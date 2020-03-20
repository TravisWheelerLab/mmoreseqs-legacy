/*******************************************************************************
 *  @file cloud_search.h
 *  @brief Bounded Forward-Backward Algorithm (QUADRATIC SPACE)
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _BOUNDED_FWDBCK_QUAD_H
#define _BOUNDED_FWDBCK_QUAD_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/edgebound.h"
// #include "objects/hmm_profile.h"
// #include "objects/sequence.h"

/* === FUNCTIONS === */
float bound_Forward_Quad(  const SEQUENCE*    query, 
                           const HMM_PROFILE* target,
                           const int          Q, 
                           const int          T, 
                           float*             st_MX3,
                           float*             st_MX,
                           float*             sp_MX, 
                           const EDGEBOUNDS*  edg,
                           const bool         test,
                           float              *sc_final );

float bound_Backward_Quad( const SEQUENCE*    query, 
                           const HMM_PROFILE* target,
                           const int          Q, 
                           const int          T, 
                           float*             st_MX3,
                           float*             st_MX,
                           float*             sp_MX, 
                           const EDGEBOUNDS*  edg,
                           const bool         test,
                           float              *sc_final );

#endif /* BOUNDED_FWDBCK_QUAD_H */