/*******************************************************************************
 *  FILE:      bounded_fwdbck_linear.c
 *  PURPOSE:   Bounded Forward/Backward Algorithm (LINEAR SPACE)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _BOUNDED_FWDBCK_LINEAR_H
#define _BOUNDED_FWDBCK_LINEAR_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/edgebound.h"
// #include "objects/hmm_profile.h"
// #include "objects/sequence.h"

/* === FUNCTIONS === */
/* ------------------------------------------------------------------------------------------ *
 *
 *  FUNCTION: bound_Forward_Linear()
 *  SYNOPSIS: Perform Edge-Bounded Forward step of Cloud Search Algorithm.
 *            Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of 
 *             (left-bound, right-bound) pairs sorted by row.  
 *            Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_().
 *            Final score produced by Forward is stored in <sc_final>.
 *
 *  ARGS:      <query>        Query sequence, 
 *             <target>       HMM model,
 *             <Q>            Query length, 
 *             <T>            Target length,
 *             <st_MX3>       Normal State (Match, Insert, Delete) Matrix (Linear Space),
 *             <st_MX>        Normal State (Match, Insert, Delete) Matrix (Quadratic Space),
 *             <sp_MX>        Special State (J,N,B,C,E) Matrix,
 *             <edg>          Edgebounds (stored row-wise)
 *             <test>         Toggles DEBUG output
 *             <sc_final>     (OUTPUT) Final Score 
 *
 *  RETURN:    Returns the final score of the Forward Algorithm.
 *
 * ------------------------------------------------------------------------------------------- */
float bound_Forward_Linear(const SEQUENCE*    query, 
                           const HMM_PROFILE* target,
                           const int          Q, 
                           const int          T, 
                           float*             st_MX3,
                           float*             st_MX,
                           float*             sp_MX, 
                           const EDGEBOUNDS*  edg,
                           const bool         test,
                           float              *sc_final );

/* ------------------------------------------------------------------------------------------ *
 *
 *  FUNCTION: bound_Backward_Linear()
 *  SYNOPSIS: Perform Edge-Bounded Backward step of Cloud Search Algorithm.
 *            Runs traditional Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of 
 *             (left-bound, right-bound) pairs sorted by row.  
 *            Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_Run3().
 *            Final score produced by Forward is stored in <sc_final>.
 *
 *  ARGS:      <query>        Query sequence, 
 *             <target>       HMM model,
 *             <Q>            Query length, 
 *             <T>            Target length,
 *             <st_MX3>       Normal State (Match, Insert, Delete) Matrix (Linear Space),
 *             <st_MX>        Normal State (Match, Insert, Delete) Matrix (Quadratic Space),
 *             <sp_MX>        Special State (J,N,B,C,E) Matrix,
 *             <edg>          Edgebounds (stored row-wise)
 *             <test>         Toggles DEBUG output
 *             <sc_final>     Final Score (OUTPUT)
 *
 *  RETURN:    Returns the final score of the Backward Algorithm.
 *
 * ------------------------------------------------------------------------------------------- */
float bound_Backward_Linear(  const SEQUENCE*    query, 
                              const HMM_PROFILE* target,
                              const int          Q, 
                              const int          T, 
                              float*             st_MX3,
                              float*             st_MX,
                              float*             sp_MX, 
                              const EDGEBOUNDS*  edg,
                              const bool         test,
                              float              *sc_final );

#endif /* _BOUNDED_FWDBCK_LINEAR_H */