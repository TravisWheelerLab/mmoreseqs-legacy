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

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/hmm_profile.h"
// #include "objects/sequence.h"

/* === FUNCTIONS === */
float forward_Linear(const SEQUENCE*   query, 
                     const HMM_PROFILE* target, 
                     const int          Q, 
                     const int          T, 
                     float*             st_MX3,
                     float*             st_MX, 
                     float*             sp_MX,
                     float*             sc_final);

float backward_Linear(const SEQUENCE*    query, 
                     const HMM_PROFILE* target, 
                     const int          Q, 
                     const int          T, 
                     float*             st_MX3,
                     float*             st_MX, 
                     float*             sp_MX,
                     float*             sc_final);

#endif /* _FWDBACK_LIN_H */
