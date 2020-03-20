/*******************************************************************************
 *
 *  FILE:    viterbi.h
 *  PURPOSE: The Viterbi Algorithm and Traceback for Sequence Alignment Search.
 *
 *  AUTHOR:  Dave Rich
 *
 *******************************************************************************/

#ifndef _VITERBI_QUAD_H
#define _VITERBI_QUAD_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/hmm_profile.h"
// #include "objects/sequence.h"

/* === FUNCTIONS === */
float viterbi_Quad(const SEQUENCE*    query,
                   const HMM_PROFILE* target,
                   const int          Q, 
                   const int          T,
                   float*             st_MX,
                   float*             sp_MX,
                   float*             sc_final);

#endif /* _VITERBI_QUAD_H */