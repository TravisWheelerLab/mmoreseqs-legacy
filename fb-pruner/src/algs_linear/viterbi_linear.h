/*******************************************************************************
 *
 *  FILE:    viterbi.h
 *  PURPOSE: The Viterbi Algorithm and Traceback for Sequence Alignment Search.
 *
 *  AUTHOR:  Dave Rich
 *
 *******************************************************************************/

#ifndef _VITERBI_LINEAR_H
#define _VITERBI_LINEAR_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/hmm_profile.h"
// #include "objects/sequence.h"

/* === FUNCTIONS === */
int viterbi_Linear(	const SEQUENCE*    query,
               		const HMM_PROFILE* target,
                    const int          Q, 
                    const int          T, 
                    MATRIX_3D*         st_MX,
                    MATRIX_2D*         sp_MX,
                    float*             sc_final);

#endif /* _VITERBI_QUAD_H */