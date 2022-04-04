/*******************************************************************************
 *  - FILE:      fwdback_quad.h
 *  - DESC:    The Forward-Backward Algorithm for Sequence Alignment Search.
 *             (Quadratic Space Alg)
 *******************************************************************************/

#ifndef _FWDBACK_QUAD_H
#define _FWDBACK_QUAD_H

/*! FUNCTION:  run_Forward_Quad()
 *  SYNOPSIS:  Perform Forward part of Forward-Backward Algorithm.
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Forward_Quad(const SEQUENCE* query,     /* query sequence */
                     const HMM_PROFILE* target, /* target hmm model */
                     const int Q,               /* query length */
                     const int T,               /* target length */
                     MATRIX_3D* st_MX,          /* normal state matrix, dim: (
                                                   NUM_NORMAL_STATES, Q+1, T+1 ) */
                     MATRIX_2D* sp_MX,          /* special state matrix, dim: (
                                                   NUM_SPECIAL_STATES, Q+1 ) */
                     float* sc_final);          /* OUTPUT: final score */

/*! FUNCTION:   run_Backward_Quad()
 *  SYNOPSIS:   Perform Backward part of Forward-Backward Algorithm.
 *  RETURN:     Return <STATUS_SUCCESS> if no errors.
 */
int run_Backward_Quad(const SEQUENCE* query,     /* query sequence */
                      const HMM_PROFILE* target, /* target hmm model */
                      const int Q,               /* query length */
                      const int T,               /* target length */
                      MATRIX_3D* st_MX,          /* normal state matrix, dim: (
                                                    NUM_NORMAL_STATES, Q+1, T+1 ) */
                      MATRIX_2D* sp_MX,          /* special state matrix, dim: (
                                                    NUM_SPECIAL_STATES, Q+1 ) */
                      float* sc_final);          /* OUTPUT: final score */

#endif /* _FWDBACK_QUAD_H */
