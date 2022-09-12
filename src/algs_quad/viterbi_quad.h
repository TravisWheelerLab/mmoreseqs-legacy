/*******************************************************************************
 *  - FILE:  viterbi_quad.h
 *  - DESC:  The Viterbi Algorithm and Traceback for Sequence Alignment Search.
 *******************************************************************************/

#ifndef _VITERBI_QUAD_H
#define _VITERBI_QUAD_H

/*! FUNCTION:  run_Viterbi_Quad()
 *  SYNOPSIS:  Run Viterbi Algorithm
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Viterbi_Quad(const SEQUENCE* query,     /* query sequence */
                     const HMM_PROFILE* target, /* target hmm model */
                     const int Q,               /* query length */
                     const int T,               /* target length */
                     MATRIX_3D* st_MX,          /* normal matrix */
                     MATRIX_2D* sp_MX,          /* special matrix */
                     float* sc_final);          /* final max score */

/*! FUNCTION:  run_Viterbi_Reverse_Quad()
 *  SYNOPSIS:  Run Backward version of Viterbi Algorithm (Seq-to-Profile,
 * general unoptimized) RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Viterbi_Reverse_Quad(const SEQUENCE* query,     /* query sequence */
                             const HMM_PROFILE* target, /* target hmm model */
                             const int Q,               /* query length */
                             const int T,               /* target length */
                             MATRIX_3D* st_MX,          /* normal matrix */
                             MATRIX_2D* sp_MX,          /* special matrix */
                             float* sc_final);          /* final max score */

#endif /* _VITERBI_QUAD_H */
