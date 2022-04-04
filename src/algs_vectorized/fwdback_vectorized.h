/*******************************************************************************
 *  - FILE:       fwdback_vectorized.c
 *  SYNOPSIS:   The Forward-Backward Algorithm for Sequence Alignment Search.
 *              (Linear O(Q), SIMD Vectorized)
 *  NOTES:
 *******************************************************************************/

#ifndef _FWDBACK_VEC_H
#define _FWDBACK_VEC_H

/*!  FUNCTION:    run_Forward_Vectorized()
 *   SYNOPSIS:    Perform Forward step of Forward-Backward Algorithm.
 *                Vectorized Implementation.
 *   RETURN:      <STATUS_SUCCESS> if no errors
 */
STATUS_FLAG
run_Forward_Vectorized(const SEQUENCE* query,
                       const HMM_PROFILE* target,
                       const int Q,
                       const int T,
                       MATRIX_3D* st_MX3,
                       MATRIX_2D* sp_MX,
                       float* sc_final);

/*!  FUNCTION:    run_Backward_Vectorized()
 *   SYNOPSIS:    Perform Backward step of Forward-Backward Algorithm.
 *                Vectorized Implementation.
 *   RETURN:      <STATUS_SUCCESS> if no errors
 */
STATUS_FLAG
run_Backward_Vectorized(const SEQUENCE* query,
                        const HMM_PROFILE* target,
                        const int Q,
                        const int T,
                        MATRIX_3D* st_MX3,
                        MATRIX_2D* sp_MX,
                        float* sc_final);

#endif /* _FWDBACK_VEC_H */
