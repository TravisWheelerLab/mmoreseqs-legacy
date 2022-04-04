/*******************************************************************************
 *  - FILE:      fwdback_linear.c
 *  - DESC:    The Forward-Backward Algorithm for Sequence Alignment Search.
 *             (Linear Space)
 *******************************************************************************/

#ifndef _FWDBACK_LIN_H
#define _FWDBACK_LIN_H

/** FUNCTION: run_Forward_Linear()
 *  SYNOPSIS: Perform Forward part of Forward-Backward Algorithm. (Linear Space
 * Implementation)
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence,
 *             <target>    HMM model,
 *             <Q>         query length,
 *             <T>         target length,
 *             <st_MX3>    Normal State (Match, Insert, Delete) 3D-Matrix,
 *                         [ NUM_NORMAL_STATES x 3 x (Q+T+1) ]
 *             <st_MX>     Normal State (Match, Insert, Delete) 3D-Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *                         [ NUM_SPECIAL_STATES x (Q+1) ]
 *             <sc_final>  OUTPUT: Score in Log-Bits
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors
 */
int run_Forward_Linear(const SEQUENCE* query,
                       const HMM_PROFILE* target,
                       const int Q,
                       const int T,
                       MATRIX_3D* st_MX3,
                       MATRIX_2D* sp_MX,
                       float* sc_final);

/** FUNCTION:  run_Backward_Linear()
 *  SYNOPSIS:  Perform Backward part of Forward-Backward Algorithm. (LINEAR ALG)
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence,
 *             <target>    HMM model,
 *             <Q>         query length,
 *             <T>         target length,
 *             <st_MX3>    Normal State (Match, Insert, Delete) Matrix,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <sc_final>  OUTPUT: Score in Log-Bits
 *
 *  RETURN:     <STATUS_SUCCESS> if no errors.
 */
int run_Backward_Linear(const SEQUENCE* query,
                        const HMM_PROFILE* target,
                        const int Q,
                        const int T,
                        MATRIX_3D* st_MX3,
                        MATRIX_2D* sp_MX,
                        float* sc_final);

#endif /* _FWDBACK_LIN_H */
