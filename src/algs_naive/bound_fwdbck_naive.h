/*******************************************************************************
 *  - FILE:     forward_backward.h
 *  - DESC:   Function prototypes for the Forward-Backward algorithm.
 *******************************************************************************/

#ifndef _BOUNDED_FWDBCK_NAIVE_H
#define _BOUNDED_FWDBCK_NAIVE_H

/*! FUNCTION: run_Bound_Forward_Naive()
 *  SYNOPSIS: Perform Forward part of Forward-Backward Algorithm.
 *
 *  ARGS:      <query>        query sequence,
 *             <target>       HMM model,
 *             <Q>            query length,
 *             <T>            target length,
 *             <st_MX>        Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>        Special State (J,N,B,C,E) Matrix,
 *             <st_MX_cloud>  Boolean Matrix of Cloud Bounds,
 *             <sc_final>     (OUTPUT) Final Score
 *
 *  RETURN:    Returns the Forward Score.
 */
float run_Bound_Forward_Naive(const SEQUENCE* query,
                              const HMM_PROFILE* target,
                              const int Q,
                              const int T,
                              MATRIX_3D* st_MX,
                              MATRIX_2D* sp_MX,
                              MATRIX_2D* st_MX_cloud,
                              float* sc_final);

/*! FUNCTION: run_Bound_Backward_Naive()
 *  SYNOPSIS: Perform Backward part of Forward-Backward Algorithm.
 *
 *  ARGS:      <query>        query sequence,
 *             <target>       HMM model,
 *             <Q>            query length,
 *             <T>            target length,
 *             <st_MX>        Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>        Special State (J,N,B,C,E) Matrix,
 *             <st_MX_cloud>  Boolean Matrix of Cloud Bounds,
 *             <sc_final>     (OUTPUT) Final Score
 *
 * RETURN:     Returns the Forward Score.
 */
float run_Bound_Backward_Naive(const SEQUENCE* query,
                               const HMM_PROFILE* target,
                               const int Q,
                               const int T,
                               MATRIX_3D* st_MX,
                               MATRIX_2D* sp_MX,
                               MATRIX_2D* st_MX_cloud,
                               float* sc_final);

#endif /* _BOUNDED_FWDBCK_NAIVE_H */
