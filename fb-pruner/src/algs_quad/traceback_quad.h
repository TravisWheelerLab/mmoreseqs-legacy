/*******************************************************************************
 *    FILE:       traceback_quad.h
 *    PURPOSE:    Traceback for Viterbi Algorithm.
 *
 *    AUTHOR:     Dave Rich
 *******************************************************************************/

#ifndef _TRACEBACK_QUAD_H
#define _TRACEBACK_QUAD_H

/*
 *  FUNCTION:  run_Traceback_Quad()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad(    const SEQUENCE*     query,       /* query sequence */
                           const HMM_PROFILE*  target,      /* HMM model */
                           const int           Q,           /* query/seq length */
                           const int           T,           /* target/model length */
                           MATRIX_3D*          st_MX,       /* Normal State (Match, Insert, Delete) Matrix */
                           MATRIX_2D*          sp_MX,       /* Special State (J,N,B,C,E) Matrix */
                           ALIGNMENT*          aln );       /* OUTPUT: Traceback Alignment */

/*
 *  FUNCTION:  run_Traceback_Quad_2()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Version 2: My implementation. Verifies that Alignment agrees with Matrix data.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad_2(     const SEQUENCE*     query,       /* query sequence */
                              const HMM_PROFILE*  target,      /* HMM model */
                              const int           Q,           /* query/seq length */
                              const int           T,           /* target/model length */
                              MATRIX_3D*          st_MX,       /* Normal State (Match, Insert, Delete) Matrix */
                              MATRIX_2D*          sp_MX,       /* Special State (J,N,B,C,E) Matrix */
                              ALIGNMENT*          aln );       /* OUTPUT: Traceback Alignment */

/*
 *  FUNCTION:  run_Traceback_Quad_3()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Version 3: My implementation. Verifies that Alignment agrees with Matrix data.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad_3(     const SEQUENCE*     query,       /* query sequence */
                              const HMM_PROFILE*  target,      /* HMM model */
                              const int           Q,           /* query/seq length */
                              const int           T,           /* target/model length */
                              MATRIX_3D*          st_MX,       /* Normal State (Match, Insert, Delete) Matrix */
                              MATRIX_2D*          sp_MX,       /* Special State (J,N,B,C,E) Matrix */
                              ALIGNMENT*          aln );       /* OUTPUT: Traceback Alignment */

#endif /* _TRACEBACK_QUAD_H */