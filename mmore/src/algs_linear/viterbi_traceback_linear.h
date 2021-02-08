/*******************************************************************************
 *  FILE:    traceback_lin.c
 *  PURPOSE: Traceback for Viterbi Algorithm (LINEAR SPACE).
 *
 *  AUTHOR:  Dave Rich   
 *******************************************************************************/

#ifndef _TRACEBACK_LIN_H
#define _TRACEBACK_LIN_H

/*
 *  FUNCTION:  run_Traceback_Linear()
 *  SYNOPSIS:  Selects the default method of run_Traceback_Lin() from the available methods.
 *             Requires that <st_MX> and <sp_MX> are still completely filled.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Linear(   const SEQUENCE*     query,       /* query sequence */
                            const HMM_PROFILE*  target,      /* HMM model */
                            const int           Q,           /* query/seq length */
                            const int           T,           /* target/model length */
                            MATRIX_3D*          st_MX3,      /* Normal State (Match, Insert, Delete) Matrix */
                            MATRIX_2D*          sp_MX,       /* Special State (J,N,B,C,E) Matrix */
                            ALIGNMENT*          aln );       /* OUTPUT: Traceback Alignment */

/*
 *  FUNCTION:  run_Traceback_Linear_via_cmp()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment. (Linear Space)
 *             Version 2: My implementation. Takes maximum next step by finding the state that fulfills equation ( <previous state> + <transition> + <score> == <current state> ).
 *             Verifies that Alignment agrees with Matrix data.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Linear_via_cmp(   const SEQUENCE*     query,       /* query sequence */
                                    const HMM_PROFILE*  target,      /* HMM model */
                                    const int           Q,           /* query/seq length */
                                    const int           T,           /* target/model length */
                                    MATRIX_3D*          st_MX,       /* Normal State (Match, Insert, Delete) Matrix */
                                    MATRIX_2D*          sp_MX,       /* Special State (J,N,B,C,E) Matrix */
                                    ALIGNMENT*          aln );       /* OUTPUT: Traceback Alignment */


#endif /* _TRACEBACK_LIN_H */