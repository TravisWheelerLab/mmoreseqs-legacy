/*******************************************************************************
 *    - FILE:       traceback_quad.h
 *    - DESC:     Traceback for Viterbi Algorithm.
 *******************************************************************************/

#ifndef _TRACEBACK_QUAD_H
#define _TRACEBACK_QUAD_H

/*! FUNCTION:  run_Traceback_Quad()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Requires that <st_MX> and <sp_MX> are still completely filled.
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad(
    const SEQUENCE* query,     /* query sequence */
    const HMM_PROFILE* target, /* HMM model */
    const int Q,               /* query/seq length */
    const int T,               /* target/model length */
    MATRIX_3D* st_MX,          /* Normal State (Match, Insert, Delete) Matrix */
    MATRIX_2D* sp_MX,          /* Special State (J,N,B,C,E) Matrix */
    ALIGNMENT* aln);           /* OUTPUT: Traceback Alignment */

/*! FUNCTION:  run_Traceback_Quad()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Version 1: based on HMMER version. Verifies that Alignment agrees
 * with Matrix data. RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad_via_hmmer(
    const SEQUENCE* query,     /* query sequence */
    const HMM_PROFILE* target, /* HMM model */
    const int Q,               /* query/seq length */
    const int T,               /* target/model length */
    MATRIX_3D* st_MX,          /* Normal State (Match, Insert, Delete) Matrix */
    MATRIX_2D* sp_MX,          /* Special State (J,N,B,C,E) Matrix */
    ALIGNMENT* aln);           /* OUTPUT: Traceback Alignment */

/*! FUNCTION:  run_Traceback_Quad_via_cmp()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Version 2: My implementation. Takes maximum next step by finding
 * the state that fulfills equation ( <previous state> + <transition> + <score>
 * == <current state> ). Verifies that Alignment agrees with Matrix data.
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad_via_cmp(
    const SEQUENCE* query,     /* query sequence */
    const HMM_PROFILE* target, /* HMM model */
    const int Q,               /* query/seq length */
    const int T,               /* target/model length */
    MATRIX_3D* st_MX,          /* Normal State (Match, Insert, Delete) Matrix */
    MATRIX_2D* sp_MX,          /* Special State (J,N,B,C,E) Matrix */
    ALIGNMENT* aln);           /* OUTPUT: Traceback Alignment */

/*! FUNCTION:  run_Traceback_Quad_via_max()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Version 3: My implementation.  Takes maximum next step from
 * current state to previous state. Warning: No verification step. RETURN:
 * Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad_via_max(
    const SEQUENCE* query,     /* query sequence */
    const HMM_PROFILE* target, /* HMM model */
    const int Q,               /* query/seq length */
    const int T,               /* target/model length */
    MATRIX_3D* st_MX,          /* Normal State (Match, Insert, Delete) Matrix */
    MATRIX_2D* sp_MX,          /* Special State (J,N,B,C,E) Matrix */
    ALIGNMENT* aln);           /* OUTPUT: Traceback Alignment */

#endif /* _TRACEBACK_QUAD_H */
