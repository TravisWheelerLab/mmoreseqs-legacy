/*******************************************************************************
 *  FILE:    traceback_sparse.h
 *  PURPOSE: Traceback for Maximum Posterior
 *
 *  AUTHOR:  Dave Rich
 *******************************************************************************/

#ifndef _MAXPOST_TRACEBACK_SPARSE_H
#define _MAXPOST_TRACEBACK_SPARSE_H


/*  FUNCTION:  run_MaxExp_Traceback_Sparse_2()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Version 2: My implementation. Verifies that Alignment agrees with Matrix data.
 *
 *    RETURN:  Return <STATUS_SUCCESS> if no errors.
 */
int run_MaxPost_Traceback_Sparse(   const SEQUENCE*      query,      /* query sequence */
                                    const HMM_PROFILE*   target,     /* HMM model */
                                    const int            Q,          /* query/seq length */
                                    const int            T,          /* target/model length */
                                    MATRIX_3D_SPARSE*    st_SMX,     /* Normal State (Match, Insert, Delete) Matrix */
                                    MATRIX_2D*           sp_MX,      /* Special State (J,N,B,C,E) Matrix */
                                    EDGEBOUNDS*          edg,        /* edgebounds of sparse matrix */
                                    ALIGNMENT*           aln );      /* OUTPUT: Traceback Alignment */


#endif /*  _MAXPOST_TRACEBACK_SPARSE_H */