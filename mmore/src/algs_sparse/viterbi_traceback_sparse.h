/*******************************************************************************
 *  FILE:    viterbi_traceback_sparse.h
 *  PURPOSE: Traceback for Viterbi Algorithm.
 *
 *  AUTHOR:  Dave Rich
 *******************************************************************************/

#ifndef _VITERBI_TRACEBACK_SPARSE_H
#define _VITERBI_TRACEBACK_SPARSE_H


/*  FUNCTION:  run_Viterbi_Traceback_Sparse()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *
 *    RETURN:  Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Viterbi_Traceback_Sparse(   const SEQUENCE*      query,      /* query sequence */
                                const HMM_PROFILE*   target,     /* HMM model */
                                const int            Q,          /* query/seq length */
                                const int            T,          /* target/model length */
                                MATRIX_3D_SPARSE*    st_SMX,     /* Normal State (Match, Insert, Delete) Matrix */
                                MATRIX_2D*           sp_MX,      /* Special State (J,N,B,C,E) Matrix */
                                EDGEBOUNDS*          edg,        /* edgebounds of sparse matrix */
                                ALIGNMENT*           aln );      /* OUTPUT: Traceback Alignment */


#endif /* _VITERBI_TRACEBACK_SPARSE_H */