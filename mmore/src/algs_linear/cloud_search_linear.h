/*******************************************************************************
 *    FILE:       cloud_search_linear.h
 *    PURPOSE:    Cloud Search for Forward-Backward Pruning Algorithm
 *                (Linear Space Alg)
 *
 *    AUTHOR:     Dave Rich
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_LINEAR_H
#define _CLOUD_SEARCH_LINEAR_H

/*! FUNCTION: run_Cloud_Forward_Linear()
 *  SYNOPSIS: Perform Forward part of Cloud Search Algorithm.
 *            Traverses the dynamic programming matrix antidiagonally, running the
 *            Forward algorithm, starting at the Viterbi alignment beginning.  
 *            At the end of each antidiagonal, compares each cell against the current maximum
 *            scoring cell.  If cell falls below (MAX_SCORE * alpha), then cell is removed 
 *            from search space.  Terminates when reaches the final cell in dp matrix or 
 *            all cells in current antidiag have been pruned.  
 *            Stores final edgebound data in <edg>.
 *  RETURN:   Returns <STATUS_SUCCESS> if no errors.
 */
STATUS_FLAG 
run_Cloud_Forward_Linear(  const SEQUENCE*      query,        /* query sequence */
                           const HMM_PROFILE*   target,       /* target hmm model */
                           const int            Q,            /* query length */
                           const int            T,            /* target length */
                           MATRIX_3D*           st_MX3,       /* normal state matrix */
                           MATRIX_2D*           sp_MX,        /* special state matrix */
                           const ALIGNMENT*     tr,           /* viterbi traceback */
                           EDGEROWS*      rows,         /* temporary edgebounds by-row vector */
                           EDGEBOUNDS*          edg,          /* OUTPUT: edgebounds of cloud search space */
                           CLOUD_PARAMS*        params,       /* pruning parameters */
                           float*               inner_sc,     /* OUTPUT: maximum score inside viterbi bounds */
                           float*               maxsc );      /* OUTPUT: highest score found during search */

/*
 *  FUNCTION: run_Cloud_Backward_Linear()
 *  SYNOPSIS: Perform Backward part of Cloud Search Algorithm.
 *            Traverses the dynamic programming matrix antidiagonally, running the
 *            Forward algorithm, starting at the Viterbi alignment ending.  
 *            At the end of each antidiagonal, compares each cell against the current maximum
 *            scoring cell.  If cell falls below (MAX_SCORE * alpha), then cell is removed 
 *            from search space.  Terminates when reaches the final cell in dp matrix or 
 *            all cells in current antidiag have been pruned.  
 *            Stores final edgebound data in <edg>.
 *  RETURN:   Maximum score.
 */
STATUS_FLAG 
run_Cloud_Backward_Linear(    const SEQUENCE*      query,        /* query sequence */
                              const HMM_PROFILE*   target,       /* target hmm model */
                              const int            Q,            /* query length */
                              const int            T,            /* target length */
                              MATRIX_3D*           st_MX3,       /* normal state matrix */
                              MATRIX_2D*           sp_MX,        /* special state matrix */
                              const ALIGNMENT*     tr,           /* viterbi traceback */
                              EDGEROWS*      rows,         /* temporary edgebounds by-row */
                              EDGEBOUNDS*          edg,          /* (OUTPUT) */
                              CLOUD_PARAMS*        params,       /* pruning parameters */
                              float*               inner_sc,     /* OUTPUT: maximum score inside viterbi bounds */
                              float*               max_sc );     /* highest score found during search */

#endif /* _CLOUD_SEARCH_LINEAR_H */