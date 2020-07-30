/* ****************************************************************************************** *
 *  FILE:      cloud_search_linear.c
 *  PURPOSE:   Cloud Search for Forward-Backward Pruning Alg.
 *             (Linear Space Alg)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
/* ****************************************************************************************** */


#ifndef _CLOUD_SEARCH_LINEAR_H
#define _CLOUD_SEARCH_LINEAR_H

/* ****************************************************************************************** *
 *  FUNCTION: run_Cloud_Forward_Linear()
 *  SYNOPSIS: Perform Forward part of Cloud Search Algorithm.
 *            Traverses the dynamic programming matrix antidiagonally, running the
 *            Forward algorithm, starting at the Viterbi alignment beginning.  
 *            At the end of each antidiagonal, compares each cell against the current maximum
 *            scoring cell.  If cell falls below (MAX_SCORE * alpha), then cell is removed 
 *            from search space.  Terminates when reaches the final cell in dp matrix or 
 *            all cells in current antidiag have been pruned.  
 *            Stores final edgebound data in <edg>.
 *  RETURN:   Maximum score.
/* ****************************************************************************************** */
int run_Cloud_Forward_Linear(  const SEQUENCE*    query,        /* query sequence */
                           const HMM_PROFILE* target,       /* target hmm model */
                           const int          Q,            /* query length */
                           const int          T,            /* target length */
                           MATRIX_3D*         st_MX3,       /* normal state matrix */
                           MATRIX_2D*         sp_MX,        /* special state matrix */
                           const ALIGNMENT*   tr,           /* viterbi traceback */
                           EDGEBOUND_ROWS*    rows,         /* temporary edgebounds by-row */
                           EDGEBOUNDS*        edg,          /* (OUTPUT) */
                           CLOUD_PARAMS*      params );     /* pruning parameters */

/* ****************************************************************************************** *
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
/* ****************************************************************************************** */
int run_Cloud_Backward_Linear(    const SEQUENCE*    query,        /* query sequence */
                              const HMM_PROFILE* target,       /* target hmm model */
                              const int          Q,            /* query length */
                              const int          T,            /* target length */
                              MATRIX_3D*         st_MX3,       /* normal state matrix */
                              MATRIX_2D*         sp_MX,        /* special state matrix */
                              const ALIGNMENT*   tr,           /* viterbi traceback */
                              EDGEBOUND_ROWS*    rows,         /* temporary edgebounds by-row */
                              EDGEBOUNDS*        edg,          /* (OUTPUT) */
                              CLOUD_PARAMS*      params );     /* pruning parameters */

#endif /* _CLOUD_SEARCH_LINEAR_H */