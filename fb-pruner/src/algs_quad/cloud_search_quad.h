/*******************************************************************************
 *  FILE:       cloud_search_quad.c
 *  PURPOSE:    Cloud Search for Forward-Backward Pruning Alg. (QUADRATIC SPACE)
 *
 *  AUTHOR:     Dave Rich
 *  BUG:        Lots.
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_QUAD_H
#define _CLOUD_SEARCH_QUAD_H

/*  
 *  FUNCTION:  run_Cloud_Forward_Run()
 *  SYNOPSIS:  Perform Forward part of Cloud Search Algorithm.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Cloud_Forward_Quad(   const SEQUENCE*      query,         /* query sequence */
                              const HMM_PROFILE*   target,        /* target hmm model */
                              const int            Q,             /* query length */
                              const int            T,             /* target length */
                              MATRIX_3D*           st_MX,         /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
                              MATRIX_2D*           sp_MX,         /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
                              const ALIGNMENT*     tr,            /* viterbi traceback */ 
                              EDGEBOUND_ROWS*      rows,          /* temporary edgebounds by-row vector */
                              EDGEBOUNDS*          edg,           /* OUTPUT: edgebounds of cloud search space */
                              CLOUD_PARAMS*        params );      /* pruning parameters */

/*  
 *  FUNCTION:  cloud_backward_Run()
 *  SYNOPSIS:  Perform Backward part of Cloud Search Algorithm.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Cloud_Backward_Quad(     const SEQUENCE*      query,         /* query sequence */
                                 const HMM_PROFILE*   target,        /* target hmm model */
                                 const int            Q,             /* query length */
                                 const int            T,             /* target length */
                                 MATRIX_3D*           st_MX,         /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
                                 MATRIX_2D*           sp_MX,         /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
                                 const ALIGNMENT*     tr,            /* viterbi traceback */ 
                                 EDGEBOUND_ROWS*      rows,          /* temporary edgebounds by-row vector */
                                 EDGEBOUNDS*          edg,           /* OUTPUT: edgebounds of cloud search space */
                                 CLOUD_PARAMS*        params );      /* pruning parameters */

#endif /* _CLOUD_SEARCH_QUAD_H */