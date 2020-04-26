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
 *  FUNCTION:  cloud_Forward_Run()
 *  SYNOPSIS:  Perform Forward part of Cloud Search Algorithm.
 */
float cloud_Forward_Quad(const SEQUENCE*      query,         /* query sequence */
                        const HMM_PROFILE*   target,        /* target model */
                        const int            Q,             /* query length */
                        const int            T,             /* target length */
                        MATRIX_3D*           st_MX,         /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
                        MATRIX_2D*           sp_MX,         /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
                        const ALIGNMENT*     tr,            /* viterbi traceback */ 
                        EDGEBOUNDS*          edg,           /* OUTPUT: edgebounds of cloud search space */
                        float                alpha,         /* PARAM: x-drop threshold */
                        int                  beta );        /* PARAM: free passes */

/*  
 *  FUNCTION: cloud_backward_Run()
 *  SYNOPSIS: Perform Backward part of Cloud Search Algorithm.
 */
float cloud_Backward_Quad(const SEQUENCE*     query,         /* query sequence */
                        const HMM_PROFILE*   target,        /* target model */
                        const int            Q,             /* query length */
                        const int            T,             /* target length */
                        MATRIX_3D*           st_MX,         /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
                        MATRIX_2D*           sp_MX,         /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
                        const ALIGNMENT*     tr,            /* viterbi traceback */ 
                        EDGEBOUNDS*          edg,           /* OUTPUT: edgebounds of cloud search space */
                        float                alpha,         /* PARAM: x-drop threshold */
                        int                  beta );        /* PARAM: free passes */

#endif /* _CLOUD_SEARCH_QUAD_H */