/*******************************************************************************
 *  - FILE:  cloud_search_quad.c
 *  - DESC:  Cloud Search for Forward-Backward Pruning Alg. (Quadratic Space)
 *******************************************************************************/

#ifndef _CLOUD_SEARCH_QUAD_H
#define _CLOUD_SEARCH_QUAD_H

/*! FUNCTION:  run_Cloud_Forward_Run()
 *  SYNOPSIS:  Perform Forward part of Cloud Search Algorithm.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Cloud_Forward_Quad(
    const SEQUENCE* query,     /* query sequence */
    const HMM_PROFILE* target, /* target hmm model */
    const int Q,               /* query length */
    const int T,               /* target length */
    MATRIX_3D*
        st_MX, /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
    MATRIX_2D*
        sp_MX,             /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
    const ALIGNMENT* tr,   /* viterbi traceback */
    EDGEBOUND_ROWS* rows,  /* temporary edgebounds by-row vector */
    VECTOR_INT* lb_vec[3], /* temporary left-bound vectors for pruning */
    VECTOR_INT* rb_vec[3], /* temporary right-bound vectors for pruning */
    EDGEBOUNDS* edg,       /* OUTPUT: edgebounds of cloud search space */
    CLOUD_PARAMS* params,  /* pruning parameters */
    float* inner_sc,       /* OUTPUT: maximum score inside viterbi bounds */
    float* maxsc);         /* OUTPUT: maximum score */

/*! FUNCTION:  run_Cloud_Backward_Run()
 *  SYNOPSIS:  Perform Backward part of Cloud Search Algorithm.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Cloud_Backward_Quad(
    const SEQUENCE* query,     /* query sequence */
    const HMM_PROFILE* target, /* target hmm model */
    const int Q,               /* query length */
    const int T,               /* target length */
    MATRIX_3D*
        st_MX, /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
    MATRIX_2D*
        sp_MX,             /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
    const ALIGNMENT* tr,   /* viterbi traceback */
    EDGEBOUND_ROWS* rows,  /* temporary edgebounds by-row vector */
    VECTOR_INT* lb_vec[3], /* temporary left-bound vectors for pruning */
    VECTOR_INT* rb_vec[3], /* temporary right-bound vectors for pruning */
    EDGEBOUNDS* edg,       /* OUTPUT: edgebounds of cloud search space */
    CLOUD_PARAMS* params,  /* pruning parameters */
    float* inner_sc,       /* OUTPUT: maximum score inside viterbi bounds */
    float* maxsc);         /* OUTPUT: maximum score */

#endif /* _CLOUD_SEARCH_QUAD_H */
