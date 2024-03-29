/*******************************************************************************
 *  - FILE:  posterior_quad.h
 *  - DESC:  The Posterior Probability and Optimal Alignment.
 *******************************************************************************/

#ifndef _POSTERIOR_QUAD_H
#define _POSTERIOR_QUAD_H

/*! FUNCTION:  run_Posterior_Quad()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward
 * <st_MX_bck>. Compute the Posterior Probability by multiplying probabilities
 * (added in log space) of Forward and Backward. Results stored in supplied
 * <st_MX_post> (can override input matrices). NOTE: Modeled after
 * <p7_Decoding()> and <p7_DomainDecoding()>. RETURN:    Return <STATUS_SUCCESS>
 * if no errors.
 */
int run_Posterior_Quad(
    SEQUENCE* q_seq,       /* query sequence */
    HMM_PROFILE* t_prof,   /* target hmm model */
    int Q,                 /* query length */
    int T,                 /* target length */
    HMM_BG* bg,            /* hmm background model */
    EDGEBOUNDS* edg,       /* edgebounds */
    MATRIX_3D* st_MX_fwd,  /* normal state matrix for forward */
    MATRIX_2D* sp_MX_fwd,  /* special state matrix for forward */
    MATRIX_3D* st_MX_bck,  /* normal state matrix for backward */
    MATRIX_2D* sp_MX_bck,  /* special state matrix for backward */
    MATRIX_3D* st_MX_post, /* OUTPUT: normal state matrix for posterior */
    MATRIX_2D* sp_MX_post, /* OUTPUT: special state matrix for posterior */
    DOMAIN_DEF* dom_def);  /* OUTPUT: domain data */

/*! FUNCTION:  run_Decode_Normal_Posterior_Quad()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices (log-space) to create
 * special state posterior into <...post> (normal-space). Can store <...post>
 * matrix in <...bck>. Again, input should be in log-space, output matrix is in
 * normal-space (needed for computing null2 score). NOTE: Modeled after
 * <p7_Decoding()> and <p7_DomainDecoding()> RETURN:    Return <STATUS_SUCCESS>
 * if no errors.
 */
int run_Decode_Posterior_Quad(
    SEQUENCE* q_seq,        /* query sequence */
    HMM_PROFILE* t_prof,    /* target hmm model */
    int Q,                  /* full query size */
    int T,                  /* full target size */
    RANGE* Q_range,         /* query range */
    RANGE* T_range,         /* target range */
    MATRIX_3D* st_MX_fwd,   /* normal state matrix for forward */
    MATRIX_2D* sp_MX_fwd,   /* special state matrix for forward */
    MATRIX_3D* st_MX_bck,   /* normal state matrix for backward */
    MATRIX_2D* sp_MX_bck,   /* special state matrix for backward */
    MATRIX_3D* st_MX_post,  /* OUTPUT: normal state matrix for posterior */
    MATRIX_2D* sp_MX_post); /* OUTPUT: normal state matrix for posterior */

/*! FUNCTION:  run_Decode_Posterior_Quad()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state
 * posterior into <...post>. Can store matrix in <...fwd> or <...bck>. NOTE:
 * Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Decode_Special_Posterior_Quad(
    SEQUENCE* q_seq,        /* query sequence */
    HMM_PROFILE* t_prof,    /* target hmm model */
    int Q,                  /* query length */
    int T,                  /* target length */
    MATRIX_2D* sp_MX_fwd,   /* special state matrix for forward */
    MATRIX_2D* sp_MX_bck,   /* special state matrix for backward */
    MATRIX_2D* sp_MX_post); /* OUTPUT: special state matrix for posterior */

/*! FUNCTION:  run_Null2_ByExpectation_Quad()
 *  SYNOPSIS:  Modeled after HMMER p7_GNull2_ByExpectation().
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Null2_ByExpectation_Quad(
    SEQUENCE* query,       /* query sequence */
    HMM_PROFILE* target,   /* target hmm model */
    int Q,                 /* query length */
    int T,                 /* target length */
    RANGE* Q_range,        /* query range */
    RANGE* T_range,        /* target range */
    MATRIX_3D* st_MX_post, /* posterior normal matrix */
    MATRIX_2D* sp_MX_post, /* posterior special matrix */
    DOMAIN_DEF* dom_def);  /* OUTPUT: domain def's null2_sc vector */

/*! FUNCTION:  run_Null2_ByExpectation_Old()
 *  SYNOPSIS:  Modeled after HMMER p7_GNull2_ByExpectation().
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Null2_ByExpectation_Quad_Old(
    SEQUENCE* query,       /* query sequence */
    HMM_PROFILE* target,   /* target hmm model */
    int Q,                 /* query length */
    int T,                 /* target length */
    MATRIX_3D* st_MX_post, /* posterior normal matrix */
    MATRIX_2D* sp_MX_post, /* posterior special matrix */
    DOMAIN_DEF* dom_def);  /* OUTPUT: domain def's null2_sc vector */

#endif /* _POSTERIOR_QUAD_H */
