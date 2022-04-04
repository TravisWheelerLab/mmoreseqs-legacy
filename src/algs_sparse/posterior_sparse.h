/*******************************************************************************
 *  - FILE:   posterior_sparse.h
 *  - DESC:    The Posterior Probability and Optimal Alignment.
 *******************************************************************************/

#ifndef _POSTERIOR_SPARSE_H
#define _POSTERIOR_SPARSE_H

/*! FUNCTION:  run_Posterior_Sparse()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward
 * <st_MX_bck>. Compute the Posterior Probability by multiplying probabilities
 * (added in log space) of Forward and Backward. Results stored in supplied
 * <st_MX_post> (can override input matrices). NOTE: Modeled after
 * <p7_ByPosteriorHeuristics()>. RETURN:    Return <STATUS_SUCCESS> if no
 * errors.
 */
STATUS_FLAG
run_Posterior_Sparse(
    SEQUENCE* q_seq,              /* query sequence */
    HMM_PROFILE* t_prof,          /* target hmm model */
    int Q,                        /* query length */
    int T,                        /* target length */
    HMM_BG* hmm_bg,               /* hmm background model */
    EDGEBOUNDS* edg,              /* edgebounds */
    MATRIX_3D_SPARSE* st_SMX_fwd, /* normal state matrix for forward */
    MATRIX_2D* sp_MX_fwd,         /* special state matrix for forward */
    MATRIX_3D_SPARSE* st_SMX_bck, /* normal state matrix for backward */
    MATRIX_2D* sp_MX_bck,         /* special state matrix for backward */
    MATRIX_3D_SPARSE*
        st_SMX_post,       /* OUTPUT: normal state matrix for posterior */
    MATRIX_2D* sp_MX_post, /* OUTPUT: special state matrix for posterior */
    MATRIX_3D_SPARSE*
        st_SMX_opt, /* OUTPUT: normal state matrix for optimal accuracy */
    MATRIX_2D*
        sp_MX_opt,       /* OUTPUT: special state matrix for optimal accuracy */
    RESULT* result,      /* OUPUT: full cloud results */
    DOMAIN_DEF* dom_def, /* OUTPUT: domain data */
    CLOCK* timer,
    TIMES* times,
    bool is_run_domains); /* if run posterior on domains */

/*! FUNCTION:  run_Decode_Normal_Posterior_Sparse()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state
 * posterior into <...post>. Can store matrix in <...fwd> or <...bck>. NOTE:
 * Modeled after <p7_Decoding()>. RETURN:    Return <STATUS_SUCCESS> if no
 * errors.
 */
STATUS_FLAG
run_Decode_Posterior_Sparse(
    SEQUENCE* q_seq,             /* query sequence */
    HMM_PROFILE* t_prof,         /* target hmm model */
    int Q,                       /* query length */
    int T,                       /* target length */
    EDGEBOUNDS* edg,             /* edgebounds */
    RANGE* dom_range,            /* domain range */
    MATRIX_3D_SPARSE* st_MX_fwd, /* normal state matrix for forward */
    MATRIX_2D* sp_MX_fwd,        /* special state matrix for forward */
    MATRIX_3D_SPARSE* st_MX_bck, /* normal state matrix for backward */
    MATRIX_2D* sp_MX_bck,        /* special state matrix for backward */
    MATRIX_3D_SPARSE*
        st_MX_post,         /* OUTPUT: normal state matrix for posterior */
    MATRIX_2D* sp_MX_post); /* OUTPUT: normal state matrix for posterior */

/*! FUNCTION:  run_Decode_Domains()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward
 * <st_MX_bck>. Compute the Posterior Probability by multiplying probabilities
 * (added in log space) of Forward and Backward. Results stored in supplied
 * <st_MX_post> (can override input matrices). NOTE: Modeled after
 * <p7_DomainDecoding()>. RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
STATUS_FLAG
run_Decode_Domains(SEQUENCE* q_seq,      /* query sequence */
                   HMM_PROFILE* t_prof,  /* target hmm model */
                   int Q,                /* query length */
                   int T,                /* target length */
                   EDGEBOUNDS* edg,      /* edgebounds */
                   MATRIX_2D* sp_MX_fwd, /* special state matrix for forward */
                   MATRIX_2D* sp_MX_bck, /* special state matrix for backward */
                   DOMAIN_DEF* dom_def); /* OUTPUT: domain data */

#endif /* _POSTERIOR_SPARSE_H */
