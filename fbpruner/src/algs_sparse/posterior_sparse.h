/*******************************************************************************
 *     FILE:   posterior_sparse.h
 *  PURPOSE:   The Posterior Probability and Optimal Alignment.
 *
 *   AUTHOR:   Dave Rich
 *******************************************************************************/

#ifndef _POSTERIOR_SPARSE_H
#define _POSTERIOR_SPARSE_H

/*! FUNCTION:  run_Posterior_Sparse()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward <st_MX_bck>.
 *             Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in supplied <st_MX_post> (can override input matrices).
 *             NOTE: Modeled after <p7_ByPosteriorHeuristics()>.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Posterior_Sparse(   SEQUENCE*               q_seq,            /* query sequence */
                        HMM_PROFILE*            t_prof,           /* target hmm model */
                        int                     Q,                /* query length */
                        int                     T,                /* target length */
                        HMM_BG*                 bg,               /* hmm background model */
                        EDGEBOUNDS*             edg,              /* edgebounds */
                        MATRIX_3D_SPARSE*       st_MX_fwd,        /* normal state matrix for forward */
                        MATRIX_2D*              sp_MX_fwd,        /* special state matrix for forward */
                        MATRIX_3D_SPARSE*       st_MX_bck,        /* normal state matrix for backward */
                        MATRIX_2D*              sp_MX_bck,        /* special state matrix for backward */
                        MATRIX_3D_SPARSE*       st_MX_post,       /* OUTPUT: normal state matrix for posterior */
                        MATRIX_2D*              sp_MX_post,       /* OUTPUT: special state matrix for posterior */     
                        DOMAIN_DEF*             dom_def );        /* OUTPUT: domain data */


/*! FUNCTION:  run_Decode_Domains()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward <st_MX_bck>.
 *             Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in supplied <st_MX_post> (can override input matrices).
 *             NOTE: Modeled after <p7_DomainDecoding()>.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Domains(     SEQUENCE*               q_seq,            /* query sequence */
                        HMM_PROFILE*            t_prof,           /* target hmm model */
                        int                     Q,                /* query length */
                        int                     T,                /* target length */
                        MATRIX_2D*              sp_MX_fwd,        /* special state matrix for forward */
                        MATRIX_2D*              sp_MX_bck,        /* special state matrix for backward */
                        MATRIX_2D*              sp_MX_post,       /* OUTPUT: special state matrix for posterior */     
                        DOMAIN_DEF*             dom_def );        /* OUTPUT: domain data */


/*! FUNCTION:  run_Decode_Normal_Posterior_Sparse()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store matrix in <...fwd> or <...bck>.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Posterior_Sparse(  SEQUENCE*            q_seq,            /* query sequence */
                              HMM_PROFILE*         t_prof,           /* target hmm model */
                              int                  Q,                /* query length */
                              int                  T,                /* target length */
                              EDGEBOUNDS*          edg,              /* edgebounds */
                              RANGE*               in_dom_range,     /* domain range */
                              MATRIX_3D_SPARSE*    st_MX_fwd,        /* normal state matrix for forward */
                              MATRIX_2D*           sp_MX_fwd,        /* special state matrix for forward */
                              MATRIX_3D_SPARSE*    st_MX_bck,        /* normal state matrix for backward */
                              MATRIX_2D*           sp_MX_bck,        /* special state matrix for backward */
                              MATRIX_3D_SPARSE*    st_MX_post,       /* OUTPUT: normal state matrix for posterior */
                              MATRIX_2D*           sp_MX_post );     /* OUTPUT: normal state matrix for posterior */


/*! FUNCTION:  run_CompositionBias_Sparse()
 *  SYNOPSIS:  Computes the Null2 Composition.
 *             Modeled after HMMER p7_GNull2_ByExpectation().
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Null2_ByExpectation_Sparse(  SEQUENCE*            query,            /* query sequence */
                                 HMM_PROFILE*         target,           /* target hmm model */
                                 int                  Q,                /* query length */
                                 int                  T,                /* target length */
                                 RANGE*               Q_range,          /* query span of bounds */
                                 RANGE*               T_range,          /* target span of bounds */
                                 EDGEBOUNDS*          edg,              /* edgebounds */
                                 MATRIX_3D_SPARSE*    st_SMX_post,      /* posterior normal matrix */
                                 MATRIX_2D*           sp_MX_post,       /* posterior special matrix */
                                 DOMAIN_DEF*          dom_def,          /* OUTPUT: domain def's null2_sc vector */
                                 float*               compo_bias );     /* OUTPUT: Null2 composition bias */

#endif /* _POSTERIOR_SPARSE_H */
