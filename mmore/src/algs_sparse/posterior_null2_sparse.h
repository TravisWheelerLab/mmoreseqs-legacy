/*******************************************************************************
 *  FILE:      posterior_null2_sparse.h
 *  PURPOSE:   Computer the Null2 Composition Bias of the Posterior Probability.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _POSTERIOR_NULL2_SPARSE_H
#define _POSTERIOR_NULL2_SPARSE_H

/*! FUNCTION:  run_Null2_ByExpectation_Sparse()
 *  SYNOPSIS:  Computes the Null2 Bias Composition.
 *             Modeled after HMMER p7_GNull2_ByExpectation().
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
STATUS_FLAG
run_Null2_ByExpectation_Sparse(  SEQUENCE*            query,            /* query sequence */
                                 HMM_PROFILE*         target,           /* target hmm model */
                                 int                  Q,                /* query length */
                                 int                  T,                /* target length */
                                 EDGEBOUNDS*          edg,              /* edgebounds */
                                 RANGE*               q_range,          /* OPTIONAL: query range of edgebound cells */
                                 RANGE*               t_range,          /* OPTIONAL: target range of edgebound cells */
                                 RANGE*               dom_range,        /* OPTIONAL: domain range in query */
                                 MATRIX_3D_SPARSE*    st_SMX_post,      /* posterior normal matrix */
                                 MATRIX_2D*           sp_MX_post,       /* posterior special matrix */
                                 DOMAIN_DEF*          dom_def,          /* OUTPUT: domain def's null2_sc vector */
                                 float*               compo_bias );     /* OUTPUT: Null2 composition bias */

#endif /* _POSTERIOR_NULL2_SPARSE_H */