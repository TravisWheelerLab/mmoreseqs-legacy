/*******************************************************************************
 *  - FILE:  bounded_fwdbck_sparse_test.c
 *  - DESC:  Bounded Forward/Backward Algorithm
 *            (Sparse Space)
 *******************************************************************************/

#ifndef _BOUND_FWDBCK_SPARSE_TEST_H
#define _BOUND_FWDBCK_SPARSE_TEST_H

/*! FUNCTION:  run_Bound_Forward_Sparse()
 *  SYNOPSIS:  Perform Edge-Bounded Forward step of Cloud Search Algorithm.
 *             Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of
 *             (left-bound, right-bound) pairs sorted by row.
 *             Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T +
 * 1)], but is reused from cloud_forward_(). Final score produced by Forward is
 * stored in <sc_final>.
 *
 *    RETURN:  Returns the final score of the Forward Algorithm.
 */
STATUS_FLAG
run_Bound_Forward_Sparse_TEST(
    const SEQUENCE* query,                 /* query sequence */
    const HMM_PROFILE* target,             /* target HMM model */
    const int Q,                           /* query length */
    const int T,                           /* target length */
    MATRIX_3D_SPARSE* restrict st_SMX_fwd, /* normal state matrix */
    MATRIX_2D* restrict sp_MX_fwd,         /* special state matrix */
    const EDGEBOUNDS* edg,                 /* edgebounds */
    const RANGE*
        dom_range,    /* (OPTIONAL) domain range for computing fwd/bck on specific
                      domain. If NULL, computes complete fwd/bck. */
    float* sc_final); /* (OUTPUT) final score */

/*! FUNCTION:  run_Bound_Backward_Sparse()
 *  SYNOPSIS:  Perform Edge-Bounded Backward step of Cloud Search Algorithm.
 *             Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of
 *             (left-bound, right-bound) pairs sorted by row.
 *             Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T +
 * 1)], but is reused from cloud_forward_(). Final score produced by Backward is
 * stored in <sc_final>.
 *
 *    RETURN:  Returns the final score of the Backward Algorithm.
 */
STATUS_FLAG
run_Bound_Backward_Sparse_TEST(
    const SEQUENCE* query,                 /* query sequence */
    const HMM_PROFILE* target,             /* target HMM model */
    const int Q,                           /* query length */
    const int T,                           /* target length */
    MATRIX_3D_SPARSE* restrict st_SMX_bck, /* normal state matrix */
    MATRIX_2D* restrict sp_MX_bck,         /* special state matrix */
    const EDGEBOUNDS* edg,                 /* edgebounds */
    const RANGE*
        dom_range,    /* (OPTIONAL) domain range for computing fwd/bck on specific
                      domain. If NULL, computes complete fwd/bck. */
    float* sc_final); /* (OUTPUT) final score */

#endif /* _BOUND_FWDBCK_SPARSE_TEST_H */
