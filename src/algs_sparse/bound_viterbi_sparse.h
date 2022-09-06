/*******************************************************************************
 *  - FILE:  bounded_viterbi_sparse.c
 *  - DESC:  Bounded Viterbi Algorithm
 *           (Sparse Matrix)
 *******************************************************************************/

#ifndef _BOUND_VITERBI_SPARSE_H
#define _BOUND_VITERBI_SPARSE_H

/*! FUNCTION: run_Bound_Viterbi_Sparse()
 *  SYNOPSIS: Perform Edge-Bounded Viterbi.
 *            Runs traditional Viterbi Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of
 *             (left-bound, right-bound) pairs sorted by row.
 *            Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T +
 * 1)], but is reused from cloud_forward_(). Final score produced by Forward is
 * stored in <sc_final>.
 *
 *  RETURN:   Returns the final score of the Forward Algorithm.
 */
STATUS_FLAG
run_Bound_Viterbi_Sparse(
    const SEQUENCE* query,                 /* query sequence */
    const HMM_PROFILE* target,             /* target HMM model */
    const int Q,                           /* query length */
    const int T,                           /* target length */
    const EDGEBOUNDS* edg,                 /* edgebounds */
    const RANGE* dom_range,                /* (OPTIONAL) domain range for computing fwd/bck on specific
                                              domain. If NULL, computes complete fwd/bck. */
    MATRIX_3D_SPARSE* restrict st_SMX_vit, /* normal state matrix */
    MATRIX_2D* restrict sp_MX_vit,         /* special state matrix */
    float* sc_final);                      /* (OUTPUT) final score */

#endif /* _BOUND_VITERBI_SPARSE_H */
