/*******************************************************************************
 *  - FILE:  bound_fwdbck_linear.c
 *  - DESC:  Bounded Forward/Backward Algorithm (Linear Space)
 *******************************************************************************/

#ifndef _BOUNDED_FWDBCK_LINEAR_H
#define _BOUNDED_FWDBCK_LINEAR_H

/*! FUNCTION: run_Bound_Forward_Linear()
 *  SYNOPSIS: Perform Edge-Bounded Forward step of Cloud Search Algorithm.
 *            Runs traditional Forward-Backward Algorithm, but only performs
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
int run_Bound_Forward_Linear(
    const SEQUENCE* query,      /* query sequence */
    const HMM_PROFILE* target,  /* target HMM model */
    const int Q,                /* query length */
    const int T,                /* target length */
    MATRIX_3D* restrict st_MX3, /* normal state matrix */
    MATRIX_2D* restrict sp_MX,  /* special state matrix */
    EDGEBOUNDS* edg,            /* edgebounds */
    float* sc_final);           /* (OUTPUT) final score */

/*! FUNCTION: run_Bound_Backward_Linear()
 *  SYNOPSIS: Perform Edge-Bounded Backward step of Cloud Search Algorithm.
 *            Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of
 *             (left-bound, right-bound) pairs sorted by row.
 *            Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T +
 * 1)], but is reused from cloud_forward_(). Final score produced by Backward is
 * stored in <sc_final>.
 *
 *  RETURN:   Returns the final score of the Backward Algorithm.
 */
int run_Bound_Backward_Linear(
    const SEQUENCE* query,      /* query sequence */
    const HMM_PROFILE* target,  /* target HMM model */
    const int Q,                /* query length */
    const int T,                /* target length */
    MATRIX_3D* restrict st_MX3, /* normal state matrix */
    MATRIX_2D* restrict sp_MX,  /* special state matrix */
    EDGEBOUNDS* edg,            /* edgebounds */
    float* sc_final);           /* (OUTPUT) final score */

#endif /* _BOUNDED_FWDBCK_LINEAR_H */
