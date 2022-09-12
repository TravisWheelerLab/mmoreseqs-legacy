/*******************************************************************************
 *  - FILE:  bounded_fwdbck_quad.c
 *  - DESC:  Bounded Forward-Backward Algorithm (Quadratic Space)
 *******************************************************************************/

#ifndef _BOUNDED_VITERBI_QUAD_H
#define _BOUNDED_VITERBI_QUAD_H

/*
 *  FUNCTION: 	run_Bound_Viterbi_Quad()
 *  SYNOPSIS: 	Perform Edge-Bounded Viterbi, for recovering alignment.
 *
 *  RETURN: 	Return <STATUS_SUCCESS> if no errors.
 */
int run_Bound_Viterbi_Quad(const SEQUENCE* query,     /* query sequence */
                           const HMM_PROFILE* target, /* target HMM model */
                           const int Q,               /* query length */
                           const int T,               /* target length */
                           MATRIX_3D* st_MX,          /* normal state matrix */
                           MATRIX_2D* sp_MX,          /* special state matrix */
                           EDGEBOUNDS* edg,           /* edgebounds */
                           float* sc_final);          /* OUTPUT: final score */

#endif /* BOUNDED_VITERBI_QUAD_H */
