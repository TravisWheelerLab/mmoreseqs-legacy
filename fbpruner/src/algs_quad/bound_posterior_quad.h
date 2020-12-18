/*******************************************************************************
 *      FILE:   bound_posterior_quad.c
 *   PURPOSE:   Bounded Forward-Backward Algorithm (QUADRATIC SPACE)
 *
 *    AUTHOR:   Dave Rich
 *******************************************************************************/

#ifndef _BOUND_POSTERIOR_QUAD_H
#define _BOUND_POSTERIOR_QUAD_H


/** FUNCTION:   run_Bound_Posterior_Quad()
 *  SYNOPSIS:   Perform Edge-Bounded Forward part of Cloud Search Algorithm.
 *
 *    RETURN:   Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Bound_Posterior_Quad(  const SEQUENCE*      query,         /* query sequence */
                           const HMM_PROFILE*   target,        /* target HMM model */
                           const int            Q,             /* query length */
                           const int            T,             /* target length */
                           EDGEBOUNDS*          edg,           /* edgebounds */
                           MATRIX_3D*           st_MX_fwd,     /* normal state matrix, forward */
                           MATRIX_2D*           sp_MX_fwd,     /* special state matrix, forward */
                           MATRIX_3D*           st_MX_bck,     /* normal state matrix, backward */
                           MATRIX_2D*           sp_MX_bck,     /* special state matrix, backward */
                           MATRIX_3D*           st_MX_post,    /* OUTPUT: normal state matrix, posterior */
                           MATRIX_2D*           sp_MX_post );  /* OUTPUT: special state matrix, posterior */

#endif /* BOUND_POSTERIOR_QUAD_H */