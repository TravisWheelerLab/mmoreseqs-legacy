/*******************************************************************************
 *     FILE:  bounded_fwdbck_linear.c
 *  PURPOSE:  Bounded Forward/Backward Algorithm 
 *            (Linear Space Alg)
 *
 *   AUTHOR:  Dave Rich
 *******************************************************************************/

#ifndef _BOUND_FWDBCK_SPARSE_H
#define _BOUND_FWDBCK_SPARSE_H

/** FUNCTION:  run_Bound_Forward_Sparse()
 *  SYNOPSIS:  Perform Edge-Bounded Forward step of Cloud Search Algorithm.
 *             Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of 
 *             (left-bound, right-bound) pairs sorted by row.  
 *             Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_().
 *             Final score produced by Forward is stored in <sc_final>.
 *
 *    RETURN:  Returns the final score of the Forward Algorithm.
 */
int run_Bound_Forward_Sparse(    const SEQUENCE*      query,         /* query sequence */
                                 const HMM_PROFILE*   target,        /* target HMM model */
                                 const int            Q,             /* query length */
                                 const int            T,             /* target length */
                                 MATRIX_3D_SPARSE*    st_SMX,         /* normal state matrix */
                                 MATRIX_2D*           sp_MX,         /* special state matrix */
                                 EDGEBOUNDS*          edg,           /* edgebounds */
                                 float*               sc_final );    /* (OUTPUT) final score */

/** FUNCTION:  run_Bound_Forward_Sparse()
 *  SYNOPSIS:  Perform Edge-Bounded Forward step of Cloud Search Algorithm.
 *             Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of 
 *             (left-bound, right-bound) pairs sorted by row.  
 *             Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_().
 *             Final score produced by Forward is stored in <sc_final>.
 *
 *    RETURN:  Returns the final score of the Forward Algorithm.
 */
int run_Bound_Backward_Sparse (  const SEQUENCE*      query,         /* query sequence */
                                 const HMM_PROFILE*   target,        /* target HMM model */
                                 const int            Q,             /* query length */
                                 const int            T,             /* target length */
                                 MATRIX_3D_SPARSE*    st_SMX,         /* normal state matrix */
                                 MATRIX_2D*           sp_MX,         /* special state matrix */
                                 EDGEBOUNDS*          edg,           /* edgebounds */
                                 float*               sc_final );    /* (OUTPUT) final score */

#endif /* _BOUND_FWDBCK_LINEAR_H */