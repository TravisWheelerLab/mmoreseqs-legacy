/*******************************************************************************
 *  FILE:      ranged_fwdback_quad.h
 *  PURPOSE:   The Forward-Backward Algorithm for Sequence Alignment Search.
 *             (Quadratic Space Alg)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _RANGED_FWDBACK_QUAD_H
#define _RANGED_FWDBACK_QUAD_H

/**   FUNCTION:   run_Ranged_Forward_Quad()
 *    SYNOPSIS:   Perform Forward part of Forward-Backward Algorithm.
 *
 *      RETURN:   Return <STATUS_SUCCESS> if no errors.
 */
int run_Ranged_Forward_Quad(  const SEQUENCE*    query,        /* query sequence */
                              const HMM_PROFILE* target,       /* target hmm model */
                              const RANGE*       Q,            /* query length */
                              const RANGE*       T,            /* target length */
                              MATRIX_3D*         st_MX,        /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
                              MATRIX_2D*         sp_MX,        /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
                              float*             sc_final );   /* OUTPUT: final score */


/**   FUNCTION:   run_Ranged_Backward_Quad()
 *    SYNOPSIS:   Perform Backward part of Forward-Backward Algorithm.
 *
 *      RETURN:   Return <STATUS_SUCCESS> if no errors.
*/
int run_Ranged_Backward_Quad(    const SEQUENCE*    query,        /* query sequence */
                                 const HMM_PROFILE* target,       /* target hmm model */
                                 const RANGE*       Q,            /* query length */
                                 const RANGE*       T,            /* target length */
                                 MATRIX_3D*         st_MX,        /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
                                 MATRIX_2D*         sp_MX,        /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
                                 float*             sc_final );   /* OUTPUT: final score */

#endif /* _RANGED_FWDBACK_QUAD_H */
