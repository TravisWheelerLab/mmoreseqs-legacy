/*******************************************************************************
 *  FILE:      bounded_fwdbck_quad.c
 *  PURPOSE:   Bounded Forward-Backward Algorithm (QUADRATIC SPACE)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _BOUNDED_VITERBI_QUAD_H
#define _BOUNDED_VITERBI_QUAD_H

/*  
 *  FUNCTION: run_Bound_Viterbi_Quad()
 *  SYNOPSIS: Perform Edge-Bounded Viterbi, for recovering alignment.
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <edg>       Bounds Data (row-wise)
 *             <sc_final>  Final Score
 *
 *  RETURN: 
 */
int run_Bound_Viterbi_Quad(      const SEQUENCE*      query,         /* query sequence */
                                 const HMM_PROFILE*   target,        /* target HMM model */
                                 const int            Q,             /* query length */
                                 const int            T,             /* target length */
                                 MATRIX_3D*           st_MX,         /* normal state matrix */
                                 MATRIX_2D*           sp_MX,         /* special state matrix */
                                 EDGEBOUNDS*          edg,           /* edgebounds */
                                 float*               sc_final );    /* (OUTPUT) final score */

#endif /* BOUNDED_VITERBI_QUAD_H */