/*******************************************************************************
 *  FILE:      expmax_quad.h
 *  PURPOSE:   The Baum-Welch Algorithm for Expectation Maximization.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _EXPMAX_QUAD_H
#define _EXPMAX_QUAD_H


/*  
 *  FUNCTION: run_ExpMax_Quad()
 *  SYNOPSIS: Perform the Baum-Welch algorithm for Expectation Maximization.
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <res>       Results Data
 *
 *  RETURN: 
 */
int run_ExpMax_Quad(    const SEQUENCE*      query, 
                        const HMM_PROFILE*   target, 
                        const int            Q, 
                        const int            T, 
                        MATRIX_3D*           st_MX_fwd, 
                        MATRIX_2D*           sp_MX_fwd,
                        MATRIX_3D*           st_MX_bck,
                        MATRIX_2D*           sp_MX_bck,
                        ALIGNMENT*           aln,
                        float*               sc_final );

int run_ExpMax_Forward_Quad(     const SEQUENCE*    query, 
                                 const HMM_PROFILE* target, 
                                 const int          Q, 
                                 const int          T, 
                                 MATRIX_3D*         st_MX, 
                                 MATRIX_2D*         sp_MX,
                                 float*             sc_final );

int run_ExpMax_Backward_Quad(    const SEQUENCE*    query, 
                                 const HMM_PROFILE* target, 
                                 const int          Q, 
                                 const int          T, 
                                 MATRIX_3D*         st_MX, 
                                 MATRIX_2D*         sp_MX,
                                 float*             sc_final );

int run_ExpMax_Viterbi_Quad(     const SEQUENCE*    query,
                                 const HMM_PROFILE* target,
                                 const int          Q, 
                                 const int          T, 
                                 MATRIX_3D*         st_MX,
                                 MATRIX_2D*         sp_MX,
                                 float*             sc_final );

#endif /* _EXPMAX_QUAD_H */
