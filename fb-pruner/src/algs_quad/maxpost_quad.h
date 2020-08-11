/*******************************************************************************
 *  FILE:      maxpost_quad.h
 *  PURPOSE:   The Max Posterior Probability and Optimal Alignment.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _MAXPOST_QUAD_H
#define _MAXPOST_QUAD_H

/*  
 *  FUNCTION:  run_MaxPost_Quad()
 *  SYNOPSIS:  Compute the Max Posterior Probability and obtain optimal alignment.
 *             Computes the Forward and Backward.  Max Posterior computed by Viterbi.
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int run_MaxPost_Quad(   const SEQUENCE*      query,            /* query sequence */
                        const HMM_PROFILE*   target,           /* target HMM model */
                        const int            Q,                /* query length */
                        const int            T,                /* target length */
                        MATRIX_3D*           st_MX_fwd,        /* normal matrix for forward */
                        MATRIX_2D*           sp_MX_fwd,        /* special matrix for forward */
                        MATRIX_3D*           st_MX_bck,        /* normal matrix for backward */
                        MATRIX_2D*           sp_MX_bck,        /* special matrix for backward */
                        MATRIX_3D*           st_MX_post,       /* OUTPUT: normal matrix for posterior (can overwrite fwd and bck data) */
                        MATRIX_2D*           sp_MX_post,       /* OUTPUT: special matrix for posterior (can overwrite fwd and bck data) */
                        MATRIX_3D*           st_MX_max,        /* OUTPUT: normal matrix for posterior (can overwrite fwd and bck data) */
                        MATRIX_2D*           sp_MX_max,        /* OUTPUT: special matrix for posterior (can overwrite fwd and bck data) */
                        ALIGNMENT*           aln,              /* OUTPUT: alignment */
                        float*               sc_final );       /* OUTPUT: final max score */

/*
 *  FUNCTION:  run_MaxPost_Posterior_Quad()
 *  SYNOPSIS:  Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in Forward Matrices.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_MaxPost_Posterior_Quad(  const SEQUENCE*      query,            /* query sequence */
                                 const HMM_PROFILE*   target,           /* target hmm model */
                                 const int            Q,                /* query length */
                                 const int            T,                /* target length */
                                 MATRIX_3D*           st_MX_fwd,        /* normal matrix for forward */
                                 MATRIX_2D*           sp_MX_fwd,        /* special matrix for forward */
                                 MATRIX_3D*           st_MX_bck,        /* normal matrix for backward */
                                 MATRIX_2D*           sp_MX_bck,        /* special matrix for backward */
                                 MATRIX_3D*           st_MX_post,       /* OUTPUT: normal matrix for posterior */
                                 MATRIX_2D*           sp_MX_post,       /* OUTPUT: special matrix for posterior */
                                 float*               sc_final );       /* OUTPUT: final max score for posterior */

/*
 *  FUNCTION:  run_MaxPost_Viterbi_Quad()
 *  SYNOPSIS:  Run Viterbi step of Maximum Posterior Probability (no transition states).
 *             Matrices are the Maximum Poster.
 *             Computes the optimal alignment through HHM.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_MaxPost_Viterbi_Quad(    const SEQUENCE*      query,            /* query sequence */
                                 const HMM_PROFILE*   target,           /* target hmm model */
                                 const int            Q,                /* query length */
                                 const int            T,                /* target length */
                                 MATRIX_3D*           st_MX_post,       /* normal matrix */
                                 MATRIX_2D*           sp_MX_post,       /* special matrix */
                                 MATRIX_3D*           st_MX_max,        /* OUTPUT: normal matrix for max posterior */
                                 MATRIX_2D*           sp_MX_max,        /* OUTPUT: special matrix for max posterior */
                                 float*               sc_final );       /* OUTPUT: final max score */

/*
 *  FUNCTION:  run_MaxPost_Traceback_Quad()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_MaxPost_Traceback_Quad(     const SEQUENCE*     query,       /* query sequence */
                                    const HMM_PROFILE*  target,      /* HMM model */
                                    const int           Q,           /* query/seq length */
                                    const int           T,           /* target/model length */
                                    MATRIX_3D*          st_MX,       /* Normal State (Match, Insert, Delete) Matrix */
                                    MATRIX_2D*          sp_MX,       /* Special State (J,N,B,C,E) Matrix */
                                    ALIGNMENT*          aln );       /* OUTPUT: Traceback Alignment */


#endif /* _MAXPOST_QUAD_H */
