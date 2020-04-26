/*******************************************************************************
 *  FILE:      merge_reorient_quad.c
 *  PURPOSE:   Functions for merging multiple EDGEBOUND objects. 
 *             Reorients from antidiagonal-wise to row-wise.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _MERGE_REORIENT_QUAD_H
#define _MERGE_REORIENT_QUAD_H

/*
 *  FUNCTION:  EDGEBOUNDS_Merge_Reorient_Quad()
 *  SYNOPSIS:  Merge two sets of EDGEBOUNDS.
 *             Reorient EDGEBOUNDS from diagonal-wise to row-wise.
 */
int EDGEBOUNDS_Merge_Reorient_Naive(const int         Q,
                                    const int         T,
                                    EDGEBOUNDS*       edg_fwd,
                                    EDGEBOUNDS*       edg_bck,
                                    EDGEBOUNDS*       edg_diag,
                                    EDGEBOUNDS*       edg_row,
                                    MATRIX_2D*        cloud_MX );

/*
 *  FUNCTION: EDGEBOUNDS_Build_From_Cloud()
 *  SYNOPSIS: Create edgebounds from a Matrix Cloud.
 */
void EDGEBOUNDS_Build_From_Cloud(const int      Q, 
                                 const int      T,
                                 EDGEBOUNDS*    edg,
                                 MATRIX_2D*     cloud_MX,
                                 int            mode );

/*
 *  FUNCTION: EDGEBOUNDS_Compare_by_Cloud()
 *  SYNOPSIS: Compare two EDGEBOUNDS by filling cloud matrices.
 */
int EDGEBOUNDS_Compare_by_Cloud( EDGEBOUNDS*    edg_a,
                                 MATRIX_2D*     mx_a,
                                 EDGEBOUNDS*    edg_b,
                                 MATRIX_2D*     mx_b );

#endif /* _MERGE_REORIENT_QUAD_H */