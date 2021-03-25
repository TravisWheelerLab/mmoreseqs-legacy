/*******************************************************************************
 *  FILE:      matrix_3d_sparse_build.h
 *  PURPOSE:   MATRIX_3D_SPARSE FLOAT object.
 *             Functions for building sparse matrix.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

#ifndef _MATRIX_3D_SPARSE_BUILD_H
#define _MATRIX_3D_SPARSE_BUILD_H

/*! FUNCTION:  MATRIX_3D_SPARSE_Shape_Like_Edgebounds()
 *  SYNOPSIS:  Creates MATRIX_3D_SPARSE object that can contain the matrix needed for 
 *             computing Bounded Forward/Backward and Bounded Viterbi algorithms. 
 *             Uses <edg_inner> as a template.
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int 
MATRIX_3D_SPARSE_Shape_Like_Edgebounds(     MATRIX_3D_SPARSE*    smx,              /* MATRIX_3D_SPARSE object */
                                            EDGEBOUNDS*          edg_inner );      /* EDGEBOUNDS of the inner (active) cells */


/*! FUNCTION:     EDGEBOUNDS_Create_Padded_Edgebounds()
 *  SYNOPSIS:     Create new EDGEBOUNDS <edg_outer> from given EDGEBOUNDS <edg_inner>.
 *                <edg_outer> contains all cells contained in <edg_inner> and pads with every cell adjacent to <edg_outer>, 
 *                If <edg_outer> already created, reuses data struct.
 *                Requires <edg_inner> is sorted.
 *
 *  RETURN:       Returns <edg_outer>.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Create_Padded_Edgebounds(    EDGEBOUNDS*    edg_inner,     /* EDGEBOUNDS of the active cells */
                                        EDGEBOUNDS*    edg_outer );   /* EDGEBOUNDS of the total cells */

/*! FUNCTION:     EDGEBOUNDS_Create_Padded_Edgebounds_Naive()
 *  SYNOPSIS:     Create new EDGEBOUNDS <edg_outer> from given EDGEBOUNDS <edg_inner>.
 *                <edg_outer> contains all cells contained in <edg_inner> and pads with every cell adjacent to <edg_outer>, 
 *                If <edg_outer> already created, reuses data struct.
 *                Requires <edg_inner> is sorted.
 *
 *  RETURN:       Returns <edg_outer>.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Create_Padded_Edgebounds_Naive(  EDGEBOUNDS*    edg_inner,     /* EDGEBOUNDS of the active cells */
                                            EDGEBOUNDS*    edg_outer );   /* EDGEBOUNDS of the total cells */

/*! FUNCTION:     MATRIX_3D_SPARSE_Map_to_Outer_Edgebounds()
 *  SYNOPSIS:     Maps <smx> data to <edg>.
 *                Keeps count of data ranges in <edg> in <count>.  
 *                For the start of each row in <edg>, the offset into <edg>'s bound list is stored in <smx->rows>.
 *                For each bound in <edg>, the offset into <smx> data is stored in <smx->offsets>.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int 
MATRIX_3D_SPARSE_Map_to_Outer_Edgebounds(   MATRIX_3D_SPARSE*    smx,        /* MATRIX_3D_SPARSE object */
                                            EDGEBOUNDS*          edg );      /* EDGEBOUNDS to map */

/*! FUNCTION:     MATRIX_3D_SPARSE_Map_to_Inner_Edgebounds()
 *  SYNOPSIS:     Maps <smx> data to <edg>.
 *                Keeps count of data ranges in <edg> in <count>.  
 *                For the start of each row in <edg>, the offset into <edg>'s bound list is stored in <smx->rows>.
 *                For each bound in <edg>, the offset into <smx> data is stored in <smx->offsets>.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int 
MATRIX_3D_SPARSE_Map_to_Inner_Edgebounds(   MATRIX_3D_SPARSE*    smx,           /* MATRIX_3D_SPARSE object */
                                            EDGEBOUNDS*          edg_inner,     /* inner EDGEBOUNDS to be mapped */
                                            EDGEBOUNDS*          edg_outer );   /* outer EDGEBOUNDS that describe data shape */

/*! FUNCTION:     MATRIX_3D_SPARSE_Map_to_Outer_Dump()
 *  SYNOPSIS:     Output map to screen. Shows bounds and data offsets.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int 
MATRIX_3D_SPARSE_Map_to_Outer_Dump(     MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                                        EDGEBOUNDS*          edg,     /* EDGEBOUNDS to map */
                                        FILE*                fp );    /* FILE pointer to output to */

/*! FUNCTION:     MATRIX_3D_SPARSE_Map_to_Inner_Dump()
 *  SYNOPSIS:     Output map to screen. Shows bounds and data offsets.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int 
MATRIX_3D_SPARSE_Map_to_Inner_Dump(     MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                                        EDGEBOUNDS*          edg,     /* EDGEBOUNDS to map */
                                        FILE*                fp );    /* FILE pointer to output to */


#endif /* _MATRIX_3D_SPARSE_H */