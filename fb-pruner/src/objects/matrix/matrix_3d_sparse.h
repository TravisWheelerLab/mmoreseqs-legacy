/*******************************************************************************
 *  FILE:      matrix_3d_sparse.h
 *  PURPOSE:   MATRIX_3D_SPARSE Float object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

#ifndef _SPARSE_MATRIX_3D_H
#define _SPARSE_MATRIX_2D_H

/* 
 *  FUNCTION: 	MATRIX_3D_SPARSE_Create_To_Contain_Edgebounds()
 *  SYNOPSIS: 	Creates sparse matrix.
 *
 *  ARGS:      <edg>     EDGEBOUNDS object.
 *
 *  RETURN:    <MATRIX_3D_SPARSE*> if no errors. Otherwise, NULL.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Create();


/* 
 *  FUNCTION: 	MATRIX_3D_SPARSE_Shape_Like_Edgebounds()
 *  SYNOPSIS: 	Creates MATRIX_3D_SPARSE object that can contain the matrix needed for 
 * 				computing Bounded Forward/Backward and Bounded Viterbi algorithms. 
 * 				Uses <in_edg> as a template
 *
 *  ARGS:      <smx>          MATRIX_3D_SPARSE object
 *             <edg_inner>    EDGEBOUNDS of the active cells
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Shape_Like_Edgebounds(  MATRIX_3D_SPARSE*    smx,
                                             EDGEBOUNDS*          edg_inner );


/*
 *  FUNCTION:     EDGEBOUNDS_Create_Outer_Edgebounds()
 *  SYNOPSIS:  	Create new EDGEBOUNDS <edg_outer> from given EDGEBOUNDS <edg_inner>.
 *             	<edg_outer> contains all cells contained in <edg_inner> plus every cell adjacent to <edg_outer>, 
 *             	If <edg_outer> already created, reuses data struct.
 * 				   Requires <edg_inner> is sorted.
 *
 *  ARGS:         <edg_inner>       EDGEBOUNDS of the active cells
 *                <edg_outer>       EDGEBOUNDS of the total cells
 *
 *  RETURN:       Returns <edg_outer>.
 */
EDGEBOUNDS* EDGEBOUNDS_Create_Outer_Edgebounds(    EDGEBOUNDS*    edg_inner,
                                                   EDGEBOUNDS*    edg_outer );


/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Map_to_Edgebounds()
 *  SYNOPSIS:  Maps data to Edgebounds.
 *
 *  ARGS:      <smx>          MATRIX_3D_SPARSE object
 *             <edg>          EDGEBOUNDS to map
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Map_to_Edgebounds(   MATRIX_3D_SPARSE*    smx,
                                          EDGEBOUNDS*          edg );


#endif /* _MATRIX_2D_H */