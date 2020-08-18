/*******************************************************************************
 *  FILE:      matrix_3d_sparse.h
 *  PURPOSE:   MATRIX_3D_SPARSE Float object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

#ifndef _SPARSE_MATRIX_3D_H
#define _SPARSE_MATRIX_3D_H

/* 
 *  FUNCTION: 	MATRIX_3D_SPARSE_Create()
 *  SYNOPSIS: 	Creates sparse matrix <smx>.
 *
 *  RETURN:    Pointer to <smx> if no errors. If errors, NULL.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Create();

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Destroy()
 *  SYNOPSIS:  Destroys <smx> and frees all memory.
 *
 *  RETURN:    NULL pointer.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Destroy( MATRIX_3D_SPARSE* smx );


/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Reuse()
 *  SYNOPSIS:  Reuses <smx> by clearing previous data (no realloc).
 *
 *  ARGS:      <edg>     EDGEBOUNDS object.
 *
 *  RETURN:    NULL pointer.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Reuse( MATRIX_3D_SPARSE* smx );


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
EDGEBOUNDS* EDGEBOUNDS_Create_Padded_Edgebounds(  	EDGEBOUNDS*    edg_inner,
                                                   EDGEBOUNDS*    edg_outer );

/*
 *  FUNCTION:     EDGEBOUNDS_Create_Padded_Edgebounds_Naive()
 *  SYNOPSIS:     Create new EDGEBOUNDS <edg_outer> from given EDGEBOUNDS <edg_inner>.
 *                <edg_outer> contains all cells contained in <edg_inner> and pads with every cell adjacent to <edg_outer>, 
 *                If <edg_outer> already created, reuses data struct.
 *                Requires <edg_inner> is sorted.
 *
 *  ARGS:         <edg_inner>       EDGEBOUNDS of the active cells
 *                <edg_outer>       EDGEBOUNDS of the total cells
 *
 *  RETURN:       Returns <edg_outer>.
 */
EDGEBOUNDS* EDGEBOUNDS_Create_Padded_Edgebounds_Naive(   EDGEBOUNDS*    edg_inner,
                                                         EDGEBOUNDS*    edg_outer );


/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Map_to_Outer_Edgebounds()
 *  SYNOPSIS:  Maps <smx> data to <edg>.
 *             Keeps count of data ranges in <edg> in <count>.  
 *             For the start of each row in <edg>, the offset into <edg>'s bound list is stored in <smx->rows>.
 *             For each bound in <edg>, the offset into <smx> data is stored in <smx->offsets>.
 *
 *  ARGS:      <smx>          MATRIX_3D_SPARSE object
 *             <edg>          EDGEBOUNDS to map
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Map_to_Outer_Edgebounds(   MATRIX_3D_SPARSE*    smx,
                                                EDGEBOUNDS*          edg );

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Map_to_Inner_Edgebounds()
 *  SYNOPSIS:  Maps <smx> data to <edg>.
 *             Keeps count of data ranges in <edg> in <count>.  
 *             For the start of each row in <edg>, the offset into <edg>'s bound list is stored in <smx->rows>.
 *             For each bound in <edg>, the offset into <smx> data is stored in <smx->offsets>.
 *
 *  ARGS:      <smx>          MATRIX_3D_SPARSE object
 *             <edg_inner>    inner EDGEBOUNDS to be mapped
 *             <edg_outer>    outer EDGEBOUNDS that describe data shape
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Map_to_Inner_Edgebounds(   MATRIX_3D_SPARSE*    smx,
                                                EDGEBOUNDS*          edg_inner,
                                                EDGEBOUNDS*          edg_outer );

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Map_to_Outer_Dump()
 *  SYNOPSIS:  Output map to screen. Shows bounds and data offsets.
 *
 *  ARGS:      <smx>          MATRIX_3D_SPARSE object
 *             <edg>          EDGEBOUNDS to map
 *             <fp>           FILE to output to
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Map_to_Outer_Dump(   MATRIX_3D_SPARSE*    smx,
                                          EDGEBOUNDS*          edg,
                                          FILE*                fp );

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Map_to_Inner_Dump()
 *  SYNOPSIS:  Output map to screen. Shows bounds and data offsets.
 *
 *  ARGS:      <smx>          MATRIX_3D_SPARSE object
 *             <edg>          EDGEBOUNDS to map
 *             <fp>           FILE to output to
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Map_to_Inner_Dump(   MATRIX_3D_SPARSE*    smx,
                                          EDGEBOUNDS*          edg,
                                          FILE*                fp );

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Get_Ref()
 *  SYNOPSIS:  Get reference to cell in corresponding to given (x,y) coordinates in complete matrix.
 *
 *  ARGS:      <smx>          MATRIX_3D_SPARSE object
 *             <q_0>          x : row/diag id
 *             <t_0>          y : offset in row/diag
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if search fails.
 */
FLT* MATRIX_3D_SPARSE_Get_Ref(   MATRIX_3D_SPARSE*    smx,
                                 int                  q_0,
                                 int                  t_0 );

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Test()
 *  SYNOPSIS:  Unit Test for MATRIX_3D_SPARSE:
 *             Fills boolean matrix and builds an EDGEBOUND data from it.
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if search fails.
 */
int MATRIX_3D_SPARSE_Test();


#endif /* _MATRIX_3D_SPARSE_H */