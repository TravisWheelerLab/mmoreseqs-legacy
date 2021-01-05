/*******************************************************************************
 *  FILE:      matrix_3d_sparse.h
 *  PURPOSE:   MATRIX_3D_SPARSE Float object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

#ifndef _SPARSE_MATRIX_3D_H
#define _SPARSE_MATRIX_3D_H

/*  FUNCTION:  MATRIX_3D_SPARSE_Create()
 *  SYNOPSIS:  Creates sparse matrix <smx>.
 *
 *  RETURN:    Pointer to <smx> if no errors. If errors, NULL.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Create();

/*  FUNCTION:  MATRIX_3D_SPARSE_Destroy()
 *  SYNOPSIS:  Destroys <smx> and frees all memory.
 *
 *  RETURN:    NULL pointer.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Destroy( MATRIX_3D_SPARSE* smx );

/*  FUNCTION:  MATRIX_3D_SPARSE_Reuse()
 *  SYNOPSIS:  Reuses <smx> by clearing previous data (no realloc).
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Reuse( MATRIX_3D_SPARSE* smx );


/*  FUNCTION:  MATRIX_3D_SPARSE_Reuse()
 *  SYNOPSIS:  Reuses <smx> by clearing previous data (no realloc).
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Reuse_Clean( MATRIX_3D_SPARSE* smx );

/*  FUNCTION:  MATRIX_3D_SPARSE_Copy()
 *  SYNOPSIS:  Copies data from <src> to <dest>.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Copy(     MATRIX_3D_SPARSE*          mx_dest,
                                             const MATRIX_3D_SPARSE*    mx_src );

/*  FUNCTION:  MATRIX_3D_SPARSE_Shape_Like_Edgebounds()
 *  SYNOPSIS:  Creates MATRIX_3D_SPARSE object that can contain the matrix needed for 
 *             computing Bounded Forward/Backward and Bounded Viterbi algorithms. 
 *             Uses <edg_inner> as a template.
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Shape_Like_Edgebounds(  MATRIX_3D_SPARSE*    smx,              /* MATRIX_3D_SPARSE object */
                                             EDGEBOUNDS*          edg_inner );      /* EDGEBOUNDS of the inner (active) cells */


/*  FUNCTION:     EDGEBOUNDS_Create_Padded_Edgebounds()
 *  SYNOPSIS:     Create new EDGEBOUNDS <edg_outer> from given EDGEBOUNDS <edg_inner>.
 *                <edg_outer> contains all cells contained in <edg_inner> and pads with every cell adjacent to <edg_outer>, 
 *                If <edg_outer> already created, reuses data struct.
 *                Requires <edg_inner> is sorted.
 *
 *  RETURN:       Returns <edg_outer>.
 */
EDGEBOUNDS* EDGEBOUNDS_Create_Padded_Edgebounds(   EDGEBOUNDS*    edg_inner,     /* EDGEBOUNDS of the active cells */
                                                   EDGEBOUNDS*    edg_outer );   /* EDGEBOUNDS of the total cells */

/*  FUNCTION:     EDGEBOUNDS_Create_Padded_Edgebounds_Naive()
 *  SYNOPSIS:     Create new EDGEBOUNDS <edg_outer> from given EDGEBOUNDS <edg_inner>.
 *                <edg_outer> contains all cells contained in <edg_inner> and pads with every cell adjacent to <edg_outer>, 
 *                If <edg_outer> already created, reuses data struct.
 *                Requires <edg_inner> is sorted.
 *
 *  RETURN:       Returns <edg_outer>.
 */
EDGEBOUNDS* EDGEBOUNDS_Create_Padded_Edgebounds_Naive(   EDGEBOUNDS*    edg_inner,     /* EDGEBOUNDS of the active cells */
                                                         EDGEBOUNDS*    edg_outer );   /* EDGEBOUNDS of the total cells */

/*  FUNCTION:     MATRIX_3D_SPARSE_Map_to_Outer_Edgebounds()
 *  SYNOPSIS:     Maps <smx> data to <edg>.
 *                Keeps count of data ranges in <edg> in <count>.  
 *                For the start of each row in <edg>, the offset into <edg>'s bound list is stored in <smx->rows>.
 *                For each bound in <edg>, the offset into <smx> data is stored in <smx->offsets>.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Map_to_Outer_Edgebounds(   MATRIX_3D_SPARSE*    smx,        /* MATRIX_3D_SPARSE object */
                                                EDGEBOUNDS*          edg );      /* EDGEBOUNDS to map */

/*  FUNCTION:     MATRIX_3D_SPARSE_Map_to_Inner_Edgebounds()
 *  SYNOPSIS:     Maps <smx> data to <edg>.
 *                Keeps count of data ranges in <edg> in <count>.  
 *                For the start of each row in <edg>, the offset into <edg>'s bound list is stored in <smx->rows>.
 *                For each bound in <edg>, the offset into <smx> data is stored in <smx->offsets>.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Map_to_Inner_Edgebounds(   MATRIX_3D_SPARSE*    smx,           /* MATRIX_3D_SPARSE object */
                                                EDGEBOUNDS*          edg_inner,     /* inner EDGEBOUNDS to be mapped */
                                                EDGEBOUNDS*          edg_outer );   /* outer EDGEBOUNDS that describe data shape */

/*  FUNCTION:     MATRIX_3D_SPARSE_Map_to_Outer_Dump()
 *  SYNOPSIS:     Output map to screen. Shows bounds and data offsets.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Map_to_Outer_Dump(   MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                                          EDGEBOUNDS*          edg,     /* EDGEBOUNDS to map */
                                          FILE*                fp );    /* FILE pointer to output to */

/*  FUNCTION:     MATRIX_3D_SPARSE_Map_to_Inner_Dump()
 *  SYNOPSIS:     Output map to screen. Shows bounds and data offsets.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Map_to_Inner_Dump(   MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                                          EDGEBOUNDS*          edg,     /* EDGEBOUNDS to map */
                                          FILE*                fp );    /* FILE pointer to output to */

/*  FUNCTION:  MATRIX_3D_SPARSE_Add()
 *  SYNOPSIS:  Returns <smx_A> as the sum of <smx_A> and <smx_B>
 *
 *  RETURN:    Returns reference to <smx_A>.
 */
int
MATRIX_3D_SPARSE_Add(   MATRIX_3D_SPARSE*    smx_A,      /* IN: matrix addend */
                        MATRIX_3D_SPARSE*    smx_B,      /* IN: matrix addend */
                        MATRIX_3D_SPARSE*    smx_res );  /* OUT: sum matrix (can be an input) */

/*  FUNCTION:  MATRIX_3D_SPARSE_Get_X()
 *  SYNOPSIS:  Get reference to cell in corresponding to given (q_0,t_0) coordinates in complete matrix.
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if search fails.
 */
int MATRIX_3D_SPARSE_Get_X(   MATRIX_3D_SPARSE*    smx,        /* MATRIX_3D_SPARSE object */
                              int                  q_0,        /* x : row/diag id */
                              int                  t_0,        /* y : offset in row/diag */
                              int*                 off_prv,    /* */
                              int*                 off_cur,    /* */
                              int*                 off_nxt );   /* */

/*  FUNCTION:  MATRIX_3D_SPARSE_Embed()
 *  SYNOPSIS:  Embed sparse matrix <smx> into matrix <mx>. 
 *
 *  RETURN:    Pointer to <mx> if success.
 *             Returns NULL if fails.
 */
MATRIX_3D* MATRIX_3D_SPARSE_Embed(  int                  Q,
                                    int                  T,
                                    MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                                    MATRIX_3D*           mx);     /* matrix */
                                    
/*  FUNCTION:  MATRIX_3D_SPARSE_Log_Embed()
 *  SYNOPSIS:  Embed sparse matrix <smx> into matrix <mx>. 
 *
 *  RETURN:    Pointer to <mx> if success.
 *             Returns NULL if fails.
 */
MATRIX_3D* MATRIX_3D_SPARSE_Log_Embed(    int                  Q,
                                          int                  T,
                                          MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                                          MATRIX_3D*           mx );    /* matrix */

/*  FUNCTION:  MATRIX_3D_SPARSE_Dump()
 *  SYNOPSIS:  Dump <smx> to file pointer <fp>.
 */
void MATRIX_3D_SPARSE_Bounds_Dump(  MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                                    FILE*                fp );    /* file pointer to be written to */


/*  FUNCTION:  MATRIX_3D_SPARSE_Bounds_Dump()
 *  SYNOPSIS:  Dump <smx> to file pointer <fp>.
 */
void MATRIX_3D_SPARSE_Dump(   MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                              FILE*                fp );    /* file pointer to be written to */

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Get_Ref()
 *  SYNOPSIS:  Get reference to cell in corresponding to given (x,y) coordinates in complete matrix.
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if search fails.
 */
FLT* MATRIX_3D_SPARSE_Get_Ref(   MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                                 int                  q_0,     /* row/diag id */
                                 int                  t_0 );   /* offset in row/diag */

/*  FUNCTION:  MATRIX_3D_SPARSE_Fill_Outer()
 *  SYNOPSIS:  Fill all cells in sparse matrix <smx> with <val>.
 */
void MATRIX_3D_SPARSE_Fill_Outer(   MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                                    float                val );

/*  FUNCTION:  MATRIX_3D_SPARSE_Exp()
 *  SYNOPSIS:  Convert matrix from log-to-normal space with the exp() function.
 */
void MATRIX_3D_SPARSE_Exp(  MATRIX_3D_SPARSE*    smx );    /* sparse matrix */


/*  FUNCTION:  MATRIX_3D_SPARSE_Log()
 *  SYNOPSIS:  Convert matrix from normal-to-log space with the log() function.
 */
void MATRIX_3D_SPARSE_Log(  MATRIX_3D_SPARSE*    smx );    /* sparse matrix */

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