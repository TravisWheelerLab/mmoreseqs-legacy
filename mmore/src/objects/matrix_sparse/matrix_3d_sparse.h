/*******************************************************************************
 *  FILE:      matrix_3d_sparse.h
 *  PURPOSE:   MATRIX_3D_SPARSE Float object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

#ifndef _SPARSE_MATRIX_3D_H
#define _SPARSE_MATRIX_3D_H

/*! FUNCTION:  MATRIX_3D_SPARSE_Create()
 *  SYNOPSIS:  Creates sparse matrix <smx>.
 *
 *  RETURN:    Pointer to <smx> if no errors. If errors, NULL.
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Create();

/*! FUNCTION:  MATRIX_3D_SPARSE_Destroy()
 *  SYNOPSIS:  Destroys <smx> and frees all memory.
 *
 *  RETURN:    NULL pointer.
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Destroy( MATRIX_3D_SPARSE* smx );

/*! FUNCTION:  MATRIX_3D_SPARSE_Reuse()
 *  SYNOPSIS:  Reuses <smx> by clearing previous data (no realloc).
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Reuse( MATRIX_3D_SPARSE* smx );

/*! FUNCTION:  MATRIX_3D_SPARSE_Reuse()
 *  SYNOPSIS:  Reuses <smx> by clearing previous data (no realloc).
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Reuse_Clean( MATRIX_3D_SPARSE* smx );

/*! FUNCTION:  MATRIX_3D_SPARSE_Copy()
 *  SYNOPSIS:  Copies data from <src> to <dest>.
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Copy(  MATRIX_3D_SPARSE*          mx_dest,
                        const MATRIX_3D_SPARSE*    mx_src );

/** FUNCTION:  MATRIX_3D_SPARSE_Shape_Like_Matrix()
 *  SYNOPSIS:  Shapes <src> like <dest>, but does not copy data.
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Shape_Like_Matrix(  MATRIX_3D_SPARSE*          mx_dest,
                                     const MATRIX_3D_SPARSE*    mx_src );

/*! FUNCTION:  MATRIX_3D_SPARSE_GetX_by1D()
 *  SYNOPSIS:  Gets reference to cell in <data> from matrix <smx> at flat index <i_0>.  
 *             Not designed for public function. This is generally computed by _GetX_byOffset(). 
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if out-of-bounds of sparse matrix.
 */
FLT* 
MATRIX_3D_SPARSE_GetX_by1D(   const MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                              const int                  i_0 );   /* flat index */

/*! FUNCTION:  MATRIX_3D_SPARSE_GetX_byOffset()
 *  SYNOPSIS:  Gets referenct to cell in <data> from matrix <smx> at flat index <i_0>.  
 *             Not designed for public function. This is generally computed by _GetX_byOffset(). 
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if out-of-bounds of sparse matrix.
 */
FLT* 
MATRIX_3D_SPARSE_GetX_byOffset(  const MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                                 const int                  qx0,     /* offset to start of data block */
                                 const int                  tx0,     /* offset into data block */
                                 const int                  st_0 );  /* state of data (3rd dimension) */

/*! FUNCTION:  MATRIX_3D_SPARSE_GetX()
 *  SYNOPSIS:  Get reference to cell in corresponding to given (q_0,t_0) coordinates in complete, embedding matrix.
 *             Search is linear in the number of discrete data blocks (bounds) per row. Not the fastest for matrix traversal, but good for spot lookups.
 *             Caller must call _Index() before calling this.
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if out-of-bounds of sparse matrix.
 */
FLT* 
MATRIX_3D_SPARSE_GetX(     const MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                           const int                  q_0,     /* query/row index */
                           const int                  t_0,     /* target/column index */
                           const int                  st_0 );  /* state index */

/*! FUNCTION:  EDGEBOUNDS_GetNumberBounds()
 *  SYNOPSIS:  Gets the number of <bounds> that are on row <id_0>.
 */
int
MATRIX_3D_SPARSE_GetNumberBounds(   const MATRIX_3D_SPARSE*    smx );   /* sparse matrix */

/*! FUNCTION:  EDGEBOUNDS_GetNumberBounds_byRow()
 *  SYNOPSIS:  Gets the number of <bounds> that are on row <id_0>.
 */
int
MATRIX_3D_SPARSE_GetNumberBounds_byRow(   const MATRIX_3D_SPARSE*    smx,           /* sparse matrix */
                                          const int                  q_0 );         /* IN: row/diag id in matrix <smx> */

/*! FUNCTION:  MATRIX_3D_SPARSE_GetBound_ByIndex_Fwd()
 *  SYNOPSIS:  Gets the offsets for the <r_0>th bound on entire matrix.
 *             WARNING: Doesn't do safety checks, so make sure to run _GetNumberBounds() to get maximum <r_0> first.
 *             Undefined behavior when requesting bound index <r_0> is greater than or equal to number in matrix.
 */
BOUND 
MATRIX_3D_SPARSE_GetBound_byIndex(  const MATRIX_3D_SPARSE*    smx,              /* sparse matrix */        
                                    const int                  r_0 );            /* IN: bound index in matrix <smx> */

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 */
int 
MATRIX_3D_SPARSE_GetOffset_ByIndex_Prv(   const MATRIX_3D_SPARSE*    smx,        
                                          const int                  r_0 );   /* index for given bound */

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 */
int 
MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(   const MATRIX_3D_SPARSE*    smx,        
                                          const int                  r_0 );   /* index for given bound */

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 */
int 
MATRIX_3D_SPARSE_GetOffset_ByIndex_Nxt(   const MATRIX_3D_SPARSE*    smx,        
                                          const int                  r_0 );   /* index for given bound */

/*! FUNCTION:  MATRIX_3D_SPARSE_GetOffset_byCoords()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 *             Gets the offsets for the current cell <qx0_prv_out>, 
 *             the cell immediately above it on the previous row <qx0_cur_out>, 
 *             and the cell immediately below it on the next row <qx0_nxt_out>.
 */
STATUS_FLAG 
MATRIX_3D_SPARSE_GetOffset_ByCoords(   const MATRIX_3D_SPARSE*    smx,        
                                       const int                  q_0,              /* IN: q_0 position in matrix */
                                       const int                  t_0,              /* IN: t_0 position in matrix */
                                       int*                       qx0_prv_out,      /* OUT: offset to position in previous row */
                                       int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                       int*                       qx0_nxt_out,      /* OUT: offset to position in next row */
                                       int*                       tx0_out );        /* OUT: offset to position in data block */


/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset ind <data> to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 *             Gets the offsets for the current cell <qx0_prv_out>, 
 *             the cell immediately above it on the previous row <qx0_cur_out>, 
 *             and the cell immediately below it on the next row <qx0_nxt_out>.
 */
STATUS_FLAG 
MATRIX_3D_SPARSE_GetOffset_ByIndex(    const MATRIX_3D_SPARSE*    smx,              /* sparse matrix */
                                       const int                  r_0,              /* index for given bound */
                                       int*                       qx0_prv_out,      /* OUT: offset to position in previous row in <data> */
                                       int*                       qx0_cur_out,      /* OUT: offset to position in current row in <data> */
                                       int*                       qx0_nxt_out );    /* OUT: offset to position in next row in <data> */

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 *             Gets the offsets for the current cell <qx0_prv_out>, 
 *             the cell immediately above it on the previous row <qx0_cur_out>, 
 *             and the cell immediately below it on the next row <qx0_nxt_out>.
 */
STATUS_FLAG 
MATRIX_3D_SPARSE_GetOffset_ByIndex_Fwd(   const MATRIX_3D_SPARSE*    smx,        
                                          const int                  r_0,              /* index for given bound */
                                          int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                          int*                       qx0_prv_out );    /* OUT: offset to position in previous row */

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 *             Gets the offsets for the current cell <qx0_prv_out>, 
 *             the cell immediately above it on the previous row <qx0_cur_out>, 
 *             and the cell immediately below it on the next row <qx0_nxt_out>.
 */
STATUS_FLAG 
MATRIX_3D_SPARSE_GetOffset_ByIndex_Bck(   const MATRIX_3D_SPARSE*    smx,        
                                          const int                  r_0,              /* index for given bound */
                                          int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                          int*                       qx0_nxt_out );    /* OUT: offset to position in next row */

/*! FUNCTION:  MATRIX_3D_SPARSE_GetBound_byRow_Fwd()
 *  SYNOPSIS:  Gets the offsets for the <r_0>th bound on <q_0>th row.
 *             Doesn't do safety checks, so make sure to run _GetNumberBounds_byRow to get maximum <r_0> first.
 *             Undefined behavior when requesting bound index <r_0> greater than exists on row.
 */
STATUS_FLAG 
MATRIX_3D_SPARSE_GetData_byRow_Fwd(    const MATRIX_3D_SPARSE*    smx,              /* sparse matrix */
                                       const int                  q_0,              /* IN: row/diag id in matrix <smx> */
                                       const int                  r_0,              /* IN: index of bounds of row <q_0> in matrix <smx> */
                                       BOUND*                     bnd_out,          /* OUT: pointer to bound in array */
                                       int*                       qx0_prv_out,      /* OUT: offset to position in previous row */
                                       int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                       int*                       qx0_nxt_out );    /* OUT: offset to position in next row */

/*! FUNCTION:  MATRIX_3D_SPARSE_GetBound_byRow_Bck()
 *  SYNOPSIS:  Gets the offsets for the <r_0>th bound on <q_0>th row.
 *             Doesn't do safety checks, so make sure to run _GetNumberBounds_byRow to get maximum <r_0> first.
 *             Undefined behavior when requesting bound index <r_0> greater than exists on row.
 */
STATUS_FLAG 
MATRIX_3D_SPARSE_GetData_byRow_Bck(    const MATRIX_3D_SPARSE*    smx,              /* sparse matrix */
                                       const int                  q_0,              /* IN: row/diag id in matrix <smx> */
                                       const int                  r_0,              /* IN: index of bound on row <q_0> in matrix <smx> */
                                       BOUND*                     bnd_out,          /* OUT: pointer to bound in array */
                                       int*                       qx0_prv_out,      /* OUT: offset to position in previous row */
                                       int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                       int*                       qx0_nxt_out );    /* OUT: offset to position in next row */

/*! FUNCTION:  MATRIX_3D_SPARSE_Embed()
 *  SYNOPSIS:  Embed sparse matrix <smx> into matrix <mx>. 
 *
 *  RETURN:    Pointer to <mx> if success.
 *             Returns NULL if fails.
 */
MATRIX_3D* 
MATRIX_3D_SPARSE_Embed(     int                  Q,
                            int                  T,
                            MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                            MATRIX_3D*           mx);     /* quadratic matrix */
                                    
/*! FUNCTION:  MATRIX_3D_SPARSE_Log_Embed()
 *  SYNOPSIS:  Embed sparse matrix <smx> into matrix <mx>. 
 *
 *  RETURN:    Pointer to <mx> if success.
 *             Returns NULL if fails.
 */
MATRIX_3D* 
MATRIX_3D_SPARSE_Log_Embed(     int                  Q,
                                int                  T,
                                MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                                MATRIX_3D*           mx );    /* matrix */

/*! FUNCTION:  MATRIX_3D_SPARSE_Dump()
 *  SYNOPSIS:  Dump <smx> to file pointer <fp>.
 */
void 
MATRIX_3D_SPARSE_Bounds_Dump(   MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                                FILE*                fp );    /* file pointer to be written to */


/*! FUNCTION:  MATRIX_3D_SPARSE_Bounds_Dump()
 *  SYNOPSIS:  Dump <smx> to file pointer <fp>.
 */
void 
MATRIX_3D_SPARSE_Dump(  MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                        FILE*                fp );    /* file pointer to be written to */

/*! FUNCTION:  MATRIX_3D_SPARSE_Fill()
 *  SYNOPSIS:  Fill all cells (including padding cells) in sparse matrix <smx> with <val>.
 */
void 
MATRIX_3D_SPARSE_Fill(  MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                        float                val );   /* value to fill data cells with  */

/*! FUNCTION:  MATRIX_3D_SPARSE_Fill()
 *  SYNOPSIS:  Fill all outer padding cells in sparse matrix <smx> with <val>.
 */
void 
MATRIX_3D_SPARSE_Fill_Outer(  MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                              float                val );   /* value to fill padding cells with */

/*! FUNCTION:  MATRIX_3D_SPARSE_Fill_Inner()
 *  SYNOPSIS:  Fill active inner cells in sparse matrix <smx> with <val>.
 */
void 
MATRIX_3D_SPARSE_Fill_Inner(  MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                              float                val );   /* val to fill active cells with */

/*! FUNCTION:  MATRIX_3D_SPARSE_Operation()
 *  SYNOPSIS:  Perform elementwise operation <op> to each cell in matrix.
 *             Matrices must be the same shape (not checked). 
 *             <mx_in> and <mx_out> can be the same matrix.
 */
void 
MATRIX_3D_SPARSE_Op(    MATRIX_3D_SPARSE*    mx_out,                 /* input sparse matrix */
                        MATRIX_3D_SPARSE*    mx_in,                  /* output sparse matrix */
                        FLT                  (*op)(FLT data) );      /* cell operation to be performed */ 

/*! FUNCTION:  MATRIX_3D_SPARSE_BinOp()
 *  SYNOPSIS:  Perform elementwise operation <op> to each cell in matrix.
 *             Matrices must be the same shape (not checked). 
 *             <mx_in> and <mx_out> can be the same matrix.
 */
void 
MATRIX_3D_SPARSE_BinOp(    MATRIX_3D_SPARSE*    mx_out,                                /* input sparse matrix */
                           MATRIX_3D_SPARSE*    mx_in_1,                               /* output sparse matrix */
                           MATRIX_3D_SPARSE*    mx_in_2,                               /* output sparse matrix */
                           FLT                  (*op)(FLT data_1, FLT data_2) );       /* cell operation to be performed */

/*! FUNCTION:  MATRIX_3D_SPARSE_Exp()
 *  SYNOPSIS:  Convert matrix from log-to-normal space with the exp() function.
 */
void 
MATRIX_3D_SPARSE_Exp(  MATRIX_3D_SPARSE*    smx );    /* sparse matrix */

/*! FUNCTION:  MATRIX_3D_SPARSE_Log()
 *  SYNOPSIS:  Convert matrix from normal-to-log space with the log() function.
 */
void 
MATRIX_3D_SPARSE_Log(  MATRIX_3D_SPARSE*    smx );    /* sparse matrix */

/*! FUNCTION:  MATRIX_3D_SPARSE_Add()
 *  SYNOPSIS:  Returns <smx_A> as the sum of <smx_A> and <smx_B>
 *
 *  RETURN:    Returns reference to <smx_A>.
 */
int
MATRIX_3D_SPARSE_Add(   MATRIX_3D_SPARSE*    smx_A,      /* IN: matrix addend */
                        MATRIX_3D_SPARSE*    smx_B,      /* IN: matrix addend */
                        MATRIX_3D_SPARSE*    smx_res );  /* OUT: sum matrix (can be an input) */

/*! FUNCTION:  MATRIX_3D_SPARSE_Test()
 *  SYNOPSIS:  Unit Test for MATRIX_3D_SPARSE.
 */
int 
MATRIX_3D_SPARSE_Test();


#endif /* _MATRIX_3D_SPARSE_H */