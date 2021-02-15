/*******************************************************************************
 *     FILE:   matrix_3d_sparse.c
 *  PURPOSE:   MATRIX_3D_SPARSE Float object.
 *
 *  AUTHOR:    Dave Rich
 *     BUG:
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

/* local imports */
#include "../structs.h"
#include "../../utilities/_utilities.h"
#include "../_objects.h"

/* header */
#include "_matrix_sparse.h"
#include "matrix_3d_sparse.h"

static int D3 = NUM_NORMAL_STATES;

/** FUNCTION: 	 MATRIX_3D_SPARSE_Create()
 *  SYNOPSIS: 	 Creates sparse matrix <smx>.
 *
 *    RETURN:   Pointer to <smx> if no errors. If errors, NULL.
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Create()
{
	MATRIX_3D_SPARSE* smx = NULL;

	smx = ERROR_malloc( sizeof(MATRIX_3D_SPARSE) );

   smx->D1           = 0;
   smx->D2           = 0;
   smx->D3           = 0;
   /* ranges */
   smx->Q_range      = (RANGE) {0,0};
   smx->T_range      = (RANGE) {0,0};
   smx->is_domain    = false;
   smx->D_range      = (RANGE) {0,0};
   /* data size */
   smx->N            = 0;
   smx->Nalloc       = 0;
   /* edgebounds */
   smx->edg_inner    = EDGEBOUNDS_Create();
   smx->edg_outer    = EDGEBOUNDS_Create();
   /* id index */
   smx->id_inner     = VECTOR_INT_Create();
   smx->id_outer     = VECTOR_INT_Create();
   /* inner map */
   smx->imap_prv     = VECTOR_INT_Create();
   smx->imap_cur     = VECTOR_INT_Create();
   smx->imap_nxt     = VECTOR_INT_Create();
   /* outer map */
   smx->omap_cur     = VECTOR_INT_Create();
   /* data */
   smx->data         = VECTOR_FLT_Create();
   smx->clean        = false;

   /* iterators */
   smx->q_0          = 0;
   smx->t_0          = 0;
   smx->r_0          = (RANGE) {0, 0};

	return smx;
}

/** FUNCTION:   MATRIX_3D_SPARSE_Destroy()
 *  SYNOPSIS:   Destroys <smx> and frees all memory.
 *
 *    RETURN:   NULL pointer.
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Destroy( MATRIX_3D_SPARSE* smx )
{
   /* edgebounds */
   smx->edg_inner    = EDGEBOUNDS_Destroy( smx->edg_inner );
   smx->edg_outer    = EDGEBOUNDS_Destroy( smx->edg_outer );
   /* id index */
   smx->id_inner     = VECTOR_INT_Destroy( smx->id_inner );
   smx->id_outer     = VECTOR_INT_Destroy( smx->id_outer );
   /* inner map */
   smx->imap_prv     = VECTOR_INT_Destroy( smx->imap_prv );
   smx->imap_cur     = VECTOR_INT_Destroy( smx->imap_cur );
   smx->imap_nxt     = VECTOR_INT_Destroy( smx->imap_nxt );
   /* outer map */
   smx->omap_cur     = VECTOR_INT_Destroy( smx->omap_cur );
   /* data */
   smx->data         = VECTOR_FLT_Destroy( smx->data );

   ERROR_free(smx);
   smx = NULL;

   return smx;
}

/** FUNCTION:   MATRIX_3D_SPARSE_Reuse()
 *  SYNOPSIS:   Reuses <smx> by clearing previous data (no realloc).
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Reuse( MATRIX_3D_SPARSE* smx )
{
   /* data */
   VECTOR_FLT_Reuse( smx->data );
   /* edgebounds */
   EDGEBOUNDS_Reuse( smx->edg_inner, 0, 0 );
   EDGEBOUNDS_Reuse( smx->edg_outer, 0, 0 );
   /* id index */
   VECTOR_INT_Reuse( smx->id_inner );
   VECTOR_INT_Reuse( smx->id_outer );
   /* inner map */
   VECTOR_INT_Reuse( smx->imap_prv );
   VECTOR_INT_Reuse( smx->imap_cur );
   VECTOR_INT_Reuse( smx->imap_nxt );
   /* outer map */
   VECTOR_INT_Reuse( smx->omap_cur );

   return smx;
}

/** FUNCTION:  MATRIX_3D_SPARSE_Reuse()
 *  SYNOPSIS:  Reuses <smx> by clearing previous data (no realloc).
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Reuse_Clean( MATRIX_3D_SPARSE* smx )
{
   MATRIX_3D_SPARSE_Reuse( smx );
   if ( smx->clean == false ) {
      VECTOR_FLT_Fill( smx->data, -INF );
   }
   smx->clean = true;

   return smx;
}

/** FUNCTION:  MATRIX_3D_SPARSE_Copy()
 *  SYNOPSIS:  Copies data from <src> to <dest>.
 *             If <dest> is null, new matrix is created.
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Copy(  MATRIX_3D_SPARSE*          mx_dest,
                        const MATRIX_3D_SPARSE*    mx_src )
{
   /* if dest and src point to the same matrix, do nothing */
   if ( mx_dest == mx_src ) {
      return mx_dest;
   }
   /* if dest is null, create a matrix */
   if ( mx_dest == NULL ) {
      mx_dest = MATRIX_3D_SPARSE_Create();
   }

   /* copies shape of matrix, but doesn't copy data */
   MATRIX_3D_SPARSE_Shape_Like_Matrix( mx_dest, mx_src );
   /* copy data */
   VECTOR_FLT_Copy( mx_dest->data, mx_src->data );   

   return mx_dest;
}

/** FUNCTION:  MATRIX_3D_SPARSE_Shape_Like_Matrix()
 *  SYNOPSIS:  Shapes <src> like <dest>, but does not copy data.
 */
MATRIX_3D_SPARSE* 
MATRIX_3D_SPARSE_Shape_Like_Matrix(  MATRIX_3D_SPARSE*          mx_dest,
                                     const MATRIX_3D_SPARSE*    mx_src )
{
   /* if dest and src point to the same matrix, do nothing */
   if ( mx_dest == mx_src ) {
      return mx_dest;
   }
   /* if dest has not been created, do it now */
   if ( mx_dest == NULL ) {
      mx_dest = MATRIX_3D_SPARSE_Create();
   }

   /* dimensions */
   mx_dest->D1       = mx_src->D1;
   mx_dest->D2       = mx_src->D2;
   mx_dest->D3       = mx_src->D3;
   /* min/max range */
   mx_dest->Q_range  = mx_src->Q_range;
   mx_dest->T_range  = mx_src->T_range;
   /* data */
   mx_dest->N        = mx_src->N;
   mx_dest->Nalloc   = mx_src->Nalloc;
   mx_dest->clean    = mx_src->clean;
   /* min/max */
   mx_dest->Q_range  = mx_src->Q_range;
   mx_dest->T_range  = mx_src->T_range;
   /* domain */
   mx_dest->is_domain   = mx_src->is_domain;
   mx_dest->D_range     = mx_src->D_range;
   /* edgebounds */
   EDGEBOUNDS_Copy( mx_dest->edg_inner, mx_src->edg_inner );
   EDGEBOUNDS_Copy( mx_dest->edg_outer, mx_src->edg_outer );
   /* id index */
   VECTOR_INT_Copy( mx_dest->imap_prv, mx_src->imap_prv );
   VECTOR_INT_Copy( mx_dest->imap_cur, mx_src->imap_cur );
   /* inner map */
   VECTOR_INT_Copy( mx_dest->imap_prv, mx_src->imap_prv );
   VECTOR_INT_Copy( mx_dest->imap_cur, mx_src->imap_cur );
   VECTOR_INT_Copy( mx_dest->imap_nxt, mx_src->imap_nxt );
   /* outer map */
   VECTOR_INT_Copy( mx_dest->omap_cur, mx_src->omap_cur );

   /* data (no copying, just ensure that they are the proper size to store data) */
   VECTOR_FLT_GrowTo( mx_dest->data, mx_src->data->N );

   return mx_dest;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_GetX_by1D()
 *  SYNOPSIS:  Gets reference to cell in <data> from matrix <smx> at flat index <i_0>.  
 *             Not designed for public function. This is generally computed by _GetX_byOffset(). 
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if out-of-bounds of sparse matrix.
 */
FLT* 
MATRIX_3D_SPARSE_GetX_by1D(   const MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                              const int                  i_0 )    /* flat index */
{
   return VECTOR_FLT_GetX( smx->data, i_0 );
}

/*! FUNCTION:  MATRIX_3D_SPARSE_GetX_byOffset()
 *  SYNOPSIS:  Gets reference to cell in <data> from matrix <smx> by offsets <qx0,tx0> at state <st_0>.
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if out-of-bounds of sparse matrix.
 */
FLT* 
MATRIX_3D_SPARSE_GetX_byOffset(  const MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                                 const int                  qx0,     /* offset to start of data block */
                                 const int                  tx0,     /* offset into data block */
                                 const int                  st_0 )   /* state of data (3rd dimension) */
{
   int flat_idx = qx0 + ( tx0 * NUM_NORMAL_STATES ) + ( st_0 );
   return MATRIX_3D_SPARSE_GetX_by1D( smx, flat_idx );
}

/*! FUNCTION:  MATRIX_3D_SPARSE_GetX()
 *  SYNOPSIS:  Get reference to cell in corresponding to given (q_0,t_0) coordinates in complete, embedding matrix.
 *             Search is linear in the number of discrete data blocks (bounds) per row. Not the fastest for matrix traversal, but good for spot lookups.
 *             Caller must call _Index() before calling this.
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if out-of-bounds of sparse matrix.
 */
float* 
MATRIX_3D_SPARSE_GetX(     const MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                           const int                  q_0,     /* query/row index */
                           const int                  t_0,     /* target/column index */
                           const int                  st_0 )   /* state index */
{
   /* get start of row from index */
   int idx_beg = EDGEBOUNDS_GetIndex_byRow( smx->edg_inner, q_0 );
   int idx_end = EDGEBOUNDS_GetIndex_byRow( smx->edg_inner, q_0 + 1 );
   /* linear search of column */
   for (int i_0 = idx_beg; i_0 < idx_end; i_0++ ) {
      BOUND bnd      = MATRIX_3D_SPARSE_GetBound_byIndex( smx, i_0 );
      /* find if bound in data block */
      bool in_range  = IS_IN_RANGE( bnd.lb, bnd.rb - 1, t_0 );
      if ( in_range == true ) {
         /* build flat index */
         int qx0 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur( smx, i_0 );
         int tx0 = t_0 - bnd.lb;
         return MATRIX_3D_SPARSE_GetX_byOffset( smx, qx0, tx0, st_0 );
      } 
   }

   /* if it reaches the end of the bound list, then (q_0,t_0) is not in matrix. */
   return NULL;
}

/*! FUNCTION:  EDGEBOUNDS_GetNumberBounds()
 *  SYNOPSIS:  Gets the total number of <bounds>.
 */
int
MATRIX_3D_SPARSE_GetNumberBounds(   const MATRIX_3D_SPARSE*    smx ) 
{
   int num_bounds;
   num_bounds = EDGEBOUNDS_GetSize( smx->edg_inner );
   return num_bounds;
}

/*! FUNCTION:  EDGEBOUNDS_GetNumberBounds_byRow()
 *  SYNOPSIS:  Gets the number of <bounds> that are on row <id_0>.
 */
int
MATRIX_3D_SPARSE_GetNumberBounds_byRow(   const MATRIX_3D_SPARSE*    smx,    
                                          const int                  q_0 )    /* IN: row/diag id in matrix <smx> */
{
   int num_bounds;
   num_bounds = EDGEBOUNDS_GetNumberBounds_byRow( smx->edg_inner, q_0 );
   return num_bounds;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_GetBound_ByIndex_Fwd()
 *  SYNOPSIS:  Gets the offsets for the <r_0>th bound on entire matrix.
 *             WARNING: Doesn't do safety checks, so make sure to run _GetNumberBounds() to get maximum <r_0> first.
 *             Undefined behavior when requesting bound index <r_0> is greater than or equal to number in matrix.
 */
inline
BOUND 
MATRIX_3D_SPARSE_GetBound_byIndex(  const MATRIX_3D_SPARSE*    smx,             
                                    const int                  r_0 )             /* IN: bound index in matrix <smx> */
{
   /* get bound */
   BOUND bnd;
   bnd = *EDGEBOUNDS_GetX( smx->edg_inner, r_0 );
   return bnd;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 *             Gets the offsets for the current cell <qx0_prv_out>, 
 *             the cell immediately above it on the previous row <qx0_cur_out>, 
 *             and the cell immediately below it on the next row <qx0_nxt_out>.
 */
inline
STATUS_FLAG 
MATRIX_3D_SPARSE_GetOffset_ByIndex(    const MATRIX_3D_SPARSE*    smx,        
                                       const int                  r_0,              /* index for given bound */
                                       int*                       qx0_prv_out,      /* OUT: offset to position in previous row */
                                       int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                       int*                       qx0_nxt_out )     /* OUT: offset to position in next row */
{
   /* get the offset to each row */
   int qx0_prv  = VEC_X( smx->imap_prv, r_0 );
   int qx0_cur  = VEC_X( smx->imap_cur, r_0 );
   int qx0_nxt  = VEC_X( smx->imap_nxt, r_0 );

   /* update values */
   *qx0_prv_out = qx0_prv;
   *qx0_cur_out = qx0_cur;
   *qx0_nxt_out = qx0_nxt;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 */
inline
int 
MATRIX_3D_SPARSE_GetOffset_ByIndex_Prv(   const MATRIX_3D_SPARSE*    smx,        
                                          const int                  r_0 )             /* index for given bound */
{
   int qx0_prv;
   qx0_prv  = VEC_X( smx->imap_prv, r_0 );
   return qx0_prv;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 */
inline
int 
MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(   const MATRIX_3D_SPARSE*    smx,        
                                          const int                  r_0 )             /* index for given bound */
{
   int qx0_cur;
   qx0_cur  = VEC_X( smx->imap_cur, r_0 );
   return qx0_cur;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 */
inline
int 
MATRIX_3D_SPARSE_GetOffset_ByIndex_Nxt(   const MATRIX_3D_SPARSE*    smx,        
                                          const int                  r_0 )             /* index for given bound */
{
   int qx0_nxt;
   qx0_nxt  = VEC_X( smx->imap_nxt, r_0 );
   return qx0_nxt;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 *             Gets the offsets for the current cell <qx0_prv_out>, 
 *             the cell immediately above it on the previous row <qx0_cur_out>, 
 *             and the cell immediately below it on the next row <qx0_nxt_out>.
 */
inline
STATUS_FLAG 
MATRIX_3D_SPARSE_GetOffset_ByIndex_Fwd(   const MATRIX_3D_SPARSE*    smx,        
                                          const int                  r_0,              /* index for given bound */
                                          int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                          int*                       qx0_prv_out )     /* OUT: offset to position in previous row */
{
   /* get the offset to each row */
   int qx0_prv  = VEC_X( smx->imap_prv, r_0 );
   int qx0_cur  = VEC_X( smx->imap_cur, r_0 );

   /* update values */
   *qx0_cur_out = qx0_cur;
   *qx0_prv_out = qx0_prv;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Get_CurrentRow_Offset()
 *  SYNOPSIS:  Gets the offset to start of datablock corresponding to <r_0>th bound in <edg_inner>.
 *             Gets the offsets for the current cell <qx0_prv_out>, 
 *             the cell immediately above it on the previous row <qx0_cur_out>, 
 *             and the cell immediately below it on the next row <qx0_nxt_out>.
 */
inline
STATUS_FLAG 
MATRIX_3D_SPARSE_GetOffset_ByIndex_Bck(   const MATRIX_3D_SPARSE*    smx,        
                                          const int                  r_0,              /* index for given bound */
                                          int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                          int*                       qx0_nxt_out )     /* OUT: offset to position in next row */
{
   /* get the offset to each row */
   int qx0_cur  = VEC_X( smx->imap_cur, r_0 );
   int qx0_nxt  = VEC_X( smx->imap_nxt, r_0 );

   /* update values */
   *qx0_cur_out = qx0_cur;
   *qx0_nxt_out = qx0_nxt;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_GeData_byIndex()
 *  SYNOPSIS:  Gets the offsets for the <r_0>th bound on entire matrix.
 *             WARNING: Doesn't do safety checks, so make sure to run _GetNumberBounds() to get maximum <r_0> first.
 *             Undefined behavior when requesting bound index <r_0> is greater than or equal to number in matrix.
 */
inline
STATUS_FLAG 
MATRIX_3D_SPARSE_GetData_byIndex(  const MATRIX_3D_SPARSE*     smx,             
                                    const int                  r_0,              /* IN: bound index in matrix <smx> */
                                    BOUND*                     bnd_out,          /* OUT: pointer to bound in array */
                                    int*                       qx0_prv_out,      /* OUT: offset to position in previous row */
                                    int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                    int*                       qx0_nxt_out )     /* OUT: offset to position in next row */
{
   /* get bound */
   bnd_out  = EDGEBOUNDS_GetX( smx->edg_inner, r_0 );
   /* get offsets */
   MATRIX_3D_SPARSE_GetOffset_ByIndex( smx, r_0, qx0_prv_out, qx0_cur_out, qx0_nxt_out );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_GetData_byRow_Fwd()
 *  SYNOPSIS:  Gets the offsets for the <r_0>th bound on <q_0>th row.
 *             WARNING: Doesn't do safety checks, so make sure to run _GetNumberBounds_byRow to get maximum <r_0> first.
 *             Undefined behavior when requesting bound index <r_0> is greater than or equal to number on row.
 */
inline
STATUS_FLAG 
MATRIX_3D_SPARSE_GetData_byRow_Fwd(    const MATRIX_3D_SPARSE*    smx,          
                                       const int                  q_0,              /* IN: row/diag id in matrix <smx> */
                                       const int                  r_0,              /* IN: bound index of row <q_0> in matrix <smx> */
                                       BOUND*                     bnd_out,          /* OUT: pointer to bound in array */
                                       int*                       qx0_prv_out,      /* OUT: offset to position in previous row */
                                       int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                       int*                       qx0_nxt_out )     /* OUT: offset to position in next row */
{
   int index;
   /* look up edgebound index for start of row */
   index    = EDGEBOUNDS_GetIndex_byRow_Fwd( smx->edg_inner, q_0 );
   /* get bound */
   bnd_out  = EDGEBOUNDS_GetX( smx->edg_inner, index + r_0 );
   /* get offsets */
   MATRIX_3D_SPARSE_GetOffset_ByIndex( smx, index + r_0, qx0_prv_out, qx0_cur_out, qx0_nxt_out );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_GetData_byRow_Bck()
 *  SYNOPSIS:  Gets the offsets for the <r_0>th bound on <q_0>th row.
 *             WARNING: Doesn't do safety checks, so make sure to run _GetNumberBounds_byRow to get maximum <r_0> first.
 *             Undefined behavior when requesting bound index <r_0> is greater than or equal to number on row.
 */
inline
STATUS_FLAG 
MATRIX_3D_SPARSE_GetData_byRow_Bck(    const MATRIX_3D_SPARSE*    smx,        
                                       const int                  q_0,              /* IN: row/diag id in matrix <smx> */
                                       const int                  r_0,              /* IN: bound index of row <q_0> in matrix <smx> */
                                       BOUND*                     bnd_out,          /* OUT: pointer to bound in array */
                                       int*                       qx0_prv_out,      /* OUT: offset to position in previous row */
                                       int*                       qx0_cur_out,      /* OUT: offset to position in current row */
                                       int*                       qx0_nxt_out )     /* OUT: offset to position in next row */
{
   int      index;
   /* look up edgebound index for start of row */
   index    = EDGEBOUNDS_GetIndex_byRow_Bck( smx->edg_inner, q_0 );
   /* get bound */
   bnd_out  = EDGEBOUNDS_GetX( smx->edg_inner, index + r_0 );
   /* get offsets */
   MATRIX_3D_SPARSE_GetOffset_ByIndex( smx, index + r_0, qx0_prv_out, qx0_cur_out, qx0_nxt_out );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Embed()
 *  SYNOPSIS:  Embed sparse matrix <smx> into matrix <mx>. 
 *
 *  RETURN:    Pointer to <mx> if success.
 *             Returns NULL if fails.
 */
MATRIX_3D* 
MATRIX_3D_SPARSE_Embed(    int                  Q,
                           int                  T,
                           MATRIX_3D_SPARSE*    smx,     /* sparse matrix to embed */
                           MATRIX_3D*           mx )     /* matrix to be embedded into */
{
   int         id_0, lb_0, rb_0;
   int         q_0, t_0;
   int         qx0, tx0;
   int         r_0b, r_0e, r_0;
   EDGEBOUNDS* edg;
   BOUND*      bnd;

   /* resize embedding matrix to contain sparse matrix */
   MATRIX_3D_Reuse( mx, NUM_NORMAL_STATES, Q+1, T+1 );
   MATRIX_3D_Fill( mx, -INF );
 
   /* edgebounds */ 
   edg  = smx->edg_inner;
   r_0b = 0;
   r_0e = EDGEBOUNDS_GetSize( edg );

   /* iterate through edgebounds */
   for ( int r_0 = r_0b; r_0 < r_0e; r_0++ ) 
   {
      /* get bound */
      bnd   = EDGEBOUNDS_GetX( edg, r_0 );
      q_0   = bnd->id;
      lb_0  = bnd->lb;
      rb_0  = bnd->rb;

      /* fetch data mapping bound start location to data block in sparse matrix */
      qx0 = VECTOR_INT_Get( smx->imap_cur, r_0 );    /* (q_0, t_0) location offset */

      for (t_0 = lb_0; t_0 < rb_0; t_0++)
      {
         tx0 = t_0 - lb_0;

         /* embed linear row into quadratic test matrix */
         float val = SMX_X(smx, MAT_ST, qx0, tx0);
         MX_3D(mx, MAT_ST, q_0, t_0) = val;
         MX_3D(mx, INS_ST, q_0, t_0) = SMX_X(smx, INS_ST, qx0, tx0);
         MX_3D(mx, DEL_ST, q_0, t_0) = SMX_X(smx, DEL_ST, qx0, tx0);
      }
   }
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Log_Embed()
 *  SYNOPSIS:  Embed sparse matrix <smx> into matrix <mx> in log scale.
 *             Does not effect sparse matrix. 
 *
 *  RETURN:    Pointer to <mx> if success.
 *             Returns NULL if fails.
 */
MATRIX_3D* 
MATRIX_3D_SPARSE_Log_Embed(   int                  Q,
                              int                  T,
                              MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                              MATRIX_3D*           mx )     /* matrix */
{
   int         id_0, lb_0, rb_0;
   int         q_0, t_0;
   int         qx0, tx0;
   int         r_0b, r_0e, r_0;
   EDGEBOUNDS* edg;
   BOUND*      bnd;

   /* resize embedding matrix to contain sparse matrix */
   MATRIX_3D_Reuse( mx, NUM_NORMAL_STATES, Q+1, T+1 );
   MATRIX_3D_Fill( mx, -INF );
 
   /* edgebounds */ 
   edg  = smx->edg_inner;
   r_0b = 0;
   r_0e = EDGEBOUNDS_GetSize( edg );

   /* iterate through edgebounds */
   for ( int r_0 = r_0b; r_0 < r_0e; r_0++ ) 
   {
      /* get bound */
      bnd   = EDGEBOUNDS_GetX(edg, r_0);
      q_0   = id_0;

      /* fetch data mapping bound start location to data block in sparse matrix */
      qx0 = VECTOR_INT_Get( smx->imap_cur, r_0 );    /* (q_0, t_0) location offset */

      for (t_0 = bnd->lb; t_0 < bnd->rb; t_0++)
      {
         tx0 = t_0 - bnd->lb;

         /* embed linear row into quadratic test matrix */
         MX_3D(mx, MAT_ST, q_0, t_0) = logf( SMX_X(smx, MAT_ST, qx0, tx0));
         MX_3D(mx, INS_ST, q_0, t_0) = logf( SMX_X(smx, INS_ST, qx0, tx0));
         MX_3D(mx, DEL_ST, q_0, t_0) = logf( SMX_X(smx, DEL_ST, qx0, tx0));
      }
   }
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Bounds_Dump()
 *  SYNOPSIS:  Dump <smx> bound data to file pointer <fp>.
 */
void 
MATRIX_3D_SPARSE_Bounds_Dump(    MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                                 FILE*                fp )     /* file pointer to be written to */
{
   int   Q, T;
   int   prv, cur, nxt, len;
   int   width = 11;
   int   hwidth = 4;
   int   swidth = 3;

   Q = smx->edg_inner->Q;
   T = smx->edg_inner->T;

   fprintf( fp, "%*s :: %*s :: %*s :: %*s ::\n",
      14, "BOUNDS", width, "PREVIOUS", width, "CURRENT", width, "NEXT");
   for ( int i = 0; i < smx->imap_cur->N; i++ ) 
   {
      BOUND* bnd = EDGEBOUNDS_GetX( smx->edg_inner, i );
      prv = smx->imap_prv->data[i];
      cur = smx->imap_cur->data[i];
      nxt = smx->imap_nxt->data[i];
      // len = bnd->rb - bnd->lb;
      len = 0;
      fprintf( fp, "[%*d]{%*d-%*d} :: %*d - %*d :: %*d - %*d :: %*d - %*d ::\n",
         swidth, bnd->id, swidth, bnd->lb, swidth, bnd->rb,
         hwidth, prv, hwidth, prv + (len * 3),
         hwidth, cur, hwidth, cur + (len * 3),
         hwidth, nxt, hwidth, nxt + (len * 3) );
   }
   fprintf(fp, "\n\n");

   BOUND*   bnd_0;
   int      q_0, t_0;
   int      lb_0, rb_0;
   int      r_0b, r_0e, r_0;
   int      cnt;
   int      pad;

   cnt = 0;
   for (int r_0 = 0; r_0 < EDGEBOUNDS_GetSize( smx->edg_outer ); r_0++)
   {
      bnd_0 = EDGEBOUNDS_GetX( smx->edg_outer, r_0 );
      fprintf(fp, "[%3d]{%3d: %3d, %3d}  ", r_0, bnd_0->id, bnd_0->lb, bnd_0->rb);

      for (int i = -1; i < Q; i++)
      {
         if (i >= bnd_0->lb && i < bnd_0->rb) {
            fprintf(fp, "%*.4f\t", pad, MSMX_X(smx, 0, cnt) );
            cnt++;
         }
         else {
            fprintf(fp, "%*s\t", pad, "****" );
         }
      }
      fprintf(fp, "\n");
   }
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Bounds_Dump()
 *  SYNOPSIS:  Dump <smx> to file pointer <fp>.
 */
void 
MATRIX_3D_SPARSE_Dump(     MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                           FILE*                fp )     /* file pointer to be written to */
{
   int      Q, T;
   BOUND*   bnd_0;
   int      q_0, t_0;
   int      lb_0, rb_0;
   int      r_0b, r_0e, r_0;
   int      cnt;
   int      pad;

   Q = smx->edg_inner->Q;
   T = smx->edg_inner->T;

   cnt = 0;
   for (int r_0 = 0; r_0 < EDGEBOUNDS_GetSize( smx->edg_outer ); r_0++)
   {
      bnd_0 = EDGEBOUNDS_GetX( smx->edg_outer, r_0 );
      fprintf(fp, "[%3d]{%3d: %3d, %3d}  ", r_0, bnd_0->id, bnd_0->lb, bnd_0->rb);

      for (int i = -1; i < T+1; i++)
      {
         if (i >= bnd_0->lb && i < bnd_0->rb) {
            fprintf(fp, "%*.4f\t", pad, MSMX_X(smx, 0, cnt) );
            cnt++;
         }
         else {
            fprintf(fp, "%*s\t", pad, "****" );
         }
      }
      fprintf(fp, "\n");
   }
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Fill()
 *  SYNOPSIS:  Fill all cells in sparse matrix <smx> with <val>.
 */
void 
MATRIX_3D_SPARSE_Fill(  MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                        float                val )    /* value to fill data cells with */
{
   VECTOR_FLT_Fill( smx->data, val );
}

/** TODO: Implement this */ 
/*! FUNCTION:  MATRIX_3D_SPARSE_Fill()
 *  SYNOPSIS:  Fill all outer padding cells in sparse matrix <smx> with <val>.
 */
void 
MATRIX_3D_SPARSE_Fill_Outer(  MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                              float                val )    /* value to fill padding cells with */
{
   VECTOR_FLT_Fill( smx->data, val );
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Fill_Inner()
 *  SYNOPSIS:  Fill active inner cells in sparse matrix <smx> with <val>.
 */
void 
MATRIX_3D_SPARSE_Fill_Inner(  MATRIX_3D_SPARSE*    smx,     /* sparse matrix */
                              float                val )    /* val to fill active cells with */
{
   VECTOR_FLT_Fill( smx->data, val );
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Operation()
 *  SYNOPSIS:  Perform elementwise operation <op> to each cell in matrix.
 *             Matrices must be the same shape (not checked). 
 *             <mx_in> and <mx_out> can be the same matrix.
 */
void 
MATRIX_3D_SPARSE_Op(    MATRIX_3D_SPARSE*    mx_out,                 /* input sparse matrix */
                        MATRIX_3D_SPARSE*    mx_in,                  /* output sparse matrix */
                        FLT                  (*op)(FLT data) )       /* cell operation to be performed */ 
{
   /* TODO: Maybe add a quick check if matrices are the same type? */


   if ( mx_out == NULL ) {
      MATRIX_3D_SPARSE_Copy( mx_out, mx_in );
   }
   mx_out->data = VECTOR_FLT_Op( mx_in->data, mx_out->data, op );
}

/*! FUNCTION:  MATRIX_3D_SPARSE_BinOp()
 *  SYNOPSIS:  Perform elementwise operation <op> to each cell in matrix.
 *             Matrices must be the same shape (not checked). 
 *             <mx_in> and <mx_out> can be the same matrix.
 */
void 
MATRIX_3D_SPARSE_BinOp(    MATRIX_3D_SPARSE*    mx_out,                                /* input sparse matrix */
                           MATRIX_3D_SPARSE*    mx_in_1,                               /* output sparse matrix */
                           MATRIX_3D_SPARSE*    mx_in_2,                               /* output sparse matrix */
                           FLT                  (*op)(FLT data_1, FLT data_2) )        /* cell operation to be performed */ 
{
   if ( mx_out == NULL ) {
      MATRIX_3D_SPARSE_Shape_Like_Matrix( mx_out, mx_in_1 );
   }
   mx_out->data = VECTOR_FLT_BinOp( mx_out->data, mx_in_1->data, mx_in_1->data, op );
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Exp()
 *  SYNOPSIS:  Convert matrix from log-to-normal space with the exp() function.
 */
void 
MATRIX_3D_SPARSE_Exp(  MATRIX_3D_SPARSE*    smx )     /* sparse matrix */
{
   VECTOR_FLT_Op( smx->data, smx->data, MATH_Exp );
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Log()
 *  SYNOPSIS:  Convert matrix from normal-to-log space with the log() function.
 */
void 
MATRIX_3D_SPARSE_Log(  MATRIX_3D_SPARSE*    smx )     /* sparse matrix */
{
   VECTOR_FLT_Op( smx->data, smx->data, MATH_Log );
}

/*! FUNCTION:  MATRIX_3D_SPARSE_Iter_Reset()
 *  SYNOPSIS:  Reset iterators for matrix.
 */
void 
MATRIX_3D_SPARSE_Iter_Reset(  MATRIX_3D_SPARSE*    smx )     /* sparse matrix */
{
   smx->r_0 = (RANGE) {0, 0};
   smx->q_0 = 0;
   smx->t_0 = 0;
}


/*! FUNCTION:  MATRIX_3D_SPARSE_Add()
 *  SYNOPSIS:  Returns <smx_A> as the sum of <smx_A> and <smx_B>
 *
 *  RETURN:    Returns reference to <smx_A>
 */
int 
MATRIX_3D_SPARSE_Add(   MATRIX_3D_SPARSE*    smx_A,      /* IN: addent matrix */
                        MATRIX_3D_SPARSE*    smx_B,      /* IN: addend matrix */
                        MATRIX_3D_SPARSE*    smx_res )   /* OUT: sum matrix (can be an input) */
{
   /* if 

   /* sum contents of matrix */
   for ( int i = 0; i < smx_A->data->N; i++ ) 
   {
      if ( smx_A->data->data[i] == -INF || smx_B->data->data[i] == -INF )  {
         smx_res->data->data[i] = -INF;
      }
      else {
         smx_res->data->data[i] = smx_A->data->data[i] + smx_B->data->data[i];
      }
      
   }
   return STATUS_SUCCESS;
}


/*! FUNCTION:  MATRIX_3D_SPARSE_Test()
 *  SYNOPSIS:  Unit Test for MATRIX_3D_SPARSE:
 *             Fills boolean matrix and builds an EDGEBOUND data from it.
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if search fails.
 */
int 
MATRIX_3D_SPARSE_Test()
{
   /* matrix dim */
   int Q = 10;
   int T = 20;

   /* transition states */
   float trans[2][2];
   trans[0][0] = 0.80;   /* 0 -> 0 */
   trans[0][1] = 0.20;   /* 0 -> 1 */
   trans[1][0] = 0.25;   /* 1 -> 0 */
   trans[1][1] = 0.75;   /* 1 -> 1 */
}