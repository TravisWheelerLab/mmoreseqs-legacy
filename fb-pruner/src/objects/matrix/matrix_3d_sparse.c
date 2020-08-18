/*******************************************************************************
 *  FILE:      matrix_3d_sparse.c
 *  PURPOSE:   MATRIX_3D_SPARSE Float object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "matrix_3d_sparse.h"

/* 
 *  FUNCTION: 	MATRIX_3D_SPARSE_Create()
 *  SYNOPSIS: 	Creates sparse matrix <smx>.
 *
 *  RETURN:    Pointer to <smx> if no errors. If errors, NULL.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Create()
{
	MATRIX_3D_SPARSE* smx = NULL;

	smx = (MATRIX_3D_SPARSE*) malloc( sizeof(MATRIX_3D_SPARSE) );
	if (smx == NULL) {
		printf("ERROR: Memory failure for allocating MATRIX_3D_SPARSE.\n");
		exit(EXIT_FAILURE);
	}

   smx->D1           = 0;
   smx->D2           = 0;
   smx->D3           = 0;

   smx->N            = 0;
   smx->Nalloc       = 0;

   smx->edg_inner    = EDGEBOUNDS_Create();
   smx->edg_outer    = EDGEBOUNDS_Create();

   smx->imap_prv     = VECTOR_INT_Create();
   smx->imap_cur     = VECTOR_INT_Create();
   smx->imap_nxt     = VECTOR_INT_Create();

   smx->omap_cur     = VECTOR_INT_Create();

   smx->data         = VECTOR_FLT_Create();
   smx->clean        = false;

	return smx;
}

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Destroy()
 *  SYNOPSIS:  Destroys <smx> and frees all memory.
 *
 *  ARGS:      <edg>     EDGEBOUNDS object.
 *
 *  RETURN:    NULL pointer.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Destroy( MATRIX_3D_SPARSE* smx )
{
   EDGEBOUNDS_Destroy( smx->edg_inner );
   EDGEBOUNDS_Destroy( smx->edg_outer );

   VECTOR_INT_Destroy( smx->imap_prv );
   VECTOR_INT_Destroy( smx->imap_cur );
   VECTOR_INT_Destroy( smx->imap_nxt );

   VECTOR_INT_Destroy( smx->omap_cur );

   VECTOR_FLT_Destroy( smx->data );

   ERRORCHECK_free(smx);
   smx = NULL;

   return smx;
}

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Reuse()
 *  SYNOPSIS:  Reuses <smx> by clearing previous data (no realloc).
 *
 *  ARGS:      <edg>     EDGEBOUNDS object.
 *
 *  RETURN:    NULL pointer.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Reuse( MATRIX_3D_SPARSE* smx )
{
   EDGEBOUNDS_Reuse( smx->edg_inner, 0, 0 );
   EDGEBOUNDS_Reuse( smx->edg_outer, 0, 0 );

   VECTOR_INT_Reuse( smx->imap_prv );
   VECTOR_INT_Reuse( smx->imap_cur );
   VECTOR_INT_Reuse( smx->imap_nxt );

   VECTOR_INT_Reuse( smx->omap_cur );

   VECTOR_FLT_Reuse( smx->data );
   smx->clean = false;

   return smx;
}


/* 
 *  FUNCTION: 	MATRIX_3D_SPARSE_Shape_Like_Edgebounds()
 *  SYNOPSIS: 	Creates MATRIX_3D_SPARSE object that can contain the matrix needed for 
 * 				computing Bounded Forward/Backward and Bounded Viterbi algorithms. 
 * 				Uses <edg_inner> as a template.
 *
 *  ARGS:      <smx>          MATRIX_3D_SPARSE object
 *             <edg_inner>    EDGEBOUNDS of the active cells
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int MATRIX_3D_SPARSE_Shape_Like_Edgebounds(  MATRIX_3D_SPARSE*    smx,
                                             EDGEBOUNDS*          edg_inner )
{
   /* verify that edgebounds are in row form */
   // if ( edg_inner->edg_mode != EDG_ROW ) {
   //    printf("ERROR: edg_inner must be in EDGE_ROW form.\n");
   //    exit(EXIT_FAILURE);
   // }

   /* get full embedding matrix dimensions */ 
   smx->D1 = edg_inner->Q + 1;
   smx->D2 = edg_inner->T + 1;
   smx->D3 = NUM_NORMAL_STATES;
   /* create edgebounds */
   smx->edg_inner = EDGEBOUNDS_Create();
   /* add inner edgebounds to sparse matrix via deep copy.  NOTE: could use reference? */
   smx->edg_inner = EDGEBOUNDS_Copy( smx->edg_inner, edg_inner );
   printf("creating padded edgebounds...\n");
	smx->edg_outer = EDGEBOUNDS_Create_Padded_Edgebounds_Naive( smx->edg_inner, smx->edg_outer );
   /* index edgebound rows */
   printf("indexing...\n");
   EDGEBOUNDS_Index( smx->edg_inner );
   EDGEBOUNDS_Index( smx->edg_outer );

   printf("INNER:\n");
   EDGEBOUNDS_Dump( smx->edg_inner, stdout );
   printf("OUTER:\n");
   EDGEBOUNDS_Dump( smx->edg_outer, stdout );

   /* map edgebounds to matrix data */
   printf("mapping outer...\n");
   MATRIX_3D_SPARSE_Map_to_Outer_Edgebounds( smx, smx->edg_outer );
   printf("mapping inner...\n");
   MATRIX_3D_SPARSE_Map_to_Inner_Edgebounds( smx, smx->edg_inner, smx->edg_outer );
   /* create matrix data */
   printf("growing...\n");
   VECTOR_FLT_GrowTo( smx->data, smx->N );
   smx->data->N = smx->N;
   /* clear data to all -INF */
   printf("filling...\n");
   VECTOR_FLT_Fill( smx->data, -INF );
   smx->clean = true;
   printf("done.\n");
}

/*
 *  FUNCTION:     EDGEBOUNDS_Create_Padded_Edgebounds()
 *  SYNOPSIS:  	Create new EDGEBOUNDS <edg_outer> from given EDGEBOUNDS <edg_inner>.
 *             	<edg_outer> contains all cells contained in <edg_inner> and pads with every cell adjacent to <edg_outer>, 
 *             	If <edg_outer> already created, reuses data struct.
 * 				   Requires <edg_inner> is sorted.
 *
 *  ARGS:         <edg_inner>       EDGEBOUNDS of the active cells
 *                <edg_outer>       EDGEBOUNDS of the total cells
 *
 *  RETURN:       Returns <edg_outer>.
 */
EDGEBOUNDS* EDGEBOUNDS_Create_Padded_Edgebounds(   EDGEBOUNDS*    edg_inner,
                                                   EDGEBOUNDS*    edg_outer )
{
   BOUND          bnd;                          /* temporary bound for adding to edgebound */
   BOUND          *bnd_1;                       /* bounds from edgebound list, for unioning */   
   int            Q, T, N;                      /* query, target, and edgebounds length */ 
   int            q_min, q_max;                 /* min and max rows in edgebounds */
   int            r_0, r_0b, r_0e;              /* (start,end] range in upper edgebound row */
   int            r_1, r_1b, r_1e;              /* (start,end] range in center edgebound row */
   int            r_2, r_2b, r_2e;              /* (start,end] range in lower edgebound row */
   int            q_0, q_1, q_2;                /* row index: prev, current, next */
   int            qx0, qx1, qx2;                /* row index (mod 3) */
   int 				n_[3];								/* counts of range arrays */
   RANGE          ranges[3][MAX_BOUNDS_PER_ROW * 3];   /* bound ranges on prev, cur, and next row */
   RANGE          rng; 
   bool           merging, is_Q; 

   /* embedded matrix dimension */
   Q = edg_inner->Q;
   T = edg_inner->T;
   /* embedded matrix length */
   N = edg_inner->N;
   /* edgebounds initialize */
   r_1b = r_1e = r_1 = 0;

   /* first and last rows in edgebounds */
   q_min = EDG_X( edg_inner, 0 ).id;
   q_max = EDG_X( edg_inner, edg_inner->N-1 ).id;

   /* set ranges to empty */
   for (int i = 0; i < 3; i++) {
   	n_[i] = 0;
   }

   /* create same size so likely no reallocation necessary */
   if ( edg_outer != NULL ) {
      edg_outer = EDGEBOUNDS_Create_by_Size( edg_inner->Nalloc );
   }
   EDGEBOUNDS_Reuse( edg_outer, edg_inner->Q, edg_inner->T );
   
   edg_outer->Q        = edg_inner->Q;
   edg_outer->T        = edg_inner->T;
   edg_outer->edg_mode = edg_inner->edg_mode;

   /* add internal bounds row-by-row  */
   for ( q_1 = q_min; q_1 <= q_max; q_1++ ) 
   {
      q_2 = q_1 - 1;
      q_0 = q_1 + 1;
      qx2 = (((q_2 % 3) + 3) % 3);   /* ensures that all mods are in the positive range */
      qx1 = (qx1 + 1) % 3;
      qx0 = (q_0 + 2) % 3;

      /* add ranges to rows */
      r_1b = r_1;

      /* get all bounds in current row */
      while ( r_1 < edg_inner->N ) 
      {
         bnd_1 = &EDG_X( edg_inner, r_1 );
         if ( bnd_1->id != q_1 ) break;

         ranges[qx2][n_[qx2]] = (RANGE){ (bnd_1->lb)-1, (bnd_1->rb)+1 };       /* add to its previous neighbors row (for insert and match state) */
         ranges[qx1][n_[qx1]] = (RANGE){ (bnd_1->lb)-1, (bnd_1->rb)+1 };       /* add to its own row (+/-1 for delete state in fwd/bck) */
         ranges[qx0][n_[qx0]] = (RANGE){ (bnd_1->lb)-1, (bnd_1->rb)+1 };       /* add to its next neighbors row (for insert and match state) */

         r_1++; 
         n_[qx0]++; n_[qx1]++; n_[qx2]++;
      }
      r_1e = r_1;    

      /* sort and merge if more than one element */
      if ( n_[qx2] > 1 )
      {
         /* sort previous range */
         for ( int i = 0; i < n_[qx2]; i++ ) 
         {
            RANGE*   cur_r    = NULL;
            RANGE*   min_r    = &ranges[qx2][i]; 
            int      min_i    = i;
            int      cmp      = 0;
            for ( int j = i; j < n_[qx2]; j++ ) 
            {
               cur_r = &ranges[qx2][j];
               cmp = RANGE_Compare( *min_r, *cur_r );
               if ( cmp > 0 ) { min_r = cur_r; }
               min_i = j;
            }
            RANGE swap           = ranges[qx2][i];
            ranges[qx2][i]       = ranges[qx2][min_i];
            ranges[qx2][min_i]   = swap;
         }

         /* merge previous range */
         for ( int i = 1; i < n_[qx2]; i++ ) 
         {
            /* if ranges overlap, merge them */
            if ( ranges[qx2][i-1].end >= ranges[qx2][i].beg ) 
            {
               ranges[qx2][i].beg = ranges[qx2][i-1].beg;
            }
            /* otherwise, add previous to edgebounds */
            else 
            {
               BOUND bnd = (BOUND) { q_2, ranges[qx2][i-1].beg, ranges[qx2][i-1].end };
               EDGEBOUNDS_Pushback( edg_outer, &bnd );
            }
         }
      }

      if ( n_[qx2] > 0 ) {
         BOUND bnd = (BOUND) { q_2, ranges[qx2][n_[qx2]-1].beg, ranges[qx2][n_[qx2]-1].end };
         EDGEBOUNDS_Pushback( edg_outer, &bnd );
      }

      n_[qx2] = 0;
   }

   q_1 = q_max;
   qx1 = q_1 % 3;
   q_0 = q_max + 1;
   qx0 = q_0 % 3;

   /* push last remaining rows */
   bnd = (BOUND) { q_1, ranges[qx1][n_[qx1]-1].beg, ranges[qx1][n_[qx1]-1].end };
   EDGEBOUNDS_Pushback( edg_outer, &bnd );
   bnd = (BOUND) { q_0, ranges[qx0][n_[qx0]-1].beg, ranges[qx0][n_[qx0]-1].end };
   EDGEBOUNDS_Pushback( edg_outer, &bnd );

   return edg_outer;
}

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
                                                         EDGEBOUNDS*    edg_outer )
{
   BOUND    bnd;
   int      id;
   int      N;

   /* create same size so likely no reallocation necessary */
   if ( edg_outer != NULL ) {
      edg_outer = EDGEBOUNDS_Create();
   }
   EDGEBOUNDS_Reuse( edg_outer, edg_inner->Q, edg_inner->T );

   N = EDGEBOUNDS_Get_Size( edg_inner );
   for ( int i = 0; i < N; i++ ) 
   {
      bnd = *EDGEBOUNDS_Get( edg_inner, i );
      bnd.lb -= 1;
      bnd.rb += 1;

      bnd.id -= 1;
      EDGEBOUNDS_Pushback( edg_outer, &bnd );
      bnd.id += 1;
      EDGEBOUNDS_Pushback( edg_outer, &bnd );
      bnd.id += 1;
      EDGEBOUNDS_Pushback( edg_outer, &bnd );
   }
   printf("sorting...\n");
   EDGEBOUNDS_Sort( edg_outer );
   EDGEBOUNDS_Dump( edg_outer, stdout );
   printf("merging...\n");
   EDGEBOUNDS_Merge( edg_outer );
   EDGEBOUNDS_Dump( edg_outer, stdout );
   exit(EXIT_SUCCESS);
}


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
                                                EDGEBOUNDS*          edg )
{
   int            N;          /* total edgebounds */
   int            offset;     /* cell count */
   int            add_offset; /* increased cell count */
   BOUND*         b_0;        /* current bound */
   VECTOR_INT*    map;        /* offset map to outer edgebounds */

   map      = smx->omap_cur;
   N        = edg->N;
   b_0      = edg->bounds;
   offset   = 0;

   VECTOR_INT_Pushback( map, offset );

   /* iterate over bounds */
   for (int i = 0; i < N; i++) 
   {
      b_0 = &(edg->bounds[i]);
      /* add edgebound cells to total offset count */
      add_offset = ( ( b_0->rb - b_0->lb ) * smx->D3 );
      offset += add_offset;
      /* map offset to bound */
      VECTOR_INT_Pushback( map, offset );
   }
   smx->N = offset;
}

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
                                                EDGEBOUNDS*          edg_outer )
{  
   int      N;             /* total edgebounds */
   int      offset;        /* cell count */
   int      offset_diff;   /* stores difference between start of outer and inner span */
   /* outer edgebounds */
   int      bo_prv_idx;    /* index of bo_prv */
   BOUND*   bo_prv;        /* previous row bound */
   int      bo_cur_idx;    /* index of bo_cur */
   BOUND*   bo_cur;        /* current row bound */
   int      bo_nxt_idx;    /* index of bo_nxt */   
   BOUND*   bo_nxt;        /* next row bound */
   /* inner edgebounds */
   int      bi_cur_idx;    /* index of bi_cur */
   BOUND*   bi_cur;        /* current row bound */
   /* row index */
   int      q_prv;         /* previous row id */
   int      q_cur;         /* current row id */
   int      q_nxt;         /* next row id */

   /* length of inner edgebound list */
   N = edg_inner->N;
   /* starting points for inner edgebound */
   bi_cur_idx = 0;
   bi_cur = &(edg_inner->bounds[0]);
   /* starting points for outer edgebounds */
   bo_prv_idx = 0;
   bo_prv = &(edg_outer->bounds[0]);
   bo_cur_idx = 0;
   bo_cur = &(edg_outer->bounds[0]);
   bo_nxt_idx = 0;
   bo_nxt = &(edg_outer->bounds[0]);

   /* iterate over edg_inner rows */
   for (int bi_cur_idx = 0; bi_cur_idx < N; bi_cur_idx++, bi_cur++) 
   {
      /* PREVIOUS */
      /* find outer edgebound span that overlaps inner edgebound span */
      /* find previous row */
      while ( (bo_prv->id == bi_cur->id - 1) == false ) {
         bo_prv++; 
         bo_prv_idx++;
      }
      /* find correct bound in row (inside range) */
      while ( (bi_cur->lb >= bo_prv->lb && bi_cur->lb < bo_prv->rb) == false ) {
         bo_prv++;
         bo_prv_idx++;
      }
      /* get offset from start of outer bound */
      offset = VECTOR_INT_Get( smx->omap_cur, bo_prv_idx );
      /* compute difference from start of outer and inner */
      offset_diff = (bi_cur->lb - bo_prv->lb) * smx->D3;
      /* add to current row inner map */
      VECTOR_INT_Pushback( smx->imap_prv, offset + offset_diff );

      /* CURRENT */
      /* find outer edgebound span that overlaps inner edgebound span */
      /* find current row */
      while ( (bo_cur->id == bi_cur->id) == false ) {
         bo_cur++; 
         bo_cur_idx++;
      }
      /* find correct bound in row (inside range) */
      while ( (bi_cur->lb >= bo_cur->lb && bi_cur->lb < bo_cur->rb) == false ) {
         bo_cur++;
         bo_cur_idx++;
      }
      /* get offset from start of outer bound */
      offset = VECTOR_INT_Get( smx->omap_cur,  bo_cur_idx );
      /* compute difference from start of outer and inner */
      offset_diff = (bi_cur->lb - bo_cur->lb) * smx->D3;
      /* add to current row inner map */
      VECTOR_INT_Pushback( smx->imap_cur, offset + offset_diff );

      /* NEXT */
      /* find outer edgebound span that overlaps inner edgebound span */
      /* find previous row */
      while ( (bo_nxt->id == bi_cur->id + 1) == false ) {
         bo_nxt++; 
         bo_nxt_idx++;
      }
      /* find correct bound in row (inside range) */
      while ( (bi_cur->lb >= bo_nxt->lb && bi_cur->lb < bo_nxt->rb) == false ) {
         bo_nxt++;
         bo_nxt_idx++;
      }
      /* get offset from start of outer bound */
      offset = VECTOR_INT_Get( smx->omap_cur, bo_nxt_idx );
      /* compute difference from start of outer and inner */
      offset_diff = (bi_cur->lb - bo_prv->lb) * smx->D3;
      /* add to current row inner map */
      VECTOR_INT_Pushback( smx->imap_nxt, offset + offset_diff );
   }
}

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
                                          FILE*                fp )
{
   int         edg_N;
   int         smx_N;
   int         offset;
   int         span;
   BOUND*      bnd;
   VECTOR_INT* map;

   edg_N = edg->N;
   smx_N = smx->data->N;
   map   = smx->omap_cur;

   for (int i = 0; i < edg_N; i++ ) 
   {
      bnd      = &EDG_X( edg, i );
      offset   = map->data[i];
      span     = (bnd->rb - bnd->lb) * smx->D3;

      fprintf(fp, "[%d/%d] {%d: %d, %d} <-> (%d-%d)/%d\n", 
         i, edg_N, bnd->id, bnd->lb, bnd->rb, offset, offset+span, smx_N );
   }
}

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
                                          FILE*                fp )
{
   int         edg_N;
   int         smx_N;
   int         offset;
   int         span;
   BOUND*      bnd;
   VECTOR_INT* map;

   edg_N = edg->N;
   smx_N = smx->data->N;

   printf("## MATRIX_3D_SPARSE\n");
   printf("## Size: %d\n", edg->N );

   printf("=> PREVIOUS:\n");
   map   = smx->imap_prv;
   for (int i = 0; i < edg_N; i++ ) 
   {
      bnd      = &EDG_X( edg, i );
      offset   = map->data[i];
      span     = (bnd->rb - bnd->lb) * smx->D3;

      fprintf(fp, "[%d/%d] {%d: %d, %d} <-> (%d-%d)/%d\n", 
         i, edg_N, bnd->id, bnd->lb, bnd->rb, offset, offset+span, smx_N );
   }

   printf("=> CURRENT:\n");
   map   = smx->imap_cur;
   for (int i = 0; i < edg_N; i++ ) 
   {
      bnd      = &EDG_X( edg, i );
      offset   = map->data[i];
      span     = (bnd->rb - bnd->lb) * smx->D3;

      fprintf(fp, "[%d/%d] {%d: %d, %d} <-> (%d-%d)/%d\n", 
         i, edg_N, bnd->id, bnd->lb, bnd->rb, offset, offset+span, smx_N );
   }

   printf("=> NEXT:\n");
   map   = smx->imap_nxt;
   for (int i = 0; i < edg_N; i++ ) 
   {
      bnd      = &EDG_X( edg, i );
      offset   = map->data[i];
      span     = (bnd->rb - bnd->lb) * smx->D3;

      fprintf(fp, "[%d/%d] {%d: %d, %d} <-> (%d-%d)/%d\n", 
         i, edg_N, bnd->id, bnd->lb, bnd->rb, offset, offset+span, smx_N );
   }
}

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Get_X()
 *  SYNOPSIS:  Get reference to cell in corresponding to given (x,y) coordinates in complete matrix.
 *             Also gets reference cells to 
 *
 *  ARGS:      <smx>          MATRIX_3D_SPARSE object
 *             <q_0>          x : row/diag id
 *             <t_0>          y : offset in row/diag
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if search fails.
 */
int MATRIX_3D_SPARSE_Get_X(   MATRIX_3D_SPARSE*    smx,        /* */
                              int                  q_0,        /* */
                              int                  t_0,        /* */
                              int*                 off_prv,    /* */
                              int*                 off_cur,    /* */
                              int*                 off_nxt )   /* */
{
   /* find bound containing (q_0,t_0) */
   int idx = EDGEBOUNDS_Search( smx->edg_inner, q_0, t_0 );
   /* if bound not found, then done */
   if (idx == -1) {
      return -1;
   }
   /* get edgebound */
   BOUND* bnd = &(smx->edg_inner->bounds[idx]);
   /* get mapped locations to start of data block */
   *off_prv = smx->imap_prv->data[idx];
   *off_cur = smx->imap_cur->data[idx];
   *off_nxt = smx->imap_nxt->data[idx];

   return idx;
}

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Embed()
 *  SYNOPSIS:  Embed 
 *
 *  ARGS:      <smx>          MATRIX_3D_SPARSE
 *             <mx>           MATRIX_3D
 *
 *  RETURN:    Pointer to <mx> if success.
 *             Returns NULL if fails.
 */
MATRIX_3D* MATRIX_3D_SPARSE_Embed(  MATRIX_3D_SPARSE*    smx,
                                    MATRIX_3D*           mx )
{
   /* TODO */
   if (mx == NULL) {
      mx = MATRIX_3D_Create( smx->D1, smx->D2, smx->D3 );
   }
   MATRIX_3D_Reuse( mx, smx->D1, smx->D2, smx->D3 );
}

/* 
 *  FUNCTION:  MATRIX_3D_SPARSE_Test()
 *  SYNOPSIS:  Unit Test for MATRIX_3D_SPARSE:
 *             Fills boolean matrix and builds an EDGEBOUND data from it.
 *
 *  RETURN:    Reference to <smx> data array location if success.
 *             Returns NULL if search fails.
 */
int MATRIX_3D_SPARSE_Test()
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