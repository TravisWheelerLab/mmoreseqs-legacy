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
#include "matrix_3d_sparse_build.h"

/** FUNCTION:  MATRIX_3D_SPARSE_Shape_Like_Edgebounds()
 *  SYNOPSIS:  Creates MATRIX_3D_SPARSE object that can contain the matrix needed for 
 *             computing Bounded Forward/Backward and Bounded Viterbi algorithms. 
 *             Uses <edg_inner> as a template.
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int 
MATRIX_3D_SPARSE_Shape_Like_Edgebounds(   MATRIX_3D_SPARSE*    smx,              /* MATRIX_3D_SPARSE object */
                                          EDGEBOUNDS*          edg_inner )       /* EDGEBOUNDS of the inner (active) cells */
{
   /* get full embedding matrix dimensions */ 
   smx->D1 = edg_inner->Q + 1;
   smx->D2 = edg_inner->T + 1;
   smx->D3 = NUM_NORMAL_STATES;

   /* add inner edgebounds to sparse matrix via deep copy.  NOTE: could use reference? */
   smx->edg_inner = EDGEBOUNDS_Copy( smx->edg_inner, edg_inner );
	
   /* create outer edgebounds that includes all cells adjacent to inner edgebounds */
   smx->edg_outer = EDGEBOUNDS_Create_Padded_Edgebounds( smx->edg_inner, smx->edg_outer );

   /* map edgebounds to matrix data */
   MATRIX_3D_SPARSE_Map_to_Outer_Edgebounds( smx, smx->edg_outer );
   MATRIX_3D_SPARSE_Map_to_Inner_Edgebounds( smx, smx->edg_inner, smx->edg_outer );

   /* find minimum and maximum Q and T values */
   // EDGEBOUNDS_Find_BoundingBox( smx->edg_inner, &smx->Q_range, &smx->T_range );
   // EDGEBOUNDS_Find_BoundingBox( smx->edg_outer, NULL, NULL );

   /* edgebounds are self-indexing, they just need to be called */
   EDGEBOUNDS_Index( smx->edg_inner );
   EDGEBOUNDS_Index( smx->edg_outer );

   /* create matrix data */
   int old_N = smx->data->N;
   VECTOR_FLT_GrowTo( smx->data, smx->N );
   smx->data->N = smx->N;

   /* if current data is clean, only clean new data */
   VECTOR_FLT_Fill( smx->data, -INF );
   if ( smx->clean = false ) {

   } else {

   }
   
   smx->edg_inner->edg_mode = EDG_ROW;
   smx->edg_outer->edg_mode = EDG_ROW;
   smx->clean = true;
}

/** FUNCTION:   EDGEBOUNDS_Create_Padded_Edgebounds()
 *  SYNOPSIS:   Create new EDGEBOUNDS <edg_outer> from given EDGEBOUNDS <edg_inner>.
 *              Method selector.
 *
 *    RETURN:   Returns <edg_outer>.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Create_Padded_Edgebounds(   EDGEBOUNDS*    edg_inner,     /* EDGEBOUNDS of the active cells */
                                       EDGEBOUNDS*    edg_outer )    /* EDGEBOUNDS of the total cells */
{
   edg_outer = EDGEBOUNDS_Create_Padded_Edgebounds_Naive( edg_inner, edg_outer );
   return edg_outer;
}

/* TODO: WIP */
/** FUNCTION:   EDGEBOUNDS_Create_Padded_Edgebounds()
 *  SYNOPSIS:   Create new EDGEBOUNDS <edg_outer> from given EDGEBOUNDS <edg_inner>.
 *              <edg_outer> contains all cells contained in <edg_inner> and pads with every cell adjacent to <edg_outer>, 
 *              If <edg_outer> already created, reuses data struct.
 *              Requires <edg_inner> to be sorted.
 *
 *    RETURN:   Returns <edg_outer>.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Create_Padded_Edgebounds_Optimal(    EDGEBOUNDS*    edg_inner,     /* EDGEBOUNDS of the active cells */
                                                EDGEBOUNDS*    edg_outer )    /* EDGEBOUNDS of the total cells */
{
   BOUND          bnd;                          /* temporary bound for adding to edgebound */
   BOUND*         bnd_1;                        /* bounds from edgebound list, for unioning */   
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
   N = EDGEBOUNDS_GetSize( edg_inner );
   /* edgebounds initialize */
   r_1b = r_1e = r_1 = 0;

   /* first and last rows in edgebounds */
   q_min = EDG_X( edg_inner, 0 ).id;
   q_max = EDG_X( edg_inner, N-1 ).id;

   /* set ranges to empty */
   for (int i = 0; i < 3; i++) {
   	n_[i] = 0;
   }

   /* create same size so likely no reallocation necessary */
   if ( edg_outer != NULL ) {
      edg_outer = EDGEBOUNDS_Create_by_Size( N );
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
      while ( r_1 < EDGEBOUNDS_GetSize( edg_inner ) ) 
      {
         bnd_1 = EDGEBOUNDS_GetX( edg_inner, r_1 );
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
               EDGEBOUNDS_Pushback( edg_outer, bnd );
            }
         }
      }

      if ( n_[qx2] > 0 ) {
         BOUND bnd = (BOUND) { q_2, ranges[qx2][n_[qx2]-1].beg, ranges[qx2][n_[qx2]-1].end };
         EDGEBOUNDS_Pushback( edg_outer, bnd );
      }

      n_[qx2] = 0;
   }

   q_1 = q_max;
   qx1 = q_1 % 3;
   q_0 = q_max + 1;
   qx0 = q_0 % 3;

   /* push last remaining rows */
   bnd = (BOUND) { q_1, ranges[qx1][n_[qx1]-1].beg, ranges[qx1][n_[qx1]-1].end };
   EDGEBOUNDS_Pushback( edg_outer, bnd );
   bnd = (BOUND) { q_0, ranges[qx0][n_[qx0]-1].beg, ranges[qx0][n_[qx0]-1].end };
   EDGEBOUNDS_Pushback( edg_outer, bnd );

   return edg_outer;
}

/*! FUNCTION:     EDGEBOUNDS_Create_Padded_Edgebounds_Naive()
 *  SYNOPSIS:     Create new EDGEBOUNDS <edg_outer> from given EDGEBOUNDS <edg_inner>.
 *                <edg_outer> contains all cells contained in <edg_inner> and pads with every cell adjacent to <edg_outer>, 
 *                If <edg_outer> already created, reuses data struct.
 *                Requires <edg_inner> is sorted.
 *                Requires 3x as much memory as non-naive version. In practice, this is just ~O(Q).
 *
 *    RETURN:     Returns <edg_outer>.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Create_Padded_Edgebounds_Naive(   EDGEBOUNDS*    edg_inner,     /* EDGEBOUNDS of the active cells */
                                             EDGEBOUNDS*    edg_outer )    /* EDGEBOUNDS of the total cells */
{
   FILE*    fp;
   BOUND    bnd;
   int      id;
   int      N;

   /* create same size so likely no reallocation necessary */
   if ( edg_outer == NULL ) {
      edg_outer = EDGEBOUNDS_Create();
   }
   /* empty edgebound list */
   EDGEBOUNDS_Reuse( edg_outer, edg_inner->Q, edg_inner->T );

   /* fill metadata */
   edg_outer->Q = edg_inner->Q;
   edg_outer->T = edg_inner->T;
   edg_outer->edg_mode = edg_inner->edg_mode;

   /* iterate over each bound in inner edgebound list */
   N = EDGEBOUNDS_GetSize( edg_inner );
   for ( int i = 0; i < N; i++ ) 
   {
      bnd = *EDGEBOUNDS_GetX( edg_inner, i );
      /* add padding to both sides of range */
      bnd.lb -= 1;
      bnd.rb += 1;

      /* add upper row */
      bnd.id -= 1;
      EDGEBOUNDS_Pushback( edg_outer, bnd );
      /* add center row */
      bnd.id += 1;
      EDGEBOUNDS_Pushback( edg_outer, bnd );
      /* add lower row */
      bnd.id += 1;
      EDGEBOUNDS_Pushback( edg_outer, bnd );
   }

   EDGEBOUNDS_Sort( edg_outer );

   #if DEBUG
   {
      fp = fopen("test_output/my.outer_premerge.edg", "w");
      EDGEBOUNDS_Dump( edg_outer, fp );
      fclose(fp);
   }
   #endif
   
   EDGEBOUNDS_Merge( edg_outer );

   return edg_outer;
}

/*! FUNCTION:     MATRIX_3D_SPARSE_Map_to_Outer_Edgebounds()
 *  SYNOPSIS:     Maps <smx> data to <edg>.
 *                Keeps count of data ranges in <edg> in <count>.  
 *                For the start of each row in <edg>, the offset into <edg>'s bound list is stored in <smx->rows>.
 *                For each bound in <edg>, the offset into <smx> data is stored in <smx->offsets>.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int 
MATRIX_3D_SPARSE_Map_to_Outer_Edgebounds(    MATRIX_3D_SPARSE*    smx,        /* MATRIX_3D_SPARSE object */
                                             EDGEBOUNDS*          edg )       /* EDGEBOUNDS to map */
{
   int            N;             /* total edgebounds */
   int            offset;        /* cell count */
   int            add_offset;    /* increased cell count */
   BOUND*         b_0;           /* current bound */
   VECTOR_INT*    map;           /* offset map to outer edgebounds */

   map      = smx->omap_cur;
   N        = EDGEBOUNDS_GetSize( edg );
   b_0      = EDGEBOUNDS_GetX( edg, 0 );
   offset   = 0;

   /* clear previous data from vector */
   VECTOR_INT_Reuse( map );

   /* initial row */
   VECTOR_INT_Pushback( map, offset );

   /* iterate over bounds */
   for (int i = 0; i < N; i++) 
   {
      b_0 = EDGEBOUNDS_GetX( edg, i );
      /* add edgebound cells to total offset count */
      add_offset = ( ( b_0->rb - b_0->lb ) * smx->D3 );
      offset += add_offset;
      /* map offset to bound */
      VECTOR_INT_Pushback( map, offset );
   }

   /* total cells in outer edgebounds */
   smx->N = offset;
}

/*! FUNCTION:     MATRIX_3D_SPARSE_Map_to_Inner_Edgebounds()
 *  SYNOPSIS:     Maps <smx> data to <edg>.
 *                Keeps count of data ranges in <edg> in <count>.  
 *                For the start of each row in <edg>, the offset into <edg>'s bound list is stored in <smx->rows>.
 *                For each bound in <edg>, the offset into <smx> data is stored in <smx->offsets>.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int 
MATRIX_3D_SPARSE_Map_to_Inner_Edgebounds(    MATRIX_3D_SPARSE*    smx,           /* MATRIX_3D_SPARSE object */
                                             EDGEBOUNDS*          edg_inner,     /* inner EDGEBOUNDS to be mapped */
                                             EDGEBOUNDS*          edg_outer )    /* outer EDGEBOUNDS that describe data shape */
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
   N = EDGEBOUNDS_GetSize( edg_inner );
   /* starting points for inner edgebound */
   bi_cur_idx = 0;
   bi_cur = EDGEBOUNDS_GetX( edg_inner, 0 );
   /* starting points for outer edgebounds */
   bo_prv_idx = 0;
   bo_prv = EDGEBOUNDS_GetX( edg_outer, 0 );
   bo_cur_idx = 0;
   bo_cur = EDGEBOUNDS_GetX( edg_outer, 0 );
   bo_nxt_idx = 0;
   bo_nxt = EDGEBOUNDS_GetX( edg_outer, 0 );

   /* clear previous data */
   VECTOR_INT_Reuse( smx->imap_prv );
   VECTOR_INT_Reuse( smx->imap_cur );
   VECTOR_INT_Reuse( smx->imap_nxt );

   /* iterate over edg_inner rows */
   for (int bi_cur_idx = 0; bi_cur_idx < N; bi_cur_idx++, bi_cur++) 
   {
      /* PREVIOUS */
      /* find outer edgebound span that overlaps inner edgebound span */
      /* find outer bound on same row as inner bound */
      while ( (bo_prv->id == bi_cur->id - 1) == false ) {
         bo_prv++; 
         bo_prv_idx++;
      }
      /* find outer bound on same column range as inner bound */
      while ( (bi_cur->lb >= bo_prv->lb && bi_cur->lb < bo_prv->rb) == false ) {
         bo_prv++;
         bo_prv_idx++;
      }
      /* get offset from start of outer bound */
      offset = VECTOR_INT_Get( smx->omap_cur, bo_prv_idx );
      /* compute difference from start of outer and inner bounds */
      offset_diff = (bi_cur->lb - bo_prv->lb) * smx->D3;
      /* add to previous row inner map */
      VECTOR_INT_Pushback( smx->imap_prv, offset + offset_diff );

      /* CURRENT */
      /* find outer edgebound span that overlaps inner edgebound span */
      /* find outer bound on same row as inner bound */
      while ( (bo_cur->id == bi_cur->id) == false ) {
         bo_cur++; 
         bo_cur_idx++;
      }
      /* find outer bound on same column range as inner bound */
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
      /* find outer bound on same row as inner bound */
      while ( (bo_nxt->id == bi_cur->id + 1) == false ) {
         bo_nxt++; 
         bo_nxt_idx++;
      }
      /* find outer bound on same column range as inner bound */
      while ( (bi_cur->lb >= bo_nxt->lb && bi_cur->lb < bo_nxt->rb) == false ) {
         bo_nxt++;
         bo_nxt_idx++;
      }
      /* get offset from start of outer bound */
      offset = VECTOR_INT_Get( smx->omap_cur, bo_nxt_idx );
      /* compute difference from start of outer and inner */
      offset_diff = (bi_cur->lb - bo_nxt->lb) * smx->D3;
      /* add to current row inner map */
      VECTOR_INT_Pushback( smx->imap_nxt, offset + offset_diff );
   }
}
   
/*! FUNCTION:     MATRIX_3D_SPARSE_Map_to_Outer_Dump()
 *  SYNOPSIS:     Output map to screen. Shows bounds and data offsets.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int 
MATRIX_3D_SPARSE_Map_to_Outer_Dump(    MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                                       EDGEBOUNDS*          edg,     /* EDGEBOUNDS to map */
                                       FILE*                fp )     /* FILE pointer to output to */
{
   int         edg_N;
   int         smx_N;
   int         offset;
   int         span;
   BOUND*      bnd;
   VECTOR_INT* map;

   edg_N = EDGEBOUNDS_GetSize( edg );
   smx_N = smx->data->N;
   map   = smx->omap_cur;

   for (int i = 0; i < edg_N; i++ ) 
   {
      bnd      = EDGEBOUNDS_GetX( edg, i );
      offset   = map->data[i];
      span     = (bnd->rb - bnd->lb) * smx->D3;

      fprintf(fp, "[%d/%d] {%d: %d, %d} <-> (%d-%d)/%d\n", 
         i, edg_N, bnd->id, bnd->lb, bnd->rb, offset, offset+span, smx_N );
   }
}

/*! FUNCTION:     MATRIX_3D_SPARSE_Map_to_Inner_Dump()
 *  SYNOPSIS:     Output map to screen. Shows bounds and data offsets.
 *
 *    RETURN:     <STATUS_SUCCESS> if no errors.
 */
int 
MATRIX_3D_SPARSE_Map_to_Inner_Dump(    MATRIX_3D_SPARSE*    smx,     /* MATRIX_3D_SPARSE object */
                                       EDGEBOUNDS*          edg,     /* EDGEBOUNDS to map */
                                       FILE*                fp )     /* FILE pointer to output to */
{
   int         edg_N;
   int         smx_N;
   int         offset;
   int         span;
   BOUND*      bnd;
   VECTOR_INT* map;
   int         q_0, t_0;

   map   = smx->omap_cur;
   edg   = smx->edg_outer;
   edg_N = EDGEBOUNDS_GetSize( edg );

   fprintf( fp, "## MATRIX_3D_SPARSE\n");
   fprintf( fp, "## Size: %d\n", edg_N );

   offset = 0;
   for ( int i = 0; i < edg_N; i++ )
   {
      bnd = EDGEBOUNDS_GetX( edg, i );
      q_0 = bnd->id;
      printf("\n");

      for ( int t_0 = bnd->lb; t_0 < bnd->rb; t_0++ ) 
      {
         fprintf( fp, "m(%d,%d): %d\n", 
            q_0, t_0, offset);

         offset += 3;
      }
   }
}

