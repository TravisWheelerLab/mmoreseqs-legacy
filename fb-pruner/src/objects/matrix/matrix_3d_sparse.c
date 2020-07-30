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
 *  FUNCTION: 	MATRIX_3D_SPARSE_Create_To_Contain_Edgebounds()
 *  SYNOPSIS: 	Creates sparse matrix.
 *
 *  ARGS:      <edg>     EDGEBOUNDS object.
 *
 *  RETURN:    <MATRIX_3D_SPARSE*> if no errors. Otherwise, NULL.
 */
MATRIX_3D_SPARSE* MATRIX_3D_SPARSE_Create()
{
	MATRIX_3D_SPARSE* smx = NULL;

	smx = (MATRIX_3D_SPARSE*) malloc( sizeof(MATRIX_3D_SPARSE) );
	if (smx == NULL) {
		printf("ERROR: Memory failure for allocating MATRIX_3D_SPARSE.\n");
		exit(EXIT_FAILURE);
	}

   smx->D1 = 0;
   smx->D2 = 0;
   smx->D3 = 0;

   smx->N         = 0;
   smx->Nalloc    = 0;

   smx->edg_inner = EDGEBOUNDS_Create();
   smx->edg_outer = EDGEBOUNDS_Create();

   smx->data   = VECTOR_FLT_Create();
   smx->clean  = false;

	return smx;
}


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
                                             EDGEBOUNDS*          edg_inner )
{
   smx->edg_inner = EDGEBOUNDS_Copy( smx->edg_outer, edg_inner );
	smx->edg_outer = EDGEBOUNDS_Create_Outer_Edgebounds( edg_inner, smx->edg_outer );
   /* */
   MATRIX_3D_SPARSE_Map_to_Edgebounds( smx, smx->edg_outer );
   /* */
   VECTOR_FLT_Resize( smx->data, smx->N );
}

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
                                                   EDGEBOUNDS*    edg_outer )
{
   BOUND          *bnd_0, *bnd_1, *bnd_2;       /* bounds from edgebound list, for unioning */   
   int            Q, T, N;                      /* query, target, and edgebounds length */ 
   int            q_min, q_max;                 /* min and max rows in edgebounds */
   int            r_0, r_0b, r_0e;              /* (start,end] range in upper edgebound row */
   int            r_1, r_1b, r_1e;              /* (start,end] range in center edgebound row */
   int            r_2, r_2b, r_2e;              /* (start,end] range in lower edgebound row */
   int            q_0, q_1, q_2;                /* row index */
   int            qx0, qx1, qx2;                /* row index (mod 3) */
   int 				n_[3];								/* counts of range arrays */
   RANGE          ranges[3][MAX_BOUNDS_PER_ROW * 3];   /* bound ranges on prev, cur, and next row */
   RANGE          rng; 
   bool           merging, is_Q; 

   /* embedded matrix dimension */
   Q = edg_inner->Q;
   T = edg_inner->T;
   N = edg_inner->N;
   /* first and last rows in edgebounds */
   q_min = EDG_X( edg_inner, 0 ).id;
   q_max = calc_Min( EDG_X( edg_inner, edg_inner->N ).id, Q-1 );

   /* set ranges to empty */
   for (int i = 0; i < 3; i++)
   	n_[i] = 0;

   /* create same size so likely no reallocation necessary */
   if ( edg_outer != NULL ) {
      edg_outer = EDGEBOUNDS_Create_by_Size( edg_inner->Nalloc );
   } else {
      EDGEBOUNDS_Reuse( edg_outer, edg_inner->Q, edg_inner->T );
   }
   edg_outer->Q        = edg_inner->Q;
   edg_outer->T        = edg_inner->T;
   edg_outer->edg_mode = edg_inner->edg_mode;

   /* add zeroth row if in edgebounds */
   if ( q_min == 0 )
   {
      q_1 = 0;
      qx1 = q_1 % 3;
      q_0 = 1;
      qx0 = q_0 % 3;
      /* add ranges to rows */
      r_1b = r_1 = 0;
      while ( bnd_1 = &EDG_X( edg_inner, r_1 ), (r_1 < edg_inner->N) && (bnd_1->id == q_1) ) 
      {
         ranges[qx1][n_[qx1]] = (RANGE){ bnd_1->lb-1, bnd_1->rb+1 };   /* add to its own row (+/-1 for delete state in fwd/bck) */
         ranges[qx0][n_[qx0]] = (RANGE){ bnd_1->lb, bnd_1->rb };       /* add to its next neighbors row (for insert and match state) */
         r_1++; 
         n_[qx1]++; n_[qx0]++;
      } 
      r_1e = r_1;
   }

   /* add internal rows  */
   for ( q_1; q_1 < edg_inner->Q - 1; q_1 ) 
   {
      qx1 = q_1 % 3;
      q_2 = q_1 - 1;
      qx2 = q_2 % 3;
      q_0 = q_1 + 1;
      qx0 = q_0 % 3;

      /* add ranges to rows */
      r_1b = r_1;
      while ( bnd_1 = &EDG_X( edg_inner, r_1 ), (r_1 < edg_inner->N) && (bnd_1->id == q_1) ) 
      {
         ranges[qx2][n_[qx2]] = (RANGE){ bnd_1->lb, bnd_1->rb };       /* add to its previous neighbors row (for insert and match state) */
         ranges[qx1][n_[qx1]] = (RANGE){ bnd_1->lb-1, bnd_1->rb+1 };   /* add to its own row (+/-1 for delete state in fwd/bck) */
         ranges[qx0][n_[qx0]] = (RANGE){ bnd_1->lb, bnd_1->rb };       /* add to its next neighbors row (for insert and match state) */
         r_1++; 
         n_[qx0]++; n_[qx1]++; n_[qx2]++;
      }
      r_1e = r_1;

      /* sort previous range */
      for ( int i = 0; i < n_[qx2]; i++ ) 
      {
         RANGE*   cur_r    = NULL;
         RANGE*   min_r    = &ranges[qx2][i]; 
         int      min_i    = i;
         int      cmp      = 0;
         for ( int j = i; j < n_[qx2]; i++ ) 
         {
            cur_r = &ranges[qx2][j];
            cmp = RANGE_Compare( min_r, cur_r );
            if ( cmp > 0 ) { min_r = cur_r; }
            min_i = j;
         }
         RANGE swap           = ranges[qx2][i];
         ranges[qx2][i]       = ranges[qx2][min_i];
         ranges[qx2][min_i]   = swap;
      }

      /* merge and add to edgebounds */
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
      BOUND bnd = (BOUND) { q_2, ranges[qx2][n_[qx2]-1].beg, ranges[qx2][n_[qx2] -1].end };
      EDGEBOUNDS_Pushback( edg_outer, &bnd );

      n_[qx2] = 0;
   }

	/* add zeroth row if in edgebounds */
   q_1 = Q;
   qx1 = q_1 % 3;
   q_2 = 1;
   qx2 = q_2 % 3;
   /* add ranges to rows */
   r_1b = r_1 = 0;
   while ( bnd_1 = &EDG_X( edg_inner, r_1 ), (r_1 < edg_inner->N) && (bnd_1->id == q_1) ) 
   {
      ranges[qx2][n_[qx2]] = (RANGE){ bnd_1->lb, bnd_1->rb };   		/* add to its prev neighbors row (for insert and match state) */ 
      ranges[qx1][n_[qx1]] = (RANGE){ bnd_1->lb-1, bnd_1->rb+1 };     /* add to its own row (+/-1 for delete state in fwd/bck) */
      r_1++; 
      n_[qx0]++; n_[qx1]++;
   } 
   r_1e = r_1;

   return edg_outer;
}


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
                                          EDGEBOUNDS*          edg )
{
   int      N;          /* total edgebounds */
   int      count;      /* cell count */
   int      offset;     /* offset into row */
   int      r_0;        /* current row */
   BOUND*   b_0;        /* current bound */

   N      = edg->N;
   count  = 0;
   offset = 0;

   b_0 = &(edg->bounds[0]);
   r_0 = b_0->id; /* first rows */

   VECTOR_INT_Pushback( smx->rows, 0 );
   VECTOR_INT_Pushback( smx->offsets, count );

   for (int i = 0; i < N; i++, b_0++) 
   {
      /* check if in new row */
      if ( b_0->id == r_0 ) {
         VECTOR_INT_Pushback( smx->rows, i );
         r_0 = b_0->id;
      }
      /* add to offsets */
      count += ( b_0->rb - b_0->lb );
      VECTOR_INT_Pushback( smx->offsets, count );
   }

   VECTOR_INT_Pushback( smx->rows, N );
   VECTOR_INT_Pushback( smx->offsets, count );

   smx->N = count;
}