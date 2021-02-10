/*******************************************************************************
 *  FILE:      edgebound_rows.c
 *  PURPOSE:   EDGEBOUND_ROWS Object.
 *             For building EDGEBOUNDS on the fly.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "structs.h"
#include "../utilities/_utilities.h"
#include "_objects.h"

/* header */
#include "edgebound_rows.h"

/*! FUNCTION:  EDGEBOUND_ROWS_Create()
 *  SYNOPSIS:  Create new EDGEBOUND_ROWS object and returns pointer.
 */
EDGEBOUND_ROWS* 
EDGEBOUND_ROWS_Create()
{
   EDGEBOUND_ROWS* edg = EDGEBOUND_ROWS_Create_by_Size( 1, 1 );
   return edg;
}

/*! FUNCTION:  EDGEBOUNDS_Create_by_Size()
 *  SYNOPSIS:  Create new EDGEBOUNDS object with chosen size and returns pointer.
 */
EDGEBOUND_ROWS* 
EDGEBOUND_ROWS_Create_by_Size(   int   Q, 
                                 int   T )
{
   EDGEBOUND_ROWS*   edg      = NULL;

   edg = ERROR_malloc( sizeof(EDGEBOUND_ROWS) );

   edg->Q         = Q;
   edg->T         = T;

   edg->N         = 0;
   edg->Nalloc    = 0;
   edg->rows_N    = NULL;
   edg->rows      = NULL;
   edg->row_max   = MAX_BOUNDS_PER_ROW;

   EDGEBOUND_ROWS_Resize( edg, Q+1 );
   return edg;
}

/*! FUNCTION: EDGEBOUND_ROWS_Destroy()
 *  SYNOPSIS: Frees all memory from EDGEBOUND_ROWS object.
 */
EDGEBOUND_ROWS* 
EDGEBOUND_ROWS_Destroy( EDGEBOUND_ROWS*  edg )
{
   if ( edg == NULL ) return edg;

   ERROR_free( edg->rows );
   ERROR_free( edg->rows_N );
   
   ERROR_free( edg );
   edg = NULL;
   return edg;
}

/*! FUNCTION: EDGEBOUND_ROWS_Reuse()
 *  SYNOPSIS: Reuses EDGEBOUND_ROWS by resizing if too a"clearing" edgebound list (does not downsize).
 */
void 
EDGEBOUND_ROWS_Reuse(   EDGEBOUND_ROWS*   edg,
                        int               Q,
                        int               T )
{
   edg->Q = Q;
   edg->T = T;

   EDGEBOUND_ROWS_Resize( edg, Q+1 );
   EDGEBOUND_ROWS_Clear( edg );
}

/*! FUNCTION: EDGEBOUND_ROWS_Clear()
 *  SYNOPSIS: Reuses EDGEBOUND_ROWS by "clearing" edgebound list (does not downsize).
 */
void EDGEBOUND_ROWS_Clear( EDGEBOUND_ROWS*   edg )
{
   for ( int i = 0; i < edg->N; i++ ) {
      edg->rows_N[i] = 0;
   }
}

/*! FUNCTION: EDGEBOUND_ROWS_GrowTo()
 *  SYNOPSIS: Resizes EDGEBOUND_ROWS if new size exceeds old size
 */
void EDGEBOUND_ROWS_GrowTo( EDGEBOUND_ROWS*  edg,
                            int              size )
{
   if ( size > edg->Nalloc ) {
      EDGEBOUND_ROWS_Resize( edg, size );
   }
   edg->N = size;
}

/*! FUNCTION: EDGEBOUND_ROWS_Resize()
 *  SYNOPSIS: Resizes EDGEBOUND_ROWS if new size exceeds old size
 */
void 
EDGEBOUND_ROWS_Resize(     EDGEBOUND_ROWS*  edg,
                           int              size )
{
   if ( size > edg->Nalloc ) 
   {
      edg->rows_N = ERROR_realloc( edg->rows_N, sizeof(int) * size );
      edg->rows   = ERROR_realloc( edg->rows, sizeof(BOUND) * ( size * edg->row_max ) );
      edg->Nalloc = size;
   }
   /* rows_N holds the number of rows */
   edg->N = size;
}

/*! FUNCTION: EDGEBOUND_ROWS_Get_RowSize()
 *  SYNOPSIS: Get the size of row <q_0>.
 */
inline
int 
EDGEBOUND_ROWS_Get_RowSize(   EDGEBOUND_ROWS*   edg,
                              int               q_0 )
{
   return edg->rows_N[q_0];
}

/*! FUNCTION: EDGEBOUND_ROWS_Get()
 *  SYNOPSIS: Return pointer to EDGEBOUND for absolute index <i>.
 */
inline
BOUND* 
EDGEBOUND_ROWS_Get(  EDGEBOUND_ROWS*   edg,
                     int               i )
{
   /* if debugging, do edgebound checks */
   #if DEBUG
      if ( i >= edg->N || i < 0 ) {
         fprintf(stderr, "ERROR: EDGEBOUND_ROWS Access Out-of-Bounds\n");
         fprintf(stderr, "dim: (%d/%d), access: (%d)\n", edg->N, edg->Nalloc, i);
         exit(EXIT_FAILURE);
      }
   #endif

   return &(edg->rows[i]);
}

/*! FUNCTION: EDGEBOUND_ROWS_Get_byRow()
 *  SYNOPSIS: Get pointer to <i_0>th bound on <q_0>th row.
 */
inline
BOUND* 
EDGEBOUND_ROWS_Get_byRow(  EDGEBOUND_ROWS*      edg,
                           int                  q_0,
                           int                  i_0 )
{
   /* convert to flat index to access proper row */
   int idx     = (q_0 * edg->row_max) + i_0;
   BOUND* bnd  = EDGEBOUND_ROWS_Get( edg, idx );
   return bnd;
}

/*! FUNCTION: EDGEBOUND_ROWS_GetLast_byRow()
 *  SYNOPSIS: Gets pointer to the last bound on <q_0>th row.
 */
inline
BOUND* 
EDGEBOUND_ROWS_GetLast_byRow(    EDGEBOUND_ROWS*      edg,
                                 int                  q_0 )
{
   BOUND*   bnd;
   int      size;

   size  = EDGEBOUND_ROWS_Get_RowSize( edg, q_0 );
   if ( size == 0 ) {
      return NULL;
   }
   bnd   = EDGEBOUND_ROWS_Get_byRow( edg, q_0, size - 1 );
   return bnd;
}

/*! FUNCTION: EDGEBOUND_ROWS_Pushback()
 *  SYNOPSIS: Add BOUND <bnd> to EDGEBOUND_ROWS list at row index <row_id>.
 */
void 
EDGEBOUND_ROWS_Pushback(   EDGEBOUND_ROWS*   edg,
                           int               row_id,
                           BOUND*            bnd )
{
   /* if bound exceeds limit for bounds per row, throw flag */
   if ( edg->rows_N[row_id] >= edg->row_max - 1 ) 
   { 
      fprintf(stderr, "ERROR: Number of edgebounds for given row (%d) has been exceeded to (%d).\n", row_id, edg->rows_N[row_id] );
      fprintf(stderr, "HELP: To increase the allowed number of edgebounds in row, increase value of compiler flag '%s'.\n", "MAX_BOUNDS_PER_ROW_SUPPORTED" );

      #if DEBUG
      {
         char s[100];
         fprintf(stderr, "BOUND ATTEMPTING TO INSERT:\n");
         fprintf(stderr, "%s\n", BOUND_To_String( *bnd, s ) );
         fprintf(stderr, "BOUND LIST AT ROW_ID:\n");
         for ( int i = 0; i < edg->rows_N[row_id]; i++ ) 
         {
            fprintf(stderr, "{%d}:\t", i );
            fprintf(stderr, "%s,\n", BOUND_To_String( edg->rows[i], s ) );
         }
      }
      #endif
      
      exit(EXIT_FAILURE);
   }

   BOUND*   edg_bnd;
   edg_bnd  = EDGEBOUND_ROWS_Get_byRow( edg, row_id, edg->rows_N[row_id] );
   *edg_bnd = *bnd;
   edg->rows_N[row_id]++;

   // printf("Added index: %d to row_id: %d\n", edg->rows_N[row_id], row_id);
   // BOUND_Dump( edg_bnd, stdout );
   // BOUND_Dump( bnd, stdout );
}

/*! FUNCTION: EDGEBOUND_ROWS_IntegrateDiag_Fwd()
 *  SYNOPSIS: Add antidiagonal bound into row-wise bounds, for the Forward Cloud Search.
 *            Looks at each cell individually in the antidiagonal.
 *            If it is right-side adjacent to the current open bound (within a tolerance value), it extends it.
 *            Otherwise, it creates a new edgebound and adds it to the list.
 *            Edgebound lists sizes are determined at compile time and do not resize.
 *            If size is exceeded, program terminates with error.
 */
void 
EDGEBOUND_ROWS_IntegrateDiag_Fwd(   EDGEBOUND_ROWS*   edg,
                                    BOUND*            bnd )
{
   BOUND*   last_bnd;
   BOUND    new_bnd;

   /* max distance between edgebounds that will be merged */
   const int tol = 0;   
   /* antidiagonal of bound */
   int d_0 = bnd->id;

   /* all cells which antidiag intersects */
   for ( int k_0 = bnd->lb; k_0 < bnd->rb; k_0++ ) 
   {
      /* row-wise coords */
      int q_0 = k_0;
      int t_0 = d_0 - q_0;

      /* get latest bound in requested row */
      last_bnd = EDGEBOUND_ROWS_GetLast_byRow( edg, q_0 );

      /* if row is not empty AND row bounds is adjacent (within tolerance) to new cell, merge them */
      if ( last_bnd != NULL && t_0 <= last_bnd->rb + tol ) 
      {
         // fprintf(stderr, "Merging, %d is in range %d:<%d-%d>\n", j, row->id, row->lb, row->rb );
         last_bnd->rb = t_0+1;
      }
      /* otherwise, create new bound and to row */
      else
      {
         // fprintf(stderr, "Creating new, %d is outside range %d:<%d-%d>\n", j, row->id, row->lb, row->rb );
         new_bnd = (BOUND){ q_0, t_0, t_0+1};
         EDGEBOUND_ROWS_Pushback( edg, q_0, &new_bnd );
      }
   }
}

/*! FUNCTION: EDGEBOUND_ROWS_IntegrateDiag_Bck()
 *  SYNOPSIS: Add antidiagonal bound into row-wise bounds, for the Forward Cloud Search.
 *            Looks at each cell individually in the antidiagonal.
 *            If it is left-side adjacent to the current open bound (within a tolerance value), it extends it.
 *            Otherwise, it creates a new edgebound and adds it to the list.
 *            Edgebound lists sizes are determined at compile time and do not resize.
 *            If size is exceeded, program terminates with error.
 */
void 
EDGEBOUND_ROWS_IntegrateDiag_Bck(   EDGEBOUND_ROWS*   edg,
                                    BOUND*            bnd )
{
   BOUND*   last_bnd;
   BOUND    new_bnd;

   /* max distance between edgebounds that will be merged */
   const int tol = 0;   
   /* antidiagonal of bound */
   int d_0 = bnd->id;

   /* all cells which antidiag intersects */
   for ( int k_0 = bnd->lb; k_0 < bnd->rb; k_0++ ) 
   {
      /* row-wise coords */
      int q_0 = k_0;
      int t_0 = d_0 - q_0;

      /* get latest bound in requested row */
      last_bnd = EDGEBOUND_ROWS_GetLast_byRow( edg, q_0 );

      /* if row is not empty AND row bounds is adjacent (within tolerance) to new cell, merge them */
      if ( ( last_bnd != NULL ) && ( t_0 >= last_bnd->lb - tol - 1 ) ) 
      {
         // fprintf(stderr, "Merging, %d is in range %d:<%d-%d>\n", j, row->id, row->lb, row->rb );
         last_bnd->rb = t_0+1;
      }
      /* otherwise, create new bound and to row */
      else
      {
         // fprintf(stderr, "Creating new, %d is outside range %d:<%d-%d>\n", j, row->id, row->lb, row->rb );
         new_bnd = (BOUND){ q_0, t_0, t_0+1};
         EDGEBOUND_ROWS_Pushback( edg, q_0, &new_bnd );
      }
   }
}

/*! FUNCTION: EDGEBOUND_ROWS_Convert()
 *  SYNOPSIS: Convert EDGEBOUNDS_ROWS <edg_in> to EDGEBOUND <edg_out>.
 */
void 
EDGEBOUND_ROWS_Convert(    EDGEBOUND_ROWS*   edg_in,
                           EDGEBOUNDS*       edg_out )
{
   EDGEBOUNDS_Reuse( edg_out, edg_in->Q, edg_in->T );
   edg_out->edg_mode = EDG_ROW;

   for ( int i = 0; i < edg_in->N; i++ ) 
   {
      for ( int j = 0; j < edg_in->rows_N[i]; j++ ) 
      {
         BOUND* bnd = EDGEBOUND_ROWS_Get_byRow( edg_in, i, j );
         EDGEBOUNDS_Pushback( edg_out, bnd );
      }
   }
}

/*! FUNCTION: EDGEBOUND_ROWS_Dump()
 *  SYNOPSIS: Print EDGEBOUND object to file.
 */
void 
EDGEBOUND_ROWS_Dump(    EDGEBOUND_ROWS*   edg,
                        FILE*             fp )
{
   /* test for bad file pointer */
   if (fp == NULL) {
      const char* obj_name = "EDGEBOUND_ROWS";
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing %s.\n", obj_name);
      exit(EXIT_FAILURE);
      return;
   }

   fprintf(fp, "\n");
   fprintf(fp, "N: %d, Nalloc: %d\n", edg->N, edg->Nalloc);
   for ( int i = 0; i < edg->N; i++ )
   {
      for ( int j = 0; j < edg->rows_N[i]; j++ ) 
      {
         BOUND* bnd = EDGEBOUND_ROWS_Get_byRow( edg, i, j );
         fprintf(fp, "[%d] ", i);
         fprintf(fp, "{ id: %d, lb: %d, rb: %d }\n", bnd->id, bnd->lb, bnd->rb);
      }
   }
   fprintf(fp, "\n");
}

/*! FUNCTION: EDGEBOUND_ROWS_Compare()
 *  SYNOPSIS: Compare two EDGEBOUND_ROWS objects.  Return 0 if equal.
 */
int 
EDGEBOUND_ROWS_Compare(    EDGEBOUND_ROWS*    edg_a,
                           EDGEBOUND_ROWS*    edg_b )
{
   if ( edg_a->N != edg_b->N ) return -1;

   for ( int i = 0 ; i < edg_a->N; i++ ) 
   {
      for ( int j = 0; j < edg_a->rows_N[i]; j++ )
      {
         BOUND* bnd_a = EDGEBOUND_ROWS_Get_byRow( edg_a, i, j );
         BOUND* bnd_b = EDGEBOUND_ROWS_Get_byRow( edg_b, i, j );
         int cmp = BOUND_Compare( *bnd_a, *bnd_b );
         if ( cmp != 0 ) {
            return cmp;
         }
      }
   }
}
