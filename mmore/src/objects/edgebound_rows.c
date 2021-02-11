/*******************************************************************************
 *  FILE:      edgebound_rows.c
 *  PURPOSE:   EDGEROWS Object.
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

/* header */
#include "_objects.h"
#include "edgebound_rows.h"

/*! FUNCTION:  EDGEROWS_Create()
 *  SYNOPSIS:  Create new EDGEROWS object and returns pointer.
 */
EDGEROWS* 
EDGEROWS_Create()
{
   EDGEROWS* edg = EDGEROWS_Create_by_Size( 1, 1 );
   return edg;
}

/*! FUNCTION:  EDGEBOUNDS_Create_by_Size()
 *  SYNOPSIS:  Create new EDGEBOUNDS object with chosen size and returns pointer.
 *             Caller must call EDGEROWS_Reuse() before use.
 */
EDGEROWS* 
EDGEROWS_Create_by_Size(   int      Q,
                                 int      T )
{
   EDGEROWS*   edg      = NULL;

   edg = ERROR_malloc( sizeof(EDGEROWS) );

   edg->Q         = Q;
   edg->T         = T;
   edg->Q_range   = (RANGE){0,0};

   edg->N         = 0;
   edg->Nalloc    = 0;
   edg->rows_N    = NULL;
   edg->rows      = NULL;
   edg->row_max   = MAX_BOUNDS_PER_ROW;

   EDGEROWS_Resize( edg, Q+1 );
   return edg;
}

/*! FUNCTION: EDGEROWS_Destroy()
 *  SYNOPSIS: Frees all memory from EDGEROWS object.
 */
EDGEROWS* 
EDGEROWS_Destroy( EDGEROWS*  edg )
{
   if ( edg == NULL ) return edg;

   ERROR_free( edg->rows );
   ERROR_free( edg->rows_N );
   
   ERROR_free( edg );
   edg = NULL;
   return edg;
}

/*! FUNCTION: EDGEROWS_Reuse()
 *  SYNOPSIS: Reuses EDGEROWS by resizing if too a"clearing" edgebound list (does not downsize).
 */
void 
EDGEROWS_Reuse(   EDGEROWS*   edg,
                        int               Q,
                        int               T,
                        RANGE             Q_range )
{
   edg->Q         = Q;
   edg->T         = T;
   edg->Q_range   = Q_range;
   int Q_size     = Q_range.end - Q_range.beg + 1;

   EDGEROWS_GrowTo( edg, Q_size );
   EDGEROWS_Clear( edg );
}

/*! FUNCTION: EDGEROWS_Clear()
 *  SYNOPSIS: Reuses EDGEROWS by "clearing" edgebound list (does not downsize).
 */
void EDGEROWS_Clear( EDGEROWS*   edg )
{
   for ( int i = 0; i < edg->N; i++ ) {
      edg->rows_N[i] = 0;
   }
}

/*! FUNCTION: EDGEROWS_GrowTo()
 *  SYNOPSIS: Resizes EDGEROWS if new size exceeds old size
 */
void EDGEROWS_GrowTo( EDGEROWS*  edg,
                            int              size )
{
   if ( size > edg->Nalloc ) {
      EDGEROWS_Resize( edg, size );
   }
   edg->N = size;
}

/*! FUNCTION: EDGEROWS_Resize()
 *  SYNOPSIS: Resizes EDGEROWS if new size exceeds old size
 */
void 
EDGEROWS_Resize(     EDGEROWS*  edg,
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

/*! FUNCTION: EDGEROWS_Get_RowSize()
 *  SYNOPSIS: Get the size of row <q_0>.
 */
inline
int 
EDGEROWS_Get_RowSize(   EDGEROWS*   edg,
                              int               q_0 )
{
   int qx0 = q_0 - edg->Q_range.beg;
   return edg->rows_N[qx0];
}

/*! FUNCTION: EDGEROWS_Get()
 *  SYNOPSIS: Return pointer to EDGEBOUND for absolute index <i>.
 *            Should not be called directly.
 */
inline
BOUND* 
EDGEROWS_Get(  EDGEROWS*   edg,
                     int               i )
{
   /* if debugging, do edgebound checks */
   #if DEBUG
      if ( i >= edg->N || i < 0 ) {
         fprintf(stderr, "ERROR: EDGEROWS Access Out-of-Bounds\n");
         fprintf(stderr, "dim: (%d/%d), access: (%d)\n", edg->N, edg->Nalloc, i);
         exit(EXIT_FAILURE);
      }
   #endif

   return &(edg->rows[i]);
}

/*! FUNCTION: EDGEROWS_Get_byRow()
 *  SYNOPSIS: Get pointer to <i_0>th bound on <q_0>th row.
 */
inline
BOUND* 
EDGEROWS_Get_byRow(  EDGEROWS*      edg,
                           int                  q_0,
                           int                  i_0 )
{
   BOUND*   bnd;
   int      idx;
   int      qx0;

   /* if boundchecking, make sure q_0 is a valid request */
   #if SAFE
   {
      if ( IS_IN_RANGE( edg->Q_range.beg, edg->Q_range.end, q_0 ) == false ) {
         fprintf( stderr, "ERROR: EDGEBOUND requested row %d, when range is: (%d,%d)\n", 
            q_0, edg->Q_range.beg, edg->Q_range.end );
      }
   }
   #endif 

   /* find map of index in full Q -> index in Q_range */
   qx0   = q_0 - edg->Q_range.beg;

   /* convert to flat index to access proper row */
   idx   = (qx0 * edg->row_max) + i_0;
   /* get row at flat index */
   bnd   = EDGEROWS_Get( edg, idx );
   return bnd;
}

/*! FUNCTION: EDGEROWS_GetLast_byRow()
 *  SYNOPSIS: Gets pointer to the last bound on <q_0>th row.
 */
inline
BOUND* 
EDGEROWS_GetLast_byRow(    EDGEROWS*      edg,
                                 int                  q_0 )
{
   BOUND*   bnd;
   int      size;

   size  = EDGEROWS_Get_RowSize( edg, q_0 );
   if ( size == 0 ) {
      return NULL;
   }
   bnd   = EDGEROWS_Get_byRow( edg, q_0, size - 1 );
   return bnd;
}

/*! FUNCTION: EDGEROWS_Pushback()
 *  SYNOPSIS: Add BOUND <bnd> to EDGEROWS list at row index <row_id>.
 */
void 
EDGEROWS_Pushback(   EDGEROWS*   edg,
                           int               q_0,
                           BOUND*            bnd )
{
   BOUND*   edg_bnd;
   int      qx0;
   int      last_idx;
   
   /* find map of index in full Q -> index in Q_range */
   qx0      = q_0 - edg->Q_range.beg;
   /* last index points to the next free bound in list */
   last_idx = EDGEROWS_Get_RowSize(edg, q_0);

   /* if bound exceeds limit for bounds per row, throw flag */
   if ( last_idx >= edg->row_max ) 
   { 
      fprintf(stderr, "ERROR: Number of edgebounds for given row (%d) has been exceeded to (%d).\n", q_0, edg->row_max );
      fprintf(stderr, "HELP: To increase the allowed number of edgebounds in row, increase value of compiler flag '%s'.\n", "MAX_BOUNDS_PER_ROW_SUPPORTED" );

      #if DEBUG
      {
         char s[100];
         fprintf(stderr, "BOUND ATTEMPTING TO INSERT:\n");
         fprintf(stderr, "%s\n", BOUND_To_String( *bnd, s ) );
         fprintf(stderr, "BOUND LIST AT ROW_ID:\n");
         for ( int i = 0; i < edg->rows_N[qx0]; i++ ) 
         {
            fprintf(stderr, "{%d}:\t", i );
            fprintf(stderr, "%s,\n", BOUND_To_String( edg->rows[i], s ) );
         }
      }
      #endif
      
      /* Optionally, we could force user to recompile with larger edgebound lists */
      // exit(EXIT_FAILURE);
      /* Otherwise, we could just bridge the last two spans together, by updating right bound */
      last_idx -= edg->row_max - 1;
      edg_bnd  = EDGEROWS_Get_byRow( edg, q_0, last_idx );
      edg_bnd->lb = MIN( edg_bnd->lb, bnd->lb );
      edg_bnd->rb = MAX( edg_bnd->rb, bnd->rb );
   }
   else {
       /* get reference to bound on row at last index and update it with new info */
      edg_bnd  = EDGEROWS_Get_byRow( edg, q_0, last_idx );
      edg_bnd->id = bnd->id;
      edg_bnd->rb = bnd->rb;
      edg_bnd->lb = bnd->lb;
      /* increment the number of bounds in row */
      edg->rows_N[qx0]++;
   }
}

/*! FUNCTION: EDGEROWS_IntegrateDiag_Fwd()
 *  SYNOPSIS: Add antidiagonal bound into row-wise bounds, for the Forward Cloud Search.
 *            Looks at each cell individually in the antidiagonal.
 *            If it is right-side adjacent to the current open bound (within a tolerance value), it extends it.
 *            Otherwise, it creates a new edgebound and adds it to the list.
 *            Edgebound lists sizes are determined at compile time and do not resize.
 *            If size is exceeded, program terminates with error.
 */
void 
EDGEROWS_IntegrateDiag_Fwd(   EDGEROWS*   edg,
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
      last_bnd = EDGEROWS_GetLast_byRow( edg, q_0 );

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
         EDGEROWS_Pushback( edg, q_0, &new_bnd );
      }
   }
}

/*! FUNCTION: EDGEROWS_IntegrateDiag_Bck()
 *  SYNOPSIS: Add antidiagonal bound into row-wise bounds, for the Forward Cloud Search.
 *            Looks at each cell individually in the antidiagonal.
 *            If it is left-side adjacent to the current open bound (within a tolerance value), it extends it.
 *            Otherwise, it creates a new edgebound and adds it to the list.
 *            Edgebound lists sizes are determined at compile time and do not resize.
 *            If size is exceeded, program terminates with error.
 */
void 
EDGEROWS_IntegrateDiag_Bck(   EDGEROWS*   edg,
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
      last_bnd = EDGEROWS_GetLast_byRow( edg, q_0 );

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
         EDGEROWS_Pushback( edg, q_0, &new_bnd );
      }
   }
}

/*! FUNCTION: EDGEROWS_Convert()
 *  SYNOPSIS: Convert EDGEROWS <edg_in> to EDGEBOUNDS <edg_out>.
 */
void 
EDGEROWS_Convert(    EDGEROWS*   edg_in,
                           EDGEBOUNDS*       edg_out )
{
   EDGEBOUNDS_Reuse( edg_out, edg_in->Q, edg_in->T );
   edg_out->edg_mode = EDG_ROW;

   /* for every row in <edg_in> */
   for ( int q_0 = edg_in->Q_range.beg; q_0 < edg_in->Q_range.end; q_0++ ) 
   {
      int row_size = EDGEROWS_Get_RowSize( edg_in, q_0 );
      /* for every bound in row */
      for ( int i_0 = 0; i_0 < row_size; i_0++ ) 
      {
         BOUND* bnd = EDGEROWS_Get_byRow( edg_in, q_0, i_0 );
         EDGEBOUNDS_Pushback( edg_out, bnd );
      }
   }
}

/*! FUNCTION: EDGEROWS_Dump()
 *  SYNOPSIS: Print EDGEBOUND object to file.
 */
void 
EDGEROWS_Dump(    EDGEROWS*   edg,
                        FILE*             fp )
{
   /* test for bad file pointer */
   if (fp == NULL) {
      const char* obj_name = "EDGEROWS";
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
         BOUND* bnd = EDGEROWS_Get_byRow( edg, i, j );
         fprintf(fp, "[%d] ", i);
         fprintf(fp, "{ id: %d, lb: %d, rb: %d }\n", bnd->id, bnd->lb, bnd->rb);
      }
   }
   fprintf(fp, "\n");
}

/*! FUNCTION: EDGEROWS_Compare()
 *  SYNOPSIS: Compare two EDGEROWS objects.  Return 0 if equal.
 */
int 
EDGEROWS_Compare(    EDGEROWS*    edg_a,
                           EDGEROWS*    edg_b )
{
   if ( edg_a->N != edg_b->N ) return -1;

   for ( int i = 0 ; i < edg_a->N; i++ ) 
   {
      for ( int j = 0; j < edg_a->rows_N[i]; j++ )
      {
         BOUND* bnd_a = EDGEROWS_Get_byRow( edg_a, i, j );
         BOUND* bnd_b = EDGEROWS_Get_byRow( edg_b, i, j );
         int cmp = BOUND_Compare( *bnd_a, *bnd_b );
         if ( cmp != 0 ) {
            return cmp;
         }
      }
   }
}
