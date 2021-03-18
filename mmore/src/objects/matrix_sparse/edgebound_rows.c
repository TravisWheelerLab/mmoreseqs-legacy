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
#include "../structs.h"
#include "../../utilities/_utilities.h"
#include "../_objects.h"

/* header */
#include "_matrix_sparse.h"
#include "edgebound_rows.h"

/*! FUNCTION:  EDGEBOUND_ROWS_Create()
 *  SYNOPSIS:  Create new EDGEBOUND_ROWS object and returns pointer.
 */
EDGEBOUND_ROWS* 
EDGEBOUND_ROWS_Create()
{
   EDGEBOUND_ROWS* edg;
   edg = EDGEBOUND_ROWS_Create_by_Size( 1, 1 );
   return edg;
}

/*! FUNCTION:  EDGEBOUNDS_Create_by_Size()
 *  SYNOPSIS:  Create new EDGEBOUNDS object with chosen size and returns pointer.
 *             Caller must call EDGEBOUND_ROWS_Reuse() before use.
 */
EDGEBOUND_ROWS* 
EDGEBOUND_ROWS_Create_by_Size(   int      Q,
                                 int      T )
{
   EDGEBOUND_ROWS*   edg      = NULL;

   edg = ERROR_malloc( sizeof(EDGEBOUND_ROWS) );

   edg->Q         = Q;
   edg->T         = T;
   edg->Q_range   = (RANGE){0,0};

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
                        int               T,
                        RANGE             Q_range )
{
   edg->Q         = Q;
   edg->T         = T;
   edg->Q_range   = Q_range;
   int Q_size     = Q_range.end - Q_range.beg + 1;

   EDGEBOUND_ROWS_GrowTo( edg, Q_size );
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
 *  SYNOPSIS: Resizes EDGEBOUND_ROWS to <size>.
 */
void 
EDGEBOUND_ROWS_Resize(     EDGEBOUND_ROWS*  edg,
                           int              size )
{
   /* allocate enough space for (SIZE * ROW_MAX) bounds */
   edg->rows_N = ERROR_realloc( edg->rows_N, sizeof(int) * size );
   edg->rows   = ERROR_realloc( edg->rows, sizeof(BOUND) * ( size * edg->row_max ) );
   edg->Nalloc = size;
   /* rows_N holds the number of rows */
   edg->N = size;
}

/*! FUNCTION: EDGEBOUND_ROWS_GetRowSize()
 *  SYNOPSIS: Get the size of row <q_0>.
 */
inline
int 
EDGEBOUND_ROWS_GetRowSize(   EDGEBOUND_ROWS*   edg,
                              int               q_0 )
{
   int qx0 = q_0 - edg->Q_range.beg;
   return edg->rows_N[qx0];
}

/*! FUNCTION: EDGEBOUND_ROWS_Get()
 *  SYNOPSIS: Return pointer to EDGEBOUND for absolute index <i>.
 *            Should not be called directly.
 */
inline
BOUND 
EDGEBOUND_ROWS_Get(  EDGEBOUND_ROWS*   edg,
                     int               i )
{
   /* if debugging, do edgebound checks */
   #if DEBUG
   {
      int N = edg->N * edg->row_max;
      if ( i >= N || i < 0 ) {
         fprintf(stderr, "ERROR: EDGEBOUND_ROWS Access Out-of-Bounds\n");
         fprintf(stderr, "dim: (%d), access: (%d)\n", N, i);
         ERRORCHECK_exit(EXIT_FAILURE);
      }
   }
   #endif

   return (edg->rows[i]);
}

/*! FUNCTION: EDGEBOUND_ROWS_Get()
 *  SYNOPSIS: Return pointer to EDGEBOUND for absolute index <i>.
 *            Should not be called directly.
 */
inline
BOUND*
EDGEBOUND_ROWS_GetX(    EDGEBOUND_ROWS*   edg,
                        int               i )
{
   /* if debugging, do edgebound checks */
   #if DEBUG
   {
      int N = edg->N * edg->row_max;
      if ( i >= N || i < 0 ) {
         fprintf(stderr, "ERROR: EDGEBOUND_ROWS Access Out-of-Bounds\n");
         fprintf(stderr, "dim: (%d), access: (%d)\n", N, i);
         ERRORCHECK_exit(EXIT_FAILURE);
      }
   }
   #endif

   return &(edg->rows[i]);
}

/*! FUNCTION: EDGEBOUND_ROWS_Get_byRow()
 *  SYNOPSIS: Get pointer to <i_0>th bound on <q_0>th row.
 */
inline
BOUND 
EDGEBOUND_ROWS_Get_byRow(  EDGEBOUND_ROWS*      edg,
                           int                  q_0,
                           int                  i_0 )
{
   BOUND    bnd;
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
   bnd   = EDGEBOUND_ROWS_Get( edg, idx );
   return bnd;
}

/*! FUNCTION: EDGEBOUND_ROWS_Get_byRow()
 *  SYNOPSIS: Get pointer to <i_0>th bound on <q_0>th row.
 */
inline
BOUND*
EDGEBOUND_ROWS_GetX_byRow(    EDGEBOUND_ROWS*      edg,
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
   bnd   = EDGEBOUND_ROWS_GetX( edg, idx );
   return bnd;
}

/*! FUNCTION: EDGEBOUND_ROWS_GetXLast_byRow()
 *  SYNOPSIS: Gets pointer to the last bound on <q_0>th row.
 */
inline
BOUND*
EDGEBOUND_ROWS_GetXLast_byRow(   EDGEBOUND_ROWS*      edg,
                                 int                  q_0 )
{
   BOUND*    bnd;
   int      size;

   size  = EDGEBOUND_ROWS_GetRowSize( edg, q_0 );
   if (size == 0) {
      return NULL;
   }
   bnd   = EDGEBOUND_ROWS_GetX_byRow( edg, q_0, size - 1 );
   return bnd;
}



/*! FUNCTION:  EDGEBOUND_ROWS_Pushback()
 *  SYNOPSIS:  Add BOUND <bnd> to EDGEBOUND_ROWS list at row index <row_id>.
 *             <bnd> need to added in sorted order.
 */
void 
EDGEBOUND_ROWS_Pushback(   EDGEBOUND_ROWS*   edg,
                           int               q_0,
                           BOUND             bnd )
{
   BOUND*   edg_bnd;
   int      qx0;
   int      last_idx;
   
   /* find map of index in full Q -> index in Q_range */
   qx0      = q_0 - edg->Q_range.beg;
   /* last index points to the next free bound in list */
   last_idx = EDGEBOUND_ROWS_GetRowSize(edg, q_0);

   /* if bound exceeds limit for bounds per row, throw flag */
   if ( last_idx >= edg->row_max ) 
   { 
      fprintf(stderr, "ERROR: Number of edgebounds for given row (%d) has been exceeded to (%d).\n", q_0, edg->row_max );
      fprintf(stderr, "HELP: To increase the allowed number of edgebounds in row, increase value of compiler flag '%s'.\n", "MAX_BOUNDS_PER_ROW_SUPPORTED" );

      #if DEBUG
      {
         char s[100];
         fprintf(stderr, "BOUND ATTEMPTING TO INSERT:\n");
         fprintf(stderr, "%s\n", BOUND_ToString( bnd, s ) );
         fprintf(stderr, "BOUND LIST AT ROW_ID:\n");
         for ( int i = 0; i < edg->rows_N[qx0]; i++ ) 
         {
            fprintf(stderr, "{%d}:\t", i );
            fprintf(stderr, "%s,\n", BOUND_ToString( edg->rows[i], s ) );
         }
      }
      #endif
      
      /* Optionally, we could force user to recompile with larger edgebound lists */
      // ERRORCHECK_exit(EXIT_FAILURE);
      /* Otherwise, we could just bridge the last two spans together, by updating right bound */
      last_idx       -= edg->row_max - 1;
      edg_bnd        = EDGEBOUND_ROWS_GetX_byRow( edg, q_0, last_idx );
      edg_bnd->lb    = MIN( edg_bnd->lb, bnd.lb );
      edg_bnd->rb    = MAX( edg_bnd->rb, bnd.rb );
   }
   else {
       /* get reference to bound on row at last index and update it with new info */
      edg_bnd        = EDGEBOUND_ROWS_GetX_byRow( edg, q_0, last_idx );
      edg_bnd->id    = bnd.id;
      edg_bnd->rb    = bnd.rb;
      edg_bnd->lb    = bnd.lb;
      /* increment the number of bounds in row */
      ARR_X( edg->rows_N, qx0 ) += 1;
   }
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
      last_bnd = EDGEBOUND_ROWS_GetXLast_byRow( edg, q_0 );

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
         EDGEBOUND_ROWS_Pushback( edg, q_0, new_bnd );
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
      last_bnd = EDGEBOUND_ROWS_GetXLast_byRow( edg, q_0 );

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
         EDGEBOUND_ROWS_Pushback( edg, q_0, new_bnd );
      }
   }
}

/** TODO: WIP */
/*! FUNCTION: EDGEBOUND_ROWS_Count()
 *  SYNOPSIS: Count number of cells covered by <edg>.
 */
int 
EDGEBOUND_ROWS_Count(    EDGEBOUND_ROWS*   edg )
{
   int count;

   return count;
}

/*! FUNCTION: EDGEBOUND_ROWS_Convert()
 *  SYNOPSIS: Convert EDGEBOUND_ROWS <edg_in> to EDGEBOUNDS <edg_out>.
 */
void 
EDGEBOUND_ROWS_Convert(    EDGEBOUND_ROWS*   edg_in,
                           EDGEBOUNDS*       edg_out )
{
   /* for every row in <edg_in> */
   for ( int q_0 = edg_in->Q_range.beg; q_0 < edg_in->Q_range.end; q_0++ ) 
   {
      int row_size = EDGEBOUND_ROWS_GetRowSize( edg_in, q_0 );
      /* for every bound in row */
      for ( int i_0 = 0; i_0 < row_size; i_0++ ) 
      {
         BOUND bnd = EDGEBOUND_ROWS_Get_byRow( edg_in, q_0, i_0 );
         // BOUND_Validate( edg_in, bnd );
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
      ERRORCHECK_exit(EXIT_FAILURE);
      return;
   }

   fprintf(fp, "\n");
   fprintf(fp, "N: %d, Nalloc: %d\n", edg->N, edg->Nalloc);
   for ( int i = 0; i < edg->N; i++ )
   {
      for ( int j = 0; j < edg->rows_N[i]; j++ ) 
      {
         BOUND* bnd  = EDGEBOUND_ROWS_GetX_byRow( edg, i, j );
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
         BOUND bnd_a = EDGEBOUND_ROWS_Get_byRow( edg_a, i, j );
         BOUND bnd_b = EDGEBOUND_ROWS_Get_byRow( edg_b, i, j );
         int cmp = BOUND_Compare( bnd_a, bnd_b );
         if ( cmp != 0 ) {
            return cmp;
         }
      }
   }

   return 0;
}

/*! FUNCTION: BOUND_Validate()
 *  SYNOPSIS: Verifies that <bnd> is a valid entry in <edg>.
 */
int
BOUND_Validate(   EDGEBOUND_ROWS*      edg,
                  BOUND*               bnd )
{
   bool id_check, q_check, lb_check, rb_check, b_check, passed;

   id_check = IS_IN_RANGE( 0, edg->Q, bnd->id );
   q_check  = IS_IN_RANGE( edg->Q_range.beg, edg->Q_range.end, bnd->id );
   lb_check = IS_IN_RANGE( 0, edg->T + 1, bnd->lb );
   rb_check = IS_IN_RANGE( 0, edg->T + 1, bnd->rb );
   b_check  = bnd->lb <= bnd->rb;
   passed   = ( id_check && q_check && lb_check && rb_check && b_check );

   if (!passed) {
      fprintf( stderr, "FAILED: bnd failed test(%d|%d|%d|%d|%d) :: bnd(%d,%d,%d) in edg(%d,%d)=>Q(%d,%d).\n", 
      id_check, q_check, lb_check, rb_check, b_check,
      bnd->id, bnd->lb, bnd->rb, edg->Q, edg->T, edg->Q_range.beg, edg->Q_range.end );
   }
   else {
      // fprintf( stderr, "PASSED: bnd passed test(%d|%d|%d|%d|%d) -> bnd(%d,%d,%d) in edg(%d,%d)=>Q(%d,%d).\n", 
      // id_check, q_check, lb_check, rb_check, b_check,
      // bnd->id, bnd->lb, bnd->rb, edg->Q, edg->T, edg->Q_range.beg, edg->Q_range.end );
   }

   return passed;
}

/*! FUNCTION: EDGEBOUND_ROWS_Stats()
 *  SYNOPSIS: Examine the number of bounds in each row, number of cells, and aggregate.
 */
int 
EDGEBOUND_ROWS_Stats(  EDGEBOUND_ROWS*    edg )
{
   /* total bounds */
   int bnd_total  = 0;
   int test_total = 0;
   /* total cells */
   int cell_total = 0; 
   /* occupancy of each row */
   int occ[MAX_BOUNDS_PER_ROW + 2];

   for (int i = 0; i < MAX_BOUNDS_PER_ROW; i++) {
      occ[i] = 0;
   }

   /* for every row in <edg_in> */
   for ( int q_0 = edg->Q_range.beg; q_0 < edg->Q_range.end; q_0++ ) 
   {
      int row_size = EDGEBOUND_ROWS_GetRowSize( edg, q_0 );
      if (row_size > MAX_BOUNDS_PER_ROW) {
         occ[row_size+1] += 1;
      }
      else {
         occ[row_size] += 1;
      }
      
      bnd_total += row_size;
      
      for (int i_0 = 0; i_0 < row_size; i_0++) {
         BOUND* bnd = EDGEBOUND_ROWS_GetX_byRow( edg, q_0, i_0 );
         cell_total += (bnd->rb - bnd->lb);
      }
   }

   /* output stats */
   printf("EDGEBOUND_ROWS => bnd_total: %d, cell_total: %d\n", bnd_total, cell_total);
   for (int i = 0; i < MAX_BOUNDS_PER_ROW; i++) {
      test_total += (occ[i] * i);
      printf("OCC[%d]: %d\n", i, occ[i]);
   }
   printf("OCC[OVER]: %d\n", occ[MAX_BOUNDS_PER_ROW+1]);
   if (test_total != bnd_total) printf("ERROR: sum_total != total\n");
}