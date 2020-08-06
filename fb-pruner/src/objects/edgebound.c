/*******************************************************************************
 *  FILE:      edgebound.c
 *  PURPOSE:   EDGEBOUNDS Object
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
#include "utilities.h"
#include "objects.h"

/* header */
#include "edgebound.h"

/*
 *  FUNCTION:  EDGEBOUNDS_Create()
 *  SYNOPSIS:  Create new EDGEBOUNDS object and returns pointer.
 */
EDGEBOUNDS* EDGEBOUNDS_Create( void )
{
   EDGEBOUNDS* edg      = NULL;
   const int   min_size = 8;

   edg = EDGEBOUNDS_Create_by_Size( min_size );

   return edg;
}

/*
 *  FUNCTION:  EDGEBOUNDS_Create()
 *  SYNOPSIS:  Create new EDGEBOUNDS object with chosen size and returns pointer.
 */
EDGEBOUNDS* EDGEBOUNDS_Create_by_Size( const int size )
{
   EDGEBOUNDS* edg      = NULL;

   edg = (EDGEBOUNDS*) malloc( sizeof(EDGEBOUNDS) );
   if (edg == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc for EDGEBOUNDS.\n");
   }

   edg->N         = 0;
   edg->Nalloc    = 0;

   edg->Q         = 0;
   edg->T         = 0;

   edg->ids       = VECTOR_INT_Create();
   edg->ids_idx   = VECTOR_INT_Create();

   edg->ids       = VECTOR_INT_Create();
   edg->ids_idx   = VECTOR_INT_Create();
   edg->edg_mode  = EDG_NONE;

   edg->bounds    = NULL;

   EDGEBOUNDS_Resize(edg, size);

   return edg;
}

/*
 *  FUNCTION: EDGEBOUNDS_Destroy()
 *  SYNOPSIS: Frees all memory from EDGEBOUNDS object.
 */
void* EDGEBOUNDS_Destroy( EDGEBOUNDS*  edg )
{
   if ( edg == NULL ) return edg;

   VECTOR_INT_Destroy( edg->ids );
   VECTOR_INT_Destroy( edg->ids_idx );

   free( edg->bounds );
   edg->bounds = NULL;

   free( edg );
   edg = NULL;

   return edg;
}

/*
 *  FUNCTION: EDGEBOUNDS_Reuse()
 *  SYNOPSIS: Reuses EDGEBOUNDS by "clearing" edgebound list (does not realloc).
 */
void EDGEBOUNDS_Reuse( EDGEBOUNDS*   edg, 
                       int           Q,
                       int           T )
{
   edg->N = 0;
   edg->Q = Q;
   edg->T = T;
}

/*
 *  FUNCTION: EDGEBOUNDS_Set_Dim()
 *  SYNOPSIS: Set the dimensions of the embedding matrix.
 */
void EDGEBOUNDS_Set_Dim(   EDGEBOUNDS*   edg, 
                           int           Q,
                           int           T )
{
   edg->Q = Q;
   edg->T = T;
}

/*
 *  FUNCTION: EDGEBOUNDS_Copy()
 *  SYNOPSIS: Create a deep copy of <edg_src> and store it in <edg_dest>.
 */
EDGEBOUNDS* EDGEBOUNDS_Copy(  EDGEBOUNDS*          edg_dest,
                              const EDGEBOUNDS*    edg_src )
{
   BOUND*   bnd;

   if ( edg_dest == NULL ) {
      edg_dest = EDGEBOUNDS_Create_by_Size( edg_src->Nalloc );
   } 
   edg_dest->Q          = edg_src->Q;
   edg_dest->T          = edg_src->T;
   edg_dest->edg_mode   = edg_src->edg_mode;


   edg_dest->edg_mode = edg_src->edg_mode;

   for ( int i = 0; i < edg_src->N; i++ ) 
   {
      bnd = &(edg_src->bounds[i]);
      EDGEBOUNDS_Pushback( edg_dest, bnd );
   }

   return edg_dest;
}

/*
 *  FUNCTION: EDGEBOUNDS_Get()
 *  SYNOPSIS: Return pointer to BOUND at index <i>.
 */
inline
BOUND* EDGEBOUNDS_Get( EDGEBOUNDS*   edg,
                       int           i )
{
   /* if debugging, do edgebound checks */
   #if DEBUG
      if ( i >= edg->N || i < 0 ) {
         fprintf(stderr, "ERROR: EDGEBOUNDS Access Out-of-Bounds\n");
         fprintf(stderr, "dim: (%d/%d), access: (%d)\n", edg->N, edg->Nalloc, i);
         exit(EXIT_FAILURE);
      }
   #endif

   return &(edg->bounds[i]);
}

/*
 *  FUNCTION: EDGEBOUNDS_Pushback()
 *  SYNOPSIS: Add BOUND to EDGEBOUNDS list.
 */
void EDGEBOUNDS_Pushback( EDGEBOUNDS*  edg,
                          BOUND*       bnd )
{
   edg->bounds[edg->N] = *bnd;
   edg->N++;

   /* resize if necessary */
   if (edg->N >= edg->Nalloc - 1) {
      EDGEBOUNDS_Resize(edg, edg->Nalloc * 2);
   }
}

/*
 *  FUNCTION: EDGEBOUNDS_Insert()
 *  SYNOPSIS: Insert/Overwrite bound into <i> index of Edgebound list.
 */
void EDGEBOUNDS_Insert( EDGEBOUNDS*    edg,
                        int            i,
                        BOUND*         bnd )
{
   BOUND*   edg_bnd;

   edg_bnd = &(edg->bounds[i]);
   *edg_bnd = *bnd;
}

/*
 *  FUNCTION: EDGEBOUNDS_Delete()
 *  SYNOPSIS: Delete BOUND at <i> index and fill from end of list <N-1>, then decrement list size.
 */
void EDGEBOUNDS_Delete( EDGEBOUNDS*    edg,
                        int            i )
{
   int N = edg->N;
   edg->bounds[i] = edg->bounds[N-1];
   edg->N -= 1;
}

/*
 *  FUNCTION: EDGEBOUNDS_Clear()
 *  SYNOPSIS: Remove all BOUNDS from EDGEBOUND list (no realloc).
 */
void EDGEBOUNDS_Clear( EDGEBOUNDS*  edg )
{
   edg->N = 0;
}

/*
 *  FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of BOUNDS allocated in EDGEBOUND object (does not downsize).
 */
void EDGEBOUNDS_GrowTo( EDGEBOUNDS*    edg,
                        int            size )
{
   if ( edg->Nalloc < size )
      EDGEBOUNDS_Resize( edg, size );
}

/*
 *  FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of BOUNDS allocated in EDGEBOUND object .
 */
void EDGEBOUNDS_Resize( EDGEBOUNDS*    edg,
                        int            size )
{
   edg->Nalloc = size;
   edg->bounds = (BOUND*) realloc( edg->bounds, sizeof(BOUND) * size );
}

/*
 *  FUNCTION:  EDGEBOUNDS_Reverse()
 *  SYNOPSIS:  Reverse order of edgebound list.
 */
void EDGEBOUNDS_Reverse( EDGEBOUNDS*   edg )
{
   BOUND tmp;
   /* iterate over half of list */
   for (int i = 0; i <= (edg->N / 2) - 1; ++i)
   {
      /* swap bound from front and back of list. */
      tmp = edg->bounds[i];
      edg->bounds[i] = edg->bounds[edg->N - i - 1];
      edg->bounds[edg->N - i - 1] = tmp;
   }
}

/*
 *  FUNCTION:  EDGEBOUNDS_Index()
 *  SYNOPSIS:  Index locations in EDGEBOUND list that start each unique BOUND id.
 *             Assumes <edg> is sorted.
 */
void EDGEBOUNDS_Index(EDGEBOUNDS *edg)
{
   int      i;       /* index of edgebound list */
   int      id_0;    /* current id */
   BOUND*   b_0;     /* pointer to current bound in list */

   VECTOR_INT_Reuse( edg->ids );
   VECTOR_INT_Reuse( edg->ids_idx );

   i     = 0;
   id_0  = edg->bounds[0].id;
   b_0   = &(edg->bounds[0]);

   VECTOR_INT_Pushback( edg->ids, id_0 );
   VECTOR_INT_Pushback( edg->ids_idx, i );

   for (i; i < edg->N; i++, b_0++)
   {
      if ( b_0->id != id_0 ) {
         id_0 = b_0->id;
         VECTOR_INT_Pushback( edg->ids, id_0 );
         VECTOR_INT_Pushback( edg->ids_idx, i );
      }
   }

   id_0 = b_0->id;
   VECTOR_INT_Pushback( edg->ids, id_0 );
   VECTOR_INT_Pushback( edg->ids_idx, edg->N );
}


/*
 *  FUNCTION:  EDGEBOUNDS_Sort()
 *  SYNOPSIS:  Sort <edg> bound list ascending: by id, lb, rb. Sorts in place.
 */
void EDGEBOUNDS_Sort( EDGEBOUNDS*   edg )
{
   EDGEBOUNDS_Sort_Sub( edg, 0, edg->N );
}

/*
 *  FUNCTION:  EDGEBOUNDS_Sort_Sub()
 *  SYNOPSIS:  Subcall to sort the edgebounds on range (beg, end]. Sorts in place.
 */
void EDGEBOUNDS_Sort_Sub(  EDGEBOUNDS*    edg,
                           int            beg,
                           int            end )
{
   /* TODO */
}

/*
 *  FUNCTION:  EDGEBOUNDS_Merge()
 *  SYNOPSIS:  Merge <edg> bound list by combining overlapping ranges.
 *             Assumes that <edg> is sorted.

 */
void EDGEBOUNDS_Merge( EDGEBOUNDS*   edg )
{
   /* if one or less edgebounds, its already sorted */
   if ( edg->N <= 1 ) return;

   BOUND*      b_0;           /* current bound */
   BOUND*      b_1;           /* previous bound */
   int         num_merges;    /* number of holes caused by  */

   b_0 = &(edg->bounds[1]);
   b_1 = &(edg->bounds[0]);

   /* iterate through all adjacent pairs of bounds. bc they are sorted, they can only merge with neighbors */
   for (int i = 1; i < edg->N; i++, b_0++, b_1++)
   {
      /* if on same row/diag and overlap, then merge */
      if ( b_1->id != b_0->id && b_1->rb >= b_0->lb ) {
         b_0->lb = b_1->lb;
         num_merges++;
      } 
      /* if not, then move previous bound to proper position */
      else 
      {
         /* to fill holes created by merging, move bound left by that amount */
         edg->bounds[ (i-1) - num_merges ] = *b_1;
      }
   }
   /* every merge step removes one bound from list */
   edg->N -= num_merges;
}

/*
 *  FUNCTION:  EDGEBOUNDS_Count_Cells()
 *  SYNOPSIS:  Count the number of cells covered by <edg>.
 *             Assumes <edg> is sorted and merged.
 */
int EDGEBOUNDS_Count_Cells( EDGEBOUNDS*   edg )
{
   BOUND* b_0;    /* current bound */
   int count;     /* total cells to current bound */

   count = 0;
   b_0   = &(edg->bounds[0]);

   for (int i = 0; i < edg->N; i++, b_0++) {
      count += (b_0->rb - b_0->lb);
   }

   return count;
}

/*
 *  FUNCTION: EDGEBOUNDS_Print()
 *  SYNOPSIS: Print EDGEBOUND object to file.
 */
void EDGEBOUNDS_Dump( EDGEBOUNDS*   edg,
                      FILE*         fp )
{
   /* test for bad file pointer */
   if (fp == NULL) {
      const char* obj_name = "EDGEBOUNDS";
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing %s.\n", obj_name);
      exit(EXIT_FAILURE);
      return;
   }

   char* edg_mode_text[] = { "EDGE_NONE", "EDGE_ANTIDIAG", "EDGE_ROW" };

   fprintf(fp, "# N: %d, Nalloc: %d\n", edg->N, edg->Nalloc);
   fprintf(fp, "# ORIENTATION: %s\n", edg_mode_text[edg->edg_mode] );
   fprintf(fp, "# Q=%d, T=%d\n", edg->Q, edg->T );
   for (unsigned int i = 0; i < edg->N; ++i)
   {
      BOUND bnd = edg->bounds[i];
      fprintf(fp, "[%d] ", i);
      fprintf(fp, "{ id: %d, lb: %d, rb: %d }\n", bnd.id, bnd.lb, bnd.rb);
   }
   fprintf(fp, "\n");
}

/*
 *  FUNCTION: EDGEBOUNDS_Dump()
 *  SYNOPSIS: Output EDGEBOUND object to file.
 */
void EDGEBOUNDS_Save(EDGEBOUNDS*  edg,
                     const char*  _filename_)
{
   FILE *fp;
   fp = fopen(_filename_, "w");
   EDGEBOUNDS_Dump(edg, fp);
   printf("Saved EDGEBOUNDS to: '%s;'\n", _filename_);
   fclose(fp);
}

/*
 *  FUNCTION: EDGEBOUNDS_Compare()
 *  SYNOPSIS: Compare two EDGEBOUNDS objects.  Return 0 if equal.
 */
int EDGEBOUNDS_Compare( EDGEBOUNDS*    edg_a,
                        EDGEBOUNDS*    edg_b )
{
   BOUND* bnd_a;
   BOUND* bnd_b;

   if ( edg_a->edg_mode != edg_b->edg_mode ) {
      printf("EDGEBOUNDS ARE NOT SAME ORIENTATION: edg_a = %d, edg_b = %d\n", 
         edg_a->edg_mode, edg_b->edg_mode );
   }

   for (int i = 0; i < edg_a->N; i++)
   {
      bnd_a = &(edg_a->bounds[i]);
      bnd_b = &(edg_b->bounds[i]);
      if ( bnd_a->id != bnd_b->id || bnd_a->lb != bnd_b->lb || bnd_a->rb != bnd_b->rb ) {
         #if DEBUG
         {
            printf("EDGEBOUND INEQUALITY at %d: edg_a = (%d,%d,%d), edg_b = (%d,%d,%d)\n", i, 
               bnd_a->id, bnd_a->lb, bnd_a->rb, 
               bnd_b->id, bnd_b->lb, bnd_b->rb);
         } 
         #endif
         return -1;
      }
   }
   return 0;
}

/*
 *  FUNCTION: EDGEBOUNDS_Count()
 *  SYNOPSIS: Count the number of cells in edgebound.
 */
int EDGEBOUNDS_Count(EDGEBOUNDS*    edg)
{
   int sum = 0;
   for (int i = 0; i < edg->N; i++)
   {
      sum += edg->bounds[i].rb - edg->bounds[i].lb;
   }
   return sum;
}

/*
 *  FUNCTION: EDGEBOUNDS_Validate()
 *  SYNOPSIS: Verifies that edgebound ranges don't go out-of-bounds of containing matrix dimensions.
 */
int EDGEBOUNDS_Validate(EDGEBOUNDS *edg)
{
   bool     valid = true;
   BOUND*   bnd   = NULL;
   int      T     = edg->T;
   int      Q     = edg->Q;

   /* if bounds are stored as rows */
   if ( edg->edg_mode == EDG_ROW ) 
   {
      for ( int i = 0; i < edg->N; i++ ) 
      {
         bnd = &(edg->bounds[i]);
         /* check that all values are non-negative */
         if ( bnd->id < 0 || bnd->lb < 0 || bnd->rb < 0 ) {
            valid = false;
            break;
         }
         /* check that left edge is less than right edge of bounds */
         if ( bnd->lb > bnd->rb ) {
            valid = false;
            break;
         }
         /* check row does not exceed row bounds */
         if ( bnd->id >= edg->Q ) {
            valid = false;
            break;
         }
         /* check that columns range is within matrix bounds  */
         if ( bnd->lb < 0 || bnd->rb > edg->T ) {
            valid = false;
            break;
         }
      }
   }
   /* if bounds are stored as anti-diagonals */
   else if ( edg->edg_mode == EDG_DIAG ) 
   {
      /* max number of anti-diagonals based on matrix dimension */
      int max_diag = (edg->Q+1) + (edg->T+1) - 1; 
      for ( int i = 0; i < edg->N; i++ ) 
      {
         bnd = &(edg->bounds[i]);
         /* valid min and max range of anti-diagonal */
         int min_rng = MAX( 0, bnd->id - (edg->T - 1) );
         int max_rng = MIN( bnd->id, edg->Q - 1 );

         /* verify all values are non-negative */
         if ( bnd->id < 0 || bnd->lb < 0 || bnd->rb < 0 ) {
            valid = false;
            break;
         }
         /* check that left edge is less than right edge of bounds */
         if ( bnd->lb > bnd->rb ) {
            valid = false;
            break;
         } 
         /* check that anti-diagonal is within bounds */
         if ( bnd->id >= max_diag ) {
            valid = false;
            break;
         }
         /* check that column range is within matrix bounds */
         if ( bnd->lb < min_rng || bnd->rb > max_rng ) {
            valid = false;
            break;
         }
      }
   }

   if ( valid == false ) {
      fprintf(stderr, "ERROR: Edgebounds are invalid for given matrix.\n");
      fprintf(stderr, "matrix dim: (%d,%d), bounds: %d,(%d,%d)\n", edg->Q, edg->T, bnd->id, bnd->lb, bnd->rb);
      exit(EXIT_FAILURE);
   }
}