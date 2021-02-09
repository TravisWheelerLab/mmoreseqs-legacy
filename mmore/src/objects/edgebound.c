/*******************************************************************************
 *  FILE:      edgebound.c
 *  PURPOSE:   EDGEBOUNDS Object
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

/* local imports */
#include "structs.h"
#include "../utilities/_utilities.h"
#include "basic/bound.h"
#include "_objects.h"

/* header */
#include "edgebound.h"

/*! FUNCTION:  EDGEBOUNDS_Create()
 *  SYNOPSIS:  Create new EDGEBOUNDS object and returns pointer.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Create()
{
   EDGEBOUNDS* edg      = NULL;
   const int   min_size = 8;

   edg = EDGEBOUNDS_Create_by_Size( min_size );

   return edg;
}

/*! FUNCTION:  EDGEBOUNDS_Create_by_Size()
 *  SYNOPSIS:  Create new EDGEBOUNDS object with chosen size and returns pointer.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Create_by_Size( const int size )
{
   EDGEBOUNDS* edg = NULL;
   edg = (EDGEBOUNDS*) ERROR_malloc( sizeof(EDGEBOUNDS) );

   edg->N         = 0;
   edg->Nalloc    = 0;

   edg->Q         = 0;
   edg->T         = 0;

   edg->ids       = VECTOR_INT_Create();
   edg->ids_idx   = VECTOR_INT_Create();
   edg->edg_mode  = EDG_NONE;

   edg->bounds    = NULL;

   EDGEBOUNDS_Resize(edg, size);

   return edg;
}

/*! FUNCTION: EDGEBOUNDS_Destroy()
 *  SYNOPSIS: Frees all memory from EDGEBOUNDS object.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Destroy( EDGEBOUNDS*  edg )
{
   if ( edg == NULL ) return edg;

   VECTOR_INT_Destroy( edg->ids );
   VECTOR_INT_Destroy( edg->ids_idx );

   ERROR_free( edg->bounds );
   edg->bounds = NULL;

   ERROR_free( edg );
   edg = NULL;

   return edg;
}

/*! FUNCTION: EDGEBOUNDS_Reuse()
 *  SYNOPSIS: Reuses EDGEBOUNDS by "clearing" edgebound list (does not realloc).
 */
void 
EDGEBOUNDS_Reuse(    EDGEBOUNDS*   edg, 
                     int           Q,
                     int           T )
{
   edg->N = 0;
   edg->Q = Q;
   edg->T = T;
}

/*! FUNCTION: EDGEBOUNDS_SetDim()
 *  SYNOPSIS: Set the dimensions of the embedding matrix.
 */
void 
EDGEBOUNDS_SetDim(  EDGEBOUNDS*   edg, 
                     int           Q,
                     int           T )
{
   edg->Q = Q;
   edg->T = T;
}

/*! FUNCTION: EDGEBOUNDS_Copy()
 *  SYNOPSIS: Create a deep copy of <edg_src> and store it in <edg_dest>.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Copy(  EDGEBOUNDS*          edg_dest,
                  const EDGEBOUNDS*    edg_src )
{
   BOUND*   bnd;

   /* if source and destination are the same, then do not copy. */
   if (edg_dest == edg_src) return edg_dest;

   /* if destination has not been created, do it now */
   if ( edg_dest == NULL ) {
      edg_dest = EDGEBOUNDS_Create_by_Size( edg_src->Nalloc );
   } 

   /* now copying begins */
   edg_dest->Q          = edg_src->Q;
   edg_dest->T          = edg_src->T;
   EDGEBOUNDS_Reuse( edg_dest, edg_src->Q, edg_src->T );
   edg_dest->edg_mode = edg_src->edg_mode;

   for ( int i = 0; i < edg_src->N; i++ ) 
   {
      bnd = &(edg_src->bounds[i]);
      EDGEBOUNDS_Pushback( edg_dest, bnd );
   }

   return edg_dest;
}

/*! FUNCTION: EDGEBOUNDS_Get()
 *  SYNOPSIS: Return pointer to BOUND at index <i>.
 */
inline
BOUND* 
EDGEBOUNDS_Get(   EDGEBOUNDS*   edg,
                  int           i )
{
   /* if debugging, do edgebound checks */
   // #if DEBUG 
   //    int N = EDGEBOUNDS_GetSize( edg );
   //    if ( i >= N || i < 0 ) {
   //       fprintf(stderr, "ERROR: EDGEBOUNDS Access Out-of-Bounds\n");
   //       fprintf(stderr, "dim: (%d/%d), access: (%d)\n", edg->N, edg->Nalloc, i);
   //       exit(EXIT_FAILURE);
   //    }
   // #endif

   return &(edg->bounds[i]);
}

/*! FUNCTION: EDGEBOUNDS_GetSize()
 *  SYNOPSIS: Get length of <edg>.
 */
inline
int 
EDGEBOUNDS_GetSize(  EDGEBOUNDS*   edg )
{
   return edg->N;
}

/*! FUNCTION: EDGEBOUNDS_SetSize()
 *  SYNOPSIS: Set length of <edg> to <size>.
 */
inline
STATUS_FLAG 
EDGEBOUNDS_SetSize(    EDGEBOUNDS*    edg,
                        int            size )
{
   edg->N = size;
}


/*! FUNCTION:  EDGEBOUNDS_Search()
 *  SYNOPSIS:  Binary search edgebounds for bound containing cell (q_0, t_0).
 *             Assumes edgebounds are sorted and merged.
 *  RETURN:    Return index of edgebound, or -1 if not contained.
 */
int 
EDGEBOUNDS_Search(   EDGEBOUNDS*    edg,     /* edgebounds  */
                     int            q_0,     /* row/diag index, position in query */
                     int            t_0 )    /* column index, position in target */
{
   int      N     = EDGEBOUNDS_GetSize( edg );
   int      idx   = N/2;
   int      cmp   = 0;
   BOUND*   bnd;

   /* binary search */
   for ( int i = N/4; i > 1; i = i/2 ) 
   {
      bnd = &(edg->bounds[idx]);

      /* check if on correct row */
      cmp = INT_Compare( q_0, bnd->id );
      if ( cmp > 0 ) 
      {
         idx -= i;
         continue;
      }
      else if ( cmp < 0 ) 
      {
         idx += i;
         continue;
      }
      else /* if ( cmp == 0 ) */ 
      {
         /* check if in correct column range */
         /* right of left bound? */
         cmp = INT_Compare( t_0, bnd->lb );
         if ( cmp < 0 ) {
            idx += i;
            continue;
         }
         /* left of right bound? */
         cmp = INT_Compare( t_0, bnd->rb - 1 );
         if ( cmp > 0 ) {
            idx -= i;
            continue;
         }
         /* if both, then it is inside range */
         return idx;
      }
   }

   return -1;
}

/*! FUNCTION: EDGEBOUNDS_Pushback()
 *  SYNOPSIS: Add BOUND to EDGEBOUNDS list.
 */
STATUS_FLAG 
EDGEBOUNDS_Pushback(    EDGEBOUNDS*  edg,
                        BOUND*       bnd )
{
   edg->bounds[edg->N] = *bnd;
   edg->N++;

   /* resize if necessary */
   if (edg->N >= edg->Nalloc - 1) {
      EDGEBOUNDS_Resize(edg, edg->Nalloc * 2);
   }

   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_Insert()
 *  SYNOPSIS: Insert/Overwrite bound into <i> index of Edgebound list.
 */
STATUS_FLAG 
EDGEBOUNDS_Insert(   EDGEBOUNDS*    edg,
                     int            i,
                     BOUND*         bnd )
{
   BOUND*   edg_bnd;

   edg_bnd = &(edg->bounds[i]);
   *edg_bnd = *bnd;

   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_Delete()
 *  SYNOPSIS: Delete BOUND at <i> index and fill from end of list <N-1>, then decrement list size.
 */
STATUS_FLAG 
EDGEBOUNDS_Delete(   EDGEBOUNDS*    edg,
                     int            i )
{
   int N = edg->N;
   edg->bounds[i] = edg->bounds[N-1];
   edg->N -= 1;

   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of BOUNDS allocated in EDGEBOUND object (does not downsize).
 */
STATUS_FLAG 
EDGEBOUNDS_GrowTo(   EDGEBOUNDS*    edg,
                     int            size )
{
   if ( edg->Nalloc < size ) {
      EDGEBOUNDS_Resize( edg, size );
   }
   
   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of BOUNDS allocated in EDGEBOUND object .
 */
void 
EDGEBOUNDS_Resize(   EDGEBOUNDS*    edg,
                     int            size )
{
   edg->Nalloc = size;
   edg->bounds = (BOUND*) realloc( edg->bounds, sizeof(BOUND) * size );
}

/*! FUNCTION:  EDGEBOUNDS_Reverse()
 *  SYNOPSIS:  Reverse order of edgebound list.
 */
STATUS_FLAG 
EDGEBOUNDS_Reverse( EDGEBOUNDS*   edg )
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

   return STATUS_SUCCESS;
}

/*! FUNCTION:  EDGEBOUNDS_Index()
 *  SYNOPSIS:  Index locations in EDGEBOUND list that start each unique BOUND id.
 *             Assumes <edg> is sorted.
 */
STATUS_FLAG 
EDGEBOUNDS_Index( EDGEBOUNDS*  edg )
{
   int      i_0;     /* index of edgebound list */
   int      id_0;    /* current id */
   BOUND*   b_0;     /* pointer to current bound in list */

   VECTOR_INT_Reuse( edg->ids );
   VECTOR_INT_Reuse( edg->ids_idx );

   i_0   = 0;
   id_0  = edg->bounds[0].id;
   b_0   = &(edg->bounds[0]);

   VECTOR_INT_Pushback( edg->ids, id_0 );
   VECTOR_INT_Pushback( edg->ids_idx, i_0 );

   for ( i_0; i_0 < edg->N; i_0++, b_0++)
   {
      if ( b_0->id != id_0 ) {
         id_0 = b_0->id;
         VECTOR_INT_Pushback( edg->ids, id_0 );
         VECTOR_INT_Pushback( edg->ids_idx, i_0 );
      }
   }

   id_0 = b_0->id;
   VECTOR_INT_Pushback( edg->ids, id_0 );
   VECTOR_INT_Pushback( edg->ids_idx, edg->N );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  EDGEBOUNDS_NxtRow()
 *  SYNOPSIS:  Iterating from front to back, gets the row index range <r_0> that are on <q_0> position, starting from <r_0e>.
 *             Skips over rows less than <q_0>.  Presumes that edgebounds is sorted and <r_e> precedes <q_0> row start.  
 */
STATUS_FLAG 
EDGEBOUNDS_NxtRow(   EDGEBOUNDS*          edg,     /* edgebounds */
                     int*                 r_0b,    /* row range begin */
                     int*                 r_0e,    /* row range end */
                     int                  id_0 )    /* query sequence position */
{
   int r_0;

   /* skip over rows before <q_0> */
   r_0 = *r_0e;
   while ( (r_0 < edg->N ) && (EDG_X(edg, r_0).id < id_0) ) {
      r_0++;
   }
   /* capture rows on <q_0> */
   *r_0b = r_0;
   while ( (r_0 < edg->N ) && (EDG_X(edg, r_0).id == id_0) ) {
      r_0++;
   }
   *r_0e = r_0;

   return STATUS_SUCCESS;
}

/*! FUNCTION:  EDGEBOUNDS_PrvRow()
 *  SYNOPSIS:  Iterating from back to front, gets the row index range <r_0> that are on <q_0> position, starting from .
 *             Skips over rows greater than <q_0>.  Presumes that edgebounds is sorted and <r_0e> precedes <q_0> row start.  
 */
STATUS_FLAG 
EDGEBOUNDS_PrvRow(   EDGEBOUNDS*          edg,     /* edgebounds */
                     int*                 r_0b,    /* row range begin */
                     int*                 r_0e,    /* row range end */
                     int                  id_0 )    /* query sequence position */
{
   int r_0;
   
   /* skip over rows before <q_0> */
   r_0 = *r_0e;
   while ( (r_0 > 0) && (EDG_X(edg, r_0).id > id_0) ) {
      r_0--;
   }
   /* capture rows on <q_0> */
   *r_0b = r_0;
   while ( (r_0 > 0) && (EDG_X(edg, r_0).id == id_0) ) {
      r_0--;
   }
   *r_0e = r_0;

   return STATUS_SUCCESS;
}

/*! FUNCTION:  EDGEBOUNDS_Sort()
 *  SYNOPSIS:  Sort <edg> bound list ascending: by id, lb, rb. Sorts in place.
 */
void 
EDGEBOUNDS_Sort( EDGEBOUNDS*   edg )
{
   int N = EDGEBOUNDS_GetSize(edg);
   EDGEBOUNDS_Sort_Sub( edg, 0, N );

   #if DEBUG 
   {
      for ( int i = 0; i < N-1; i++ ) 
      {
         BOUND* cur = &edg->bounds[i];
         BOUND* nxt = &edg->bounds[i+1];
         int cmp = BOUND_Compare( *cur, *nxt );
         if ( (cmp <= 0) == false ) {
            printf("ERROR: bad sort. %d, %d v %d: { %d %d %d } vs { %d %d %d }\n",
               cmp, i, i+1, cur->id, cur->lb, cur->rb, nxt->id, nxt->lb, nxt->rb );
         }
      }
   }
   #endif
}

/*! FUNCTION:  EDGEBOUNDS_Sort_Sub()
 *  SYNOPSIS:  Subcall to sort the edgebounds on range (beg, end]. Sorts in place.
 */
void 
EDGEBOUNDS_Sort_Sub(    EDGEBOUNDS*    edg,
                        int            beg,
                        int            end )
{
   const int begin_select_sort = 4;
   int N = EDGEBOUNDS_GetSize(edg);
   if (N <= 1) {
      return;
   } 

   int size = end - beg;
   /* run selection sort if below threshold */
   if ( size <= begin_select_sort ) {
      EDGEBOUNDS_Sort_Sub_Selectsort( edg, beg, end );
   } 
   /* otherwise run quiksort */
   else {
      EDGEBOUNDS_Sort_Sub_Quicksort( edg, beg, end );
   }
}

/*! FUNCTION:  EDGEBOUNDS_Sort_Sub_Selectsort()
 *  SYNOPSIS:  Selection Sorts subarray of <edg> data in ascending order on range (beg,end].  
 */
void 
EDGEBOUNDS_Sort_Sub_Selectsort(  EDGEBOUNDS*    edg,
                                 int            beg,
                                 int            end )
{
   for ( int i = beg; i < end; i++ ) 
   {
      /* initial minimum value found */
      int   min_idx = i;
      BOUND min_val = edg->bounds[i];
      for ( int j = i+1; j < end; j++ ) {
         /* if new minimum found, update value and index */
         int cmp = BOUND_Compare( min_val, edg->bounds[j] );
         if ( cmp > 0 ) {
            min_idx = j;
            min_val = edg->bounds[j];
         }
      }
      /* swap new minimum to left-most position */
      EDGEBOUNDS_Swap( edg, i, min_idx );
   }
}

/*! FUNCTION:  EDGEBOUNDS_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Quick Sorts subarray of <edg> data in ascending order on range (beg,end].  
 */
void 
EDGEBOUNDS_Sort_Sub_Quicksort(   EDGEBOUNDS*    edg,
                                 int            beg,
                                 int            end )
{
   /* partition pointers */
   int   r_idx    = end - 1;
   int   l_idx    = beg + 1;
   BOUND*  rhs    = &(edg->bounds[beg + 1]);
   BOUND*  lhs    = &(edg->bounds[end - 1]);

   /* select random pivot value */
   int      range       = end - beg;
   int      pivot_idx   = RNG_INT_Range( beg, end );
   BOUND    pivot_val   = edg->bounds[pivot_idx];
   EDGEBOUNDS_Swap( edg, pivot_idx, beg );

   /* partition on pivot */
   while ( l_idx <= r_idx )
   {
      /* find next right partition element that is less than pivot element */
      while ( (l_idx <= r_idx) && (BOUND_Compare( pivot_val, edg->bounds[r_idx] ) < 0) ) {
         r_idx--;
      }
      /* find next left partition element that is greater than pivot element */
      while ( (l_idx <= r_idx) && (BOUND_Compare( pivot_val, edg->bounds[l_idx] ) >= 0) ) {
         l_idx++;
      }
      /* if left and right index have not crossed, then swap elements */
      if ( l_idx <= r_idx ) {
         EDGEBOUNDS_Swap( edg, l_idx, r_idx );
      }
   }
   /* move partition element to barrier between left and right index */
   EDGEBOUNDS_Swap( edg, beg, r_idx );
   /* sort both partitions (omit partition element) */
   EDGEBOUNDS_Sort_Sub( edg, beg, r_idx );
   EDGEBOUNDS_Sort_Sub( edg, r_idx+1, end );
}

/*! FUNCTION:  EDGEBOUNDS_Swap()
 *  SYNOPSIS:  Swaps the values of <edg> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
inline
void 
EDGEBOUNDS_Swap(  EDGEBOUNDS*    edg,
                  int            i,
                  int            j )
{
   BOUND swap = edg->bounds[i];
   edg->bounds[i] = edg->bounds[j];
   edg->bounds[j] = swap;
}

/*! FUNCTION:  EDGEBOUNDS_Merge()
 *  SYNOPSIS:  Merge <edg>'s bound list by combining overlapping ranges. In-place.
 *             Assumes that <edg> is sorted.
 */
void 
EDGEBOUNDS_Merge( EDGEBOUNDS*    edg )
{
   /* if one or less edgebounds, its already merged */
   int N = EDGEBOUNDS_GetSize(edg);
   if ( N <= 1 ) return;

   EDGEBOUNDS_Merge_Sub( edg, 0, N );
}

/*! FUNCTION:  EDGEBOUNDS_Merge_Sub()
 *  SYNOPSIS:  Merge <edg>'s subarray of bound list by combining overlapping ranges. In-place.
 *             Assumes that <edg> is sorted.
 */
void 
EDGEBOUNDS_Merge_Sub(   EDGEBOUNDS*    edg,
                        int            beg,
                        int            end )
{
   BOUND*      b_0;           /* current bound */
   BOUND*      b_1;           /* previous bound */
   BOUND*      b_x;           /* new bound to fill */
   int         num_merges;    /* number of holes caused by  */
   int         num_fills;     /* number of bounds filled in by merge */

   /* initialize number of merges and fills */
   num_merges  = 0;
   num_fills   = 0;

   /* iterate through all adjacent pairs of bounds. bc they are sorted, they can only merge with neighbors */
   for ( int i = beg+1; i < end; i++ )
   {
      /* get next two adjacent bounds */
      b_0  = EDGEBOUNDS_Get( edg, i );
      b_1  = EDGEBOUNDS_Get( edg, i-1 );
      /* if adjacent bounds are on same row/diag and ranges overlap, then merge */
      if ( b_1->id == b_0->id && b_1->rb >= b_0->lb ) 
      {
         /* update current bounds by adding previous bounds to it */
         b_0->lb = b_1->lb;
         b_0->rb = MAX( b_0->rb, b_1->rb );
         num_merges++;
      } 
      /* if not, then move previous bound to proper position */
      else 
      {
         /* to fill holes created by merging, move bound by that amount */
         b_x   = EDGEBOUNDS_Get( edg, num_fills );
         /* add previous edgebound last unmerged to next open place in list */
         *b_x  = *b_1;
         /* number of fills determines resulting edgebound size */
         num_fills++;
      }
   }
   /* add final edgebound to list */
   b_0   = EDGEBOUNDS_Get( edg, end-1 );
   /* to fill holes created by merging, move bound by that amount */
   b_x   = EDGEBOUNDS_Get( edg, num_fills );
   /* add previous edgebound last unmerged to next open place in list */
   *b_x  = *b_0;
   /* number of fills determines resulting edgebound size */
   num_fills++;

   /* every merge step removes one bound from list */
   EDGEBOUNDS_SetSize( edg, num_fills );
}

/*! FUNCTION: EDGEBOUNDS_Print()
 *  SYNOPSIS: Print EDGEBOUND object to file.
 */
void 
EDGEBOUNDS_Dump(  EDGEBOUNDS*    edg,
                  FILE*          fp )
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

   int N = EDGEBOUNDS_GetSize( edg );
   for (int i = 0; i < edg->N; ++i)
   {
      BOUND bnd = edg->bounds[i];
      fprintf(fp, "[%d] ", i);
      fprintf(fp, "{ id: %d, lb: %d, rb: %d }\n", bnd.id, bnd.lb, bnd.rb);
   }
   fprintf(fp, "\n");
}

/*! FUNCTION: EDGEBOUNDS_Print()
 *  SYNOPSIS: Print EDGEBOUND object to file.
 */
void 
EDGEBOUNDS_Sub_Dump(    EDGEBOUNDS*    edg,
                        FILE*          fp,
                        int            beg, 
                        int            end )
{
   /* test for bad file pointer */
   if (fp == NULL) {
      const char* obj_name = "EDGEBOUNDS";
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing %s.\n", obj_name);
      exit(EXIT_FAILURE);
      return;
   }

   for (unsigned int i = beg; i < end; ++i)
   {
      BOUND bnd = edg->bounds[i];
      fprintf(fp, "[%d] ", i);
      fprintf(fp, "{ id: %d, lb: %d, rb: %d }\n", bnd.id, bnd.lb, bnd.rb);
   }
   printf("\n");
}

/*! FUNCTION: EDGEBOUNDS_Dump()
 *  SYNOPSIS: Output EDGEBOUND object to file.
 */
void 
EDGEBOUNDS_Save(  EDGEBOUNDS*   edg,
                  const char*   filename )
{
   FILE *fp;
   fp = fopen(filename, "w");
   EDGEBOUNDS_Dump(edg, fp);
   printf("Saved EDGEBOUNDS to: '%s'\n", filename);
   fclose(fp);
}

/*! FUNCTION: EDGEBOUNDS_Compare()
 *  SYNOPSIS: Compare two EDGEBOUNDS objects.  Return 0 if equal.
 */
int 
EDGEBOUNDS_Compare(  EDGEBOUNDS*    edg_a,
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

/*! FUNCTION: EDGEBOUNDS_Count()
 *  SYNOPSIS: Count the number of cells in edgebound.
 */
int 
EDGEBOUNDS_Count( EDGEBOUNDS*  edg )
{
   int sum = 0;
   for (int i = 0; i < edg->N; i++)
   {
      sum += edg->bounds[i].rb - edg->bounds[i].lb;
   }
   return sum;
}

/*! FUNCTION: EDGEBOUNDS_Validate()
 *  SYNOPSIS: Verifies that edgebound ranges don't go out-of-bounds of containing matrix dimensions.
 *            For testing. 
 */
int 
EDGEBOUNDS_Validate( EDGEBOUNDS*  edg )
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

/*! FUNCTION: EDGEBOUNDS_Cover_Matrix()
 *  SYNOPSIS: Creates an edgebounds that covers every cell in DP Matrix with dimensions {Q x T}.
 *            For testing.
 */
void
EDGEBOUNDS_Cover_Matrix(   EDGEBOUNDS*    edg, 
                           int            Q,
                           int            T )
{
   EDGEBOUNDS_Clear(edg);
   for (int q_0 = 0; q_0 <= Q; q_0++) {
      EDGEBOUNDS_Pushback(edg, &(BOUND){ q_0, 0, T+1 });
   }
}

/*! FUNCTION:  EDGEBOUNDS_Cover_Range()
 *  SYNOPSIS:  Creates edgebound space that fills square with Q_range in Query and T_range in Target.
 *             For testing.
 */
STATUS_FLAG
EDGEBOUNDS_Cover_Range(    EDGEBOUNDS*    edg,
                           RANGE          Q_range,
                           RANGE          T_range )
{
   int Q_beg, Q_end;
   int T_beg, T_end;

   Q_beg = Q_range.beg;
   Q_end = Q_range.end;
   T_beg = T_range.beg;
   T_end = T_range.end;

   EDGEBOUNDS_Clear(edg);
   for (int q_0 = Q_beg; q_0 < Q_end; q_0++) {
      EDGEBOUNDS_Pushback(edg, &(BOUND){q_0, T_beg, T_end});
   }
}


/*! FUNCTION: EDGEBOUNDS_Clear()
 *  SYNOPSIS: Remove all BOUNDS from EDGEBOUND list.
 */
void 
EDGEBOUNDS_Clear( EDGEBOUNDS* edg ) 
{
   edg->N = 0;
}

/*! FUNCTION: EDGEBOUNDS_Find_BoundingBox()
 *  SYNOPSIS: Find the min/max range of values contained in the edgebounds.
 *            Assumes edgebounds have been sorted.
 */
int 
EDGEBOUNDS_Find_BoundingBox(  EDGEBOUNDS*   edg,
                              RANGE*        Q_range,
                              RANGE*        T_range )
{
   int Q, T;
   Q = edg->Q;
   T = edg->T;

   /* find the target and query range */
   Q_range->beg = Q + 1;
   Q_range->end = 0;
   T_range->beg = T + 1;
   T_range->end = 0;

   /* if list is empty, return */
   if (edg->N == 0) {
      return STATUS_SUCCESS;
   }

   /* create bounding box */
   for (int i = 0; i < edg->N; i++) {
      if ( T_range->beg > edg->bounds[i].lb ) {
         T_range->beg = edg->bounds[i].lb;
      }
      if ( T_range->end < edg->bounds[i].rb ) {
         T_range->end = edg->bounds[i].rb;
      }
   }
   Q_range->beg = edg->bounds[0].id;
   Q_range->end = edg->bounds[edg->N - 1].id;
   /* edge checks */
   T_range->beg = MAX(T_range->beg, 0);
   T_range->end = MIN(T_range->end, T + 1);
   Q_range->beg = MAX(Q_range->beg, 0);
   Q_range->end = MIN(Q_range->end, Q + 1);

   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_SetDomain()
 *  SYNOPSIS: Build an EDGEBOUND <edg_out> from QxT EDGEBOUNDS <edg_in>
 *            and constraining the query range to <in_dom_range> => <q_beg, q_end>.
 *            Simply shifts all query indexes by <q_beg>.
 *            NOTE: Assumes <edg_in> has been sorted.
 */
int 
EDGEBOUNDS_SetDomain(  EDGEBOUNDS*    edg_in,
                        EDGEBOUNDS*    edg_out,
                        RANGE          Q_range )
{
   BOUND bnd;
   int Q_len;
   Q_len = Q_range.end - Q_range.beg + 1; 

   /* clears old data */
   EDGEBOUNDS_Reuse( edg_out, Q_len, edg_in->T );
   /* shift all query indexes to set <q_beg> to 0. */
   for (int i = 0; i < edg_in->N; i++)
   {
      bnd = EDG_X(edg_in, i);
      bnd.id -= Q_range.beg;
      EDGEBOUNDS_Pushback( edg_out, &bnd );
   }

   return STATUS_SUCCESS;
}
