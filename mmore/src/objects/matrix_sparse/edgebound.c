/*******************************************************************************
 *  FILE:      edgebound.c
 *  PURPOSE:   EDGEBOUNDS Object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *  TODO:
 *    -  Need to replace BOUNDS* <bounds> with a VECTOR_BOUND*.
 *       Right now, with the way things are referenced in other algorithms it is too much of a chore.
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
#include "edgebound.h"

/* index padding */
const int index_pad = 1;

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
   EDGEBOUNDS* edg;
   edg = ERROR_malloc( sizeof(EDGEBOUNDS) );

   edg->Q         = 0;
   edg->T         = 0;
   /* index */
   edg->id_index  = VECTOR_INT_Create();
   edg->edg_mode  = EDG_NONE;
   /* data */
   edg->bounds      = VECTOR_BOUND_Create();
   edg->is_indexed   = false;
   edg->is_sorted    = false;
   edg->is_merged    = false;

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

   edg->id_index  = VECTOR_INT_Destroy( edg->id_index );
   edg->bounds    = VECTOR_BOUND_Destroy( edg->bounds );
   edg            = ERROR_free( edg );

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
   edg->Q            = Q;
   edg->T            = T;
   edg->is_indexed   = false;
   edg->is_merged    = false;
   edg->is_sorted    = false;
   VECTOR_INT_Reuse( edg->id_index );
   VECTOR_BOUND_Reuse( edg->bounds );
}

/*! FUNCTION: EDGEBOUNDS_Clear()
 *  SYNOPSIS: Reuses EDGEBOUNDS by "clearing" edgebound list (does not realloc).
 */
void 
EDGEBOUNDS_Clear(    EDGEBOUNDS*   edg )
{
   edg->Q            = 0;
   edg->T            = 0;
   edg->is_indexed   = false;
   edg->is_merged    = false;
   edg->is_sorted    = false;
   VECTOR_INT_Reuse( edg->id_index );
   VECTOR_BOUND_Reuse( edg->bounds );
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
   BOUND   bnd;

   /* if source and destination are the same, then do not copy. */
   if (edg_dest == edg_src) {
      return edg_dest;
   } 
   /* if destination has not been created, do it now */
   if ( edg_dest == NULL ) {
      size_t size    = EDGEBOUNDS_GetSize( edg_src );
      edg_dest       = EDGEBOUNDS_Create_by_Size( size );
   } 

   /* now copying begins */
   edg_dest->Q          = edg_src->Q;
   edg_dest->T          = edg_src->T;
   edg_dest->Q_range    = edg_src->Q_range;
   edg_dest->T_range    = edg_src->T_range;
   edg_dest->is_indexed = edg_src->is_indexed;
   edg_dest->is_merged  = edg_src->is_merged;
   edg_dest->is_sorted  = edg_src->is_sorted;
   edg_dest->edg_mode   = edg_src->edg_mode;

   EDGEBOUNDS_Reuse(    edg_dest, edg_src->Q,   edg_src->T );
   VECTOR_BOUND_Copy(   edg_dest->bounds,       edg_src->bounds );
   VECTOR_INT_Copy(     edg_dest->id_index,     edg_src->id_index );

   return edg_dest;
}

/*! FUNCTION: EDGEBOUNDS_Set()
 *  SYNOPSIS: Gets BOUND at index <i>.
 */
inline
STATUS_FLAG
EDGEBOUNDS_Set(   EDGEBOUNDS*    edg,
                  int            i,
                  BOUND          val )
{
   *EDGEBOUNDS_GetX( edg, i ) = val;
   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_GetX()
 *  SYNOPSIS: Get reference to BOUND at index <i>.
 */
inline
BOUND* 
EDGEBOUNDS_GetX(  const EDGEBOUNDS*    edg,
                  int                  i )
{
   BOUND* bnd = VECTOR_BOUND_GetX( edg->bounds, i );
   return bnd;
}

/*! FUNCTION: EDGEBOUNDS_Get()
 *  SYNOPSIS: Gets BOUND at index <i>.
 */
inline
BOUND
EDGEBOUNDS_Get(   const EDGEBOUNDS*    edg,
                  int                  i )
{
   BOUND bnd = VECTOR_BOUND_Get( edg->bounds, i );
   return bnd;
}

/*! FUNCTION:  EDGEBOUNDS_GetNumberBounds_byRow()
 *  SYNOPSIS:  Gets the number of <bounds> that are on row <id_0>.
 *             Caller must have already run _Index().
 */
inline
int
EDGEBOUNDS_GetNumberBounds_byRow(   const EDGEBOUNDS*   edg,
                                    const int           id_0 )
{
   int id_offset;
   int num_bounds;
   /* if outside of range, then zero. */
   bool in_range = IS_IN_RANGE( edg->Q_range.beg, edg->Q_range.end, id_0 );
   if ( in_range == false ) {
      num_bounds = 0;
      return num_bounds;
   }
   /* if in range, table lookup */
   id_offset   = id_0 - edg->Q_range.beg;
   num_bounds  = VEC_X( edg->id_index, id_offset + 1 ) - VEC_X( edg->id_index, id_offset );
   return num_bounds;
}

/*! FUNCTION:  EDGEBOUNDS_GetIndex_byRow()
 *  SYNOPSIS:  Gets the index into <bounds> at the start of row <id_0>.
 *             Caller must have already run _Index().
 *             Returns -1 if row not in <edg> range.
 */
inline
int
EDGEBOUNDS_GetIndex_byRow(    const EDGEBOUNDS*   edg,
                              const int           id_0 )
{
   int id_offset;
   int index;
   int N;
   
   id_offset   = id_0 - edg->Q_range.beg + index_pad;
   N           = VECTOR_INT_GetSize( edg->id_index );

   /* if below range, then zero. */
   if ( id_offset < 0 + index_pad ) {
      id_offset = 0 + index_pad;
   }
   /* if above range, then max */
   elif ( id_offset >= (N-1) ) {
      id_offset = (N-1);
   }
   /* if in range, table lookup */
   index = VECTOR_INT_Get( edg->id_index, id_offset );
   return index;
}

/*! FUNCTION:  EDGEBOUNDS_GetIndex_byRow_Fwd()
 *  SYNOPSIS:  Gets the index into <bounds> at the start of row <id_0>.
 *             Edgechecks specific to forward direction.
 *             Caller must have already run _Index().
 *             Returns -1 if row not in <edg> range.
 */
inline
int
EDGEBOUNDS_GetIndex_byRow_Fwd(   const EDGEBOUNDS*   edg,
                                 const int           id_0 )
{
   int id_offset;
   int index;
   int N;
   
   id_offset   = id_0 - edg->Q_range.beg + index_pad;
   N           = VECTOR_INT_GetSize( edg->id_index );

   /* if below range, then zero. */
   if ( id_offset < 0 + index_pad ) {
      id_offset = 0 + index_pad;
   }
   /* if above range, then max */
   if ( id_offset >= (N-1) ) {
      id_offset = (N-1);
   }
   /* if in range, table lookup */
   index = VECTOR_INT_Get( edg->id_index, id_offset );
   return index;
}

/*! FUNCTION:  EDGEBOUNDS_GetIndex_byRow_Bck()
 *  SYNOPSIS:  Gets the index into <bounds> at the start of row <id_0>.
 *             Edgechecks specific to backward direction.
 *             Caller must have already run _Index().
 *             Returns -1 if row not in <edg> range.
 */
inline
int
EDGEBOUNDS_GetIndex_byRow_Bck(   const EDGEBOUNDS*   edg,
                                 const int           id_0 )
{
   int id_offset;
   int index;
   int N;
   
   id_offset   = id_0 - edg->Q_range.beg + index_pad;
   N           = VECTOR_INT_GetSize( edg->id_index );

   /* if below range, then zero. */
   if ( id_offset < 0 + index_pad ) {
      id_offset = 0 + index_pad;
   }
   /* if above range, then max */
   if ( id_offset >= (N-1) ) {
      id_offset = (N-1);
   }
   /* if in range, table lookup */
   index = VEC_X( edg->id_index, id_offset ) - 1;
   return index;
}

/*! FUNCTION:  EDGEBOUNDS_GetBoundRange_byRow()
 *  SYNOPSIS:  Gets the index range <r_0b,r_0e> of bound row <id_0>.
 *             Returns pointer to head of row in <bounds>.
 *             Caller must have already run _Index().
 */
STATUS_FLAG
EDGEBOUNDS_GetIndexRange_byRow(  EDGEBOUNDS*    edg,
                                 int            id_0,
                                 RANGE*         index_out )
{
   int      id_offset;
   RANGE    index;
   int      num_bounds;
   BOUND*   row_head;
   
   /* otherwise, we look it up in the index table */
   index.beg = EDGEBOUNDS_GetIndex_byRow_Fwd( edg, id_0 );
   index.end = EDGEBOUNDS_GetIndex_byRow_Fwd( edg, id_0 + 1 );
   if ( index_out != NULL ) {
      *index_out = index;
   }
   return STATUS_SUCCESS;
}

/*! FUNCTION:  EDGEBOUNDS_GetBounds_byRow()
 *  SYNOPSIS:  Gets <i_0>th bound of row <id_0>.
 *             Caller must have already run _Index().
 */
inline
BOUND
EDGEBOUNDS_GetRow(   EDGEBOUNDS*    edg,
                     int            id_0,
                     int            i_0 )
{
   /* edgebounds check: if outside range or range set is less than <i_0>, then there is no edgebound to get */
   #if SAFE
   {
      
   }
   #endif
   
   /* offset id by Q_range */
   int      id_offset   = id_0 - edg->Q_range.beg;
   /* get index to start of row */
   int      index       = VEC_X( edg->id_index, id_offset );
   /* get bound from data */
   BOUND    bnd         = EDG_X( edg, index + i_0 );
   return bnd;
}

/*! FUNCTION: EDGEBOUNDS_GetSize()
 *  SYNOPSIS: Get length of <edg>.
 */
inline
size_t 
EDGEBOUNDS_GetSize(  const EDGEBOUNDS*   edg )
{
   size_t size = VECTOR_BOUND_GetSize( edg->bounds );
   return size;
}

/*! FUNCTION: EDGEBOUNDS_SetSize()
 *  SYNOPSIS: Set length of <edg> to <size>.
 */
inline
STATUS_FLAG 
EDGEBOUNDS_SetSize(     EDGEBOUNDS*    edg,
                        size_t         size )
{
   VECTOR_BOUND_SetSize( edg->bounds, size );
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
      bnd = EDGEBOUNDS_GetX( edg, idx );

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
                        BOUND        bnd )
{
   STATUS_FLAG status;
   status = VECTOR_BOUND_Pushback( edg->bounds, bnd );
   return status;
}

/*! FUNCTION:  EDGEBOUNDS_Insert()
 *  SYNOPSIS:  Delete BOUND at <i> index and fill from end of list <N-1>, then decrement list size.
 *             WARNING: This will break a sorted list.
 */
STATUS_FLAG 
EDGEBOUNDS_Insert(   EDGEBOUNDS*    edg,
                     int            i )
{
   // VECTOR_BOUND_Insert( edg->bounds, i );
   return STATUS_SUCCESS;
}

/*! FUNCTION:  EDGEBOUNDS_Delete()
 *  SYNOPSIS:  Delete BOUND at <i> index and fill from end of list <N-1>, then decrement list size.
 *             WARNING: This will break a sorted list.
 */
STATUS_FLAG 
EDGEBOUNDS_Delete(   EDGEBOUNDS*    edg,
                     int            i )
{
   VECTOR_BOUND_Delete( edg->bounds, i );
   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of BOUNDS allocated in EDGEBOUND object (does not downsize).
 */
STATUS_FLAG 
EDGEBOUNDS_GrowTo(   EDGEBOUNDS*    edg,
                     size_t         size )
{
   VECTOR_BOUND_GrowTo( edg->bounds, size );
   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of BOUNDS allocated in EDGEBOUND object .
 */
void 
EDGEBOUNDS_Resize(   EDGEBOUNDS*    edg,
                     size_t         size )
{
   VECTOR_BOUND_Resize( edg->bounds, size );
}

/*! FUNCTION:  EDGEBOUNDS_Reverse()
 *  SYNOPSIS:  Reverse order of edgebound list.
 */
STATUS_FLAG 
EDGEBOUNDS_Reverse( EDGEBOUNDS*   edg )
{
   VECTOR_BOUND_Reverse( edg->bounds );
   return STATUS_SUCCESS;
}

/*! FUNCTION:  EDGEBOUNDS_Index()
 *  SYNOPSIS:  Index locations in EDGEBOUND list that start each unique BOUND id.
 *             Assumes <edg> is sorted.
 */
STATUS_FLAG 
EDGEBOUNDS_Index( EDGEBOUNDS*  edg )
{
   RANGE    Q_range;
   int      Q_size;
   int      N;

   /* number of entries in edgebounds map */
   N = EDGEBOUNDS_GetSize( edg );
   /* find min/max range of Q and T */
   EDGEBOUNDS_Find_BoundingBox( edg, NULL, NULL );
   /* we need an index to span all reachable Q values */
   Q_range  = edg->Q_range;
   Q_size   = Q_range.end - Q_range.beg + 1;
   // printf("Len: %d, Q_range: (%d,%d), Q_size: %d\n", N, Q_range.beg, Q_range.end, Q_size);

   /* clear old data */
   VECTOR_INT_Reuse( edg->id_index );
   /* set size to cover Q_range (plus leftside/rightside pad) */
   VECTOR_INT_SetSize( edg->id_index, Q_size + (2 * index_pad) );

   /* initial id size is minimum Q value */
   int      id_0;
   int      id_offset;
   BOUND    b_0;

   /* initial pad index for look backward safety */
   for (int i = 1; i < index_pad + 1; i++) {
      id_0 = Q_range.beg - i;
      id_offset = id_0 - Q_range.beg + index_pad;
      VEC_X( edg->id_index, id_offset ) = -1;
   }

   id_0 = Q_range.beg;
   /* go through edgebound list and map <id_0>'s to position in <bounds>. */
   for ( int i_0 = 0; i_0 < N; i_0++)
   {
      b_0 = EDGEBOUNDS_Get( edg, i_0 );
      /* assign <id_0> index to its first instance in <bounds>. */
      /* if there are no <bounds> with <id_0>, then it is assigned to the edgebound proceed where it would be */
      while ( b_0.id >= id_0 ) 
      {
         id_offset = id_0 - Q_range.beg + index_pad;
         VEC_X( edg->id_index, id_offset ) = i_0;
         id_0 += 1;
      }
   }
   
   /* final pad index for look forward safety */
   for (int i = 1; i < index_pad + 1; i++) {
      id_0        = Q_range.end + i;
      id_offset   = id_0 - Q_range.beg + index_pad;
      VEC_X( edg->id_index, id_offset ) = N; /* should be (Q_size - 1) index */
   }

   edg->is_indexed = true;

   return STATUS_SUCCESS;
}

/*! FUNCTION:  EDGEBOUNDS_NxtRow()
 *  SYNOPSIS:  Iterating from front to back, gets the row index range <r_0> that are on <q_0> position, starting from <r_0e>.
 *             Skips over rows less than <q_0>.  Presumes that edgebounds is sorted and <r_e> precedes <q_0> row start.  
 */
inline
STATUS_FLAG 
EDGEBOUNDS_NxtRow(   EDGEBOUNDS* restrict       edg,     /* edgebounds */
                     int* restrict              r_0b,    /* row range begin */
                     int* restrict              r_0e,    /* row range end */
                     int                        id_0 )   /* query sequence position */
{
   int   N = EDGEBOUNDS_GetSize( edg );
   int   r_0b_old = *r_0b; 
   int   r_0e_old = *r_0e;
   int   r_0;

   /* skip over rows before <q_0> */
   r_0 = *r_0e;
   while ( (r_0 < N ) && (EDG_X(edg, r_0).id < id_0) ) {
      r_0++;
   }
   /* capture rows on <q_0> */
   *r_0b = r_0;
   while ( (r_0 < N ) && (EDG_X(edg, r_0).id == id_0) ) {
      r_0++;
   }
   *r_0e = r_0;

   // printf("[%d] (%d,%d) -> (%d,%d)\n", id_0, r_0b_old, r_0e_old, *r_0b, *r_0e );

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
   VECTOR_BOUND_Sort( edg->bounds );
   edg->is_sorted = true;
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
   VECTOR_BOUND_Swap( edg->bounds, i, j );
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
      b_0  = EDGEBOUNDS_GetX( edg, i );
      b_1  = EDGEBOUNDS_GetX( edg, i-1 );
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
         b_x   = EDGEBOUNDS_GetX( edg, num_fills );
         /* add previous edgebound last unmerged to next open place in list */
         *b_x  = *b_1;
         /* number of fills determines resulting edgebound size */
         num_fills++;
      }
   }
   /* add final edgebound to list */
   b_0   = EDGEBOUNDS_GetX( edg, end-1 );
   /* to fill holes created by merging, move bound by that amount */
   b_x   = EDGEBOUNDS_GetX( edg, num_fills );
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
EDGEBOUNDS_Dump(  const EDGEBOUNDS*    edg,
                  FILE*                fp )
{
   /* test for bad file pointer */
   if (fp == NULL) {
      const char* obj_name = "EDGEBOUNDS";
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing %s.\n", obj_name);
      ERRORCHECK_exit(EXIT_FAILURE);
      return;
   }

   int   N        = EDGEBOUNDS_GetSize( edg );
   int   Nalloc   = 0; /* TODO: add _getSizeAlloc() function */
   char* edg_mode_text[]   = { "EDGE_NONE", "EDGE_ANTIDIAG", "EDGE_ROW" };

   // fprintf(fp, "# N: %d, Nalloc: %d\n", N, edg->Nalloc );
   fprintf(fp, "# ORIENTATION: %s\n", edg_mode_text[edg->edg_mode] );
   fprintf(fp, "# Q=%d, T=%d\n", edg->Q, edg->T );

   for (int i = 0; i < N; ++i)
   {
      BOUND bnd = EDGEBOUNDS_Get( edg, i );
      fprintf(fp, "[%d] ", i);
      fprintf(fp, "{ id: %d, lb: %d, rb: %d }\n", bnd.id, bnd.lb, bnd.rb);
   }
   fprintf(fp, "\n");
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
   if ( edg_a->edg_mode != edg_b->edg_mode ) {
      printf("EDGEBOUNDS ARE NOT SAME ORIENTATION: edg_a = %d, edg_b = %d\n", 
         edg_a->edg_mode, edg_b->edg_mode );
   }

   int N_a  = EDGEBOUNDS_GetSize( edg_a );
   int N_b  = EDGEBOUNDS_GetSize( edg_b );
   if ( N_a - N_b != 0 ) {
      return -1;
   }

   for (int i = 0; i < N_b; i++)
   {
      BOUND bnd_a    = EDGEBOUNDS_Get( edg_a, i );
      BOUND bnd_b    = EDGEBOUNDS_Get( edg_b, i );
      if ( bnd_a.id != bnd_b.id || bnd_a.lb != bnd_b.lb || bnd_a.rb != bnd_b.rb ) {
         #if DEBUG
         {
            printf("EDGEBOUND INEQUALITY at %d: edg_a = (%d,%d,%d), edg_b = (%d,%d,%d)\n", i, 
               bnd_a.id, bnd_a.lb, bnd_a.rb, 
               bnd_b.id, bnd_b.lb, bnd_b.rb);
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
   int sum  = 0;
   int N    = EDGEBOUNDS_GetSize( edg );
   for (int i = 0; i < N; i++)
   {
      sum += EDGEBOUNDS_Get( edg, i ).rb - EDGEBOUNDS_Get( edg, i ).lb;
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
   bool     valid    = true;
   BOUND*   bnd      = NULL;
   int      T        = edg->T;
   int      Q        = edg->Q;
   int      N        = EDGEBOUNDS_GetSize( edg );

   /* if bounds are stored as rows */
   if ( edg->edg_mode == EDG_ROW ) 
   {
      for ( int i = 0; i < N; i++ ) 
      {
         bnd = EDGEBOUNDS_GetX( edg, i );
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
      for ( int i = 0; i < N; i++ ) 
      {
         bnd = EDGEBOUNDS_GetX( edg, i );
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
      ERRORCHECK_exit(EXIT_FAILURE);
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
   EDGEBOUNDS_Reuse(edg, Q, T);
   for (int q_0 = 0; q_0 <= Q; q_0++) 
   {
      BOUND bnd = (BOUND){ q_0, 0, T+1 };
      EDGEBOUNDS_Pushback(edg, bnd);
   }
   EDGEBOUNDS_Index(edg);
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
      BOUND bnd = (BOUND){ q_0, T_beg, T_end };
      EDGEBOUNDS_Pushback(edg, bnd);
   }
}

/*! FUNCTION: EDGEBOUNDS_Find_BoundingBox()
 *  SYNOPSIS: Find the min/max range of values contained in the edgebounds.
 *            Assumes edgebounds have been sorted.
 */
int 
EDGEBOUNDS_Find_BoundingBox(  EDGEBOUNDS*   edg,
                              RANGE*        Q_range_out,
                              RANGE*        T_range_out )
{
   RANGE*   Q_range;
   RANGE*   T_range; 
   int      N;
   int      Q, T;

   N        = EDGEBOUNDS_GetSize( edg );
   Q_range  = &edg->Q_range;
   T_range  = &edg->T_range; 
   Q        = edg->Q;
   T        = edg->T;

   /* find the target and query range */
   Q_range->beg = Q + 1;
   Q_range->end = 0;
   T_range->beg = T + 1;
   T_range->end = 0;

   /* if list is empty, return */
   if ( N == 0 ) {
      return STATUS_SUCCESS;
   }

   /* create bounding box */
   for (int i = 0; i < N; i++) 
   {
      if ( T_range->beg > EDGEBOUNDS_Get( edg, i ).lb ) {
         T_range->beg = EDGEBOUNDS_Get( edg, i ).lb;
      }
      if ( T_range->end < EDGEBOUNDS_Get( edg, i ).rb ) {
         T_range->end = EDGEBOUNDS_Get( edg, i ).rb;
      }
   }
   Q_range->beg = EDGEBOUNDS_Get( edg, 0 ).id;
   Q_range->end = EDGEBOUNDS_Get( edg, N-1 ).id;

   /* edge checks */
   T_range->beg = MAX(T_range->beg, 0);
   T_range->end = MIN(T_range->end, T + 1);
   Q_range->beg = MAX(Q_range->beg, 0);
   Q_range->end = MIN(Q_range->end, Q + 1);

   /* output if not null */
   if ( Q_range_out != NULL ) {
      *Q_range_out = *Q_range;
   }
   if ( T_range_out != NULL ) {
      *T_range_out = *T_range;
   }

   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_SetDomain()
 *  SYNOPSIS: Build an EDGEBOUND <edg_out> from QxT EDGEBOUNDS <edg_in>
 *            and constraining the query range to <in_dom_range> => <q_beg, q_end>.
 *            Simply shifts all query indexes by <q_beg>.
 *            NOTE: Assumes <edg_in> has been sorted.
 */
int 
EDGEBOUNDS_SetDomain(   EDGEBOUNDS*    edg_in,
                        EDGEBOUNDS*    edg_out,
                        RANGE          Q_range )
{
   BOUND    bnd;
   int      Q_len;
   int      N;

   Q_len = Q_range.end - Q_range.beg + 1; 
   N     = EDGEBOUNDS_GetSize( edg_in );

   /* clears old data */
   EDGEBOUNDS_Reuse( edg_out, Q_len, edg_in->T );
   /* shift all query indexes to set <q_beg> to 0. */
   for (int i = 0; i < N; i++)
   {
      bnd = EDG_X(edg_in, i);
      bnd.id -= Q_range.beg;
      EDGEBOUNDS_Pushback( edg_out, bnd );
   }

   return STATUS_SUCCESS;
}

