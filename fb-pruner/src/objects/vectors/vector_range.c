/*******************************************************************************
 *  FILE:      vector_range.c
 *  PURPOSE:   VECTOR_RANGE Object Functions
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

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "vector_template.h"

/*
 *  FUNCTION:  VECTOR_RANGE_Create()
 *  SYNOPSIS:  Create new VECTOR_RANGE object and returns pointer.
 */
VECTOR_RANGE* VECTOR_RANGE_Create()
{
   const int init_size = VECTOR_INIT_SIZE;
   return VECTOR_RANGE_Create_by_Size( init_size );
}

/*
 *  FUNCTION:  VECTOR_RANGE_Create()
 *  SYNOPSIS:  Create new VECTOR_RANGE object at specific size and returns pointer.
 */
VECTOR_RANGE* VECTOR_RANGE_Create_by_Size( int    size )
{
   VECTOR_RANGE *vec = NULL;
   vec = (VECTOR_RANGE *) malloc( sizeof(VECTOR_RANGE) );
   if ( vec == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n");
      exit(EXIT_FAILURE);
   }

   vec->data   = NULL;
   vec->N      = 0;
   vec->Nalloc = 0;

   VECTOR_RANGE_Resize( vec, size );

   return vec;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_RANGE.
 */
void* VECTOR_RANGE_Destroy( VECTOR_RANGE*   vec )
{
   if ( vec == NULL ) return vec;
   free(vec->data);
   free(vec);

   vec = NULL;
   return vec;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_RANGE object by resetting size counter (no realloc) .
 */
void VECTOR_RANGE_Reuse( VECTOR_RANGE*   vec )
{
   vec->N = 0;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Fill()
 *  SYNOPSIS:  Fill VECTOR_RANGE object with val.
 */
void VECTOR_RANGE_Fill(   VECTOR_RANGE*   vec, 
                        RANGE           val )
{
   for ( int i = 0; i < vec->N; i++ ) {
      vec->data[i] = val;
   }
}

/*
 *  FUNCTION:  VECTOR_RANGE_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_RANGE for <dest> if <dest> is NULL.
 */
VECTOR_RANGE* VECTOR_RANGE_Copy(  VECTOR_RANGE*   src, 
                              VECTOR_RANGE*   dest )
{

   if ( dest == NULL ) {
      dest = VECTOR_RANGE_Create();
   }
   /* allocate variable-sized data */
   VECTOR_RANGE_Resize( dest, src->Nalloc );
   /* copy variable-sized data */
   memcpy( dest->data, src->data, sizeof(RANGE) * src->N );
   /* copy base data */
   dest->N = src->N;

   return dest;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_RANGE_Resize(    VECTOR_RANGE*   vec, 
                           int           size )
{
   vec->data = (RANGE*) realloc( vec->data, sizeof(RANGE) * size );
   if ( vec->data == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n" );
      exit(EXIT_FAILURE);
   }
   vec->Nalloc = size;
}

/*
 *  FUNCTION:  VECTOR_RANGE_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_RANGE_GrowTo( VECTOR_RANGE*   vec, 
                        int           size )
{
   if ( vec->Nalloc < size ) {
      VECTOR_RANGE_Resize( vec, size );
   }
}

/*
 *  FUNCTION:  VECTOR_RANGE_Push()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array. 
 *             Warning: Does not handle resizing or check for out-of-bounds. For that, use Pushback().
 */
inline
void VECTOR_RANGE_Push(   VECTOR_RANGE*   vec, 
                        RANGE           val )
{
   vec->data[vec->N] = val;
   vec->N++;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
inline
void VECTOR_RANGE_Pushback(  VECTOR_RANGE*   vec, 
                           RANGE           val )
{
   VECTOR_RANGE_Push( vec, val );

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_RANGE_Resize( vec, vec->N * 2 );
   }
}

/*
 *  FUNCTION:  VECTOR_RANGE_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
inline
RANGE VECTOR_RANGE_Pop( VECTOR_RANGE*   vec )
{
   RANGE data = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_RANGE_Resize( vec, vec->N / 2 );
   }

   return data;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_RANGE_Set(    VECTOR_RANGE*   vec, 
                        int           idx, 
                        RANGE           val )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_RANGE access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   vec->data[idx] = val;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return value of data.
 *  RETURN:    Return data at <idx>.
 */
inline
RANGE VECTOR_RANGE_Get(  VECTOR_RANGE*   vec, 
                     int           idx )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_RANGE access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return &(vec->data[idx]);
}

/*
 *  FUNCTION:  VECTOR_RANGE_Get_X()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
inline
RANGE* VECTOR_RANGE_Get_X(  VECTOR_RANGE*   vec, 
                        int           idx )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_RANGE access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return &(vec->data[idx]);
}

/*
 *  FUNCTION:  VECTOR_RANGE_Get_Size()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
inline
RANGE VECTOR_RANGE_Get_Size(   VECTOR_RANGE*   vec )
{
   return vec->N;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Set_Size()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
inline
void VECTOR_RANGE_Set_Size(  VECTOR_RANGE*   vec, 
                           int           size )
{
   vec->N = size;
}


/*
 *  FUNCTION:  VECTOR_RANGE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_RANGE_Search(  VECTOR_RANGE*   vec, 
                        RANGE           val )
{
   int N       = VECTOR_RANGE_Get_Size( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = RANGE_Compare( val, vec->data[i] );

      if ( cmp > 0) {
         idx += i;
      }
      else if ( cmp < 0 ) {
         idx -= i;
      }
      else {
         found = idx;
         return found;
      }
   }
   return found;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_RANGE_Search_First(  VECTOR_RANGE*   vec, 
                              RANGE           val )
{
   int N       = VECTOR_RANGE_Get_Size( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = RANGE_Compare( val, vec->data[i] );

      if ( cmp > 0) {
         idx += i;
      }
      else if ( cmp < 0 ) {
         idx -= i;
      }
      else {
         found = idx;
         idx -= i;
      }
   }
   return found;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_RANGE_Search_Last(   VECTOR_RANGE*   vec, 
                              RANGE           val )
{
   int N       = VECTOR_RANGE_Get_Size( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = RANGE_Compare( val, vec->data[i] );

      if ( cmp > 0) {
         idx += i;
      }
      else if ( cmp < 0 ) {
         idx -= i;
      }
      else {
         found = idx;
         idx += i;
      }
   }
   return found;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_RANGE_Compare(    VECTOR_RANGE*   vec_A, 
                           VECTOR_RANGE*   vec_B )
{
   for (int i = 0; i < vec_A->N; i++) 
   {
      if ( RANGE_Compare( vec_A->data[i], vec_B->data[i] ) != 0 ) 
      {
         return RANGE_Compare( vec_A->data[i], vec_B->data[i] ) != 0 );       
      }
   }
   return 0;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
void VECTOR_RANGE_Sort( VECTOR_RANGE*    vec )
{
   int N = VECTOR_RANGE_Get_Size(vec);
   VECTOR_RANGE_Sort_Sub( vec, 0, N );
}

/*
 *  FUNCTION:  VECTOR_RANGE_Sort_Sub()
 *  SYNOPSIS:  Sorts subarray of <vec> data in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
void VECTOR_RANGE_Sort_Sub(  VECTOR_RANGE*    vec,
                           int            beg,
                           int            end )
{
   const int begin_select_sort = 16;
   
   /* run selection sort if below threshold */
   if ( end - beg > begin_select_sort ) {
      VECTOR_RANGE_Sort_Sub_Selectsort( vec, beg, end );
   }
   VECTOR_RANGE_Sort_Sub_Quicksort( vec, beg, end );
}

/*
 *  FUNCTION:  VECTOR_RANGE_Sort_Sub_Selectsort()
 *  SYNOPSIS:  Selection Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
void VECTOR_RANGE_Sort_Sub_Selectsort(   VECTOR_RANGE*    vec,
                                       int            beg,
                                       int            end )
{
   for (int i = beg; i < end; i++) 
   {
      /* initial minimum value found */
      int min_idx = i;
      RANGE min_val = vec->data[i];
      for (int j = i+1; j < end; j++) {
         /* if new minimum found, update value and index */
         int cmp = RANGE_Compare( min_val, vec->data[j] );
         if ( cmp < 0 ) {
            min_idx = j;
            min_val = vec->data[j];
         }
      }
      /* swap new minimum to left-most position */
      VECTOR_RANGE_Swap( vec, i, min_idx );
   }
   return;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Quick Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
void VECTOR_RANGE_Sort_Sub_Quicksort( VECTOR_RANGE*    vec,
                                    int            beg,
                                    int            end )
{
   /* partition pointers */
   int   r_idx = beg + 1;
   int   l_idx = end - 1;
   RANGE*  rhs   = &(vec->data[beg + 1]);
   RANGE*  lhs   = &(vec->data[end - 1]);

   /* select random pivot value */
   int   range = end - beg;
   int pivot_idx = (rand() % range) + beg;
   RANGE pivot_val = vec->data[pivot_idx];
   VECTOR_RANGE_Swap( vec, pivot, beg );

   /* partition on pivot */
   while ( l_idx <= r_idx )
   {
      while ( (l_idx <= r_idx) && (RANGE_Compare( pivot_val, vec->data[r_idx] ) < 0) ) {
         r_idx--;
      }
      while ( (l_idx <= r_idx) && (RANGE_Compare( pivot_val, vec->data[l_idx] ) >= 0) ) {
         l_idx++;
      }
      if ( l_idx <= r_idx ) {
         VECTOR_RANGE_Swap( vec, l_idx, r_idx );
         r_idx--;
         l_idx++;
      }
   }
   VECTOR_RANGE_Swap( vec, beg, r_idx );
}

/*
 *  FUNCTION:  VECTOR_RANGE_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
inline
void VECTOR_RANGE_Swap(   VECTOR_RANGE*    vec,
                        int            i,
                        int            j )
{
   RANGE swap = vec->data[i];
   vec->data[i] = vec->data[j];
   vec->data[j] = swap;
}

/*
 *  FUNCTION:  VECTOR_RANGE_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
inline
void VECTOR_RANGE_Reverse(   VECTOR_RANGE*    vec )
{
   int N = VECTOR_RANGE_Get_Size( vec );

   for (int i = 0; i < (N/2); i++) {
      VECTOR_RANGE_Swap( vec, i, (N-1)-i );
   }
}

/*
 *  FUNCTION:  VECTOR_RANGE_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_RANGE_Dump(   VECTOR_RANGE*    vec,
                        FILE*          fp )
{
   fprintf(fp, "%s: ", "VECTOR_RANGE");
   fprintf(fp, "[ ");
   for ( int i = 0; i < vec->N; i++ ) {
      fprintf(fp, "%3d ", vec->data[i] );
   }
   fprintf(fp, "]\n" );
}
