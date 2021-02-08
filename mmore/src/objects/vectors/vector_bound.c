/*******************************************************************************
 *  FILE:      vector_bound.c
 *  PURPOSE:   VECTOR_BOUND Object Functions.
 *             Template for building vector classes.
 *             Run "scripts/builder-helper/build_vector_classes_from_template" to update.
 *             Requires data primitive to have BOUND_Compare().
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

/* local imports */
#include "../../objects/structs.h"
#include "../../utilities/_utilities.h"
#include "../../objects/_objects.h"

/* header */
#include "_vectors.h"
#include "vector_bound.h"

/*! FUNCTION:  VECTOR_BOUND_Create()
 *  SYNOPSIS:  Create new VECTOR_BOUND object and returns pointer.
 */
VECTOR_BOUND* 
VECTOR_BOUND_Create()
{
   const int init_size = VECTOR_INIT_SIZE;
   return VECTOR_BOUND_Create_by_Size( init_size );
}

/*! FUNCTION:  VECTOR_BOUND_Create()
 *  SYNOPSIS:  Create new VECTOR_BOUND object at specific size and returns pointer.
 */
VECTOR_BOUND* 
VECTOR_BOUND_Create_by_Size( size_t    size )
{
   VECTOR_BOUND* vec;
   vec = ERROR_malloc( sizeof(VECTOR_BOUND) );

   vec->data   = NULL;
   vec->N      = 0;
   vec->Nalloc = 0;

   VECTOR_BOUND_Resize( vec, size );

   return vec;
}

/*! FUNCTION:  VECTOR_BOUND_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_BOUND.
 */
VECTOR_BOUND* 
VECTOR_BOUND_Destroy( VECTOR_BOUND*   vec )
{
   if ( vec == NULL ) return NULL;

   /* NOTE: This loop prevents memory leaks with data types which allocate dynamic memory 
    *       Hopefully, this can be safely used in template, as this loop should be 
    *       inlined and optimized out by non-dynamic data types(?)
    */
   // for (int i = 0; i < vec->N; i++) {
   //    VEC_X( vec, i ) = BOUND_Destroy( VEC_X( vec, i ) );
   // }

   vec->data   = ERROR_free(vec->data);
   vec         = ERROR_free(vec);

   return NULL;
}

/*! FUNCTION:  VECTOR_BOUND_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_BOUND object by resetting size counter (no realloc) .
 */
STATUS_FLAG 
VECTOR_BOUND_Reuse( VECTOR_BOUND*   vec )
{
   /* NOTE: This loop prevents memory leaks with data types which allocate dynamic memory 
    *       Hopefully, this can be safely used in template, as this loop should be 
    *       inlined and optimized out by non-dynamic data types(?)
    */
   // for (int i = 0; i < vec->N; i++) {
   //    VEC_X( vec, i ) = BOUND_Destroy( VEC_X( vec, i ) );
   // }

   vec->N = 0;
}

/*! FUNCTION:  VECTOR_BOUND_Fill()
 *  SYNOPSIS:  Fill VECTOR_BOUND object with val.
 */
STATUS_FLAG 
VECTOR_BOUND_Fill(  VECTOR_BOUND*   vec, 
                  BOUND           val )
{
   for ( int i = 0; i < vec->N; i++ ) {
      VEC_X( vec, i ) = BOUND_Create( val );
   }
}

/*! FUNCTION:  VECTOR_BOUND_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_BOUND for <dest> if <dest> is NULL.
 */
VECTOR_BOUND* 
VECTOR_BOUND_Copy(  VECTOR_BOUND*   dest, 
                  VECTOR_BOUND*   src )
{

   if ( dest == NULL ) {
      dest = VECTOR_BOUND_Create();
   }
   /* allocate variable-sized data */
   VECTOR_BOUND_Resize( dest, src->N );
   /* copy variable-sized data */
   for (int i = 0; i < src->N; i++ ) {
      VEC_X( dest, i ) = BOUND_Create( VEC_X( src, i ) );
   }
   /* copy base data */
   dest->N = src->N;

   return dest;
}

/*! FUNCTION:  VECTOR_BOUND_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
STATUS_FLAG 
VECTOR_BOUND_Resize(   VECTOR_BOUND*   vec, 
                     size_t        size )
{
   vec->data = ERROR_realloc( vec->data, sizeof(BOUND) * size );
   vec->Nalloc = size;
}

/*! FUNCTION:  VECTOR_BOUND_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
STATUS_FLAG 
VECTOR_BOUND_GrowTo(   VECTOR_BOUND*   vec, 
                     size_t        size )
{
   if ( vec->Nalloc < size ) {
      VECTOR_BOUND_Resize( vec, size );
   }
}

/*! FUNCTION:  VECTOR_BOUND_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
inline
BOUND* 
VECTOR_BOUND_GetArray(   VECTOR_BOUND*   vec )
{
   return vec->data;
}

/*! FUNCTION:  VECTOR_BOUND_Push()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array. 
 *             Warning: Does not handle resizing or check for out-of-bounds. For that, use Pushback().
 */
inline
STATUS_FLAG 
VECTOR_BOUND_Push(  VECTOR_BOUND*   vec, 
                  BOUND           val )
{
   /* NOTE: This push() creates another copy of the data to store in vector (in the case of dynamically allocated data) */
   VEC_X( vec, vec->N ) = BOUND_Create( val );
   vec->N++;
}

/*! FUNCTION:  VECTOR_BOUND_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
inline
STATUS_FLAG 
VECTOR_BOUND_Pushback(    VECTOR_BOUND*   vec, 
                        BOUND           val )
{
   VECTOR_BOUND_Push( vec, val );

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_BOUND_Resize( vec, vec->N * 2 );
   }
}

/*! FUNCTION:  VECTOR_BOUND_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
inline
BOUND 
VECTOR_BOUND_Pop( VECTOR_BOUND*   vec )
{
   BOUND data = VEC_X( vec, vec->N - 1 );
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_BOUND_Resize( vec, vec->N / 2 );
   }

   return data;
}

/*! FUNCTION:  VECTOR_BOUND_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data array. 
 */
inline
STATUS_FLAG 
VECTOR_BOUND_Append(   VECTOR_BOUND*   vec, 
                     BOUND*          append,
                     size_t        L )
{
   size_t N_new;
   N_new = vec->N + L;

   /* resize array */
   if ( vec->Nalloc < N_new ) {
      VECTOR_BOUND_Resize( vec, N_new );
   }

   /* copy data over */
   for (int i = 0; i < L; i++) {
      VECTOR_BOUND_Push( vec, append[i] );
   }
}

/*! FUNCTION:  VECTOR_BOUND_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG 
VECTOR_BOUND_Set(   VECTOR_BOUND*   vec, 
                  int           idx, 
                  BOUND           val )
{
   /* if debugging, do edgebound checks */
   #if SAFE 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_BOUND access out-of-bounds.\n");
      fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   VEC_X( vec, idx ) = BOUND_Destroy( VEC_X( vec, idx ) );
   VEC_X( vec, idx ) = BOUND_Create( val );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_BOUND_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return value of data.
 *  RETURN:    Return data at <idx>.
 */
inline
BOUND 
VECTOR_BOUND_Get(   VECTOR_BOUND*   vec, 
                  int           idx )
{
   /* if debugging, do edgebound checks */
   #if SAFE 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_BOUND access out-of-bounds.\n");
      fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return (vec->data[idx]);
}

/*! FUNCTION:  VECTOR_BOUND_GetX()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
inline
BOUND* 
VECTOR_BOUND_GetX( VECTOR_BOUND*   vec, 
                  int           idx )
{
   /* if debugging, do edgebound checks */
   #if SAFE 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_BOUND access out-of-bounds.\n");
      fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return &(vec->data[idx]);
}

/*! FUNCTION:  VECTOR_BOUND_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
inline
int 
VECTOR_BOUND_GetSize(   VECTOR_BOUND*   vec )
{
   return vec->N;
}

/*! FUNCTION:  VECTOR_BOUND_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>. 
 *             Will allocate memory if necessary.
 */
inline
STATUS_FLAG 
VECTOR_BOUND_SetSize(     VECTOR_BOUND*   vec, 
                        size_t        size )
{
   VECTOR_BOUND_GrowTo( vec, size );
   vec->N = size;
}


/*! FUNCTION:  VECTOR_BOUND_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_BOUND_Search(   VECTOR_BOUND*   vec, 
                     BOUND           val )
{
   int N       = VECTOR_BOUND_GetSize( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = BOUND_Compare( val, vec->data[i] );

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

/*! FUNCTION:  VECTOR_BOUND_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_BOUND_Search_First(  VECTOR_BOUND*   vec, 
                           BOUND           val )
{
   int N       = VECTOR_BOUND_GetSize( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = BOUND_Compare( val, vec->data[i] );

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

/*! FUNCTION:  VECTOR_BOUND_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_BOUND_Search_Last(    VECTOR_BOUND*   vec, 
                           BOUND           val )
{
   int N       = VECTOR_BOUND_GetSize( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = BOUND_Compare( val, vec->data[i] );

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

/*! FUNCTION:  VECTOR_BOUND_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int 
VECTOR_BOUND_Compare(  VECTOR_BOUND*   vec_A, 
                     VECTOR_BOUND*   vec_B )
{
   for (int i = 0; i < vec_A->N; i++) {
      if ( BOUND_Compare( vec_A->data[i], vec_B->data[i] ) != 0 ) {
         return BOUND_Compare( vec_A->data[i], vec_B->data[i] );       
      }
   }
   return 0;
}

/*! FUNCTION:  VECTOR_BOUND_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG 
VECTOR_BOUND_Sort( VECTOR_BOUND*    vec )
{
   int N = VECTOR_BOUND_GetSize( vec );
   VECTOR_BOUND_Sort_Sub( vec, 0, N );

   #if DEBUG 
   {
      for ( int i = 0; i < N-1; i++ ) 
      {
         BOUND cur = vec->data[i];
         BOUND nxt = vec->data[i+1];
         char s_cur[50];
         char s_nxt[50];
         int cmp = BOUND_Compare( cur, nxt );
         if ( (cmp <= 0) == false ) {
            fprintf(stderr, "ERROR: bad sort. %d, %d v %d: %s vs %s\n",
               cmp, i, i+1, BOUND_To_String(cur, s_cur), BOUND_To_String(nxt, s_nxt) );
         }
      }
   }
   #endif
}

/*! FUNCTION:  VECTOR_BOUND_Sort_Sub()
 *  SYNOPSIS:  Sorts subarray of <vec> data in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
STATUS_FLAG 
VECTOR_BOUND_Sort_Sub(    VECTOR_BOUND*    vec,
                        int            beg,
                        int            end )
{
   const int begin_select_sort = 4;
   int N = VECTOR_BOUND_GetSize(vec);
   if (N <= 1) {
      return STATUS_SUCCESS;
   } 

   int size = end - beg;
   /* run selection sort if below threshold */
   if ( size <= begin_select_sort ) {
      VECTOR_BOUND_Sort_Sub_Selectsort( vec, beg, end );
   } 
   /* otherwise run quiksort */
   else {
      VECTOR_BOUND_Sort_Sub_Quicksort( vec, beg, end );
   }

   return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_BOUND_Sort_Sub_Selectsort()
 *  SYNOPSIS:  Selection Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
STATUS_FLAG 
VECTOR_BOUND_Sort_Sub_Selectsort(     VECTOR_BOUND*    vec,
                                    int            beg,
                                    int            end )
{
   /* find the minimum element of remaining unsorted list */
   for (int i = beg; i < end; i++) 
   {
      /* initial minimum value found */
      int min_idx = i;
      BOUND min_val = vec->data[i];
      for (int j = i+1; j < end; j++) {
         /* if new minimum found, update value and index */
         int cmp = BOUND_Compare( min_val, vec->data[j] );
         if ( cmp > 0 ) {
            min_idx = j;
            min_val = vec->data[j];
         }
      }
      /* swap new minimum to left-most position */
      VECTOR_BOUND_Swap( vec, i, min_idx );
   }
}

/*! FUNCTION:  VECTOR_BOUND_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Quick Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
STATUS_FLAG 
VECTOR_BOUND_Sort_Sub_Quicksort(   VECTOR_BOUND*    vec,
                                 int            beg,
                                 int            end )
{
   /* partition pointers */
   int   r_idx = beg + 1;
   int   l_idx = end - 1;
   BOUND*  rhs   = &(vec->data[beg + 1]);
   BOUND*  lhs   = &(vec->data[end - 1]);

   /* select random pivot value */
   int   range       = end - beg;
   int   pivot_idx   = RNG_INT_Range( beg, end );
   BOUND   pivot_val   = vec->data[pivot_idx];
   VECTOR_BOUND_Swap( vec, pivot_idx, beg );

   /* partition on pivot */
   while ( l_idx <= r_idx )
   {
      /* find next right partition element that is less than pivot element */
      while ( (l_idx <= r_idx) && (BOUND_Compare( pivot_val, vec->data[r_idx] ) < 0) ) {
         r_idx--;
      }
      /* find next left partition element that is greater than pivot element */
      while ( (l_idx <= r_idx) && (BOUND_Compare( pivot_val, vec->data[l_idx] ) >= 0) ) {
         l_idx++;
      }
      /* if left and right index have not crossed, then swap elements */
      if ( l_idx <= r_idx ) {
         VECTOR_BOUND_Swap( vec, l_idx, r_idx );
      }
   }
   /* move partition element to barrier between left and right index */
   VECTOR_BOUND_Swap( vec, beg, r_idx );
   /* sort both partitions (omit partition element) */
   VECTOR_BOUND_Sort_Sub( vec, beg, r_idx );
   VECTOR_BOUND_Sort_Sub( vec, r_idx+1, end );
}

/*! FUNCTION:  VECTOR_BOUND_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
inline
STATUS_FLAG 
VECTOR_BOUND_Swap(  VECTOR_BOUND*    vec,
                  int            i,
                  int            j )
{
   BOUND swap = vec->data[i];
   vec->data[i] = vec->data[j];
   vec->data[j] = swap;
}

/*! FUNCTION:  VECTOR_BOUND_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
inline
STATUS_FLAG VECTOR_BOUND_Reverse(   VECTOR_BOUND*    vec )
{
   int N = VECTOR_BOUND_GetSize( vec );

   for (int i = 0; i < (N/2); i++) {
      VECTOR_BOUND_Swap( vec, i, (N-1)-i );
   }
}

/*! FUNCTION:  VECTOR_BOUND_LoadTSV()
 *  SYNOPSIS:  Load tsv from <filename> and store in <vec>.
 */
STATUS_FLAG 
VECTOR_BOUND_LoadTSV(  VECTOR_BOUND*    vec,
                     char*          filename )
{
   FILE* fp  = fopen(filename, "r+");
   char buf[1024]; 

   fclose(fp);
}

/*! FUNCTION:  VECTOR_BOUND_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
STATUS_FLAG 
VECTOR_BOUND_Dump(  VECTOR_BOUND*    vec,
                  FILE*          fp )
{
   /* stringification of template object */
   char s[50];

   fprintf(fp, "%s: ", "VECTOR_BOUND");
   fprintf(fp, "[ ");
   for ( int i = 0; i < vec->N; i++ ) {
      fprintf(fp, "%s, ", BOUND_To_String(vec->data[i], s) );
   }
   fprintf(fp, "]\n" );
}


/*! FUNCTION:  VECTOR_BOUND_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_BOUND.
 */
STATUS_FLAG 
VECTOR_BOUND_Unit_Test()
{
   VECTOR_BOUND* vec = VECTOR_BOUND_Create();

   VECTOR_BOUND_Destroy( vec );
}