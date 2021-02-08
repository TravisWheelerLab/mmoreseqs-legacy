/*******************************************************************************
 *  FILE:      vector_int.c
 *  PURPOSE:   VECTOR_INT Object Functions.
 *             Template for building vector classes.
 *             Run "scripts/builder-helper/build_vector_classes_from_template" to update.
 *             Requires data primitive to have INT_Compare().
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
#include "vector_int.h"

/*! FUNCTION:  VECTOR_INT_Create()
 *  SYNOPSIS:  Create new VECTOR_INT object and returns pointer.
 */
VECTOR_INT* 
VECTOR_INT_Create()
{
   const int init_size = VECTOR_INIT_SIZE;
   return VECTOR_INT_Create_by_Size( init_size );
}

/*! FUNCTION:  VECTOR_INT_Create()
 *  SYNOPSIS:  Create new VECTOR_INT object at specific size and returns pointer.
 */
VECTOR_INT* 
VECTOR_INT_Create_by_Size( size_t    size )
{
   VECTOR_INT* vec;
   vec = ERROR_malloc( sizeof(VECTOR_INT) );

   vec->data   = NULL;
   vec->N      = 0;
   vec->Nalloc = 0;

   VECTOR_INT_Resize( vec, size );

   return vec;
}

/*! FUNCTION:  VECTOR_INT_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_INT.
 */
VECTOR_INT* 
VECTOR_INT_Destroy( VECTOR_INT*   vec )
{
   if ( vec == NULL ) return NULL;

   /* NOTE: This loop prevents memory leaks with data types which allocate dynamic memory 
    *       Hopefully, this can be safely used in template, as this loop should be 
    *       inlined and optimized out by non-dynamic data types(?)
    */
   // for (int i = 0; i < vec->N; i++) {
   //    VEC_X( vec, i ) = INT_Destroy( VEC_X( vec, i ) );
   // }

   vec->data   = ERROR_free(vec->data);
   vec         = ERROR_free(vec);

   return NULL;
}

/*! FUNCTION:  VECTOR_INT_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_INT object by resetting size counter (no realloc) .
 */
STATUS_FLAG 
VECTOR_INT_Reuse( VECTOR_INT*   vec )
{
   /* NOTE: This loop prevents memory leaks with data types which allocate dynamic memory 
    *       Hopefully, this can be safely used in template, as this loop should be 
    *       inlined and optimized out by non-dynamic data types(?)
    */
   // for (int i = 0; i < vec->N; i++) {
   //    VEC_X( vec, i ) = INT_Destroy( VEC_X( vec, i ) );
   // }

   vec->N = 0;
}

/*! FUNCTION:  VECTOR_INT_Fill()
 *  SYNOPSIS:  Fill VECTOR_INT object with val.
 */
STATUS_FLAG 
VECTOR_INT_Fill(  VECTOR_INT*   vec, 
                  INT           val )
{
   for ( int i = 0; i < vec->N; i++ ) {
      VEC_X( vec, i ) = INT_Create( val );
   }
}

/*! FUNCTION:  VECTOR_INT_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_INT for <dest> if <dest> is NULL.
 */
VECTOR_INT* 
VECTOR_INT_Copy(  VECTOR_INT*   dest, 
                  VECTOR_INT*   src )
{

   if ( dest == NULL ) {
      dest = VECTOR_INT_Create();
   }
   /* allocate variable-sized data */
   VECTOR_INT_Resize( dest, src->N );
   /* copy variable-sized data */
   for (int i = 0; i < src->N; i++ ) {
      VEC_X( dest, i ) = INT_Create( VEC_X( src, i ) );
   }
   /* copy base data */
   dest->N = src->N;

   return dest;
}

/*! FUNCTION:  VECTOR_INT_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
STATUS_FLAG 
VECTOR_INT_Resize(   VECTOR_INT*   vec, 
                     size_t        size )
{
   vec->data = ERROR_realloc( vec->data, sizeof(INT) * size );
   vec->Nalloc = size;
}

/*! FUNCTION:  VECTOR_INT_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
STATUS_FLAG 
VECTOR_INT_GrowTo(   VECTOR_INT*   vec, 
                     size_t        size )
{
   if ( vec->Nalloc < size ) {
      VECTOR_INT_Resize( vec, size );
   }
}

/*! FUNCTION:  VECTOR_INT_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
inline
INT* 
VECTOR_INT_GetArray(   VECTOR_INT*   vec )
{
   return vec->data;
}

/*! FUNCTION:  VECTOR_INT_Push()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array. 
 *             Warning: Does not handle resizing or check for out-of-bounds. For that, use Pushback().
 */
inline
STATUS_FLAG 
VECTOR_INT_Push(  VECTOR_INT*   vec, 
                  INT           val )
{
   /* NOTE: This push() creates another copy of the data to store in vector (in the case of dynamically allocated data) */
   VEC_X( vec, vec->N ) = INT_Create( val );
   vec->N++;
}

/*! FUNCTION:  VECTOR_INT_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
inline
STATUS_FLAG 
VECTOR_INT_Pushback(    VECTOR_INT*   vec, 
                        INT           val )
{
   VECTOR_INT_Push( vec, val );

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_INT_Resize( vec, vec->N * 2 );
   }
}

/*! FUNCTION:  VECTOR_INT_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
inline
INT 
VECTOR_INT_Pop( VECTOR_INT*   vec )
{
   INT data = VEC_X( vec, vec->N - 1 );
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_INT_Resize( vec, vec->N / 2 );
   }

   return data;
}

/*! FUNCTION:  VECTOR_INT_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data array. 
 */
inline
STATUS_FLAG 
VECTOR_INT_Append(   VECTOR_INT*   vec, 
                     INT*          append,
                     size_t        L )
{
   size_t N_new;
   N_new = vec->N + L;

   /* resize array */
   if ( vec->Nalloc < N_new ) {
      VECTOR_INT_Resize( vec, N_new );
   }

   /* copy data over */
   for (int i = 0; i < L; i++) {
      VECTOR_INT_Push( vec, append[i] );
   }
}

/*! FUNCTION:  VECTOR_INT_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG 
VECTOR_INT_Set(   VECTOR_INT*   vec, 
                  int           idx, 
                  INT           val )
{
   /* if debugging, do edgebound checks */
   #if SAFE 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_INT access out-of-bounds.\n");
      fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   VEC_X( vec, idx ) = INT_Destroy( VEC_X( vec, idx ) );
   VEC_X( vec, idx ) = INT_Create( val );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_INT_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return value of data.
 *  RETURN:    Return data at <idx>.
 */
inline
INT 
VECTOR_INT_Get(   VECTOR_INT*   vec, 
                  int           idx )
{
   /* if debugging, do edgebound checks */
   #if SAFE 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_INT access out-of-bounds.\n");
      fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return (vec->data[idx]);
}

/*! FUNCTION:  VECTOR_INT_GetX()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
inline
INT* 
VECTOR_INT_GetX( VECTOR_INT*   vec, 
                  int           idx )
{
   /* if debugging, do edgebound checks */
   #if SAFE 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_INT access out-of-bounds.\n");
      fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return &(vec->data[idx]);
}

/*! FUNCTION:  VECTOR_INT_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
inline
int 
VECTOR_INT_GetSize(   VECTOR_INT*   vec )
{
   return vec->N;
}

/*! FUNCTION:  VECTOR_INT_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>. 
 *             Will allocate memory if necessary.
 */
inline
STATUS_FLAG 
VECTOR_INT_SetSize(     VECTOR_INT*   vec, 
                        size_t        size )
{
   VECTOR_INT_GrowTo( vec, size );
   vec->N = size;
}


/*! FUNCTION:  VECTOR_INT_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_INT_Search(   VECTOR_INT*   vec, 
                     INT           val )
{
   int N       = VECTOR_INT_GetSize( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = INT_Compare( val, vec->data[i] );

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

/*! FUNCTION:  VECTOR_INT_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_INT_Search_First(  VECTOR_INT*   vec, 
                           INT           val )
{
   int N       = VECTOR_INT_GetSize( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = INT_Compare( val, vec->data[i] );

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

/*! FUNCTION:  VECTOR_INT_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_INT_Search_Last(    VECTOR_INT*   vec, 
                           INT           val )
{
   int N       = VECTOR_INT_GetSize( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = INT_Compare( val, vec->data[i] );

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

/*! FUNCTION:  VECTOR_INT_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int 
VECTOR_INT_Compare(  VECTOR_INT*   vec_A, 
                     VECTOR_INT*   vec_B )
{
   for (int i = 0; i < vec_A->N; i++) {
      if ( INT_Compare( vec_A->data[i], vec_B->data[i] ) != 0 ) {
         return INT_Compare( vec_A->data[i], vec_B->data[i] );       
      }
   }
   return 0;
}

/*! FUNCTION:  VECTOR_INT_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG 
VECTOR_INT_Sort( VECTOR_INT*    vec )
{
   int N = VECTOR_INT_GetSize( vec );
   VECTOR_INT_Sort_Sub( vec, 0, N );

   #if DEBUG 
   {
      for ( int i = 0; i < N-1; i++ ) 
      {
         INT cur = vec->data[i];
         INT nxt = vec->data[i+1];
         char s_cur[50];
         char s_nxt[50];
         int cmp = INT_Compare( cur, nxt );
         if ( (cmp <= 0) == false ) {
            fprintf(stderr, "ERROR: bad sort. %d, %d v %d: %s vs %s\n",
               cmp, i, i+1, INT_To_String(cur, s_cur), INT_To_String(nxt, s_nxt) );
         }
      }
   }
   #endif
}

/*! FUNCTION:  VECTOR_INT_Sort_Sub()
 *  SYNOPSIS:  Sorts subarray of <vec> data in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
STATUS_FLAG 
VECTOR_INT_Sort_Sub(    VECTOR_INT*    vec,
                        int            beg,
                        int            end )
{
   const int begin_select_sort = 4;
   int N = VECTOR_INT_GetSize(vec);
   if (N <= 1) {
      return STATUS_SUCCESS;
   } 

   int size = end - beg;
   /* run selection sort if below threshold */
   if ( size <= begin_select_sort ) {
      VECTOR_INT_Sort_Sub_Selectsort( vec, beg, end );
   } 
   /* otherwise run quiksort */
   else {
      VECTOR_INT_Sort_Sub_Quicksort( vec, beg, end );
   }

   return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_INT_Sort_Sub_Selectsort()
 *  SYNOPSIS:  Selection Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
STATUS_FLAG 
VECTOR_INT_Sort_Sub_Selectsort(     VECTOR_INT*    vec,
                                    int            beg,
                                    int            end )
{
   /* find the minimum element of remaining unsorted list */
   for (int i = beg; i < end; i++) 
   {
      /* initial minimum value found */
      int min_idx = i;
      INT min_val = vec->data[i];
      for (int j = i+1; j < end; j++) {
         /* if new minimum found, update value and index */
         int cmp = INT_Compare( min_val, vec->data[j] );
         if ( cmp > 0 ) {
            min_idx = j;
            min_val = vec->data[j];
         }
      }
      /* swap new minimum to left-most position */
      VECTOR_INT_Swap( vec, i, min_idx );
   }
}

/*! FUNCTION:  VECTOR_INT_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Quick Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
STATUS_FLAG 
VECTOR_INT_Sort_Sub_Quicksort(   VECTOR_INT*    vec,
                                 int            beg,
                                 int            end )
{
   /* partition pointers */
   int   r_idx = beg + 1;
   int   l_idx = end - 1;
   INT*  rhs   = &(vec->data[beg + 1]);
   INT*  lhs   = &(vec->data[end - 1]);

   /* select random pivot value */
   int   range       = end - beg;
   int   pivot_idx   = RNG_INT_Range( beg, end );
   INT   pivot_val   = vec->data[pivot_idx];
   VECTOR_INT_Swap( vec, pivot_idx, beg );

   /* partition on pivot */
   while ( l_idx <= r_idx )
   {
      /* find next right partition element that is less than pivot element */
      while ( (l_idx <= r_idx) && (INT_Compare( pivot_val, vec->data[r_idx] ) < 0) ) {
         r_idx--;
      }
      /* find next left partition element that is greater than pivot element */
      while ( (l_idx <= r_idx) && (INT_Compare( pivot_val, vec->data[l_idx] ) >= 0) ) {
         l_idx++;
      }
      /* if left and right index have not crossed, then swap elements */
      if ( l_idx <= r_idx ) {
         VECTOR_INT_Swap( vec, l_idx, r_idx );
      }
   }
   /* move partition element to barrier between left and right index */
   VECTOR_INT_Swap( vec, beg, r_idx );
   /* sort both partitions (omit partition element) */
   VECTOR_INT_Sort_Sub( vec, beg, r_idx );
   VECTOR_INT_Sort_Sub( vec, r_idx+1, end );
}

/*! FUNCTION:  VECTOR_INT_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
inline
STATUS_FLAG 
VECTOR_INT_Swap(  VECTOR_INT*    vec,
                  int            i,
                  int            j )
{
   INT swap = vec->data[i];
   vec->data[i] = vec->data[j];
   vec->data[j] = swap;
}

/*! FUNCTION:  VECTOR_INT_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
inline
STATUS_FLAG VECTOR_INT_Reverse(   VECTOR_INT*    vec )
{
   int N = VECTOR_INT_GetSize( vec );

   for (int i = 0; i < (N/2); i++) {
      VECTOR_INT_Swap( vec, i, (N-1)-i );
   }
}

/*! FUNCTION:  VECTOR_INT_LoadTSV()
 *  SYNOPSIS:  Load tsv from <filename> and store in <vec>.
 */
STATUS_FLAG 
VECTOR_INT_LoadTSV(  VECTOR_INT*    vec,
                     char*          filename )
{
   FILE* fp  = fopen(filename, "r+");
   char buf[1024]; 

   fclose(fp);
}

/*! FUNCTION:  VECTOR_INT_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
STATUS_FLAG 
VECTOR_INT_Dump(  VECTOR_INT*    vec,
                  FILE*          fp )
{
   /* stringification of template object */
   char s[50];

   fprintf(fp, "%s: ", "VECTOR_INT");
   fprintf(fp, "[ ");
   for ( int i = 0; i < vec->N; i++ ) {
      fprintf(fp, "%s, ", INT_To_String(vec->data[i], s) );
   }
   fprintf(fp, "]\n" );
}


/*! FUNCTION:  VECTOR_INT_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_INT.
 */
STATUS_FLAG 
VECTOR_INT_Unit_Test()
{
   VECTOR_INT* vec = VECTOR_INT_Create();

   VECTOR_INT_Destroy( vec );
}