/*******************************************************************************
 *  FILE:      vector_char.c
 *  PURPOSE:   VECTOR_CHAR Object Functions.
 *             Template for building vector classes.
 *             Run "scripts/builder-helper/build_vector_classes_from_template" to update.
 *             Requires data primitive to have CHAR_Compare().
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
#include "vector_char.h"

/** FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Create new VECTOR_CHAR object and returns pointer.
 */
VECTOR_CHAR* VECTOR_CHAR_Create()
{
   const int init_size = VECTOR_INIT_SIZE;
   return VECTOR_CHAR_Create_by_Size( init_size );
}

/** FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Create new VECTOR_CHAR object at specific size and returns pointer.
 */
VECTOR_CHAR* VECTOR_CHAR_Create_by_Size( int    size )
{
   VECTOR_CHAR *vec = NULL;
   vec = (VECTOR_CHAR *) malloc( sizeof(VECTOR_CHAR) );
   if ( vec == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n");
      exit(EXIT_FAILURE);
   }

   vec->data   = NULL;
   vec->N      = 0;
   vec->Nalloc = 0;

   VECTOR_CHAR_Resize( vec, size );

   return vec;
}

/*
 *  FUNCTION:  VECTOR_CHAR_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_CHAR.
 */
void* VECTOR_CHAR_Destroy( VECTOR_CHAR*   vec )
{
   if ( vec == NULL ) return NULL;

   free(vec->data);
   free(vec);

   return NULL;
}

/*
 *  FUNCTION:  VECTOR_CHAR_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_CHAR object by resetting size counter (no realloc) .
 */
void VECTOR_CHAR_Reuse( VECTOR_CHAR*   vec )
{
   vec->N = 0;
}

/*
 *  FUNCTION:  VECTOR_CHAR_Fill()
 *  SYNOPSIS:  Fill VECTOR_CHAR object with val.
 */
void VECTOR_CHAR_Fill(   VECTOR_CHAR*   vec, 
                        CHAR           val )
{
   for ( int i = 0; i < vec->N; i++ ) {
      vec->data[i] = val;
   }
}

/** FUNCTION:  VECTOR_CHAR_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_CHAR for <dest> if <dest> is NULL.
 */
VECTOR_CHAR* VECTOR_CHAR_Copy(   VECTOR_CHAR*   src, 
                                 VECTOR_CHAR*   dest )
{

   if ( dest == NULL ) {
      dest = VECTOR_CHAR_Create();
   }
   /* allocate variable-sized data */
   VECTOR_CHAR_Resize( dest, src->Nalloc );
   /* copy variable-sized data */
   memcpy( dest->data, src->data, sizeof(CHAR) * src->N );
   /* copy base data */
   dest->N = src->N;

   return dest;
}

/*
 *  FUNCTION:  VECTOR_CHAR_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_CHAR_Resize(   VECTOR_CHAR*   vec, 
                           int            size )
{
   vec->data = (CHAR*) realloc( vec->data, sizeof(CHAR) * size );
   if ( vec->data == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n" );
      exit(EXIT_FAILURE);
   }
   vec->Nalloc = size;
}

/*
 *  FUNCTION:  VECTOR_CHAR_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_CHAR_GrowTo(   VECTOR_CHAR*   vec, 
                           int            size )
{
   if ( vec->Nalloc < size ) {
      VECTOR_CHAR_Resize( vec, size );
   }
}

/*
 *  FUNCTION:  VECTOR_CHAR_Push()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array. 
 *             Warning: Does not handle resizing or check for out-of-bounds. For that, use Pushback().
 */
inline
void VECTOR_CHAR_Push(   VECTOR_CHAR*   vec, 
                        CHAR           val )
{
   vec->data[vec->N] = val;
   vec->N++;
}

/*
 *  FUNCTION:  VECTOR_CHAR_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
inline
void VECTOR_CHAR_Pushback(  VECTOR_CHAR*   vec, 
                           CHAR           val )
{
   VECTOR_CHAR_Push( vec, val );

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_CHAR_Resize( vec, vec->N * 2 );
   }
}

/*
 *  FUNCTION:  VECTOR_CHAR_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
inline
CHAR VECTOR_CHAR_Pop( VECTOR_CHAR*   vec )
{
   CHAR data = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_CHAR_Resize( vec, vec->N / 2 );
   }

   return data;
}

/*
 *  FUNCTION:  VECTOR_CHAR_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_CHAR_Set(    VECTOR_CHAR*   vec, 
                        int           idx, 
                        CHAR           val )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_CHAR access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   vec->data[idx] = val;
}

/*
 *  FUNCTION:  VECTOR_CHAR_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return value of data.
 *  RETURN:    Return data at <idx>.
 */
inline
CHAR VECTOR_CHAR_Get(  VECTOR_CHAR*   vec, 
                     int           idx )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_CHAR access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return (vec->data[idx]);
}

/*
 *  FUNCTION:  VECTOR_CHAR_Get_X()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
inline
CHAR* VECTOR_CHAR_Get_X(  VECTOR_CHAR*   vec, 
                        int           idx )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_CHAR access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return &(vec->data[idx]);
}

/*
 *  FUNCTION:  VECTOR_CHAR_Get_Size()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
inline
int VECTOR_CHAR_Get_Size(   VECTOR_CHAR*   vec )
{
   return vec->N;
}

/*
 *  FUNCTION:  VECTOR_CHAR_Set_Size()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
inline
void VECTOR_CHAR_Set_Size(  VECTOR_CHAR*   vec, 
                           int           size )
{
   vec->N = size;
}


/*
 *  FUNCTION:  VECTOR_CHAR_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_CHAR_Search(  VECTOR_CHAR*   vec, 
                        CHAR           val )
{
   int N       = VECTOR_CHAR_Get_Size( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = CHAR_Compare( val, vec->data[i] );

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
 *  FUNCTION:  VECTOR_CHAR_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_CHAR_Search_First(  VECTOR_CHAR*   vec, 
                              CHAR           val )
{
   int N       = VECTOR_CHAR_Get_Size( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = CHAR_Compare( val, vec->data[i] );

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
 *  FUNCTION:  VECTOR_CHAR_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_CHAR_Search_Last(   VECTOR_CHAR*   vec, 
                              CHAR           val )
{
   int N       = VECTOR_CHAR_Get_Size( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = CHAR_Compare( val, vec->data[i] );

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
 *  FUNCTION:  VECTOR_CHAR_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_CHAR_Compare(    VECTOR_CHAR*   vec_A, 
                           VECTOR_CHAR*   vec_B )
{
   for (int i = 0; i < vec_A->N; i++) 
   {
      if ( CHAR_Compare( vec_A->data[i], vec_B->data[i] ) != 0 ) 
      {
         return CHAR_Compare( vec_A->data[i], vec_B->data[i] );       
      }
   }
   return 0;
}

/*
 *  FUNCTION:  VECTOR_CHAR_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
void VECTOR_CHAR_Sort( VECTOR_CHAR*    vec )
{
   int N = VECTOR_CHAR_Get_Size( vec );
   VECTOR_CHAR_Sort_Sub( vec, 0, N );

   #if DEBUG 
   {
      for ( int i = 0; i < N-1; i++ ) 
      {
         CHAR cur = vec->data[i];
         CHAR nxt = vec->data[i+1];
         char s_cur[50];
         char s_nxt[50];
         int cmp = CHAR_Compare( cur, nxt );
         if ( (cmp <= 0) == false ) {
            fprintf(stderr, "ERROR: bad sort. %d, %d v %d: %s vs %s\n",
               cmp, i, i+1, CHAR_To_String(cur, s_cur), CHAR_To_String(nxt, s_nxt) );
         }
      }
   }
   #endif
}

/*
 *  FUNCTION:  VECTOR_CHAR_Sort_Sub()
 *  SYNOPSIS:  Sorts subarray of <vec> data in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
void VECTOR_CHAR_Sort_Sub(  VECTOR_CHAR*    vec,
                           int            beg,
                           int            end )
{
   const int begin_select_sort = 4;
   int N = VECTOR_CHAR_Get_Size(vec);
   if (N <= 1) return;

   int size = end - beg;
   /* run selection sort if below threshold */
   if ( size <= begin_select_sort ) {
      VECTOR_CHAR_Sort_Sub_Selectsort( vec, beg, end );
   } 
   /* otherwise run quiksort */
   else {
      VECTOR_CHAR_Sort_Sub_Quicksort( vec, beg, end );
   }
}

/*
 *  FUNCTION:  VECTOR_CHAR_Sort_Sub_Selectsort()
 *  SYNOPSIS:  Selection Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
void VECTOR_CHAR_Sort_Sub_Selectsort(   VECTOR_CHAR*    vec,
                                       int            beg,
                                       int            end )
{
   /* find the minimum element of remaining unsorted list */
   for (int i = beg; i < end; i++) 
   {
      /* initial minimum value found */
      int min_idx = i;
      CHAR min_val = vec->data[i];
      for (int j = i+1; j < end; j++) {
         /* if new minimum found, update value and index */
         int cmp = CHAR_Compare( min_val, vec->data[j] );
         if ( cmp > 0 ) {
            min_idx = j;
            min_val = vec->data[j];
         }
      }
      /* swap new minimum to left-most position */
      VECTOR_CHAR_Swap( vec, i, min_idx );
   }
}

/*
 *  FUNCTION:  VECTOR_CHAR_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Quick Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
void VECTOR_CHAR_Sort_Sub_Quicksort( VECTOR_CHAR*    vec,
                                    int            beg,
                                    int            end )
{
   /* partition pointers */
   int   r_idx = beg + 1;
   int   l_idx = end - 1;
   CHAR*  rhs   = &(vec->data[beg + 1]);
   CHAR*  lhs   = &(vec->data[end - 1]);

   /* select random pivot value */
   int   range       = end - beg;
   int   pivot_idx   = RNG_INT_Range( beg, end );
   CHAR   pivot_val   = vec->data[pivot_idx];
   VECTOR_CHAR_Swap( vec, pivot_idx, beg );

   /* partition on pivot */
   while ( l_idx <= r_idx )
   {
      /* find next right partition element that is less than pivot element */
      while ( (l_idx <= r_idx) && (CHAR_Compare( pivot_val, vec->data[r_idx] ) < 0) ) {
         r_idx--;
      }
      /* find next left partition element that is greater than pivot element */
      while ( (l_idx <= r_idx) && (CHAR_Compare( pivot_val, vec->data[l_idx] ) >= 0) ) {
         l_idx++;
      }
      /* if left and right index have not crossed, then swap elements */
      if ( l_idx <= r_idx ) {
         VECTOR_CHAR_Swap( vec, l_idx, r_idx );
      }
   }
   /* move partition element to barrier between left and right index */
   VECTOR_CHAR_Swap( vec, beg, r_idx );
   /* sort both partitions (omit partition element) */
   VECTOR_CHAR_Sort_Sub( vec, beg, r_idx );
   VECTOR_CHAR_Sort_Sub( vec, r_idx+1, end );
}

/*
 *  FUNCTION:  VECTOR_CHAR_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
inline
void VECTOR_CHAR_Swap(   VECTOR_CHAR*    vec,
                        int            i,
                        int            j )
{
   CHAR swap = vec->data[i];
   vec->data[i] = vec->data[j];
   vec->data[j] = swap;
}

/*
 *  FUNCTION:  VECTOR_CHAR_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
inline
void VECTOR_CHAR_Reverse(   VECTOR_CHAR*    vec )
{
   int N = VECTOR_CHAR_Get_Size( vec );

   for (int i = 0; i < (N/2); i++) {
      VECTOR_CHAR_Swap( vec, i, (N-1)-i );
   }
}

/*
 *  FUNCTION:  VECTOR_CHAR_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
void VECTOR_CHAR_Dump(   VECTOR_CHAR*    vec,
                        FILE*          fp )
{
   /* stringification of template object */
   char s[50];

   fprintf(fp, "%s: ", "VECTOR_CHAR");
   fprintf(fp, "[ ");
   for ( int i = 0; i < vec->N; i++ ) {
      fprintf(fp, "%s, ", CHAR_To_String(vec->data[i], s) );
   }
   fprintf(fp, "]\n" );
}


/*
 *  FUNCTION:  VECTOR_CHAR_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_CHAR.
 */
void VECTOR_CHAR_Unit_Test()
{
   VECTOR_CHAR* vec = VECTOR_CHAR_Create();

   VECTOR_CHAR_Destroy( vec );
}