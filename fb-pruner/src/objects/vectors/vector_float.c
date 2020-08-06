/*******************************************************************************
 *  FILE:      vector_float.c
 *  PURPOSE:   VECTOR_FLT Object Functions
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
 *  FUNCTION:  VECTOR_FLT_Create()
 *  SYNOPSIS:  Create new VECTOR_FLT object and returns pointer.
 */
VECTOR_FLT* VECTOR_FLT_Create()
{
   const int init_size = VECTOR_INIT_SIZE;
   return VECTOR_FLT_Create_by_Size( init_size );
}

/*
 *  FUNCTION:  VECTOR_FLT_Create()
 *  SYNOPSIS:  Create new VECTOR_FLT object at specific size and returns pointer.
 */
VECTOR_FLT* VECTOR_FLT_Create_by_Size( int    size )
{
   VECTOR_FLT *vec = NULL;
   vec = (VECTOR_FLT *) malloc( sizeof(VECTOR_FLT) );
   if ( vec == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n");
      exit(EXIT_FAILURE);
   }

   vec->data   = NULL;
   vec->N      = 0;
   vec->Nalloc = 0;

   VECTOR_FLT_Resize( vec, size );

   return vec;
}

/*
 *  FUNCTION:  VECTOR_FLT_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_FLT.
 */
void* VECTOR_FLT_Destroy( VECTOR_FLT*   vec )
{
   if ( vec == NULL ) return vec;
   free(vec->data);
   free(vec);

   vec = NULL;
   return vec;
}

/*
 *  FUNCTION:  VECTOR_FLT_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_FLT object by resetting size counter (no realloc) .
 */
void VECTOR_FLT_Reuse( VECTOR_FLT*   vec )
{
   vec->N = 0;
}

/*
 *  FUNCTION:  VECTOR_FLT_Fill()
 *  SYNOPSIS:  Fill VECTOR_FLT object with val.
 */
void VECTOR_FLT_Fill(   VECTOR_FLT*   vec, 
                        FLT           val )
{
   for ( int i = 0; i < vec->N; i++ ) {
      vec->data[i] = val;
   }
}

/*
 *  FUNCTION:  VECTOR_FLT_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_FLT for <dest> if <dest> is NULL.
 */
VECTOR_FLT* VECTOR_FLT_Copy(  VECTOR_FLT*   src, 
                              VECTOR_FLT*   dest )
{

   if ( dest == NULL ) {
      dest = VECTOR_FLT_Create();
   }
   /* allocate variable-sized data */
   VECTOR_FLT_Resize( dest, src->Nalloc );
   /* copy variable-sized data */
   memcpy( dest->data, src->data, sizeof(FLT) * src->N );
   /* copy base data */
   dest->N = src->N;

   return dest;
}

/*
 *  FUNCTION:  VECTOR_FLT_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_FLT_Resize(    VECTOR_FLT*   vec, 
                           int           size )
{
   vec->data = (FLT*) realloc( vec->data, sizeof(FLT) * size );
   if ( vec->data == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n" );
      exit(EXIT_FAILURE);
   }
   vec->Nalloc = size;
}

/*
 *  FUNCTION:  VECTOR_FLT_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_FLT_GrowTo( VECTOR_FLT*   vec, 
                        int           size )
{
   if ( vec->Nalloc < size ) {
      VECTOR_FLT_Resize( vec, size );
   }
}

/*
 *  FUNCTION:  VECTOR_FLT_Push()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array. 
 *             Warning: Does not handle resizing or check for out-of-bounds. For that, use Pushback().
 */
inline
void VECTOR_FLT_Push(   VECTOR_FLT*   vec, 
                        FLT           val )
{
   vec->data[vec->N] = val;
   vec->N++;
}

/*
 *  FUNCTION:  VECTOR_FLT_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
inline
void VECTOR_FLT_Pushback(  VECTOR_FLT*   vec, 
                           FLT           val )
{
   VECTOR_FLT_Push( vec, val );

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_FLT_Resize( vec, vec->N * 2 );
   }
}

/*
 *  FUNCTION:  VECTOR_FLT_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
inline
FLT VECTOR_FLT_Pop( VECTOR_FLT*   vec )
{
   FLT data = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_FLT_Resize( vec, vec->N / 2 );
   }

   return data;
}

/*
 *  FUNCTION:  VECTOR_FLT_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_FLT_Set(    VECTOR_FLT*   vec, 
                        int           idx, 
                        FLT           val )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_FLT access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   vec->data[idx] = val;
}

/*
 *  FUNCTION:  VECTOR_FLT_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return value of data.
 *  RETURN:    Return data at <idx>.
 */
inline
FLT VECTOR_FLT_Get(  VECTOR_FLT*   vec, 
                     int           idx )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_FLT access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return &(vec->data[idx]);
}

/*
 *  FUNCTION:  VECTOR_FLT_Get_X()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
inline
FLT* VECTOR_FLT_Get_X(  VECTOR_FLT*   vec, 
                        int           idx )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_FLT access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return &(vec->data[idx]);
}

/*
 *  FUNCTION:  VECTOR_FLT_Get_Size()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
inline
FLT VECTOR_FLT_Get_Size(   VECTOR_FLT*   vec )
{
   return vec->N;
}

/*
 *  FUNCTION:  VECTOR_FLT_Set_Size()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
inline
void VECTOR_FLT_Set_Size(  VECTOR_FLT*   vec, 
                           int           size )
{
   vec->N = size;
}


/*
 *  FUNCTION:  VECTOR_FLT_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_FLT_Search(  VECTOR_FLT*   vec, 
                        FLT           val )
{
   int N       = VECTOR_FLT_Get_Size( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = FLT_Compare( val, vec->data[i] );

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
 *  FUNCTION:  VECTOR_FLT_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_FLT_Search_First(  VECTOR_FLT*   vec, 
                              FLT           val )
{
   int N       = VECTOR_FLT_Get_Size( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = FLT_Compare( val, vec->data[i] );

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
 *  FUNCTION:  VECTOR_FLT_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_FLT_Search_Last(   VECTOR_FLT*   vec, 
                              FLT           val )
{
   int N       = VECTOR_FLT_Get_Size( vec );
   int idx     = (N/2);
   int found   = -1; 

   for (int i = N/4; i >= 1; i /= 2)
   {
      int cmp = FLT_Compare( val, vec->data[i] );

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
 *  FUNCTION:  VECTOR_FLT_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_FLT_Compare(    VECTOR_FLT*   vec_A, 
                           VECTOR_FLT*   vec_B )
{
   for (int i = 0; i < vec_A->N; i++) 
   {
      if ( FLT_Compare( vec_A->data[i], vec_B->data[i] ) != 0 ) 
      {
         return FLT_Compare( vec_A->data[i], vec_B->data[i] ) != 0 );       
      }
   }
   return 0;
}

/*
 *  FUNCTION:  VECTOR_FLT_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
void VECTOR_FLT_Sort( VECTOR_FLT*    vec )
{
   int N = VECTOR_FLT_Get_Size(vec);
   VECTOR_FLT_Sort_Sub( vec, 0, N );
}

/*
 *  FUNCTION:  VECTOR_FLT_Sort_Sub()
 *  SYNOPSIS:  Sorts subarray of <vec> data in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
void VECTOR_FLT_Sort_Sub(  VECTOR_FLT*    vec,
                           int            beg,
                           int            end )
{
   const int begin_select_sort = 16;
   
   /* run selection sort if below threshold */
   if ( end - beg > begin_select_sort ) {
      VECTOR_FLT_Sort_Sub_Selectsort( vec, beg, end );
   }
   VECTOR_FLT_Sort_Sub_Quicksort( vec, beg, end );
}

/*
 *  FUNCTION:  VECTOR_FLT_Sort_Sub_Selectsort()
 *  SYNOPSIS:  Selection Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
void VECTOR_FLT_Sort_Sub_Selectsort(   VECTOR_FLT*    vec,
                                       int            beg,
                                       int            end )
{
   for (int i = beg; i < end; i++) 
   {
      /* initial minimum value found */
      int min_idx = i;
      FLT min_val = vec->data[i];
      for (int j = i+1; j < end; j++) {
         /* if new minimum found, update value and index */
         int cmp = FLT_Compare( min_val, vec->data[j] );
         if ( cmp < 0 ) {
            min_idx = j;
            min_val = vec->data[j];
         }
      }
      /* swap new minimum to left-most position */
      VECTOR_FLT_Swap( vec, i, min_idx );
   }
   return;
}

/*
 *  FUNCTION:  VECTOR_FLT_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Quick Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
void VECTOR_FLT_Sort_Sub_Quicksort( VECTOR_FLT*    vec,
                                    int            beg,
                                    int            end )
{
   /* partition pointers */
   int   r_idx = beg + 1;
   int   l_idx = end - 1;
   FLT*  rhs   = &(vec->data[beg + 1]);
   FLT*  lhs   = &(vec->data[end - 1]);

   /* select random pivot value */
   int   range = end - beg;
   int pivot_idx = (rand() % range) + beg;
   FLT pivot_val = vec->data[pivot_idx];
   VECTOR_FLT_Swap( vec, pivot, beg );

   /* partition on pivot */
   while ( l_idx <= r_idx )
   {
      while ( (l_idx <= r_idx) && (FLT_Compare( pivot_val, vec->data[r_idx] ) < 0) ) {
         r_idx--;
      }
      while ( (l_idx <= r_idx) && (FLT_Compare( pivot_val, vec->data[l_idx] ) >= 0) ) {
         l_idx++;
      }
      if ( l_idx <= r_idx ) {
         VECTOR_FLT_Swap( vec, l_idx, r_idx );
         r_idx--;
         l_idx++;
      }
   }
   VECTOR_FLT_Swap( vec, beg, r_idx );
}

/*
 *  FUNCTION:  VECTOR_FLT_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
inline
void VECTOR_FLT_Swap(   VECTOR_FLT*    vec,
                        int            i,
                        int            j )
{
   FLT swap = vec->data[i];
   vec->data[i] = vec->data[j];
   vec->data[j] = swap;
}

/*
 *  FUNCTION:  VECTOR_FLT_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
inline
void VECTOR_FLT_Reverse(   VECTOR_FLT*    vec )
{
   int N = VECTOR_FLT_Get_Size( vec );

   for (int i = 0; i < (N/2); i++) {
      VECTOR_FLT_Swap( vec, i, (N-1)-i );
   }
}

/*
 *  FUNCTION:  VECTOR_FLT_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_FLT_Dump(   VECTOR_FLT*    vec,
                        FILE*          fp )
{
   fprintf(fp, "%s: ", "VECTOR_FLT");
   fprintf(fp, "[ ");
   for ( int i = 0; i < vec->N; i++ ) {
      fprintf(fp, "%3d ", vec->data[i] );
   }
   fprintf(fp, "]\n" );
}
