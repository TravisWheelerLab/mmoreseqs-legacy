/*******************************************************************************
 *  FILE:      vector_template.c
 *  PURPOSE:   VECTOR_TMP Object Functions
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
 *  FUNCTION:  VECTOR_TMP_Create()
 *  SYNOPSIS:  Create new VECTOR_TMP object and returns pointer.
 */
VECTOR_TMP* VECTOR_TMP_Create()
{
   const int init_size = VECTOR_INIT_SIZE;
   return VECTOR_TMP_Create_by_Size( init_size );
}

/*
 *  FUNCTION:  VECTOR_TMP_Create()
 *  SYNOPSIS:  Create new VECTOR_TMP object at specific size and returns pointer.
 */
VECTOR_TMP* VECTOR_TMP_Create_by_Size( int    size )
{
   VECTOR_TMP *vec = NULL;
   vec         = (VECTOR_TMP *) malloc( sizeof(VECTOR_TMP) );
   if ( vec == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n");
      exit(EXIT_FAILURE);
   }

   vec->data   = NULL;
   vec->N      = 0;
   vec->Nalloc = 0;

   VECTOR_TMP_Resize( vec, size );

   return vec;
}

/*
 *  FUNCTION:  VECTOR_TMP_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_TMP.
 */
void* VECTOR_TMP_Destroy( VECTOR_TMP*   vec )
{
   if ( vec == NULL ) return vec;
   free(vec->data);
   free(vec);

   vec = NULL;
   return vec;
}

/*
 *  FUNCTION:  VECTOR_TMP_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_TMP object by resetting size counter (no realloc) .
 */
void VECTOR_TMP_Reuse( VECTOR_TMP*   vec )
{
   vec->N = 0;
}

/*
 *  FUNCTION:  VECTOR_TMP_Fill()
 *  SYNOPSIS:  Fill VECTOR_TMP object with val.
 */
void VECTOR_TMP_Fill(   VECTOR_TMP*   vec, 
                        TMP           val )
{
   for ( int i = 0; i < vec->N; i++ ) {
      vec->data[i] = val;
   }
}

/*
 *  FUNCTION:  VECTOR_TMP_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_TMP for <dest> if <dest> is NULL.
 */
VECTOR_TMP* VECTOR_TMP_Copy(  VECTOR_TMP*   src, 
                              VECTOR_TMP*   dest )
{

   if ( dest == NULL ) {
      dest = VECTOR_TMP_Create();
   }
   /* allocate variable-sized data */
   VECTOR_TMP_Resize( dest, src->Nalloc );
   /* copy variable-sized data */
   memcpy( dest->data, src->data, sizeof(TMP) * src->N );
   /* copy base data */
   dest->N = src->N;

   return dest;
}

/*
 *  FUNCTION:  VECTOR_TMP_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_TMP_Resize(    VECTOR_TMP*   vec, 
                           int           size )
{
   vec->data = (TMP*) realloc( vec->data, sizeof(TMP) * size );
   if ( vec->data == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n" );
      exit(EXIT_FAILURE);
   }
   vec->Nalloc = size;
}

/*
 *  FUNCTION:  VECTOR_TMP_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_TMP_GrowTo( VECTOR_TMP*   vec, 
                        int           size )
{
   if ( vec->Nalloc < size ) {
      VECTOR_TMP_Resize( vec, size );
   }
}

/*
 *  FUNCTION:  VECTOR_TMP_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_TMP_Pushback(  VECTOR_TMP*   vec, 
                           TMP           val )
{
   vec->data[vec->N] = val;
   vec->N++;

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_TMP_Resize( vec, vec->N * 2 );
   }
}

/*
 *  FUNCTION:  VECTOR_TMP_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
inline
TMP VECTOR_TMP_Pop( VECTOR_TMP*   vec )
{
   TMP data = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_TMP_Resize( vec, vec->N / 2 );
   }

   return data;
}

/*
 *  FUNCTION:  VECTOR_TMP_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_TMP_Set(    VECTOR_TMP*   vec, 
                        int           idx, 
                        TMP           val )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_TMP access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   vec->data[idx] = val;
}

/*
 *  FUNCTION:  VECTOR_TMP_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
inline
TMP VECTOR_TMP_Get(  VECTOR_TMP*   vec, 
                     int           idx )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_TMP access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return vec->data[idx];
}


/*
 *  FUNCTION:  VECTOR_TMP_Get_Ref()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
inline
TMP* VECTOR_TMP_Get_Ref(   VECTOR_TMP*   vec, 
                           int           idx )
{
   /* if debugging, do edgebound checks */
   #if DEBUG 
   if ( idx >= vec->N || idx < 0 ) {
      fprintf(stderr, "ERROR: VECTOR_TMP access out-of-bounds.\n");
      fprintf(stderr, "dim: (%d/%d), access: %d\n", vec->N, vec->Nalloc, idx);
      exit(EXIT_FAILURE);
   }
   #endif

   return &(vec->data[idx]);
}


/*
 *  FUNCTION:  VECTOR_TMP_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_TMP_Compare(    VECTOR_TMP*   vec_A, 
                           VECTOR_TMP*   vec_B )
{
   for (int i = 0; i < vec_A->N; i++) {
      if ( vec_A->data[i] != vec_B->data[i] ) {
         if ( vec_A->data[i] > vec_B->data[i] ) {
            return 1;
         } else {
            return -1;
         }
      }
   }
   return 0;
}


/*
 *  FUNCTION:  VECTOR_TMP_Sort()
 *  SYNOPSIS:  Sort <vec> data array in ascending order.
 */
void VECTOR_TMP_Sort( VECTOR_TMP*    vec )
{
   /* TODO */
}


/*
 *  FUNCTION:  VECTOR_TMP_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_TMP_Dump(   VECTOR_TMP*    vec,
                        FILE*          fp )
{
   fprintf(fp, "%s: ", "VECTOR_TMP");
   fprintf(fp, "[ ");
   for ( int i = 0; i < vec->N; i++ ) {
      fprintf(fp, "%3d ", vec->data[i] );
   }
   fprintf(fp, "]\n" );
}