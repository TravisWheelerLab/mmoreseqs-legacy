/*******************************************************************************
 *  FILE:      vector_char.c
 *  PURPOSE:   VECTOR_CHAR Object Functions
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
 *  FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Create new VECTOR_CHAR object and returns pointer.
 */
VECTOR_CHAR* VECTOR_CHAR_Create()
{
   const int init_size = VECTOR_INIT_SIZE;
   return VECTOR_CHAR_Create_by_Size( init_size );
}

/*
 *  FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Create new VECTOR_CHAR object at specific size and returns pointer.
 */
VECTOR_CHAR* VECTOR_CHAR_Create_by_Size( int    size )
{
   VECTOR_CHAR *vec = NULL;
   vec         = (VECTOR_CHAR *) malloc( sizeof(VECTOR_CHAR) );
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
   if ( vec == NULL ) return vec;
   free(vec->data);
   free(vec);

   vec = NULL;
   return vec;
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

/*
 *  FUNCTION:  VECTOR_CHAR_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_CHAR for <dest> if <dest> is NULL.
 */
VECTOR_CHAR* VECTOR_CHAR_Copy(  VECTOR_CHAR*   src, 
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
void VECTOR_CHAR_Resize(    VECTOR_CHAR*   vec, 
                           int           size )
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
void VECTOR_CHAR_GrowTo( VECTOR_CHAR*   vec, 
                        int           size )
{
   if ( vec->Nalloc < size ) {
      VECTOR_CHAR_Resize( vec, size );
   }
}

/*
 *  FUNCTION:  VECTOR_CHAR_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_CHAR_Pushback(  VECTOR_CHAR*   vec, 
                           CHAR           val )
{
   vec->data[vec->N] = val;
   vec->N++;

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
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
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

   return vec->data[idx];
}


/*
 *  FUNCTION:  VECTOR_CHAR_Get_Ref()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
inline
CHAR* VECTOR_CHAR_Get_Ref(   VECTOR_CHAR*   vec, 
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
 *  FUNCTION:  VECTOR_CHAR_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_CHAR_Compare(    VECTOR_CHAR*   vec_A, 
                           VECTOR_CHAR*   vec_B )
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
 *  FUNCTION:  VECTOR_CHAR_Sort()
 *  SYNOPSIS:  Sort <vec> data array in ascending order.
 */
void VECTOR_CHAR_Sort( VECTOR_CHAR*    vec )
{
   /* TODO */
}


/*
 *  FUNCTION:  VECTOR_CHAR_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_CHAR_Dump(   VECTOR_CHAR*    vec,
                        FILE*          fp )
{
   fprintf(fp, "%s: ", "VECTOR_CHAR");
   fprintf(fp, "[ ");
   for ( int i = 0; i < vec->N; i++ ) {
      fprintf(fp, "%3d ", vec->data[i] );
   }
   fprintf(fp, "]\n" );
}