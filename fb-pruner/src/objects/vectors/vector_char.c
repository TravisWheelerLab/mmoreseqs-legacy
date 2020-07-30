/*******************************************************************************
 *  FILE:      vector_int.c
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
#include "vector_char.h"

/*
 *  FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Create new VECTOR_CHAR object and returns pointer.
 */
VECTOR_CHAR* VECTOR_CHAR_Create()
{
   const int init_size = 8;
   return VECTOR_CHAR_Create_by_Size( init_size );
}

/*
 *  FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Create new VECTOR_CHAR object at specific size and returns pointer.
 */
VECTOR_CHAR* VECTOR_CHAR_Create_by_Size( int size )
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

/* destructor */
void* VECTOR_CHAR_Destroy( VECTOR_CHAR* vec )
{
   if ( vec == NULL ) return vec;
   free(vec->data);
   // vec->data = NULL;

   free(vec);
   vec = NULL;
   return vec;
}

/* reuse by resetting counter*/
void VECTOR_CHAR_Reuse( VECTOR_CHAR* vec )
{
   vec->N = 0;
}

/* set all active indexes to zero */
void VECTOR_CHAR_Fill(  VECTOR_CHAR*   vec, 
                        char           val )
{
   for ( int i = 0; i < vec->N; i++ ) {
      vec->data[i] = val;
   }
}

/* deep copy */
VECTOR_CHAR* VECTOR_CHAR_Copy( VECTOR_CHAR* src )
{
   VECTOR_CHAR* vec = VECTOR_CHAR_Create();
   VECTOR_CHAR_Resize( vec, src->Nalloc );
   /* copy base data */
   memcpy( vec, src, sizeof(VECTOR_CHAR) );
   /* copy variable-sized data */
   memcpy( vec->data, src->data, sizeof(int) * src->N );

   return vec;
}

/* resize the array */
void VECTOR_CHAR_Resize( VECTOR_CHAR*  vec, 
                         int           size )
{
   vec->data = (char*) realloc( vec->data, sizeof(int) * size );
   if ( vec->data == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n" );
      exit(EXIT_FAILURE);
   }
   vec->Nalloc = size;
}

/* push element onto end of array */
void VECTOR_CHAR_Pushback(    VECTOR_CHAR*   vec, 
                              char           val )
{
   vec->data[vec->N] = val;
   vec->N++;

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_CHAR_Resize( vec, vec->N * 2 );
   }
}

/* pop element from end of array */
char VECTOR_CHAR_Pop( VECTOR_CHAR* vec )
{
   int tmp = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_CHAR_Resize( vec, vec->N / 2 );
   }

   return tmp;
}
/* get data at index (no checks) */
char* VECTOR_CHAR_Get(  VECTOR_CHAR*   vec, 
                        int            idx )
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

/* compare two VECTOR_CHAR objects */
int VECTOR_CHAR_Compare(   VECTOR_CHAR*   vecA, 
                           VECTOR_CHAR*   vecB )
{
   for (int i = 0; i < vecA->N; i++) {
      if ( vecA->data[i] != vecB->data[i] ) {
         if ( vecA->data[i] > vecB->data[i] ) {
            return 1;
         } else {
            return -1;
         }
      }
   }
   return 0;
}

/* output VECTOR_CHAR to file */
void VECTOR_CHAR_Dump(     VECTOR_CHAR*   vec,
                           FILE*          fp )
{
   fprintf(fp, "INTEGER VECTOR:\n");
   fprintf(fp, "[ ");
   for ( int i = 0; i < vec->N; i++ ) {
      fprintf(fp, "%3c ", vec->data[i] );
   }
   fprintf(fp, "]\n" );
}