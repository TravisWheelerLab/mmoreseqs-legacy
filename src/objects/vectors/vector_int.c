/*******************************************************************************
 *  @file VECTOR_INT.c
 *  @brief INT VECTOR Objects
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

/* local imports */
#include "objects/structs.h"

/* header */
#include "vector_int.h"

/* constructor */
VECTOR_INT* VECTOR_INT_Create()
{
   const int init_size = 8;
   VECTOR_INT *vec = NULL;
   vec         = (VECTOR_INT *) malloc( sizeof(VECTOR_INT) );
   if ( vec == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n");
      exit(EXIT_FAILURE);
   }
   vec->data   = NULL;
   vec->N      = 0;
   vec->Nalloc = 0;
   VECTOR_INT_Resize( vec, init_size );

   return vec;
}

/* destructor */
void VECTOR_INT_Destroy( VECTOR_INT *vec )
{
   free(vec->data);
   free(vec);
}

/* reuse */
void VECTOR_INT_Reuse( VECTOR_INT *vec )
{
   vec->N = 0;
}

/* deep copy */
VECTOR_INT* VECTOR_INT_Copy( VECTOR_INT *src )
{
   VECTOR_INT* vec = VECTOR_INT_Create();
   VECTOR_INT_Resize( vec, src->Nalloc );
   /* copy base data */
   memcpy( vec, src, sizeof(VECTOR_INT) );
   /* copy variable-sized data */
   memcpy( vec->data, src->data, sizeof(int) * src->N );

   return vec;
}

/* resize the array */
void VECTOR_INT_Resize( VECTOR_INT *vec, int size )
{
   vec->data = (int *) realloc( vec->data, sizeof(int) * size );
   if ( vec->data == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n" );
      exit(EXIT_FAILURE);
   }
   vec->Nalloc = size;
}

/* push element onto end of array */
void VECTOR_INT_Pushback( VECTOR_INT *vec, int val )
{
   vec->data[vec->N] = val;
   vec->N++;

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_INT_Resize( vec, vec->N * 2 );
   }
}

/* pop element from end of array */
int VECTOR_INT_Pop( VECTOR_INT *vec )
{
   int tmp = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_INT_Resize( vec, vec->N / 2 );
   }

   return tmp;
}

/* set data at index (no bound checks) */
void VECTOR_INT_Set( VECTOR_INT *vec, int idx, int val )
{
   vec->data[idx] = val;
}

/* get data at index (no checks) */
int VECTOR_INT_Get( VECTOR_INT *vec, int idx )
{
   return vec->data[idx];
}

/* compare two VECTOR_INT objects */
int VECTOR_INT_Compare( VECTOR_INT *vecA, VECTOR_INT *vecB )
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