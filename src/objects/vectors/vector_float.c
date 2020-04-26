/*******************************************************************************
 *  @file VECTOR_FLOAT.c
 *  @brief FLOAT VECTOR objects
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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "vector_float.h"

/* constructor */
VECTOR_FLOAT* VECTOR_FLOAT_Create()
{
   VECTOR_FLOAT *vec;
   const int init_size = 8;
   vec = (VECTOR_FLOAT *) malloc( sizeof(VECTOR_FLOAT) );
   vec->data = (float *) malloc( sizeof(float) * init_size );
   vec->N = 0;
   vec->Nalloc = init_size;
   return vec;
}

/* destructor */
void VECTOR_FLOAT_Destroy( VECTOR_FLOAT *vec )
{
   free(vec->data);
   free(vec);
}

/* deep copy */
VECTOR_FLOAT* VECTOR_FLOAT_Copy( VECTOR_FLOAT *src )
{
   VECTOR_FLOAT *vec;
   vec = (VECTOR_FLOAT *) malloc( sizeof(VECTOR_FLOAT) );
   /* copy base data */
   memcpy( vec, src, sizeof(VECTOR_FLOAT) );
   /* copy variable-sized data */
   vec->data = (float *) malloc( sizeof(float) * src->Nalloc );
   memcpy( vec->data, src->data, sizeof(float) * src->N );

   return vec;
}

/* resize the array */
void VECTOR_FLOAT_Resize( VECTOR_FLOAT *vec, const float growth_factor )
{
   vec->data = (float *) realloc( vec->data, sizeof(float) * vec->Nalloc * growth_factor );
   vec->Nalloc *= growth_factor;
}

/* push element onto end of array */
void VECTOR_FLOAT_Pushback( VECTOR_FLOAT *vec, const float val )
{
   vec->data[vec->N] = val;
   vec->N++;

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_FLOAT_Resize( vec, 2 );
   }
}

/* pop element from end of array */
float VECTOR_FLOAT_Pop( VECTOR_FLOAT *vec )
{
   float tmp = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_FLOAT_Resize( vec, 0.5 );
   }

   return tmp;
}

/* set data at index (no bound checks) */
void VECTOR_FLOAT_Set( VECTOR_FLOAT *vec, const int idx, const float val )
{
   vec->data[idx] = val;
}

/* get data at index (no checks) */
float VECTOR_FLOAT_Get( VECTOR_FLOAT *vec, const int idx )
{
   return vec->data[idx];
}

/* compare two VECTOR_FLOAT objects */
int VECTOR_FLOAT_Compare( const VECTOR_FLOAT *vecA, const VECTOR_FLOAT *vecB )
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