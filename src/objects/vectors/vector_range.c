/*******************************************************************************
 *  @file VECTOR_RANGE.c
 *  @brief RANGE VECTORS Objects
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
#include "../structs.h"

/* header */
#include "vector_range.h"

/* constructor */
VECTOR_RANGE* VECTOR_RANGE_Create()
{
   VECTOR_RANGE *vec;
   const int init_size = 8;
   vec = (VECTOR_RANGE *) malloc( sizeof(VECTOR_RANGE) );
   vec->data = (RANGE *) malloc( sizeof(RANGE) * init_size );
   vec->N = 0;
   vec->Nalloc = init_size;
   return vec;
}

/* destructor */
void VECTOR_RANGE_Destroy( VECTOR_RANGE *vec )
{
   free(vec->data);
   free(vec);
}

/* deep copy */
VECTOR_RANGE* VECTOR_RANGE_Copy( VECTOR_RANGE *src )
{
   VECTOR_RANGE *vec;
   vec = (VECTOR_RANGE *) malloc( sizeof(VECTOR_RANGE) );
   /* copy base data */
   memcpy( vec, src, sizeof(VECTOR_RANGE) );
   /* copy variable-sized data */
   vec->data = (RANGE *) malloc( sizeof(RANGE) * src->Nalloc );
   memcpy( vec->data, src->data, sizeof(RANGE)  * src->N );

   return vec;
}

/* resize the array */
void VECTOR_RANGE_Resize( VECTOR_RANGE *vec, const float growth_factor )
{
   vec->data = (RANGE *) realloc( vec->data, sizeof(RANGE) * vec->Nalloc * growth_factor );
   vec->Nalloc *= growth_factor;
}

/* push element onto end of array */
void VECTOR_RANGE_Pushback( VECTOR_RANGE *vec, const RANGE val )
{
   vec->data[vec->N] = val;
   vec->N++;

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_RANGE_Resize( vec, 2 );
   }
}

/* pop element from end of array */
RANGE VECTOR_RANGE_Pop( VECTOR_RANGE *vec )
{
   RANGE tmp = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_RANGE_Resize( vec, 0.5 );
   }

   return tmp;
}

/* set data at index (no bound checks) */
void VECTOR_RANGE_Set( VECTOR_RANGE *vec, const int idx, const RANGE val )
{
   vec->data[idx] = val;
}

/* get data at index (no checks) */
RANGE VECTOR_RANGE_Get( VECTOR_RANGE *vec, const int idx )
{
   return vec->data[idx];
}

/* compare two VECTOR_RANGE objects */
int VECTOR_RANGE_Compare( const VECTOR_RANGE *vecA, const VECTOR_RANGE *vecB )
{
   // for (int i = 0; i < vecA->N; i++) {
   //    if ( vecA->data[i] != vecB->data[i] ) {
   //       if ( vecA->data[i] > vecB->data[i] ) {
   //          return 1;
   //       } else {
   //          return -1;
   //       }
   //    }
   // }
   return 0;
}

void VECTOR_RANGE_Dump(FILE *fp, VECTOR_RANGE *vec)
{
   for (int i = 0; i < vec->N; i++)
   {
      
   }
}