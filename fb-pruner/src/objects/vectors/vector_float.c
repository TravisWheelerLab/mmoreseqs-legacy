/*******************************************************************************
 *  @file VECTOR_FLT.c
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
VECTOR_FLT* VECTOR_FLT_Create()
{
   const int init_size = VECTOR_INIT_SIZE;

   VECTOR_FLT* vec;
   vec = (VECTOR_FLT*) malloc( sizeof(VECTOR_FLT) );
   if ( vec == NULL ) {
      fprintf(stderr, "ERROR: Memory Allocation error.\n");
      exit(EXIT_FAILURE);
   }

   vec->N = 0;
   vec->Nalloc = init_size;
   VECTOR_FLT_Resize( vec, init_size );

   return vec;
}

/* destructor */
void VECTOR_FLT_Destroy( VECTOR_FLT*   vec )
{
   free(vec->data);
   free(vec);
}

/* empty vector */
void VECTOR_FLT_Reuse( VECTOR_FLT*   vec )
{
   vec->N = 0;
}

/* set all active indexes to val */
void VECTOR_FLT_Fill(  VECTOR_FLT*   vec, 
                       FLT           val )
{
   for ( int i = 0; i < vec->N; i++ ) {
      vec->data[i] = val;
   }
}

/* deep copy */
VECTOR_FLT* VECTOR_FLT_Copy(  VECTOR_FLT*    dest,
                              VECTOR_FLT*    src )
{
   if ( dest == NULL ) {
      dest = (VECTOR_FLT*) malloc( sizeof(VECTOR_FLT) );
   }
   
   /* copy base data */
   dest->N = src->N;
   dest->Nalloc = src->Nalloc;

   /* copy variable-sized data */
   dest->data = (FLT*) malloc( sizeof(FLT) * src->Nalloc );
   memcpy( dest->data, src->data, sizeof(FLT) * src->N );

   return dest;
}

/* resize the array */
void VECTOR_FLT_Resize( VECTOR_FLT*    vec, 
                        const int      size )
{
   if ( size > vec->Nalloc ) {
      vec->data = (FLT*) realloc( vec->data, sizeof(FLT) * size );
      vec->Nalloc = size;
   }
}

/* push element onto end of array */
void VECTOR_FLT_Pushback(  VECTOR_FLT*    vec, 
                           const FLT      val )
{
   vec->data[vec->N] = val;
   vec->N++;

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_FLT_Resize( vec, 2 );
   }
}

/* pop element from end of array */
FLT VECTOR_FLT_Pop( VECTOR_FLT*   vec )
{
   float tmp = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_FLT_Resize( vec, 0.5 );
   }

   return tmp;
}

/* set data at index (no bound checks) */
void VECTOR_FLT_Set( VECTOR_FLT*    vec, 
                     const int      idx, 
                     const float    val )
{
   vec->data[idx] = val;
}

/* get data at index (no checks) */
FLT VECTOR_FLT_Get(  VECTOR_FLT*    vec, 
                     const int      idx )
{
   return vec->data[idx];
}

/* compare two VECTOR_FLT objects */
int VECTOR_FLT_Compare( const VECTOR_FLT*    vec_A, 
                        const VECTOR_FLT*    vec_B )
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