/*******************************************************************************
 *  FILE:      vector_trace.c
 *  PURPOSE:   VECTOR_TRACE objects
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
#include "vector_trace.h"

/* constructor */
VECTOR_TRACE* VECTOR_TRACE_Create()
{
   const int init_size = VECTOR_INIT_SIZE;   
   return VECTOR_TRACE_Create_by_Size( init_size );
}

/* constructor */
VECTOR_TRACE* VECTOR_TRACE_Create_by_Size( int size )
{
   VECTOR_TRACE*  vec = NULL;
   
   vec = (VECTOR_TRACE *) malloc( sizeof(VECTOR_TRACE) );
   if ( vec == NULL ) {
      fprintf(stderr, "ERROR: Failure to malloc.\n");
      exit(EXIT_FAILURE);
   }

   vec->N      = 0;
   vec->Nalloc = 0;
   vec->data   = NULL;

   VECTOR_TRACE_Resize( vec, size );

   return vec;
}

/* destructor */
void VECTOR_TRACE_Destroy( VECTOR_TRACE *vec )
{
   free(vec->data);
   free(vec);
}

/* deep copy */
VECTOR_TRACE* VECTOR_TRACE_Copy( VECTOR_TRACE *src )
{
   VECTOR_TRACE *vec;
   vec = (VECTOR_TRACE *) malloc( sizeof(VECTOR_TRACE) );
   /* copy base data */
   memcpy( vec, src, sizeof(VECTOR_TRACE) );
   /* copy variable-sized data */
   vec->data = (TRACE *) malloc( sizeof(TRACE) * src->Nalloc );
   memcpy( vec->data, src->data, sizeof(TRACE) * src->N );

   return vec;
}

/* resize the array */
void VECTOR_TRACE_Resize(  VECTOR_TRACE*  vec, 
                           int            size )
{
   vec->data = (TRACE*) realloc( vec->data, sizeof(TRACE) * size );
   vec->Nalloc = size;
}

/* resize the array */
void VECTOR_TRACE_GrowTo(  VECTOR_TRACE*  vec, 
                           int            size )
{
   if ( size > vec->N ) {
      VECTOR_TRACE_Resize( vec, size );
   }
}

/* push element onto end of array */
void VECTOR_TRACE_Pushback( VECTOR_TRACE *vec, TRACE val )
{
   vec->data[vec->N] = val;
   vec->N++;

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_TRACE_Resize( vec, 2 );
   }
}

/* pop element from end of array */
TRACE VECTOR_TRACE_Pop( VECTOR_TRACE *vec )
{
   TRACE tmp = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_TRACE_Resize( vec, 0.5 );
   }

   return tmp;
}

/* set data at index (no bound checks) */
void VECTOR_TRACE_Set( VECTOR_TRACE *vec, int idx, TRACE val )
{
   vec->data[idx] = val;
}

/* get data at index (no checks) */
TRACE VECTOR_TRACE_Get( VECTOR_TRACE *vec, int idx )
{
   return vec->data[idx];
}

/* compare two VECTOR_TRACE objects */
int VECTOR_TRACE_Compare( VECTOR_TRACE *vecA, VECTOR_TRACE *vecB )
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