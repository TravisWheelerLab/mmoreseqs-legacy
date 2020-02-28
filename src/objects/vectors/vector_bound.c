/*******************************************************************************
 *  @file VECTOR_BOUND.c
 *  @brief VECTOR BOUND objects
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
#include "vector_bound.h"

/* constructor */
VECTOR_BOUND* VECTOR_BOUND_Create()
{
   VECTOR_BOUND *vec;
   const int init_size = 8;
   vec = (VECTOR_BOUND *) malloc( sizeof(VECTOR_BOUND) );
   if (vec == NULL) {
      perror("ERROR: Error while malloc'ing VECTOR_BOUND.\n");
      exit(0);
   }
   vec->data = (BOUND *) malloc( sizeof(BOUND) * init_size );
   if (vec == NULL) {
      perror("ERROR: Error while malloc'ing VECTOR_BOUND.\n");
      exit(0);
   }
   vec->N = 0;
   vec->Nalloc = init_size;
   return vec;
}

/* destructor */
void VECTOR_BOUND_Destroy( VECTOR_BOUND  *vec )
{
   free(vec->data);
   free(vec);
}

/* deep copy */
VECTOR_BOUND* VECTOR_BOUND_Copy( VECTOR_BOUND  *src )
{
   VECTOR_BOUND *vec;
   vec = (VECTOR_BOUND *) malloc( sizeof(VECTOR_BOUND) );
   if (vec == NULL) {
      perror("ERROR: Unable to malloc VECTOR_BOUND.\n");
   }
   /* copy base data */
   memcpy( vec, src, sizeof(VECTOR_BOUND) );
   /* copy variable-sized data */
   vec->data = (BOUND *) malloc( sizeof(BOUND) * src->Nalloc );
   memcpy( vec->data, src->data, sizeof(BOUND) * src->N );

   return vec;
}

/* resize the array */
void VECTOR_BOUND_Resize( VECTOR_BOUND  *vec, 
                          float         growth_factor )
{
   vec->data = (BOUND *) realloc( vec->data, sizeof(BOUND) * vec->Nalloc * growth_factor );
   vec->Nalloc *= growth_factor;
}

/* push element onto end of array */
void VECTOR_BOUND_Pushback( VECTOR_BOUND *vec, BOUND val )
{
   vec->data[vec->N] = val;
   vec->N++;

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_BOUND_Resize( vec, 2 );
   }
}

/* pop element from end of array */
BOUND VECTOR_BOUND_Pop( VECTOR_BOUND  *vec )
{
   BOUND tmp = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_BOUND_Resize( vec, 0.5 );
   }

   return tmp;
}

/* set data at index (no bound checks) */
void VECTOR_BOUND_Set( VECTOR_BOUND  *vec, 
                       int           idx, 
                       BOUND         val )
{
   vec->data[idx] = val;
}

/* get data at index (no checks) */
BOUND VECTOR_BOUND_Get( VECTOR_BOUND  *vec, 
                        int           idx )
{
   return vec->data[idx];
}

/* compare two VECTOR_BOUND objects */
int VECTOR_BOUND_Compare( VECTOR_BOUND  *vecA, 
                          VECTOR_BOUND  *vecB )
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