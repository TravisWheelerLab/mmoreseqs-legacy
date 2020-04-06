/*******************************************************************************
 *  @file VECTOR_CHAR.c
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
#include "objects/structs.h"

/* header */
#include "vector_char.h"

/* constructor */
VECTOR_CHAR* VECTOR_CHAR_Create()
{
   VECTOR_CHAR* vec;
   const int init_size = 8;
   vec = (VECTOR_CHAR *) malloc( sizeof(VECTOR_CHAR) );
   vec->data = (char *) malloc( sizeof(char) * init_size );
   vec->N = 0;
   vec->Nalloc = init_size;
   return vec;
}

/* destructor */
void VECTOR_CHAR_Destroy( VECTOR_CHAR*  vec )
{
   free(vec->data);
   free(vec);
}

/* deep copy */
VECTOR_CHAR* VECTOR_CHAR_Copy( VECTOR_CHAR*  src )
{
   VECTOR_CHAR* vec;
   vec = (VECTOR_CHAR*)malloc( sizeof(VECTOR_CHAR) );
   /* copy base data */
   memcpy( vec, src, sizeof(VECTOR_CHAR) );
   /* copy variable-sized data */
   vec->data = (char*) malloc( sizeof(char) * src->Nalloc );
   memcpy( vec->data, src->data, sizeof(char) * src->N );

   return vec;
}

/* resize the array */
void VECTOR_CHAR_Resize( VECTOR_CHAR*  vec, 
                         const float   growth_factor )
{
   vec->data = (char*) realloc( vec->data, sizeof(char) * vec->Nalloc * growth_factor );
   vec->Nalloc *= growth_factor;
}

/* push element onto end of array */
void VECTOR_CHAR_Pushback( VECTOR_CHAR*  vec, 
                           const char    val )
{
   vec->data[vec->N] = val;
   vec->N++;

   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_CHAR_Resize( vec, 2 );
   }
}

/* pop element from end of array */
char VECTOR_CHAR_Pop( VECTOR_CHAR*  vec )
{
   char tmp = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_CHAR_Resize( vec, 0.5 );
   }

   return tmp;
}

/* set data at index (no bound checks) */
void VECTOR_CHAR_Set( VECTOR_CHAR*  vec, 
                      const int     idx, 
                      const char    val )
{
   vec->data[idx] = val;
}

/* get data at index (no checks) */
char VECTOR_CHAR_Get( VECTOR_CHAR*  vec, 
                      const int     idx )
{
   return vec->data[idx];
}

/* compare two VECTOR_CHAR objects */
int VECTOR_CHAR_Compare( const VECTOR_CHAR*  vecA, 
                         const VECTOR_CHAR*  vecB )
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