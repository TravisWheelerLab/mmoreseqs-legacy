/*******************************************************************************
 *  @file vector.c
 *  @brief OOP Vectors 
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* macros */
#define DATATYPE float
#define DEFAULT_SIZE 8

/* macros for concatenating tokens into single token */
#define M_CONC(A, B) M_CONC_(A, B)
#define M_CONC_(A, B) A ## B

#define FUNCTION_NAME(operation) M_CONC(STRUCT_NAME, M_CONC(_, operation))
#define STRUCT_NAME M_CONC(DATATYPE, Vec)

/* Vector struct */
typedef struct {
   DATATYPE *data;   /* array of data type */
   int N;            /* current length of array in use */
   int alloc;        /* current length of array allocated */
}  STRUCT_NAME;

/* constructor */
void FUNCTION_NAME(create) ( STRUCT_NAME *ptr )
{
   ptr = (STRUCT_NAME *) malloc( sizeof(STRUCT_NAME) );
   ptr->data = (DATATYPE *) malloc( sizeof(DATATYPE) * DEFAULT_SIZE );
   ptr->N = 0;
   ptr->alloc = DEFAULT_SIZE;
}

/* resize the array */
void FUNCTION_NAME(resize) ( STRUCT_NAME *ptr, int new_alloc )
{
   ptr->data = (DATATYPE *) realloc( ptr->data, sizeof(DATATYPE) * new_alloc );
   ptr->alloc = new_alloc;
}

/* push element onto end of array */
void FUNCTION_NAME(pushback) ( STRUCT_NAME *ptr, DATATYPE val )
{
   ptr->data[ptr->N] = val;
   ptr->N += 1;

   /* if array is full, resize */
   if (ptr->N >= ptr->alloc - 1) {
      FUNCTION_NAME(resize) ( ptr, ptr->alloc * 2 );
   }
}

/* pop element from end of array */
DATATYPE FUNCTION_NAME(pop) ( STRUCT_NAME *ptr, DATATYPE val )
{
   DATATYPE tmp = ptr->data[ptr->N-1];
   ptr->N -= 1;

   // /* if array is less than half used, resize */
   // if (ptr->N >= ptr->alloc / 2) {
   //    FUNCTION_NAME(resize) ( ptr, ptr->alloc / 2 );
   // }

   return tmp;
}

/* set data at index (no checks) */
void FUNCTION_NAME(set) ( STRUCT_NAME *ptr, DATATYPE val, int idx )
{
   ptr->data[idx] = val;
}

/* get data at index (no checks) */
DATATYPE FUNCTION_NAME(get) ( STRUCT_NAME *ptr, int idx )
{
   return ptr->data[idx];
}

/* destructor */
void FUNCTION_NAME(destroy) ( STRUCT_NAME *ptr )
{
   free(ptr->data);
   free(ptr);
}

// int main() 
// {
//    STRUCT_NAME *vec;

//    FUNCTION_NAME(create) (vec);

//    printf("VEC: N=%d, alloc=%d\n", vec->N, vec->alloc);

//    for (int i = 0; i < 32; ++i) {
//       printf("adding %d...\n", i);
//       FUNCTION_NAME(pushback) (vec, i);
//    }
// }