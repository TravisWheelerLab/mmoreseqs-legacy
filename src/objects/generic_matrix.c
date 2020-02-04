/*******************************************************************************
 *  @file vector.c
 *  @brief Generic OOP Generic N-D Matrices 
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
#include <math.h>
#include <ctype.h>

/* macros */
#define DATATYPE      float
#define DEFAULT_SIZE  8
#define GROWTH_RATE   2

/* macros for concatenating tokens into single token */
#define M_CONC(A, B) M_CONC_(A, B)
#define M_CONC_(A, B) A ## B

/* generate function name tokens by data type and method */
#define FUNCTION_NAME(operation) M_CONC(STRUCT_NAME, M_CONC(_, operation))
#define STRUCT_NAME M_CONC(DATATYPE, Matrix)

/* Vector struct */
typedef struct {
   int N_dim;           /* number of dimensions */
   int *dims;           /* dimension sizes */
   DATATYPE *data;      /* array of data type */
}  STRUCT_NAME;

/* constructor */
STRUCT_NAME *FUNCTION_NAME(create) ( int num_dim, int *dims )
{
   STRUCT_NAME *ptr;
   ptr = (STRUCT_NAME *) malloc( sizeof(STRUCT_NAME) );
   ptr->N_dim = num_dim;
   ptr->dims = (int *) malloc( sizeof(int) * num_dim );

   int area = 1;
   for (int i = 0; i < num_dim; i++) {
      area *= dims[i];
   }

   ptr->data = (DATATYPE *) malloc( sizeof(DATATYPE) * area );
}

/* destructor */
void FUNCTION_NAME(destroy) ( STRUCT_NAME *ptr )
{
   free(ptr->data);
   free(ptr->dims);
   free(ptr);
}