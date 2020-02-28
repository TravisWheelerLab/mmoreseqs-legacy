/*******************************************************************************
 *  @file VECTOR_FLOAT.c
 *  @brief CHARACTER VECTOR objects
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _VECTOR_FLOAT_H
#define _VECTOR_FLOAT_H

/* import datatypes */
#include "../structs.h"

/* Vector struct */
typedef struct {
   float *data;      /* array of data type */
   int N;            /* current length of array in use */
   int Nalloc;       /* current length of array allocated */
}  VECTOR_FLOAT;


/* constructor */
VECTOR_FLOAT* VECTOR_FLOAT_Create();

/* destructor */
void VECTOR_FLOAT_Destroy( VECTOR_FLOAT *vec );

/* deep copy */
VECTOR_FLOAT* VECTOR_FLOAT_Copy();

/* resize the array */
void VECTOR_FLOAT_Resize( VECTOR_FLOAT *vec, const float growth_factor );

/* push element onto end of array */
void VECTOR_FLOAT_Pushback( VECTOR_FLOAT *vec, const float val );

/* pop element from end of array */
float VECTOR_FLOAT_Pop( VECTOR_FLOAT *vec );

/* set data at index (no bound checks) */
void VECTOR_FLOAT_Set( VECTOR_FLOAT *vec, const int idx, const float val );

/* get data at index (no bound checks) */
float VECTOR_FLOAT_Get( VECTOR_FLOAT *vec, const int idx );

/* equality test */
int VECTOR_FLOAT_Compare( const VECTOR_FLOAT *vecA, const VECTOR_FLOAT *vecB );

#endif /* _VECTOR_FLOAT_H */