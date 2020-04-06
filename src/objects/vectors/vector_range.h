/*******************************************************************************
 *  @file VECTOR_RANGE.c
 *  @brief RANGE VECTORS Objects
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _VECTOR_RANGE_H
#define _VECTOR_RANGE_H

/* import datatypes */
#include "objects/structs.h"

/* VECTOR struct */
typedef struct {
   RANGE *data;      /* array of data type */
   int N;            /* current length of array in use */
   int Nalloc;       /* current length of array allocated */
}  VECTOR_RANGE;


/* constructor */
VECTOR_RANGE* VECTOR_RANGE_Create();

/* destructor */
void VECTOR_RANGE_Destroy( VECTOR_RANGE *vec );

/* deep copy */
VECTOR_RANGE* VECTOR_RANGE_Copy();

/* resize the array */
void VECTOR_RANGE_Resize( VECTOR_RANGE *vec, const float growth_factor );

/* push element onto end of array */
void VECTOR_RANGE_Pushback( VECTOR_RANGE *vec, const RANGE val );

/* pop element from end of array */
RANGE VECTOR_RANGE_Pop( VECTOR_RANGE *vec );

/* set data at index (no bound checks) */
void VECTOR_RANGE_Set( VECTOR_RANGE *vec, const int idx, const RANGE val );

/* get data at index (no bound checks) */
RANGE VECTOR_RANGE_Get( VECTOR_RANGE *vec, const int idx );

/* equality test */
int VECTOR_RANGE_Compare( const VECTOR_RANGE *vecA, const VECTOR_RANGE *vecB );

#endif /* _VECTOR_RANGE_H */