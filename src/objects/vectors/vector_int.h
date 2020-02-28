/*******************************************************************************
 *  @file VECTOR_INT.c
 *  @brief INT VECTOR Objects
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _VECTOR_INT_H
#define _VECTOR_INT_H

/* import datatypes */
#include "../structs.h"

// /* Vector struct */
// typedef struct {
//    int *data;        /* array of data type */
//    int N;             current length of array in use 
//    int Nalloc;       /* current length of array allocated */
// }  VECTOR_INT;


/* constructor */
VECTOR_INT* VECTOR_INT_Create();

/* destructor */
void VECTOR_INT_Destroy( VECTOR_INT *vec );

/* deep copy */
VECTOR_INT* VECTOR_INT_Copy();

/* resize the array */
void VECTOR_INT_Resize( VECTOR_INT *vec, float growth_factor );

/* push element onto end of array */
void VECTOR_INT_Pushback( VECTOR_INT *vec, int val );

/* pop element from end of array */
int VECTOR_INT_Pop( VECTOR_INT *vec );

/* set data at index (no bound checks) */
void VECTOR_INT_Set( VECTOR_INT *vec, int idx, int val );

/* get data at index (no bound checks) */
int VECTOR_INT_Get( VECTOR_INT *vec, int idx );

/* equality test */
int VECTOR_INT_Compare( VECTOR_INT *vecA, VECTOR_INT *vecB );

#endif 