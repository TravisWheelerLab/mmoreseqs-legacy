/*******************************************************************************
 *  @file vector_range_2d.c
 *  @brief 2D RANGE VECTOR objects
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _VECTOR_RANGE_2D_H
#define _VECTOR_RANGE_2D_H

/* import datatype dependencies */
#include "../structs.h"
#include "../edgebound.h"
#include "vector_int.h"
#include "vector_range.h"

/* VECTOR_2D struct */
typedef struct {
   VECTOR_INT     *id;          /* the ID for each range data */
   VECTOR_RANGE   *data;        /* array of data type */
   int            N;            /* current length of array in use */
   int            Nalloc;       /* current length of array allocated */
}  VECTOR_RANGE_2D;

/* constructor */
VECTOR_RANGE_2D* VECTOR_RANGE_2D_Create();

/* constructor with initial size and fill  */
VECTOR_RANGE_2D* VECTOR_RANGE_2D_Create_Size( const int init_size );

/* destructor */
void VECTOR_RANGE_2D_Destroy( VECTOR_RANGE_2D *vec );

/* deep copy */
VECTOR_RANGE_2D* VECTOR_RANGE_2D_Copy();

/* resize the array */
void VECTOR_RANGE_2D_Resize( VECTOR_RANGE_2D *vec, float growth_factor );

/* push element onto end of array */
void VECTOR_RANGE_2D_Pushback( VECTOR_RANGE_2D *vec, VECTOR_RANGE val );

/* pop element from end of array */
VECTOR_RANGE VECTOR_RANGE_2D_Pop( VECTOR_RANGE_2D *vec );

/* set data at index (no bound checks) */
void VECTOR_RANGE_2D_Set( VECTOR_RANGE_2D *vec, int idx, VECTOR_RANGE val );

/* get data at index (no bound checks) */
VECTOR_RANGE VECTOR_RANGE_2D_Get( VECTOR_RANGE_2D *vec, int idx );

/* equality test */
int VECTOR_RANGE_2D_Compare( VECTOR_RANGE_2D *vecA, VECTOR_RANGE_2D *vecB );

/* merge current diagonal bound into vectors */
void VECTOR_RANGE_2D_MergeFwd( VECTOR_RANGE_2D *vec, BOUND bnd );

/* convert to edgebound */
EDGEBOUNDS* VECTOR_RANGE_2D_Convert_to_Edgebound( VECTOR_RANGE_2D *vec );

/* unit test */
void VECTOR_RANGE_2D_UnitTest();

#endif /* _VECTOR_RANGE_2D_H */