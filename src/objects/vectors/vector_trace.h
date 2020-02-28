/*******************************************************************************
 *  @file VECTOR_TRACE.c
 *  @brief TRACE VECTOR objects
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _VECTOR_TRACE_H
#define _VECTOR_TRACE_H

/* import datatypes */
#include "../structs.h"

// typedef struct {
//    int         i;            /* index in query */
//    int         j;            /* index in target */
//    int         st;           /* state at index */
// } TRACE;

// /* VECTOR struct */
// typedef struct {
//    TRACE       *data;        /* array of data type */
//    int         N;            /* current length of array in use */
//    int         Nalloc;       /* current length of array allocated */
// }  VECTOR_TRACE;


/* constructor */
VECTOR_TRACE* VECTOR_TRACE_Create();

/* destructor */
void VECTOR_TRACE_Destroy( VECTOR_TRACE *vec );

/* deep copy */
VECTOR_TRACE* VECTOR_TRACE_Copy();

/* resize the array */
void VECTOR_TRACE_Resize( VECTOR_TRACE *vec, float growth_factor );

/* push element onto end of array */
void VECTOR_TRACE_Pushback( VECTOR_TRACE *vec, TRACE val );

/* pop element from end of array */
TRACE VECTOR_TRACE_Pop( VECTOR_TRACE *vec );

/* set data at index (no bound checks) */
void VECTOR_TRACE_Set( VECTOR_TRACE *vec, int idx, TRACE val );

/* get data at index (no bound checks) */
TRACE VECTOR_TRACE_Get( VECTOR_TRACE *vec, int idx );

/* equality test */
int VECTOR_TRACE_Compare( VECTOR_TRACE *vecA, VECTOR_TRACE *vecB );

#endif /* _VECTOR_TRACE_H */