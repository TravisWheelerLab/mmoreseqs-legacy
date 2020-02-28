/*******************************************************************************
 *  @file VECTOR_CHAR.c
 *  @brief CHARACTER VECTOR objects
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _VECTOR_CHAR_H
#define _VECTOR_CHAR_H

/* import datatypes */
#include "../structs.h"

/* Vector struct */
typedef struct {
   char *data;      /* array of data type */
   int N;            /* current length of array in use */
   int Nalloc;       /* current length of array allocated */
}  VECTOR_CHAR;


/* constructor */
VECTOR_CHAR* VECTOR_CHAR_Create();

/* destructor */
void VECTOR_CHAR_Destroy( VECTOR_CHAR *vec );

/* deep copy */
VECTOR_CHAR* VECTOR_CHAR_Copy();

/* resize the array */
void VECTOR_CHAR_Resize( VECTOR_CHAR *vec, const float growth_factor );

/* push element onto end of array */
void VECTOR_CHAR_Pushback( VECTOR_CHAR *vec, const char val );

/* pop element from end of array */
char VECTOR_CHAR_Pop( VECTOR_CHAR *vec );

/* set data at index (no bound checks) */
void VECTOR_CHAR_Set( VECTOR_CHAR *vec, const int idx, const char val );

/* get data at index (no bound checks) */
char VECTOR_CHAR_Get( VECTOR_CHAR *vec, const int idx );

/* equality test */
int VECTOR_CHAR_Compare( const VECTOR_CHAR *vecA, const VECTOR_CHAR *vecB );

#endif /* _VECTOR_CHAR_H */