/*******************************************************************************
 *  FILE:      vector_int.c
 *  PURPOSE:   VECTOR_INT Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_INT_H
#define _VECTOR_INT_H

/*
 *  FUNCTION:  VECTOR_INT_Create()
 *  SYNOPSIS:  Create new VECTOR_INT object and returns pointer.
 */
VECTOR_INT* VECTOR_INT_Create();

/*
 *  FUNCTION:  VECTOR_INT_Create()
 *  SYNOPSIS:  Create new VECTOR_INT object at specific size and returns pointer.
 */
VECTOR_INT* VECTOR_INT_Create_by_Size( int size );

/* destructor */
void* VECTOR_INT_Destroy( VECTOR_INT* 	vec );

/* reuse by resetting counter*/
void VECTOR_INT_Reuse( VECTOR_INT* 		vec );

/* set all active indexes to zero */
void VECTOR_INT_Clear( VECTOR_INT* 		vec );

/* deep copy */
VECTOR_INT* VECTOR_INT_Copy();

/* resize the array */
void VECTOR_INT_Resize( VECTOR_INT* 	vec, 
								int  		   	size );

/* push element onto end of array */
void VECTOR_INT_Pushback( 	VECTOR_INT* 	vec, 
									int 				val );

/* pop element from end of array */
int VECTOR_INT_Pop( VECTOR_INT* 	vec );

/* set data at index (no bound checks) */
void VECTOR_INT_Set(	VECTOR_INT* 	vec, 
							int 				idx, 
							int 				val );

/* get data at index (no bound checks) */
int* VECTOR_INT_Get( VECTOR_INT* 	vec, 
							int 				idx );

/* equality test */
int VECTOR_INT_Compare( VECTOR_INT* 	vecA, 
								VECTOR_INT* 	vecB );

/* output VECTOR_INT to file */
void VECTOR_INT_Dump(   VECTOR_INT* 	vec,
                        FILE*       	fp );

#endif 