/*******************************************************************************
 *  FILE:      vector_int.c
 *  PURPOSE:   VECTOR_CHAR Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_CHAR_H
#define _VECTOR_CHAR_H

/*
 *  FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Create new VECTOR_CHAR object and returns pointer.
 */
VECTOR_CHAR* VECTOR_CHAR_Create();

/*
 *  FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Create new VECTOR_CHAR object at specific size and returns pointer.
 */
VECTOR_CHAR* VECTOR_CHAR_Create_by_Size( int size );

/* destructor */
void* VECTOR_CHAR_Destroy( VECTOR_CHAR* 	vec );

/* reuse by resetting counter*/
void VECTOR_CHAR_Reuse( VECTOR_CHAR* 		vec );

/* set all active indexes to zero */
void VECTOR_CHAR_Clear( VECTOR_CHAR* 		vec );

/* deep copy */
VECTOR_CHAR* VECTOR_CHAR_Copy();

/* resize the array */
void VECTOR_CHAR_Resize( 	VECTOR_CHAR* 	vec, 
						 			int  		   	size );

/* push element onto end of array */
void VECTOR_CHAR_Pushback( 	VECTOR_CHAR* 	vec, 
										char 				val );

/* pop element from end of array */
char VECTOR_CHAR_Pop( VECTOR_CHAR* 	vec );

/* set data at index (no bound checks) */
void VECTOR_CHAR_Set(	VECTOR_CHAR* 	vec, 
								int 				idx, 
								char 				val );

/* get data at index (no bound checks) */
char* VECTOR_CHAR_Get( 	VECTOR_CHAR* 	vec, 
								int 				idx );

/* equality test */
int VECTOR_CHAR_Compare( 	VECTOR_CHAR* 	vecA, 
									VECTOR_CHAR* 	vecB );

/* output VECTOR_CHAR to file */
void VECTOR_CHAR_Dump(   	VECTOR_CHAR* 	vec,
                        	FILE*       	fp );

#endif 