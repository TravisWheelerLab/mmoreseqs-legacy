/*******************************************************************************
 *  @file VECTOR_FLT.c
 *  @brief CHARACTER VECTOR objects
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _VECTOR_FLT_H
#define _VECTOR_FLT_H

/* constructor */
VECTOR_FLT* VECTOR_FLT_Create();

/* destructor */
void VECTOR_FLT_Destroy( VECTOR_FLT* 	vec );

/* deep copy */
VECTOR_FLT* VECTOR_FLT_Copy();

/* resize the array */
void VECTOR_FLT_Resize( VECTOR_FLT* 	vec, 
						const int 		size );

/* push element onto end of array */
void VECTOR_FLT_Pushback( VECTOR_FLT* 	vec, 
						  const float 	val );

/* pop element from end of array */
float VECTOR_FLT_Pop( VECTOR_FLT* 	vec );

/* set data at index (no bound checks) */
void VECTOR_FLT_Set( VECTOR_FLT* 	vec, 
					 const int 		idx, 
					 const float 	val );

/* get data at index (no bound checks) */
float VECTOR_FLT_Get( 	VECTOR_FLT* 	vec, 
						const int 		idx );

/* equality test */
int VECTOR_FLT_Compare( const VECTOR_FLT* 	vec_A, 
						const VECTOR_FLT* 	vec_B );

#endif /* _VECTOR_FLT_H */