/*******************************************************************************
 *  @file VECTOR_BOUND.c
 *  @brief VECTOR BOUND Vectors 
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _VECTOR_BOUND_H
#define _VECTOR_BOUND_H

/* constructor */
VECTOR_BOUND* VECTOR_BOUND_Create();

/* destructor */
void VECTOR_BOUND_Destroy( VECTOR_BOUND*  vec );

/* deep copy */
VECTOR_BOUND* VECTOR_BOUND_Copy();

/* resize the array */
void VECTOR_BOUND_Resize( VECTOR_BOUND*   vec, 
                          int             size);

/* resize the array if smaller than new size */
void VECTOR_BOUND_GrowTo( VECTOR_BOUND*   vec, 
                          int             size );

/* push element onto end of array */
void VECTOR_BOUND_Pushback( VECTOR_BOUND*  vec, 
                            BOUND          val );

/* pop element from end of array */
BOUND VECTOR_BOUND_Pop( VECTOR_BOUND*  vec );

/* set data at index (no bound checks) */
void VECTOR_BOUND_Set( VECTOR_BOUND*  vec, 
                       int            idx, 
                       BOUND          val );

/* get data at index (no bound checks) */
BOUND VECTOR_BOUND_Get( VECTOR_BOUND*  vec, 
                        int            idx );

/* equality test */
int VECTOR_BOUND_Compare( VECTOR_BOUND*  vecA, 
                          VECTOR_BOUND*  vecB );

/* sort vector ascending */
int VECTOR_BOUND_sort( VECTOR_BOUND*   vec );


#endif /* _VECTOR_BOUND_H */