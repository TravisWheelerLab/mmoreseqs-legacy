/*******************************************************************************
 *  FILE:      bound.h
 *  PURPOSE:   BOUND Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _BOUND_H
#define _BOUND_H

/* === DATA STRUCT === */

/* === FUNCTIONS === */

/*! FUNCTION:  BOUND_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
BOUND
BOUND_Create( const BOUND   data );

/*! FUNCTION:  BOUND_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
BOUND
BOUND_Destroy( BOUND   data );

/*! FUNCTION:  BOUND_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing. 
 */
BOUND
BOUND_Clear( BOUND   data );

/*! FUNCTION:  BOUND_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* 
BOUND_To_String(  const BOUND    d,
                  char*          buf );

/*! FUNCTION:  BOUND_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int 
BOUND_Compare(    const BOUND   a, 
                  const BOUND   b );

/*! FUNCTION:  BOUND_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
int 
BOUND_CompareTo(    const void*   a, 
                  const void*   b );

#endif /* _BOUND_H */
