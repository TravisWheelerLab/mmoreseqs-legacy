/******************************************************************************
 *  FILE:      trace.c
 *  PURPOSE:   TRACE Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _TRACE_H
#define _TRACE_H

/*! FUNCTION:  TRACE_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
TRACE
TRACE_Create( const TRACE   data );

/*! FUNCTION:  TRACE_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
TRACE
TRACE_Destroy( TRACE   data );

/*! FUNCTION:  TRACE_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing. 
 */
TRACE
TRACE_Clear( TRACE   data );

/*! FUNCTION:  TRACE_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* 
TRACE_To_String(  const TRACE   	d,
                  char*       	buf );

/*! FUNCTION:  TRACE_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int 
TRACE_Compare( const TRACE  a, 
               const TRACE  b );

/*! FUNCTION:  TRACE_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
int 
TRACE_CompareTo(    const void*   a, 
                  const void*   b );

#endif /* _TRACE_H */