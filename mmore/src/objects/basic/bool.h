/*******************************************************************************
 *  FILE:      bool.h
 *  PURPOSE:   BOOL Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _BOOL_H
#define _BOOL_H

/*! FUNCTION:  BOOL_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
BOOL
BOOL_Create( const BOOL   data );

/*! FUNCTION:  BOOL_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
BOOL
BOOL_Destroy( BOOL   data );

/*! FUNCTION:  BOOL_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *    RETURN:  Pointer to <buf>
 */
char* 
BOOL_To_String( const BOOL   data,
                char*        buf );

/*! FUNCTION:  BOOL_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int 
BOOL_Compare(  const BOOL   a, 
               const BOOL   b );

#endif /* _BOOL_H */
