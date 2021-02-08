/*******************************************************************************
 *  FILE:      float.h
 *  PURPOSE:   FLT Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _FLT_H
#define _FLT_H

/*! FUNCTION:  FLT_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
FLT
FLT_Create( const FLT   data );

/*! FUNCTION:  FLT_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
FLT
FLT_Destroy( FLT   data );

/*! FUNCTION:  FLT_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* 
FLT_To_String( 	const FLT   	d,
						char*       	buf );

/*! FUNCTION:  FLT_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int 
FLT_Compare(  	const FLT   a, 
					const FLT   b );

#endif /* _FLT_H */
