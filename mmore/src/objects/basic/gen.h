/*******************************************************************************
 *  FILE:      gen.h
 *  PURPOSE:   GEN Object.  A union which can hold most primitive datatypes.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _GEN_H
#define _GEN_H

#include "../structs.h"

/*! FUNCTION:  GEN_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
GEN
GEN_Create( const GEN   data );

/*! FUNCTION:  GEN_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
GEN
GEN_Destroy( GEN   data );

/*! FUNCTION:  GEN_To_String()
 *  SYNOPSIS:  Create a string representation of <data>.
 *             If it is of float-like type, formats with <sig_digits> as number of significant digits.
 *             Stores it in a char* buffer <buf>. Caller must have preallocated buffer. 
 *             Strings are truncated to buf_size.
 *    RETURN:  Pointer to <buf>.
 */
char* 
GEN_To_String( 	const GEN      data,
                  char*          buf,
                  const int      buf_size,
                  const int      sig_digits );

/*! FUNCTION:  BOUND_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int 
GEN_Compare(   const GEN   a, 
               const GEN   b );

#endif /* _GEN_H */