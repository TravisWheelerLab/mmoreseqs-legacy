/*******************************************************************************
 *  FILE:      gen_data.h
 *  PURPOSE:   GEN_DATA Object.  A union which can hold most primitive datatypes.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _GEN_DATA_H
#define _GEN_DATA_H

#include "../structs.h"

/*! FUNCTION:  GEN_DATA_Create()
 *  SYNOPSIS:  Create a GEN_DATA struct. 
 *             <data> should be reference to the data to be stored.
 *             <type> should be one of the enumerated datatypes.
 *             <size> should be sizeof(<type>). 
 *    RETURN:  Pointer to GEN_DATA struct.
 */
GEN_DATA 
GEN_DATA_Create( 	const void*       data,
                  const DATATYPE    type,
                  const size_t      size ); 

/*! FUNCTION:  GEN_DATA_To_String()
 *  SYNOPSIS:  Create a string representation of <data>.
 *             If it is of float-like type, formats with <sig_digits> as number of significant digits.
 *             Stores it in a char* buffer <buf>. Caller must have preallocated buffer. 
 *             Strings are truncated to buf_size.
 *    RETURN:  Pointer to <buf>.
 */
char* 
GEN_DATA_To_String( 	const GEN_DATA    data,
                     char*             buf,
                     const int         buf_size,
                     const int         sig_digits );

/*! FUNCTION:  BOUND_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int 
GEN_DATA_Compare(    const GEN_DATA   a, 
                     const GEN_DATA   b );

#endif /* _GEN_DATA_H */
