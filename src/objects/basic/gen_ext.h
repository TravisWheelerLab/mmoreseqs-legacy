/*******************************************************************************
 *  - FILE:      gen_ext.h
 *  - DESC:    GEN Object.  A union which can hold most primitive datatypes.
 *             Extended functionality.
 *******************************************************************************/

#ifndef _GEN_EXT_H
#define _GEN_EXT_H

#include "../structs.h"

/*! FUNCTION:  GEN_Wrap()
 *  SYNOPSIS:  Create a GEN struct.
 *             <data> should be reference to the data to be stored.
 *             <type> should be one of the enumerated datatypes.
 *             <size> should be sizeof(<type>).
 *    RETURN:  Pointer to GEN struct.
 */
GEN GEN_Wrap(const void* data, const DATATYPE type, const size_t size);

/*! FUNCTION:  GEN_GetSize()
 *  SYNOPSIS:  Get the size of an enumerated datatype.
 *    RETURN:  Pointer to <buf>, or NULL if error.
 */
size_t GEN_GetSize(const DATATYPE type);

#endif /* _GEN_EXT_H */
