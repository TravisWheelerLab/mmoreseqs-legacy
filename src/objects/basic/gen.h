/*******************************************************************************
 *  - FILE:      gen.h
 *  - DESC:    GEN Object.  A union which can hold most primitive datatypes.
 *******************************************************************************/

#ifndef _GEN_H
#define _GEN_H

#include "../structs.h"

/*! FUNCTION:  GEN_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
GEN GEN_Create(const GEN data);

/*! FUNCTION:  GEN_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
GEN GEN_Destroy(GEN data);

/*! FUNCTION:  GEN_ToString()
 *  SYNOPSIS:  Create a string representation of <data>.
 *             If it is of float-like type, formats with <sig_digits> as number
 * of significant digits. Stores it in a char* buffer <buf>. Caller must have
 * preallocated buffer. Strings are truncated to buf_size. RETURN:  Pointer to
 * <buf>.
 */
char* GEN_ToString(const GEN data,
                   char* buf,
                   const int buf_size,
                   const int sig_digits);

/*! FUNCTION:  GEN_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *    RETURN:  Pointer to <buf>
 */
GEN GEN_FromString(char* str);

/*! FUNCTION:  BOUND_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  POSITIVE if (a > b),
 *             ZERO if equal,
 *             NEGATIVE if (a < b)
 */
int GEN_Compare(const GEN a, const GEN b);

/*! FUNCTION:  GEN_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POSITIVE if (a > b),
 *             ZERO if equal,
 *             NEGATIVE if (a < b)
 */
int GEN_CompareTo(const void* a, const void* b);

#endif /* _GEN_H */
