/*******************************************************************************
 *  - FILE:      bool.h
 *  - DESC:    BOOL Object
 *******************************************************************************/

#ifndef _BOOL_H
#define _BOOL_H

/*! FUNCTION:  BOOL_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
BOOL BOOL_Create(const BOOL data);

/*! FUNCTION:  BOOL_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
BOOL BOOL_Destroy(BOOL data);

/*! FUNCTION:  BOOL_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do
 * nothing.
 */
BOOL BOOL_Clear(BOOL data);

/*! FUNCTION:  BOOL_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *    RETURN:  Pointer to <buf>
 */
char* BOOL_ToString(const BOOL data, char* buf);

/*! FUNCTION:  BOOL_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
BOOL BOOL_FromString(char* str);

/*! FUNCTION:  BOOL_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  pos if (a > b),
 *             0 if equal,
 *             neg if (a < b)
 */
int BOOL_Compare(const BOOL a, const BOOL b);

/*! FUNCTION:  BOOL_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b),
 *             0 if equal,
 *             NEG if (a < b)
 */
int BOOL_CompareTo(const void* a, const void* b);

#endif /* _BOOL_H */
