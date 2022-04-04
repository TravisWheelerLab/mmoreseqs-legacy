/*******************************************************************************
 *  - FILE:      char.c
 *  - DESC:    CHAR Object
 *******************************************************************************/

#ifndef _CHAR_H
#define _CHAR_H

/*! FUNCTION:  CHAR_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
CHAR CHAR_Create(const CHAR data);

/*! FUNCTION:  CHAR_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
CHAR CHAR_Destroy(CHAR data);

/*! FUNCTION:  CHAR_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do
 * nothing.
 */
CHAR CHAR_Clear(CHAR data);

/*! FUNCTION:  CHAR_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* CHAR_ToString(const CHAR data, char* buf);

/*! FUNCTION:  CHAR_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
CHAR CHAR_FromString(char* str);

/*! FUNCTION:  CHAR_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b),
 *             0 if equal,
 *             neg if (a < b)
 */
int CHAR_Compare(const CHAR a, const CHAR b);

/*! FUNCTION:  CHAR_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b),
 *             0 if equal,
 *             NEG if (a < b)
 */
int CHAR_CompareTo(const void* a, const void* b);

#endif /* _CHAR_H */
