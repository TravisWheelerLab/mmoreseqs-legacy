/*******************************************************************************
 *  - FILE:      int.h
 *  - DESC:    INT Object
 *******************************************************************************/

#ifndef _INT_H
#define _INT_H

/*! FUNCTION:  INT_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
INT INT_Create(const INT data);

/*! FUNCTION:  INT_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
INT INT_Destroy(INT data);

/*! FUNCTION:  INT_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do
 * nothing.
 */
INT INT_Clear(INT data);

/*! FUNCTION:  INT_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *    RETURN:  Pointer to <buf>
 */
char* INT_ToString(const INT data, char* buf);

/*! FUNCTION:  INT_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
INT INT_FromString(char* str);

/*! FUNCTION:  INT_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  pos if (a > b),
 *             0 if equal,
 *             neg if (a < b)
 */
int INT_Compare(const INT a, const INT b);

/*! FUNCTION:  INT_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b),
 *             0 if equal,
 *             NEG if (a < b)
 */
int INT_CompareTo(const void* a, const void* b);

/*! FUNCTION:  INT_Swap()
 *  SYNOPSIS:  Swap values of <a> and <b>
 */
void INT_Swap(INT* a, INT* b);

#endif /* _INT_H */
