/*******************************************************************************
 *  - FILE:      float.h
 *  - DESC:    FLT Object
 *******************************************************************************/

#ifndef _FLT_H
#define _FLT_H

/*! FUNCTION:  FLT_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
FLT FLT_Create(const FLT data);

/*! FUNCTION:  FLT_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
FLT FLT_Destroy(FLT data);

/*! FUNCTION:  FLT_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do
 * nothing.
 */
FLT FLT_Clear(FLT data);

/*! FUNCTION:  FLT_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* FLT_ToString(const FLT data, char* buf);

/*! FUNCTION:  FLT_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* FLT_ToExpString(const FLT data, char* buf);

/*! FUNCTION:  FLT_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
FLT FLT_FromString(char* str);

/*! FUNCTION:  FLT_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b),
 *             0 if equal,
 *             neg if (a < b)
 */
int FLT_Compare(const FLT a, const FLT b);

/*! FUNCTION:  FLT_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b),
 *             0 if equal,
 *             NEG if (a < b)
 */
int FLT_CompareTo(const void* a, const void* b);

#endif /* _FLT_H */
