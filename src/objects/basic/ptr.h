/*******************************************************************************
 *  - FILE:      ptr.h
 *  - DESC:    PTR Object.
 *             Void pointer wrapper.
 *******************************************************************************/

#ifndef _PTR_H
#define _PTR_H

/* include where PTR is defined */
#include "../structs.h"

/*! FUNCTION:  PTR_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
PTR PTR_Create(const PTR data);

/*! FUNCTION:  PTR_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
PTR PTR_Destroy(PTR data);

/*! FUNCTION:  PTR_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do
 * nothing.
 */
PTR PTR_Clear(PTR data);

/*! FUNCTION:  PTR_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a preallocated char* buffer <buf> of length
 * <buf_size>. RETURN:  Pointer to <buf>.
 */
char* PTR_ToString(const PTR data, char* buf);

/*! FUNCTION:  PTR_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
PTR PTR_FromString(char* str);

/*! FUNCTION:  PTR_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  pos if (a > b),
 *             0 if equal,
 *             neg if (a < b)
 */
int PTR_Compare(const PTR a, const PTR b);

/*! FUNCTION:  PTR_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b),
 *             0 if equal,
 *             NEG if (a < b)
 */
int PTR_CompareTo(const void* a, const void* b);

/*! FUNCTION:  PTR_Swap()
 *  SYNOPSIS:  Swap values of <a> and <b>
 */
void PTR_Swap(PTR* a, PTR* b);

#endif /* _PTR_H */
