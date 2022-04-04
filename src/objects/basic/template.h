/*******************************************************************************
 *  - FILE: template.h
 *  - DESC:  XXX Object
 *******************************************************************************/

#ifndef _XXX_H
#define _XXX_H

/* include where XXX is defined */
#include "../structs.h"

/*! FUNCTION:  XXX_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
XXX XXX_Create(const XXX data);

/*! FUNCTION:  XXX_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
XXX XXX_Destroy(XXX data);

/*! FUNCTION:  XXX_Clear()
 *  SYNOPSIS:  Initialize or Clear <data>.  If pointer data, sets to null.
 * Otherwise, do nothing.
 */
XXX XXX_Clear(XXX data);

/*! FUNCTION:  XXX_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a preallocated char* buffer <buf> of length
 * <buf_size>. RETURN:  Pointer to <buf>.
 */
char* XXX_ToString(const XXX data, char* buf);

/*! FUNCTION:  XXX_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
XXX XXX_FromString(char* str);

/*! FUNCTION:  XXX_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  POSITIVE if (a > b),
 *             ZERO if equal,
 *             NEGATIVE if (a < b)
 */
int XXX_Compare(const XXX a, const XXX b);

/*! FUNCTION:  XXX_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POSITIVE if (a > b),
 *             ZERO if equal,
 *             NEGATIVE if (a < b)
 */
int XXX_CompareTo(const void* a, const void* b);

/*! FUNCTION:  XXX_Equals()
 *  SYNOPSIS:  Checks if <a> and <b> are equal.
 *    RETURN:  TRUE if equal, FALSE otherwise
 */
bool XXX_Equals(const XXX a, const XXX b);

/*! FUNCTION:  XXX_Swap()
 *  SYNOPSIS:  Swap values of <a> and <b>
 */
void XXX_Swap(XXX* a, XXX* b);

#endif /* _XXX_H */
