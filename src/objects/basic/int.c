/*******************************************************************************
 *  - FILE:      int.c
 *  - DESC:    INT Object
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "../structs.h"
#include "../../utilities/_utilities.h"
#include "../_objects.h"

/* header */
#include "int.h"

/*! FUNCTION:  INT_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
inline INT
INT_Create(const INT data) {
  return data;
}

/*! FUNCTION:  INT_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
inline INT
INT_Destroy(INT data) {
  return data;
}

/*! FUNCTION:  INT_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing.
 */
inline INT
INT_Clear(INT data) {
  return data;
}

/*! FUNCTION:  INT_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline char*
INT_ToString(const INT data,
             char* buf) {
  sprintf(buf, "%d", data);
  return buf;
}

/*! FUNCTION:  INT_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
inline INT
INT_FromString(char* str) {
  INT data;
  data = atoi(str);
  return data;
}

/*! FUNCTION:  INT_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  POSITIVE if (a > b),
 *             ZERO if equal,
 *             NEGATIVE if (a < b)
 */
inline int
INT_Compare(const INT a,
            const INT b) {
  return (a - b);
}

/*! FUNCTION:  INT_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POSITIVE if (a > b),
 *             ZERO if equal,
 *             NEGATIVE if (a < b)
 */
inline int
INT_CompareTo(const void* a,
              const void* b) {
  INT* x = (INT*)a;
  INT* y = (INT*)b;

  return INT_Compare(*x, *y);
}

/*! FUNCTION:  INT_Swap()
 *  SYNOPSIS:  Swap values of <a> and <b>.
 */
inline void
INT_Swap(INT* a,
         INT* b) {
  INT tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}
