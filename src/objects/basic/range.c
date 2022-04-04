/*******************************************************************************
 *  - FILE:      trace.c
 *  - DESC:    RANGE Object
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
#include "_basic.h"
#include "range.h"

/*! FUNCTION:  RANGE_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
inline RANGE
RANGE_Create(const RANGE data) {
  return data;
}

/*! FUNCTION:  RANGE_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
inline RANGE
RANGE_Destroy(RANGE data) {
  return data;
}

/*! FUNCTION:  RANGE_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing.
 */
inline RANGE
RANGE_Clear(RANGE data) {
  return data;
}

/*! FUNCTION:  RANGE_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline char*
RANGE_ToString(const RANGE data,
               char* buf) {
  sprintf(buf, "(%d,%d)", data.beg, data.end);
  return buf;
}

/*! FUNCTION:  RANGE_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
inline RANGE
RANGE_FromString(char* str) {
  RANGE data;
  return data;
}

/*! FUNCTION:  RANGE_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b),
 *             0 if equal,
 *             neg if (a < b)
 */
inline int
RANGE_Compare(const RANGE a,
              const RANGE b) {
  if ((a.beg - b.beg) != 0) {
    return (a.beg - b.beg);
  }

  return (a.end - b.end);
}

/*! FUNCTION:  RANGE_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b),
 *             0 if equal,
 *             NEG if (a < b)
 */
inline int
RANGE_CompareTo(const void* a,
                const void* b) {
  RANGE* x = (RANGE*)a;
  RANGE* y = (RANGE*)b;

  return RANGE_Compare(*x, *y);
}
