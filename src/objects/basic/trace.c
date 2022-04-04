/*******************************************************************************
 *  - FILE:      trace.c
 *  - DESC:    TRACE Object
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
#include "trace.h"

/*! FUNCTION:  TRACE_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
inline TRACE
TRACE_Create(const TRACE data) {
  return data;
}

/*! FUNCTION:  TRACE_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
inline TRACE
TRACE_Destroy(TRACE data) {
  return data;
}

/*! FUNCTION:  TRACE_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing.
 */
inline TRACE
TRACE_Clear(TRACE data) {
  return data;
}

/*! FUNCTION:  TRACE_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline char*
TRACE_ToString(const TRACE data,
               char* buf) {
  sprintf(buf, "{ [%d] (%d,%d) }", data.st, data.q_0, data.t_0);
  return buf;
}

/*! FUNCTION:  TRACE_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
inline TRACE
TRACE_FromString(char* str) {
  TRACE data;
  return data;
}

/*! FUNCTION:  TRACE_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  POS if (a > b),
 *             0 if equal,
 *             NEG if (a < b)
 */
inline int
TRACE_Compare(const TRACE a,
              const TRACE b) {
  if ((a.q_0 - b.q_0) != 0) {
    return (a.q_0 - b.q_0);
  }

  if ((a.t_0 - b.t_0) != 0) {
    return (a.t_0 - b.t_0);
  }

  return (a.st - b.st);
}

/*! FUNCTION:  TRACE_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b),
 *             0 if equal,
 *             NEG if (a < b)
 */
inline int
TRACE_CompareTo(const void* a,
                const void* b) {
  TRACE* x = (TRACE*)a;
  TRACE* y = (TRACE*)b;

  return TRACE_Compare(*x, *y);
}
