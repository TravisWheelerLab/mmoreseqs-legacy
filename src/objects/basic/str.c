/*******************************************************************************
 *  - FILE:  str.c
 *  - DESC:   STR Object.
 *             Mostly a simple wrapper for char*.
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

/* header */
#include "../_objects.h"
#include "_basic.h"
#include "str.h"

/*! FUNCTION:  STR_Create()
 *  SYNOPSIS:  Create <str>, fills with input <str_chrs> and return pointer to <str>.
 */
inline STR
STR_Create(const char* chrs) {
  /* if string points to null, return null */
  if (chrs == NULL)
    return NULL;

  STR str;
  str = strdup(chrs);
  return str;
}

/*! FUNCTION:  STR_Destroy()
 *  SYNOPSIS:  Destroys <str>, frees memory and returns NULL pointer.
 */
inline STR
STR_Destroy(STR str) {
  if (str == NULL)
    return str;

  str = ERROR_free(str);
  return str;
}

/*! FUNCTION:  STR_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing.
 */
inline STR
STR_Clear(STR data) {
  return NULL;
}

/*! FUNCTION:  STR_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *    RETURN:  Pointer to <buf>
 */
inline char*
STR_ToString(const STR data,
             char* buf) {
  if (data == NULL) {
    sprintf(buf, "%s", "(null)");
  } else {
    sprintf(buf, "%s", data);
  }

  return buf;
}

/*! FUNCTION:  STR_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
inline STR
STR_FromString(char* str) {
  STR data;
  data = STR_Create(str);
  return data;
}

/*! FUNCTION:  STR_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  pos if (a > b),
 *             0 if equal,
 *             neg if (a < b)
 */
inline int
STR_Compare(const STR a,
            const STR b) {
  return strcmp(a, b);
}

/*! FUNCTION:  STR_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b),
 *             0 if equal,
 *             NEG if (a < b)
 */
inline int
STR_CompareTo(const void* a,
              const void* b) {
  STR* x = (STR*)a;
  STR* y = (STR*)b;

  return STR_Compare(*x, *y);
}
