/*******************************************************************************
 *  - FILE:   vector_float_ext.c
 *  - DESC:    VECTOR_FLT Object (Extended) Functions.
 *             Float type-specific vector operations.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>

/* local imports */
#include "../../objects/structs.h"
#include "../../utilities/_utilities.h"
#include "../../objects/_objects.h"

/* header */
#include "_vectors.h"
#include "vector_float.h"

/*! FUNCTION:  VECTOR_FLT_Logf()
 *  SYNOPSIS:  Perform logrithmic function on <vec>
 */
STATUS_FLAG
VECTOR_FLT_Logf(VECTOR_FLT* vec) {
  int L = VECTOR_FLT_GetSize(vec);
  for (int i = 0; i < L; i++) {
    if (VEC_X(vec, i) > 0.0) {
      VEC_X(vec, i) = logf(VEC_X(vec, i));
    } else {
      VEC_X(vec, i) = -INF;
    }
  }
  return STATUS_SUCCESS;
}
