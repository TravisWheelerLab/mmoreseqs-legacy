/*******************************************************************************
 *  - FILE:      logmath.h
 *  - DESC:    Miscellaneous Math Functions.
 *             Primary related to logrithmic math.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* objects */
#include "../objects/structs.h"
#include "../objects/_objects.h"

/* header */
#include "_utilities.h"
#include "mymath.h"

/* macros */
#define LOGSUM_SCALE 1000.0f
#define LOGSUM_TBL 16000

/* table of logsum values */
static float LOGSUM_LOOKUP[LOGSUM_TBL];
bool LOGSUM_INITIALIZED = false;

/*! FUNCTION:  MATH_Max()
 *  SYNOPSIS:  Returns maximum value of two floats <x> and <y>.
 */
inline float
MATH_Max(const float x,
         const float y) {
  if (x > y) {
    return x;
  }
  return y;
}

/*! FUNCTION:  MATH_Min()
 *  SYNOPSIS:  Returns mimumum value of two floats <x> and <y>.
 */
inline float
MATH_Min(const float x,
         const float y) {
  if (x < y) {
    return x;
  }
  return y;
}

/*! FUNCTION:  MATH_Identity()
 *  SYNOPSIS:  Identity function.  Simple passthrough function that returns the value of <x>.
 */
inline float
MATH_Identity(const float x) {
  return x;
}

/*! FUNCTION:  MATH_Log()
 *  SYNOPSIS:  Computes the natural logrithm of <x>.
 */
inline float
MATH_Log(const float x) {
  return logf(x);
}

/*! FUNCTION:  MATH_Exp()
 *  SYNOPSIS:  Computes the exponential of <x>.
 */
inline float
MATH_Exp(const float x) {
  return expf(x);
}

/*! FUNCTION:  MATH_One()
 *  SYNOPSIS:  Method Selector for computing zero.
 */
float MATH_One() {
  return MATH_LogOne();
}

/*! FUNCTION:  MATH_NormalOne()
 *  SYNOPSIS:  Return normal space one = 1.
 */
float MATH_NormalOne() {
  return 1.0f;
}

/*! FUNCTION:  MATH_LogOne()
 *  SYNOPSIS:  Return the log(0) = -INF.
 */
inline float
MATH_LogOne() {
  return 0.0f;
}

/*! FUNCTION:  MATH_Zero()
 *  SYNOPSIS:  Method Selector for computing zero.
 */
float MATH_Zero() {
  return MATH_LogZero();
}

/*! FUNCTION:  MATH_NormalZero()
 *  SYNOPSIS:  Return normal space zero = 0.
 */
float MATH_NormalZero() {
  return 0.0f;
}

/*! FUNCTION:  MATH_LogZero()
 *  SYNOPSIS:  Reteurn the log(0) = -INF.
 */
inline float
MATH_LogZero() {
  return -INF;
}

/*! FUNCTION:  MATH_Sum()
 *  SYNOPSIS:  Method Selector.
 */
inline float
MATH_Sum(const float x,
         const float y) {
  return MATH_LogSum(x, y);
}

/*! FUNCTION:  MATH_NormalSum()
 *  SYNOPSIS:  Method Selector.
 */
inline float
MATH_NormalSum(const float x,
               const float y) {
  return x + y;
}

/*! FUNCTION:  MATH_Logsum_Init()
 *  SYNOPSIS:  Initializes values for the Logsum lookup table, LOGSUM_LOOKUP.
 *             Computes the natural logarithm for real numbers, scaled down by constant LOGSUM_SCALE (to prevent underflow of real numbers).
 *             Must be initialized before calling MATH_Sum();
 */
void MATH_Logsum_Init() {
  if (LOGSUM_INITIALIZED == false) {
    LOGSUM_INITIALIZED = true;
    int i;
    for (i = 0; i < LOGSUM_TBL; i++) {
      LOGSUM_LOOKUP[i] = log(1. + exp((double)-i / LOGSUM_SCALE));
    }
  }
}

/*! FUNCTION:  MATH_LogSum()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (approximation).
 *             Speedup using LOGSUM_LOOKUP table means no exp() or log() operations are performed.
 *             LOGSUM_LOOKUP must have been initialized before use.
 *             NOTE: This one of the vals x,y = -inf, then other value overrides it.
 */
inline float
MATH_LogSum(const float x,
            const float y) {
  float max, min;

  if (x > y) {
    max = x;
    min = y;
  } else {
    max = y;
    min = x;
  }

  return (min == -INF || (max - min) >= 15.7f) ? max : max + LOGSUM_LOOKUP[(int)((max - min) * LOGSUM_SCALE)];
}

/*! FUNCTION:  MATH_Logsum_explicit()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (exact).
 *             Slower than MATH_Logsum().
 */
inline float
MATH_Logsum_explicit(const float x,
                     const float y) {
  float sum;
  sum = exp(x) + exp(y);

  if (sum < 0) {
    fprintf(stderr, "ERROR: Log of negative value attempted: exp(%.3f) + exp(%.3e) = %.3e + %.3e = %.3f.\n",
            x, y, exp(x), exp(y), sum);
    ERRORCHECK_exit(EXIT_FAILURE);
    return -INF;
  }
  elif (sum == 0) {
    return -INF;
  }

  return log(sum);
}

/*! FUNCTION:  MATH_Logdiff_explicit()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real difference (exact).
 *             Slower than MATH_Logsum().
 */
inline float
MATH_Logdiff_explicit(const float x,
                      const float y) {
  float sum;
  sum = exp(x) - exp(y);

  if (sum < 0) {
    fprintf(stderr, "ERROR: Attempted to compute log of negative value: exp(%.3f) + exp(%.3e) = %.3e + %.3e = %.3f.\n",
            x, y, exp(x), exp(y), sum);
    // ERRORCHECK_exit(EXIT_FAILURE);
    return -INF;
  }
  elif (sum == 0) {
    return -INF;
  }

  return log(sum);
}

/*! FUNCTION:  MATH_Prod()
 *  SYNOPSIS:  Takes the product of two numbers.
 */
inline float
MATH_Prod(const float x,
          const float y) {
  return MATH_LogProd(x, y);
}

/*! FUNCTION:  MATH_NormalProd()
 *  SYNOPSIS:  Takes the product of two numbers.
 */
inline float
MATH_NormalProd(const float x,
                const float y) {
  return x * y;
}

/*! FUNCTION:  MATH_LogProd()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real product.
 */
inline float
MATH_LogProd(const float x,
             const float y) {
  return x + y;
}

/*! FUNCTION:  MATH_Logsum_Dump()
 *  SYNOPSIS:  Output LOGSUM_LOOKUP table to file.
 */
void MATH_Logsum_Dump(FILE* fp) {
  fprintf(fp, "=== LOGSUM TABLE ===\n");
  for (int i = 0; i < 16000; i += 160) {
    fprintf(fp, "[%d] %.2f\n", i, LOGSUM_LOOKUP[i]);
  }
  fprintf(fp, "\n\n");
}

/*! FUNCTION:  MATH_NegLn2Real()
 *  SYNOPSIS:  Convert negative natural logarithm probability to real probability.
 */
float MATH_NegLn2Real(const float negln_prob) {
  float real_prob = 1.0f / (exp(negln_prob));
  return real_prob;
}

/*! FUNCTION:  MATH_Real2NegLn()
 *  SYNOPSIS:  Convert real probabilities to negative natural logarithm probability.
 */
float MATH_Real2NegLn(const float real_prob) {
  float negln_prob = -1.0f * log(real_prob);
  return negln_prob;
}

/*! FUNCTION:  MATH_CmpTol()
 *  SYNOPSIS:  Checks whether two floats <a> and <b> are within static tolerance <tol> of eachother.
 *    RETURN:  Returns true if numbers absolute difference is less than tolerance <tol>.
 */
inline bool
MATH_CmpTol(const float a,
            const float b) {
  /* acceptable tolerance range for "equality tests" */
  const float tol = 1e-5;
  bool is_diff_within_tol = (fabs(a - b) < tol);
  return is_diff_within_tol;
}

/*! FUNCTION:  MATH_CmpTol_byTol()
 *  SYNOPSIS:  Checks whether two floats <a> and <b> are within custom tolerance <tol> of eachother.
 *    RETURN:  Returns true if numbers absolute difference is less than tolerance <tol>.
 */
inline bool
MATH_CmpTol_byTol(const float a,
                  const float b,
                  const float tol) {
  /* acceptable tolerance range for "equality tests" */
  bool is_diff_within_tol = (fabs(a - b) < tol);
  return is_diff_within_tol;
}
