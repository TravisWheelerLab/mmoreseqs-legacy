/*******************************************************************************
 *  FILE:      logmath.h
 *  PURPOSE:   Miscellaneous Math Functions.
 *             Primary related to logrithmic math.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _LOGMATH_H
#define _LOGMATH_H

/* === imports === */
/* NONE */

/* === macros === */
#define LOGSUM_SCALE    1000.0f
#define LOGSUM_TBL      16000

/* === public functions === */

/*! FUNCTION:  calc_Max()
 *  SYNOPSIS:  Returns maximum value of two floats.
 */
float 
calc_Max(   float   x, 
            float   y );

/*! FUNCTION:  calc_Min()
 *  SYNOPSIS:  Returns mimumum value of two floats.
 */
float 
calc_Min(   float   x, 
            float   y );


/*! FUNCTION:  logsum_Init()
 *  SYNOPSIS:  Initializes values for the Logsum lookup table, LOGSUM_LOOKUP. 
 *             Computes the natural logarithm for real numbers, scaled down by constant LOGSUM_SCALE (to prevent underflow of real numbers).
 *             Must be initialized before calling logsum();
 */
void 
logsum_Init();

/*! FUNCTION:  logsum()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (approximation).
 *             Speedup using LOGSUM_LOOKUP table means no exp() or log() operations are performed.
 *             LOGSUM_LOOKUP must have been initialized before use.
 *             Note: This one of the vals x,y = -inf, then other value overrides it.
 */
float 
logsum(  float     x,
         float     y );

/*! FUNCTION:  logsum_post()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (approximation).
 *             Speedup using LOGSUM_LOOKUP table means no exp() or log() operations are performed.
 *             LOGSUM_LOOKUP must have been initialized before use.
 *             Note: If one of the vals x,y = -inf, then other value overrides it.
 *             Note: Handles
 */
float 
logsum_post(   float    x,
               float    y );

/*! FUNCTION:  logsum_exact()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (exact).
 *             Slower than logsum().
 */
float 
logsum_exact(  float  x, 
               float  y );

/*! FUNCTION:  logprod()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real product.
 */
float 
logprod(    float  x,
            float  y );

/*! FUNCTION:  logsum_Dump()
 *  SYNOPSIS:  Output LOGSUM_LOOKUP table to file.
 */
void 
logsum_Dump( FILE *fp );

/*! FUNCTION:  negln2real()
 *  SYNOPSIS:  Convert negative natural logarithm probability to real probability.
 */
float 
negln2real( float    negln_prob );

/*! FUNCTION:  real2negln()
 *  SYNOPSIS:  Convert real probabilities to negative natural logarithm probability.
 */
float 
real2negln( float    real_prob );

/*! FUNCTION:  cmp_tol()
 *  SYNOPSIS:  Checks whether two floats <a> and <b> are within static tolerance <tol> of eachother.  
 *    RETURN:  Returns true if numbers absolute difference is less than tolerance <tol>.
 */
bool 
cmp_tol(    const float    a, 
            const float    b);

/*! FUNCTION:  cmp_tol()
 *  SYNOPSIS:  Checks whether two floats <a> and <b> are within custom tolerance <tol> of eachother.  
 *    RETURN:  Returns true if numbers absolute difference is less than tolerance <tol>.
 */
bool 
cmp_tol_cust(  const float    a, 
               const float    b,
               const float    tol );

/*! FUNCTION:  my_delay()
 *  SYNOPSIS:  Hold for given number of <secs>
 */
void 
my_delay( int  secs );

#endif /* _LOGMATH_H */