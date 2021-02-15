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

/* === public functions === */

/*! FUNCTION:  MATH_Max()
 *  SYNOPSIS:  Returns maximum value of two floats.
 */
float 
MATH_Max(   const float   x, 
            const float   y );

/*! FUNCTION:  MATH_Min()
 *  SYNOPSIS:  Returns mimumum value of two floats.
 */
float 
MATH_Min(   const float   x, 
            const float   y );

/*! FUNCTION:  MATH_Identity()
 *  SYNOPSIS:  Identity function.  Simple passthrough function that returns the value of <x>.
 */
float 
MATH_Identity(   const float   x );

/*! FUNCTION:  MATH_Exp()
 *  SYNOPSIS:  Computes the exponential of <x>.
 */
float 
MATH_Exp(   const float   x );

/*! FUNCTION:  MATH_Log()
 *  SYNOPSIS:  Computes the natural logrithm of <x>.
 */
float 
MATH_Log(   const float   x );

/*! FUNCTION:  MATH_Zero()
 *  SYNOPSIS:  Method Selector for computing zero.
 */
float 
MATH_Zero();

/*! FUNCTION:  MATH_LogZero()
 *  SYNOPSIS:  Reteurn the log(0) = -INF.
 */
float 
MATH_LogZero();

/*! FUNCTION:  MATH_Zero()
 *  SYNOPSIS:  Method Selector for computing zero.
 */
float 
MATH_One();

/*! FUNCTION:  MATH_LogZero()
 *  SYNOPSIS:  Return the log(1) = 0.0f.
 */
float 
MATH_LogOne();

/*! FUNCTION:  MATH_Sum()
 *  SYNOPSIS:  Method Selector for computing sum of two numbers.
 */
float 
MATH_Sum(   const float     x,
            const float     y );

/*! FUNCTION:  MATH_Logsum_Init()
 *  SYNOPSIS:  Initializes values for the Logsum lookup table, LOGSUM_LOOKUP. 
 *             Computes the natural logarithm for real numbers, scaled down by constant LOGSUM_SCALE (to prevent underflow of real numbers).
 *             Must be initialized before calling MATH_Sum();
 */
void 
MATH_Logsum_Init();

/*! FUNCTION:  MATH_Logsum_exact()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (exact).
 *             Slower than MATH_Sum().
 */
float 
MATH_LogSum(   const float  x, 
               const float  y );

/*! FUNCTION:  MATH_Logsum_exact()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (exact).
 *             Slower than MATH_Sum().
 */
float 
MATH_Logsum_exact(   const float  x, 
                     const float  y );

/*! FUNCTION:  MATH_Prod()
 *  SYNOPSIS:  Method selector for computing product of two numbers.
 */
float 
MATH_Prod(  const float  x,
            const float  y );

/*! FUNCTION:  MATH_LogProd()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real product.
 */
float 
MATH_LogProd(  const float  x,
               const float  y );

/*! FUNCTION:  MATH_Logsum_Dump()
 *  SYNOPSIS:  Output LOGSUM_LOOKUP table to file.
 */
void 
MATH_Logsum_Dump( FILE *fp );

/*! FUNCTION:  MATH_NegLn2Real()
 *  SYNOPSIS:  Convert negative natural logarithm probability to real probability.
 */
float 
MATH_NegLn2Real( const float    negln_prob );

/*! FUNCTION:  MATH_Real2NegLn()
 *  SYNOPSIS:  Convert real probabilities to negative natural logarithm probability.
 */
float 
MATH_Real2NegLn( const float    real_prob );

/*! FUNCTION:  MATH_CmpTol()
 *  SYNOPSIS:  Checks whether two floats <a> and <b> are within static tolerance <tol> of eachother.  
 *    RETURN:  Returns true if numbers absolute difference is less than tolerance <tol>.
 */
bool 
MATH_CmpTol(   const float    a, 
               const float    b );

/*! FUNCTION:  MATH_CmpTol()
 *  SYNOPSIS:  Checks whether two floats <a> and <b> are within custom tolerance <tol> of eachother.  
 *    RETURN:  Returns true if numbers absolute difference is less than tolerance <tol>.
 */
bool 
MATH_CmpTol_byTol(   const float    a, 
                     const float    b,
                     const float    tol );

#endif /* _LOGMATH_H */