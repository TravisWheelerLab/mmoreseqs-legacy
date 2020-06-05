/*******************************************************************************
 *  @file misc.h
 *  @brief Miscellaneous Helper Functions.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* objects */
#include "objects/structs.h"
#include "objects/edgebound.h"

/* local imports */
#include "structs.h"
#include "objects.h"

/* header */
#include "utilities.h"

/* GLOBAL STATIC VARIABLES */
/* table of logsum values */
static float LOGSUM_LOOKUP[LOGSUM_TBL];
bool LOGSUM_INITIALIZED = false;

/*
 *  FUNCTION:  calc_Max()
 *  SYNOPSIS:  Returns maximum value of two floats.
 */
inline
float calc_Max( float   x,
                float   y )
{
   if (x > y)
   {
      return x;
   }
   return y;
}

/*
 *  FUNCTION:  calc_Min()
 *  SYNOPSIS:  Returns mimumum value of two floats.
 */
inline
float calc_Min( float   x,
                float   y )
{
   if (x < y)
   {
      return x;
   }
   return y;
}

/*
 *  FUNCTION:  logsum_Init()
 *  SYNOPSIS:  Initializes values for the Logsum lookup table, LOGSUM_LOOKUP.
 *             Computes the natural logarithm for real numbers, scaled down by constant LOGSUM_SCALE (to prevent underflow of real numbers).
 *             Must be initialized before calling logsum();
 */
void logsum_Init()
{
   if (LOGSUM_INITIALIZED == false)
   {
      LOGSUM_INITIALIZED = true;
      int i;
      for (i = 0; i < LOGSUM_TBL; i++)
      {
         LOGSUM_LOOKUP[i] = log(1. + exp((double) - i / LOGSUM_SCALE));
      }
   }
}

/*
 *  FUNCTION:  logsum()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (approximation).
 *             Speedup using LOGSUM_LOOKUP table means no exp() or log() operations are performed.
 *             LOGSUM_LOOKUP must have been initialized before use.
 */
inline
float logsum( float     x,
              float     y )
{
   float max, min;

   if (x > y)
   {
      max = x;
      min = y;
   }
   else
   {
      max = y;
      min = x;
   }

   return (min == -INF || (max - min) >= 15.7f) ?
          max : max + LOGSUM_LOOKUP[ (int)((max - min) * LOGSUM_SCALE) ];
}

/*
 *  FUNCTION:  logsum_exact()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (exact).
 *             Slower than logsum().
 */
inline
float logsum_exact( float  x,
                    float  y )
{
   return log( exp(x) + exp(y) );
}

/*
 *  FUNCTION:  logsum_Dump()
 *  SYNOPSIS:  Output LOGSUM_LOOKUP table to file.
 */
void logsum_Dump( FILE *fp )
{
   fprintf(fp, "=== LOGSUM TABLE ===\n");
   for (int i = 0; i < 16000; i += 160)
   {
      fprintf(fp, "[%d] %.2f\n", i, LOGSUM_LOOKUP[i]);
   }
   fprintf(fp, "\n\n");
}

/*
 *  FUNCTION:  negln2real()
 *  SYNOPSIS:  Convert negative natural logarithm probability to real probability.
 */
float negln2real( float negln_prob )
{
   float real_prob = 1.0f / ( exp( negln_prob ) );
   return real_prob;
}

/*
 *  FUNCTION:  real2negln()
 *  SYNOPSIS:  Convert real probabilities to negative natural logarithm probability.
 */
float real2negln(float real_prob)
{
   float negln_prob = -1.0f * log(real_prob);
   return negln_prob;
}

/*
 *  FUNCTION:  cmp_tol()
 *  SYNOPSIS:  Checks whether two floats <a> and <b> are within static tolerance <tol> of eachother.
 *    RETURN:  Returns true if numbers absolute difference is less than tolerance <tol>.
 */
inline
bool cmp_tol( const float a,
              const float b)
{
   /* acceptable tolerance range for "equality tests" */
   const float tol = 1e-5;
   return ( fabs( a - b ) < tol );
}

/*
 *  FUNCTION:  cmp_tol()
 *  SYNOPSIS:  Checks whether two floats <a> and <b> are within custom tolerance <tol> of eachother.
 *    RETURN:  Returns true if numbers absolute difference is less than tolerance <tol>.
 */
inline
bool cmp_tol_cust(const float a,
                  const float b,
                  const float tol )
{
   /* acceptable tolerance range for "equality tests" */
   return ( fabs( a - b ) < tol );
}

/*
 *  FUNCTION:  cmp_tol()
 *  SYNOPSIS:  Checks whether two floats <a> and <b> are within custom tolerance <tol> of eachother.
 *    RETURN:  Returns <1> on success, <0> on failure
 */
int my_mkdir( const char* folderpath )
{
  
}