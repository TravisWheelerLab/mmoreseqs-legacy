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
#include <time.h>

/* objects */
#include "../objects/structs.h"
#include "../objects/structs.h"
#include "../objects/edgebound.h"

/* header */
#include "misc.h"

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

/*! FUNCTION:  logsum()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (approximation).
 *             Speedup using LOGSUM_LOOKUP table means no exp() or log() operations are performed.
 *             LOGSUM_LOOKUP must have been initialized before use.
 *             Note: This one of the vals x,y = -inf, then other value overrides it.
 */
inline
float logsum( float     x,
              float     y )
{
   float max, min;

   if (x > y) {
      max = x;
      min = y;
   }
   else {
      max = y;
      min = x;
   }

   return (min == -INF || (max - min) >= 15.7f) ?
          max : max + LOGSUM_LOOKUP[ (int)((max - min) * LOGSUM_SCALE) ];
}

/*! FUNCTION:  logsum_post()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (approximation).
 *             Speedup using LOGSUM_LOOKUP table means no exp() or log() operations are performed.
 *             LOGSUM_LOOKUP must have been initialized before use.
 *             Note: If one of the vals x,y = -inf, then other value overrides it.
 *             Note: Handles
 */
inline
float logsum_post(   float     x,
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

/*! FUNCTION:  logsum_exact()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real sum (exact).
 *             Slower than logsum().
 */
inline
float logsum_exact( float  x,
                    float  y )
{
   return log( exp(x) + exp(y) );
}

/*! FUNCTION:  logprod()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real product.
 *             This correctly handles -INF and INF.
 */
inline
float logprod( float  x,
               float  y )
{
   // float max, min, ret;

   return x + y;

   // if (x > y) {
   //    max = x;
   //    min = y;
   // }
   // else {
   //    max = y;
   //    min = x;
   // }
   
   // // if ( max == INF ) {
   // //    return min;
   // // }
   // if ( max == -INF ) {
   //    return -INF;
   // } 
   // if ( min == -INF ) {
   //    return max;
   // }
   // else {
   //    return min + max;
   // }
}

/*! NOTE: This is a WIP.  Currently does not work as intended.
 *  FUNCTION:  logdiff()
 *  SYNOPSIS:  Takes two log-scaled numbers and returns the log-scale of their real diff (approximation).
 *             Speedup using LOGSUM_LOOKUP table means no exp() or log() operations are performed.
 *             LOGSUM_LOOKUP must have been initialized before use.
 */
inline
float logdiff( float  x,
               float  y )
{
   float max; 
   float min;

   y = -y;
   
   if ( x == -INF && y == -INF ) return -INF;
   if ( x == INF && y == INF ) return INF;

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

/*! FUNCTION:  logsum_Dump()
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

/*! FUNCTION:  negln2real()
 *  SYNOPSIS:  Convert negative natural logarithm probability to real probability.
 */
float negln2real( float negln_prob )
{
   float real_prob = 1.0f / ( exp( negln_prob ) );
   return real_prob;
}

/*! FUNCTION:  real2negln()
 *  SYNOPSIS:  Convert real probabilities to negative natural logarithm probability.
 */
float real2negln(float real_prob)
{
   float negln_prob = -1.0f * log(real_prob);
   return negln_prob;
}

/*! FUNCTION:  cmp_tol()
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

/*! FUNCTION:  cmp_tol()
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

/*! FUNCTION:  my_mkdir()
 *  SYNOPSIS:  Make directory.
 *    RETURN:  Returns <STATUS_SUCCESS> on success.
 */
int my_mkdir( const char* folderpath )
{
    return STATUS_SUCCESS;
}
  
/*! FUNCTION:  my_delay()
 *  SYNOPSIS:  Hold for given number of <secs>
 */
void my_delay( int milli_seconds ) 
{ 
    // Converting time into milli_seconds 
    int nano_seconds =  1000 * milli_seconds; 
  
    // Storing start time 
    clock_t start_time = clock(); 
  
    // looping till required time is not achieved 
    while (clock() < start_time + nano_seconds) {}; 
} 