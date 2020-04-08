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
#include "testing.h"
#include "hmm_parser.h"

/* header */
#include "utilities/utility.h"

/* GLOBAL STATIC VARIABLES */
/* table of logsum values */
static float LOGSUM_LOOKUP[LOGSUM_TBL];
bool LOGSUM_INITIALIZED = false;

/* allocates memory and checks for size */
int alloc_pointer( void**  ptr,
                   int     size )
{
   *ptr = realloc( *ptr, size );

   return ( *ptr == NULL );
}

/* Get max value of two floats */
inline
float calc_Max (float x, float y)
{
   if (x > y)
   {
      return x;
   }
   return y;
}

/* Get min value of two floats */
inline
float calc_Min (float x, float y)
{
   if (x < y)
   {
      return x;
   }
   return y;
}

/* initialize the logsum table */
void init_Logsum ()
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

/* Takes two logscale numbers and returns the log of their real sum (approx) */
// static inline
float calc_Logsum (float x, float y)
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

/* Takes two log numbers and returns the log of their real sum (exact) */
float calc_Logsum_exact(float x, float y)
{
   return log( exp(x) + exp(y) );
}

/* Outputs the LOGSUM_TBL to FILE pointer */
void print_Logsum(FILE *fp)
{
   fprintf(fp, "=== LOGSUM TABLE ===\n");
   for (int i = 0; i < 16000; i += 160)
   {
      fprintf(fp, "[%d] %.2f\n", i, LOGSUM_LOOKUP[i]);
   }
   fprintf(fp, "\n\n");
}

/* check for NULL pointer */
void check_null_ptr(void *ptr, char *str)
{
   if (ptr == NULL)
   {
      perror("ERROR: issuing while malloc'ing %s.\n");
      exit(EXIT_SUCCESS);
   }
}

/* convert negative natural logs to real probabilities */
float negln2real(float negln_prob)
{
   float real_prob = 1.0f / ( exp( negln_prob ) );
   return real_prob;
}

/* convert real probabilities to negative natural logs */
float real2negln(float real_prob)
{
   float negln_prob = -1.0f * log(real_prob);
   return negln_prob;
}
