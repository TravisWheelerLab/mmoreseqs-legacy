/******************************************************************************
 *  @file misc.h
 *  @brief Miscellaneous helper functions.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

#ifndef _MISC_H
#define _MISC_H

/* constants */
#define LOGSUM_SCALE 1000.f
#define LOGSUM_TBL 16000

/* CHECK FOR NULL POINTER */
int alloc_pointer( void**  ptr, 
                   int     size );

/* MINIMUM & MAXIMUM FUNCTIONS */
float calc_Max ( float  x, 
                 float  y );
float calc_Min ( float  x, 
                 float  y );

/* LOGRITHMIC SUM FUNCTIONS */
void init_Logsum ();

/* Takes two logscale numbers and returns the log of their real sum (approx) */
float calc_Logsum( float   x, 
                   float   y );

/* Takes two log numbers and returns the log of their real sum (exact) */
float calc_Logsum_exact( float   x, 
                         float   y);
/* Outputs the LOGSUM_TBL to FILE pointer */
void print_Logsum( FILE*   fp );

/* CONVERT REAL NUMBERS to NEGATIVE LOG NUMBERS */
/* convert negative natural logs to real probabilities */
float negln2real( float  negln_prob );

/* convert real probabilities to negative natural logs */
float real2negln( float  real_prob );

/* check if two numbers are within a tolerance of eachother */
bool cmp_tol(const float a, 
			 const float b);

/* frees memory at pointer and sets to null */
void myfree(void* ptr);

#endif /* _MISC_H */