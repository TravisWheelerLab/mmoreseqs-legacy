/*******************************************************************************
 *  FILE:      statistics.c
 *  PURPOSE:   Statistics for converting bit-score, P-value, and E-value
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
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

#define SMALLX1 0.001

/*
 *  FUNCTION:  	gumbel_pdf()
 *  SYNOPSIS:  	Return the right-tail mass about Gumbel Probability Density Function, G(x). 
 * 				P( G(X) > x)
 */
double gumbel_pdf(double x, double mu, double lambda)
{
  double y;
  y = lambda * (x - mu);
  return (lambda * exp(-y - exp(-y)));
}


/*
 *  FUNCTION:  	gumbel_log_pdf()
 *  SYNOPSIS:  	
 */
double gumbel_log_pdf(double x, double mu, double lambda)
{
  double y;
  y = lambda * (x - mu);
  return (log(lambda) -y - exp(-y));
}


/*
 *  FUNCTION:  	gumbel_cdf()
 *  SYNOPSIS:  	
 */
double gumbel_cdf(double x, double mu, double lambda)
{
  double y;
  y = lambda*(x-mu);
  return exp(-exp(-y));
}


/*
 *  FUNCTION:  	gumbel_log_cdf()
 *  SYNOPSIS:  	
 */
double gumbel_log_cdf(double x, double mu, double lambda)
{
  double y;
  y = lambda*(x-mu);
  return (-exp(-y));
}

/*
 *  FUNCTION:  	gumbel_surv()
 *  SYNOPSIS:  	
 */
double gumbel_surv(double x, double mu, double lambda)
{
	double y  = lambda*(x-mu);
	double ey = -exp(-y);

	/* Use 1-e^x ~ -x approximation here when e^-y is small. */
	if (fabs(ey) < SMALLX1) 
		return -ey;
	else                       
		return 1 - exp(ey);
}


/*
 *  FUNCTION:  	gumbel_log_surv()
 *  SYNOPSIS:  	
 */
double gumbel_log_surv(double x, double mu, double lambda)
{
	double y  = lambda*(x-mu);
	double ey = -exp(-y);

	/* The real calculation is log(1-exp(-exp(-y))).
	* For "large" y, -exp(-y) is small, so 1-exp(-exp(-y) ~ exp(-y),
	* and log of that gives us -y.
	* For "small y", exp(-exp(-y) is small, and we can use log(1-x) ~ -x. 
	*/
	if      (fabs(ey)      < SMALLX1) 
		return -y;
	else if (fabs(exp(ey)) < SMALLX1) 
		return -exp(ey);
	else
		return log(1-exp(ey));
}