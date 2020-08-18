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
#include "statistics.h"

#define SMALLX1 0.001

/*
 *  FUNCTION:  	GUMBEL_pdf()
 *  SYNOPSIS:  	Return the right-tail mass about Gumbel Probability Density Function, G(mu, lambda). 
 * 				<mu> and <lambda> are parameters of the Gumbel distribution.
 * 				y = lambda * (x - mu)
 * 				Pr( G(mu,lambda) > x ) = lamda * exp(-(y) - exp(-(y))).
 */
double GUMBEL_pdf(	double x, 
							double mu, 
							double lambda )
{
  double y;
  y = lambda * (x - mu);
  return (lambda * exp(-y - exp(-y)));
}


/*
 *  FUNCTION:  	GUMBEL_log_pdf()
 *  SYNOPSIS:  	
 */
double GUMBEL_log_pdf( 	double x, 
								double mu, 
								double lambda )
{
  double y;
  y = lambda * (x - mu);
  return (log(lambda) -y - exp(-y));
}


/*
 *  FUNCTION:  	GUMBEL_cdf()
 *  SYNOPSIS:  	
 */
double GUMBEL_cdf(	double x, 
							double mu, 
							double lambda )
{
  double y;
  y = lambda * (x - mu);
  return (exp(-exp(-y)));
}


/*
 *  FUNCTION:  	GUMBEL_log_cdf()
 *  SYNOPSIS:  	
 */
double GUMBEL_log_cdf( 	double x, 
								double mu, 
								double lambda )
{
  double y;
  y = lambda * (x - mu);
  return (-exp(-y));
}

/*
 *  FUNCTION:  	GUMBEL_survivor()
 *  SYNOPSIS:  	
 */
double GUMBEL_survivor( 	double x, 
									double mu, 
									double lambda )
{
	double y  = lambda*(x-mu);
	double ey = -exp(-y);

	/* Use 1-e^x ~ -x approximation here when e^-y is small. */
	if (fabs(ey) < SMALLX1) {
		return -ey;
	}
	else {
		return (1 - exp(ey));
	}
}


/*
 *  FUNCTION:  	GUMBEL_log_survivor()
 *  SYNOPSIS:  	
 */
double GUMBEL_log_survivor( 	double x, 
										double mu, 
										double lambda)
{
	double y  = lambda*(x-mu);
	double ey = -exp(-y);

	/* The real calculation is log(1-exp(-exp(-y))).
	* For "large" y, -exp(-y) is small, so 1-exp(-exp(-y) ~ exp(-y),
	* and log of that gives us -y.
	* For "small y", exp(-exp(-y) is small, and we can use log(1-x) ~ -x. 
	*/
	if ( fabs(ey) < SMALLX1 ) {
		return -y;
	} 
	else if ( fabs(exp(ey)) < SMALLX1 ) {
		return -exp(ey);
	} 
	else {
		return log(1 - exp(ey));
	}
}