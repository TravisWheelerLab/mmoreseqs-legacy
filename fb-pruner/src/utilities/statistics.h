/*******************************************************************************
 *  FILE:      statistics.c
 *  PURPOSE:   Statistics for converting bit-score, P-value, and E-value.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

#ifndef _STATISTICS_H
#define _STATISTICS_H

/*
 *  FUNCTION:     GUMBEL_pdf()
 *  SYNOPSIS:     Return the right-tail mass about Gumbel Probability Density Function, G(mu, lambda). 
 *             <mu> and <lambda> are parameters of the Gumbel distribution.
 *             y = lambda * (x - mu)
 *             Pr( G(mu,lambda) > x ) = lamda * exp(-(y) - exp(-(y))).
 */
double GUMBEL_pdf(   double x, 
                     double mu, 
                     double lambda );


/*
 *  FUNCTION:     GUMBEL_log_pdf()
 *  SYNOPSIS:     
 */
double GUMBEL_log_pdf(  double x, 
                        double mu, 
                        double lambda );


/*
 *  FUNCTION:     GUMBEL_cdf()
 *  SYNOPSIS:     
 */
double GUMBEL_cdf(   double x, 
                     double mu, 
                     double lambda );


/*
 *  FUNCTION:     GUMBEL_log_cdf()
 *  SYNOPSIS:     
 */
double GUMBEL_log_cdf(  double x, 
                        double mu, 
                        double lambda );

/*
 *  FUNCTION:     GUMBEL_survivor()
 *  SYNOPSIS:     
 */
double GUMBEL_survivor(    double x, 
                           double mu, 
                           double lambda );


/*
 *  FUNCTION:     GUMBEL_log_survivor()
 *  SYNOPSIS:     
 */
double GUMBEL_log_survivor(   double x, 
                              double mu, 
                              double lambda);

#endif /* _STATISTICS_H */