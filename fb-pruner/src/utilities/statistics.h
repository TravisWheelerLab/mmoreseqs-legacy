/*******************************************************************************
 *  FILE:      statistics.c
 *  PURPOSE:   Statistics for converting bit-score, P-value, and E-value.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

#ifndef _STATISTICS_H
#define _STATISTICS_H

#include "../objects/structs.h"

/* 
 * FUNCTION:  STATS_Nats_to_Bitscore()
 * SYNOPSIS:  Given a score in Nats <nat_sc> and 
 */
int STATS_Nats_to_Bitscore(   float    nats_sc,
                              float    filter_sc,
                              float*   bit_sc );

/* 
 * FUNCTION:   BG_FilterScore()
 * SYNOPSIS:   Calculates the filter null model score.
 *
 * PURPOSE:    Modeled after p7_bg_FilterScore() in HMMER.
 *             Calculates the filter null model <bg> score for sequence
 *             <dsq> of length <L>, and return it in 
 *             <*ret_sc>.
 *            
 *             The score is calculated as an HMM Forward score using
 *             the two-state filter null model. It is a log-odds ratio,
 *             relative to the iid background frequencies, in nats:
 *             same as main model Forward scores.
 *
 *             The filter null model has no length distribution of its
 *             own; the same geometric length distribution (controlled
 *             by <bg->p1>) that the null1 model uses is imposed.
 */
int BG_FilterScore(  HMM_PROFILE*   t_prof, 
                     SEQUENCE*      q_seq, 
                     MATRIX_3D*     st_MX,
                     MATRIX_2D*     sp_MX,
                     float*         filter_sc );

/* 
 * FUNCTION:  p7_bg_FilterScore()
 * SYNOPSIS:  Calculates the filter null model score.
 */
int NULL_MODEL_Forward(    HMM_PROFILE*   t_prof, 
                           SEQUENCE*      q_seq, 
                           MATRIX_3D*     st_MX,
                           MATRIX_2D*     sp_MX,
                           float*         null_sc );

/* 
 * FUNCTION:  EXPONENTIAL_survivor()
 * SYNOPSIS:  Calculates the survivor function, $P(X>x)$ (that is, 1-CDF,
 *            the right tail probability mass) for an exponential distribution,
 *            given value <x>, offset <mu>, and decay parameter <lambda>.
 */
double EXPONENTIAL_survivor(  double x, 
                              double mu, 
                              double lambda );


/*
 *  FUNCTION:     GUMBEL_pdf()
 *  SYNOPSIS:     Return the right-tail mass about Gumbel Probability Density Function, G(mu, lambda). 
 *                <mu> and <lambda> are parameters of the Gumbel distribution.
 *                y = lambda * (x - mu)
 *                Pr( G(mu,lambda) > x ) = lamda * exp(-(y) - exp(-(y))).
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