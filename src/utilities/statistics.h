/*******************************************************************************
 *  - FILE:      statistics.h
 *  - DESC:    Statistical functions for converting between scores in NATS or
 *BITS, P-value, and E-value.
 *******************************************************************************/

#ifndef _STATISTICS_H
#define _STATISTICS_H

#include "../objects/structs.h"

/*! FUNCTION:   STATS_Nats_to_Bitscore()
 *  SYNOPSIS:   Given a score in Nats <nat_sc> and
 */
float STATS_Nats_to_Bits(float nats_sc);

/*! FUNCTION:  	STATS_Bits_to_Nats()
 *  SYNOPSIS:  	Converts a score in BITS to score in NATS.
 */
float STATS_Bits_to_Nats(float bit_sc);

/*! FUNCTION:  	STATS_Pval_to_Eval()
 *  SYNOPSIS:  	Converts a P-VALUE <pval> to and E-VALUE <eval>.
 * 					Requires <db_size>
 */
float STATS_Pval_to_Eval(float pval, int db_size);

/*! FUNCTION:  	STATS_Eval_to_Pval()
 *  SYNOPSIS:  	Converts a E-VALUE <pval> to and P-VALUE <eval>.
 * 					Requires <db_size> (number of sequences in database).
 */
float STATS_Eval_to_Pval(float eval, int db_size);

/*! FUNCTION:  	STATS_Viterbi_Natsc_to_Eval()
 *  SYNOPSIS:  	Get an E-VALUE <eval> from a Viterbi NAT score <natsc>.
 * 					Requires <gumbel_params> = <mu,tau> parameters which fits Gumbel
 * distribution for the model. Requires <db_size>, which is the number of
 * sequence queries in the database. Optionally, can include the <null1_bias>
 * model and <null2_bias> sequence biases, or left as 0.0f. Optionally, other
 * scores <presc>, <seqsc>, <pval> can be captured or left NULL. RETURN:
 * E-value.
 */
float STATS_Viterbi_Nats_to_Eval(
    float natsc, /* INPUT: Viterbi output score (in NATS) */
    float*
        presc_p,              /* OPT OUTPUT: non-sequence bias corrected score (in BITS) */
    float* seqsc_p,           /* OPT OUTPUT: sequence bias corrected score (in BITS) */
    float* pval_p,            /* OPT OUTPUT: P-value (probably to a match given a random
                       sequence) */
    float* eval_p,            /* OPT OUTPUT: E-value (expected number of matches in given
                       <db_size> of random sequences) */
    DIST_PARAM gumbel_params, /* parameters for fitting gumbel model */
    int db_size,              /* number of query sequences in database */
    float null1_hmm_bias,     /* null1 model bias (also called null_sc) */
    float null2_seq_bias);    /* null2 sequence bias (also called seq_bias) */

/*! FUNCTION:  	STATS_Fwdback_Natsc_to_Eval()
 *  SYNOPSIS:  	Get an E-VALUE <eval> from a Forward-Backward NAT score <natsc>.
 * 					Requires <gumbel_params> = <mu,tau> parameters which fits Gumbel
 * distribution for the model. Requires <db_size>, which is the number of
 * sequence queries in the database. Optionally, can include the <null1_bias>
 * model and <null2_bias> sequence biases, or left as 0.0f. Optionally, other
 * scores <presc>, <seqsc>, <pval> can be captured or left NULL. RETURN:
 * E-value.
 */
float STATS_Fwdback_Nats_to_Eval(
    float natsc, /* INPUT: Viterbi output score (in NATS) */
    float*
        presc_p,           /* OPT OUTPUT: non-sequence bias corrected score (in BITS) */
    float* seqsc_p,        /* OPT OUTPUT: sequence bias corrected score (in BITS) */
    float* pval_p,         /* OPT OUTPUT: P-value (probably to a match given a random
                       sequence) */
    float* eval_p,         /* OPT OUTPUT: E-value (expected number of matches in given
                       <db_size> of random sequences) */
    DIST_PARAM exp_params, /* parameters for fitting exponential model */
    int db_size,           /* number of query sequences in database */
    float null1_hmm_bias,  /* null1 model bias (also called null_sc) */
    float null2_seq_bias); /* null2 sequence bias (also called seq_bias) */

#endif /* _STATISTICS_H */
