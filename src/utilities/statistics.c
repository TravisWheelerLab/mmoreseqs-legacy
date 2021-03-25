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

/* library imports */
#include "easel.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"

/* local imports */
#include "../objects/structs.h"
#include "../objects/_objects.h"

/* header */
#include "_utilities.h"
#include "statistics.h"

/*! FUNCTION:  	STATS_Nats_to_Bits()
 *  SYNOPSIS:  	Converts a score in NATS to score in BITS.
 */
float 
STATS_Nats_to_Bits( 	float 	nat_sc )
{
	float bit_sc;
	bit_sc = nat_sc / (float)CONST_LOG2;
	return bit_sc;
}

/*! FUNCTION:  	STATS_Bits_to_Nats()
 *  SYNOPSIS:  	Converts a score in BITS to score in NATS.
 */
float 
STATS_Bits_to_Nats( 	float 	bit_sc )
{
	float nat_sc;
	nat_sc = bit_sc * (float)CONST_LOG2;
	return nat_sc;
}

/*! FUNCTION:  	STATS_Pval_to_Eval()
 *  SYNOPSIS:  	Converts a P-VALUE <pval> to and E-VALUE <eval>.
 * 					Requires <db_size>
 */
float 
STATS_Pval_to_Eval( 	float 	pval,
							int 		db_size )
{
	float eval;
	eval = pval * (float)db_size;
	return eval;
}

/*! FUNCTION:  	STATS_Eval_to_Pval()
 *  SYNOPSIS:  	Converts a E-VALUE <pval> to and P-VALUE <eval>.
 * 					Requires <db_size> (number of sequences in database).
 */
float 
STATS_Eval_to_Pval( 	float 	eval,
							int 		db_size )
{
	float pval;
	pval = eval / (float)db_size;
	return pval;
}

/*! FUNCTION:  	STATS_Viterbi_Nats_to_Eval()
 *  SYNOPSIS:  	Get an E-VALUE <eval> from a Viterbi NAT score <natsc>.
 * 					Requires <gumbel_params> = <mu,tau> parameters which fits Gumbel distribution for the model.
 * 					Requires <db_size>, which is the number of sequence queries in the database.
 * 					Optionally, can include the <null1_bias> model and <null2_bias> sequence biases, or left as 0.0f.
 * 					Optionally, other scores <presc>, <seqsc>, <pval> can be captured or left NULL.
 *  RETURN: 		E-value.
 */
inline
float 
STATS_Viterbi_Nats_to_Eval( 	float 			natsc, 				/* INPUT: Viterbi output score (in NATS) */
										float* 			presc_p, 			/* OPT OUTPUT: non-sequence bias corrected score (in BITS) */
										float* 			seqsc_p, 			/* OPT OUTPUT: sequence bias corrected score (in BITS) */
										float* 			pval_p, 				/* OPT OUTPUT: P-value (probably to a match given a random sequence) */
										float* 			eval_p, 				/* OPT OUTPUT: E-value (expected number of matches in given <db_size> of random sequences) */
										DIST_PARAM		gumbel_params, 	/* parameters for fitting gumbel model */
										int 				db_size, 			/* number of query sequences in database */  	
										float				null1_hmm_bias, 	/* null1 model bias (also called null_sc) (in NATS) */
										float				null2_seq_bias ) 	/* null2 sequence bias (also called seq_bias) (in NATS) */
{
	float presc, seqsc, ln_pval, pval, eval;
	float mu, lambda;
	mu 		= gumbel_params.param1;
	lambda 	= gumbel_params.param2;

	presc 	= STATS_Nats_to_Bits( natsc - null1_hmm_bias );
	seqsc 	= STATS_Nats_to_Bits( natsc - (null1_hmm_bias + null2_seq_bias) );
	ln_pval 	= esl_gumbel_logsurv( seqsc, mu, lambda );
	pval 		= exp(ln_pval);
	eval 		= STATS_Pval_to_Eval( pval, db_size );

	/* output scores we want to keep */
	if (presc_p != NULL) *presc_p = presc;
	if (seqsc_p != NULL) *seqsc_p = seqsc;
	if (pval_p != NULL) 	*pval_p 	= pval;
	if (eval_p != NULL) 	*eval_p 	= eval;
	
	return eval;
}

/*! FUNCTION:  	STATS_Fwdback_Nats_to_Eval()
 *  SYNOPSIS:  	Get an E-VALUE <eval> from a Forward-Backward NAT score <natsc>.
 * 					Requires <gumbel_params> = <mu,tau> parameters which fits Gumbel distribution for the model.
 * 					Requires <db_size>, which is the number of sequence queries in the database.
 * 					Optionally, can include the <null1_bias> model and <null2_bias> sequence biases, or left as 0.0f.
 * 					Optionally, other scores <presc>, <seqsc>, <pval> can be captured or left NULL.
 *  RETURN: 		E-value.
 */
inline
float 
STATS_Fwdback_Nats_to_Eval( 	float 			natsc, 				/* INPUT: Viterbi output score (in NATS) */
										float* 			presc_p, 			/* OPT OUTPUT: non-sequence bias corrected score (in BITS) */
										float* 			seqsc_p, 			/* OPT OUTPUT: sequence bias corrected score (in BITS) */
										float* 			pval_p, 				/* OPT OUTPUT: P-value (probably to a match given a random sequence) */
										float* 			eval_p, 				/* OPT OUTPUT: E-value (expected number of matches in given <db_size> of random sequences) */
										DIST_PARAM		exp_params, 		/* parameters for fitting exponential model */
										int 				db_size, 			/* number of query sequences in database */  	
										float				null1_hmm_bias, 	/* null1 model bias (also called null_sc) (in NATS) */
										float				null2_seq_bias ) 	/* null2 sequence bias (also called seq_bias) (in NATS) */
{
	float presc, seqsc, ln_pval, pval, eval;
	float mu, lambda;
	mu 		= exp_params.param1;
	lambda 	= exp_params.param2;

	presc 	= STATS_Nats_to_Bits( natsc - null1_hmm_bias );
	seqsc 	= STATS_Nats_to_Bits( natsc - (null1_hmm_bias + null2_seq_bias) );
	ln_pval 	= esl_exp_logsurv( seqsc, mu, lambda );
	pval 		= exp(ln_pval);
	eval 		= STATS_Pval_to_Eval( pval, db_size );

	/* output scores to keep */
	if (presc_p != NULL) *presc_p = presc;
	if (seqsc_p != NULL) *seqsc_p = seqsc;
	if (pval_p != NULL) *pval_p = pval;
	if (eval_p != NULL) *eval_p = eval;
	
	return eval;
}

/* TODO: WIP */
/*! FUNCTION:  	STATS_Fwdback_Eval_to_Nats()
 *  SYNOPSIS:  	Get an NAT score <natsc> from a Forward-Backward NAT score <natsc>.
 * 					Requires <gumbel_params> = <mu,tau> parameters which fits Gumbel distribution for the model.
 * 					Requires <db_size>, which is the number of sequence queries in the database.
 * 					Optionally, can include the <null1_bias> model and <null2_bias> sequence biases, or left as 0.0f.
 * 					Optionally, other scores <presc>, <seqsc>, <pval> can be captured or left NULL.
 *  RETURN: 		Natscore.
 */
inline
float 
STATS_Fwdback_Eval_to_Nats( 	float 			eval, 				/* INPUT: E-value (expected number of matches in given <db_size> of random sequences) */
										float* 			pval_p, 				/* OPT OUTPUT: P-value (probably to a match given a random sequence) */
										float* 			seqsc_p, 			/* OPT OUTPUT: sequence bias corrected score (in BITS) */
										float* 			presc_p, 			/* OPT OUTPUT: non-sequence bias corrected score (in BITS) */
										float* 			natsc_p, 			/* OPT OUTPUT: Score (in NATS)  */
										DIST_PARAM		exp_params, 		/* parameters for fitting exponential model */
										int 				db_size, 			/* number of query sequences in database */  	
										float				null1_hmm_bias, 	/* null1 model bias (also called null_sc) */
										float				null2_seq_bias ) 	/* null2 sequence bias (also called seq_bias) */
{
	float natsc, presc, seqsc, ln_pval, pval;
	float mu, lambda;
	mu 		= exp_params.param1;
	lambda 	= exp_params.param2;

	pval 		= STATS_Eval_to_Pval( eval, db_size );
	ln_pval 	= log(ln_pval);
	/* TODO: Need the inverse to this exponential function to get score (user must supply accurate bias scores) */
	// seqsc 	= esl_exp_logsurv( ln_pval, mu, lambda );
	presc 	= STATS_Nats_to_Bits( natsc - null1_hmm_bias );
	seqsc 	= STATS_Nats_to_Bits( natsc - (null1_hmm_bias + null2_seq_bias) );

	/* output scores to keep (non-NULL pointers) */
	if (pval_p != NULL) 	*pval_p 	= pval;
	if (seqsc_p != NULL) *seqsc_p = seqsc;
	if (presc_p != NULL) *presc_p = presc;
	if (natsc_p != NULL) *natsc_p	= natsc;
	
	return natsc;
}
