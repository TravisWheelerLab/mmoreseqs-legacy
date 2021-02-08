/*******************************************************************************
 *  FILE:      hmm_bg.c
 *  PURPOSE:   HMM_BG Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUGS:
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

/* Easel imports */
#include "easel.h"
#include "esl_hmm.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"
#include "esl_sq.h"

/* local imports */
#include "structs.h"
#include "../utilities/_utilities.h"
#include "_objects.h"

/* header */
#include "hmm_bg.h"

/* NOTE: Modeled after P7_BG: p7_bg.c */

/*!   FUNCTION:  HMM_BG_Create()
 *    SYNOPSIS:  Create <hmm_bg>, allocate memory and return pointer.
 */
HMM_BG* 
HMM_BG_Create()
{
	ESL_ALPHABET* abc = esl_alphabet_Create(eslAMINO);
	if ( abc == NULL ) {
		printf("ERROR: malloc failed.\n");
		exit(EXIT_FAILURE);
	}

	HMM_BG* bg 	= (HMM_BG*) ERROR_malloc( sizeof(HMM_BG) );
	bg->f 		= NULL;
	bg->fhmm		= NULL;
	bg->abc 		= NULL;
	bg->sq 		= NULL;

	bg->f 		= (float*) ERROR_malloc( sizeof(float) * abc->K );
	bg->fhmm 	= esl_hmm_Create(abc, 2);
	bg->abc  	= abc;

	/* fill hardcoded amino frequencies */
	for ( int i = 0; i < abc->K; i++ ) {
		bg->f[i] = BG_MODEL[i];
	}

	bg->p1    = 350.0 / 351.0;
	bg->omega = 1.0 / 256.0;

	return bg;
}

/*!  FUNCTION:  HMM_BG_Destroy()
 *   SYNOPSIS:  Destroy <hmm_bg>, free memory, and return NULL pointer.
 */
HMM_BG* 
HMM_BG_Destroy( HMM_BG* bg )
{
	if (bg != NULL)
	{
		if (bg->fhmm != NULL) esl_alphabet_Destroy( bg->abc );
		bg->f = ERROR_free( bg->f );
		if (bg->fhmm != NULL) esl_hmm_Destroy( bg->fhmm );
		if (bg->sq != NULL) {
			bg->sq->dsq = ERROR_free( bg->sq->dsq ); 	/* sq->dsq is not handled by sq_Destroy() */
			esl_sq_Destroy( bg->sq );
		} 
		bg = ERROR_free( bg );
	}
	return NULL;
}

/*!  FUNCTION:  HMM_BG_SetSequence()
 *   SYNOPSIS:  Set the sequence to create digitized sequence.
 */
void HMM_BG_SetSequence( 	HMM_BG*		bg,
                           SEQUENCE* 	seq )
{
	/* if holding old sequence, get rid of it */
	if ( bg->sq != NULL ) esl_sq_Destroy( bg->sq );
	/* create digitized sequence */
	bg->sq 	= esl_sq_CreateDigital(bg->abc);
	esl_abc_CreateDsq(bg->abc, seq->seq, &bg->sq->dsq);
	// esl_sq_CreateDigitalFrom(abc, name, dsq, n, desc, acc, ss); /* Not necessary since we only use this struct for this function */
}

/*!  FUNCTION:  HMM_BG_UnsetSequence()
 *   SYNOPSIS:  Set the sequence to create digitized sequence.
 */
void HMM_BG_UnsetSequence( 	HMM_BG*		bg,
                              SEQUENCE* 	seq )
{
	// free( &bg->sq->dsq );
	esl_sq_Destroy(bg->sq);
}

/*!  FUNCTION:  HMM_BG_SetLength()
 *   SYNOPSIS:  Set the null model length distribution.
 *
 *   PURPOSE:  Sets the geometric null model length
 *             distribution in <bg> to a mean of <L> residues.
 */
void HMM_BG_SetLength( 	HMM_BG*		bg,
                        int 			L )
{
	bg->p1 = (float) L / (float) (L + 1);

	bg->fhmm->t[0][0] = bg->p1;
	bg->fhmm->t[0][1] = 1.0f - bg->p1;
}

/*!  FUNCTION: p7_bg_NullOne()
 *
 *   PURPOSE:  Calculate the null1 log score, for sequence <dsq>
 *             of length <L> "aligned" to the base null model <bg>.
 *
 *   NOTE:     Because the residue composition in null1 <bg> is the
 *             same as the background used to calculate residue
 *             scores in profiles and null models, all we have to
 *             do here is score null model transitions.
 *
 *             Can accept a NULL for *dsq, in which case the returned
 *             value will be (float) L * log(bg->p1) + log(1.-bg->p1);
 */
void 
HMM_BG_NullOne( 	const HMM_BG* 		bg,
                  int 					L,
                  float*				null1_sc )
{
	*null1_sc = (float) L * log( bg->p1 ) + log( 1.0 - bg->p1 );
}

/*!  FUNCTION:  p7_bg_SetFilter()
 *   SYNOPSIS:  Configure filter HMM with new model composition.
 *
 *   PURPOSE: The "filter HMM" is an experimental filter in the
 *            acceleration pipeline for avoiding biased composition
 *            sequences. It has no effect on final scoring, if a
 *            sequence passes all steps of the pipeline; it is only
 *            used to eliminate biased sequences from further
 *            consideration early in the pipeline, before the big guns
 *            of domain postprocessing are applied.
 *
 *            At least at present, it doesn't actually work as well as
 *            one would hope.  This will be an area of future work.
 *            What we really want to do is make a better null model of
 *            real protein sequences (and their biases), and incorporate
 *            that model into the flanks (NCJ states) of the profile.
 *
 *            <compo> is the average model residue composition, from
 *            either the HMM or the copy in a profile or optimized
 *            profile. <M> is the length of the model in nodes.
 *
 *            The expected length of the filter HMM's generated
 *            sequence is set to a default (about 400). You need a
 *            subsequent call to <p7_bg_SetLength()> to set it to the
 *            target sequence length. In hmmscan, this requires a
 *            call after every new model is read and <p7_pli_NewModel()>
 *            is called, because <NewModel()> is calling <p7_bg_SetFilter()>
 *            to copy the new model's composition <compo>. [Failure to
 *            do this properly was bug #h85, 14 Dec 2010.]
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      SRE:J4/25: generalized to use composition vector, not
 *                       specifically an HMM.
 *
 * Note:      This looks like a two-state HMM, but if you start thinking
 *            about its length distribution ("oh my god, L0 assumes a
 *            fixed L=400 expectation, it's all wrong, it's not
 *            conditional on the target sequence length and length
 *            modeling's messed up!"), don't panic. It's set up as a
 *            conditional-on-L model that generates according to P(x |
 *            model, L) P(L); the P(L) term is added in
 *            p7_bg_FilterScore() below.
 *
 *            Additionally, and not to confuse you further, but the
 *            t[0][0] transition is dependent on L.  The initial
 *            setting here is just a dummy. When p7_bg_SetLength()
 *            sets p1 for the null1 model length distribution, it sets
 *            t[0][0] to the same thing. This is controlling the
 *            relative expected balance of background sequence to
 *            biased sequence, not the overall length distribution.
 *
 *            All of this is ad hoc, and little of it has been
 *            optimized against data.
 */
void 
HMM_BG_SetFilter(	   HMM_BG* 			bg,
                     int 				M,
                     const float*	compo )
{
	float L0 = 400.0;						/* mean length in state 0 of filter HMM (normal background) */
	float L1 = (float) M / 8.0; 		/* mean length in state 1 of filter HMM (biased segment) */

	/* State 0 is the normal iid model. */
	bg->fhmm->t[0][0] =   L0 / (L0 + 1.0f);
	bg->fhmm->t[0][1] = 1.0f / (L0 + 1.0f);
	bg->fhmm->t[0][2] = 1.0f;          	/* 1.0 transition to E means we'll set length distribution externally. */
	esl_vec_FCopy(bg->f, bg->abc->K, bg->fhmm->e[0]);

	/* State 1 is the potentially biased model composition. */
	bg->fhmm->t[1][0] = 1.0f / (L1 + 1.0f);
	bg->fhmm->t[1][1] =   L1 / (L1 + 1.0f);
	bg->fhmm->t[1][2] = 1.0f;         	/* 1.0 transition to E means we'll set length distribution externally. */
	esl_vec_FCopy(compo, bg->abc->K, bg->fhmm->e[1]);

	bg->fhmm->pi[0] = 0.999;
	bg->fhmm->pi[1] = 0.001;

	esl_hmm_Configure(bg->fhmm, bg->f);
}

/*! FUNCTION: 	BG_FilterScore()
 *  SYNOPSIS:  Calculates the filter null model score.
 *
 *  PURPOSE:   Modeled after p7_bg_FilterScore() in HMMER.
 * 				Calculates the filter null model <bg> score for sequence
 *            	<dsq> of length <L>, and return it in
 *            	<*ret_sc>.
 *
 *            	The score is calculated as an HMM Forward score using
 *            	the two-state filter null model. It is a log-odds ratio,
 *            	relative to the iid background frequencies, in nats:
 *            	same as main model Forward scores.
 *
 *            	The filter null model has no length distribution of its
 *            	own; the same geometric length distribution (controlled
 *            	by <bg->p1>) that the null1 model uses is imposed.
 */
void 
HMM_BG_FilterScore( 	HMM_BG* 			bg,
                     SEQUENCE* 		seq,
                     float*			result )
{
	ESL_DSQ* 	dsq 	= bg->sq->dsq;
	int 			L 		= seq->N;
	int 			M 		= bg->fhmm->M;

	/* build matrix for computation */
	ESL_HMX* hmx 	= esl_hmx_Create(L, M); /* optimization target: this can be a 2-row matrix, and it can be stored in <bg>. */

	/* compute null_sc (composition bias) */
	float null_sc;
	esl_hmm_Forward(bg->sq->dsq, L, bg->fhmm, hmx, &null_sc);
	
	/* impose the length distribution */
	*result = null_sc + (float) L * logf(bg->p1) + logf(1.0 - bg->p1);

	esl_hmx_Destroy(hmx);
}