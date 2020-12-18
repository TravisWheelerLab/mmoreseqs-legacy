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
#include "../utilities/utilities.h"
#include "objects.h"

/* header */
#include "hmm_bg.h"

/* NOTE: Modeled after P7_BG: p7_bg.c */

/*
 *  FUNCTION:  HMM_BG_Create()
 *  SYNOPSIS:
 */
HMM_BG* HMM_BG_Create()
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

/*
 *  FUNCTION:  HMM_BG_Destroy()
 *  SYNOPSIS:
 */
void* HMM_BG_Destroy( HMM_BG* bg )
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

/* NOTE: Will be moved to SEQUENCE */
/* FUNCTION:  HMM_BG_SetSequence()
 * SYNOPSIS:  Set the sequence to create digitized sequence.
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

/* FUNCTION:  HMM_BG_UnsetSequence()
 * SYNOPSIS:  Set the sequence to create digitized sequence.
 */
void HMM_BG_UnsetSequence( 	HMM_BG*		bg,
                           SEQUENCE* 	seq )
{
	// free( &bg->sq->dsq );
	esl_sq_Destroy(bg->sq);
}


/* FUNCTION:  HMM_BG_SetLength()
 * SYNOPSIS:  Set the null model length distribution.
 *
 * 	PURPOSE:  Sets the geometric null model length
 *            distribution in <bg> to a mean of <L> residues.
 */
void HMM_BG_SetLength( 	HMM_BG*		bg,
                        int 			L )
{
	bg->p1 = (float) L / (float) (L + 1);

	bg->fhmm->t[0][0] = bg->p1;
	bg->fhmm->t[0][1] = 1.0f - bg->p1;
}

/* Function:  p7_bg_NullOne()
 *
 * Purpose:   Calculate the null1 lod score, for sequence <dsq>
 *            of length <L> "aligned" to the base null model <bg>.
 *
 * Note:      Because the residue composition in null1 <bg> is the
 *            same as the background used to calculate residue
 *            scores in profiles and null models, all we have to
 *            do here is score null model transitions.
 *
 *            Can accept a NULL for *dsq, in which case the returned
 *            value will be (float) L * log(bg->p1) + log(1.-bg->p1);
 */
void HMM_BG_NullOne( 	const HMM_BG* 		bg,
                        int 					L,
                        float*				null_sc )
{
	*null_sc = (float) L * log( bg->p1 ) + log( 1.0 - bg->p1 );
}

/* Function:  p7_bg_SetFilter()
 * Synopsis:  Configure filter HMM with new model composition.
 *
 * Purpose:   The "filter HMM" is an experimental filter in the
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
 *            conditional-on-L model that generates accordint to P(x |
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
void HMM_BG_SetFilter(	HMM_BG* 			bg,
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

/*
 * FUNCTION: 	BG_FilterScore()
 * SYNOPSIS:  	Calculates the filter null model score.
 *
 * PURPOSE:   	Modeled after p7_bg_FilterScore() in HMMER.
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
void HMM_BG_FilterScore( 	HMM_BG* 			bg,
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


/* Function:  p7_domaindef_ByPosteriorHeuristics()
 * Synopsis:  Define domains in a sequence using posterior probs.
 * Incept:    SRE, Sat Feb 23 08:17:44 2008 [Janelia]
 *
 * Purpose:   Given a sequence <sq> and model <om> for which we have
 *            already calculated a Forward and Backward parsing
 *            matrices <oxf> and <oxb>; use posterior probability
 *            heuristics to determine an annotated domain structure;
 *            and for each domain found, score it (with null2
 *            calculations) and obtain an optimal accuracy alignment,
 *            using <fwd> and <bck> matrices as workspace for the
 *            necessary full-matrix DP calculations. Caller provides a
 *            new or reused <ddef> object to hold these results.
 *            A <bg> is provided for (possible) use in biased-composition
 *            score correction (used in nhmmer), and a boolean
 *            <long_target> argument is provided to allow nhmmer-
 *            specific modifications to the behavior of this function
 *            (TRUE -> from nhmmer).
 *
 *            Upon return, <ddef> contains the definitions of all the
 *            domains: their bounds, their null-corrected Forward
 *            scores, and their optimal posterior accuracy alignments.
 *
 * Returns:   <eslOK> on success.
 *
 *            <eslERANGE> on numeric overflow in posterior
 *            decoding. This should not be possible for multihit
 *            models.
 */
// int
// p7_domaindef_ByPosteriorHeuristics( const ESL_SQ *sq, 
// 												const ESL_SQ *ntsq, 
// 												P7_OPROFILE *om,
//                                    	P7_OMX *oxf, P7_OMX *oxb, P7_OMX *fwd, P7_OMX *bck,
//                                    P7_DOMAINDEF *ddef, P7_BG *bg, int long_target,
//                                    P7_BG *bg_tmp, float *scores_arr, float *fwd_emissions_arr) 
// {
// 	int i, j;
//    int triggered;
//    int d;
//    int i2, j2;
//    int last_j2;
//    int nc;
//    int saveL     = om->L;   /* Save the length config of <om>; will restore upon return */
//    int save_mode = om->mode;   /* Likewise for the mode. */
//    int status;

//    // p7_domaindef_DumpPosteriors( stdout, ddef);

//    if ((status = p7_domaindef_GrowTo(ddef, sq->n))      != eslOK) return status;  /* ddef's btot,etot,mocc now ready for seq of length n */
//    if ((status = p7_DomainDecoding(om, oxf, oxb, ddef)) != eslOK) return status;  /* ddef->{btot,etot,mocc} now made.                    */

//    esl_vec_FSet(ddef->n2sc, sq->n + 1, 0.0);        /* ddef->n2sc null2 scores are initialized                        */
//    ddef->nexpected = ddef->btot[sq->n];             /* posterior expectation for # of domains (same as etot[sq->n])   */

//    p7_oprofile_ReconfigUnihit(om, saveL);     /* process each domain in unihit mode, regardless of om->mode     */
//    i     = -1;
//    triggered = FALSE;
// }

/*
 *      NOTE: HOW TO CONVERT row-coords to diag-coords
 *       MMX3(i-1,j-1) => MMX3(, d_2)
 *       MMX3(i,  j-1) => MMX3(, d_1)
 *       MMX3(i,  j  ) => MMX3(, d_1)
 */

/*  FUNCTION:  COMPUTE_Bias_Compo()
 *  SYNOPSIS:
 *
 *    RETURN:
 */
int COMPUTE_Bias_Compo( WORKER* worker )
{
   SEQUENCE*      q_seq    = worker->q_seq;
   HMM_PROFILE*   t_prof   = worker->t_prof;
   VECTOR_FLT*    null2sc  = VECTOR_FLT_Create();

   VECTOR_FLT_Set_Size( null2sc, q_seq->N + 1 );
   VECTOR_FLT_Fill( null2sc, 0.0 );

}


/* Function:  p7_GNull2_ByExpectation()
 * Synopsis:  Calculate null2 model from posterior probabilities.
 * Incept:    SRE, Thu Feb 28 09:52:28 2008 [Janelia]
 *
 * Purpose:   Calculate the "null2" model for the envelope encompassed
 *            by a posterior probability calculation <pp> for model
 *            <gm>.  Return the null2 odds emission probabilities
 *            $\frac{f'{x}}{f{x}}$ in <null2>, which caller
 *            provides as space for at least <alphabet->Kp> residues.
 *
 *            The expectation method is applied to envelopes in
 *            simple, well resolved regions (regions containing just a
 *            single envelope, where no stochastic traceback
 *            clustering was required).
 *
 *            Make sure that the posterior probability matrix <pp> has
 *            been calculated by the caller for only the envelope; thus
 *            its rows are numbered <1..Ld>, for envelope <ienv..jenv>
 *            of length <Ld=jenv-ienv+1>.
 *
 * Args:      gm    - profile, in any mode, target length model set to <L>
 *            pp    - posterior prob matrix, for <gm> against domain envelope <dsq+i-1> (offset)
 *            null2 - RETURN: null2 odds ratios per residue; <0..Kp-1>; caller allocated space
 *
 * Returns:   <eslOK> on success; <null2> contains the null2 scores. The 0
 *            row of <pp> has been used as temp space, and happens to contain
 *            the expected frequency that each M,I,N,C,J state is used in this
 *            <pp> matrix to generate residues.
 *
 * Throws:    (no abnormal error conditions)
 */
int
NULL2_ByExpectation(    HMM_PROFILE*      prof,
                        int               i,
                        int               j,
                        VECTOR_FLT*       null2 )
{
   // int      M      = gm->M;
   // int      Ld     = pp->L;
   // float  **dp     = pp->dp;
   // float   *xmx    = pp->xmx;
   // float    xfactor;
   // int      x;        /* over symbols 0..K-1                       */
   // int      i;        /* over offset envelope dsq positions 1..Ld  */
   // int      k;        /* over model M states 1..M, I states 1..M-1 */

   // /* Calculate expected # of times that each emitting state was used
   //  * in generating the Ld residues in this domain.
   //  * The 0 row in <wrk> is used to hold these numbers.
   //  */
   // esl_vec_FCopy(pp->dp[1],            (M + 1)*p7G_NSCELLS, pp->dp[0]);
   // esl_vec_FCopy(pp->xmx + p7G_NXCELLS,  p7G_NXCELLS,       pp->xmx);
   // for (i = 2; i <= Ld; i++)
   // {
   //    esl_vec_FAdd(pp->dp[0], pp->dp[i],             (M + 1)*p7G_NSCELLS);
   //    esl_vec_FAdd(pp->xmx,   pp->xmx + i * p7G_NXCELLS, p7G_NXCELLS);
   // }

   // /* Convert those expected #'s to log frequencies; these we'll use as
   //  * the log posterior weights.
   //  */
   // esl_vec_FLog(pp->dp[0], (M + 1)*p7G_NSCELLS);
   // esl_vec_FLog(pp->xmx,   p7G_NXCELLS);

   // esl_vec_FIncrement(pp->dp[0], (M + 1)*p7G_NSCELLS, -log((float)Ld));
   // esl_vec_FIncrement(pp->xmx,   p7G_NXCELLS,       -log((float)Ld));

   // /* Calculate null2's log odds emission probabilities, by taking
   //  * posterior weighted sum over all emission vectors used in paths
   //  * explaining the domain.
   //  * This is dog-slow; a point for future optimization.
   //  */
   // xfactor = XMX(0, p7G_N);
   // xfactor = p7_FLogsum(xfactor, XMX(0, p7G_C));
   // xfactor = p7_FLogsum(xfactor, XMX(0, p7G_J));
   // esl_vec_FSet(null2, gm->abc->K, -eslINFINITY);
   // for (x = 0; x < gm->abc->K; x++)
   // {
   //    for (k = 1; k < M; k++)
   //    {
   //       null2[x] = p7_FLogsum(null2[x], MMX(0, k) + p7P_MSC(gm, k, x));
   //       null2[x] = p7_FLogsum(null2[x], IMX(0, k) + p7P_ISC(gm, k, x));
   //    }
   //    null2[x] = p7_FLogsum(null2[x], MMX(0, M) + p7P_MSC(gm, k, x));
   //    null2[x] = p7_FLogsum(null2[x], xfactor);
   // }

   // esl_vec_FExp (null2, gm->abc->K);
   // /* now null2[x] = \frac{f_d(x)}{f_0(x)} for all x in alphabet,
   //  * 0..K-1, where f_d(x) are the ad hoc "null2" residue frequencies
   //  * for this envelope.
   //  */

   // /* make valid scores for all degeneracies, by averaging the odds ratios. */
   // esl_abc_FAvgScVec(gm->abc, null2); /* does not set gap, nonres, missing  */
   // null2[gm->abc->K]    = 1.0;        /* gap character    */
   // null2[gm->abc->Kp - 2] = 1.0;    /* nonresidue "*"   */
   // null2[gm->abc->Kp - 1] = 1.0;    /* missing data "~" */

   // return eslOK;
}