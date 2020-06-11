/*******************************************************************************
 *  FILE:      p7_funcs.c
 *  PURPOSE:   Imports HMMER and Easel library functions.
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
#include <ctype.h>

/* import vector intrinsics */
#include <xmmintrin.h>    /* SSE  */
#include <emmintrin.h>    /* SSE2 */
#include <pmmintrin.h>   /* DENORMAL_MODE */

/* import easel */
#include "../easel/easel.h"
#include "../easel/esl_dirichlet.h"
#include "../easel/esl_hmm.h"
#include "../easel/esl_sq.h"
#include "../easel/esl_sqio.h"
#include "../easel/esl_scorematrix.h"
#include "../easel/esl_getopts.h"
#include "../easel/esl_vectorops.h"

/* import hmmer */
#include "hmmer/p7_structs.h"

/* self header */
#include "hmmer/p7_funcs.h"

/* Function:  p7_bg_Create()
 * Synopsis:  Create a <P7_BG> null model object.
 *
 * Purpose:   Allocate a <P7_BG> object for digital alphabet <abc>,
 *            initializes it to appropriate default values, and
 *            returns a pointer to it.
 *            
 *            For protein models, default iid background frequencies
 *            are set (by <p7_AminoFrequencies()>) to average
 *            Swiss-Prot residue composition. For DNA, RNA and other
 *            alphabets, default frequencies are set to a uniform
 *            distribution.
 *            
 *            The model composition <bg->mcomp[]> is not initialized
 *            here; neither is the filter null model <bg->fhmm>.  To
 *            use the filter null model, caller will want to
 *            initialize these fields by calling
 *            <p7_bg_SetFilter()>.
 *
 * Throws:    <NULL> on allocation failure.
 *
 * Xref:      STL11/125.
 */
P7_BG *
p7_bg_Create(const ESL_ALPHABET *abc)
{
  P7_BG *bg = NULL;
  int    status;

  ESL_ALLOC(bg, sizeof(P7_BG));
  bg->f     = NULL;
  bg->fhmm  = NULL;

  ESL_ALLOC(bg->f,     sizeof(float) * abc->K);
  if ((bg->fhmm = esl_hmm_Create(abc, 2)) == NULL) goto ERROR;

  if       (abc->type == eslAMINO)
    {
      if (p7_AminoFrequencies(bg->f) != eslOK) goto ERROR;
    }
  else
    esl_vec_FSet(bg->f, abc->K, 1. / (float) abc->K);

  bg->p1    = 350./351.;
  bg->omega = 1./256.;
  bg->abc   = abc;
  return bg;

 ERROR:
  p7_bg_Destroy(bg);
  return NULL;
}


/* Function:  p7_bg_Destroy()
 *
 * Purpose:   Frees a <P7_BG> object.
 *
 * Returns:   (void)
 *
 * Xref:      SRE:STL11/125.
 */
void
p7_bg_Destroy(P7_BG *bg)
{
  if (bg != NULL) {
    if (bg->f     != NULL) free(bg->f);
    if (bg->fhmm  != NULL) esl_hmm_Destroy(bg->fhmm);
    free(bg);
  }
  return;
}


/* Function:  p7_builder_Create()
 * Synopsis:  Create a default HMM construction configuration.
 *
 * Purpose:   Create a construction configuration for building
 *            HMMs in alphabet <abc>, and return a pointer to it.
 *            
 *            An application configuration <go> may optionally be
 *            provided. If <go> is <NULL>, default parameters are
 *            used. If <go> is non-<NULL>, it must include appropriate
 *            settings for all of the following ``standard build options'':
 *            
 *            Model construction:   --fast --hand --symfrac --fragthresh
 *            Relative weighting:   --wgsc --wblosum --wpb --wgiven --wid
 *            Effective seq #:      --eent --eclust --enone --eset --ere --esigma --eid
 *            Prior scheme:         --pnone --plaplace
 *            E-val calibration:    --EmL --EmN --EvL --EvN --EfL --EfN --Eft
 *            run-to-run variation: --seed
 *            
 *            See <hmmbuild.c> or other big users of the build
 *            pipeline for an example of appropriate <ESL_GETOPTS>
 *            initializations of these 24 options.
 */
P7_BUILDER *
p7_builder_Create(const ESL_GETOPTS*   go, 
                  const ESL_ALPHABET*  abc)
{
  P7_BUILDER *bld = NULL;
  int         seed;
  int         status;


  ESL_ALLOC(bld, sizeof(P7_BUILDER));
  bld->prior        = NULL;
  bld->r            = NULL;
  bld->S            = NULL;
  bld->Q            = NULL;
  bld->eset         = -1.0;   /* -1.0 = unset; must be set if effn_strategy is p7_EFFN_SET */
  bld->re_target    = -1.0;

  if (go == NULL) 
    {
      bld->arch_strategy = p7_ARCH_FAST;
      bld->wgt_strategy  = p7_WGT_PB;
      bld->effn_strategy = p7_EFFN_ENTROPY;
      seed               = 42;
    }
  else 
    {
      if      (esl_opt_GetBoolean(go, "--fast"))    bld->arch_strategy = p7_ARCH_FAST;
      else if (esl_opt_GetBoolean(go, "--hand"))    bld->arch_strategy = p7_ARCH_HAND;

    
      if      (esl_opt_GetBoolean(go, "--wpb"))     bld->wgt_strategy = p7_WGT_PB;
      else if (esl_opt_GetBoolean(go, "--wgsc"))    bld->wgt_strategy = p7_WGT_GSC;
      else if (esl_opt_GetBoolean(go, "--wblosum")) bld->wgt_strategy = p7_WGT_BLOSUM;
      else if (esl_opt_GetBoolean(go, "--wnone"))   bld->wgt_strategy = p7_WGT_NONE;
      else if (esl_opt_GetBoolean(go, "--wgiven"))  bld->wgt_strategy = p7_WGT_GIVEN;

      if      (esl_opt_GetBoolean(go, "--eent"))    bld->effn_strategy = p7_EFFN_ENTROPY;
      else if (esl_opt_GetBoolean(go, "--eentexp")) bld->effn_strategy = p7_EFFN_ENTROPY_EXP;
      else if (esl_opt_GetBoolean(go, "--eclust"))  bld->effn_strategy = p7_EFFN_CLUST;
      else if (esl_opt_GetBoolean(go, "--enone"))   bld->effn_strategy = p7_EFFN_NONE;
      else if (esl_opt_IsOn      (go, "--eset"))  { bld->effn_strategy = p7_EFFN_SET;      bld->eset = esl_opt_GetReal(go, "--eset"); }

      seed = esl_opt_GetInteger(go, "--seed");
    }

  bld->max_insert_len = 0;

  /* The default RE target is alphabet dependent. */
  if (go != NULL &&  esl_opt_IsOn (go, "--ere")) 
    bld->re_target = esl_opt_GetReal(go, "--ere");
  else {
    switch (abc->type) {
    case eslAMINO:  bld->re_target = p7_ETARGET_AMINO; break;
    case eslDNA:    bld->re_target = p7_ETARGET_DNA;   break;
    case eslRNA:    bld->re_target = p7_ETARGET_DNA;   break;
    default:        bld->re_target = p7_ETARGET_OTHER; break;
    }
  }

  bld->symfrac    = (go != NULL) ?  esl_opt_GetReal   (go, "--symfrac")    : 0.5; 
  bld->fragthresh = (go != NULL) ?  esl_opt_GetReal   (go, "--fragthresh") : 0.5; 
  bld->wid        = (go != NULL) ?  esl_opt_GetReal   (go, "--wid")        : 0.62;
  bld->esigma     = (go != NULL) ?  esl_opt_GetReal   (go, "--esigma")     : 45.0;
  bld->eid        = (go != NULL) ?  esl_opt_GetReal   (go, "--eid")        : 0.62;
  bld->EmL        = (go != NULL) ?  esl_opt_GetInteger(go, "--EmL")        : 200;
  bld->EmN        = (go != NULL) ?  esl_opt_GetInteger(go, "--EmN")        : 200;
  bld->EvL        = (go != NULL) ?  esl_opt_GetInteger(go, "--EvL")        : 200;
  bld->EvN        = (go != NULL) ?  esl_opt_GetInteger(go, "--EvN")        : 200;
  bld->EfL        = (go != NULL) ?  esl_opt_GetInteger(go, "--EfL")        : 100;
  bld->EfN        = (go != NULL) ?  esl_opt_GetInteger(go, "--EfN")        : 200;
  bld->Eft        = (go != NULL) ?  esl_opt_GetReal   (go, "--Eft")        : 0.04;

  /* Normally we reinitialize the RNG to original seed before calibrating each model.
   * This eliminates run-to-run variation.
   * As a special case, seed==0 means choose an arbitrary seed and shut off the
   * reinitialization; this allows run-to-run variation.
   */

  bld->r            = esl_randomness_CreateFast(seed);
  bld->do_reseeding = (seed == 0) ? FALSE : TRUE;

  if      (go && esl_opt_GetBoolean(go, "--pnone") )     bld->prior = NULL;
  else if (go && esl_opt_GetBoolean(go, "--plaplace") )  bld->prior = p7_prior_CreateLaplace(abc);
  else
    {
      switch (abc->type) {
      case eslAMINO: bld->prior = p7_prior_CreateAmino();      break;
      case eslDNA:   bld->prior = p7_prior_CreateNucleic();    break;
      case eslRNA:   bld->prior = p7_prior_CreateNucleic();    break;
      default:       bld->prior = p7_prior_CreateLaplace(abc); break;
      }
      if (bld->prior == NULL) goto ERROR;
    }


  bld->abc       = abc;
  bld->errbuf[0] = '\0';

  bld->popen   = -1;
  bld->pextend = -1;

  return bld;
  
 ERROR:
  p7_builder_Destroy(bld);
  return NULL;
}


/* Function:  p7_SingleBuilder()
 * Synopsis:  Build a new HMM from a single sequence.
 *
 * Purpose:   Take the sequence <sq> and a build configuration <bld>, and
 *            build a new HMM.
 *            
 *            The single sequence scoring system in the <bld>
 *            configuration must have been previously initialized by
 *            <p7_builder_SetScoreSystem()>.
 *            
 * Args:      bld       - build configuration
 *            sq        - query sequence
 *            bg        - null model (needed to paramaterize insert emission probs)
 *            opt_hmm   - optRETURN: new HMM
 *            opt_gm    - optRETURN: profile corresponding to <hmm>
 *            opt_om    - optRETURN: optimized profile corresponding to <gm>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if <bld> isn't properly configured somehow.
 */
int
p7_SingleBuilder( P7_BUILDER*    bld, 
                  ESL_SQ*        sq, 
                  P7_BG*         bg, 
                  P7_HMM**       opt_hmm,
                  P7_TRACE**     opt_tr, 
                  P7_PROFILE**   opt_gm, 
                  P7_OPROFILE**  opt_om )
{
  printf("=== SINGLEBUILDER ===\n");
  P7_HMM   *hmm = NULL;
  P7_TRACE *tr  = NULL;
  int       k;
  int       status;

  
  bld->errbuf[0] = '\0';
  if (! bld->Q) ESL_XEXCEPTION(eslEINVAL, "score system not initialized");

  if ((status = p7_Seqmodel(bld->abc, sq->dsq, sq->n, sq->name, bld->Q, bg->f, bld->popen, bld->pextend, &hmm)) != eslOK) goto ERROR;
  if ((status = p7_hmm_SetComposition(hmm))                                                                     != eslOK) goto ERROR;
  if ((status = p7_hmm_SetConsensus(hmm, sq))                                                                   != eslOK) goto ERROR; 
  if ((status = calibrate(bld, hmm, bg, opt_gm, opt_om))                                                        != eslOK) goto ERROR;

  if ( bld->abc->type == eslDNA ||  bld->abc->type == eslRNA ) 
  {
    if (bld->w_len > 0)           hmm->max_length = bld->w_len;
    else if (bld->w_beta == 0.0)  hmm->max_length = hmm->M *4;
    else if ( (status =  p7_Builder_MaxLength(hmm, bld->w_beta)) != eslOK) goto ERROR;
  }

  /* DAVID RICH EDIT */
  printf("=== HMM IN SINGLEBUILDER ===\n");
  FILE *fp = NULL;
  fp = fopen("examples/hmm_singlebuilder.txt", "w+");
  p7_hmm_Dump(fp, hmm);
  fclose(fp);

  /* note that <opt_gm> and <opt_om> were already set by calibrate() call above. */
  if (opt_hmm   != NULL) *opt_hmm = hmm; else p7_hmm_Destroy(hmm);
  if (opt_tr    != NULL) *opt_tr  = tr;
  return eslOK;

 ERROR:
  p7_hmm_Destroy(hmm);
  if (tr        != NULL) p7_trace_Destroy(tr);
  if (opt_gm    != NULL) p7_profile_Destroy(*opt_gm);
  if (opt_om    != NULL) p7_oprofile_Destroy(*opt_om);
  return status;
}

/* Function:  p7_Seqmodel()
 * Synopsis:  Make a profile HMM from a single sequence.
 *
 * Purpose:   Make a profile HMM from a single sequence, for
 *            probabilistic Smith/Waterman alignment, HMMER3-style.
 *            
 *            The query is digital sequence <dsq> of length <M>
 *            residues in alphabet <abc>, named <name>. 
 *            
 *            The scoring system is given by <Q>, <f>, <popen>, and
 *            <pextend>. <Q> is a $K \times K$ matrix giving
 *            conditional residue probabilities $P(a \mid b)}$; these
 *            are typically obtained by reverse engineering a score
 *            matrix like BLOSUM62. <f> is a vector of $K$ background
 *            frequencies $p_a$. <popen> and <pextend> are the
 *            probabilities assigned to gap-open ($t_{MI}$ and
 *            $t_{MD}$) and gap-extend ($t_{II}$ and $t_{DD}$)
 *            transitions.
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success, and a newly allocated HMM is returned
 *            in <ret_hmm>. 
 *
 * Throws:    <eslEMEM> on allocation error, and <*ret_hmm> is <NULL>.
 */
int
p7_Seqmodel( const ESL_ALPHABET* abc, 
                   ESL_DSQ*      dsq, 
                   int           M, 
                   char*         name,
                   ESL_DMATRIX*  Q, 
                   float*        f, 
                   double        popen, 
                   double        pextend,
                  P7_HMM**       ret_hmm)
{
  printf("=== SEQ_MODEL ===\n");
  int     status;
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[HMM created from a query sequence]";
  int     k;

  if ((hmm = p7_hmm_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; }
  
  for (k = 0; k <= M; k++)
    {
      /* Use rows of P matrix as source of match emission vectors */
      if (k > 0) esl_vec_D2F(Q->mx[(int) dsq[k]], abc->K, hmm->mat[k]);

      /* Set inserts to background for now. This will be improved. */
      esl_vec_FCopy(f, abc->K, hmm->ins[k]);

      hmm->t[k][p7H_MM] = 1.0 - 2 * popen;
      hmm->t[k][p7H_MI] = popen;
      hmm->t[k][p7H_MD] = popen;
      hmm->t[k][p7H_IM] = 1.0 - pextend;
      hmm->t[k][p7H_II] = pextend;
      hmm->t[k][p7H_DM] = 1.0 - pextend;
      hmm->t[k][p7H_DD] = pextend;
    }

  /* Deal w/ special stuff at node M, overwriting a little of what we
   * just did. 
   */
  hmm->t[M][p7H_MM] = 1.0 - popen;
  hmm->t[M][p7H_MD] = 0.;
  hmm->t[M][p7H_DM] = 1.0;
  hmm->t[M][p7H_DD] = 0.;
  
  /* Add mandatory annotation
   */
  p7_hmm_SetName(hmm, name);
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  hmm->nseq     = 1;
  p7_hmm_SetCtime(hmm);
  hmm->checksum = 0;

  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}

/* Function:  p7_hmm_SetComposition()
 * Synopsis:  Calculate and set model composition, <hmm->compo[]>
 *
 * Purpose:   Calculates the mean residue composition emitted by
 *            model <hmm>, and set <hmm->compo[]> to it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, in which case 
 *            values in <hmm->compo[]> are unchanged.
 *            
 * Note:      In principle, you should be able to normalize
 *            hmm->compo[] by dividing thru by the sum of
 *            mocc[] and iocc[] vectors, and that's what the
 *            3.0 release version did. This allowed p7_hmm_Validate()
 *            to check compo[] for summation to 1.0, as a smoke
 *            check for bugs here. The problem with that 
 *            is numerical roundoff error accumulation, when
 *            hmm->M is large [bug #h84]. To fix #h84, we
 *            simply renormalize compo[], rather than the fancier
 *            previous version. This avoids error accumulation,
 *            but it also guarantees that compo[] will trivially
 *            pass the hmm_Validation() step; it's not really
 *            validating the SetComposition() calculation at all.                                 
 *            (For description of #h84, error analysis, and the fix,
 *            xref J7/7; SRE, Tue Nov  2 14:32:29 2010)
 */
int
p7_hmm_SetComposition(P7_HMM *hmm)
{
  float *mocc = NULL;
  float *iocc = NULL;
  int    k;
  int    status;

  ESL_ALLOC(mocc, sizeof(float) * (hmm->M+1));
  ESL_ALLOC(iocc, sizeof(float) * (hmm->M+1));

  p7_hmm_CalculateOccupancy(hmm, mocc, iocc);
  esl_vec_FSet(hmm->compo, hmm->abc->K, 0.0);
  esl_vec_FAddScaled(hmm->compo, hmm->ins[0], iocc[0], hmm->abc->K);
  for (k = 1; k <= hmm->M; k++)
    {
      esl_vec_FAddScaled(hmm->compo, hmm->mat[k], mocc[k], hmm->abc->K);
      esl_vec_FAddScaled(hmm->compo, hmm->ins[k], iocc[k], hmm->abc->K);
    }

  esl_vec_FNorm(hmm->compo, hmm->abc->K);
  hmm->flags  |= p7H_COMPO;

  free(mocc);
  free(iocc);
  return eslOK;

 ERROR:
  if (mocc != NULL) free(mocc);
  if (iocc != NULL) free(iocc);
  return status;
}
  

/* Function:  p7_hmm_SetConsensus()
 * Synopsis:  Set the consensus residue line of the HMM.
 *
 * Purpose:   Sets the consensus annotation line of the model <hmm>.
 *            
 *            Behavior differs, depending on whether this is a
 *            single-sequence model (i.e. phmmer) or a standard
 *            model of a multiple sequence alignment. If <sq> is
 *            non-<NULL> this is a single-sequence model and <sq> is
 *            the digital sequence it was built from. If <sq> is <NULL>
 *            this is a standard multiple-sequence model.
 *            
 *            In a standard model, the most likely (highest emission
 *            probability) residue is the consensus at each position.
 *            In a single-sequence model, the consensus is the
 *            sequence itself.
 *            
 *            In both cases, if the emission probability is $\geq$
 *            certain threshold, the residue is upper cased. The
 *            threshold is arbitrarily set to 0.9 for nucleic acid
 *            alphabets (<eslDNA>, <eslRNA>) and 0.5 for amino acid
 *            alphabets (<eslAMINO>) and all other alphabets.
 *            
 *            The special handling of single-sequence models avoids
 *            a counterintuitive case where the most likely residue is
 *            not the original residue. For example, under the
 *            BLOSUM62 matrix, given an observed M, the most likely
 *            aligned residue is an L, not an M. (Because L is so much
 *            more likely a priori than M.)
 *
 * Args:      hmm    - model with valid probability parameters mat[1..M][x]
 *            sq     - NULL if a standard model;
 *                     or the query sequence for a single-sequence model.
 *           
 * Returns:   <eslOK> on success. The <p7H_CONS> flag on the <hmm> is raised
 *            if it wasn't already. The <hmm->consensus> line is set.
 *
 * Throws:    <eslEMEM> on allocation error. The <p7H_CONS> is dropped, even
 *            if it was up to begin with, and the <hmm->consensus> is <NULL>,
 *            even if we had one to begin with.
 *
 * Xref:      SRE:J8/26.
 */
int
p7_hmm_SetConsensus(P7_HMM *hmm, ESL_SQ *sq)
{
  int   k, x;
  float mthresh;
  int   status;
  
  /* allocation, if needed */
  if (! hmm->consensus) ESL_ALLOC(hmm->consensus, sizeof(char) * (hmm->M+2));

  /* set our arbitrary threshold for upper/lower casing */
  if      (hmm->abc->type == eslAMINO) mthresh = 0.5;
  else if (hmm->abc->type == eslDNA)   mthresh = 0.9;
  else if (hmm->abc->type == eslRNA)   mthresh = 0.9;
  else                                 mthresh = 0.5;

  hmm->consensus[0] = ' ';
  for (k = 1; k <= hmm->M; k++) 
    {
      x = (sq ?  sq->dsq[k] : esl_vec_FArgMax(hmm->mat[k], hmm->abc->K));
      hmm->consensus[k] = ((hmm->mat[k][x] >= mthresh) ? toupper(hmm->abc->sym[x]) : tolower(hmm->abc->sym[x]));
    }
  hmm->consensus[hmm->M+1] = '\0';
  hmm->flags  |= p7H_CONS; 
  return eslOK;

 ERROR:
  if (hmm->consensus) free(hmm->consensus);
  hmm->consensus = NULL;
  hmm->flags    &= (~p7H_CONS);  
  return status;
}

/* Function:  p7_hmm_CalculateOccupancy()
 * Synopsis:  Calculate match occupancy and insert expected use count vectors.
 *
 * Purpose:   Calculate a vector <mocc[1..M]> containing probability
 *            that each match state is used in a sampled glocal path through
 *            the model. Caller provides allocated space (<M+1> floats)
 *            for <mocc>.
 *            
 *            Caller may optionally provide an array <iocc[0..M]> as
 *            well, which (if provided) will be set to contain the
 *            expected number of times that a sampled path would contain
 *            each insert state.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_hmm_CalculateOccupancy(const P7_HMM *hmm, float *mocc, float *iocc)
{
  int k;

  mocc[0] = 0.;                              /* no M_0 state */
  mocc[1] = hmm->t[0][p7H_MI] + hmm->t[0][p7H_MM];   /* initialize w/ 1 - B->D_1 */
  for (k = 2; k <= hmm->M; k++)
      mocc[k] = mocc[k-1] * (hmm->t[k-1][p7H_MM] + hmm->t[k-1][p7H_MI]) +
        (1.0-mocc[k-1]) * hmm->t[k-1][p7H_DM];
  if (iocc != NULL) {
    iocc[0] = hmm->t[0][p7H_MI] / hmm->t[0][p7H_IM];
    for (k = 1; k <= hmm->M; k++)
      iocc[k] = mocc[k] * hmm->t[k][p7H_MI] / hmm->t[k][p7H_IM];
  }

  return eslOK;
}

