/*******************************************************************************
 *  FILE:      work_posterior.h
 *  PURPOSE:   Workflow Subroutines for Posterior Algorithms.
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
#include <time.h>

/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../parsers/_parsers.h"
#include "../algs_linear/_algs_linear.h"
#include "../algs_quad/_algs_quad.h"
#include "../algs_naive/_algs_naive.h"
#include "../algs_sparse/_algs_sparse.h"
#include "../reporting/_reporting.h"

/* easel library */
#include "easel.h"
#include "esl_gumbel.h"
#include "esl_exponential.h"

/* header */
#include "_work.h"
#include "work_index.h"

/* compute correction bias and convert natscore -> bitscore -> pval -> eval */
void 
WORK_posterior( WORKER* worker )
{
   FILE*          fp       = NULL;

   HMM_BG*        bg       = worker->hmm_bg;
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      q_seq    = worker->q_seq;
   RESULT*        result   = worker->result;
   TASKS*         tasks    = worker->tasks;

   /* size */
   int      T           = worker->t_prof->N;
   int      Q           = worker->q_seq->N;
   /* alignment scores */
   float    sc;
   float    nat_sc      = 0.0f;  /* score in NATS */
   float    pre_sc      = 0.0f;  /* adjusted for compo bias / score in BITS */
   float    seq_sc      = 0.0f;  /* adjusted for sequence bias / score in BITS */
   float    ln_pval     = 0.0f;  /* natural log of p-value */
   float    pval        = 0.0f;  /* p-value */
   float    eval        = 0.0f;  /* e-value */
   /* bias correction */
   float    null_sc     = 0.0f;  /* in NATS */
   float    filter_sc   = 0.0f;  /* in NATS */
   float    seq_bias    = 0.0f;  /* in NATS */
   /* parameters for exponential distribution, for converting bitscore -> p-value */
   float    tau         = worker->t_prof->forward_dist.param1;
   float    lambda      = worker->t_prof->forward_dist.param2;
   /* number of sequences in database, for computing e-value */
   int      n_seqs      = worker->q_index->N;

   /* init scores */
   seq_bias = 0.0f;
   nat_sc   = worker->scores->lin_bound_fwd;

   /* initialize hmm_bg */
   HMM_BG_SetSequence( bg, q_seq );
   HMM_BG_SetFilter( bg, t_prof->N, t_prof->bg_model->compo );
   HMM_BG_SetLength( bg, q_seq->N );
   /* compute null one */
   HMM_BG_NullOne( bg, q_seq->N, &null_sc );
   /* compute nullscore for bias */
   HMM_BG_FilterScore( bg, q_seq, &filter_sc );
   /* free digitized sequence TODO: move to sequence */
   // HMM_BG_UnsetSequence( bg, seq );

   /* Posterior and Find Domains */
   if ( worker->args->is_compo_bias )
   {
      /*! NOTE: Find domains and assesses domain-specific correction bias
       *  See p7_domaindef_ByPosteriorHeuristics().
       *  Currently, there is no domain finding. The search cloud is considered 
       *  a single complete domain and bias is assessed for the entire region. 
       */

      /* For testing: Quadratic assesses bias for the entire tightest bounding box containing entire cloud. */ 
      if ( tasks->quad_bias_corr == true )
      {
         run_Posterior_Quad(
            worker->q_seq, worker->t_prof, Q, T, worker->hmm_bg, worker->edg_row,
            worker->st_MX_fwd, worker->sp_MX_fwd, worker->st_MX_bck, worker->sp_MX_bck, worker->st_MX_bck, worker->sp_MX_bck, 
            worker->dom_def );
      }
      /* Sparse assesses bias only for cells contained by cloud */ 
      if ( tasks->sparse_bias_corr == true )
      {
         /* if running full fwdbackward or pruned (for comparison testing) */
         if ( worker->args->is_run_full == true ) {
            /* cloud fills entire dp matrix */
            EDGEBOUNDS_Cover_Matrix(worker->edg_row, Q, T);
            result->cloud_cells  = EDGEBOUNDS_Count( worker->edg_row );
         }

         /* build sparse matrices */
         CLOCK_Start( worker->clok );
         MATRIX_3D_SPARSE_Shape_Like_Edgebounds( worker->st_SMX_fwd, worker->edg_row );
         MATRIX_3D_SPARSE_Fill_Outer( worker->st_SMX_fwd, -INF );
         MATRIX_3D_SPARSE_Copy( worker->st_SMX_bck, worker->st_SMX_fwd );
         if ( worker->st_SMX_post != worker->st_SMX_bck ) {
            MATRIX_3D_SPARSE_Copy( worker->st_SMX_post, worker->st_SMX_fwd );
         }
         CLOCK_Stop( worker->clok );
         worker->times->sp_build_mx = CLOCK_Duration( worker->clok );

         /* TODO: Need to break up posterior into WORK subroutines */
         /* run full posterior */
         bool is_run_domains = worker->args->is_run_domains;
         CLOCK_Start( worker->clok );
         /* TODO: remove timer from this */
         run_Posterior_Sparse(
            worker->q_seq, worker->t_prof, Q, T, worker->hmm_bg, worker->edg_row,
            worker->st_SMX_fwd, worker->sp_MX_fwd, worker->st_SMX_bck, worker->sp_MX_bck, 
            worker->st_SMX_post, worker->sp_MX_post, worker->st_SMX_fwd, worker->sp_MX_fwd,
            worker->result, worker->dom_def, worker->clok, worker->times, is_run_domains );
         CLOCK_Stop( worker->clok );
         worker->times->sp_posterior = CLOCK_Duration( worker->clok );
      }

      /* compute sequence bias */
      nat_sc   = result->bound_fwd_natsc;
      seq_bias = result->final_scores.seq_bias;
      seq_bias = logsum(0.0, worker->hmm_bg->omega + worker->dom_def->seq_bias);
      printf("nat_sc: %7.3f, null_sc: %7.3f, (final) seq_bias: %7.3f,\n", nat_sc, null_sc, seq_bias);
      seq_bias = worker->dom_def->seq_bias;
   }

   /* compute pre_score and sequence_score by accounting for bias and convert from nats -> bits */
   pre_sc = (nat_sc - null_sc) / eslCONST_LOG2;
   seq_sc = (nat_sc - (null_sc + seq_bias)) / eslCONST_LOG2;
   fprintf(stdout, "# nat_sc = %11.4f, null_sc = %11.4f, seq_bias = %7.4f\n",
      nat_sc, null_sc, seq_bias );
   fprintf(stdout, "# pre_sc = %7.4f, seq_sc = %7.4f\n", 
      pre_sc, seq_sc );

   /* compute log of P-value */ 
   ln_pval  = esl_exp_logsurv( seq_sc, tau, lambda );
   pval     = exp(ln_pval);
   eval     = pval * n_seqs;
   fprintf(stdout, "# ln_pval = %7.4f, pval = %9.2e, eval = %9.2e\n",
      ln_pval, pval, eval );

   /* save scores */
   result->final_scores.nat_sc      = nat_sc;
   result->final_scores.null_sc     = null_sc;
   result->final_scores.filter_sc   = filter_sc;
   result->final_scores.seq_bias    = seq_bias;
   result->final_scores.pre_sc      = pre_sc;
   result->final_scores.seq_sc      = seq_sc;
   result->final_scores.pre_sc      = pre_sc;
   result->final_scores.ln_pval     = ln_pval;
   result->final_scores.pval        = pval;
   result->final_scores.eval        = eval;
}

