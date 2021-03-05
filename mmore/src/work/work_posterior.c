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
#include "work_posterior.h"

/*! FUNCTION:  WORK_posterior()
 *  SYNOPSIS:  Run full posterior, for full sequence and domain-specific. 
 */
void 
WORK_posterior( WORKER* worker )
{
   WORK_posterior_sparse( worker );
}

/*! FUNCTION:  WORK_posterior_sparse()
 *  SYNOPSIS:  Run posterior for full sequence using sparse matrices.
 */
void 
WORK_posterior_sparse( WORKER* worker )
{
   ARGS* args = worker->args;

   /* compute hmm model bias */
   WORK_null1_hmm_bias( worker );
   /* build sparse matrix */
   WORK_build_sparse_matrix( worker );
   /* compute fwdback */
   WORK_bound_fwdback_sparse( worker );
   /* find domains */
   WORK_decode_domains( worker );
   /* compute posterior */
   WORK_decode_posterior( worker );
   /* compute sequence bias */
   WORK_null2_seq_bias( worker );

   /* compute posterior alignment */
   if ( args->is_run_postaln == true ) 
   {
      /* compute optimal accuracy from posterior */
      WORK_optimal_accuracy( worker );
      /* backtrace optimal accuracy for posterior alignment */
      WORK_optacc_traceback( worker );
   }
   
   /* compute viterbi alignment */
   if ( args->is_run_vitaln == true )
   {
      /* compute viterbi */
      WORK_viterbi_sparse( worker );
      /* backtrace viterbi */
      WORK_viterbi_traceback_sparse( worker );
   }
}

/*! FUNCTION:  WORK_null1_hmm_bias()
 *  SYNOPSIS:  Compute the correction bias for the hmm model.
 */
void 
WORK_null1_hmm_bias( WORKER* worker )
{
   FILE*          fp             = NULL;
   TASKS*         tasks          = worker->tasks;
   ARGS*          args           = worker->args;
   CLOCK*         timer          = worker->timer;
   /* input data */
   SEQUENCE*      q_seq          = worker->q_seq;
   int            Q              = q_seq->N;
   HMM_PROFILE*   t_prof         = worker->t_prof;
   int            T              = t_prof->N;
   HMM_BG*        hmm_bg         = worker->hmm_bg;
   /* output data */
   TIMES*         times          = worker->times;
   RESULT*        result         = worker->result;
   ALL_SCORES*    scores         = &result->scores;
   SCORES*        finalsc        = &result->final_scores;
   DOMAIN_DEF*    dom_def        = worker->dom_def;

   /* temp */
   scores->null1_hmm_bias     = 0.0f;
   scores->null1_filtersc     = 0.0f;
   finalsc->null_omega_natsc  = 0.0f;
   scores->null_omega         = 0.0f;
   dom_def->null_omega        = 0.0f;
   dom_def->null1_hmm_bias    = 0.0f;
   
   /* NOTE: Currently filter score is not used for anything.
    *       In addition Set/UnsetSequence process has a memory leak
    *       when building ESL_DSQ.
    */
   /* initialize hmm_bg */
   // HMM_BG_SetSequence( hmm_bg, q_seq );
   HMM_BG_SetFilter( hmm_bg, t_prof->N, t_prof->bg_model->compo );
   HMM_BG_SetLength( hmm_bg, q_seq->N );
   /* compute null one */
   HMM_BG_NullOne( hmm_bg, q_seq->N, &scores->null1_hmm_bias );
   /* compute nullscore for bias */
   // HMM_BG_FilterScore( hmm_bg, q_seq, &scores->null1_filtersc );
   /* fetch omega (prior prob of no bias) */
   finalsc->null_omega_natsc  = hmm_bg->omega;
   scores->null_omega         = hmm_bg->omega;
   dom_def->null_omega        = hmm_bg->omega;
   dom_def->null1_hmm_bias    = scores->null1_hmm_bias;
   /* free digitized sequence TODO: move to sequence */
   // HMM_BG_UnsetSequence( hmm_bg, q_seq );
}

/*! FUNCTION:  WORK_decode_domains()
 *  SYNOPSIS:  Find domain ranges in posterior.
 */
void
WORK_decode_domains( WORKER* worker )
{
   FILE*          fp             = NULL;
   TASKS*         tasks          = worker->tasks;
   ARGS*          args           = worker->args;
   CLOCK*         timer          = worker->timer;
   /* input data */
   SEQUENCE*      q_seq          = worker->q_seq;
   int            Q              = q_seq->N;
   HMM_PROFILE*   t_prof         = worker->t_prof;
   int            T              = t_prof->N;
   EDGEBOUNDS*    edg            = worker->edg_row;
   /* working data */
   MATRIX_2D*     sp_MX_fwd      = worker->sp_MX_fwd;
   MATRIX_2D*     sp_MX_bck      = worker->sp_MX_bck;
   /* output data */
   TIMES*         times          = worker->times;
   RESULT*        result         = worker->result;
   ALL_SCORES*    scores         = &result->scores;
   SCORES*        final_scores   = &result->final_scores;
   DOMAIN_DEF*    dom_def        = worker->dom_def;
   STATS*         stats          = worker->stats;

   /* find Domain ranges */
   CLOCK_Start( timer );

   run_Decode_Domains( 
      q_seq, t_prof, Q, T, edg, sp_MX_fwd, sp_MX_bck, dom_def );
   stats->n_reported_domains += dom_def->dom_ranges->N;
   CLOCK_Stop( timer );
   times->sp_decodedom = CLOCK_Duration( timer );

   if ( args->verbose_level >= VERBOSE_HIGH )
   {
      fprintf(stdout, "# DOMAINS FOUND: %ld\n", dom_def->dom_ranges->N);
      for (int i = 0; i < dom_def->dom_ranges->N; i++) 
      {
         RANGE r_0 = VEC_X( dom_def->dom_ranges, i );
         fprintf(stdout, "[%d] {%d,%d}\n", i, r_0.beg, r_0.end);
      }
   }
}

/*! FUNCTION:  WORK_decode_posterior()
 *  SYNOPSIS:  Compute posterior from forward and backward matrices.
 */
void 
WORK_decode_posterior( WORKER* worker )
{
   FILE*                fp             = NULL;
   ARGS*                args           = worker->args;
   TASKS*               tasks          = worker->tasks;
   CLOCK*               timer          = worker->timer;
   /* input data */
   SEQUENCE*            q_seq          = worker->q_seq;
   int                  Q              = q_seq->N;
   HMM_PROFILE*         t_prof         = worker->t_prof;
   int                  T              = t_prof->N;
   EDGEBOUNDS*          edg            = worker->edg_row;
   /* working data */
   MATRIX_3D_SPARSE*    st_SMX         = worker->st_SMX;
   MATRIX_3D_SPARSE*    st_SMX_fwd     = worker->st_SMX_fwd;
   MATRIX_3D_SPARSE*    st_SMX_bck     = worker->st_SMX_bck;
   MATRIX_3D_SPARSE*    st_SMX_post    = worker->st_SMX_post;
   MATRIX_2D*           sp_MX          = worker->sp_MX;
   MATRIX_2D*           sp_MX_fwd      = worker->sp_MX_fwd;
   MATRIX_2D*           sp_MX_bck      = worker->sp_MX_bck;
   MATRIX_2D*           sp_MX_post     = worker->sp_MX_post;
   /* output data */
   TIMES*               times          = worker->times;
   RESULT*              result         = worker->result;
   ALL_SCORES*          scores         = &result->scores;
   SCORES*              final_scores   = &result->final_scores;
   float                sc;

   /* compute Posterior */
   CLOCK_Start( timer );
   run_Decode_Posterior_Sparse( 
      q_seq, t_prof, Q, T, edg, NULL,
      st_SMX_fwd, sp_MX_fwd, st_SMX_bck, sp_MX_bck, st_SMX_post, sp_MX_post );
   CLOCK_Stop( timer );
   times->sp_posterior = CLOCK_Duration( timer );

   fprintf(stdout, "# ==> Posterior (full cloud)\n");
   #if DEBUG 
   {
      fp = fopen("test_output/my.posterior.sp.log.mx", "w+");
      MATRIX_3D_SPARSE_Log_Embed(Q, T, st_SMX_post, debugger->test_MX);
      MATRIX_2D_Log(sp_MX_post);
      DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_post, fp);
      MATRIX_2D_Exp(sp_MX_post);
      fclose(fp);

      fp = fopen("test_output/my.posterior.sp.exp.mx", "w+");
      MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_post, debugger->test_MX);
      DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_post, fp);
      fclose(fp);
   }
   #endif
}

/*! FUNCTION:  WORK_decode_posterior()
 *  SYNOPSIS:  Compute the correction bias for the sequence.
 */
void 
WORK_null2_seq_bias( WORKER* worker )
{
   FILE*                fp             = NULL;
   ARGS*                args           = worker->args;
   TASKS*               tasks          = worker->tasks;
   CLOCK*               timer          = worker->timer;
   /* input data */
   SEQUENCE*            q_seq          = worker->q_seq;
   int                  Q              = q_seq->N;
   HMM_PROFILE*         t_prof         = worker->t_prof;
   int                  T              = t_prof->N;
   MATRIX_3D_SPARSE*    st_SMX         = worker->st_SMX;
   MATRIX_3D_SPARSE*    st_SMX_post    = worker->st_SMX_post;
   MATRIX_2D*           sp_MX          = worker->sp_MX;
   MATRIX_2D*           sp_MX_post     = worker->sp_MX_post;
   EDGEBOUNDS*          edg            = worker->edg_row;
   /* output data */
   TIMES*               times          = worker->times;
   RESULT*              result         = worker->result;
   ALL_SCORES*          scores         = &result->scores;
   SCORES*              finalsc        = &result->final_scores;
   DOMAIN_DEF*          dom_def        = worker->dom_def;
   float                null2_seq_bias;

   /* Composition Bias */
   CLOCK_Start( timer );
   run_Null2_ByExpectation_Sparse( 
      q_seq, t_prof, Q, T, edg, NULL, NULL, NULL,
      st_SMX_post, sp_MX_post, dom_def, &null2_seq_bias );
   CLOCK_Stop( timer );
   times->sp_biascorr      = CLOCK_Duration( timer );
   scores->null2_seq_bias           = null2_seq_bias;
   finalsc->null1_hmm_bias_natsc    = null2_seq_bias;
   finalsc->null1_hmm_bias_bitsc    = STATS_Nats_to_Bits( null2_seq_bias );

   if ( args->verbose_level >= VERBOSE_HIGH )
   {
      fprintf(stdout, "# ==> Null2 Compo Bias (full cloud): %11.4f %11.4f\n", 
         null2_seq_bias, null2_seq_bias/CONST_LOG2);
   }
}