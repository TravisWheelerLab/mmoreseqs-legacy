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
#include "work_posterior_bydom.h"

/*! FUNCTION:  WORK_posterior()
 *  SYNOPSIS:  Run full posterior, for full sequence and domain-specific. 
 */
void 
WORK_posterior_bydom( WORKER* worker )
{
   WORK_posterior_sparse_bydom( worker );

   // /* compute hmm model bias */
   // WORK_null1_hmm_bias( worker );
   // /* build sparse matrix */
   // WORK_build_sparse_matrix( worker );
   // /* compute fwdback */
   // WORK_bound_fwdback_sparse( worker );
   // /* decode domains */
   // WORK_decode_domains( worker );
   // /* compute posterior */
   // WORK_decode_posterior( worker );
   // /* compute sequence bias */
   // WORK_null2_seq_bias( worker );
   // /* compute optimal accuracy from posterior */
   // // WORK_optimal_accuracy( worker );
   // /* backtrace optimal accuracy for posterior alignment */
   // // WORK_optacc_traceback( worker );
}

/*! FUNCTION:  WORK_posterior_sparse_bydom()
 *  SYNOPSIS:  Run posterior using sparse matrices, per domain.
 */
void 
WORK_posterior_sparse_bydom( WORKER* worker )
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
   DOMAIN_DEF*          dom_def        = worker->dom_def;
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
   /* loop vars */
   RANGE                D_range;
   int                  D_size;
   float                null1_hmm_bias, null2_seq_bias;
   float                pre_sc, fwd_sc, bck_sc;
   float                post_sc, opt_sc, dom_sc;

   /* run through Domains and compute score, bias correction, and optimal alignment */
   if ( args->is_run_domains == true )
   {
      dom_def->n_domains      = dom_def->dom_ranges->N;
      dom_def->dom_sumsc      = 0.0f;
      dom_def->dom_sumbias    = 0.0f;
      dom_def->n_residues     = 0;
      null1_hmm_bias          = dom_def->null1_hmm_bias;
      times->dom_start        = CLOCK_GetTime( timer );

      /* TODO: fix null1 bias computation */
      null1_hmm_bias = 0.0f;

      for (int i = 0; i < dom_def->n_domains; i++)
      {
         dom_def->idx = i;
         D_range = VEC_X( dom_def->dom_ranges, i );

         printf("Domain (%d of %d): {%d,%d}\n", 
            i+1, dom_def->n_domains, D_range.beg, D_range.end);

         // /* Reparameterize sequence to only cover domain range. */
         D_size = D_range.end - D_range.beg + 1;
         // SEQUENCE_SetDomain( q_seq, D_range );
         // EDGEBOUNDS_SetDomain( edg, edg_dom, D_range );
         // HMM_PROFILE_ReconfigLength( t_prof, q_seq->N );
         
         /* TODO: Should be able to eliminate this */
         /* clear previous data */
         MATRIX_3D_SPARSE_Fill( st_SMX_fwd, -INF );
         MATRIX_3D_SPARSE_Fill( st_SMX_bck, -INF );
         if (st_SMX_bck != st_SMX_post) {
            MATRIX_3D_SPARSE_Fill( st_SMX_post, -INF );
         }

         /* compute Forward/Backward for the domain range */
         CLOCK_Start( timer );
         run_Bound_Forward_Sparse( 
            q_seq, t_prof, Q, T, st_SMX_fwd, sp_MX_fwd, edg, &D_range, &fwd_sc );
         CLOCK_Stop( timer );
         times->dom_bound_fwd += CLOCK_Duration( timer );

         if ( args->verbose_level >= VERBOSE_HIGH )
         {
            fprintf(stdout, "# ==> Forward   (domain %d/%d): %11.4f %11.4f\n", 
               i+1, dom_def->n_domains, fwd_sc, fwd_sc/CONST_LOG2);
         }
         #if DEBUG 
         {
            fp = fopen("test_output/my.sparse_fwd.dom.mx", "w+");
            MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_fwd, debugger->test_MX);
            DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_fwd, fp);
            fclose(fp);
         }
         #endif

         /* compute Forward/Backward for the domain range */
         CLOCK_Start( timer );
         run_Bound_Backward_Sparse( 
            q_seq, t_prof, Q, T, st_SMX_bck, sp_MX_bck, edg, &D_range, &bck_sc );
         CLOCK_Stop( timer );
         times->dom_bound_bck += CLOCK_Duration( timer );
         
         if ( args->verbose_level >= VERBOSE_HIGH )
         {
            fprintf(stdout, "# ==> Backward  (domain %d/%d): %11.4f %11.4f\n", 
               i+1, dom_def->n_domains, bck_sc, bck_sc/CONST_LOG2);
         }
         #if DEBUG 
         {
            fp = fopen("test_output/my.sparse_bck.dom.mx", "w+");
            MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_bck, debugger->test_MX);
            DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_bck, fp);
            fclose(fp);
         }
         #endif

         /* compute Posterior (forward * backward) for domain range */
         CLOCK_Start( timer );
         run_Decode_Posterior_Sparse( 
            q_seq, t_prof, Q, T, edg, &D_range,
            st_SMX_fwd, sp_MX_fwd, st_SMX_bck, sp_MX_bck, st_SMX_post, sp_MX_post );
         CLOCK_Stop( timer );
         times->dom_posterior += CLOCK_Duration( timer );

         if ( args->verbose_level >= VERBOSE_HIGH )
         {
            fprintf(stdout, "# ==> Posterior (domain %d/%d)\n", 
               i+1, dom_def->n_domains );
         }
         #if DEBUG 
         {
            fp = fopen("test_output/my.sparse_post.dom.mx", "w+");
            MATRIX_3D_SPARSE_Log_Embed(Q, T, st_SMX_post, debugger->test_MX);
            MATRIX_2D_Log(sp_MX_post);
            DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_post, fp);
            MATRIX_2D_Exp(sp_MX_post);
            fclose(fp);
         }
         #endif
         
         /* run Null2 Score to compute Composition Bias */
         CLOCK_Start( timer );
         run_Null2_ByExpectation_Sparse( q_seq, t_prof, Q, T, edg, &D_range, NULL, NULL,
            st_SMX_post, sp_MX_post, dom_def, &null2_seq_bias );
         CLOCK_Stop( timer );
         times->dom_biascorr += CLOCK_Duration( timer );

          if ( args->verbose_level >= VERBOSE_HIGH )
         {
            fprintf(stdout, "# ==> Null2 Compo Bias  (domain %d/%d): %11.4f %11.4f\n", 
               i+1, dom_def->n_domains, null2_seq_bias, null2_seq_bias/CONST_LOG2);
         }

         /** TODO: Get optimal alignment */ 
         
         
         /* add domain data */
         VECTOR_FLT_Pushback( dom_def->dom_fwdsc, fwd_sc );
         VECTOR_FLT_Pushback( dom_def->dom_bias, null2_seq_bias );

         /* check if best score */
         pre_sc = (fwd_sc - (null1_hmm_bias)) / CONST_LOG2;
         dom_sc = (fwd_sc - (null1_hmm_bias + null2_seq_bias)) / CONST_LOG2;
         if ( dom_sc > dom_def->best_sc )
         {
            dom_def->best        = i;
            dom_def->best_sc     = dom_sc;
            dom_def->best_fwdsc  = fwd_sc;
            dom_def->best_presc  = pre_sc;
            dom_def->best_bias   = null2_seq_bias;
            dom_def->best_range  = D_range;
         }

         /* constructed score over all domains */
         dom_def->dom_sumsc   += fwd_sc;
         dom_def->dom_sumbias += null2_seq_bias;
         dom_def->n_residues  += (D_range.end - D_range.beg + 1);
      }

      /* final reconstructed score */
      if ( dom_def->n_domains > 0 )
      {
         /* constructed score over all domains */
         dom_def->dom_sumbias = MATH_Sum(0.0f, log(dom_def->null_omega) + dom_def->dom_sumbias);
         dom_def->dom_sumsc  += (Q - dom_def->n_residues) * log((float) Q / (float) (Q + 3));
         dom_def->dom_sumsc   = (dom_def->dom_sumsc - (dom_def->nullsc + dom_def->dom_sumbias)) / CONST_LOG2;
      }
      else {
         dom_def->dom_sumbias = -INF;
         dom_def->dom_sumsc   = -INF;
      }

      /* total dom-scoring time */
      times->dom_end    = CLOCK_GetTime( timer );
      times->dom_total  = CLOCK_GetDiff( timer, times->dom_start, times->dom_end );
   }
}

/*! FUNCTION:  	WORK_bound_fwdback_sparse_bydom()
 *  SYNOPSIS:  	Run "bound forward/backward" step of pruned forward/backward (aka adaptive-band forward/backward).
 *                Uses the sparse matrix implementation.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_bound_fwdback_sparse_bydom( WORKER* worker )
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
   DOMAIN_DEF*          dom_def        = worker->dom_def;
   EDGEBOUNDS*          edg_dom        = dom_def->edg;
   /* working data */
   MATRIX_3D_SPARSE*    st_SMX         = worker->st_SMX;
   MATRIX_3D_SPARSE*    st_SMX_fwd     = worker->st_SMX_fwd;
   MATRIX_3D_SPARSE*    st_SMX_bck     = worker->st_SMX_bck;
   MATRIX_2D*           sp_MX          = worker->sp_MX;
   MATRIX_2D*           sp_MX_fwd      = worker->sp_MX_fwd;
   MATRIX_2D*           sp_MX_bck      = worker->sp_MX_bck;
   /* output data */
   TIMES*               times          = worker->times;
   RESULT*              result         = worker->result;
   ALL_SCORES*          scores         = &result->scores;
   SCORES*              final_scores   = &result->final_scores; 
   /* domain data */
   DOMAIN_X*            domain         = &dom_def->domain;
   int                  idx            = dom_def->idx;
   RANGE                D_range        = VEC_X( dom_def->dom_ranges, idx );
   float                fwd_sc, bck_sc;

   /* compute Forward/Backward for the domain range */
   CLOCK_Start( timer );
   run_Bound_Forward_Sparse( 
      q_seq, t_prof, Q, T, st_SMX_fwd, sp_MX_fwd, edg_dom, &D_range, &fwd_sc );
   CLOCK_Stop( timer );
   times->dom_bound_fwd += CLOCK_Duration( timer );

   if ( args->verbose_level >= VERBOSE_HIGH )
   {
      fprintf(stdout, "# ==> Forward   (domain %d/%d): %11.4f %11.4f\n", 
         idx + 1, dom_def->n_domains, fwd_sc, fwd_sc/CONST_LOG2);
   }
   #if DEBUG 
   {
      fp = fopen("test_output/my.sparse_fwd.dom.mx", "w+");
      MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_fwd, debugger->test_MX);
      DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_fwd, fp);
      fclose(fp);
   }
   #endif


}

/*! FUNCTION:  WORK_null1_hmm_bias_bydom()
 *  SYNOPSIS:  Compute the correction bias for the hmm model.
 *             Only for domain region.
 */
void 
WORK_null1_hmm_bias_bydom( WORKER* worker )
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
   SCORES*        final_scores   = &result->final_scores;
   
   /* initialize hmm_bg */
   HMM_BG_SetSequence( hmm_bg, q_seq );
   
   HMM_BG_SetFilter( hmm_bg, t_prof->N, t_prof->bg_model->compo );
   HMM_BG_SetLength( hmm_bg, q_seq->N );
   /* compute null one */
   HMM_BG_NullOne( hmm_bg, q_seq->N, &scores->null1_hmm_bias );
   /* compute nullscore for bias */
   HMM_BG_FilterScore( hmm_bg, q_seq, &scores->null1_filtersc );
   /* fetch omega (prior prob of no bias) */
   scores->null_omega = hmm_bg->omega;
   /* free digitized sequence TODO: move to sequence */
   // HMM_BG_UnsetSequence( hmm_bg, q_seq );
}

/*! FUNCTION:  WORK_decode_posterior_bydom()
 *  SYNOPSIS:  Compute posterior from forward and backward matrices.
 *             Only for domain region.
 */
void 
WORK_decode_posterior_bydom( WORKER* worker )
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
      fp = fopen("test_output/my.sparse_post.mx", "w+");
      MATRIX_3D_SPARSE_Log_Embed(Q, T, st_SMX_post, debugger->test_MX);
      MATRIX_2D_Log(sp_MX_post);
      DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_post, fp);
      MATRIX_2D_Exp(sp_MX_post);
      fclose(fp);
   }
   #endif
}

/*! FUNCTION:  WORK_decode_posterior_bydom()
 *  SYNOPSIS:  Compute the correction bias for the sequence.
 *             Only for domain region.
 */
void 
WORK_null2_seq_bias_bydom( WORKER* worker )
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
   SCORES*              final_scores   = &result->final_scores;
   DOMAIN_DEF*          dom_def        = worker->dom_def;
   float                null2_seq_bias;

   /* Composition Bias */
   CLOCK_Start( timer );
   run_Null2_ByExpectation_Sparse( 
      q_seq, t_prof, Q, T, edg, NULL, NULL, NULL,
      st_SMX_post, sp_MX_post, dom_def, &null2_seq_bias );
   CLOCK_Stop( timer );
   times->sp_biascorr      = CLOCK_Duration( timer );
   scores->null2_seq_bias  = null2_seq_bias;

   fprintf(stdout, "# ==> Null2 Compo Bias (full cloud): %11.4f %11.4f\n", null2_seq_bias, null2_seq_bias/CONST_LOG2);
}