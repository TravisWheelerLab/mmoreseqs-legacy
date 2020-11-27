/*******************************************************************************
 *  FILE:      posterior_quad.h
 *  PURPOSE:   The Maximum Posterior Probability and Optimal Alignment.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "../objects/structs.h"
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../algs_quad/algs_quad.h"
#include "../algs_linear/algs_linear.h"
#include "../parsers/parsers.h"

/* header */
#include "posterior_quad.h"


/*! FUNCTION:  run_Posterior_Quad()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward <st_MX_bck>.
 *             Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in supplied <st_MX_post> (can override input matrices).
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Posterior_Quad(  SEQUENCE*      q_seq,            /* query sequence */
                     HMM_PROFILE*   t_prof,           /* target hmm model */
                     int            Q,                /* query length */
                     int            T,                /* target length */
                     HMM_BG*        bg,               /* hmm background model */
                     EDGEBOUNDS*    edg,              /* edgebounds */
                     MATRIX_3D*     st_MX_fwd,        /* normal state matrix for forward */
                     MATRIX_2D*     sp_MX_fwd,        /* special state matrix for forward */
                     MATRIX_3D*     st_MX_bck,        /* normal state matrix for backward */
                     MATRIX_2D*     sp_MX_bck,        /* special state matrix for backward */
                     MATRIX_3D*     st_MX_post,       /* OUTPUT: normal state matrix for posterior */
                     MATRIX_2D*     sp_MX_post,       /* OUTPUT: special state matrix for posterior */     
                     DOMAIN_DEF*    dom_def )         /* OUTPUT: domain data */
{
   printf("=== POSTERIOR HEURISTICS ===\n");
   // printf("==> cutoffs: rt1=%6.3f, rt2=%6.3f, rt3=%6.3f\n",
   //    dom_def->rt1, dom_def->rt2, dom_def->rt3 );
   
   RANGE    Q_range, T_range;
   int      Q_size, T_size;
   float    sc, sc1, sc2;
   
   FILE*    fp;
   float    sc_fwd_full, sc_fwd_rng, sc_bck_full, sc_bck_rng;
   
   /* find the target and query range */
   Q_range.beg = Q + 1;
   Q_range.end = 0;
   T_range.beg = T + 1;
   T_range.end = 0;

   /* create bounding box */
   // EDGEBOUNDS_Dump(edg, stdout);
   for (int i = 0; i < edg->N; i++) {
      if ( T_range.beg > edg->bounds[i].lb ) {
         T_range.beg = edg->bounds[i].lb;
      }
      if ( T_range.end < edg->bounds[i].rb ) {
         T_range.end = edg->bounds[i].rb;
      }
   }
   Q_range.beg = edg->bounds[0].id;
   Q_range.end = edg->bounds[edg->N - 1].id;
   /* edge checks */
   T_range.beg = MAX(T_range.beg, 0);
   T_range.end = MIN(T_range.end, T);
   Q_range.beg = MAX(Q_range.beg, 0);
   Q_range.end = MIN(Q_range.end, Q);
   /* resize matrix to cover bounding box */
   Q_size = Q_range.end - Q_range.beg;
   T_size = T_range.end - T_range.beg;

   printf("Q: %d {%d,%d}, T: %d {%d,%d}\n", 
      Q, Q_range.beg, Q_range.end, T, T_range.beg, T_range.end );

   /* temporary override */
   // Q_range.beg = 0;
   // Q_range.end = Q;
   // T_range.beg = 0;
   // T_range.end = T;
   // printf_vhi("Q: %d {%d,%d}, T: %d {%d,%d}\n", 
   //    Q, Q_range.beg, Q_range.end, T, T_range.beg, T_range.end );

   
   /* resize special states */
   // MATRIX_2D_Reuse( sp_MX_fwd, NUM_SPECIAL_STATES, Q + 1 ); 
   // MATRIX_2D_Reuse( sp_MX_bck, NUM_SPECIAL_STATES, Q + 1 );
   // MATRIX_2D_Reuse( sp_MX_post, NUM_SPECIAL_STATES, Q + 1 );
   /* resize normal states */
   // MATRIX_3D_Reuse( st_MX_fwd, NUM_NORMAL_STATES, Q + 1, T + 1 );
   // MATRIX_3D_Reuse( st_MX_bck, NUM_NORMAL_STATES, Q + 1, T + 1 );
   // MATRIX_3D_Reuse( st_MX_post, NUM_NORMAL_STATES, Q + 1, T + 1 );
   /* clear previous data */
   // MATRIX_3D_Fill( st_MX_fwd, -INF );
   // MATRIX_3D_Fill( st_MX_bck, -INF );
   // MATRIX_3D_Fill( st_MX_post, -INF );

   /* run forward/backward on entire of profile/sequence */
   // fprintf(stdout, "# ==> Full Forward\n");
   // run_Forward_Quad(
   //    q_seq, t_prof, Q, T, st_MX_fwd, sp_MX_fwd, &sc_fwd_full );
   // fprintf(stdout, "# ==> Full Backward\n");
   // run_Backward_Quad(
   //    q_seq, t_prof, Q, T, st_MX_bck, sp_MX_bck, &sc_bck_full );
   
   // fp = fopen("fwd_full.mx", "w+");
   // DP_MATRIX_Dump(Q, T, st_MX_fwd, sp_MX_fwd, stdout);
   // fclose(fp);
   // fp = fopen("bck_full.mx", "w+");
   // DP_MATRIX_Dump(Q, T, st_MX_bck, sp_MX_bck, stdout);
   // fclose(fp);

   /* constrain sequence and profile */
   SEQUENCE_SetSubseq( q_seq, Q_range.beg, Q_range.end );
   HMM_PROFILE_SetSubmodel( t_prof, T_range.beg, T_range.end ); 
   HMM_PROFILE_ReconfigLength( t_prof, q_seq->N );
   
   /* resize special states */
   MATRIX_2D_Reuse( sp_MX_fwd, NUM_SPECIAL_STATES, Q_size + 1 ); 
   MATRIX_2D_Reuse( sp_MX_bck, NUM_SPECIAL_STATES, Q_size + 1 );
   MATRIX_2D_Reuse( sp_MX_post, NUM_SPECIAL_STATES, Q_size + 1 );
   /* resize normal states */
   MATRIX_3D_Reuse( st_MX_fwd, NUM_NORMAL_STATES, Q_size + 1, T_size + 1 );
   MATRIX_3D_Reuse( st_MX_bck, NUM_NORMAL_STATES, Q_size + 1, T_size + 1 );
   MATRIX_3D_Reuse( st_MX_post, NUM_NORMAL_STATES, Q_size + 1, T_size + 1 );
   /* clear previous data */
   // MATRIX_3D_Fill( st_MX_fwd, -INF );
   // MATRIX_3D_Fill( st_MX_bck, -INF );
   // MATRIX_3D_Fill( st_MX_post, -INF );

   /* run forward/backward on entire of profile/sequence */
   fprintf(stdout, "# ==> Ranged Forward\n");
   run_Forward_Quad(
      q_seq, t_prof, Q_size, T_size, st_MX_fwd, sp_MX_fwd, &sc_fwd_rng );
   fprintf(stdout, "# ==> Ranged Backward\n");
   run_Backward_Quad(
      q_seq, t_prof, Q_size, T_size, st_MX_bck, sp_MX_bck, &sc_bck_rng );

   /** TODO: fix Ranged functions */
   /* run ranged forward and backward */
   // fprintf(stdout, "# ==> Ranged Forward\n");
   // run_Ranged_Forward_Quad(
   //       q_seq, t_prof, &Q_range, &T_range, st_MX_fwd, sp_MX_fwd, &sc_fwd_rng );
   // fprintf(stdout, "# ==> Ranged Backward\n");
   // run_Ranged_Backward_Quad( 
   //       q_seq, t_prof, &Q_range, &T_range, st_MX_bck, sp_MX_bck, &sc_bck_rng );

   // fp = fopen("fwd_rng.mx", "w+");
   // DP_MATRIX_Dump(Q_size, T_size, st_MX_fwd, sp_MX_fwd, stdout);
   // fclose(fp);
   // fp = fopen("bck_rng.mx", "w+");
   // DP_MATRIX_Dump(Q_size, T_size, st_MX_bck, sp_MX_bck, stdout);
   // fclose(fp);

   fprintf(stdout, "# fwd score: (full) %7.3f, (ranged) %7.3f\n", sc_fwd_full, sc_fwd_rng);
   fprintf(stdout, "# bck score: (full) %7.3f, (ranged) %7.3f\n", sc_bck_full, sc_bck_rng);

   /* compute Posterior (forward * backward) */
   fprintf(stdout, "# ==> Posterior\n");
   run_Decode_Posterior_Quad( q_seq, t_prof, Q_size, T_size, &Q_range, &T_range,
         st_MX_fwd, sp_MX_fwd, st_MX_bck, sp_MX_bck, st_MX_post, sp_MX_post );
   
   #if DEBUG 
   {
      // FILE* fpout; 
      // fpout = fopen("test_output/my.posterior.mx", "w");
      // DP_MATRIX_Dump(Q, T, worker->st_MX_post, worker->sp_MX_post, fpout );
      // fclose(fpout);
      // fpout = fopen("test_output/my.posterior.log.mx", "w");
      // DP_MATRIX_Log_Dump(Q, T, worker->st_MX_post, worker->sp_MX_post, fpout );
      // fclose(fpout);
   }
   #endif
   
   /* run Null2 Score to compute Composition Bias */
   fprintf(stdout, "# ==> Null2 Compo Bias OLD\n");
   run_Null2_ByExpectation_Quad_Old( q_seq, t_prof, Q, T, 
         st_MX_post, sp_MX_post, dom_def );
   float seq_bias_old = dom_def->seq_bias;
   fprintf(stdout, "# ==> Null2 Compo Bias\n");
   run_Null2_ByExpectation_Quad( q_seq, t_prof, Q, T, &Q_range, &T_range, 
         st_MX_post, sp_MX_post, dom_def );
   float seq_bias_new = dom_def->seq_bias;
   fprintf(stdout, "seq_bias: %f %f\n", seq_bias_old, seq_bias_new);

   fprintf(stdout, "# ==> Posterior (end)\n");
   SEQUENCE_UnsetSubseq( q_seq );
   HMM_PROFILE_UnsetSubmodel( t_prof ); 
   HMM_PROFILE_ReconfigLength( t_prof, q_seq->N );

   return STATUS_SUCCESS;
}


/*! FUNCTION:  run_Decode_Normal_Posterior_Quad()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store matrix in <...fwd> or <...bck>.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Posterior_Quad( SEQUENCE*         q_seq,            /* query sequence */
                           HMM_PROFILE*      t_prof,           /* target hmm model */
                           int               Q,                /* full query length */
                           int               T,                /* full target length */
                           RANGE*            Q_range,          /* query range */
                           RANGE*            T_range,          /* target range */
                           MATRIX_3D*        st_MX_fwd,        /* normal state matrix for forward */
                           MATRIX_2D*        sp_MX_fwd,        /* special state matrix for forward */
                           MATRIX_3D*        st_MX_bck,        /* normal state matrix for backward */
                           MATRIX_2D*        sp_MX_bck,        /* special state matrix for backward */
                           MATRIX_3D*        st_MX_post,       /* OUTPUT: normal state matrix for posterior */
                           MATRIX_2D*        sp_MX_post )      /* OUTPUT: normal state matrix for posterior */
{
   printf("=== run_Decode_Posterior_Quad ===\n");
   // printf("==> FWD:\n");
   // DP_MATRIX_Dump(Q, T, st_MX_fwd, sp_MX_fwd, stdout );
   // printf("==> BCK:\n");
   // DP_MATRIX_Dump(Q, T, st_MX_bck, sp_MX_bck, stdout );

   /* query and target index */
   int      q_0, q_1;
   int      t_0, t_1;
   /* mapped query and target index */
   int      qx0, qx1;
   int      tx0, tx1;
   /* state index */
   int      st_0;
   /* overall score */
   float    overall_sc;
   /* common scale factor denominator */
   float    denom;
   /* temp mx scores */
   float    mmx, imx, dmx, smx;
   float    mmx_, imx_, dmx_, smx_;

   overall_sc  =  XMX_X( sp_MX_fwd, SP_C, Q ) + 
                  XSC_X(t_prof, SP_C, SP_MOVE);
   printf("overall_sc: %f %f ==> %f\n", 
      XMX_X( sp_MX_fwd, SP_C, Q ), XSC_X(t_prof, SP_C, SP_MOVE), overall_sc);
   // overall_sc = 0.0;

   /* init zero row */
   q_0 = 0;
   XMX_X(sp_MX_post, SP_E, q_0) = 0.0;
   XMX_X(sp_MX_post, SP_N, q_0) = 0.0;
   XMX_X(sp_MX_post, SP_J, q_0) = 0.0;
   XMX_X(sp_MX_post, SP_B, q_0) = 0.0;
   XMX_X(sp_MX_post, SP_C, q_0) = 0.0;

   for ( t_0 = 0; t_0 <= T; t_0++ )
   {
      MMX_X(st_MX_post, q_0, t_0) = 0.0;
      IMX_X(st_MX_post, q_0, t_0) = 0.0;
      DMX_X(st_MX_post, q_0, t_0) = 0.0; 
   }
   
   /* every position in query */
   for ( q_0 = 1; q_0 <= Q; q_0++ )
   {
      q_1   = q_0 - 1;
      denom = 0.0;
      t_0   = 0;

      MMX_X(st_MX_post, q_0, t_0) = 0.0;
      IMX_X(st_MX_post, q_0, t_0) = 0.0;
      DMX_X(st_MX_post, q_0, t_0) = 0.0; 

      /* every position in target */
      for ( t_0 = 1; t_0 < T; t_0++ )
      {
         t_1 = t_0 - 1;

         /* normal states */
         mmx = MMX_X(st_MX_fwd, q_0, t_0) + 
               MMX_X(st_MX_bck, q_0, t_0) - 
               overall_sc;

         mmx_ = expf(mmx);
         MMX_X(st_MX_post, q_0, t_0) = expf(mmx);
         denom += MMX_X(st_MX_post, q_0, t_0);

         imx = IMX_X(st_MX_fwd, q_0, t_0) + 
               IMX_X(st_MX_bck, q_0, t_0) - 
               overall_sc;
         imx_ = exp(imx);
         IMX_X(st_MX_post, q_0, t_0) = exp(imx);
         denom += IMX_X(st_MX_post, q_0, t_0);
         
         DMX_X(st_MX_post, q_0, t_0) = 0.0;
      }

      /* unrolled loop */
      t_0 = T;
      t_1 = T-1;

      /* normal states */
      mmx = MMX_X(st_MX_fwd, q_0, t_0) + 
            MMX_X(st_MX_bck, q_0, t_0) -
            overall_sc;
      MMX_X(st_MX_post, q_0, t_0) = exp(mmx);
      denom += MMX_X(st_MX_post, q_0, t_0);

      IMX_X(st_MX_post, q_0, t_0) = 0.0;
      DMX_X(st_MX_post, q_0, t_0) = 0.0;

      /* special states */
      XMX_X(sp_MX_post, SP_E, q_0) = 0.0;
      XMX_X(sp_MX_post, SP_B, q_0) = 0.0;

      smx = XMX_X(sp_MX_fwd, SP_N, q_1) +
            XMX_X(sp_MX_bck, SP_N, q_0) +
            XSC_X(t_prof, SP_N, SP_LOOP) -
            overall_sc;
      XMX_X(sp_MX_post, SP_N, q_0) = exp(smx);

      smx = XMX_X(sp_MX_fwd, SP_J, q_1) +
            XMX_X(sp_MX_bck, SP_J, q_0) +
            XSC_X(t_prof, SP_J, SP_LOOP) -
            overall_sc;
      XMX_X(sp_MX_post, SP_J, q_0) = exp(smx);

      smx = XMX_X(sp_MX_fwd, SP_C, q_1) +
            XMX_X(sp_MX_bck, SP_C, q_0) +
            XSC_X(t_prof, SP_C, SP_LOOP) -
            overall_sc;
      XMX_X(sp_MX_post, SP_C, q_0) = exp(smx);

      denom += XMX_X(sp_MX_post, SP_N, q_0) + 
               XMX_X(sp_MX_post, SP_J, q_0) + 
               XMX_X(sp_MX_post, SP_C, q_0);

      /* normalize by scaling row by common factor denominator */
      denom = 1.0 / denom;
      for ( t_0 = 1; t_0 < T; t_0++ ) {
         MMX_X(st_MX_post, q_0, t_0) *= denom;
         IMX_X(st_MX_post, q_0, t_0) *= denom;
      }
      MMX_X(st_MX_post, q_0, T) *= denom;
      XMX_X(sp_MX_post, SP_N, t_0) *= denom; 
      XMX_X(sp_MX_post, SP_J, t_0) *= denom; 
      XMX_X(sp_MX_post, SP_C, t_0) *= denom;
   }
   // printf("==> POSTERIOR:\n");
   // DP_MATRIX_Log_Dump(Q, T, st_MX_post, sp_MX_post, stdout );

   return STATUS_SUCCESS;
}


/*! FUNCTION:  run_Decode_Special_Posterior_Quad()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store matrix in <...fwd> or <...bck>.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Special_Posterior_Quad(  SEQUENCE*         q_seq,            /* query sequence */
                                    HMM_PROFILE*      t_prof,           /* target hmm model */
                                    int               Q,                /* query length */
                                    int               T,                /* target length */
                                    MATRIX_2D*        sp_MX_fwd,        /* special state matrix for forward */
                                    MATRIX_2D*        sp_MX_bck,        /* special state matrix for backward */
                                    MATRIX_2D*        sp_MX_post )      /* OUTPUT: special state matrix for posterior */
{
   printf("=== run_Decode_Special_Posterior_Quad ===\n");
   /* query index */
   int q_0, q_1;
   /* target index */ 
   int t_0, t_1;

   // printf("==> SPEC FWD:\n");
   // DP_MATRIX_SpecExp_Dump( Q, T, NULL, sp_MX_fwd, stdout );
   // printf("==> SPEC BCK:\n");
   // DP_MATRIX_SpecExp_Dump( Q, T, NULL, sp_MX_bck, stdout );

   /* init zero row */
   q_0 = 0;
   XMX_X(sp_MX_post, SP_E, q_0) = -INF;
   XMX_X(sp_MX_post, SP_N, q_0) = -INF;
   XMX_X(sp_MX_post, SP_J, q_0) = -INF;
   XMX_X(sp_MX_post, SP_B, q_0) = -INF;
   XMX_X(sp_MX_post, SP_C, q_0) = -INF;

   /* every position in query */
   for ( q_0 = 1; q_0 <= Q; q_0++ )
   {
      q_1 = q_0 - 1;

      /* update special states */
      XMX_X(sp_MX_post, SP_E, q_0) = -INF;
      XMX_X(sp_MX_post, SP_N, q_0) = logsum( XMX_X(sp_MX_fwd, SP_N, q_0), 
                                             XMX_X(sp_MX_bck, SP_N, q_0) );
      XMX_X(sp_MX_post, SP_J, q_0) = logsum( XMX_X(sp_MX_fwd, SP_J, q_0), 
                                             XMX_X(sp_MX_bck, SP_J, q_0) );
      XMX_X(sp_MX_post, SP_B, q_0) = -INF;
      XMX_X(sp_MX_post, SP_C, q_0) = logsum( XMX_X(sp_MX_fwd, SP_C, q_0), 
                                             XMX_X(sp_MX_bck, SP_C, q_0) );
   }

   // printf("=== SPEC POSTERIOR ===\n");
   // DP_MATRIX_SpecExp_Dump( Q, T, NULL, sp_MX_post, stdout );

   return STATUS_SUCCESS;
}


/*! FUNCTION:  run_Null2_By_Expectation_Quad()
 *  SYNOPSIS:  Modeled after HMMER p7_GNull2_ByExpectation().
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Null2_ByExpectation_Quad(    SEQUENCE*         q_seq,            /* query sequence */
                                 HMM_PROFILE*      t_prof,           /* target hmm model */
                                 int               Q,                /* query length */
                                 int               T,                /* target length */
                                 RANGE*            Q_range,          /* query range */
                                 RANGE*            T_range,          /* target range */
                                 MATRIX_3D*        st_MX_post,       /* posterior normal matrix */
                                 MATRIX_2D*        sp_MX_post,       /* posterior special matrix */
                                 DOMAIN_DEF*       dom_def )         /* OUTPUT: domain def's null2_sc vector */
{
   // printf("=== run_Null2_ByExpectation ===\n");
   int      Q_beg, Q_end, Q_len;
   int      T_beg, T_end, T_len;
   int      q_0, t_0; 
   int      st_0, k_0;
   float    x_factor;
   float    null2sc;

   Q_beg = MAX(0, Q_range->beg);
   Q_end = MIN(Q, Q_range->end);
   Q_len = Q_range->end - Q_range->beg;
   T_beg = MAX(0, T_range->beg);
   T_end = MIN(T, T_range->end);
   T_len = T_range->end - T_range->beg;

   // printf("=== POSTERIOR ===\n");
   // DP_MATRIX_Log_Dump(Q->end, T, st_MX_post, sp_MX_post, stdout);
   // printf("Q,T=(%d,%d)\n", q_seq->N, t_prof->N );
   // printf("=================\n");
   
   VECTOR_FLT_SetSizeTo( dom_def->st_freq, (T+1) * NUM_NORMAL_STATES );
   VECTOR_FLT_SetSizeTo( dom_def->sp_freq,  NUM_SPECIAL_STATES );
   VECTOR_FLT_SetSizeTo( dom_def->null2_sc, NUM_AMINO_PLUS_SPEC );

   // printf("st_freq->N = %d %d\n", dom_def->st_freq->Nalloc, dom_def->st_freq->N);
   // printf("sp_freq->N = %d %d\n", dom_def->sp_freq->Nalloc, dom_def->sp_freq->N);
   // printf("null2_sc->N = %d %d\n", dom_def->null2_sc->Nalloc, dom_def->null2_sc->Nalloc);

   q_0 = Q_beg;
   /* sum over each position in target model into vectors  */
   /* for each position in query domain */
   for ( t_0 = T_beg; t_0 <= T_end; t_0++ ) {
      /* for each normal state emissions */
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) = MX_3D( st_MX_post, st_0, q_0, t_0);
      }
   }
   /* for each special state emissions */
   for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
      VEC_X( dom_def->sp_freq, st_0 ) = MX_2D( sp_MX_post, st_0, q_0 );
   }

   // printf("<1>\n");
   // printf("==> NORMAL STATES <=\n");
   // for (int t_0 = 0; t_0 <= T; t_0++) {
   //    printf("[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + DEL_ST ) );
   // }
   // printf("==> SPECIAL STATES <=\n");
   // for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
   //    printf("[%d]: %12.9f\n", st_0, 
   //       VEC_X( dom_def->sp_freq, st_0 ) );
   // }

   /* sum over each position in query sequence domain into vectors  */
   for ( q_0 = Q_beg + 1; q_0 <= Q_end; q_0++ ) {
      /* for each position in target model */
      for ( t_0 = T_beg; t_0 <= T_end; t_0++ ) {
         /* for each normal state emissions */
         for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
            VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) += MX_3D( st_MX_post, st_0, q_0, t_0 );
         }
      }
      /* for each special state emissions */
      for ( st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
         VEC_X( dom_def->sp_freq, st_0 ) += MX_2D( sp_MX_post, st_0, q_0 );
      }
   }

   // printf("<2>\n");
   // printf("==> NORMAL STATES <=\n");
   // for (int t_0 = 0; t_0 <= T; t_0++) {
   //    printf("[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + DEL_ST ) );
   // }
   // printf("==> SPECIAL STATES <=\n");
   // for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
   //    printf("[%d]: %12.9f\n", st_0, 
   //       VEC_X( dom_def->sp_freq, st_0 ) );
   // }

   /* convert expected numbers to log frequencies */
   /* for each position in query domain */
   for ( t_0 = T_beg; t_0 <= T_end; t_0++ ) {
      /* for each normal state emissions */
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) = log( VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) );
      }
   }
   /* for each special state emissions */
   for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
      VEC_X( dom_def->sp_freq, st_0 ) = log( VEC_X( dom_def->sp_freq, st_0 ) );
   }

   // printf("<3>\n");
   // printf("==> NORMAL STATES <=\n");
   // for (int t_0 = 0; t_0 <= T; t_0++) {
   //    printf("[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + DEL_ST ) );
   // }
   // printf("==> SPECIAL STATES <=\n");
   // for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
   //    printf("[%d]: %12.9f\n", st_0, 
   //       VEC_X( dom_def->sp_freq, st_0 ) );
   // }

   float neglog_Q = -log( (float)Q_len );
   // printf("neglog_Q = %d %f\n", q_len, neglog_Q);
   /* for each position in query domain */
   for ( t_0 = T_beg; t_0 <= T_end; t_0++ ) {
      /* for each normal state emissions */
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) += neglog_Q;
      }
   }
   /* for each special state emissions */
   for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
      VEC_X( dom_def->sp_freq, st_0 ) += neglog_Q;
   }

   // printf("<4>\n");
   // printf("==> NORMAL STATES <=\n");
   // for (int t_0 = 0; t_0 <= T; t_0++) {
   //    printf("[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + DEL_ST ) );
   // }
   // printf("==> SPECIAL STATES <=\n");
   // for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
   //    printf("[%d]: %12.9f\n", st_0, 
   //       VEC_X( dom_def->sp_freq, st_0 ) );
   // }

   /* x-factor: */
   x_factor = VEC_X( dom_def->sp_freq, SP_N);
   x_factor = logsum( x_factor,
                      VEC_X(dom_def->sp_freq, SP_C) );
   x_factor = logsum( x_factor,
                      VEC_X(dom_def->sp_freq, SP_J) );
   
   // printf("<5>\n");
   // printf("xfactor: %f, NCJ: %f %f %f\n", 
   //    x_factor, VEC_X( dom_def->sp_freq, SP_N), VEC_X(dom_def->sp_freq, SP_C), VEC_X(dom_def->sp_freq, SP_J) );

   /* initialize null2 vector */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = -INF;
   }
   for ( k_0 = NUM_AMINO; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = 0.0;
   }

   // printf("<6>\n");
   /* null2 log emissions probabilities found by summing over 
    * all emmissions used in paths explaining the domain. 
    */
   /* for each amino acid */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) 
   {
      // printf("k_0 = %d\n", k_0);
      /* for each position in model */
      for ( t_0 = T_beg + 1; t_0 < T; t_0++ ) 
      {
         // printf("t_0 = %d\n", t_0);
         VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                   VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ) + MSC_X( t_prof, t_0, k_0 ) );
         VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                   VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ) + ISC_X( t_prof, t_0, k_0 ) );

         // if ( t_0 == 14 ) {
         //    printf("(k_0,t_0)=%2d,%2d, MSC=%f, ISC=%f, NULL2=%f\n",
         //       k_0, t_0, 
         //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ) + MSC_X( t_prof, t_0, k_0 ),
         //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ) + ISC_X( t_prof, t_0, k_0 ),
         //       VEC_X( dom_def->null2_sc, k_0 )
         //    );
         // }
      }
      t_0 = T;
      VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ) + MSC_X( t_prof, t_0, k_0 ) );
      VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                x_factor );                            
   }

   // printf("==> NULL2 (pre) <=\n");
   // for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
   //    printf("[%2d] %f\n", k_0, VEC_X( dom_def->null2_sc, k_0) );
   // }

   /* convert log to normal scores */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = exp( VEC_X( dom_def->null2_sc, k_0 ) );
   }

   // printf("==> NULL2 (norm) <=\n");
   // for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
   //    printf("[%2d] %f\n", k_0, dom_def->null2_sc->data[k_0] );
   // }

   /* average the odds ratios */
   for ( k_0 = NUM_AMINO+1; k_0 < NUM_AMINO_PLUS_SPEC-3; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = exp( VEC_X( dom_def->null2_sc, k_0 ) );
   }

   // printf("==> NULL2 (avg) <=\n");
   // for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
   //    printf("[%2d] %f\n", k_0, dom_def->null2_sc->data[k_0] );
   // }

   /* set scores for all degenerencies */
   VEC_X( dom_def->null2_sc, NUM_AMINO_PLUS_SPEC - 1 ) = 1.0;
   for ( k_0 = NUM_AMINO; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = 1.0;
   }

   // printf("==> NULL2 (end) <=\n");
   // for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
   //    printf("[%2d] %f\n", k_0, dom_def->null2_sc->data[k_0] );
   // }

   /* sum over all nullscore vector to get seqbias */
   dom_def->seq_bias = 0.0;
   for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      dom_def->seq_bias += VEC_X( dom_def->null2_sc, k_0 );
   }
   // printf("seqbias: %f\n", dom_def->seq_bias);

   /* convert nullscore to logscores */
   for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = log( VEC_X( dom_def->null2_sc, k_0 ) );
   }

   /* TODO: Do I need to go in and out of logspace? */
   // printf("==> NULL2 (logsc) <=\n");
   // for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
   //    printf("[%2d] %f\n", k_0, dom_def->null2_sc->data[k_0] );
   // }

   /* multiply (add in logspace) all biases by symbol in alphabet */
   dom_def->seq_bias = VEC_X( dom_def->null2_sc, 0 );
   for ( k_0 = 1; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      dom_def->seq_bias += VEC_X( dom_def->null2_sc, k_0 );
   }
   // printf("seqbias: %f\n", dom_def->seq_bias);

   return STATUS_SUCCESS;
}

/*! FUNCTION:  run_Null2_ByExpectation_Old()
 *  SYNOPSIS:  Modeled after HMMER p7_GNull2_ByExpectation().
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Null2_ByExpectation_Quad_Old(   SEQUENCE*            q_seq,            /* query sequence */
                                    HMM_PROFILE*         t_prof,           /* target hmm model */
                                    int                  Q,                /* query length */
                                    int                  T,                /* target length */
                                    MATRIX_3D*           st_MX_post,       /* posterior normal matrix */
                                    MATRIX_2D*           sp_MX_post,       /* posterior special matrix */
                                    DOMAIN_DEF*          dom_def )         /* OUTPUT: domain def's null2_sc vector */
{
   // printf("=== run_Null2_ByExpectation_OLD ===\n");
   int      q_beg, q_end, q_len;
   int      t_beg, t_end, t_len;
   int      q_0, t_0; 
   int      st_0, k_0;
   float    x_factor;
   float    null2sc;

   q_beg = 0;
   q_end = Q;
   q_len = q_end - q_beg;
   t_beg = 0;
   t_end = T;
   t_len = t_end - t_beg;

   // printf("=== POSTERIOR ===\n");
   // // DP_MATRIX_Dump(Q, T, st_MX_post, sp_MX_post, stdout);
   // DP_MATRIX_Log_Dump(Q, T, st_MX_post, sp_MX_post, stdout);
   // printf("Q,T=(%d,%d)\n", q_seq->N, t_prof->N );
   // printf("=================\n");
   
   VECTOR_FLT_GrowTo( dom_def->st_freq, (T+1) * NUM_NORMAL_STATES );
   VECTOR_FLT_GrowTo( dom_def->sp_freq, NUM_SPECIAL_STATES );
   VECTOR_FLT_GrowTo( dom_def->null2_sc, NUM_AMINO_PLUS_SPEC );

   q_0 = q_beg;
   /* sum over each position in target model into vectors  */
   /* for each position in query domain */
   for ( t_0 = 0; t_0 <= T; t_0++ ) {
      /* for each normal state emissions */
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) = MX_3D( st_MX_post, st_0, q_0, t_0);
      }
   }
   /* for each special state emissions */
   for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
      VEC_X( dom_def->sp_freq, st_0 ) = MX_2D( sp_MX_post, st_0, q_0 );
   }

   // printf("<0>\n");
   // printf("==> NORMAL STATES <=\n");
   // for (int t_0 = 0; t_0 <= T; t_0++) {
   //    printf("[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + DEL_ST ) );
   // }
   // printf("==> SPECIAL STATES <=\n");
   // for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
   //    printf("[%d]: %12.9f\n", st_0, 
   //       VEC_X( dom_def->sp_freq, st_0 ) );
   // }

   /* sum over each position in query sequence domain into vectors  */
   for ( q_0 = q_beg + 1; q_0 <= q_end; q_0++ ) {
      /* for each position in target model */
      for ( t_0 = 0; t_0 <= T; t_0++ ) {
         /* for each normal state emissions */
         for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
            VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) += MX_3D( st_MX_post, st_0, q_0, t_0 );
         }
      }
      /* for each special state emissions */
      for ( st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
         VEC_X( dom_def->sp_freq, st_0 ) += MX_2D( sp_MX_post, st_0, q_0 );
      }
   }

   // printf("<1>\n");
   // printf("==> NORMAL STATES <=\n");
   // for (int t_0 = 0; t_0 <= T; t_0++) {
   //    printf("[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + DEL_ST ) );
   // }
   // printf("==> SPECIAL STATES <=\n");
   // for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
   //    printf("[%d]: %12.9f\n", st_0, 
   //       VEC_X( dom_def->sp_freq, st_0 ) );
   // }

   /* convert expected numbers to log frequencies */
   /* for each position in query domain */
   for ( t_0 = 0; t_0 <= T; t_0++ ) {
      /* for each normal state emissions */
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) = log( VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) );
      }
   }
   /* for each special state emissions */
   for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
      VEC_X( dom_def->sp_freq, st_0 ) = log( VEC_X( dom_def->sp_freq, st_0 ) );
   }

   // printf("<2>\n");
   // printf("==> NORMAL STATES <=\n");
   // for (int t_0 = 0; t_0 <= T; t_0++) {
   //    printf("[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + DEL_ST ) );
   // }
   // printf("==> SPECIAL STATES <=\n");
   // for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
   //    printf("[%d]: %12.9f\n", st_0, 
   //       VEC_X( dom_def->sp_freq, st_0 ) );
   // }

   float neglog_Q = -log( (float)q_len );
   // printf("neglog_Q = %d %f\n", q_len, neglog_Q);
   /* for each position in query domain */
   for ( t_0 = 0; t_0 <= T; t_0++ ) {
      /* for each normal state emissions */
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) += neglog_Q;
      }
   }
   /* for each special state emissions */
   for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
      VEC_X( dom_def->sp_freq, st_0 ) += neglog_Q;
   }

   // printf("<3>\n");
   // printf("==> NORMAL STATES <=\n");
   // for (int t_0 = 0; t_0 <= T; t_0++) {
   //    printf("[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ), 
   //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + DEL_ST ) );
   // }
   // printf("==> SPECIAL STATES <=\n");
   // for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
   //    printf("[%d]: %12.9f\n", st_0, 
   //       VEC_X( dom_def->sp_freq, st_0 ) );
   // }

   /* x-factor: */
   x_factor = VEC_X( dom_def->sp_freq, SP_N);
   x_factor = logsum( x_factor,
                      VEC_X(dom_def->sp_freq, SP_C) );
   x_factor = logsum( x_factor,
                      VEC_X(dom_def->sp_freq, SP_J) );

   // printf("<4>\n");
   // printf("xfactor: %f, NCJ: %f %f %f\n", 
   //    x_factor, VEC_X( dom_def->sp_freq, SP_N), VEC_X(dom_def->sp_freq, SP_C), VEC_X(dom_def->sp_freq, SP_J) );

   /* initialize null2 vector */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = -INF;
   }
   for ( k_0 = NUM_AMINO; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = 0.0;
   }

   /* null2 log emissions probabilities found by summing over 
    * all emmissions used in paths explaining the domain. 
    */
   /* for each amino acid */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) {
      /* for each position in model */
      for ( t_0 = 1; t_0 < T; t_0++ ) {
         VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                   VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ) + MSC_X( t_prof, t_0, k_0 ) );
         VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                   VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ) + ISC_X( t_prof, t_0, k_0 ) );

         // if ( t_0 == 14 ) {
         //    printf("(k_0,t_0)=%2d,%2d, MSC=%f, ISC=%f, NULL2=%f\n",
         //       k_0, t_0, 
         //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ) + MSC_X( t_prof, t_0, k_0 ),
         //       VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST ) + ISC_X( t_prof, t_0, k_0 ),
         //       VEC_X( dom_def->null2_sc, k_0 )
         //    );
         // }
      }
      t_0 = T;
      VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST ) + MSC_X( t_prof, t_0, k_0 ) );
      VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                x_factor );                            
   }

   // printf("==> NULL2 (pre) <=\n");
   // for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
   //    printf("[%2d] %f\n", k_0, VEC_X( dom_def->null2_sc, k_0) );
   // }

   /* convert log to normal scores */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = exp( VEC_X( dom_def->null2_sc, k_0 ) );
   }

   // printf("==> NULL2 (norm) <=\n");
   // for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
   //    printf("[%2d] %f\n", k_0, dom_def->null2_sc->data[k_0] );
   // }

   /* average the odds ratios */
   for ( k_0 = NUM_AMINO+1; k_0 < NUM_AMINO_PLUS_SPEC-3; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = exp( VEC_X( dom_def->null2_sc, k_0 ) );
   }

   // printf("==> NULL2 (avg) <=\n");
   // for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
   //    printf("[%2d] %f\n", k_0, dom_def->null2_sc->data[k_0] );
   // }

   /* set scores for all degenerencies */
   VEC_X( dom_def->null2_sc, NUM_AMINO_PLUS_SPEC - 1 ) = 1.0;
   for ( k_0 = NUM_AMINO; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = 1.0;
   }

   // printf("==> NULL2 (end) <=\n");
   // for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
   //    printf("[%2d] %f\n", k_0, dom_def->null2_sc->data[k_0] );
   // }

   /* sum over all nullscore vector to get seqbias */
   dom_def->seq_bias = 0.0;
   for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      dom_def->seq_bias += VEC_X( dom_def->null2_sc, k_0 );
   }
   // printf("seqbias: %f\n", dom_def->seq_bias);

   /* convert nullscore to logscores */
   for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = log( VEC_X( dom_def->null2_sc, k_0 ) );
   }

   /* TODO: Do I need to go in and out of logspace? */
   // printf("==> NULL2 (logsc) <=\n");
   // for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
   //    printf("[%2d] %f\n", k_0, dom_def->null2_sc->data[k_0] );
   // }

   /* multiply (add in logspace) all biases by symbol in alphabet */
   dom_def->seq_bias = VEC_X( dom_def->null2_sc, 0 );
   for ( k_0 = 1; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      dom_def->seq_bias += VEC_X( dom_def->null2_sc, k_0 );
   }
   // printf("seqbias: %f\n", dom_def->seq_bias);

   return STATUS_SUCCESS;
}