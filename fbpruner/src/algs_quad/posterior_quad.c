/*******************************************************************************
 *  FILE:      posterior_quad.h
 *  PURPOSE:   The Maximum Posterior Probability and Optimal Alignment.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *    - None known.
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
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../algs_quad/_algs_quad.h"
#include "../algs_linear/_algs_linear.h"
#include "../parsers/_parsers.h"

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
   printf("=== run_Posterior_Quad() ===\n");
   // printf("==> cutoffs: rt1=%6.3f, rt2=%6.3f, rt3=%6.3f\n",
   //    dom_def->rt1, dom_def->rt2, dom_def->rt3 );
   
   RANGE    Q_range, T_range;
   int      Q_size, T_size;
   float    sc, sc1, sc2;
   
   FILE*    fp;
   float    sc_fwd_full, sc_fwd_rng, sc_fwd_con;
   float    sc_bck_full, sc_bck_rng, sc_bck_con;
   float    seq_bias;
   
   /* find the target and query range */
   Q_range.beg = Q;
   Q_range.end = 0;
   T_range.beg = T;
   T_range.end = 0;

   /* create bounding box */
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

   printf("BOUNDING BOX => Q: %d {%d,%d} -> %d, T: %d {%d,%d} -> %d\n", 
         Q, Q_range.beg, Q_range.end, Q_size, T, T_range.beg, T_range.end, T_size );
   
   /* TESTING => run full forward/backward */
   #if DEBUG 
   {
      /* temporary override to complete sequence */
      // Q_range.beg = 0;
      // Q_range.end = Q+1;
      // T_range.beg = 0;
      // T_range.end = T+1;
      // printf("BOUNDING BOX => Q: %d {%d,%d}, T: %d {%d,%d}\n", 
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
   }
   #endif

   /* TESTING => override to explicit values */
   #if DEBUG
   {
      // Q_range.beg = 16;
      // Q_range.end = 816;
      // T_range.beg = 3;
      // T_range.end = 791;
      Q_range.beg = 0;
      Q_range.end = Q;
      T_range.beg = 0;
      T_range.end = T;
      /* resize matrix to cover bounding box */
      Q_size = Q_range.end - Q_range.beg;
      T_size = T_range.end - T_range.beg;
      printf("BOUNDING BOX => Q: %d {%d,%d} -> %d, T: %d {%d,%d} -> %d\n", 
         Q, Q_range.beg, Q_range.end, Q_size, T, T_range.beg, T_range.end, T_size );
   }
   #endif

   /* choose posterior calculation method => 0 for all */
   #define POST_METHOD 1

   /* METHOD 1: For Testing Only. Use full forward/backward, for checking Posterior computation */
   #if POST_METHOD == 0 || POST_METHOD == 1
   {
      /* resize special states */
      MATRIX_2D_Reuse( sp_MX_fwd, NUM_SPECIAL_STATES, Q + 1 ); 
      MATRIX_2D_Reuse( sp_MX_bck, NUM_SPECIAL_STATES, Q + 1 );
      MATRIX_2D_Reuse( sp_MX_post, NUM_SPECIAL_STATES, Q + 1 );
      /* resize normal states */
      MATRIX_3D_Reuse( st_MX_fwd, NUM_NORMAL_STATES, Q + 1, T + 1 );
      MATRIX_3D_Reuse( st_MX_bck, NUM_NORMAL_STATES, Q + 1, T + 1 );
      MATRIX_3D_Reuse( st_MX_post, NUM_NORMAL_STATES, Q + 1, T + 1 );
      /* clear previous data */
      MATRIX_3D_Fill( st_MX_fwd, -INF );
      MATRIX_3D_Fill( st_MX_bck, -INF );
      MATRIX_3D_Fill( st_MX_post, -INF );

      /* run full forward and backward */
      fprintf(stdout, "# ==> Full Forward\n");
      run_Forward_Quad(
            q_seq, t_prof, Q, T, st_MX_fwd, sp_MX_fwd, &sc_fwd_full );
      fprintf(stdout, "# ==> Full Backward\n");
      run_Backward_Quad( 
            q_seq, t_prof, Q, T, st_MX_bck, sp_MX_bck, &sc_bck_full );

      #if DEBUG
      {
         fprintf(stdout, "# fwd score: (full) %7.3f\n", sc_fwd_full);
         fprintf(stdout, "# bck score: (full) %7.3f\n", sc_bck_full);
         DP_MATRIX_Save(Q_size, T_size, st_MX_fwd, sp_MX_fwd, "test_output/my.full_fwd.mx");
         DP_MATRIX_Save(Q_size, T_size, st_MX_bck, sp_MX_bck, "test_output/my.full_bck.mx");
      }
      #endif

      /* compute Posterior (forward * backward) */
      fprintf(stdout, "# ==> Posterior\n");
      run_Decode_Posterior_Quad( q_seq, t_prof, Q_size, T_size, &Q_range, &T_range,
            st_MX_fwd, sp_MX_fwd, st_MX_bck, sp_MX_bck, st_MX_post, sp_MX_post );

      #if DEBUG
      {
         fp = fopen("test_output/my.full_post.mx", "w+");
         DP_MATRIX_Log_Dump(Q_size, T_size, st_MX_post, sp_MX_post, fp);
         fclose(fp);
      }
      #endif
   }
   #endif

   /* METHOD 2: Use normal forward/backward, but constrain model/sequence */
   /* ISSUE: This does not work because it only computes special states on constrained range */
   #if POST_METHOD == 0 || POST_METHOD == 2
   {
      /* constrain sequence and profile to bounding range (for use with standard forward/backward) */
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
      MATRIX_3D_Fill( st_MX_fwd, -INF );
      MATRIX_3D_Fill( st_MX_bck, -INF );
      MATRIX_3D_Fill( st_MX_post, -INF );

      /* run forward/backward on entire of profile/sequence */
      fprintf(stdout, "# ==> Constrained Forward\n");
      run_Forward_Quad(
         q_seq, t_prof, Q_size - 1, T_size - 1, st_MX_fwd, sp_MX_fwd, &sc_fwd_con );
      fprintf(stdout, "# ==> Constrained Backward\n");
      run_Backward_Quad(
         q_seq, t_prof, Q_size - 1, T_size - 1, st_MX_bck, sp_MX_bck, &sc_bck_con );

      /* release constraints */
      SEQUENCE_UnsetSubseq( q_seq );
      HMM_PROFILE_UnsetSubmodel( t_prof ); 
      HMM_PROFILE_ReconfigLength( t_prof, q_seq->N );

      #if DEBUG
      {
         fprintf(stdout, "# fwd score: (con) %7.3f\n", sc_fwd_con);
         fprintf(stdout, "# bck score: (con) %7.3f\n", sc_bck_con);
         DP_MATRIX_Save(Q_size, T_size, st_MX_fwd, sp_MX_fwd, "test_output/my.constrained_fwd.mx");
         DP_MATRIX_Save(Q_size, T_size, st_MX_bck, sp_MX_bck, "test_output/my.constrained_bck.mx");
      }
      #endif

      /* compute Posterior (forward * backward) */
      fprintf(stdout, "# ==> Posterior\n");
      run_Decode_Posterior_Quad( q_seq, t_prof, Q_size, T_size, &Q_range, &T_range,
            st_MX_fwd, sp_MX_fwd, st_MX_bck, sp_MX_bck, st_MX_post, sp_MX_post );

      #if DEBUG
      {
         fp = fopen("test_output/my.constrained_post.mx", "w+");
         DP_MATRIX_Log_Dump(Q_size, T_size, st_MX_post, sp_MX_post, fp);
         fclose(fp);
      }
      #endif
   }
   #endif
   
   /* METHOD 3: Use ranged forward/backward, which only computes for given range */
   #if POST_METHOD == 0 || POST_METHOD == 3
   {
      /* resize special states */
      MATRIX_2D_Reuse( sp_MX_fwd, NUM_SPECIAL_STATES, Q_size + 1 ); 
      MATRIX_2D_Reuse( sp_MX_bck, NUM_SPECIAL_STATES, Q_size + 1 );
      MATRIX_2D_Reuse( sp_MX_post, NUM_SPECIAL_STATES, Q_size + 1 );
      /* resize normal states */
      MATRIX_3D_Reuse( st_MX_fwd, NUM_NORMAL_STATES, Q_size + 1, T_size + 1 );
      MATRIX_3D_Reuse( st_MX_bck, NUM_NORMAL_STATES, Q_size + 1, T_size + 1 );
      MATRIX_3D_Reuse( st_MX_post, NUM_NORMAL_STATES, Q_size + 1, T_size + 1 );
      /* clear previous data */
      MATRIX_3D_Fill( st_MX_fwd, -INF );
      MATRIX_3D_Fill( st_MX_bck, -INF );
      MATRIX_3D_Fill( st_MX_post, -INF );

      /* run ranged forward and backward */
      fprintf(stdout, "# ==> Ranged Forward\n");
      run_Ranged_Forward_Quad(
            q_seq, t_prof, &Q_range, &T_range, st_MX_fwd, sp_MX_fwd, &sc_fwd_rng );
      fprintf(stdout, "# ==> Ranged Backward\n");
      run_Ranged_Backward_Quad( 
            q_seq, t_prof, &Q_range, &T_range, st_MX_bck, sp_MX_bck, &sc_bck_rng );

      #if DEBUG
      {
         fprintf(stdout, "# fwd score: (rng) %7.3f\n", sc_fwd_rng);
         fprintf(stdout, "# bck score: (rng) %7.3f\n", sc_bck_rng);
         DP_MATRIX_Save(Q_size, T_size, st_MX_fwd, sp_MX_fwd, "test_output/my.ranged_fwd.mx");
         DP_MATRIX_Save(Q_size, T_size, st_MX_bck, sp_MX_bck, "test_output/my.ranged_bck.mx");
      }
      #endif

      /* compute Posterior (forward * backward) */
      fprintf(stdout, "# ==> Posterior\n");
      run_Decode_Posterior_Quad( q_seq, t_prof, Q_size, T_size, &Q_range, &T_range,
            st_MX_fwd, sp_MX_fwd, st_MX_bck, sp_MX_bck, st_MX_post, sp_MX_post );

      #if DEBUG
      {
         fp = fopen("test_output/my.ranged_post.mx", "w+");
         DP_MATRIX_Log_Dump(Q_size, T_size, st_MX_post, sp_MX_post, fp);
         fclose(fp);
      }
      #endif
   }
   #endif

   fprintf(stdout, ">> FULL vs RANGED:\n");
   fprintf(stdout, "# fwd score: (full) %7.3f, (rng) %7.3f, (con) %7.3f\n", sc_fwd_full, sc_fwd_rng, sc_fwd_con);
   fprintf(stdout, "# bck score: (full) %7.3f, (rng) %7.3f, (con) %7.3f\n", sc_bck_full, sc_bck_rng, sc_bck_con);

   fprintf(stdout, "# ==> Null2 Compo Bias\n");
   run_Null2_ByExpectation_Quad( q_seq, t_prof, Q_size, T_size, &Q_range, &T_range, 
         st_MX_post, sp_MX_post, dom_def );
   seq_bias = dom_def->seq_bias;

   fprintf(stdout, "# seq_bias: %.9f\n", seq_bias );

   fprintf(stdout, "# ==> Posterior (end)\n");

   return STATUS_SUCCESS;
}


/*! FUNCTION:  run_Decode_Normal_Posterior_Quad()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices (log-space) to create special state posterior into <...post> (normal-space).
 *             Can store <...post> matrix in <...bck>.
 *             Again, input should be in log-space, output matrix is in normal-space (needed for computing null2 score).
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
   printf("=== run_Decode_Posterior_Quad() ===\n");

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

   /* init zero row in query */
   q_0 = 0;
   qx0 = q_0; 

   XMX_X(sp_MX_post, SP_E, q_0) = 0.0;
   XMX_X(sp_MX_post, SP_N, q_0) = 0.0;
   XMX_X(sp_MX_post, SP_J, q_0) = 0.0;
   XMX_X(sp_MX_post, SP_B, q_0) = 0.0;
   XMX_X(sp_MX_post, SP_C, q_0) = 0.0;

   for ( t_0 = 0; t_0 <= T; t_0++ )
   {
      tx0 = t_0;
      /* zero column is -inf in logspace.  We can skip this step and convert to normal space now. */
      MMX_X(st_MX_post, q_0, t_0) = 0.0;
      IMX_X(st_MX_post, q_0, t_0) = 0.0;
      DMX_X(st_MX_post, q_0, t_0) = 0.0; 
   }
   
   /* every position in query */
   for ( q_0 = 1; q_0 <= Q; q_0++ )
   {
      q_1 = q_0 - 1;
      qx0 = q_0;  
      qx1 = qx0 - 1;

      denom = 0.0;

      /* unrolled first loop: special case for left edge of range */
      {
         t_0 = 0;
         tx0 = t_0;
         /* zero column is -inf in logspace.  We can skip this step and convert to normal space now. */
         MMX_X(st_MX_post, q_0, t_0) = 0.0;
         IMX_X(st_MX_post, q_0, t_0) = 0.0;
         DMX_X(st_MX_post, q_0, t_0) = 0.0;
      }

      /* MAIN RECURSION */
      /* FOR every position in TARGET profile */
      for ( t_0 = 1; t_0 < T; t_0++ )
      {
         t_1 = t_0 - 1;
         tx0 = t_0;
         tx1 = tx0 - 1;

         /* normal states */
         mmx = MMX_X(st_MX_fwd, q_0, t_0) + 
               MMX_X(st_MX_bck, q_0, t_0) - 
               overall_sc;
         MMX_X(st_MX_post, q_0, t_0) = exp(mmx);
         denom += MMX_X(st_MX_post, q_0, t_0);

         imx = IMX_X(st_MX_fwd, q_0, t_0) + 
               IMX_X(st_MX_bck, q_0, t_0) - 
               overall_sc;
         IMX_X(st_MX_post, q_0, t_0) = exp(imx);
         denom += IMX_X(st_MX_post, q_0, t_0);
         
         DMX_X(st_MX_post, q_0, t_0) = 0.0;
      }

      /* unrolled final loop: special case for right edge of range */
      {
         t_0 = T;
         t_1 = t_0 - 1;
         tx0 = t_0;
         tx1 = tx0 - 1;

         /* normal states */
         mmx = MMX_X(st_MX_fwd, q_0, t_0) + 
               MMX_X(st_MX_bck, q_0, t_0) -
               overall_sc;
         MMX_X(st_MX_post, q_0, t_0) = exp(mmx);
         denom += MMX_X(st_MX_post, q_0, t_0);

         IMX_X(st_MX_post, q_0, t_0) = 0.0;
         DMX_X(st_MX_post, q_0, t_0) = 0.0;
      }
      
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
      // if (q_0 < 5) printf("denom(q_0 = %d): %9.4f %9.4f\n", q_0, denom, 1.0 / denom);
      denom = 1.0 / denom;

      for ( t_0 = 1; t_0 < T; t_0++ ) {
         MMX_X(st_MX_post, q_0, t_0) *= denom;
         IMX_X(st_MX_post, q_0, t_0) *= denom;
      }
      MMX_X(st_MX_post, q_0, T) *= denom;
      IMX_X(st_MX_post, q_0, T)  = 0.0;

      XMX_X(sp_MX_post, SP_N, q_0) *= denom; 
      XMX_X(sp_MX_post, SP_J, q_0) *= denom; 
      XMX_X(sp_MX_post, SP_C, q_0) *= denom;
   }

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
run_Null2_ByExpectation_Quad(    SEQUENCE*         query,            /* query sequence */
                                 HMM_PROFILE*      target,           /* target hmm model */
                                 int               Q,                /* query length */
                                 int               T,                /* target length */
                                 RANGE*            Q_range,          /* query range */
                                 RANGE*            T_range,          /* target range */
                                 MATRIX_3D*        st_MX_post,       /* posterior normal matrix */
                                 MATRIX_2D*        sp_MX_post,       /* posterior special matrix */
                                 DOMAIN_DEF*       dom_def )         /* OUTPUT: domain def's null2_sc vector */
{
   printf("=== run_Null2_ByExpectation_Quad() ===\n");
   FILE*    fp;
   int      Q_beg, Q_end, Q_len;
   int      T_beg, T_end, T_len;
   int      q_0, t_0; 
   int      qx0, tx0;
   int      st_0, k_0;
   float    mmx, imx, dmx;
   float    x_factor;
   float    null2sc;
   float    val;

   Q_beg = MAX(0, Q_range->beg);
   Q_end = MIN(Q+1, Q_range->end);
   T_beg = MAX(0, T_range->beg);
   T_end = MIN(T+1, T_range->end);
   Q_len = Q_end - Q_beg;
   T_len = T_end - T_beg;

   /* set vectors to proper length */
   MATRIX_2D_Reuse( dom_def->st_freq, T+1, NUM_NORMAL_STATES );
   VECTOR_FLT_SetSize( dom_def->sp_freq,  NUM_SPECIAL_STATES );
   VECTOR_FLT_SetSize( dom_def->null2_sc, NUM_AMINO_PLUS_SPEC );
   VECTOR_FLT_SetSize( dom_def->null2_exp, Q+1 );

   /* fill empty vectors, in normal space */
   MATRIX_2D_Fill(  dom_def->st_freq,   0.0f );
   VECTOR_FLT_Fill( dom_def->sp_freq,   0.0f );
   VECTOR_FLT_Fill( dom_def->null2_sc,  0.0f );
   VECTOR_FLT_Fill( dom_def->null2_exp, 0.0f );

   /* sum over each position in query sequence domain into vectors  */
   for ( q_0 = Q_beg; q_0 < Q_end; q_0++ ) {
      /* for each position in target model */
      qx0 = q_0 - Q_beg;
      for ( t_0 = T_beg; t_0 < T_end; t_0++ ) {
         /* for each normal state emissions */
         tx0 = t_0 - T_beg;
         for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
            // printf("%d/%d, %d/%d, %d/%d\n", q_0, Q_end, t_0, T_end, st_0, NUM_NORMAL_STATES);
            MX_2D( dom_def->st_freq, t_0, st_0 ) += MX_3D( st_MX_post, st_0, qx0, tx0 );
         }
      }
      /* for each special state emissions */
      for ( st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
         VEC_X( dom_def->sp_freq, st_0 ) += MX_2D( sp_MX_post, st_0, q_0 );
      }
   }

   #if DEBUG
   {
      fp = fopen("test_output/my.post_vec.2.quad.csv", "w+");
      fprintf(fp, "<2>\n");
      fprintf(fp, "==> NORMAL STATES (st_freq) <=\n");
      for (int t_0 = T_beg; t_0 < T_end; t_0++) {
         tx0 = t_0 - T_beg;
         fprintf(fp, "ST[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
            MX_2D( dom_def->st_freq, tx0, MAT_ST ), 
            MX_2D( dom_def->st_freq, tx0, INS_ST ), 
            MX_2D( dom_def->st_freq, tx0, DEL_ST ) );
      }
      fprintf(fp, "==> SPECIAL STATES (sp_freq) <=\n");
      for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
         fprintf(fp, "SP[%d]: %12.9f\n", st_0, 
            VEC_X( dom_def->sp_freq, st_0 ) );
      }
      fclose(fp);
   }
   #endif 

   /* convert probabilities to log frequencies */
   /* for each position in query domain */
   for ( t_0 = T_beg; t_0 < T_end; t_0++ ) {
      tx0 = t_0 - T_beg;
      /* for each normal state emissions */
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         MX_2D( dom_def->st_freq, t_0, st_0 ) = log( MX_2D( dom_def->st_freq, t_0, st_0 ));
      }
   }
   /* for each special state emissions */
   for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
      VEC_X( dom_def->sp_freq, st_0 ) = log( VEC_X( dom_def->sp_freq, st_0 ) );
   }

   #if DEBUG
   {
      fp = fopen("test_output/my.post_vec.3.quad.csv", "w+");
      fprintf(fp, "<3>\n");
      fprintf(fp, "==> NORMAL STATES (st_freq) <=\n");
      for (int t_0 = T_beg; t_0 < T_end; t_0++) {
         tx0 = t_0 - T_beg;
         fprintf(fp, "ST[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
            MX_2D( dom_def->st_freq, tx0, MAT_ST ), 
            MX_2D( dom_def->st_freq, tx0, INS_ST ), 
            MX_2D( dom_def->st_freq, tx0, DEL_ST ) );
      }
      fprintf(fp, "==> SPECIAL STATES (sp_freq) <=\n");
      for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
         fprintf(fp, "SP[%d]: %12.9f\n", st_0, 
            VEC_X( dom_def->sp_freq, st_0 ) );
      }
      fclose(fp);
   }
   #endif 

   /*  This is the cumulative score contributed by each index in the model.
    *  Divide by Q to take the average per cell (or subtract by the log(Q) to compute in log space) 
    */
   float neglog_Q = -log( (float)Q );
   printf("neglog_Q = %d %f\n", Q_len, neglog_Q);
   /* for each position in query domain */
   for ( t_0 = T_beg; t_0 < T_end; t_0++ ) {
      tx0 = t_0 - T_beg;
      /* for each normal state emissions */
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         MX_2D( dom_def->st_freq, t_0, st_0 ) += neglog_Q;
      }
   }
   /* for each special state emissions */
   for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
      VEC_X( dom_def->sp_freq, st_0 ) += neglog_Q;
   }

   #if DEBUG
   {
      fp = fopen("test_output/my.post_vec.4.quad.csv", "w+");
      fprintf(fp, "<4>\n");
      fprintf(fp, "==> NORMAL STATES (st_freq) <=\n");
      for (int t_0 = T_beg; t_0 < T_end; t_0++) {
         tx0 = t_0 - T_beg;
         fprintf(fp, "ST[%3d]: %12.9f %12.9f %12.9f\n", t_0, 
            MX_2D( dom_def->st_freq, tx0, MAT_ST ), 
            MX_2D( dom_def->st_freq, tx0, INS_ST ), 
            MX_2D( dom_def->st_freq, tx0, DEL_ST ) );
      }
      fprintf(fp, "==> SPECIAL STATES (sp_freq) <=\n");
      for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
         fprintf(fp, "SP[%d]: %12.9f\n", st_0, 
            VEC_X( dom_def->sp_freq, st_0 ) );
      }
      fclose(fp);
   }
   #endif 

   /* x-factor: */
   x_factor = VEC_X( dom_def->sp_freq, SP_N);
   x_factor = logsum( x_factor,
                      VEC_X(dom_def->sp_freq, SP_C) );
   x_factor = logsum( x_factor,
                      VEC_X(dom_def->sp_freq, SP_J) );
   
  #if DEBUG
   {
      fprintf(stdout, "xfactor: %9.4f,  NCJ: %9.4f %9.4f %9.4f\n", 
         x_factor, VEC_X( dom_def->sp_freq, SP_N), VEC_X(dom_def->sp_freq, SP_C), VEC_X(dom_def->sp_freq, SP_J) );
   }
   #endif

   /* initialize null2 vector for log-space */
   VECTOR_FLT_Fill( dom_def->null2_sc, -INF );

   /* null2 log emissions probabilities found by summing over 
    * all emissions used in paths explaining the domain. 
    */
   /* for each amino acid */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) 
   {
      /* for each position in model */
      for ( t_0 = T_beg; t_0 < T_end - 1; t_0++ ) 
      {
         tx0 = t_0 - T_beg;
         /* look at the log frequencies (weighted probability of position in model contributed to path score ) 
          * at the model position and multiply them by the score contribution
          */ 
         mmx = MX_2D( dom_def->st_freq, t_0, MAT_ST) + MSC_X(target, t_0, k_0);
         imx = MX_2D( dom_def->st_freq, t_0, INS_ST) + ISC_X(target, t_0, k_0);

         VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ), mmx );
         VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ), imx );
      }
      t_0 = T_end - 1;
      tx0 = t_0 - T_beg;
      mmx = MX_2D( dom_def->st_freq, t_0, MAT_ST) + MSC_X(target, t_0, k_0);
      
      VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ), mmx);
      VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ), x_factor );                            
   }

   #if DEBUG
   {
      // fp = stdout;
      fp = fopen("test_output/my.null2.quad.csv", "w+");
      fprintf(fp, "==> NULL2 (pre) <=\n");
      for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
         fprintf(fp, "%d %c %.7f\n", k_0, AA[k_0], VEC_X( dom_def->null2_sc, k_0) );
      }
   }
   #endif
   
   /* convert log to normal scores (for computing averages) */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = exp( VEC_X( dom_def->null2_sc, k_0 ) );
   }

   #if DEBUG
   {
      fprintf(fp, "==> NULL2 (log->normal) <=\n");
      for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
         fprintf(fp, "%d %c %.7f\n", k_0, AA[k_0], VEC_X( dom_def->null2_sc, k_0) );
      }
   }
   #endif

   /* compute the degeneracy characters by averaging the odds ratios of degeneracy matches (currently only for X character)  */
   VEC_X( dom_def->null2_sc, AMINO_X ) = 0.0f;
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) {
      VEC_X( dom_def->null2_sc, AMINO_X ) += VEC_X( dom_def->null2_sc, k_0 );
   }
   VEC_X( dom_def->null2_sc, AMINO_X ) /= NUM_AMINO;

   /* set scores for all missing, non-residue, and gap characters */
   VEC_X( dom_def->null2_sc, AMINO_MIS ) = 1.0;
   VEC_X( dom_def->null2_sc, AMINO_NON ) = 1.0;
   VEC_X( dom_def->null2_sc, AMINO_GAP ) = 1.0;

   #if DEBUG
   {
      fprintf(fp, "==> NULL2 (end) <=\n");
      for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
         fprintf(fp, "%d %c %.7f\n", k_0, AA[k_0], VEC_X( dom_def->null2_sc, k_0) );
      }
      if (fp != stdout) fclose(fp);
   }
   #endif

   /*  Now have the expected bias per symbol in the alphabet. 
    *  Insert these scores into a sequence according to there expected bias.
    *  NOTE: Should be able to remove this for main pipeline.
    */
   for ( q_0 = Q_beg; q_0 < Q_end; q_0++ ) {
      qx0 = q_0 - Q_beg;
      k_0 = AA_REV[query->seq[q_0]];
      VEC_X( dom_def->null2_exp, qx0) = logf( VEC_X( dom_def->null2_sc, k_0 ));
   }

   /* Finally, for bias correction, sum over the entire expected sequence bias. */
   dom_def->seq_bias = 0.0f;
   for ( q_0 = Q_beg; q_0 < Q_end; q_0++ ) {
      qx0 = q_0 - Q_beg;
      k_0 = AA_REV[query->seq[q_0]];
      dom_def->seq_bias += logf( VEC_X( dom_def->null2_sc, k_0 ));
   }

   #if DEBUG
   {
      // fp = stdout;
      fp = fopen("test_output/my.null2sc.quad.csv", "w+");
      fprintf(fp, "# === NULL2SC (Expectation by Amino) ===\n");
      for (int k_0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++) {
         fprintf(fp, "%d %c %.9f\n", k_0, AA[k_0], VEC_X( dom_def->null2_sc, k_0 ));
      }
      fprintf(fp, "# === NULL2SC (Expectation by Position) ===\n");
      for ( q_0 = Q_beg; q_0 < Q_end; q_0++ ) {
         qx0 = q_0 - Q_beg;
         fprintf(fp, "%d %.9f\n", q_0, VEC_X( dom_def->null2_exp, qx0 ));
      }
      fprintf(fp, "# === DOMCORRECTION: %.9f\n", dom_def->seq_bias);
      if (fp != stdout) fclose(fp);

      printf("# === DOMCORRECTION: %.9f\n", dom_def->seq_bias);
   }
   #endif

   return STATUS_SUCCESS;
}
