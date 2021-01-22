/*******************************************************************************
 *  FILE:      posterior_sparse.h
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
#include "../algs_linear/algs_linear.h"
#include "../parsers/parsers.h"

/* header */
#include "_algs_sparse.h"
#include "posterior_sparse.h"

/*! FUNCTION:  run_Posterior_Sparse()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward <st_MX_bck>.
 *             Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in supplied <st_MX_post> (can override input matrices).
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
STATUS_FLAG 
run_Posterior_Sparse(   SEQUENCE*               q_seq,            /* query sequence */
                        HMM_PROFILE*            t_prof,           /* target hmm model */
                        int                     Q,                /* query length */
                        int                     T,                /* target length */
                        HMM_BG*                 hmm_bg,           /* hmm background model */
                        EDGEBOUNDS*             edg,              /* edgebounds */
                        MATRIX_3D_SPARSE*       st_SMX_fwd,       /* normal state matrix for forward */
                        MATRIX_2D*              sp_MX_fwd,        /* special state matrix for forward */
                        MATRIX_3D_SPARSE*       st_SMX_bck,       /* normal state matrix for backward */
                        MATRIX_2D*              sp_MX_bck,        /* special state matrix for backward */
                        MATRIX_3D_SPARSE*       st_SMX_post,      /* OUTPUT: normal state matrix for posterior */
                        MATRIX_2D*              sp_MX_post,       /* OUTPUT: special state matrix for posterior */   
                        MATRIX_3D_SPARSE*       st_SMX_opt,       /* OUTPUT: normal state matrix for optimal accuracy */
                        MATRIX_2D*              sp_MX_opt,        /* OUTPUT: special state matrix for optimal accuracy */       
                        RESULT*                 result,           /* OUPUT: full cloud results */
                        DOMAIN_DEF*             dom_def,          /* OUTPUT: domain data */
                        bool                    is_run_domains )  /* if run posterior on domains */  
{
   printf("=== run_Posterior_Sparse() [BEGIN] ===\n");
   FILE*    fp;

   /* working space needed for computing optimal alignment. Can recycle <...fwd> */
   EDGEBOUNDS*    edg_dom;
   RANGE          Q_range, T_range;
   int            Q_size, T_size, D_size;
   RANGE          D_range;
   DOMAIN_X       domain;
   float          compo_bias, null_sc;
   float          pre_sc, fwd_sc, bck_sc, post_sc, opt_sc, dom_sc;
   int            N_domains;
   int            D_total;

   /* init and clear data from previous search */
   DOMAIN_DEF_Reuse( dom_def );
   edg_dom = dom_def->edg;
   
   /* find the min/max of query and target range in core model */
   EDGEBOUNDS_Find_BoundingBox( edg, &Q_range, &T_range );
   Q_size = Q_range.end - Q_range.beg;
   T_size = T_range.end - T_range.beg;
   printf("BOUNDING: { Q:(%d,%d), T:(%d,%d) }\n",
      Q_range.beg, Q_range.end, T_range.beg, T_range.end);

   /* store begin/end points */
   result->query_start  = Q_range.beg;
   result->query_end    = Q_range.end;
   result->target_start = T_range.beg;
   result->target_start = T_range.end;

   #if DEBUG
   {
      EDGEBOUNDS_Save( edg, "test_output/my.cloud.edg" );
   }
   #endif

   /* compute score and bias for entire cloud */
   {
      /* compute Forward */
      run_Bound_Forward_Sparse( 
         q_seq, t_prof, Q, T, st_SMX_fwd, sp_MX_fwd, edg, NULL, &fwd_sc );
      fprintf(stdout, "# ==> Forward  (full cloud): %11.4f %11.4f\n", 
         fwd_sc, fwd_sc/CONST_LOG2);
      result->bound_fwd_natsc = fwd_sc;
      #if DEBUG
      {
         DP_MATRIX_Clean(Q, T, debugger->test_MX, NULL);
         MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_fwd, debugger->test_MX);
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_fwd, "test_output/my.sparse_fwd.mx");
      }
      #endif

      /* compute Backward */
      run_Bound_Backward_Sparse( 
         q_seq, t_prof, Q, T, st_SMX_bck, sp_MX_bck, edg, NULL, &bck_sc );
      fprintf(stdout, "# ==> Backward (full cloud): %11.4f %11.4f\n", 
         bck_sc, bck_sc/CONST_LOG2);
      result->bound_bck_natsc = bck_sc;
      #if DEBUG && VIZ
      {
         DP_MATRIX_Clean(Q, T, debugger->test_MX, NULL);
         MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_bck, debugger->test_MX);
         DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_bck, "test_output/my.sparse_bck.mx");
      }
      #endif

      /* find the Domain ranges */
      fprintf(stdout, "# ==> Find Domains\n");
      run_Decode_Domains( 
         q_seq, t_prof, Q, T, sp_MX_fwd, sp_MX_bck, dom_def );
      N_domains = dom_def->dom_ranges->N;
      /* list found domains */
      fprintf(stdout, "DOMAINS FOUND: %d\n", dom_def->dom_ranges->N);
      for (int i = 0; i < dom_def->dom_ranges->N; i++) 
      {
         RANGE r_0 = VEC_X( dom_def->dom_ranges, i );
         fprintf(stdout, "[%d] {%d,%d}\n", i, r_0.beg, r_0.end);
      }

      /* compute Posterior */
      fprintf(stdout, "# ==> Posterior (full cloud)\n");
      run_Decode_Posterior_Sparse( 
         q_seq, t_prof, Q, T, edg, NULL,
         st_SMX_fwd, sp_MX_fwd, st_SMX_bck, sp_MX_bck, st_SMX_post, sp_MX_post );

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

      /* Composition Bias */
      fprintf(stdout, "# ==> Null2 Compo Bias (full cloud)\n");
      run_Null2_ByExpectation_Sparse( 
         q_seq, t_prof, Q, T, edg, NULL, NULL, NULL,
         st_SMX_post, sp_MX_post, dom_def, &compo_bias );
      fprintf(stdout, "# ==> Null2 Compo Bias (full cloud): %11.4f %11.4f\n", compo_bias, compo_bias/CONST_LOG2);
      result->final_scores.seq_bias;

      // /* Optimal Alignment */
      // fprintf(stdout, "# ==> Optimal Accuracy (full cloud)\n");
      // run_Posterior_Optimal_Accuracy_Sparse( 
      //    q_seq, t_prof, Q, T, edg, NULL,
      //    st_SMX_post, sp_MX_post, st_SMX_opt, sp_MX_opt, &opt_sc );
      // fprintf(stdout, "# ==> Optimal Accuracy (full cloud): %11.4f\n", opt_sc);

      #if DEBUG 
      {
         fp = fopen("test_output/my.sparse_opt.mx", "w+");
         MATRIX_3D_SPARSE_Log_Embed(Q, T, st_SMX_opt, debugger->test_MX);
         MATRIX_2D_Log(sp_MX_opt);
         DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_opt, fp);
         MATRIX_2D_Exp(sp_MX_post);
         fclose(fp);
      }
      #endif

      // /* Optimal Alignment Traceback */
      // fprintf(stdout, "# ==> Optimal Alignment (full cloud)\n");
      // run_Posterior_Optimal_Traceback_Sparse( 
      //    q_seq, t_prof, Q, T, edg, NULL, 
      //    st_SMX_post, sp_MX_post, dom_def->align );
      // ALIGNMENT_Build_MMSEQS_Style( dom_def->align, q_seq, t_prof );
      // ALIGNMENT_Build_HMMER_Style( dom_def->align, q_seq, t_prof );
      // fprintf(stdout, "# ==> Cigar Alignment: %s", dom_def->align->cigar_aln);
   }

   /* run through Domains and compute score, bias correction, and optimal alignment */
   printf("is_run_domains: %d\n", is_run_domains);
   if ( is_run_domains == true )
   {
      dom_def->n_domains      = N_domains;
      dom_def->dom_sumsc      = 0.0f;
      dom_def->dom_sumbias    = 0.0f;
      D_total                 = 0;
      null_sc                 = dom_def->nullsc;

      for (int i = 0; i < N_domains; i++)
      {
         D_range = VEC_X( dom_def->dom_ranges, i );
         printf("Domain (%d of %d): {%d,%d}\n", 
            i+1, N_domains, D_range.beg, D_range.end);

         // /* Reparameterize sequence to only cover domain range. */
         D_size = D_range.end - D_range.beg + 1;
         // SEQUENCE_Set_Domain( q_seq, D_range );
         // EDGEBOUNDS_Set_Domain( edg, edg_dom, D_range );
         // HMM_PROFILE_ReconfigLength( t_prof, q_seq->N );
         
         /* clear previous data */
         MATRIX_3D_SPARSE_Fill_Outer( st_SMX_fwd, -INF );
         MATRIX_3D_SPARSE_Fill_Outer( st_SMX_bck, -INF );
         if (st_SMX_bck != st_SMX_post) {
            MATRIX_3D_SPARSE_Fill_Outer( st_SMX_post, -INF );
         }

         /* compute Forward/Backward for the domain range */
         run_Bound_Forward_Sparse( 
            q_seq, t_prof, Q, T, st_SMX_fwd, sp_MX_fwd, edg, &D_range, &fwd_sc );
         fprintf(stdout, "# ==> Forward   (domain %d/%d): %11.4f %11.4f\n", 
            i+1, N_domains, fwd_sc, fwd_sc/CONST_LOG2);
         #if DEBUG 
         {
            fp = fopen("test_output/my.sparse_fwd.dom.mx", "w+");
            MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_fwd, debugger->test_MX);
            DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_fwd, fp);
            fclose(fp);
         }
         #endif

         /* compute Forward/Backward for the domain range */
         run_Bound_Backward_Sparse( 
            q_seq, t_prof, Q, T, st_SMX_bck, sp_MX_bck, edg, &D_range, &bck_sc );
         fprintf(stdout, "# ==> Backward  (domain %d/%d): %11.4f %11.4f\n", 
            i+1, N_domains, bck_sc, bck_sc/CONST_LOG2);
         #if DEBUG 
         {
            fp = fopen("test_output/my.sparse_bck.dom.mx", "w+");
            MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_bck, debugger->test_MX);
            DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_bck, fp);
            fclose(fp);
         }
         #endif

         /* compute Posterior (forward * backward) for domain range */
         fprintf(stdout, "# ==> Posterior (domain %d/%d)\n", 
            i+1, N_domains);
         run_Decode_Posterior_Sparse( 
            q_seq, t_prof, Q, T, edg, &D_range,
            st_SMX_fwd, sp_MX_fwd, st_SMX_bck, sp_MX_bck, st_SMX_post, sp_MX_post );
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
         fprintf(stdout, "# ==> Null2 Compo Bias\n");
         run_Null2_ByExpectation_Sparse( q_seq, t_prof, Q, T, edg, &D_range, NULL, NULL,
            st_SMX_post, sp_MX_post, dom_def, &compo_bias );
         fprintf(stdout, "# ==> Null2 Compo Bias  (domain %d/%d): %11.4f %11.4f\n", 
            i+1, N_domains, compo_bias, compo_bias/CONST_LOG2);

         /** TODO: Get optimal alignment */ 
         
         /* add domain data */
         VECTOR_FLT_Pushback( dom_def->dom_fwdsc, fwd_sc );
         VECTOR_FLT_Pushback( dom_def->dom_bias, compo_bias );

         /* check if best score */
         pre_sc = (fwd_sc - (null_sc)) / CONST_LOG2;
         dom_sc = (fwd_sc - (null_sc + compo_bias)) / CONST_LOG2;
         if ( dom_sc > dom_def->best_sc )
         {
            dom_def->best        = i;
            dom_def->best_sc     = dom_sc;
            dom_def->best_fwdsc  = fwd_sc;
            dom_def->best_presc  = pre_sc;
            dom_def->best_bias   = compo_bias;
            dom_def->best_range  = D_range;
         }

         /* constructed score over all domains */
         dom_def->dom_sumsc   += fwd_sc;
         dom_def->dom_sumbias += compo_bias;
         D_total              += (D_range.end - D_range.beg + 1);
      }

      /* final reconstructed score */
      if (N_domains > 0)
      {
         /* constructed score over all domains */
         dom_def->dom_sumbias = logsum(0.0f, log(hmm_bg->omega) + dom_def->dom_sumbias);
         dom_def->dom_sumsc  += (Q - D_total) * log((float) Q / (float) (Q + 3));
         dom_def->dom_sumsc   = (dom_def->dom_sumsc - (dom_def->nullsc + dom_def->dom_sumbias)) / CONST_LOG2;
      }
      else {
         dom_def->dom_sumbias = -INF;
         dom_def->dom_sumsc   = -INF;
      }
   }

   printf("=== run_Posterior_Sparse() [END] ===\n");
   return STATUS_SUCCESS;
}


/*! FUNCTION:  run_Decode_Normal_Posterior_Sparse()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store output <...post> matrix in input <...bck> matrix.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
STATUS_FLAG
run_Decode_Posterior_Sparse(  SEQUENCE*            q_seq,            /* query sequence */
                              HMM_PROFILE*         t_prof,           /* target hmm model */
                              int                  Q,                /* query length */
                              int                  T,                /* target length */
                              EDGEBOUNDS*          edg,              /* edgebounds */
                              RANGE*               dom_range,        /* OPTIONAL: domain range */
                              MATRIX_3D_SPARSE*    st_SMX_fwd,       /* normal state matrix for forward */
                              MATRIX_2D*           sp_MX_fwd,        /* special state matrix for forward */
                              MATRIX_3D_SPARSE*    st_SMX_bck,       /* normal state matrix for backward */
                              MATRIX_2D*           sp_MX_bck,        /* special state matrix for backward */
                              MATRIX_3D_SPARSE*    st_SMX_post,      /* OUTPUT: normal state matrix for posterior */
                              MATRIX_2D*           sp_MX_post )      /* OUTPUT: normal state matrix for posterior */
{
   printf("=== run_Decode_Posterior_Sparse ===\n");
   FILE* fp = NULL;

   /* query index */
   int      q_0, q_1;
   int      qx0, qx1;
   /* target index */ 
   int      t_0, t_1; 
   int      tx0, tx1;
   /* state index */
   int      st_0;
   /* edgebound index */
   int      i_0;
   int      r_0, r_0b, r_0e;
   BOUND*   bnd;
   RANGE    Q_range;
   RANGE    T_range;
   int      id, lb_0, rb_0;
   int      lb_T, rb_T;
   int      len;
   /* overall score */
   float    overall_sc;
   /* common scale factor denominator */
   float    denom;
   /* temp mx scores */
   float    mmx, imx, dmx, smx;
   float    mmx_, imx_, dmx_, smx_;
   /* check if query position is in the domain */
   bool     is_q_0_in_dom_range;
   bool     is_q_1_in_dom_range;

   /* --------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   overall_sc  =  XMX_X( sp_MX_fwd, SP_C, Q ) + 
                  XSC_X(t_prof, SP_C, SP_MOVE);
   // printf("overall_sc: %f %f ==> %f\n", XMX_X( sp_MX_fwd, SP_C, Q ), XSC_X(t_prof, SP_C, SP_MOVE), overall_sc);
   
   /* domain range (query sequence) */
   if (dom_range == NULL) {
      Q_range.beg = 0;
      Q_range.end = Q + 1;
   } else {
      Q_range = *dom_range;
   }
   /* target range */
   T_range.beg = 0;
   T_range.end = T + 1;

   /* COMPUTE POSTERIOR */
   if (true)
   {
      /* init index */
      q_0 = 0;
      r_0b = r_0e = 0;
      /* check if query position is in domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );
      /* get edgebound range */
      EDGEBOUNDS_NxtRow( edg, &r_0b, &r_0e, q_0 );

      /* initialize special states for zero row */
      XMX_X(sp_MX_post, SP_E, q_0) = -INF;
      XMX_X(sp_MX_post, SP_N, q_0) = -INF;
      XMX_X(sp_MX_post, SP_J, q_0) = -INF;
      XMX_X(sp_MX_post, SP_B, q_0) = -INF;
      XMX_X(sp_MX_post, SP_C, q_0) = -INF;

      /* if sequence position is in domain range */ 
      if ( is_q_0_in_dom_range == true )
      {
         /* FOR every BOUND in zero ROW */
         for (r_0 = r_0b; r_0 < r_0e; r_0++)
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0);
            id    = bnd->id;
            lb_T  = bnd->lb <= 0;
            lb_0  = MAX(bnd->lb, T_range.beg);   /* can't overflow left edge */
            rb_T  = bnd->rb >= T;
            rb_0  = MIN(bnd->rb, T_range.end);   /* can't overflow right edge */

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0 = VECTOR_INT_Get( st_SMX_fwd->imap_cur, r_0 );    /* (q_0, t_0) location offset */

            /* FOR every position in TARGET profile */
            for (t_0 = lb_0; t_0 < rb_0; t_0++)
            {
               tx0 = t_0 - bnd->lb;
               /* zero column is -inf in logspace.  We can skip this step and convert to normal space now. */
               MSMX_X(st_SMX_post, qx0, tx0) = -INF;
               ISMX_X(st_SMX_post, qx0, tx0) = -INF;
               DSMX_X(st_SMX_post, qx0, tx0) = -INF; 
            }
         }
      }

      /* FOR every position in QUERY sequence (row in matrix) */
      for (q_0 = 1; q_0 <= Q; q_0++)
      {
         q_1 = q_0 - 1;
          
         /* check if query position is in domain */
         is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );
         /* get edgebound range */
         EDGEBOUNDS_NxtRow( edg, &r_0b, &r_0e, q_0 );

         /* if sequence position is in domain range */ 
         if ( is_q_0_in_dom_range == true )
         {
            /* FOR every BOUND in current ROW */
            for (r_0 = r_0b; r_0 < r_0e; r_0++)
            {
               /* get bound data */
               bnd   = &EDG_X(edg, r_0);
               id    = bnd->id;
               lb_T  = bnd->lb <= 0;
               lb_0  = MAX(bnd->lb, T_range.beg);   /* can't overflow left edge */
               rb_T  = bnd->rb >= T;
               rb_0  = MIN(bnd->rb, T_range.end);   /* can't overflow right edge */

               /* fetch data mapping bound start location to data block in sparse matrix */
               qx0 = VECTOR_INT_Get( st_SMX_fwd->imap_cur, r_0 );    /* (q_0, t_0) location offset */
               qx1 = VECTOR_INT_Get( st_SMX_fwd->imap_prv, r_0 );    /* (q_1, t_0) location offset */

               /* unrolled first loop: special case for left edge of range */
               if ( true )
               {
                  t_0 = lb_0;
                  tx0 = t_0 - bnd->lb;

                  /* zero column is -inf in logspace.  We can skip this step and convert to normal space now. */
                  MSMX_X(st_SMX_post, qx0, tx0) = -INF;
                  ISMX_X(st_SMX_post, qx0, tx0) = -INF;
                  DSMX_X(st_SMX_post, qx0, tx0) = -INF; 
               }

               /* MAIN RECURSION */
               /* FOR every position in TARGET profile */
               for (t_0 = lb_0 + 1; t_0 < rb_0 - 1; t_0++)
               {
                  tx0 = t_0 - bnd->lb;

                  /* normal states */
                  mmx = MSMX_X(st_SMX_fwd, qx0, tx0) +
                        MSMX_X(st_SMX_bck, qx0, tx0) -
                        overall_sc;
                  MSMX_X(st_SMX_post, qx0, tx0) = mmx;

                  imx = ISMX_X(st_SMX_fwd, qx0, tx0) + 
                        ISMX_X(st_SMX_bck, qx0, tx0) -
                        overall_sc;
                  ISMX_X(st_SMX_post, qx0, tx0) = imx;

                  DSMX_X(st_SMX_post, qx0, tx0) = -INF;
               }

               /* unrolled final loop: special case for right edge of range */
               if (true && rb_0 > 1)  
               {
                  t_0 = rb_0 - 1;
                  tx0 = t_0 - bnd->lb;

                  /* normal states */
                  mmx = MSMX_X(st_SMX_fwd, qx0, tx0) + 
                        MSMX_X(st_SMX_bck, qx0, tx0) -
                        overall_sc;
                  MSMX_X(st_SMX_post, qx0, tx0) = mmx;

                  ISMX_X(st_SMX_post, qx0, tx0) = -INF;
                  DSMX_X(st_SMX_post, qx0, tx0) = -INF;
               }
            }
         }

         /* special states */
         XMX_X(sp_MX_post, SP_E, q_0) = -INF;
         XMX_X(sp_MX_post, SP_B, q_0) = -INF;

         smx = XMX_X(sp_MX_fwd, SP_N, q_1) +
               XMX_X(sp_MX_bck, SP_N, q_0) +
               XSC_X(t_prof, SP_N, SP_LOOP) -
               overall_sc;
         XMX_X(sp_MX_post, SP_N, q_0) = smx;

         smx = XMX_X(sp_MX_fwd, SP_J, q_1) +
               XMX_X(sp_MX_bck, SP_J, q_0) +
               XSC_X(t_prof, SP_J, SP_LOOP) -
               overall_sc;
         XMX_X(sp_MX_post, SP_J, q_0) = smx;

         smx = XMX_X(sp_MX_fwd, SP_C, q_1) +
               XMX_X(sp_MX_bck, SP_C, q_0) +
               XSC_X(t_prof, SP_C, SP_LOOP) -
               overall_sc;
         XMX_X(sp_MX_post, SP_C, q_0) = smx;
      }
   }

   /* CONVERT TO NORMAL SPACE */
   if ( true )
   {
      MATRIX_3D_SPARSE_Exp( st_SMX_post );
      MATRIX_2D_Exp( sp_MX_post );
   }

   /* NORMALIZE MATRIX */
   if ( true )
   {
      /* init indexes */
      q_0  = 0;
      r_0b = r_0e = 0;  
      /* check if query position is in domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );
      /* get edgebound range */
      EDGEBOUNDS_NxtRow( edg, &r_0b, &r_0e, q_0 );

      /* FOR every position in QUERY sequence (row in matrix) */
      for (q_0 = 1; q_0 <= Q; q_0++)
      {
         q_1 = q_0 - 1;
         t_0 = 0;
         /* check if query position is in domain */
         is_q_0_in_dom_range = (q_0 >= Q_range.beg && q_0 < Q_range.end);
         /* get edgebound range */
         EDGEBOUNDS_NxtRow( edg, &r_0b, &r_0e, q_0 );

         denom = 0.0;

         /* only compute if in domain range */
         if ( is_q_0_in_dom_range == true )
         {
            /* FOR every BOUND in current ROW */
            for (r_0 = r_0b; r_0 < r_0e; r_0++)
            {
               /* get bound data */
               bnd   = &EDG_X(edg, r_0);
               id    = bnd->id;
               lb_T  = bnd->lb <= 0;
               lb_0  = MAX(bnd->lb, T_range.beg);   /* can't overflow left edge */
               rb_T  = bnd->rb >= T;
               rb_0  = MIN(bnd->rb, T_range.end);   /* can't overflow right edge */

               /* fetch data mapping bound start location to data block in sparse matrix */
               qx0 = VECTOR_INT_Get( st_SMX_fwd->imap_cur, r_0 );    /* (q_0, t_0) location offset */
               qx1 = VECTOR_INT_Get( st_SMX_fwd->imap_prv, r_0 );    /* (q_1, t_0) location offset */

               /* MAIN RECURSION */
               /* FOR every position in TARGET profile */
               for (t_0 = lb_0 + 1; t_0 < rb_0 - 1; t_0++)
               {
                  t_1 = t_0 - 1; 
                  tx0 = t_0 - bnd->lb;
                  tx1 = tx0 - 1;

                  /* normal states */
                  denom += MSMX_X(st_SMX_post, qx0, tx0);
                  denom += ISMX_X(st_SMX_post, qx0, tx0);
               }

               /* unrolled final loop: special case for right edge of range */
               if ( true && (rb_0 > 1) )  
               {
                  denom += MSMX_X(st_SMX_post, qx0, tx0);
               }
            }
         }

         /* special states */
         denom += XMX_X(sp_MX_post, SP_N, q_0) + 
                  XMX_X(sp_MX_post, SP_J, q_0) + 
                  XMX_X(sp_MX_post, SP_C, q_0);

         /* normalize by scaling row by common factor denominator */
         // if (q_0 < 5) printf("denom(q_0 = %d): %9.4f %9.4f\n", q_0, denom, 1.0 / denom);
         denom = 1.0 / denom;

         /* apply denominator scaling factor to entire row */
         if ( is_q_0_in_dom_range == true )
         {
            /* FOR every BOUND in current ROW */
            for (r_0 = r_0b; r_0 < r_0e; r_0++)
            {
               /* get bound data */
               bnd   = &EDG_X(edg, r_0);
               id    = bnd->id;
               lb_T  = bnd->lb <= 0;
               lb_0  = MAX(bnd->lb, T_range.beg);   /* can't overflow left edge */
               rb_T  = bnd->rb >= T;
               rb_0  = MIN(bnd->rb, T_range.end);   /* can't overflow right edge */

               /* fetch data mapping bound start location to data block in sparse matrix */
               qx0 = VECTOR_INT_Get( st_SMX_fwd->imap_cur, r_0 );    /* (q_0, t_0) location offset */
               qx1 = VECTOR_INT_Get( st_SMX_fwd->imap_prv, r_0 );    /* (q_1, t_0) location offset */

               /* if left bound touches the left edge */
               if (true)
               {
                  t_0 = lb_0;
                  tx0 = t_0 - bnd->lb;
                  /* normalize by scaling row by common factor denominator */
                  MSMX_X(st_SMX_post, qx0, tx0) *= denom;
                  ISMX_X(st_SMX_post, qx0, tx0)  = 0.0;
                  DSMX_X(st_SMX_post, qx0, tx0)  = 0.0;
               }

               /* normalize all values to fraction of total column weight */
               /* FOR every position in TARGET profile */
               for (t_0 = lb_0 + 1; t_0 < rb_0 - 1; t_0++)
               {
                  t_1 = t_0 - 1; 
                  tx0 = t_0 - bnd->lb;
                  tx1 = tx0 - 1;

                  /* normalize by scaling row by common factor denominator */
                  MSMX_X(st_SMX_post, qx0, tx0) *= denom;
                  ISMX_X(st_SMX_post, qx0, tx0) *= denom;
                  DSMX_X(st_SMX_post, qx0, tx0)  = 0.0;
               }

               /* if right bound touches the right edge */
               if ( true && rb_0 > 1 ) 
               {
                  t_0 = rb_0 - 1;
                  t_1 = t_0 - 1;
                  tx0 = t_0 - bnd->lb;
                  tx1 = tx0 - 1;
                  /* normalize by scaling row by common factor denominator */
                  MSMX_X(st_SMX_post, qx0, tx0) *= denom;
                  ISMX_X(st_SMX_post, qx0, tx0)  = 0.0;
                  DSMX_X(st_SMX_post, qx0, tx0)  = 0.0;
               }
            }
         }

         /* normalize by scaling row by common factor denominator */
         XMX_X(sp_MX_post, SP_N, q_0) *= denom; 
         XMX_X(sp_MX_post, SP_J, q_0) *= denom; 
         XMX_X(sp_MX_post, SP_C, q_0) *= denom;
      }
   }
   
   return STATUS_SUCCESS;
}


/*! FUNCTION:  run_Decode_Domains()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward <st_MX_bck>.
 *             Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in supplied <st_MX_post> (can override input matrices).
 *             NOTE: Modeled after <p7_DomainDecoding()>.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
STATUS_FLAG
run_Decode_Domains(     SEQUENCE*               q_seq,            /* query sequence */
                        HMM_PROFILE*            t_prof,           /* target hmm model */
                        int                     Q,                /* query length */
                        int                     T,                /* target length */
                        MATRIX_2D*              sp_MX_fwd,        /* special state matrix for forward */
                        MATRIX_2D*              sp_MX_bck,        /* special state matrix for backward */   
                        DOMAIN_DEF*             dom_def )         /* OUTPUT: domain data */
{
   FILE*          fp;
   VECTOR_FLT*    idx;

   int      q_0, q_1;
   int      q_beg, q_end;
   int      t_0, t_1; 
   int      t_beg, t_end;
   double   b_0, e_0;
   double   btot_add, etot_add;
   double   njcp, np, jp, cp;
   float    mocc_0, btot_0, btot_1, etot_0, etot_1;
   float    rt1_test, rt2_btest, rt2_etest;
   float    rt1_max, rt2_bmax, rt2_emax;
   float    rt1_crit, rt2_crit;
   float    rt1, rt2, rt3;
   float    overall_logp;
   bool     is_in_domain; 
   bool     rt2_bhit, rt2_ehit;
   bool     is_multiple_domains;

   /* scalar for preventing underflow error */
   overall_logp = XMX_X( sp_MX_fwd, SP_C, Q ) + 
                  XSC_X(t_prof, SP_C, SP_MOVE);
   // printf("# overall_logp: %7.4f, Q = %d\n", overall_logp, Q );

   /* domain threshold test 1: */
   rt1 = dom_def->rt1;
   /* domain threshold test 2: */
   rt2 = dom_def->rt2;

   VECTOR_FLT_SetSize( dom_def->b_tot, Q+1 );
   VECTOR_FLT_SetSize( dom_def->e_tot, Q+1 );
   VECTOR_FLT_SetSize( dom_def->m_occ, Q+1 );

   VEC_X( dom_def->b_tot, 0 ) = 0.0f;
   VEC_X( dom_def->e_tot, 0 ) = 0.0f;
   VEC_X( dom_def->m_occ, 0 ) = 0.0f;

   /* compute posterior for B, E, and core model states (via HMMER method) */
   if (true)
   {
      for (q_0 = 1; q_0 <= Q; q_0++) 
      {
         q_1 = q_0 - 1;

         /* probability of having reached the begin state */
         btot_add  = XMX_X( sp_MX_fwd, SP_B, q_1 ) + 
                     XMX_X( sp_MX_bck, SP_B, q_1 ) - 
                     overall_logp;
         VEC_X( dom_def->b_tot, q_0 ) = VEC_X( dom_def->b_tot, q_1 ) + 
                                        exp(btot_add);

         /* probability of having reached the end state */
         etot_add  = XMX_X( sp_MX_fwd, SP_E, q_0 ) + 
                     XMX_X( sp_MX_bck, SP_E, q_0 ) - 
                     overall_logp;
         VEC_X( dom_def->e_tot, q_0 ) = VEC_X( dom_def->e_tot, q_1 ) + 
                                        exp(etot_add); 

         /* probability of being in the main model (match,insert,delete) */
         np  = XMX_X( sp_MX_fwd, SP_N, q_1 ) + 
               XMX_X( sp_MX_bck, SP_N, q_0 ) + 
               XSC_X( t_prof, SP_N, SP_LOOP ) - 
               overall_logp;
         np  = expf(np);
         jp  = XMX_X( sp_MX_fwd, SP_J, q_1 ) + 
               XMX_X( sp_MX_bck, SP_J, q_0 ) + 
               XSC_X( t_prof, SP_J, SP_LOOP ) - 
               overall_logp;
         jp  = expf(jp);
         cp  = XMX_X( sp_MX_fwd, SP_C, q_1 ) + 
               XMX_X( sp_MX_bck, SP_C, q_0 ) + 
               XSC_X( t_prof, SP_C, SP_LOOP ) - 
               overall_logp;
         cp  = expf(cp);
         njcp = np + jp + cp;
         
         VEC_X( dom_def->m_occ, q_0 ) = 1.0 - njcp;
      }
   }

   #if DEBUG
   {
      idx = VECTOR_FLT_Create();
      VECTOR_FLT_SetSize( idx, Q+1 );
      for (q_0 = 0; q_0 <= Q; q_0++) {
         VEC_X( idx, q_0 ) = q_0;
      } 

      fp = fopen("test_output/my.njcp.exp.csv", "w+");
      VECTOR_FLT_Dump( idx, fp );
      VECTOR_FLT_Dump( dom_def->b_tot, fp );
      VECTOR_FLT_Dump( dom_def->e_tot, fp );
      VECTOR_FLT_Dump( dom_def->m_occ, fp );
      fclose(fp);

      idx = VECTOR_FLT_Destroy(idx);
   }
   #endif

   dom_def->n_regions   = 0;
   dom_def->n_domains   = 0;
   
   /* find domain ranges (HMMER method) */
   if (true)
   {
      VECTOR_RANGE_Reuse( dom_def->dom_ranges );

      /* flag to check if we are currently in a domain */
      is_in_domain   = false;
      q_beg = q_end  = -1;

      for (q_0 = 1; q_0 <= Q; q_0++)
      {
         q_1 = q_0 - 1;
         mocc_0 = VEC_X(dom_def->m_occ, q_0);
         btot_0 = VEC_X(dom_def->b_tot, q_0);
         btot_1 = VEC_X(dom_def->b_tot, q_1);
         etot_0 = VEC_X(dom_def->e_tot, q_0);
         etot_1 = VEC_X(dom_def->e_tot, q_1);
         /* total probability in main model */
         rt1_test    = mocc_0;
         /* check if sufficient probability we have reached begin state for domain at q_0th position */
         rt2_btest   = mocc_0 - (btot_0 - btot_1);
         /* check if sufficient probability we have reached end state for domain at q_0th position */
         rt2_etest   = mocc_0 - (etot_0 - etot_1);

         #if DEBUG
         {
            // if ( (q_0 >= 800 && q_0 < 850) || (q_0 < 50) ) {
            //    printf("RT_SCOR[%d]: %11.4f %11.4f, %11.4f\n", q_0, rt1_test, rt2_btest, rt2_etest);
            //    printf("RT_TEST[%d]: %11d %11d, %11d\n", q_0, rt1_test >= rt1, rt2_btest < rt2, rt2_etest < rt2);
            // }
         }
         #endif

         /* if start point of new domain has not been found */
         if (is_in_domain == false)
         {
            /* check if we have found a new begin state */
            if ( rt2_btest < rt2 ) {
               q_beg = q_0;
            } 
            /* check if this is the first position in region, init begin state position */
            else if ( q_beg == -1 ) {
               q_beg = q_0;
            }
            /* check if we are in a new domain region */
            if ( rt1_test >= rt1 ) {
               is_in_domain = true;
            }
         }
         /* check if we have reached an end point, creating a valid domain range */
         else if ( rt2_etest < rt2 )
         {
            /*! TODO: Determine whether multiple domains in envelope */
            is_multiple_domains = false;

            /* if single domain */
            if ( is_multiple_domains == false )
            {
               q_end = q_0;
               VECTOR_RANGE_Pushback( dom_def->dom_ranges, (RANGE){q_beg, q_end} );
               dom_def->n_regions++;
               dom_def->n_domains++;

               printf("ADDED REGION: (%d,%d)\n", q_beg, q_end);

               /* end of domain means we are no longer in current domain range */ 
               q_beg = -1;
               is_in_domain = FALSE;
            }
            else /* if multiple domains */
            {
               fprintf(stderr, "ERROR: Multi-domain regions currently not supported.\n");
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
      }
   }

   /* best single region (? WIP) */
   if (false)
   {
      VECTOR_RANGE_Reuse( dom_def->dom_ranges );

      printf("==> DOMAIN SEARCH (best)\n");
      /* flag to check if we are in a current domain */
      is_in_domain = FALSE;
      t_beg = -1;
      t_end = -1;
      rt1_max = -INF;
      rt2_bmax = INF;
      rt2_emax = INF;
      rt2_crit = 0.8;
      rt2_bhit = false;
      rt2_ehit = false;

      for (q_0 = 1; q_0 <= Q; q_0++)
      {
         q_1 = q_0 - 1;
         mocc_0 = VEC_X(dom_def->m_occ, q_0);
         btot_0 = VEC_X(dom_def->b_tot, q_0);
         btot_1 = VEC_X(dom_def->b_tot, q_1);
         etot_0 = VEC_X(dom_def->e_tot, q_0);
         etot_1 = VEC_X(dom_def->e_tot, q_1);
         /* total probability in main model */
         rt1_test    = mocc_0;
         /* check if sufficient probability we have reached begin state for domain at q_0th position */
         rt2_btest   = btot_0;
         /* check if sufficient probability we have reached end state for domain at q_0th position */
         rt2_etest   = etot_0;

         #if DEBUG
         {
            // if ( (q_0 >= 800 && q_0 < 850) || (q_0 < 50) ) {
            //    printf("RT_SCOR[%d]: %11.4f %11.4f, %11.4f\n", q_0, rt1_test, rt2_btest, rt2_etest);
            //    printf("RT_TEST[%d]: %11d %11d, %11d\n", q_0, rt1_test >= rt1, rt2_btest < rt2, rt2_etest < rt2);
            // }
         }
         #endif

         /* find e state which contributes greatest to final score */
         if ( rt2_bhit == false && rt2_btest > rt2_crit ) {
            q_beg = q_0;
            rt2_bhit = true;
         } 
         /* find b state which contributes greatest to final score */
         if ( rt2_ehit == false && rt2_etest > rt2_crit )
         {
            q_end = q_0;
            rt2_ehit = true;
         }
      }
      printf("BEST REGION: (%d,%d)\n", q_beg, q_end);
      VECTOR_RANGE_Pushback( dom_def->dom_ranges, (RANGE){q_beg, q_end} );
   }

   return STATUS_SUCCESS;
}



