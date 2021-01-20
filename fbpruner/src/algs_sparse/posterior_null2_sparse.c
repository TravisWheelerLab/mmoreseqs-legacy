/*******************************************************************************
 *  FILE:      posterior_null2_sparse.c
 *  PURPOSE:   Computer the Null2 Composition Bias of the Posterior Probability.
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

/* local imports */
#include "../objects/structs.h"
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../algs_linear/algs_linear.h"
#include "../parsers/parsers.h"

/* header */
#include "_algs_sparse.h"
#include "posterior_null2_sparse.h"

/*! FUNCTION:  run_Null2_ByExpectation_Sparse()
 *  SYNOPSIS:  Modeled after HMMER p7_GNull2_ByExpectation().
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
STATUS_FLAG
run_Null2_ByExpectation_Sparse(  SEQUENCE*            query,            /* query sequence */
                                 HMM_PROFILE*         target,           /* target hmm model */
                                 int                  Q,                /* query length */
                                 int                  T,                /* target length */
                                 EDGEBOUNDS*          edg,              /* edgebounds */
                                 RANGE*               q_range,          /* OPTIONAL: query range of edgebound cells */
                                 RANGE*               t_range,          /* OPTIONAL: target range of edgebound cells */
                                 RANGE*               dom_range,        /* OPTIONAL: domain range */
                                 MATRIX_3D_SPARSE*    st_SMX_post,      /* posterior normal matrix */
                                 MATRIX_2D*           sp_MX_post,       /* posterior special matrix */
                                 DOMAIN_DEF*          dom_def,          /* OUTPUT: domain def's null2_sc vector */
                                 float*               compo_bias )      /* OUTPUT: Null2 composition bias */
{
   printf("=== run_Null2_ByExpectation_Sparse() [BEGIN] ===\n");
   FILE*    fp;     
   
   float    bias;
   int      Q_beg, Q_end, Q_len;
   int      T_beg, T_end, T_len;
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* maps column index into data index (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */
   int      tx0, tx1;                        /* maps target index into data index (target)  */
   int      st_0, k_0;                       /* state and amino acid index */
   int      r_0, r_0b, r_0e;                 /* edgebound list index, beginning and end */
   int      id_0, lb_0, rb_0, rb_T;          /* edgebound range indexs */
   float    mmx, imx, dmx;                   /* {MID} state */
   RANGE    Q_range;
   RANGE    T_range;
   bool     is_q_0_in_dom_range;

   BOUND*   bnd;
   float    x_factor;
   float    null2sc;
   float    neglog_Q, neglog_cnt;                         

   /* query sequence */
   Q_range.beg = 0;
   Q_range.end = Q;
   if ( q_range != NULL ) {
      Q_range = *q_range;
   }
   if ( dom_range != NULL ) {
      Q_range = *dom_range;
   }
   /* target range */
   T_range.beg = 0;
   T_range.end = T + 1;
   if ( t_range != NULL ) {
      T_range = *t_range;
   }

   Q_len = Q_range.end - Q_range.beg + 1;
   T_len = T_range.end - T_range.beg + 1;

   /* TODO: temp var */
   const int NULL2_LEN = NUM_AMINO_PLUS_SPEC;
   
   MATRIX_2D_Reuse( dom_def->st_freq, T+1, NUM_NORMAL_STATES );
   VECTOR_FLT_SetSize( dom_def->st_num, T+1 );
   VECTOR_FLT_SetSize( dom_def->sp_freq,  NUM_SPECIAL_STATES );
   VECTOR_FLT_SetSize( dom_def->null2_sc, NUM_AMINO_PLUS_SPEC );
   VECTOR_FLT_SetSize( dom_def->null2_exp, Q+1 );

   MATRIX_2D_Fill( dom_def->st_freq, 0.0f );
   VECTOR_FLT_Fill( dom_def->st_num, 0.0f );
   VECTOR_FLT_Fill( dom_def->sp_freq, 0.0f );
   VECTOR_FLT_Fill( dom_def->null2_sc, 0.0f );
   VECTOR_FLT_Fill( dom_def->null2_exp, 0.0f );

   /* init indexes */
   r_0b = r_0e = 0;

   /* sum over each position in target model into vectors  */
   /* FOR every position in QUERY sequence (row in matrix) */
   for (q_0 = 0; q_0 <= Q; q_0++)
   {
      q_1 = q_0 - 1;
      /* get edgebound range */
      EDGEBOUNDS_NxtRow( edg, &r_0b, &r_0e, q_0 );
      /* check if query position is in domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );

      if ( is_q_0_in_dom_range == true )
      {
         /* FOR every BOUND in current ROW */
         for (r_0 = r_0b; r_0 < r_0e; r_0++)
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0);
            id_0  = bnd->id;
            lb_0  = bnd->lb;
            rb_0  = bnd->rb;

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0 = VECTOR_INT_Get( st_SMX_post->imap_cur, r_0 );    /* (q_0, t_0) location offset */
            qx1 = VECTOR_INT_Get( st_SMX_post->imap_prv, r_0 );    /* (q_1, t_0) location offset */

            /* initial location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

            /* MAIN RECURSION */
            /* FOR every position in TARGET profile */
            for (t_0 = lb_0; t_0 < rb_0; t_0++)
            {
               tx0 = t_0 - bnd->lb;

               /* normal state emissions */
               mmx = MSMX_X( st_SMX_post, qx0, tx0);
               imx = ISMX_X( st_SMX_post, qx0, tx0);
               dmx = DSMX_X( st_SMX_post, qx0, tx0);

               MX_2D( dom_def->st_freq, t_0, M_ST ) += mmx; 
               MX_2D( dom_def->st_freq, t_0, I_ST ) += imx;
               MX_2D( dom_def->st_freq, t_0, D_ST ) += dmx;

               /* add to count */
               // VEC_X( dom_def->st_num, t_0 ) += 1;
            }
         }
      }

      /* special state emissions */
      for ( st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
         VEC_X( dom_def->sp_freq, st_0 ) += XMX_X( sp_MX_post, st_0, q_0 );
      }
   }

   #if DEBUG
   {
      fp = fopen("test_output/my.post_vec.2.sparse.csv", "w+");
      fprintf(fp, "<2>\n");
      fprintf(fp, "==> NORMAL STATES (st_freq) <=\n");
      for (t_0 = T_range.beg; t_0 < T_range.end; t_0++) {
         tx0 = t_0 - T_range.beg;
         fprintf(fp, "ST[%3d]: %12.7f %12.7f %12.7f\n", t_0, 
            MX_2D( dom_def->st_freq, t_0, M_ST ), 
            MX_2D( dom_def->st_freq, t_0, I_ST ), 
            MX_2D( dom_def->st_freq, t_0, D_ST ) );
      }
      fprintf(fp, "==> SPECIAL STATES (sp_freq) <=\n");
      for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
         fprintf(fp, "SP[%d]: %12.7f\n", st_0, 
            VEC_X( dom_def->sp_freq, st_0 ) );
      }
      fclose(fp);
   }
   #endif 

   /* convert probabilities to log frequencies */
   /* for each position in query domain */
   for ( t_0 = 0; t_0 <= T; t_0++ ) {
      tx0 = t_0 - T_range.beg;
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
      fp = fopen("test_output/my.post_vec.3.sparse.csv", "w+");
      fprintf(fp, "<3>\n");
      fprintf(fp, "==> NORMAL STATES (st_freq) <=\n");
      for (int t_0 = T_range.beg; t_0 < T_range.end; t_0++) {
         tx0 = t_0 - T_range.beg;
         fprintf(fp, "ST[%3d]: %12.7f %12.7f %12.7f\n", t_0, 
            MX_2D( dom_def->st_freq, tx0, MAT_ST ), 
            MX_2D( dom_def->st_freq, tx0, INS_ST ), 
            MX_2D( dom_def->st_freq, tx0, DEL_ST ) );
      }
      fprintf(fp, "==> SPECIAL STATES (sp_freq) <=\n");
      for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
         fprintf(fp, "SP[%d]: %12.7f\n", st_0, 
            VEC_X( dom_def->sp_freq, st_0 ) );
      }
      fclose(fp);
   }
   #endif 

   /*  This is the cumulative score contributed by each position in the model.
    *  Divide by Q to take the average per cell (or subtract by the log(Q) to compute in log space) 
    */
   neglog_Q = -log( (float)Q );
   // neglog_Q = -log( (float)Q_len );
   // printf("neglog_Q = %d %f\n", Q_len, neglog_Q);
   // for ( t_0 = T_beg; t_0 < T_end; t_0++ ) {
   //    tx0 = t_0 - T_beg;
   //    VEC_X( dom_def->st_num, t_0 ) = -log( (float)VEC_X( dom_def->st_num, t_0 ));
   // }

   /* for each position in query domain */
   for ( t_0 = T_range.beg; t_0 < T_range.end; t_0++ ) {
      tx0 = t_0 - T_range.beg;
      // neglog_cnt = VEC_X( dom_def->st_num, t_0 );

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
      fp = fopen("test_output/my.post_vec.4.sparse.csv", "w+");
      fprintf(fp, "<4>\n");
      fprintf(fp, "==> NORMAL STATES (st_freq) <=\n");
      for (int t_0 = T_range.beg; t_0 < T_range.end; t_0++) {
         tx0 = t_0 - T_range.beg;
         fprintf(fp, "ST[%3d]: %12.7f %12.7f %12.7f\n", t_0, 
            MX_2D( dom_def->st_freq, tx0, MAT_ST ), 
            MX_2D( dom_def->st_freq, tx0, INS_ST ), 
            MX_2D( dom_def->st_freq, tx0, DEL_ST ) );
      }
      fprintf(fp, "==> SPECIAL STATES (sp_freq) <=\n");
      for (int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++) {
         fprintf(fp, "SP[%d]: %12.7f\n", st_0, 
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
      // fprintf(stdout, "xfactor: %9.4f,  NCJ: %9.4f %9.4f %9.4f\n", 
      //    x_factor, VEC_X( dom_def->sp_freq, SP_N), VEC_X(dom_def->sp_freq, SP_C), VEC_X(dom_def->sp_freq, SP_J) );
   }
   #endif

   /* initialize null2 vector for log-space */
   VECTOR_FLT_Fill( dom_def->null2_sc, -INF );

   /* null2 log emissions probabilities found by summing over 
    * all emmissions used in paths explaining the domain. 
    */
   /* for each amino acid */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) 
   {
      /* for each position in model */
      for ( t_0 = T_range.beg; t_0 < T_range.end - 1; t_0++ ) 
      {
         tx0 = t_0 - T_range.beg;
         /* look at the log frequencies (weighted probability of position in model contributed to path score ) 
          * at the model position and multiply them by the score contribution
          */ 
         mmx = MX_2D( dom_def->st_freq, t_0, MAT_ST) + MSC_X(target, t_0, k_0);
         imx = MX_2D( dom_def->st_freq, t_0, INS_ST) + ISC_X(target, t_0, k_0);

         VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ), mmx );
         VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ), imx );
      }
      t_0 = T_range.end - 1;
      tx0 = t_0 - T_range.beg;
      mmx = MX_2D( dom_def->st_freq, t_0, MAT_ST) + MSC_X(target, t_0, k_0);
      
      VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ), mmx);
      VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ), x_factor );                            
   }

   #if DEBUG
   {
      // fp = stdout;
      fp = fopen("test_output/my.null2.sparse.csv", "w+");
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
   for ( q_0 = Q_range.beg; q_0 < Q_range.end; q_0++ ) {
      qx0 = q_0 - Q_range.beg;
      k_0 = AA_REV[query->seq[q_0]];
      float val = VEC_X( dom_def->null2_sc, k_0 );
      VEC_X( dom_def->null2_exp, qx0) = logf( val );
   }

   /* Finally, for bias correction, sum over the entire expected sequence bias. */
   bias = 0.0f;
   for ( q_0 = Q_range.beg; q_0 < Q_range.end; q_0++ ) {
      qx0 = q_0 - Q_range.beg;
      k_0 = AA_REV[query->seq[q_0]];
      bias += logf( VEC_X( dom_def->null2_sc, k_0 ));
   }

   *compo_bias = bias;
   dom_def->seq_bias = bias;

   #if DEBUG
   {
      // fp = stdout;
      fp = fopen("test_output/my.null2sc.sparse.csv", "w+");
      fprintf(fp, "# === NULL2SC (Expectation by Amino) ===\n");
      for (int k_0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++) {
         fprintf(fp, "%d %c %.9f\n", k_0, AA[k_0], VEC_X( dom_def->null2_sc, k_0 ));
      }
      fprintf(fp, "# === NULL2SC (Expectation by Position) ===\n");
      for ( q_0 = Q_range.beg; q_0 < Q_range.end; q_0++ ) {
         qx0 = q_0 - Q_range.beg;
         fprintf(fp, "%d %.9f\n", q_0, VEC_X( dom_def->null2_exp, qx0 ));
      }
      fprintf(fp, "# === DOMCORRECTION: %.9f\n", dom_def->seq_bias);
      if (fp != stdout) fclose(fp);

      // printf("# === DOMCORRECTION: %.9f\n", dom_def->seq_bias);
   }
   #endif

   return STATUS_SUCCESS;
}