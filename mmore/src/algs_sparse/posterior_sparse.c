/*******************************************************************************
 *  FILE:      posterior_sparse.h
 *  PURPOSE:   The Maximum Posterior Probability and Optimal Alignment.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *  NOTES:
 *    - There are a few suboptimal arrangements I have to clarify the posterior process.
 *    - Currently, posterior_decoding converts from log to normal space in a separate middle step.
 *      This could be integrated into the first pass through the matrix and reduce passes.
 *      The third step also takes two passes: once to get the total per row, and again 
 *      to divide by the total denom.  This summing process could also be done during the 
 *      first pass if it were in normal space
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
#include "../algs_linear/_algs_linear.h"
#include "../parsers/_parsers.h"

/* header */
#include "_algs_sparse.h"
#include "posterior_sparse.h"

/*! FUNCTION:  run_Decode_Normal_Posterior_Sparse()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store output <...post> matrix in input <...bck> matrix.
 *             Input <...fwd> and <...bck> matrices must be in log space.
 *             Output <...post> matrices are in normal space.
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
   MATH_Logsum_Init();

   overall_sc  =  XMX_X( sp_MX_fwd, SP_C, Q ) + 
                  XSC_X(t_prof, SP_C, SP_MOVE);
   
   printf("overall_sc: %f %f ==> %f\n", 
      XMX_X( sp_MX_fwd, SP_C, Q ), XSC_X(t_prof, SP_C, SP_MOVE), overall_sc);
   
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

   #if DEBUG
   {
      fp = fopen("test_output/my.raw_post.exp.mx", "w+");
      MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_post, debugger->test_MX);
      DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_post, fp);
      fclose(fp); 
   }
   #endif

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
                  t_1 = t_0 - 1; 
                  tx0 = t_0 - bnd->lb;
                  tx1 = tx0 - 1;
                  
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
                        EDGEBOUNDS*             edg,              /* edgebounds */
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
   
   /* find all domain ranges (HMMER method) */
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
               // printf("ADDED REGION: (%d,%d)\n", q_beg, q_end);

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



