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
#include "../algs_sparse/algs_sparse.h"
#include "../algs_linear/algs_linear.h"
#include "../parsers/parsers.h"

/* header */
#include "posterior_sparse.h"

/*! FUNCTION:  run_Posterior_Sparse()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward <st_MX_bck>.
 *             Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in supplied <st_MX_post> (can override input matrices).
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Posterior_Sparse(   SEQUENCE*               q_seq,            /* query sequence */
                        HMM_PROFILE*            t_prof,           /* target hmm model */
                        int                     Q,                /* query length */
                        int                     T,                /* target length */
                        HMM_BG*                 bg,               /* hmm background model */
                        EDGEBOUNDS*             edg,              /* edgebounds */
                        MATRIX_3D_SPARSE*       st_SMX_fwd,       /* normal state matrix for forward */
                        MATRIX_2D*              sp_MX_fwd,        /* special state matrix for forward */
                        MATRIX_3D_SPARSE*       st_SMX_bck,       /* normal state matrix for backward */
                        MATRIX_2D*              sp_MX_bck,        /* special state matrix for backward */
                        MATRIX_3D_SPARSE*       st_SMX_post,      /* OUTPUT: normal state matrix for posterior */
                        MATRIX_2D*              sp_MX_post,       /* OUTPUT: special state matrix for posterior */     
                        DOMAIN_DEF*             dom_def )         /* OUTPUT: domain data */
{
   printf("=== POSTERIOR HEURISTICS (sparse) ===\n");
   printf("==> cutoffs: rt1=%6.3f, rt2=%6.3f, rt3=%6.3f\n",
      dom_def->rt1, dom_def->rt2, dom_def->rt3 );
   
   FILE*    fp;
   RANGE    Q_range, T_range;
   int      Q_size, T_size;
   float    sc, sc1, sc2;
   
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
   T_range.end = MIN(T_range.end, T + 1);
   Q_range.beg = MAX(Q_range.beg, 0);
   Q_range.end = MIN(Q_range.end, Q + 1);
   /* resize matrix to cover bounding box */
   Q_size = Q_range.end - Q_range.beg;
   T_size = T_range.end - T_range.beg;

   #if DEBUG 
   {
      fp = fopen("test_output/my.sparse_bck.in.mx", "w+");
      MATRIX_3D_SPARSE_Embed(Q, T, st_SMX_bck, debugger->test_MX);
      DP_MATRIX_Dump(Q, T, debugger->test_MX, sp_MX_bck, fp);
      fclose(fp);
   }
   #endif

   /* compute Posterior (forward x backward) */
   fprintf(stdout, "# ==> Posterior\n");
   run_Decode_Posterior_Sparse( q_seq, t_prof, Q, T, st_SMX_fwd->edg_inner,
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
   
   /* run Null2 Score to compute Composition Bias */
   fprintf(stdout, "# ==> Null2 Compo Bias\n");
   run_Null2_ByExpectation_Sparse( q_seq, t_prof, Q, T, &Q_range, &T_range, st_SMX_post->edg_inner,
         st_SMX_post, sp_MX_post, dom_def );

   fprintf(stdout, "# ==> Posterior (end)\n");

   return STATUS_SUCCESS;
}


/*! FUNCTION:  run_Decode_Normal_Posterior_Sparse()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store output <...post> matrix in input <...bck> matrix.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Posterior_Sparse(  SEQUENCE*            q_seq,            /* query sequence */
                              HMM_PROFILE*         t_prof,           /* target hmm model */
                              int                  Q,                /* query length */
                              int                  T,                /* target length */
                              EDGEBOUNDS*          edg,              /* edgebounds */
                              MATRIX_3D_SPARSE*    st_SMX_fwd,       /* normal state matrix for forward */
                              MATRIX_2D*           sp_MX_fwd,        /* special state matrix for forward */
                              MATRIX_3D_SPARSE*    st_SMX_bck,       /* normal state matrix for backward */
                              MATRIX_2D*           sp_MX_bck,        /* special state matrix for backward */
                              MATRIX_3D_SPARSE*    st_SMX_post,      /* OUTPUT: normal state matrix for posterior */
                              MATRIX_2D*           sp_MX_post )      /* OUTPUT: normal state matrix for posterior */
{
   FILE* fp = NULL;

   printf("=== run_Decode_Posterior_Sparse ===\n");
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
   int      id_0, lb_0, rb_0;
   int      lb_T, rb_T;
   int      len;
   /* overall score */
   float    overall_sc;
   /* common scale factor denominator */
   float    denom;
   /* temp mx scores */
   float    mmx, imx, dmx, smx;
   float    mmx_, imx_, dmx_, smx_;

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   overall_sc  =  XMX_X( sp_MX_fwd, SP_C, Q ) + 
                  XSC_X(t_prof, SP_C, SP_MOVE);
   printf("overall_sc: %f %f ==> %f\n", 
      XMX_X( sp_MX_fwd, SP_C, Q ), XSC_X(t_prof, SP_C, SP_MOVE), overall_sc);

   /* COMPUTE POSTERIOR */
   if (true)
   {
      /* initialize special states for zero row */
      q_0 = 0;
      XMX_X(sp_MX_post, SP_E, q_0) = -INF;
      XMX_X(sp_MX_post, SP_N, q_0) = -INF;
      XMX_X(sp_MX_post, SP_J, q_0) = -INF;
      XMX_X(sp_MX_post, SP_B, q_0) = -INF;
      XMX_X(sp_MX_post, SP_C, q_0) = -INF;

      /* pass over preceding edgebounds from list */
      q_0  = 0;
      r_0b = r_0 = 0;  /* beginning index for current row in list */
      /* add every inner edgebound from current row */
      r_0b = r_0;
      while ( (r_0 < edg->N) && (EDG_X(edg, r_0).id < q_0) ) {
         r_0++;
      }
      r_0e = r_0;

      /* process zero row edgebounds from list */
      r_0b = r_0;
      while ( (r_0 < edg->N) && (EDG_X(edg, r_0).id == q_0) ) {
         r_0++;
      }
      r_0e = r_0;

      /* FOR every BOUND in zero ROW */
      for (r_0 = r_0b; r_0 < r_0e; r_0++)
      {
         /* get bound data */
         bnd   = &EDG_X(edg, r_0);
         id_0  = bnd->id; 
         lb_T  = (bnd->lb <= 0);
         lb_0  = MAX(1, bnd->lb);         /* can't overflow the left edge */
         rb_T  = (bnd->rb >= T+1);
         rb_0  = MIN(bnd->rb, T+1);       /* can't overflow the right edge */

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

      /* FOR every position in QUERY sequence (row in matrix) */
      for (q_0 = 1; q_0 <= Q; q_0++)
      {
         q_1 = q_0 - 1;
         t_0 = 0;

         /* add every inner edgebound from current row */
         r_0b = r_0;
         while ( (r_0 < edg->N) && (EDG_X(edg, r_0).id == q_0) ) {
            r_0++;
         }
         r_0e = r_0;

         /* FOR every BOUND in current ROW */
         for (r_0 = r_0b; r_0 < r_0e; r_0++)
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0);
            id_0  = bnd->id; 
            lb_T  = (bnd->lb <= 0);
            lb_0  = MAX(1, bnd->lb);         /* can't overflow the left edge */
            rb_T  = (bnd->rb >= T+1);
            rb_0  = MIN(bnd->rb, T+1);       /* can't overflow the right edge */
            len   = rb_0 - lb_0;

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0 = VECTOR_INT_Get( st_SMX_fwd->imap_cur, r_0 );    /* (q_0, t_0) location offset */
            qx1 = VECTOR_INT_Get( st_SMX_fwd->imap_prv, r_0 );    /* (q_1, t_0) location offset */

            /* unrolled first loop: special case for left edge of range */
            if (lb_T)
            {
               t_0 = bnd->lb;
               tx0 = t_0 - bnd->lb;
               /* zero column is -inf in logspace.  We can skip this step and convert to normal space now. */
               MSMX_X(st_SMX_post, qx0, tx0) = -INF;
               ISMX_X(st_SMX_post, qx0, tx0) = -INF;
               DSMX_X(st_SMX_post, qx0, tx0) = -INF; 
            }

            /* MAIN RECURSION */
            /* FOR every position in TARGET profile */
            for (t_0 = lb_0; t_0 < rb_0 - 1; t_0++)
            {
               t_1 = t_0 - 1; 
               tx0 = t_0 - bnd->lb;
               tx1 = tx0 - 1;

               /* normal states */
               /* logprod handles the case in which both are  */
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
               t_1 = t_0 - 1;
               tx0 = t_0 - bnd->lb;
               tx1 = tx0 - 1;

               /* normal states */
               mmx = MSMX_X(st_SMX_fwd, qx0, tx0) + 
                     MSMX_X(st_SMX_bck, qx0, tx0) -
                     overall_sc;
               MSMX_X(st_SMX_post, qx0, tx0) = mmx;

               ISMX_X(st_SMX_post, qx0, tx0) = -INF;
               DSMX_X(st_SMX_post, qx0, tx0) = -INF;
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
   if (true)
   {
      MATRIX_3D_SPARSE_Exp( st_SMX_post );
      MATRIX_2D_Exp( sp_MX_post );
   }
   

   /* NORMALIZE MATRIX */
   if (false)
   {
      /* pass over preceding edgebounds from list */
      q_0  = 0;
      r_0b = r_0 = 0;  /* beginning index for current row in list */
      /* add every inner edgebound from current row */
      r_0b = r_0;
      while ( (r_0 < edg->N) && (EDG_X(edg, r_0).id < q_0) ) {
         r_0++;
      }
      r_0e = r_0;

      /* process zero row edgebounds from list */
      r_0b = r_0;
      while ( (r_0 < edg->N) && (EDG_X(edg, r_0).id == q_0) ) {
         r_0++;
      }
      r_0e = r_0;

      /* FOR every position in QUERY sequence (row in matrix) */
      for (q_0 = 1; q_0 <= Q; q_0++)
      {
         q_1 = q_0 - 1;
         t_0 = 0;

         denom = 0.0;

         /* add every inner edgebound from current row */
         r_0b = r_0;
         while ( (r_0 < edg->N) && (EDG_X(edg, r_0).id == q_0) ) {
            r_0++;
         }
         r_0e = r_0;

         /* FOR every BOUND in current ROW */
         for (r_0 = r_0b; r_0 < r_0e; r_0++)
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0);
            id_0  = bnd->id; 
            lb_T  = (bnd->lb <= 0);
            lb_0  = MAX(1, bnd->lb);         /* can't overflow the left edge */
            rb_T  = (bnd->rb >= T+1);
            rb_0  = MIN(bnd->rb, T+1);       /* can't overflow the right edge */

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0 = VECTOR_INT_Get( st_SMX_fwd->imap_cur, r_0 );    /* (q_0, t_0) location offset */
            qx1 = VECTOR_INT_Get( st_SMX_fwd->imap_prv, r_0 );    /* (q_1, t_0) location offset */

            /* MAIN RECURSION */
            /* FOR every position in TARGET profile */
            for (t_0 = lb_0; t_0 < rb_0 - 1; t_0++)
            {
               t_1 = t_0 - 1; 
               tx0 = t_0 - bnd->lb;
               tx1 = tx0 - 1;

               /* normal states */
               denom += MSMX_X(st_SMX_post, qx0, tx0);
               denom += ISMX_X(st_SMX_post, qx0, tx0);
            }

            /* unrolled final loop: special case for right edge of range */
            if ( true && (rb_0 - lb_0 > 1) )  
            {
               denom += MSMX_X(st_SMX_post, qx0, tx0);
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
         /* FOR every BOUND in current ROW */
         for (r_0 = r_0b; r_0 < r_0e; r_0++)
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0);
            id_0  = bnd->id; 
            lb_T  = (bnd->lb <= 0);
            lb_0  = MAX(1, bnd->lb);         /* can't overflow the left edge */
            rb_T  = (bnd->rb >= T+1);
            rb_0  = MIN(bnd->rb, T+1);       /* can't overflow the right edge */

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0 = VECTOR_INT_Get( st_SMX_fwd->imap_cur, r_0 );    /* (q_0, t_0) location offset */
            qx1 = VECTOR_INT_Get( st_SMX_fwd->imap_prv, r_0 );    /* (q_1, t_0) location offset */

            /* if left bound touches the left edge */
            if (lb_T)
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
            for (t_0 = lb_0; t_0 < rb_0 - 1; t_0++)
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
            
            /* normalize by scaling row by common factor denominator */
            XMX_X(sp_MX_post, SP_N, q_0) *= denom; 
            XMX_X(sp_MX_post, SP_J, q_0) *= denom; 
            XMX_X(sp_MX_post, SP_C, q_0) *= denom;
         }
      }
   }
   
   return STATUS_SUCCESS;
}


/*! FUNCTION:  run_Decode_Special_Posterior_Sparse()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store matrix in <...fwd> or <...bck>.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Special_Posterior_Sparse(   SEQUENCE*         q_seq,            /* query sequence */
                                       HMM_PROFILE*      t_prof,           /* target hmm model */
                                       int               Q,                /* query length */
                                       int               T,                /* target length */
                                       MATRIX_2D*        sp_MX_fwd,        /* special state matrix for forward */
                                       MATRIX_2D*        sp_MX_bck,        /* special state matrix for backward */
                                       MATRIX_2D*        sp_MX_post )      /* OUTPUT: special state matrix for posterior */
{
   printf("=== run_Decode_Special_Posterior_Sparse ===\n");
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

   return STATUS_SUCCESS;
}


/*! FUNCTION:  run_Null2_ByExpectation_Sparse()
 *  SYNOPSIS:  Modeled after HMMER p7_GNull2_ByExpectation().
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Null2_ByExpectation_Sparse(  SEQUENCE*            query,            /* query sequence */
                                 HMM_PROFILE*         target,           /* target hmm model */
                                 int                  Q,                /* query length */
                                 int                  T,                /* target length */
                                 RANGE*               Q_range,          /* query span of bounds */
                                 RANGE*               T_range,          /* target span of bounds */
                                 EDGEBOUNDS*          edg,              /* edgebounds */
                                 MATRIX_3D_SPARSE*    st_SMX_post,      /* posterior normal matrix */
                                 MATRIX_2D*           sp_MX_post,       /* posterior special matrix */
                                 DOMAIN_DEF*          dom_def )         /* OUTPUT: domain def's null2_sc vector */
{
   printf("=== run_Null2_ByExpectation (sparse) ===\n");
   FILE*    fp;     
   int      Q_beg, Q_end, Q_len;
   int      T_beg, T_end, T_len;
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* maps column index into data index (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */
   int      tx0, tx1;                        /* maps target index into data index (target)  */
   int      st_0, k_0;                       /* state and amino acid index */
   int      r_0, r_0b, r_0e;                 /* edgebound list index, beginning and end */
   int      id_0, lb_0, rb_0, rb_T;            /* edgebound range indexs */
   float    mmx, imx, dmx;                   /* {MID} state */

   BOUND*   bnd;
   float    x_factor;
   float    null2sc;
   float    neglog_Q, neglog_cnt;                         

   Q_beg = Q_range->beg;
   Q_end = Q_range->end;
   // T_beg = MAX(1, T_range->beg);
   T_beg = T_range->beg;

   T_end = T_range->end;
   Q_len = Q_end - Q_beg;
   T_len = T_end - T_beg;

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

   /* pass over edgebounds preceeding zero row */
   r_0 = r_0b = r_0e = 0;
   r_0b = r_0;
   while ( (r_0 < edg->N) && (EDG_X(edg, r_0).id < q_0) ) {
      r_0++;
   }
   r_0e = r_0;

   /* sum over each position in target model into vectors  */
   /* FOR every position in QUERY sequence (row in matrix) */
   for (q_0 = 0; q_0 <= Q; q_0++)
   {
      q_1 = q_0 - 1;

      /* add every inner edgebound from current row */
      r_0b = r_0;
      while ( (r_0 < edg->N) && (EDG_X(edg, r_0).id == q_0) ) {
         r_0++;
      }
      r_0e = r_0;

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
            t_1 = t_0 - 1; 
            tx0 = t_0 - bnd->lb;
            tx1 = tx0 - 1;

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

      /* special state emissions */
      for ( st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
         VEC_X( dom_def->sp_freq, st_0 ) += MX_2D( sp_MX_post, st_0, q_0 );
      }
   }

   #if DEBUG
   {
      fp = fopen("test_output/my.post_vec.2.sparse.csv", "w+");
      fprintf(fp, "<2>\n");
      fprintf(fp, "==> NORMAL STATES (st_freq) <=\n");
      for (int t_0 = T_beg; t_0 < T_end; t_0++) {
         tx0 = t_0 - T_beg;
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

   /* convert probabilities to log frequencies */
   /* for each position in query domain */
   for ( t_0 = 0; t_0 <= T; t_0++ ) {
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
      fp = fopen("test_output/my.post_vec.3.sparse.csv", "w+");
      fprintf(fp, "<3>\n");
      fprintf(fp, "==> NORMAL STATES (st_freq) <=\n");
      for (int t_0 = T_beg; t_0 < T_end; t_0++) {
         tx0 = t_0 - T_beg;
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
   printf("neglog_Q = %d %f\n", Q_len, neglog_Q);
   // for ( t_0 = T_beg; t_0 < T_end; t_0++ ) {
   //    tx0 = t_0 - T_beg;
   //    VEC_X( dom_def->st_num, t_0 ) = -log( (float)VEC_X( dom_def->st_num, t_0 ));
   // }

   /* for each position in query domain */
   for ( t_0 = T_beg; t_0 < T_end; t_0++ ) {
      tx0 = t_0 - T_beg;
      neglog_cnt = VEC_X( dom_def->st_num, t_0 );
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
      for (int t_0 = T_beg; t_0 < T_end; t_0++) {
         tx0 = t_0 - T_beg;
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
      fprintf(stdout, "xfactor: %9.4f,  NCJ: %9.4f %9.4f %9.4f\n", 
         x_factor, VEC_X( dom_def->sp_freq, SP_N), VEC_X(dom_def->sp_freq, SP_C), VEC_X(dom_def->sp_freq, SP_J) );
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
      fp = fopen("test_output/my.null2sc.sparse.csv", "w+");
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

