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
   
   RANGE Q_range, T_range;
   int   Q_size, T_size;
   float sc, sc1, sc2;
   
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

   /* compute Posterior (forward * backward) */
   fprintf(stdout, "# ==> Posterior\n");
   run_Decode_Posterior_Sparse( q_seq, t_prof, Q, T, edg,
         st_SMX_fwd, sp_MX_fwd, st_SMX_bck, sp_MX_bck, st_SMX_post, sp_MX_post );
   
   #if DEBUG 
   {
      MATRIX_3D_SPARSE_Embed(st_SMX_post, debugger->test_MX);
      DP_MATRIX_Save(Q, T, debugger->test_MX, sp_MX_post, "test_output/sparse_posterior.mx");
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
   fprintf(stdout, "# ==> Null2 Compo Bias\n");
   run_Null2_ByExpectation_Sparse( q_seq, t_prof, Q, T, &Q_range, &T_range, st_SMX_post->edg_inner,
         st_SMX_post, sp_MX_post, dom_def );

   fprintf(stdout, "# ==> Posterior (end)\n");

   return STATUS_SUCCESS;
}


/*! FUNCTION:  run_Decode_Normal_Posterior_Sparse()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store matrix in <...fwd> or <...bck>.
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
   int      r_0, r_0b, r_0e;
   BOUND*   bnd;
   int      id, lb_0, lb_T, rb_0, rb_T;
   /* overall score */
   float    overall_sc;
   /* common scale factor denominator */
   float    denom;
   /* temp mx scores */
   float    mmx, imx, dmx, smx;
   float    mmx_, imx_, dmx_, smx_;

   // printf("==> POSTERIOR :: forward\n");
   // MATRIX_3D_SPARSE_Dump( st_SMX_fwd, stdout );
   // printf("==> POSTERIOR :: backward\n");
   // MATRIX_3D_SPARSE_Dump( st_SMX_bck, stdout );

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

   /* add every inner edgebound from current row */
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
      id    = bnd->id;
      lb_0  = bnd->lb;           /* can't overflow the left edge */
      rb_0  = bnd->rb;           /* can't overflow the right edge */

      /* fetch data mapping bound start location to data block in sparse matrix */
      qx0 = VECTOR_INT_Get( st_SMX_fwd->imap_cur, r_0 );    /* (q_0, t_0) location offset */
      qx1 = VECTOR_INT_Get( st_SMX_fwd->imap_prv, r_0 );    /* (q_1, t_0) location offset */

      /* initial location for square matrix and mapping to sparse matrix */
      t_0 = lb_0;
      tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

      /* FOR every position in TARGET profile */
      for (t_0 = lb_0; t_0 < rb_0; t_0++)
      {
         tx0 = t_0 - bnd->lb;
         MSMX_X(st_SMX_post, qx0, tx0) = 0.0;
         ISMX_X(st_SMX_post, qx0, tx0) = 0.0;
         DSMX_X(st_SMX_post, qx0, tx0) = 0.0; 
      }
   }

   /* MAIN RECURSION */
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
      // printf("q_0 = %d/%d, r_0 = {%d,%d}\n", q_0, Q, r_0b, r_0e);

      /* FOR every BOUND in current ROW */
      for (r_0 = r_0b; r_0 < r_0e; r_0++)
      {
         /* get bound data */
         bnd   = &EDG_X(edg, r_0);
         id    = bnd->id; 
         lb_T  = ( bnd->lb == 0 );
         lb_0  = MAX(1, bnd->lb);         /* can't overflow the left edge */
         rb_T  = ( bnd->rb > T );
         rb_0  = MIN(bnd->rb, T);         /* can't overflow the right edge */
         // printf("r_0 = {%3d:%3d,%3d} = {%3d:%3d,%3d}\n", bnd->id, bnd->lb, bnd->rb, id, lb_0, rb_0);

         /* fetch data mapping bound start location to data block in sparse matrix */
         qx0 = VECTOR_INT_Get( st_SMX_fwd->imap_cur, r_0 );    /* (q_0, t_0) location offset */
         qx1 = VECTOR_INT_Get( st_SMX_fwd->imap_prv, r_0 );    /* (q_1, t_0) location offset */
         // printf("r_0 = %d, qx0 = %d, qx1 = %d\n", r_0, qx0, qx1);

         /* initial location for square matrix and mapping to sparse matrix */
         t_0 = lb_0;
         tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

         /* special case for zero column */
         if ( lb_T ) 
         {
            t_0 = lb_0;
            tx0 = t_0 - bnd->lb;

            MSMX_X(st_SMX_post, qx0, tx0) = 0.0;
            ISMX_X(st_SMX_post, qx0, tx0) = 0.0;
            DSMX_X(st_SMX_post, qx0, tx0) = 0.0; 

            lb_0 += 1;
         }

         /* MAIN RECURSION */
         /* FOR every position in TARGET profile */
         for (t_0 = lb_0; t_0 < rb_0; t_0++)
         {
            t_1 = t_0 - 1; 
            tx0 = t_0 - bnd->lb;
            tx1 = tx0 - 1;
            // printf("t_0 = %d/%d, tx0 = %d ? %d, N = %d\n", t_0, rb_0, tx0, qx0 + tx0, st_SMX_post->data->Nalloc );

            /* normal states */
            mmx = MSMX_X(st_SMX_fwd, qx0, tx0) + 
                  MSMX_X(st_SMX_bck, qx0, tx0) -
                  overall_sc;
            mmx_ = expf(mmx);
            MSMX_X(st_SMX_post, qx0, tx0) = mmx_;
            denom += MSMX_X(st_SMX_post, qx0, tx0);

            imx = ISMX_X(st_SMX_fwd, qx0, tx0) + 
                  ISMX_X(st_SMX_bck, qx0, tx0) -
                  overall_sc;
            imx_ = expf(imx);
            ISMX_X(st_SMX_post, qx0, tx0) = imx_;
            denom += ISMX_X(st_SMX_post, qx0, tx0);

            DSMX_X(st_SMX_post, qx0, tx0) = 0.0;
         }

         /* unrolled loop */
         if ( rb_T )  
         {
            // printf("rb_T\n");
            t_0 = T;
            t_1 = t_0 - 1;
            tx0 = t_0 - bnd->lb;
            tx1 = tx0 - 1;

            /* normal states */
            mmx = MSMX_X(st_SMX_fwd, qx0, tx0) + 
                  MSMX_X(st_SMX_bck, qx0, tx0) -
                  overall_sc;
            mmx_ = expf(mmx);
            MSMX_X(st_SMX_post, qx0, tx0) = mmx_;
            denom += MSMX_X(st_SMX_post, qx0, tx0);

            ISMX_X(st_SMX_post, qx0, tx0) = 0.0;
            DSMX_X(st_SMX_post, qx0, tx0) = 0.0;

            /* embed linear row into quadratic test matrix */
            #if DEBUG
            {
               // MX_2D(cloud_MX, q_0, t_0) = 1.0;
               // MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
               // MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
               // MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
            }
            #endif
         }

      }

      /* special states */
      XMX_X(sp_MX_post, SP_E, q_0) = 0.0;
      XMX_X(sp_MX_post, SP_B, q_0) = 0.0;

      smx = XMX_X(sp_MX_fwd, SP_N, q_1) +
            XMX_X(sp_MX_bck, SP_N, q_0) +
            XSC_X(t_prof, SP_N, SP_LOOP) -
            overall_sc;
      if (q_0 >= Q - 50 && q_0 < Q) {
         printf("N(%d): %9.4f %9.4f %9.4f %9.4f\n", 
            q_0, XMX_X(sp_MX_fwd, SP_N, q_1), XMX_X(sp_MX_bck, SP_N, q_0), smx, expf(smx) );
      }
      XMX_X(sp_MX_post, SP_N, q_0) = expf(smx);

      smx = XMX_X(sp_MX_fwd, SP_J, q_1) +
            XMX_X(sp_MX_bck, SP_J, q_0) +
            XSC_X(t_prof, SP_J, SP_LOOP) -
            overall_sc;
      if (q_0 >= Q - 50 && q_0 < Q) {
         printf("J(%d): %9.4f %9.4f %9.4f %9.4f\n", 
            q_0, XMX_X(sp_MX_fwd, SP_J, q_1), XMX_X(sp_MX_bck, SP_J, q_0), smx, expf(smx) );
      }
      XMX_X(sp_MX_post, SP_J, q_0) = expf(smx);

      smx = XMX_X(sp_MX_fwd, SP_C, q_1) +
            XMX_X(sp_MX_bck, SP_C, q_0) +
            XSC_X(t_prof, SP_C, SP_LOOP) -
            overall_sc;
      if (q_0 >= Q - 50 && q_0 < Q) {
         printf("C(%d): %9.4f %9.4f %9.4f %9.4f\n", 
            q_0, XMX_X(sp_MX_fwd, SP_C, q_1), XMX_X(sp_MX_bck, SP_C, q_0), smx, expf(smx) );
      }
      XMX_X(sp_MX_post, SP_C, q_0) = expf(smx);

      denom += XMX_X(sp_MX_post, SP_N, q_0) + 
               XMX_X(sp_MX_post, SP_J, q_0) + 
               XMX_X(sp_MX_post, SP_C, q_0);

      /* apply denominator scaling factor to entire row */
      /* FOR every BOUND in current ROW */
      for (r_0 = r_0b; r_0 < r_0e; r_0++)
      {
         /* get bound data */
         bnd   = &EDG_X(edg, r_0);
         id    = bnd->id; 
         lb_0  = MAX(1, bnd->lb);         /* can't overflow the left edge */
         rb_T  = ( bnd->rb > T );
         rb_0  = MIN(bnd->rb, T);         /* can't overflow the right edge */

         /* fetch data mapping bound start location to data block in sparse matrix */
         qx0 = VECTOR_INT_Get( st_SMX_fwd->imap_cur, r_0 );    /* (q_0, t_0) location offset */
         qx1 = VECTOR_INT_Get( st_SMX_fwd->imap_prv, r_0 );    /* (q_1, t_0) location offset */

         /* initial location for square matrix and mapping to sparse matrix */
         t_0 = lb_0;
         tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

         /* FOR every position in TARGET profile */
         for (t_0 = lb_0; t_0 < rb_0; t_0++)
         {
            t_1 = t_0 - 1; 
            tx0 = t_0 - bnd->lb;
            tx1 = tx0 - 1;

            // printf("[%2d]: denom :=> %9f\n", q_0, denom);
            /* normalize by scaling row by common factor denominator */
            denom = 1.0 / denom;
            for ( t_0 = 1; t_0 < T; t_0++ ) {
               MSMX_X(st_SMX_post, qx0, tx0) *= denom;
               ISMX_X(st_SMX_post, qx0, tx0) *= denom;
            }
            MSMX_X(st_SMX_post, q_0, T) *= denom;
            XMX_X(sp_MX_post, SP_N, q_0) *= denom; 
            XMX_X(sp_MX_post, SP_J, q_0) *= denom; 
            XMX_X(sp_MX_post, SP_C, q_0) *= denom;
         }
      }
   }

   /** TODO: optimize? */
   /* convert all outer cells from logspace to normal space (-inf -> 0.0) */
   for (int i = 0; i < st_SMX_post->data->N; i++) 
   {
      if (st_SMX_post->data->data[i] == -INF) {
         st_SMX_post->data->data[i] = 0.0f;
      }
   }

   // printf("==> POSTERIOR :: posterior\n");
   // MATRIX_3D_SPARSE_Dump( st_SMX_post, stdout );
   // MATRIX_2D_Dump( sp_MX_post, stdout );

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
   int      Q_beg, Q_end, Q_len;
   int      T_beg, T_end, T_len;
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* maps column index into data index (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */
   int      tx0, tx1;                        /* maps target index into data index (target)  */
   int      st_0, k_0;                       /* state and amino acid index */
   int      r_0, r_0b, r_0e;                 /* edgebound list index, beginning and end */
   int      id, lb_0, rb_0, rb_T;            /* edgebound range indexs */
   BOUND*   bnd;
   float    x_factor;
   float    null2sc;

   Q_len = Q_range->end - Q_range->beg;
   T_len = T_range->end - T_range->beg;
   Q_beg = Q_range->beg;
   Q_end = Q_range->end;
   T_beg = T_range->beg;
   T_end = T_range->end;

   // printf("=== POSTERIOR ===\n");
   // DP_MATRIX_Log_Dump(Q->end, T, st_MX_post, sp_MX_post, stdout);
   // printf("Q,T=(%d,%d)\n", query->N, target->N );
   // printf("=================\n");
   
   VECTOR_FLT_GrowTo( dom_def->st_freq, (T+1) * NUM_NORMAL_STATES );
   VECTOR_FLT_GrowTo( dom_def->sp_freq,  NUM_SPECIAL_STATES );
   VECTOR_FLT_GrowTo( dom_def->null2_sc, NUM_AMINO_PLUS_SPEC );

   /* initialize values */
   /* for each position in query domain */
   for ( t_0 = 0; t_0 <= T; t_0++ ) {
      /* for each normal state emissions */
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0 ) = 0.0;
      }
   }
   /* for each special state emissions */
   for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
      VEC_X( dom_def->sp_freq, st_0 ) = 0.0;
   }

   /* sum over each position in target model into vectors  */
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
         id    = bnd->id;
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
            VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + M_ST ) += MSMX_X( st_SMX_post, qx0, tx0);
            VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + I_ST ) += ISMX_X( st_SMX_post, qx0, tx0);
            VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + D_ST ) += DSMX_X( st_SMX_post, qx0, tx0);
         }
      }

      /* special state emissions */
      for ( st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
         VEC_X( dom_def->sp_freq, st_0 ) += MX_2D( sp_MX_post, st_0, q_0 );
      }
   }

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

   float neglog_Q = -log( (float)Q_len );
   // printf("neglog_Q = %d %f\n", Q_len, neglog_Q);
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

   printf("<4>\n");
   printf("xfactor: %f, NCJ: %f %f %f\n", 
      x_factor, VEC_X( dom_def->sp_freq, SP_N), VEC_X(dom_def->sp_freq, SP_C), VEC_X(dom_def->sp_freq, SP_J) );

   /* x-factor: */
   x_factor = VEC_X( dom_def->sp_freq, SP_N);
   x_factor = logsum( x_factor,
                      VEC_X(dom_def->sp_freq, SP_C) );
   // x_factor = logsum( x_factor,
   //                    VEC_X(dom_def->sp_freq, SP_J) );
   
   printf("<5>\n");
   printf("xfactor: %f, NCJ: %f %f %f\n", 
      x_factor, VEC_X( dom_def->sp_freq, SP_N), VEC_X(dom_def->sp_freq, SP_C), VEC_X(dom_def->sp_freq, SP_J) );

   /* initialize null2 vector with logscore */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = -INF;
   }
   for ( k_0 = NUM_AMINO; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      VEC_X( dom_def->null2_sc, k_0 ) = 0.0;
   }

   printf("<6>\n");
   for (t_0 = 1; t_0 < T; t_0++) {
      printf("st_freq[%d]: %9.4f %9.4f %9.4f\n", 
         t_0, VEC_X(dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST), VEC_X(dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST), VEC_X(dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + DEL_ST) );
   }

   /* temp(?): convert inf and -inf to zero */
   for ( t_0 = 0; t_0 <= T; t_0++ ) {
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         float val = VEC_X(dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0);
         if ( val == INF || val == -INF ) {
             VEC_X(dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + st_0) = 0.0;
         }
      }
      
   }

   /* null2 log emissions probabilities found by summing over 
    * all emmissions used in paths explaining the domain. 
    */
   /* for each amino acid */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) 
   {
      /* for each position in model */
      for ( t_0 = 1; t_0 < T; t_0++ ) 
      {
         if ( VEC_X( dom_def->null2_sc, k_0 ) != -INF ) {
            printf("(k_0,t_0)= (%3d,%3d) MSC= %9.4f, ISC= %9.4f, NULL2= %9.4f\n",
               k_0, t_0, 
               VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST) + MSC_X(target, t_0, k_0),
               VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST) + ISC_X(target, t_0, k_0),
               VEC_X( dom_def->null2_sc, k_0 )
            );
         }

         VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                   VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST) + MSC_X(target, t_0, k_0) );
         VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                   VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + INS_ST) + ISC_X(target, t_0, k_0) );
      }
      t_0 = T;
      VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                VEC_X( dom_def->st_freq, (t_0 * NUM_NORMAL_STATES) + MAT_ST) + MSC_X(target, t_0, k_0) );
      VEC_X( dom_def->null2_sc, k_0 ) = logsum( VEC_X( dom_def->null2_sc, k_0 ),
                                                x_factor );                            
   }

   printf("==> NULL2 (pre) <=\n");
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

   return STATUS_SUCCESS;
}

