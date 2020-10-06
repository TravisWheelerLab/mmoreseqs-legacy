/*******************************************************************************
 *  FILE:      max_quad.h
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
#include "structs.h"
#include "utilities.h"
#include "objects.h"
#include "algs_quad.h"
#include "algs_linear.h"
#include "parsers.h"

/* header */
#include "posterior_quad.h"

/** FUNCTION:  run_MaxPost_Quad()
 *  SYNOPSIS:  Compute the Max Posterior Probability and obtain optimal alignment.
 *             Computes the Forward and Backward.  Max Posterior computed by Viterbi.
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int 
run_MaxPost_Quad( SEQUENCE*            query,            /* query sequence */
                  HMM_PROFILE*         target,           /* target HMM model */
                  int                  Q,                /* query length */
                  int                  T,                /* target length */
                  MATRIX_3D*           st_MX_fwd,        /* normal matrix for forward */
                  MATRIX_2D*           sp_MX_fwd,        /* special matrix for forward */
                  MATRIX_3D*           st_MX_bck,        /* normal matrix for backward */
                  MATRIX_2D*           sp_MX_bck,        /* special matrix for backward */
                  MATRIX_3D*           st_MX_post,       /* OUTPUT: normal matrix for posterior (can overwrite fwd and bck data) */
                  MATRIX_2D*           sp_MX_post,       /* OUTPUT: special matrix for posterior (can overwrite fwd and bck data) */
                  ALIGNMENT*           aln,              /* OUTPUT: alignment */
                  DOMAIN_DEF*          dom_def )         /* OUTPUT: domain definition scores */
{
   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   /* Posterior */
   run_Posterior_Quad( 
      query, target, Q, T, st_MX_fwd, sp_MX_fwd, st_MX_bck, sp_MX_bck, st_MX_post, sp_MX_post, dom_def );

   /* Viterbi for Posterior */
   // run_MaxPost_Viterbi_Quad( query, target, Q, T, st_MX_post, sp_MX_post, st_MX_max, sp_MX_max, sc_final );
}

/** FUNCTION:  run_Posterior_Quad()
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
                     MATRIX_3D*     st_MX_fwd,        /* normal state matrix for forward */
                     MATRIX_2D*     sp_MX_fwd,        /* special state matrix for forward */
                     MATRIX_3D*     st_MX_bck,        /* normal state matrix for backward */
                     MATRIX_2D*     sp_MX_bck,        /* special state matrix for backward */
                     MATRIX_3D*     st_MX_post,       /* OUTPUT: normal state matrix for posterior */
                     MATRIX_2D*     sp_MX_post,       /* OUTPUT: special state matrix for posterior */     
                     DOMAIN_DEF*    dom_def )         /* OUTPUT: domain data */
{
   printf("=== POSTERIOR HEURISTICS ===\n");
   printf("==> cutoffs: rt1=%6.3f, rt2=%6.3f, rt3=%6.3f\n",
      dom_def->rt1, dom_def->rt2, dom_def->rt3 );

   int      d_0;                          /* domain index */
   int      q_0, q_1;                     /* query index */
   int      t_0, t_1;                     /* target index */
   int      q_beg, q2_beg;                /* start points of domain */
   int      q_end, q2_end, q2_end_last;   /* end points of domain */
   int      q_max;                        /* find center of highest scoring region */
   int      q_size;                       /* size of domain in sequence */
   int      n_cluster;                    /* size of cluster */
   bool     is_triggered;                 /* specifies whether the conditions have been met for a possible alignment domain */
   bool     is_multidomain_region;        /* specifies whether the region contains multiple domains */
   float    curr_sc;                      /* current score for testing domain edge thresholds */
   RANGE    Q_range;
   int      Q_dom;                        /* size of query domain */
   float    max_beg, max_end, max_sc;     /* scores for finding best scoring alignment */
   float    fwd_sc, bck_sc;               /* domain scores */

   /* resize domain definition structure to accept  */
   DOMAIN_DEF_GrowTo( dom_def, Q );

   printf("decode_spec...\n");
   /* compute the posterior probabilties (only need special states) */
   run_Decode_Special_Posterior_Quad( 
      q_seq, t_prof, Q, T, sp_MX_fwd, sp_MX_bck, sp_MX_post );
   printf("decode_domain...\n");
   /* compute domain transition probabilities */
   run_DecodeDomain_Quad(
      q_seq, t_prof, Q, T, sp_MX_fwd, sp_MX_bck, dom_def );

   /* find the best scoring region */
   q_max = -1;
   max_sc = -INF;
   for ( q_0 = 1; q_0 <= Q; q_0++ )
   {
      q_1 = q_0 - 1;
      curr_sc = dom_def->m_occ->data[q_0];
      /* if new max found, update location */
      if ( max_sc < curr_sc ) {
         q_max = q_0;
         max_sc = curr_sc;
      }
   }

   /* find beginning of high scoring region, starting from max score location */
   q_beg = -1;
   max_beg = -INF;
   for ( q_0 = q_max; q_0 >= 1; q_0-- )
   {
      q_1 = q_0 - 1;
      curr_sc = dom_def->m_occ->data[q_0] - ( dom_def->b_tot->data[q_0] - dom_def->b_tot->data[q_1] );
      /* if score falls below cutoff threshold */
      q_beg = q_0;
      if ( curr_sc < dom_def->rt2 ) {
         break;
      }
   }

   /* find end of high scoring region, starting from max score location  */
   q_end = -1;
   max_end = -INF;
   for ( q_0 = q_max; q_0 <= Q; q_0++ )
   {
      q_1 = q_0 - 1;
      curr_sc = dom_def->m_occ->data[q_0] - ( dom_def->e_tot->data[q_0] - dom_def->e_tot->data[q_1] );
      /* if score falls below cutoff threshold */
      q_end = q_0;
      if ( curr_sc < dom_def->rt2 ) {
         break;
      }
   }

   printf("q_range=(%d,%d), q_max=%d\n", q_beg, q_end, q_max );

   /* set sequence to be the domain subrange */
   SEQUENCE_SetSubsequence( q_seq, q_beg, q_end );
   Q_dom = q_end - q_beg + 1;
   /* reconfigure model for domain */
   if ( t_prof->mode == MODE_MULTIGLOCAL|| t_prof->mode == MODE_MULTILOCAL ) {
      HMM_PROFILE_ReconfigMultihit( t_prof, q_seq->N );
   }
   else {
      HMM_PROFILE_ReconfigUnihit( t_prof, q_seq->N );
   }

   /* resize matrices to fit domain */
   MATRIX_3D_Reuse( st_MX_fwd, NUM_NORMAL_STATES, Q_dom+1, T+1 );
   MATRIX_2D_Reuse( sp_MX_fwd, NUM_SPECIAL_STATES, Q_dom+1 );
   MATRIX_3D_Reuse( st_MX_bck, NUM_NORMAL_STATES, Q_dom+1, T+1 );
   MATRIX_2D_Reuse( sp_MX_bck, NUM_SPECIAL_STATES, Q_dom+1 );
   /* NOTE: most likely, posterior dp matrixes will just recycle fwd or bck dp matrices, but do this just in case */
   MATRIX_3D_Reuse( st_MX_post, NUM_NORMAL_STATES, Q_dom+1, T+1 );
   MATRIX_2D_Reuse( sp_MX_post, NUM_SPECIAL_STATES, Q_dom+1 );

   /* capture the posterior */
   run_Forward_Quad( q_seq, t_prof, Q_dom, T, st_MX_fwd, sp_MX_fwd, &fwd_sc );
   run_Backward_Quad( q_seq, t_prof, Q_dom, T, st_MX_bck, sp_MX_bck, &bck_sc );
   run_Decode_Posterior_Quad( q_seq, t_prof, Q_dom, T, st_MX_bck, sp_MX_bck, st_MX_bck, sp_MX_bck, st_MX_post, sp_MX_post );
   /* compute correction bias */
   run_Null2_ByExpectation( q_seq, t_prof, Q_dom, T, st_MX_post, sp_MX_post, dom_def );

   /* set sequence to be full sequence */
   SEQUENCE_UnsetSubsequence( q_seq );
   /* reconfigure model to its state from beginning of function */
   if ( t_prof->mode == MODE_MULTIGLOCAL|| t_prof->mode == MODE_MULTILOCAL ) {
      HMM_PROFILE_ReconfigMultihit( t_prof, q_seq->N );
   }
   else {
      HMM_PROFILE_ReconfigUnihit( t_prof, q_seq->N );
   }
   return STATUS_SUCCESS;
}

/** FUNCTION:  run_Posterior_Quad2()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward <st_MX_bck>.
 *             Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in supplied <st_MX_post> (can override input matrices).
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Posterior_Quad2( SEQUENCE*      q_seq,            /* query sequence */
                     HMM_PROFILE*   t_prof,           /* target hmm model */
                     int            Q,                /* query length */
                     int            T,                /* target length */
                     MATRIX_3D*     st_MX_fwd,        /* normal state matrix for forward */
                     MATRIX_2D*     sp_MX_fwd,        /* special state matrix for forward */
                     MATRIX_3D*     st_MX_bck,        /* normal state matrix for backward */
                     MATRIX_2D*     sp_MX_bck,        /* special state matrix for backward */
                     MATRIX_3D*     st_MX_post,       /* OUTPUT: normal state matrix for posterior */
                     MATRIX_2D*     sp_MX_post,       /* OUTPUT: special state matrix for posterior */     
                     ALIGNMENT*     aln,
                     DOMAIN_DEF*    dom_def )         /* OUTPUT: domain data */
{
   printf("=== POSTERIOR HEURISTICS ===\n");
   printf("==> cutoffs: rt1=%6.3f, rt2=%6.3f, rt3=%6.3f\n",
      dom_def->rt1, dom_def->rt2, dom_def->rt3 );

   int      d_0;                          /* domain index */
   int      q_0, q_1;                     /* query index */
   int      t_0, t_1;                     /* target index */
   int      q_beg, q2_beg;                /* start points of domain */
   int      q_end, q2_end, q2_end_last;   /* end points of domain */
   int      q_size;                       /* size of domain in sequence */
   int      n_cluster;                    /* size of cluster */
   bool     is_triggered;                 /* specifies whether the conditions have been met for a possible alignment domain */
   bool     is_multidomain_region;        /* specifies whether the region contains multiple domains */
   float    threshold_sc;                 /* threshold for testing domain edges */
   RANGE    Q_range;

   /* resize domain definition structure to accept  */
   DOMAIN_DEF_GrowTo( dom_def, Q );

   printf("decode_spec...\n");
   /* compute the posterior probabilties */
   // run_Decode_Normal_Posterior_Quad( 
   //    q_seq, t_prof, Q, T, st_MX_fwd, st_MX_bck, st_MX_post );
   // run_Decode_Special_Posterior_Quad( 
   //    q_seq, t_prof, Q, T, sp_MX_fwd, sp_MX_bck, sp_MX_post );
   printf("decode_domain...\n");
   /* compute domain transition probabilities */
   run_DecodeDomain_Quad(
      q_seq, t_prof, Q, T, sp_MX_fwd, sp_MX_bck, dom_def );

   printf("init null score...\n");
   /* initialize null2 scores */
   for ( q_0 = 0; q_0 < Q; q_0++ ) {
      dom_def->null2_sc->data[q_0] = 0.0;
   }

   /* iterate through domain data by index in sequence */
   for ( q_0 = 1; q_0 <= Q; q_0++ )
   {
      q_1 = q_0 - 1;

      if ( is_triggered == false ) {
         threshold_sc = dom_def->m_occ->data[q_0] - ( dom_def->b_tot->data[q_0] - dom_def->b_tot->data[q_1] );
         printf("BEG(i)=%d :: mocc[j]: %6.3f, btot[j]: %6.3f, btot[j-1]: %6.3f, sc: %6.3f\n",
            q_beg, dom_def->m_occ->data[q_0], dom_def->b_tot->data[q_0], dom_def->b_tot->data[q_1], threshold_sc );
      }
      else {
         threshold_sc = dom_def->m_occ->data[q_0] - ( dom_def->e_tot->data[q_0] - dom_def->e_tot->data[q_1] );
         printf("END(i)=%d :: mocc[j]: %6.3f, etot[j]: %6.3f, etot[j-1]: %6.3f, sc: %6.3f\n",
            q_0, dom_def->m_occ->data[q_0], dom_def->e_tot->data[q_0], dom_def->e_tot->data[q_1], threshold_sc );
      }

      /* if we have not found the start of a potential new domain */
      if ( is_triggered == false ) 
      {
         /* this logic is explained in hmmer code: xref J2/101 */
         if ( threshold_sc < dom_def->rt2 ) {
            q_beg = q_0;
         }
         else if ( q_0 == -1 ) {
            q_beg = q_0;
         }
         
         /* this logic is explained in hmmer code: xref J2/101 */
         threshold_sc = dom_def->m_occ->data[q_0];

         if ( threshold_sc >= dom_def->rt1 ) {
            is_triggered = true;
         }
      }
      else /* if we have found the start of a potential new domain */
      { 
         /* this logic is explained in hmmer code: xref J2/101 */
         threshold_sc = dom_def->m_occ->data[q_0] - ( dom_def->e_tot->data[q_0] - dom_def->e_tot->data[q_1] );

         /* if we have found the end of a domain */
         if ( threshold_sc < dom_def->rt2 )
         {
            q_end = q_0;
            q_size = q_end - q_beg + 1;
            printf("REGION: (%d-%d)\n", q_beg, q_end);

            /* grow matrices to prepare for forward/backward over domain */
            MATRIX_3D_Reuse( st_MX_fwd, NUM_NORMAL_STATES, 2, Q+1 );
            MATRIX_2D_Reuse( sp_MX_fwd, NUM_SPECIAL_STATES, Q+1 );
            MATRIX_3D_Reuse( st_MX_bck, NUM_NORMAL_STATES, 2, Q+1 );
            MATRIX_2D_Reuse( sp_MX_bck, NUM_SPECIAL_STATES, Q+1 );
            dom_def->n_regions++;

            is_multidomain_region = test_Multidomain_Region( dom_def, q_beg, q_end );

            /* verifies whether is a multidomain region (domains overlap) */
            if ( is_multidomain_region == true )
            {
               dom_def->n_clustered++;
               
               q2_end_last = 0;
               for ( int d = 0; d < n_cluster; d++ ) {
                  /* === p7_spensemble_GetClusterCoords() === */

                  if ( q2_beg <= q2_end_last ) {
                     dom_def->n_overlaps++;
                  }
                  dom_def->n_envelopes++;

                  /* === rescore_isolated_domain() === */
               }
            }
            else /* single domain region (domains do not overlap) */
            {
               dom_def->n_envelopes++;
               
               Q_range = (RANGE){q_beg,q_end};
               run_Rescore_Isolated_Domain(
                  q_seq, t_prof, &Q_range, T, st_MX_fwd, sp_MX_fwd, st_MX_bck, sp_MX_bck, aln, dom_def );
            }
            q_beg          = -1;
            is_triggered   = false;
         }
      }
   }

   /* reconfigure model to its state from beginning of function */
   if ( t_prof->mode == MODE_MULTIGLOCAL|| t_prof->mode == MODE_MULTILOCAL ) {
      HMM_PROFILE_ReconfigMultihit( t_prof, q_seq->N );
   }
   else {
      HMM_PROFILE_ReconfigUnihit( t_prof, q_seq->N );
   }
   return STATUS_SUCCESS;
}


/** FUNCTION:  test_Multidomain_Region()
 *  SYNOPSIS:  Returns whether sequence range <q_beg,q_end> contains multiple domains.
 *             NOTES: More precisely: return TRUE if  \max_z [ \min (B(z), E(z)) ]  >= rt3
 *             where
 *             E(z) = expected number of E states occurring in region before z is emitted
 *                = \sum_{y=i}^{z} eocc[i]  =  etot[z] - etot[i-1]
 *             B(z) = expected number of B states occurring in region after z is emitted
 *                = \sum_{y=z}^{j} bocc[i]  =  btot[j] - btot[z-1]
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
bool
test_Multidomain_Region(   DOMAIN_DEF*    dom_def,
                           int            q_beg,
                           int            q_end )
{
   int   z;
   int   max;
   float expected_n;
   bool  is_multidomain_region;

   max = -1.0;
   for ( z = q_beg; z <= q_end; z++ ) 
   {
      expected_n  = MIN( 
                     (dom_def->e_tot->data[z] - dom_def->e_tot->data[q_beg - 1]),      /* E(z) */
                     (dom_def->b_tot->data[q_end] - dom_def->b_tot->data[z - 1]) );    /* B(z) */
      max = MAX( max, expected_n );
   }
   is_multidomain_region = ( ( max >= dom_def->rt3 ) ? true : false );
   return is_multidomain_region;
}

/** FUNCTION:  run_Decode_Posterior_Quad()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create normal state posterior into <...post>.
 *             Can store matrix in <...fwd> or <...bck>.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Posterior_Quad( SEQUENCE*         q_seq,            /* query sequence */
                           HMM_PROFILE*      t_prof,           /* target hmm model */
                           int               Q,                /* query length */
                           int               T,                /* target length */
                           MATRIX_3D*        st_MX_fwd,        /* normal state matrix for forward */
                           MATRIX_2D*        sp_MX_fwd,        /* special state matrix for forward */
                           MATRIX_3D*        st_MX_bck,        /* normal state matrix for backward */
                           MATRIX_2D*        sp_MX_bck,        /* special state matrix for backward */
                           MATRIX_3D*        st_MX_post,       /* OUTPUT: normal state matrix for posterior */
                           MATRIX_2D*        sp_MX_post )      /* OUTPUT: special state matrix for posterior */
{
   run_Decode_Normal_Posterior_Quad( 
      q_seq, t_prof, Q, T, st_MX_fwd, st_MX_bck, st_MX_post );
      
   run_Decode_Special_Posterior_Quad( 
      q_seq, t_prof, Q, T, sp_MX_fwd, sp_MX_bck, sp_MX_post );

}

/** FUNCTION:  run_Decode_Normal_Posterior_Quad()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store matrix in <...fwd> or <...bck>.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Normal_Posterior_Quad(   SEQUENCE*         q_seq,            /* query sequence */
                                    HMM_PROFILE*      t_prof,           /* target hmm model */
                                    int               Q,                /* query length */
                                    int               T,                /* target length */
                                    MATRIX_3D*        st_MX_fwd,        /* normal state matrix for forward */
                                    MATRIX_3D*        st_MX_bck,        /* normal state matrix for backward */
                                    MATRIX_3D*        st_MX_post )      /* OUTPUT: normal state matrix for posterior */
{
   /* query index */
   int q_0, q_1;
   /* target index */ 
   int t_0, t_1;

   /* init normal state */
   q_0 = 0;
   for ( t_0 = 0; t_0 <= T; t_0++ )
   {
      MMX_X(st_MX_post, q_0, t_0) = -INF;
      IMX_X(st_MX_post, q_0, t_0) = -INF;
      DMX_X(st_MX_post, q_0, t_0) = -INF; 
   }
   
   /* every position in query */
   for ( q_0 = 1; q_0 <= Q; q_0++ )
   {
      t_0 = 0;
      MMX_X(st_MX_post, q_0, t_0) = -INF;
      IMX_X(st_MX_post, q_0, t_0) = -INF;
      DMX_X(st_MX_post, q_0, t_0) = -INF; 

      /* every position in target */
      for ( t_0 = 1; t_0 < T; t_0++ )
      {
         MMX_X(st_MX_post, q_0, t_0) = logsum(  MMX_X(st_MX_fwd, q_0, t_0), 
                                                MMX_X(st_MX_bck, q_0, t_0) );
         IMX_X(st_MX_post, q_0, t_0) = logsum(  IMX_X(st_MX_fwd, q_0, t_0), 
                                                IMX_X(st_MX_bck, q_0, t_0) );
         DMX_X(st_MX_post, q_0, t_0) = logsum(  DMX_X(st_MX_fwd, q_0, t_0),
                                                DMX_X(st_MX_bck, q_0, t_0) ); 
      }

      /* unrolled loop */
      t_0 = T;
      MMX_X(st_MX_post, q_0, t_0) = logsum(  MMX_X(st_MX_fwd, q_0, t_0), 
                                             MMX_X(st_MX_bck, q_0, t_0) );
      IMX_X(st_MX_post, q_0, t_0) = -INF;
      DMX_X(st_MX_post, q_0, t_0) = -INF;
   }

   return STATUS_SUCCESS;
}

/** FUNCTION:  run_Decode_Special_Posterior_Quad()
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
   /* query index */
   int q_0, q_1;
   /* target index */ 
   int t_0, t_1;

   // printf("=== SPEC FWD ===\n");
   // DP_MATRIX_SpecExp_Dump( Q, T, NULL, sp_MX_fwd, stdout );

   // printf("=== SPEC BCK ===\n");
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

/** FUNCTION:  run_DecodeDomain_Posterior_Quad()
 *  SYNOPSIS:  Using posterior special state matrix <sp_MX_post> to compute domains.
 *             NOTE: Modeled after <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_DecodeDomain_Posterior_Quad( SEQUENCE*         q_seq,            /* query sequence */
                                 HMM_PROFILE*      t_prof,           /* target hmm model */
                                 int               Q,                /* query length */
                                 int               T,                /* target length */
                                 MATRIX_2D*        sp_MX_post,       /* special state matrix for posterior */ 
                                 DOMAIN_DEF*       dom_def )         /* OUTPUT: domain data */ 
{
   /* query index */
   int q_0, q_1;
   /* target index */ 
   int t_0, t_1;
   /* cumulative probability of core model states */
   float njc_pr;

   /* grow vector to necessary size */
   VECTOR_FLT_GrowTo(dom_def->b_tot, Q+1);
   VECTOR_FLT_GrowTo(dom_def->e_tot, Q+1);
   VECTOR_FLT_GrowTo(dom_def->m_occ, Q+1);
   VECTOR_FLT_GrowTo(dom_def->null2_sc, Q+1);

   /* init zero row */
   q_0 = 0;
   dom_def->b_tot->data[0] = 0.0;
   dom_def->e_tot->data[0] = 0.0;
   dom_def->m_occ->data[0] = 0.0;

   for ( q_0 = 1; q_0 <= Q; q_0++ )
   {
      q_1 = q_0 - 1;
      // printf("[%d] B(i-1)= %6.3f E(i)= %6.3f\n", 
      //    q_0, XMX_X(sp_MX_post, SP_B, q_1), XMX_X(sp_MX_post, SP_E, q_0) );

      /* cumulative probability of begin state */
      dom_def->b_tot->data[q_0] = logsum( dom_def->b_tot->data[q_1], 
                                          XMX_X(sp_MX_post, SP_B, q_1) ); 
      /* cumulative probability of end state */
      dom_def->e_tot->data[q_0] = logsum( dom_def->b_tot->data[q_0], 
                                          XMX_X(sp_MX_post, SP_E, q_0) ); 
      /* cumulative proabaility emitted by core model (not init or term states) */
      njc_pr  =   logsum( logsum( XMX_X(sp_MX_post, SP_N, q_1), 
                                  XMX_X(sp_MX_post, SP_N, q_0) ), 
                                  XSC_X(t_prof, SP_N, SP_LOOP) );
      njc_pr +=   logsum( logsum( XMX_X(sp_MX_post, SP_J, q_1), 
                                  XMX_X(sp_MX_post, SP_J, q_0) ), 
                                  XSC_X(t_prof, SP_J, SP_LOOP) );
      njc_pr +=   logsum( logsum( XMX_X(sp_MX_post, SP_C, q_1), 
                                  XMX_X(sp_MX_post, SP_C, q_0) ), 
                                  XSC_X(t_prof, SP_C, SP_LOOP) );
      dom_def->m_occ->data[q_0] += 1.0 - njc_pr;
   }
   return STATUS_SUCCESS;
}

/** FUNCTION:  run_DecodeDomain_Quad()
 *  SYNOPSIS:  Using posterior special state matrix <sp_MX_bck> and <sp_MX_fwd> to compute domains.
 *             NOTE: Modeled after <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_DecodeDomain_Quad(  SEQUENCE*         q_seq,            /* query sequence */
                        HMM_PROFILE*      t_prof,           /* target hmm model */
                        int               Q,                /* query length */
                        int               T,                /* target length */
                        MATRIX_2D*        sp_MX_fwd,       /* special state matrix for forward */ 
                        MATRIX_2D*        sp_MX_bck,       /* special state matrix for backward */ 
                        DOMAIN_DEF*       dom_def )         /* OUTPUT: domain data */ 
{
   /* query index */
   int q_0, q_1;
   /* target index */ 
   int t_0, t_1;
   /* cumulative probability of core model states */
   float tmp_pr;
   float njc_pr;
   float overall_logp;

   overall_logp = XMX_X(sp_MX_fwd, SP_C, Q) + XSC_X(t_prof, SP_C, SP_MOVE);
   // printf("# overall_logp: %7.4f\n", overall_logp);

   /* grow vector to necessary size */
   VECTOR_FLT_GrowTo(dom_def->b_tot, Q+1);
   VECTOR_FLT_GrowTo(dom_def->e_tot, Q+1);
   VECTOR_FLT_GrowTo(dom_def->m_occ, Q+1);
   VECTOR_FLT_GrowTo(dom_def->null2_sc, Q+1);

   /* init zero row */
   q_0 = 0;
   dom_def->b_tot->data[0] = 0.0;
   dom_def->e_tot->data[0] = 0.0;
   dom_def->m_occ->data[0] = 0.0;

   for ( q_0 = 1; q_0 <= Q; q_0++ )
   {
      q_1 = q_0 - 1;

      /* cumulative probability of begin state */
      dom_def->b_tot->data[q_0] = dom_def->b_tot->data[q_1] + 
                                  exp( XMX_X(sp_MX_fwd, SP_B, q_1) +
                                       XMX_X(sp_MX_bck, SP_B, q_1) -
                                       overall_logp ); 
      /* cumulative probability of end state */
      dom_def->e_tot->data[q_0] = dom_def->e_tot->data[q_1] +
                                  exp( XMX_X(sp_MX_fwd, SP_E, q_0) +
                                       XMX_X(sp_MX_bck, SP_E, q_0) -
                                       overall_logp ); 
      /* cumulative proabaility emitted by core model (not init or term states) */
      njc_pr  =   exp(  XMX_X(sp_MX_fwd, SP_N, q_1) +
                        XMX_X(sp_MX_bck, SP_N, q_0) + 
                        XSC_X(t_prof, SP_N, SP_LOOP) -
                        overall_logp );
      njc_pr +=   exp(  XMX_X(sp_MX_fwd, SP_J, q_1) +
                        XMX_X(sp_MX_bck, SP_J, q_0) + 
                        XSC_X(t_prof, SP_J, SP_LOOP) -
                        overall_logp );
      njc_pr +=   exp(  XMX_X(sp_MX_fwd, SP_C, q_1) + 
                        XMX_X(sp_MX_bck, SP_C, q_0) + 
                        XSC_X(t_prof, SP_C, SP_LOOP) -
                        overall_logp );
      dom_def->m_occ->data[q_0] = 1.0 - njc_pr;
   }

   for ( q_0 = 1; q_0 <= Q; q_0++ )
   {
      q_1 = q_0 - 1;
      // printf("[%2d] B(i-1)= %8.3f %8.3f %8.3f || E(i)= %8.3f %8.3f %8.3f || M_OCC(i)= %8.3f\n", 
      //    q_0, 
      //    dom_def->b_tot->data[q_0], XMX_X(sp_MX_fwd, SP_B, q_1), XMX_X(sp_MX_bck, SP_B, q_1),
      //    dom_def->e_tot->data[q_0], XMX_X(sp_MX_fwd, SP_E, q_0), XMX_X(sp_MX_bck, SP_E, q_0),
      //    dom_def->m_occ->data[q_0] );
      
      // printf("[%2d] B(i-1)= %8.3f %8.3f %8.3f || E(i)= %8.3f %8.3f %8.3f\n", 
      //    q_0, 
      //    exp(dom_def->b_tot->data[q_0]), exp(XMX_X(sp_MX_fwd, SP_B, q_1)), exp(XMX_X(sp_MX_bck, SP_B, q_1)),
      //    exp(dom_def->e_tot->data[q_0]), exp(XMX_X(sp_MX_fwd, SP_E, q_0)), exp(XMX_X(sp_MX_bck, SP_E, q_0)) );
   }

   return STATUS_SUCCESS;
}

/** FUNCTION:  run_Rescore_Isolated_Domain()
 *  SYNOPSIS:  
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Rescore_Isolated_Domain(  SEQUENCE*         q_seq,            /* query sequence */
                              HMM_PROFILE*      t_prof,           /* target hmm model */
                              RANGE*            Q_range,          /* query range */
                              int               T,                /* target length */
                              MATRIX_3D*        st_MX_fwd,        /* normal state matrix for forward */
                              MATRIX_2D*        sp_MX_fwd,        /* special state matrix for forward */
                              MATRIX_3D*        st_MX_bck,        /* normal state matrix for backward */
                              MATRIX_2D*        sp_MX_bck,        /* special state matrix for backward */
                              ALIGNMENT*        aln,              /* OUTPUT: domain alignment */
                              DOMAIN_DEF*       dom_def )         /* OUTPUT: domain data */
{
   int   Q; 
   float fwd_sc;
   float bck_sc;

   Q = Q_range->end - Q_range->beg + 1;

   // SEQUENCE_SetSubSequence( q_seq, Q_range->beg, Q_range->end );

   MATRIX_3D_Reuse( st_MX_fwd, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
   MATRIX_2D_Reuse( sp_MX_fwd, NUM_SPECIAL_STATES, Q+1 );
   MATRIX_3D_Reuse( st_MX_fwd, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
   MATRIX_2D_Reuse( sp_MX_fwd, NUM_SPECIAL_STATES, Q+1 );

   run_Forward_Linear( q_seq, t_prof, Q, T, st_MX_fwd, sp_MX_fwd, &fwd_sc );
   run_Backward_Linear( q_seq, t_prof, Q, T, st_MX_fwd, sp_MX_fwd, &bck_sc );

   return STATUS_SUCCESS;
}


/** FUNCTION:  run_Posterior_Viterbi_Quad()
 *  SYNOPSIS:  Run Viterbi step of Maximum Posterior Probability (no transition states).
 *             Matrices are the Maximum Poster.
 *             Computes the optimal alignment through HHM.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Posterior_Viterbi_Quad(   SEQUENCE*         query,            /* query sequence */
                              HMM_PROFILE*      target,           /* target hmm model */
                              int               Q,                /* query length */
                              int               T,                /* target length */
                              MATRIX_3D*        st_MX_post,       /* normal matrix */
                              MATRIX_2D*        sp_MX_post,       /* special matrix */
                              MATRIX_3D*        st_MX_max,        /* OUTPUT: normal matrix for max posterior */
                              MATRIX_2D*        sp_MX_max,        /* OUTPUT: special matrix for max posterior */
                              float*            sc_final )        /* OUTPUT: final max score */
{
   // /* vars for aliasing matrix accesses */
   // MATRIX_3D*  st_MX;
   // MATRIX_2D*  sp_MX;

   // /* vars for accessing query/target data structs */
   // char     a;                               /* store current character in sequence */
   // int      A;                               /* store int value of character */
   // char*    seq;                             /* alias for getting seq */
   // int      N;                               /* length of edgebound list */
   // bool     is_local;                        /* whether using local or global alignments */

   // /* vars for indexing into data matrices by row-col */
   // int      b, d, i, j, k;                   /* antidiagonal, row, column indices */
   // int      q_0, q_1;                        /* real index of current and previous rows (query) */
   // int      qx0, qx1;                        /* mod mapping of column index into data matrix (query) */
   // int      t_0, t_1;                        /* real index of current and previous columns (target) */

   // /* vars for indexing into data matrices by anti-diag */
   // int      d_0, d_1, d_2;                   /* real index of current and previous antidiagonals */
   // int      dx0, dx1, dx2;                   /* mod mapping of antidiagonal index into data matrix */
   // int      k_0, k_1;                        /* offset into antidiagonal */
   // int      d_st, d_end, d_cnt;              /* starting and ending diagonal indices */
   // int      dim_T, dim_Q, dim_TOT;           /* dimensions of submatrix being searched */
   // int      dim_min, dim_max;                /* diagonal index where num cells reaches highest point and diminishing point */ 
   // int      num_cells;                       /* number of cells in current diagonal */

   // /* vars for indexing into edgebound lists */
   // BOUND*   bnd;                             /* current bound */
   // int      id;                              /* id in edgebound list (row/diag) */
   // int      r_0;                             /* current index for current row */
   // int      r_0b, r_0e;                      /* begin and end indices for current row in edgebound list */
   // int      r_1;                             /* current index for previous row */
   // int      r_1b, r_1e;                      /* begin and end indices for current row in edgebound list */
   // int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   // int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   // int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   // int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   // bool     rb_T;                            /* checks if edge touches right bound of matrix */

   // /* vars for recurrance scores */
   // float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   // float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   // float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   // float    prv_loop, prv_move;            /* previous loop and move for special states */
   // float    prv_sum, prv_best;             /* temp sub_totaling vars */
   // float    sc_best;                         /* final best scores */
   // float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */
   
   // /* debugger tools */
   // FILE*       dbfp;
   // MATRIX_2D*  cloud_MX;
   // MATRIX_2D*  cloud_MX3;
   // MATRIX_3D*  test_MX;
   // MATRIX_3D*  test_MX3;
   // int         num_writes;
   // int         num_clears;

   // /* initialize debugging matrix */
   // #if DEBUG
   // {
   //    cloud_MX    = debugger->cloud_MX;
   //    cloud_MX3   = debugger->cloud_MX3;
   //    test_MX     = debugger->test_MX;
   //    test_MX3    = debugger->test_MX3;

   //    MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
   //    MATRIX_2D_Fill( cloud_MX, 0 );
   //    MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
   //    MATRIX_2D_Fill( cloud_MX3, 0 );
   //    MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
   //    MATRIX_3D_Fill( test_MX, -INF );
   //    MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
   //    MATRIX_3D_Fill( test_MX3, -INF );

   //    num_writes = 0;
   //    num_clears = 0;
   // }
   // #endif

   // /* --------------------------------------------------------------------------------- */

   // /* query sequence */
   // seq         = query->seq;
   // st_MX       = st_MX_post;
   // sp_MX       = sp_MX_post;
   // /* local or global alignments? */
   // is_local    = target->isLocal;
   // sc_E        = (is_local) ? 0 : -INF;

   // q_0 = 0;
   // qx0 = q_0;

   // /* FOR every position in QUERY seq */
   // for (q_0 = 1; q_0 <= Q; q_0++)
   // {
   //    q_1 = q_0 - 1;
   //    qx0 = q_0;
   //    qx1 = q_1;

   //    /* FOR every position in TARGET profile */
   //    for (t_0 = 1; t_0 < T; t_0++)
   //    {
   //       t_1 = t_0 - 1;

   //       /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
   //       /* best previous state transition (match takes the diag element of each prev state) */
   //       prv_M = MMX(qx1, t_1);
   //       prv_I = IMX(qx1, t_1);
   //       prv_D = DMX(qx1, t_1);
   //       prv_B = XMX(SP_B, q_1); /* from begin match state (new alignment) */
   //       /* best-to-match */
   //       prv_best = calc_Max( 
   //                      calc_Max( prv_M, prv_I ), 
   //                      calc_Max( prv_D, prv_B ) );
   //       MMX(qx0, t_0)  = prv_best;

   //       /* FIND BEST PATH TO INSERT STATE (FROM MATCH OR INSERT) */
   //       /* previous states (match takes the left element of each state) */
   //       prv_M = MMX(qx1, t_0);
   //       prv_I = IMX(qx1, t_0);
   //       /* best-to-insert */
   //       prv_best = calc_Max(prv_M, prv_I);
   //       IMX(qx0, t_0) = prv_best;

   //       /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
   //       /* previous states (match takes the left element of each state) */
   //       prv_M = MMX(qx0, t_1);
   //       prv_D = DMX(qx0, t_1);
   //       /* best-to-delete */
   //       prv_best = calc_Max(prv_M, prv_D);
   //       DMX(qx0, t_0) = prv_best;

   //       /* UPDATE E STATE */
   //       prv_E = XMX(SP_E, q_0);
   //       prv_M = MMX(qx0, t_0);
   //       /* best-to-e-state */
   //       XMX(SP_E, q_0) = calc_Max( prv_E, prv_M );

   //       /* embed linear row into quadratic test matrix */
   //       #if DEBUG
   //       {
   //          MX_2D(cloud_MX, q_0, t_0) = 1.0;
   //          MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
   //          MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
   //          MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
   //       }
   //       #endif
   //    }

   //    /* UNROLLED FINAL LOOP ITERATION */
   //    t_0 = T;
   //    t_1 = t_0 - 1;

   //    /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
   //    /* best previous state transition (match takes the diag element of each prev state) */
   //    prv_M = MMX(qx1, t_1);
   //    prv_I = IMX(qx1, t_1);
   //    prv_D = DMX(qx1, t_1);
   //    prv_B = XMX(SP_B, q_1); /* from begin match state (new alignment) */
   //    /* best-to-match */
   //    prv_best = calc_Max(
   //                   calc_Max( prv_M, prv_I ),
   //                   calc_Max( prv_D, prv_B ) );
   //    MMX(qx0, t_0) = prv_best;

   //    /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
   //    /* previous states (match takes the left element of each state) */
   //    prv_M = MMX(qx0, t_1);
   //    prv_D = DMX(qx0, t_1);
   //    /* best-to-delete */
   //    prv_best = calc_Max( prv_M, prv_D );
   //    DMX(qx0, t_0) = prv_best;

   //    /* UPDATE E STATE */
   //    prv_E = XMX(SP_E, q_0);
   //    prv_M = MMX(qx0, t_0);
   //    prv_D = DMX(qx0, t_0);
   //    XMX(SP_E, q_0) = calc_Max( prv_E, 
   //                         calc_Max( prv_M, prv_D ) );

   //    /* SPECIAL STATES */
   //    /* J state */
   //    prv_J = XMX(SP_J, q_1);     /* J->J */
   //    prv_E = XMX(SP_E, q_0);     /* E->J is E's "loop" */
   //    XMX(SP_J, q_0) = calc_Max( prv_J, prv_E );

   //    /* C state */
   //    prv_C = XMX(SP_C, q_1);
   //    prv_E = XMX(SP_E, q_0);
   //    XMX(SP_C, q_0) = calc_Max( prv_C, prv_E );

   //    /* N state */
   //    prv_N = XMX(SP_N, q_1);
   //    XMX(SP_N, q_0) = prv_N;

   //    /* B state */
   //    prv_N = XMX(SP_N, q_0);     /* N->B is N's move */
   //    prv_J = XMX(SP_J, q_0);     /* J->B is J's move */
   //    XMX(SP_B, q_0) = calc_Max( prv_N, prv_J );

   //    /* embed linear row into quadratic test matrix */
   //    #if DEBUG
   //    {
   //       MX_2D(cloud_MX, q_0, t_0) = 1.0;
   //       MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
   //       MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
   //       MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
   //    }
   //    #endif
   // }

   // /* T state (stores final state score) */
   // sc_best = XMX(SP_C, Q);
   // *sc_final = sc_best;

   // return STATUS_SUCCESS;
}

/** FUNCTION:  run_Null2_By_Expectation()
 *  SYNOPSIS:  Modeled after HMMER p7_GNull2_ByExpectation().
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Null2_ByExpectation(   SEQUENCE*         query,            /* query sequence */
                           HMM_PROFILE*      target,           /* target hmm model */
                           int               Q,                /* query length */
                           int               T,                /* target length */
                           MATRIX_3D*        st_MX_post,       /* posterior normal matrix */
                           MATRIX_2D*        sp_MX_post,       /* posterior special matrix */
                           DOMAIN_DEF*       dom_def )         /* OUTPUT: domain def's null2_sc vector */
{
   int   q_0, t_0, st_0, k_0;
   float x_factor;
   
   VECTOR_FLT_GrowTo( dom_def->st_freq, (Q+1) * NUM_NORMAL_STATES );
   VECTOR_FLT_GrowTo( dom_def->sp_freq, NUM_SPECIAL_STATES );
   VECTOR_FLT_GrowTo( dom_def->null2_sc, NUM_AMINO_PLUS_SPEC );
   
   /* initialize vectors */
   t_0 = 0;
   for ( q_0 = 0; q_0 < Q; q_0++ ) {
      for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
         dom_def->st_freq->data[(q_0*NUM_NORMAL_STATES) + st_0] = -INF;
      }
   }
   for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
      dom_def->sp_freq->data[st_0] = -INF;
   }

   /* add emissions from remaining positions into vectors  */
   for ( int t_0 = 1; t_0 <= T; t_0++ ) 
   {
      /* normal state emissions */
      for ( q_0 = 0; q_0 < Q; q_0++ ) {
         for ( st_0 = 0; st_0 < NUM_NORMAL_STATES; st_0++ ) {
            dom_def->st_freq->data[(q_0*NUM_NORMAL_STATES) + st_0] = logsum(  dom_def->st_freq->data[(q_0*NUM_NORMAL_STATES) + st_0],
                                                                              MX_3D( st_MX_post, st_0, q_0, t_0 ) );
         }
      }
      /* special state emissions */
      for ( int st_0 = 0; st_0 < NUM_SPECIAL_STATES; st_0++ ) {
         dom_def->sp_freq->data[st_0] = logsum( dom_def->sp_freq->data[st_0],
                                                MX_2D( sp_MX_post, st_0, t_0 ) );
      }
   }

   x_factor = XMX_X( sp_MX_post, SP_N, 0 );
   x_factor = logsum( x_factor,
                      XMX_X( sp_MX_post, SP_C, 0 ) );
   x_factor = logsum( x_factor,
                      XMX_X( sp_MX_post, SP_J, 0 ) );

   /* initialize null2 vector */
   for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      dom_def->null2_sc->data[k_0] = -INF;
   }

   /* for each amino acid */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) 
   {
      /* for each position in model */
      for ( q_0 = 1; q_0 < Q; q_0++ ) 
      {
         dom_def->null2_sc->data[k_0] = logsum( dom_def->null2_sc->data[k_0],
                                                MMX_X( st_MX_post, 0, q_0 ) + MSC_X( target, q_0, k_0 ) );
         dom_def->null2_sc->data[k_0] = logsum( dom_def->null2_sc->data[k_0],
                                                IMX_X( st_MX_post, 0, q_0 ) + ISC_X( target, q_0, k_0 ) );
      }
      q_0 = Q;
      dom_def->null2_sc->data[k_0] = logsum( dom_def->null2_sc->data[k_0],
                                             MMX_X( st_MX_post, 0, Q ) + MSC_X( target, q_0, k_0 ) );
      dom_def->null2_sc->data[k_0] = logsum( dom_def->null2_sc->data[k_0],
                                             x_factor );
   }

   /* convert log to normal scores */
   for ( k_0 = 0; k_0 < NUM_AMINO; k_0++ ) {
      printf("[%d] %8.3f\n", k_0, dom_def->null2_sc->data[k_0] );
      dom_def->null2_sc->data[k_0] = exp( dom_def->null2_sc->data[k_0] );
   }

   /* set scores for all degenerencies */
   for ( k_0 = NUM_AMINO; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      dom_def->null2_sc->data[k_0] = 1.0;
   }

   dom_def->seq_bias = 0.0;
   for ( k_0 = 0; k_0 < NUM_AMINO_PLUS_SPEC; k_0++ ) {
      printf("[%d] %8.3f\n", k_0, dom_def->null2_sc->data[k_0] );
      dom_def->seq_bias += dom_def->null2_sc->data[k_0];
   }

   return STATUS_SUCCESS;
}