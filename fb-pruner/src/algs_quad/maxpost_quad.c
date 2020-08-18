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

/* header */
#include "maxpost_quad.h"

/*  
 *  FUNCTION:  run_MaxPost_Quad()
 *  SYNOPSIS:  Compute the Max Posterior Probability and obtain optimal alignment.
 *             Computes the Forward and Backward.  Max Posterior computed by Viterbi.
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int run_MaxPost_Quad(   const SEQUENCE*      query,            /* query sequence */
                        const HMM_PROFILE*   target,           /* target HMM model */
                        const int            Q,                /* query length */
                        const int            T,                /* target length */
                        MATRIX_3D*           st_MX_fwd,        /* normal matrix for forward */
                        MATRIX_2D*           sp_MX_fwd,        /* special matrix for forward */
                        MATRIX_3D*           st_MX_bck,        /* normal matrix for backward */
                        MATRIX_2D*           sp_MX_bck,        /* special matrix for backward */
                        MATRIX_3D*           st_MX_post,       /* OUTPUT: normal matrix for posterior (can overwrite fwd and bck data) */
                        MATRIX_2D*           sp_MX_post,       /* OUTPUT: special matrix for posterior (can overwrite fwd and bck data) */
                        MATRIX_3D*           st_MX_max,        /* OUTPUT: normal matrix for max posterior (can overwrite fwd and bck data) */
                        MATRIX_2D*           sp_MX_max,        /* OUTPUT: special matrix for max posterior (can overwrite fwd and bck data) */
                        ALIGNMENT*           aln,              /* OUTPUT: alignment */
                        float*               sc_final )        /* OUTPUT: final max score */
{
   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   /* Forward */
   run_Forward_Quad( query, target, Q, T, st_MX_fwd, sp_MX_fwd, sc_final );

   /* Backward */
   run_Backward_Quad( query, target, Q, T, st_MX_bck, sp_MX_bck, sc_final );

   /* Posterior */
   run_MaxPost_Posterior_Quad( query, target, Q, T, st_MX_fwd, sp_MX_fwd, st_MX_bck, sp_MX_bck, st_MX_post, sp_MX_post, sc_final );

   // /* Viterbi for Posterior */
   // run_MaxPost_Viterbi_Quad( query, target, Q, T, st_MX_post, sp_MX_post, st_MX_max, sp_MX_max, sc_final );

   /* Traceback */
   run_MaxPost_Traceback_Quad( query, target, Q, T, st_MX_max, sp_MX_max, aln );
}

/*
 *  FUNCTION:  run_MaxPost_Posterior_Quad()
 *  SYNOPSIS:  Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in Forward Matrices.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_MaxPost_Posterior_Quad(  const SEQUENCE*      query,            /* query sequence */
                                 const HMM_PROFILE*   target,           /* target hmm model */
                                 const int            Q,                /* query length */
                                 const int            T,                /* target length */
                                 MATRIX_3D*           st_MX_fwd,        /* normal matrix for forward */
                                 MATRIX_2D*           sp_MX_fwd,        /* special matrix for forward */
                                 MATRIX_3D*           st_MX_bck,        /* normal matrix for backward */
                                 MATRIX_2D*           sp_MX_bck,        /* special matrix for backward */
                                 MATRIX_3D*           st_MX_post,       /* OUTPUT: normal matrix for posterior */
                                 MATRIX_2D*           sp_MX_post,       /* OUTPUT: special matrix for posterior */
                                 float*               sc_final )        /* OUTPUT: final max score for posterior */
{
   /* computes Posterior for Normal States */
   for (int st = 0; st <= NUM_NORMAL_STATES; st++) {
      for (int q_0 = 0; q_0 <= Q; q_0++) {
         for (int t_0 = 0; t_0 <= T; t_0++) {
            double fwd_sc = MX_3D( st_MX_fwd, st, q_0, t_0 );
            double bck_sc = MX_3D( st_MX_bck, st, q_0, t_0 );
            MX_3D( st_MX_post, st, q_0, t_0 ) = logsum( fwd_sc, bck_sc );
         }
      }
   }

   /* compute Posterior for Special States */
   for (int st = 0; st <= NUM_SPECIAL_STATES; st++) {
      for (int q_0 = 0; q_0 <= Q; q_0++) {
         double fwd_sc = MX_2D( sp_MX_fwd, st, q_0 );
         double bck_sc = MX_2D( sp_MX_bck, st, q_0 );
         MX_2D( sp_MX_post, st, q_0 );
      }
   }

   /* TODO: get sc_final? */
}


/*
 *  FUNCTION:  run_MaxPost_Viterbi_Quad()
 *  SYNOPSIS:  Run Viterbi step of Maximum Posterior Probability (no transition states).
 *             Matrices are the Maximum Poster.
 *             Computes the optimal alignment through HHM.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_MaxPost_Viterbi_Quad(    const SEQUENCE*      query,            /* query sequence */
                                 const HMM_PROFILE*   target,           /* target hmm model */
                                 const int            Q,                /* query length */
                                 const int            T,                /* target length */
                                 MATRIX_3D*           st_MX_post,       /* normal matrix */
                                 MATRIX_2D*           sp_MX_post,       /* special matrix */
                                 MATRIX_3D*           st_MX_max,        /* OUTPUT: normal matrix for max posterior */
                                 MATRIX_2D*           sp_MX_max,        /* OUTPUT: special matrix for max posterior */
                                 float*               sc_final )        /* OUTPUT: final max score */
{
   /* vars for aliasing matrix accesses */
   MATRIX_3D*  st_MX;
   MATRIX_2D*  sp_MX;

   /* vars for accessing query/target data structs */
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   char*    seq;                             /* alias for getting seq */
   int      N;                               /* length of edgebound list */
   bool     is_local;                        /* whether using local or global alignments */

   /* vars for indexing into data matrices by row-col */
   int      b, d, i, j, k;                   /* antidiagonal, row, column indices */
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* mod mapping of column index into data matrix (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */

   /* vars for indexing into data matrices by anti-diag */
   int      d_0, d_1, d_2;                   /* real index of current and previous antidiagonals */
   int      dx0, dx1, dx2;                   /* mod mapping of antidiagonal index into data matrix */
   int      k_0, k_1;                        /* offset into antidiagonal */
   int      d_st, d_end, d_cnt;              /* starting and ending diagonal indices */
   int      dim_T, dim_Q, dim_TOT;           /* dimensions of submatrix being searched */
   int      dim_min, dim_max;                /* diagonal index where num cells reaches highest point and diminishing point */ 
   int      num_cells;                       /* number of cells in current diagonal */

   /* vars for indexing into edgebound lists */
   BOUND*   bnd;                             /* current bound */
   int      id;                              /* id in edgebound list (row/diag) */
   int      r_0;                             /* current index for current row */
   int      r_0b, r_0e;                      /* begin and end indices for current row in edgebound list */
   int      r_1;                             /* current index for previous row */
   int      r_1b, r_1e;                      /* begin and end indices for current row in edgebound list */
   int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   bool     rb_T;                            /* checks if edge touches right bound of matrix */

   /* vars for recurrance scores */
   float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;            /* previous loop and move for special states */
   float    prv_sum, prv_best;             /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */
   
   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_2D*  cloud_MX3;
   MATRIX_3D*  test_MX;
   MATRIX_3D*  test_MX3;
   int         num_writes;
   int         num_clears;

   /* initialize debugging matrix */
   #if DEBUG
   {
      cloud_MX    = debugger->cloud_MX;
      cloud_MX3   = debugger->cloud_MX3;
      test_MX     = debugger->test_MX;
      test_MX3    = debugger->test_MX3;

      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
      MATRIX_2D_Fill( cloud_MX3, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
      MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
      MATRIX_3D_Fill( test_MX3, -INF );

      num_writes = 0;
      num_clears = 0;
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* query sequence */
   seq         = query->seq;
   st_MX       = st_MX_post;
   sp_MX       = sp_MX_post;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   q_0 = 0;
   qx0 = q_0;

   /* FOR every position in QUERY seq */
   for (q_0 = 1; q_0 <= Q; q_0++)
   {
      q_1 = q_0 - 1;
      qx0 = q_0;
      qx1 = q_1;

      /* FOR every position in TARGET profile */
      for (t_0 = 1; t_0 < T; t_0++)
      {
         t_1 = t_0 - 1;

         /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         prv_M = MMX(qx1, t_1);
         prv_I = IMX(qx1, t_1);
         prv_D = DMX(qx1, t_1);
         prv_B = XMX(SP_B, q_1); /* from begin match state (new alignment) */
         /* best-to-match */
         prv_best = calc_Max( 
                        calc_Max( prv_M, prv_I ), 
                        calc_Max( prv_D, prv_B ) );
         MMX(qx0, t_0)  = prv_best;

         /* FIND BEST PATH TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         prv_M = MMX(qx1, t_0);
         prv_I = IMX(qx1, t_0);
         /* best-to-insert */
         prv_best = calc_Max(prv_M, prv_I);
         IMX(qx0, t_0) = prv_best;

         /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         prv_M = MMX(qx0, t_1);
         prv_D = DMX(qx0, t_1);
         /* best-to-delete */
         prv_best = calc_Max(prv_M, prv_D);
         DMX(qx0, t_0) = prv_best;

         /* UPDATE E STATE */
         prv_E = XMX(SP_E, q_0);
         prv_M = MMX(qx0, t_0);
         /* best-to-e-state */
         XMX(SP_E, q_0) = calc_Max( prv_E, prv_M );

         /* embed linear row into quadratic test matrix */
         #if DEBUG
         {
            MX_2D(cloud_MX, q_0, t_0) = 1.0;
            MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
            MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
            MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
         }
         #endif
      }

      /* UNROLLED FINAL LOOP ITERATION */
      t_0 = T;
      t_1 = t_0 - 1;

      /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
      /* best previous state transition (match takes the diag element of each prev state) */
      prv_M = MMX(qx1, t_1);
      prv_I = IMX(qx1, t_1);
      prv_D = DMX(qx1, t_1);
      prv_B = XMX(SP_B, q_1); /* from begin match state (new alignment) */
      /* best-to-match */
      prv_best = calc_Max(
                     calc_Max( prv_M, prv_I ),
                     calc_Max( prv_D, prv_B ) );
      MMX(qx0, t_0) = prv_best;

      /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
      /* previous states (match takes the left element of each state) */
      prv_M = MMX(qx0, t_1);
      prv_D = DMX(qx0, t_1);
      /* best-to-delete */
      prv_best = calc_Max( prv_M, prv_D );
      DMX(qx0, t_0) = prv_best;

      /* UPDATE E STATE */
      prv_E = XMX(SP_E, q_0);
      prv_M = MMX(qx0, t_0);
      prv_D = DMX(qx0, t_0);
      XMX(SP_E, q_0) = calc_Max( prv_E, 
                           calc_Max( prv_M, prv_D ) );

      /* SPECIAL STATES */
      /* J state */
      prv_J = XMX(SP_J, q_1);     /* J->J */
      prv_E = XMX(SP_E, q_0);     /* E->J is E's "loop" */
      XMX(SP_J, q_0) = calc_Max( prv_J, prv_E );

      /* C state */
      prv_C = XMX(SP_C, q_1);
      prv_E = XMX(SP_E, q_0);
      XMX(SP_C, q_0) = calc_Max( prv_C, prv_E );

      /* N state */
      prv_N = XMX(SP_N, q_1);
      XMX(SP_N, q_0) = prv_N;

      /* B state */
      prv_N = XMX(SP_N, q_0);     /* N->B is N's move */
      prv_J = XMX(SP_J, q_0);     /* J->B is J's move */
      XMX(SP_B, q_0) = calc_Max( prv_N, prv_J );

      /* embed linear row into quadratic test matrix */
      #if DEBUG
      {
         MX_2D(cloud_MX, q_0, t_0) = 1.0;
         MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
         MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
         MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
      }
      #endif
   }

   /* T state (stores final state score) */
   sc_best = XMX(SP_C, Q);
   *sc_final = sc_best;

   return STATUS_SUCCESS;
}


/*
 *  FUNCTION:  run_MaxPost_Traceback_Quad()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Greedy algorithm that chooses the maximum next state.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_MaxPost_Traceback_Quad(     const SEQUENCE*     query,       /* query sequence */
                                    const HMM_PROFILE*  target,      /* HMM model */
                                    const int           Q,           /* query/seq length */
                                    const int           T,           /* target/model length */
                                    MATRIX_3D*          st_MX,       /* Normal State (Match, Insert, Delete) Matrix */
                                    MATRIX_2D*          sp_MX,       /* Special State (J,N,B,C,E) Matrix */
                                    ALIGNMENT*          aln )        /* OUTPUT: Traceback Alignment */
{
   /* vars for accessing query/target data structs */
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   char*    seq;                             /* alias for getting seq */
   bool     is_local;                        /* whether */

   /* vars for indexing into data matrices */
   int      i, j;                            /* antidiagonal, row, column indices */
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* map row index into data matrix (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */
   int      tx0, tx1;                        /* map column index into data matrix (target) */

   /* vars for recurrance scores */
   float    cur;                             /* current state score */
   float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* vars for alignment traceback */
   int            st_cur, st_prv;            /* current, previous state in traceback */
   int            st_max;                    /* previous state with maximum score */
   int            t_max;                     /* t_0 with maximum score (for B->M) */
   float          sc_cur, sc_prv;            /* current, previous score of state */
   float          sc_max;                    /* maximum score for previous state */
   TRACE*         tr;                        /* trace object for appending */
   const float    tol = 1e-5;                /* acceptable tolerance range for "equality tests" */

   int passes = 0;

   /* --------------------------------------------------------------------------- */

   /* real sequence indexes */
   q_0 = Q;
   q_1 = q_0 - 1;
   t_0 = 0;
   t_1 = t_0 - 1;
   /* map sequence indexes to data */
   qx0 = q_0;
   qx1 = q_1;
   tx0 = t_0;
   tx1 = t_1;

   seq = query->seq; 
   tr  = aln->traces->data; 

   /* local or global? */
   is_local = target->isLocal;

   /* allocate memory for trace */
   ALIGNMENT_Reuse( aln, Q, T );

   /* Backalnacing, so C is end state */
   ALIGNMENT_Append( aln, tr, T_ST, q_0, t_0 );
   ALIGNMENT_Append( aln, tr, C_ST, q_0, t_0 );
   st_prv = C_ST;

   /* End of trace is S state */
   while (st_prv != S_ST)
   {
      /* if we have reached the end of the query sequence */
      if (q_0 == 0) {
         ALIGNMENT_Append( aln, tr, S_ST, q_0, t_0 );
         break;
      } 

      /* reset maximums of previous scores */
      st_max = X_ST; /* unknown state */
      sc_max = -INF;

      /* previous query and target index */
      q_1 = q_0 - 1;
      t_1 = t_0 - 1;
      /* map sequence indexes to data */
      qx0 = q_0;
      qx1 = q_1;
      tx0 = t_0;
      tx1 = t_1;
      
      /* get next sequence character */
      a = seq[q_1];
      A = AA_REV[a];

      /* jump from current state to the prev state */
      switch ( st_prv )
      {
         /* C STATE to {C,E} */
         case C_ST:  /* C(q_0) comes from C(q_1) or E(q_0) */
         {
            /* current state */
            cur = XMX(SP_C, q_0);
            
            /* possible previous states before transition scores */
            prv_C = XMX(SP_C, q_1) + XSC(SP_C, SP_LOOP);
            prv_E = XMX(SP_E, q_0) + XSC(SP_E, SP_MOVE);

            /* find maximum of previous states */
            if ( sc_max < prv_C ) {
               st_max = C_ST;
               sc_max = prv_C;
            }
            if ( sc_max < prv_E ) {
               st_max = E_ST;
               sc_max = prv_E;
            }

            /* update next state to max */
            st_cur = st_max;
            sc_cur = sc_max;
         } break;

         /* E STATE to {M,D} */
         case E_ST:  /* E connects from any M state. t_0 is set here. */
         {
            /* current state */
            cur = XMX(SP_E, q_0);

            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible E_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            if ( is_local ) /* local mode: ends in M */
            {
               /* can't come from D, in a *local* Viterbi alignment. */
               st_cur = M_ST;

               /* next state must be M */
               st_max = M_ST;

               /* find maximum of previous states (all possible M(q,t) for current q=q_0 ) */
               for ( t_0 = T; t_0 >= 1; t_0-- ) 
               {
                  /* possible previous state */
                  prv_M = MMX(q_0, t_0);

                  /* find maximum of possible previous states */
                  if ( sc_max < prv_M ) {
                     sc_max = prv_M;
                     t_max = t_0;
                  }
               }
            }
            else /* glocal mode: we either come from D_M or M_M */
            {
               /* possible previous states */
               prv_M = MMX(q_0, T);
               prv_D = DMX(q_0, T);

               /* find maximum of possible previous states */
               if ( sc_max < prv_M ) {
                  sc_max = prv_M;
                  t_max = T;
               }
               if ( sc_max < prv_D ) {
                  sc_max = prv_D;
                  t_max = T;
               }
            }
            
            /* update next state to max */
            st_cur = st_max;
            sc_cur = sc_max;
            /* set t_0 of end of alignment */
            t_0 = t_max;
         } break;

         /* M STATE to {B,M,I,D} */
         case M_ST:  /* M connects from {M,I,D}(q_1, t_1), or B(q_1) */
         {
            /* current state */
            cur = MMX(q_0, t_0);    

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible M_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M) + MSC(t_0, A);
            prv_M = MMX(q_1, t_1) + TSC(t_1, M2M) + MSC(t_0, A);
            prv_I = IMX(q_1, t_1) + TSC(t_1, I2M) + MSC(t_0, A);
            prv_D = DMX(q_1, t_1) + TSC(t_1, D2M) + MSC(t_0, A);

            /* find maximum of previous states */
            if ( sc_max < prv_B ) {
               st_max = B_ST;
               sc_max = prv_B;
            }
            if ( sc_max < prv_M ) {
               st_max = M_ST;
               sc_max = prv_M;
            }
            if ( sc_max < prv_I ) {
               st_max = I_ST;
               sc_max = prv_I;
            }
            if ( sc_max < prv_D ) {
               st_max = D_ST;
               sc_max = prv_D;
            }

            /* update next state to max */
            st_cur = st_max;
            sc_cur = sc_max;
            /* update index of previous state */
            t_0--; q_0--;
         } break;

         /* D STATE to {M,D} */
         case D_ST:  /* D connects from M,D at (q_0, t_1) */
         {
            /* current state */
            cur = DMX(q_0, t_0);

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible D_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_M = MMX(q_0, t_1) + TSC(t_1, M2D);
            prv_D = DMX(q_0, t_1) + TSC(t_1, D2D);

            /* find maximum of previous states */
            if ( sc_max < prv_M ) {
               st_max = M_ST;
               sc_max = prv_M;
            }
            if ( sc_max < prv_D ) {
               st_max = D_ST;
               sc_max = prv_D;
            }

            /* update next state to max */
            st_cur = st_max;
            sc_cur = sc_max;
            /* update index to previous state */
            t_0--;
         } break;

         /* I STATE to {M,I} */
         case I_ST:  /* I connects from M,I at (q_1, t_0) */
         {
            /* current state */
            cur = IMX(q_0, t_0);

            /* No valid alignment goes to -INF */
            if (IMX(q_0, t_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible I_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_M = MMX(q_1, t_0) + TSC(t_0, M2I) + ISC(t_0, A);
            prv_I = IMX(q_1, t_0) + TSC(t_0, I2I) + ISC(t_0, A);

            /* find maximum of previous states */
            if ( sc_max < prv_M ) {
               st_max = M_ST;
               sc_max = prv_M;
            }
            if ( sc_max < prv_I ) {
               st_max = I_ST;
               sc_max = prv_I;
            }

            /* update next state to max */
            st_cur = st_max;
            sc_cur = sc_max;
            /* update index to previous state */
            q_0--;
         } break;

         /* N STATE to {N,S} */
         case N_ST:  /* N connects from S,N */
         {
            /* current state */
            cur = XMX(SP_N, q_0);

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible N_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* if at beginning of query sequence, then alignment completes at S state, else N state */
            if ( q_0 <= 0 ) {
               st_cur = S_ST;
            } 
            else {
               st_cur = N_ST;
            }
         } break;

         /* B STATE to {N,J} */
         case B_ST:  /* B connects from N, J */
         {
            /* current state */
            cur = XMX(SP_B, q_0);

            /* possible previous states */
            prv_N = XMX(SP_N, q_0) + XSC(SP_N, SP_MOVE);
            prv_J = XMX(SP_J, q_0) + XSC(SP_J, SP_MOVE);

            /* find maximum of previous states */
            if ( sc_max < prv_N ) {
               st_max = N_ST;
               sc_max = prv_N;
            }
            if ( sc_max < prv_J ) {
               st_max = J_ST;
               sc_max = prv_J;
            }
            /* update next state to max */
            st_cur = st_max;
            sc_cur = sc_max;
         } break;

         /* J STATE to {J,E} */
         case J_ST:  /* J connects from E(q_0) or J(q_1) */
         {
            /* current state */
            cur = XMX(SP_J, q_0);

            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible J_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_J = XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP);
            prv_E = XMX(SP_E, q_0) + XSC(SP_E, SP_LOOP);

            /* find maximum of previous states */
            if ( sc_max < prv_J ) {
               st_max = J_ST;
               sc_max = prv_J;
            }
            if ( sc_max < prv_E ) {
               st_max = E_ST;
               sc_max = prv_E;
            }
            /* update next state to max */
            st_cur = st_max;
            sc_cur = sc_max;
         } break;

         default:
         {
            fprintf( stderr, "ERROR: Hit Bogus State!!!\n");
            exit(EXIT_FAILURE);
         }
      }

      ALIGNMENT_Append( aln, tr, st_cur, q_0, t_0 );

      /* For NCJ, we deferred i decrement. (only if looping within state) */
      if ( (st_cur == N_ST || st_cur == J_ST || st_cur == C_ST) && (st_cur == st_prv) ) {
         q_0--;
      }

      /* Update previous state */
      st_prv = st_cur;
   }

   /* reverse order of alnaceback */
   ALIGNMENT_Reverse( aln );

   /* find end and begin alignment points (first and last match state) */
   int N  = aln->traces->N;
   for (int i = 0; i < N; ++i) {
      if ( tr[i].st == B_ST ) {
         VECTOR_INT_Pushback( aln->tr_beg, i + 1 );
      }
      if ( tr[i].st == E_ST ) {
         VECTOR_INT_Pushback( aln->tr_end, i - 1 );
      }
   }
   aln->beg = aln->tr_beg->data[0];
   aln->end = aln->tr_end->data[0];

   #if DEBUG
   {
      MATRIX_2D* cloud_MX = debugger->cloud_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      for ( int i = 0; i < N; i++ ) {
         if ( tr[i].st == M_ST || tr[i].st == I_ST || tr[i].st == D_ST )
            MX_2D( cloud_MX, tr[i].i, tr[i].j ) = -1.0;
      }
   }
   #endif

   return STATUS_SUCCESS;
}