/*******************************************************************************
 *  FILE:    traceback_quad.c
 *  PURPOSE: Traceback for Viterbi Algorithm (QUADRATIC SPACE).
 *
 *  AUTHOR:  Dave Rich
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
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"

/* header */
#include "_algs_quad.h"
#include "viterbi_traceback_quad.h"


/** FUNCTION:  run_Traceback_Quad()
 *  SYNOPSIS:  Selects the default method of run_Traceback_Quad() from the available methods.
 *             Requires that <st_MX> and <sp_MX> are still completely filled.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad(    const SEQUENCE*     query,       /* query sequence */
                           const HMM_PROFILE*  target,      /* HMM model */
                           const int           Q,           /* query/seq length */
                           const int           T,           /* target/model length */
                           MATRIX_3D*          st_MX,       /* Normal State (Match, Insert, Delete) Matrix */
                           MATRIX_2D*          sp_MX,       /* Special State (J,N,B,C,E) Matrix */
                           ALIGNMENT*          aln )        /* OUTPUT: Traceback Alignment */
{
   run_Traceback_Quad_via_cmp( 
      query, target, Q, T, st_MX, sp_MX, aln );
}

/** FUNCTION:  run_Traceback_Quad()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment. 
 *             Version 1: based on HMMER version. Verifies that Alignment agrees with Matrix data.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad_via_hmmer(   const SEQUENCE*     query,       /* query sequence */
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
   int      qx0, qx1;                        /* mapping of row index into data matrix (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */
   int      tx0, tx1;                        /* mapping of col index into data matrix (query) */

   /* temporary stores */
   float    cur_sc;                          /* current score */

   /* vars for alignment traceback */
   int            st_cur, st_prv;            /* current, previous state in traceback */
   TRACE*         tr;                        /* trace object for appending */
   const float    tol = 1e-5;                /* acceptable tolerance range for "equality tests" */

   /* --------------------------------------------------------------------------- */

   q_0 = Q;
   q_1 = q_0 - 1;

   t_0 = 0;
   t_1 = t_0 - 1;

   seq   = query->seq; 
   tr    = aln->traces->data; 

   /* local or global? */
   is_local = target->isLocal;

   /* allocate memory for trace */
   ALIGNMENT_Reuse( aln, Q, T );

   /* Backalnacing, so C is end state */
   ALIGNMENT_Append( aln, T_ST, q_0, t_0 );
   ALIGNMENT_Append( aln, C_ST, q_0, t_0 );
   st_prv = C_ST;

   /* End of alnace is S state */
   while (st_prv != S_ST)
   {
      if (q_0 == 0) 
      {
         ALIGNMENT_Append( aln, S_ST, q_0, t_0 );
         break;
      } 
      q_1 = q_0 - 1;
      t_1 = t_0 - 1;
      
      /* get current sequence character */
      a = seq[q_1];
      A = AA_REV[a];

      /* jump from current state to the prev state */
      switch (st_prv)
      {
         /* C STATE to {C,E} */
         case C_ST:  /* C(i) comes from C(i-1) or E(i) */
            if (XMX(SP_C, q_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible C_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            if ( CMP_TOL( XMX(SP_C, q_0), XMX(SP_C, q_1) + XSC(SP_C, SP_LOOP) ) )
               st_cur = C_ST;
            else if ( CMP_TOL( XMX(SP_C, q_0), XMX(SP_E, q_0) + XSC(SP_E, SP_MOVE) ) )
               st_cur = E_ST;
            else {
               fprintf( stderr, "ERROR: Failed to trace from B_ST at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
            break;

         /* E STATE to {M,D} */
         case E_ST:  /* E connects from any M state. k set here */
            if (XMX(SP_E, q_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible E_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            if ( is_local )  /* local mode: ends in M */
            {
               st_cur = M_ST;    /* can't come from D, in a *local* Viterbi alnace. */
               for (t_0 = T; t_0 >= 1; t_0--) {
                  // fprintf( stderr, "testing E at (%d, %d) => (%.2f v. %.2f) \n", q_0, t_0, XMX(SP_E, q_0), MMX(q_0, t_0) );
                  if ( CMP_TOL( XMX(SP_E, q_0), MMX(q_0, t_0) ) )
                     break;
               }
               if (t_0 == 0) {
                  fprintf( stderr, "ERROR: Failed to trace from E_ST at (%d,%d)\n", q_0, t_0);
                  ERRORCHECK_exit(EXIT_FAILURE);
               }
            }
            else     /* glocal mode: we either come from D_M or M_M */
            {
               if ( CMP_TOL( XMX(SP_E, q_0), MMX(q_0, T) ) ) {
                  st_cur = M_ST;
                  t_0 = T;
               }
               else if ( CMP_TOL( XMX(SP_E, q_0), DMX(q_0, T) ) ) {
                  st_cur = D_ST;
                  t_0 = T;
               }
               else {
                  fprintf( stderr, "ERROR: Failed to trace from E_ST at (%d,%d)\n", q_0, t_0);
                  ERRORCHECK_exit(EXIT_FAILURE);
               }
            }
            break;

         /* M STATE to {B,M,I,D} */
         case M_ST:  /* M connects from i-1,k-1, or B */
            if (MMX(q_0, t_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible M_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            if ( CMP_TOL( MMX(q_0, t_0), XMX(SP_B, q_1) + TSC(t_1, B2M) + MSC(t_0, A) ) )
               st_cur = B_ST;
            else if ( CMP_TOL( MMX(q_0, t_0), MMX(q_1, t_1) + TSC(t_1, M2M) + MSC(t_0, A) ) )
               st_cur = M_ST;
            else if ( CMP_TOL( MMX(q_0, t_0), IMX(q_1, t_1) + TSC(t_1, I2M) + MSC(t_0, A) ) )
               st_cur = I_ST;
            else if ( CMP_TOL( MMX(q_0, t_0), DMX(q_1, t_1) + TSC(t_1, TM) + MSC(t_0, A) ) )
               st_cur = D_ST;
            else {
               fprintf( stderr, "ERROR: Failed to trace from M_ST at (%d,%d)\n", t_0, q_0);
               fprintf( stderr, "TOL: %f vs %f\n", MMX(q_0, t_0), MMX(q_1, t_1) + TSC(t_1, TM) + MSC(t_0, A) );
               ERRORCHECK_exit(EXIT_FAILURE);
            }
            t_0--; q_0--;
            break;

         /* D STATE to {M,D} */
         case D_ST:  /* D connects from M,D at i,k-1 */
            if (DMX(q_0, t_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible D_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            t_1 = t_0-1;
            if ( CMP_TOL( DMX(q_0, t_0), MMX(q_0, t_1) + TSC(t_1, M2D) ) )
               st_cur = M_ST;
            else if ( CMP_TOL( DMX(q_0, t_0), DMX(q_0, t_1) + TSC(t_1, TD) ) )
               st_cur = D_ST;
            else {
               fprintf( stderr, "ERROR: Failed to alnace from D_ST at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
            t_0--;
            break;

         /* I STATE to {M,I} */
         case I_ST:  /* I connects from M,I at i-1,k */
            if (IMX(q_0, t_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible I_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            if ( CMP_TOL( IMX(q_0, t_0), MMX(q_1, t_0) + TSC(t_0, M2I) + ISC(t_0, A) ) )
               st_cur = M_ST;
            else if ( CMP_TOL( IMX(q_0, t_0), IMX(q_1, t_0) + TSC(t_0, I2I) + ISC(t_0, A) ) )
               st_cur = I_ST;
            else {
               fprintf( stderr, "ERROR: Failed to alnace from I_ST at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
            q_0--;
            break;

         /* N STATE to {N,S} */
         case N_ST:  /* N connects from S, N */
            if (XMX(SP_N, q_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible N_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            st_cur = ( (q_0 <= 0) ? S_ST : N_ST );
            break;

         /* B STATE to {N,J} */
         case B_ST:  /* B connects from N, J */
            if ( CMP_TOL( XMX(SP_B, q_0), XMX(SP_N, q_0) + XSC(SP_N, SP_MOVE) ) )
               st_cur = N_ST;
            else if ( CMP_TOL( XMX(SP_B, q_0), XMX(SP_J, q_0) + XSC(SP_J, SP_MOVE) ) )
               st_cur = J_ST;
            else {
               fprintf( stderr, "ERROR: Failed to alnace from B_ST at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
            break;

         /* J STATE to {J,E} */
         case J_ST:  /* J connects from E(i) or J(i-1) */
            if (XMX(SP_J, q_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible J_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            if ( CMP_TOL( XMX(SP_J, q_0), XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP) ) )
               st_cur = J_ST;
            else if ( CMP_TOL( XMX(SP_J, q_0), XMX(SP_E, q_0) + XSC(SP_E, SP_LOOP) ) )
               st_cur = E_ST;
            else {
               fprintf( stderr, "ERROR: Failed to alnace from J_ST at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
            break;

         default:
            fprintf( stderr, "ERROR: Hit Bogus State!!!\n");
            ERRORCHECK_exit(EXIT_FAILURE);
      }

      ALIGNMENT_Append( aln, st_cur, q_0, t_0 );

      /* For NCJ, we deferred i decrement. */
      if ( (st_cur == N_ST || st_cur == J_ST || st_cur == C_ST) && st_cur == st_prv) {
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
            MX_2D( cloud_MX, tr[i].q_0, tr[i].t_0 ) = -1.0;
      }
   }
   #endif

   return STATUS_SUCCESS;
}

/** FUNCTION:  run_Traceback_Quad_via_cmp()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Version 2: My implementation. Takes maximum next step by finding the state that fulfills equation ( <previous state> + <transition> + <score> == <current state> ).
 *             Verifies that Alignment agrees with Matrix data.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad_via_cmp(     const SEQUENCE*     query,       /* query sequence */
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
   int      qx0, qx1;                        /* maps row index into data matrix (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */
   int      tx0, tx1;                        /* maps column index into data matrix (target) */

   /* vars for recurrance scores */
   float    cur;                             /* current state */
   float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;              /* previous loop and move for special states */
   float    prv_sum, prv_best;               /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* vars for alignment traceback */
   int            st_cur, st_prv;            /* current, previous state in traceback */
   TRACE*         tr;                        /* trace object for appending */
   const float    tol = 1e-5;                /* acceptable tolerance range for "equality tests" */

   /* --------------------------------------------------------------------------- */

   /* real sequence indexes */
   q_0 = Q;
   q_1 = q_0 - 1;
   t_0 = 0;
   t_1 = t_0 - 1;
   /* map sequence index to data index */
   qx0 = q_0;
   qx1 = qx1;
   tx0 = t_0;
   tx1 = tx1;

   seq   = query->seq; 
   tr    = aln->traces->data; 

   /* local or global? */
   is_local = target->isLocal;

   /* allocate memory for trace */
   ALIGNMENT_Reuse( aln, Q, T );

   /* Backalnacing, so C is end state */
   ALIGNMENT_Append( aln, T_ST, q_0, t_0 );
   ALIGNMENT_Append( aln, C_ST, q_0, t_0 );
   st_prv = C_ST;

   /* End of alnace is S state */
   while (st_prv != S_ST)
   {
      /* if we have reached the end of the query sequence */
      if (q_0 == 0) {
         ALIGNMENT_Append( aln, S_ST, q_0, t_0 );
         break;
      } 

      /* real sequence indexes */
      q_1 = q_0 - 1;
      t_1 = t_0 - 1;
      /* map sequence index to data index */
      qx0 = q_0;
      qx1 = qx1;
      tx0 = t_0;
      tx1 = tx1;
      
      /* get next sequence character */
      a = seq[q_1];
      A = AA_REV[a];

      /* jump from current state to the previous state */
      switch ( st_prv )
      {
         /* C STATE to {C,E} */
         case C_ST:  /* C(q_0) comes from C(q_1) or E(q_0) */
         {
            /* current score */
            cur   = XMX(SP_C, q_0);

            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible C_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
            
            /* possible previous states */
            prv_C = XMX(SP_C, q_1) + XSC(SP_C, SP_LOOP);
            prv_E = XMX(SP_E, q_0) + XSC(SP_E, SP_MOVE);

            if ( CMP_TOL( cur, prv_C ) ) {
               st_cur = C_ST;
            }
            else if ( CMP_TOL( cur, prv_E ) ) {
               st_cur = E_ST;
            }
            else {
               fprintf( stderr, "ERROR: Failed to trace from B_ST at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         } break;

         /* E STATE to {M,D} */
         case E_ST:  /* E connects from any M state. t_0 is set here. */
         {
            /* current score */
            cur = XMX(SP_E, q_0);

            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible E_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            if ( is_local )  /* local mode: ends in M */
            {
               /* can't come from D, in a *local* Viterbi alignment. */
               st_cur = M_ST;

               /* possible previous states (any M state) */
               for ( t_0 = T; t_0 >= 1; t_0-- ) 
               {
                  /* possible previous state */
                  prv_M = MMX(q_0, t_0);

                  /* verifies if scores agree with true previous state in alignment */
                  if ( CMP_TOL( cur, prv_M ) ) {
                     break;
                  }
               }
               /* if no entry point into M found */
               if ( t_0 == 0 ) {
                  fprintf( stderr, "ERROR: Failed to trace from E_ST at (%d,%d)\n", q_0, t_0);
                  ERRORCHECK_exit(EXIT_FAILURE);
               }
            }
            else     /* glocal mode: we either come from D_M or M_M */
            {
               /* possible previous states */
               prv_M = MMX(q_0, T);
               prv_D = DMX(q_0, T);

               /* verifies if scores agree with true previous state in alignment */
               if ( CMP_TOL( cur, prv_M ) ) {
                  st_cur = M_ST;
                  t_0 = T;
               }
               else if ( CMP_TOL( cur, prv_D ) ) {
                  st_cur = D_ST;
                  t_0 = T;
               }
               else {
                  fprintf( stderr, "ERROR: Failed to trace from E_ST at (%d,%d)\n", q_0, t_0);
                  ERRORCHECK_exit(EXIT_FAILURE);
               }
            }
         } break;

         /* M STATE to {B,M,I,D} */
         case M_ST:  /* M connects from (q_1, t_1), or B */
         {
            /* current score */
            cur = MMX(q_0, t_0);      

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible M_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M) + MSC(t_0, A);
            prv_M = MMX(q_1, t_1) + TSC(t_1, M2M) + MSC(t_0, A);
            prv_I = IMX(q_1, t_1) + TSC(t_1, I2M) + MSC(t_0, A);
            prv_D = DMX(q_1, t_1) + TSC(t_1, TM) + MSC(t_0, A);

            /* verifies if scores agree with true previous state in alignment */
            if ( CMP_TOL( cur, prv_B ) ) {
               st_cur = B_ST;
            }
            else if ( CMP_TOL( cur, prv_M ) ) {
               st_cur = M_ST;
            }
            else if ( CMP_TOL( cur, prv_I ) ) {
               st_cur = I_ST;
            }
            else if ( CMP_TOL( cur, prv_D ) ) {
               st_cur = D_ST;
            }
            else {
               fprintf( stderr, "ERROR: Failed to trace from M_ST at (%d,%d)\n", t_0, q_0);
               fprintf( stderr, "TOL: %f vs %f\n", MMX(q_0, t_0), MMX(q_1, t_1) + TSC(t_1, TM) + MSC(t_0, A) );
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            /* update index to previous state */
            t_0--; q_0--;
         } break;

         /* D STATE to {M,D} */
         case D_ST:  /* D connects from M,D at (q_0, t_1) */
         {
            /* current score */
            cur = DMX(q_0, t_0);

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible D_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_M = MMX(q_0, t_1) + TSC(t_1, M2D);
            prv_D = DMX(q_0, t_1) + TSC(t_1, TD);

            /* verifies if scores agree with true previous state in alignment */
            if ( CMP_TOL( cur, prv_M ) ) {
               st_cur = M_ST;
            }
            else if ( CMP_TOL( cur, prv_D ) ) {
               st_cur = D_ST;
            }
            else {
               fprintf( stderr, "ERROR: Failed to trace from D_ST at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            /* update index to previous state */
            t_0--;
         } break;

         /* I STATE to {M,I} */
         case I_ST:  /* I connects from M,I at (q_1, t_0) */
         {
            /* current score */
            cur = IMX(q_0, t_0);

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible I_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_M = MMX(q_1, t_0) + TSC(t_0, M2I) + ISC(t_0, A);
            prv_I = IMX(q_1, t_0) + TSC(t_0, I2I) + ISC(t_0, A);

            /* verifies if scores agree with true previous state in alignment */
            if ( CMP_TOL( cur, prv_M ) ) {
               st_cur = M_ST;
            }
            else if ( CMP_TOL( cur, prv_I ) ) {
               st_cur = I_ST;
            }
            else {
               fprintf( stderr, "ERROR: Failed to trace from I_ST at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
            q_0--;
         } break;

         /* N STATE to {N,S} */
         case N_ST:  /* N connects from S,N */
         {
            /* current score */
            cur = XMX(SP_N, q_0);

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible N_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
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
            /* current score */
            cur = XMX(SP_B, q_0);

            /* possible previous states */
            prv_N = XMX(SP_N, q_0) + XSC(SP_N, SP_MOVE);
            prv_J = XMX(SP_J, q_0) + XSC(SP_J, SP_MOVE);

            /* verifies that scores agree with true previous state in alignment */
            if ( CMP_TOL( cur, prv_N ) ) {
               st_cur = N_ST;
            }
            else if ( CMP_TOL( cur, prv_J ) ) {
               st_cur = J_ST;
            }
            else {
               fprintf( stderr, "ERROR: Failed to alnace from B_ST at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         } break;

         /* J STATE to {J,E} */
         case J_ST:  /* J connects from E(q_0) or J(q_1) */
         {
            /* current score */
            cur = XMX(SP_J, q_0);

            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible J_ST reached at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_J = XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP);
            prv_E = XMX(SP_E, q_0) + XSC(SP_E, SP_LOOP);

            /* verifies that scores agree with true previous state in alignment */
            if ( CMP_TOL( cur, prv_J ) ) {
               st_cur = J_ST;
            }
            else if ( CMP_TOL( cur, prv_E ) ) {
               st_cur = E_ST;
            }
            else {
               fprintf( stderr, "ERROR: Failed to alnace from J_ST at (%d,%d)\n", q_0, t_0);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         } break;

         default:
         {
            fprintf( stderr, "ERROR: Hit Bogus State!!!\n");
            ERRORCHECK_exit(EXIT_FAILURE);
         }
      }
      
      ALIGNMENT_Append( aln, st_cur, q_0, t_0 );

      /* For {N,C,J}, we deferred i decrement. */
      if ( (st_cur == N_ST || st_cur == J_ST || st_cur == C_ST) && (st_cur == st_prv) ) {
         q_0--;
      }

      /* Update previous state */
      st_prv = st_cur;
   }

   /* reverse order of traceback */
   ALIGNMENT_Reverse( aln );
   /* scan traceback for all begin, end states */
   ALIGNMENT_Find_Length( aln );

   #if DEBUG
   {
      MATRIX_2D* cloud_MX = debugger->cloud_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      for ( int i = 0; i < aln->traces->N; i++ ) 
      {
         tr = aln->traces->data;
         if ( tr[i].st == M_ST || tr[i].st == I_ST || tr[i].st == D_ST ) {
            MX_2D( cloud_MX, tr[i].q_0, tr[i].t_0 ) = -1.0;
         }
      }
   }
   #endif

   return STATUS_SUCCESS;
}


/* TODO: Finish implementataion */
/** FUNCTION:  run_Traceback_Quad_via_max()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Version 3: My implementation.  Takes maximum next step from current state to previous state.  
 *             Warning: No verification step.  
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Quad_via_max(     const SEQUENCE*     query,       /* query sequence */
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

   seq   = query->seq; 
   tr    = aln->traces->data; 

   /* local or global? */
   is_local = target->isLocal;

   /* allocate memory for trace */
   ALIGNMENT_Reuse( aln, Q, T );

   /* Backalnacing, so C is end state */
   ALIGNMENT_Append( aln, T_ST, q_0, t_0 );
   ALIGNMENT_Append( aln, C_ST, q_0, t_0 );
   st_prv = C_ST;

   /* End of trace is S state */
   while (st_prv != S_ST)
   {
      /* if we have reached the end of the query sequence */
      if (q_0 == 0) {
         ALIGNMENT_Append( aln, S_ST, q_0, t_0 );
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
               ERRORCHECK_exit(EXIT_FAILURE);
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
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M) + MSC(t_0, A);
            prv_M = MMX(q_1, t_1) + TSC(t_1, M2M) + MSC(t_0, A);
            prv_I = IMX(q_1, t_1) + TSC(t_1, I2M) + MSC(t_0, A);
            prv_D = DMX(q_1, t_1) + TSC(t_1, TM) + MSC(t_0, A);

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
               ERRORCHECK_exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_M = MMX(q_0, t_1) + TSC(t_1, M2D);
            prv_D = DMX(q_0, t_1) + TSC(t_1, TD);

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
               ERRORCHECK_exit(EXIT_FAILURE);
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
               ERRORCHECK_exit(EXIT_FAILURE);
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
               ERRORCHECK_exit(EXIT_FAILURE);
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
            ERRORCHECK_exit(EXIT_FAILURE);
         }
      }

      ALIGNMENT_Append( aln, st_cur, q_0, t_0 );

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
            MX_2D( cloud_MX, tr[i].q_0, tr[i].t_0 ) = -1.0;
      }
   }
   #endif

   return STATUS_SUCCESS;
}