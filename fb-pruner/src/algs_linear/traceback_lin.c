/*******************************************************************************
 *  FILE:    traceback.c
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
#include "structs.h"
#include "utilities.h"
#include "objects.h"
#include "algs_quad.h"

/* header */
#include "traceback_lin.h"

/*
 *  FUNCTION:  run_Traceback_Quad()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Lin(     const SEQUENCE*     query,       /* query sequence */
                           const HMM_PROFILE*  target,      /* HMM model */
                           const int           Q,           /* query/seq length */
                           const int           T,           /* target/model length */
                           MATRIX_3D*          st_MX,       /* Normal State (Match, Insert, Delete) Maalnix */
                           MATRIX_2D*          sp_MX,       /* Special State (J,N,B,C,E) Maalnix */
                           ALIGNMENT*          aln )        /* Traceback Alignment */
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
   ALIGNMENT_Append( aln, tr, T_ST, q_0, t_0 );
   ALIGNMENT_Append( aln, tr, C_ST, q_0, t_0 );
   st_prv = C_ST;

   /* End of alnace is S state */
   while (st_prv != S_ST)
   {
      if (q_0 == 0) 
      {
         ALIGNMENT_Append( aln, tr, S_ST, q_0, t_0 );
         break;
      } 
      q_1 = q_0-1;
      t_1 = t_0-1;
      
      a = seq[q_1];
      A = AA_REV[a];

      /* jump from current state to the prev state */
      switch (st_prv)
      {
         /* C STATE to {C,E} */
         case C_ST:  /* C(i) comes from C(i-1) or E(i) */
            if (XMX(SP_C, q_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible C_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            if ( CMP_TOL( XMX(SP_C, q_0), XMX(SP_C, q_1) + XSC(SP_C, SP_LOOP) ) )
               st_cur = C_ST;
            else if ( CMP_TOL( XMX(SP_C, q_0), XMX(SP_E, q_0) + XSC(SP_E, SP_MOVE) ) )
               st_cur = E_ST;
            else {
               fprintf( stderr, "ERROR: Failed to trace from B_ST at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }
            break;

         /* E STATE to {M,D} */
         case E_ST:  /* E connects from any M state. k set here */
            if (XMX(SP_E, q_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible E_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
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
                  exit(EXIT_FAILURE);
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
                  exit(EXIT_FAILURE);
               }
            }
            break;

         /* M STATE to {B,M,I,D} */
         case M_ST:  /* M connects from i-1,k-1, or B */
            if (MMX(q_0, t_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible M_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            if ( CMP_TOL( MMX(q_0, t_0), XMX(SP_B, q_1) + TSC(t_1, B2M) + MSC(t_0, A) ) )
               st_cur = B_ST;
            else if ( CMP_TOL( MMX(q_0, t_0), MMX(q_1, t_1) + TSC(t_1, M2M) + MSC(t_0, A) ) )
               st_cur = M_ST;
            else if ( CMP_TOL( MMX(q_0, t_0), IMX(q_1, t_1) + TSC(t_1, I2M) + MSC(t_0, A) ) )
               st_cur = I_ST;
            else if ( CMP_TOL( MMX(q_0, t_0), DMX(q_1, t_1) + TSC(t_1, D2M) + MSC(t_0, A) ) )
               st_cur = D_ST;
            else {
               fprintf( stderr, "ERROR: Failed to trace from M_ST at (%d,%d)\n", t_0, q_0);
               fprintf( stderr, "TOL: %f vs %f\n", MMX(q_0, t_0), MMX(q_1, t_1) + TSC(t_1, D2M) + MSC(t_0, A) );
               exit(EXIT_FAILURE);
            }
            t_0--; q_0--;
            break;

         /* D STATE to {M,D} */
         case D_ST:  /* D connects from M,D at i,k-1 */
            if (DMX(q_0, t_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible D_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            t_1 = t_0-1;
            if ( CMP_TOL( DMX(q_0, t_0), MMX(q_0, t_1) + TSC(t_1, M2D) ) )
               st_cur = M_ST;
            else if ( CMP_TOL( DMX(q_0, t_0), DMX(q_0, t_1) + TSC(t_1, D2D) ) )
               st_cur = D_ST;
            else {
               fprintf( stderr, "ERROR: Failed to alnace from D_ST at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }
            t_0--;
            break;

         /* I STATE to {M,I} */
         case I_ST:  /* I connects from M,I at i-1,k */
            if (IMX(q_0, t_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible I_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            if ( CMP_TOL( IMX(q_0, t_0), MMX(q_1, t_0) + TSC(t_0, M2I) + ISC(t_0, A) ) )
               st_cur = M_ST;
            else if ( CMP_TOL( IMX(q_0, t_0), IMX(q_1, t_0) + TSC(t_0, I2I) + ISC(t_0, A) ) )
               st_cur = I_ST;
            else {
               fprintf( stderr, "ERROR: Failed to alnace from I_ST at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }
            q_0--;
            break;

         /* N STATE to {N,S} */
         case N_ST:  /* N connects from S, N */
            if (XMX(SP_N, q_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible N_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
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
               exit(EXIT_FAILURE);
            }
            break;

         /* J STATE to {J,E} */
         case J_ST:  /* J connects from E(i) or J(i-1) */
            if (XMX(SP_J, q_0) == -INF ) {
               fprintf( stderr, "ERROR: Impossible J_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            if ( CMP_TOL( XMX(SP_J, q_0), XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP) ) )
               st_cur = J_ST;
            else if ( CMP_TOL( XMX(SP_J, q_0), XMX(SP_E, q_0) + XSC(SP_E, SP_LOOP) ) )
               st_cur = E_ST;
            else {
               fprintf( stderr, "ERROR: Failed to alnace from J_ST at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }
            break;

         default:
            fprintf( stderr, "ERROR: Hit Bogus State!!!\n");
            exit(EXIT_FAILURE);
      }

      ALIGNMENT_Append( aln, tr, st_cur, q_0, t_0 );

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
         if ( tr[i].st == M_ST || tr[i].st == I_ST || tr[i].st == D_ST ) {
            MX_2D( cloud_MX, tr[i].i, tr[i].j ) = -1.0;
         }
      }
   }
   #endif

   return STATUS_SUCCESS;
}
