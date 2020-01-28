/*******************************************************************************
 *  @file viterbi.c
 *  @brief The Viterbi Algorithm for Sequence Alignment Search.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

/* macros */
#define getName(var) #var
#define SCALE_FACTOR 1000

/* external imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports (after struct declarations) */
#include "structs.h"
#include "misc.h"
#include "hmm_parser.h"
#include "viterbi.h"

/*
 *  FUNCTION:  viterbi_Run()
 *  SYNOPSIS:  Run Viterbi Algorithm (Seq-to-Profile, general unoptimized)
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence,
 *             <target>    HMM model,
 *             <res>       Results Data
 *             <tr>        Traceback Data
 *
 *  RETURN:
 */
float viterbi_Run (const SEQ* query,
                   const HMM_PROFILE* target,
                   int Q, int T,
                   float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                   float sp_MX[ NUM_SPECIAL_STATES * (Q + 1) ],
                   RESULTS* res,
                   TRACEBACK* tr)
{
   char   a;           /* store current character in sequence */
   int    A;           /* store int value of character */
   int    i, j, k = 0; /* row, column indices */
   char   *seq  = query->seq; /* alias for getting seq */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_best;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc1, sc2, sc3, sc4;

   /* local or global (multiple alignments) */
   bool   isLocal = target->isLocal;
   float  sc_E = (isLocal) ? 0 : -INF;

   /* initialize special states (?) */
   XMX(SP_N, 0) = 0;                                        /* S->N, p=1             */
   XMX(SP_B, 0) = XSC(SP_N, SP_MOVE);                       /* S->N->B, no N-tail    */
   XMX(SP_E, 0) = XMX(SP_C, 0) = XMX(SP_J, 0) = -INF;       /* need seq to get here  */

   /* initialize zero row (top-edge) */
   for (j = 0; j <= T; j++)
   { MMX(0, j) = IMX(0, j) = DMX(0, j) = -INF; }         /* need seq to get here  */

   /* FOR every position in QUERY seq */
   for (i = 1; i <= Q; i++)
   {
      /* Get next character in Query */
      a = seq[i - 1];
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX(i, 0) = IMX(i, 0) = DMX(i, 0) = -INF;
      XMX(SP_E, i) = -INF;

      /* FOR every position in TARGET profile */
      for (j = 1; j < T; j++)
      {
         /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         sc1 = prev_mat = MMX(i - 1, j - 1)  + TSC(j - 1, M2M);
         sc2 = prev_ins = IMX(i - 1, j - 1)  + TSC(j - 1, I2M);
         sc3 = prev_del = DMX(i - 1, j - 1)  + TSC(j - 1, D2M);
         sc4 = prev_beg = XMX(SP_B, i - 1) + TSC(j - 1, B2M); /* from begin match state (new alignment) */
         /* best-to-match */
         prev_best = calc_Max( calc_Max( prev_mat, prev_ins ) , calc_Max( prev_del, prev_beg ) );
         MMX(i, j)  = prev_best + MSC(j, A);

         // printf("(%d, %d)  \tMM: %.3f\t IM: %.3f\t DM: %.3f\t BM: %.3f\t MSC: %.3f\n", i, j, sc1, sc2, sc3, sc4, MSC(j, A));

         /* FIND BEST PATH TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX(i - 1, j) + TSC(j, M2I);
         prev_ins = IMX(i - 1, j) + TSC(j, I2I);
         /* best-to-insert */
         prev_best = calc_Max(prev_mat, prev_ins);
         IMX(i, j) = prev_best + ISC(j, A);

         /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX(i, j - 1) + TSC(j - 1, M2D);
         prev_del = DMX(i, j - 1) + TSC(j - 1, D2D);
         /* best-to-delete */
         prev_best = calc_Max(prev_mat, prev_del);
         DMX(i, j) = prev_best;

         /* UPDATE E STATE */
         sc1 = XMX(SP_E, i);
         sc2 = MMX(i, j) + sc_E;
         sc4 = XMX(SP_E, i) = calc_Max( XMX(SP_E, i), MMX(i, j) + sc_E );

         // printf("(%d, %d)  \tEE: %.3f\t ME: %.3f\t E_sc: %.3f\t E_MAX: %.3f\n", i,j, sc1, sc2, sc_E, sc4);

         /* (?!) check for best score seen so far (best score should end in a match state) */
         if (sc_max < MMX(i, j))
         {
            tr->sc_last_m = MMX(i, j);
            tr->last_m.i = i;
            tr->last_m.j = j;
         }

      }

      /* UNROLLED FINAL LOOP ITERATION */
      j = T;

      /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
      /* best previous state transition (match takes the diag element of each prev state) */
      sc1 = prev_mat = MMX(i - 1, j - 1)  + TSC(j - 1, M2M);
      sc2 = prev_ins = IMX(i - 1, j - 1)  + TSC(j - 1, I2M);
      sc3 = prev_del = DMX(i - 1, j - 1)  + TSC(j - 1, D2M);
      sc4 = prev_beg = XMX(SP_B, i - 1) + TSC(j - 1, B2M); /* from begin match state (new alignment) */
      /* best-to-match */
      prev_best = calc_Max(
                     calc_Max( prev_mat, prev_ins ),
                     calc_Max( prev_del, prev_beg )
                  );
      MMX(i, j) = prev_best + MSC(j, A);

      // DAVID RICH EDIT
      // printf("(%d, %d)  \tMM: %.3f\t IM: %.3f\t DM: %.3f\t BM: %.3f\t MSC: %.3f\n", i, j, sc1, sc2, sc3, sc4, MSC(j, A));

      /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
      /* previous states (match takes the left element of each state) */
      prev_mat = MMX(i, j - 1) + TSC(j - 1, M2D);
      prev_del = DMX(i, j - 1) + TSC(j - 1, D2D);
      /* best-to-delete */
      prev_best = calc_Max( prev_mat, prev_del );
      DMX(i, j) = prev_best;

      /* UPDATE E STATE */
      sc1 = XMX(SP_E, i);
      sc2 = MMX(i, j);
      sc3 = DMX(i, j);
      sc4 = XMX(SP_E, i) = calc_Max( sc1, calc_Max( sc2, sc3 ) );

      // printf("(%d, %d) \t EE: %.3f\t ME: %.3f\t DE: %.3f\t E_MAX: %.3f\n", i, j, sc1, sc2, sc3, sc4);

      /* SPECIAL STATES */
      /* J state */
      sc_1 = XMX(SP_J, i - 1) + XSC(SP_J, SP_LOOP); /* J->J */
      sc_2 = XMX(SP_E, i)   + XSC(SP_E, SP_LOOP);   /* E->J is E's "loop" */
      XMX(SP_J, i) = calc_Max( sc_1, sc_2 );

      /* C state */
      sc_1 = XMX(SP_C, i - 1) + XSC(SP_C, SP_LOOP);
      sc_2 = XMX(SP_E, i)   + XSC(SP_E, SP_MOVE);
      XMX(SP_C, i) = calc_Max( sc_1, sc_2 );

      /* N state */
      XMX(SP_N, i) = XMX(SP_N, i - 1) + XSC(SP_N, SP_LOOP);

      /* B state */
      sc_1 = XMX(SP_N, i) + XSC(SP_N, SP_MOVE);     /* N->B is N's move */
      sc_2 = XMX(SP_J, i) + XSC(SP_J, SP_MOVE);     /* J->B is J's move */
      XMX(SP_B, i) = calc_Max( sc_1, sc_2 );
   }

   /* T state */
   sc_best = XMX(SP_C, Q) + XSC(SP_C, SP_MOVE);

   tr->sc_max = sc_best;
   return sc_best;
}

/*
 *  FUNCTION:  traceback_Build()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence,
 *             <target>    HMM model,
 *             <Q>         query/seq length,
 *             <T>         target/model length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <tr>        Traceback Alignment
 *
 *  RETURN:
 */
void traceback_Build (const SEQ* query,
                      const HMM_PROFILE* target,
                      int Q, int T,
                      float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                      float sp_MX[ NUM_SPECIAL_STATES * (Q + 1) ],
                      TRACEBACK* tr)
{
   int    i = Q;              /* position in query seq (row) (1...L) */
   int    j = 0;              /* position in target model (col) (1...M) */
   char   a;                  /* store current character in sequence */
   int    A;                  /* store int value of character */
   int    st_cur, st_prv;    /* current, previous state in trace */
   float  tol = 1e-5;         /* acceptable tolerance for "equality" tests */
   char   *seq = query->seq;  /* alias for getting seq */

   /* local or global? */
   bool is_local = target->isLocal;

   /* allocate memory for trace */
   const int min_size = 256;
   tr->N = 0;
   tr->size = min_size;
   tr->traces = malloc(min_size * sizeof(TRACE));

   /* Backtracing, so C is end state */
   traceback_Append(tr, T_ST, i, j);
   traceback_Append(tr, C_ST, i, j);
   st_prv = C_ST;

   /* End of trace is S state */
   while (st_prv != S_ST)
   {
      a = seq[i - 1];
      A = AA_REV[a];

      /* jump from current state to the prev state */
      switch (st_prv)
      {
         /* C STATE to {C,E} */
         case C_ST:  /* C(i) comes from C(i-1) or E(i) */
            if (XMX(SP_C, i) == -INF ) {
               printf("ERROR: Impossible C_ST reached at (%d,%d)\n", i, j);
               exit(1);
            }

            if ( CMP_TOL( XMX(SP_C, i), XMX(SP_C, i - 1) + XSC(SP_C, SP_LOOP) ) )
               st_cur = C_ST;
            else if ( CMP_TOL( XMX(SP_C, i), XMX(SP_E, i) + XSC(SP_E, SP_MOVE) ) )
               st_cur = E_ST;
            else {
               printf("ERROR: Failed to trace from B_ST at (%d,%d)\n", i, j);
               exit(1);
            }
            break;

         /* E STATE to {M,D} */
         case E_ST:  /* E connects from any M state. k set here */
            if (XMX(SP_E, i) == -INF ) {
               printf("ERROR: Impossible E_ST reached at (%d,%d)\n", i, j);
               exit(1);
            }

            if ( is_local )  /* local mode: ends in M */
            {
               st_cur = M_ST;    /* can't come from D, in a *local* Viterbi trace. */
               for (j = T; j >= 1; j--) {
                  // printf("testing E at (%d, %d) => (%.2f v. %.2f) \n", i, j, XMX(SP_E, i), MMX(i, j) );
                  if ( CMP_TOL( XMX(SP_E, i), MMX(i, j) ) )
                     break;
               }
               if (j == 0) {
                  printf("ERROR: Failed to trace from E_ST at (%d,%d)\n", i, j);
                  exit(1);
               }
            }
            else     /* glocal mode: we either come from D_M or M_M */
            {
               if ( CMP_TOL( XMX(SP_E, i), MMX(i, T) ) ) {
                  st_cur = M_ST;
                  j = T;
               }
               else if ( CMP_TOL( XMX(SP_E, i), DMX(i, T) ) ) {
                  st_cur = D_ST;
                  j = T;
               }
               else {
                  printf("ERROR: Failed to trace from E_ST at (%d,%d)\n", i, j);
                  exit(1);
               }
            }
            break;

         /* M STATE to {B,M,I,D} */
         case M_ST:  /* M connects from i-1,k-1, or B */
            if (MMX(i, j) == -INF ) {
               printf("ERROR: Impossible M_ST reached at (%d,%d)\n", i, j);
               exit(1);
            }

            if ( CMP_TOL( MMX(i, j), XMX(SP_B, i - 1) + TSC(j - 1, B2M) + MSC(j, A) ) )
               st_cur = B_ST;
            else if ( CMP_TOL( MMX(i, j), MMX(i - 1, j - 1) + TSC(j - 1, M2M) + MSC(j, A) ) )
               st_cur = M_ST;
            else if ( CMP_TOL( MMX(i, j), IMX(i - 1, j - 1) + TSC(j - 1, I2M) + MSC(j, A) ) )
               st_cur = I_ST;
            else if ( CMP_TOL( MMX(i, j), DMX(i - 1, j - 1) + TSC(j - 1, D2M) + MSC(j, A) ) )
               st_cur = D_ST;
            else {
               printf("ERROR: Failed to trace from M_ST at (%d,%d)\n", i, j);
               printf("TOL: %f vs %f\n", MMX(i, j), MMX(i - 1, j - 1) + TSC(j - 1, D2M) + MSC(j, A) );
               exit(1);
            }
            j--; i--;
            break;

         /* D STATE to {M,D} */
         case D_ST:  /* D connects from M,D at i,k-1 */
            if (DMX(i, j) == -INF ) {
               printf("ERROR: Impossible D_ST reached at (%d,%d)\n", i, j);
               exit(1);
            }

            if ( CMP_TOL( DMX(i, j), MMX(i, j - 1) + TSC(j - 1, M2D) ) )
               st_cur = M_ST;
            else if ( CMP_TOL( DMX(i, j), DMX(i, j - 1) + TSC(j - 1, D2D) ) )
               st_cur = D_ST;
            else {
               printf("ERROR: Failed to trace from D_ST at (%d,%d)\n", i, j);
               exit(1);
            }
            j--;
            break;

         /* I STATE to {M,I} */
         case I_ST:  /* I connects from M,I at i-1,k */
            if (IMX(i, j) == -INF ) {
               printf("ERROR: Impossible I_ST reached at (%d,%d)\n", i, j);
               exit(1);
            }

            if ( CMP_TOL( IMX(i, j), MMX(i - 1, j) + TSC(j, M2I) + ISC(j, A) ) )
               st_cur = M_ST;
            else if ( CMP_TOL( IMX(i, j), IMX(i - 1, j) + TSC(j, I2I) + ISC(j, A) ) )
               st_cur = I_ST;
            else {
               printf("ERROR: Failed to trace from I_ST at (%d,%d)\n", i, j);
               exit(1);
            }
            i--;
            break;

         /* N STATE to {N,S} */
         case N_ST:  /* N connects from S, N */
            if (XMX(SP_N, i) == -INF ) {
               printf("ERROR: Impossible N_ST reached at (%d,%d)\n", i, j);
               exit(1);
            }

            st_cur = ( (i <= 0) ? S_ST : N_ST );
            break;

         /* B STATE to {N,J} */
         case B_ST:  /* B connects from N, J */
            if ( CMP_TOL( XMX(SP_B, i), XMX(SP_N, i) + XSC(SP_N, SP_MOVE) ) )
               st_cur = N_ST;
            else if ( CMP_TOL( XMX(SP_B, i), XMX(SP_J, i) + XSC(SP_J, SP_MOVE) ) )
               st_cur = J_ST;
            else {
               printf("ERROR: Failed to trace from B_ST at (%d,%d)\n", i, j);
               exit(1);
            }
            break;

         /* J STATE to {J,E} */
         case J_ST:  /* J connects from E(i) or J(i-1) */
            if (XMX(SP_J, i) == -INF ) {
               printf("ERROR: Impossible J_ST reached at (%d,%d)\n", i, j);
               exit(1);
            }

            if ( CMP_TOL( XMX(SP_J, i), XMX(SP_J, i - 1) + XSC(SP_J, SP_LOOP) ) )
               st_cur = J_ST;
            else if ( CMP_TOL( XMX(SP_J, i), XMX(SP_E, i) + XSC(SP_E, SP_LOOP) ) )
               st_cur = E_ST;
            else {
               printf("ERROR: Failed to trace from J_ST at (%d,%d)\n", i, j);
               exit(1);
            }
            break;

         default:
            printf("ERROR: Hit Bogus State!!!\n");
            exit(1);
      }

      /* Add new state and (i,j) to trace */
      int state_num[] = {MAT_ST, INS_ST, DEL_ST, SP_E, SP_N, SP_J, SP_C, SP_B, -1, -1};
      // static char * states[] = {"ST_M",
      //                           "ST_I",
      //                           "ST_D",
      //                           "ST_E",
      //                           "ST_N",
      //                           "ST_J",
      //                           "ST_C",
      //                           "ST_B",
      //                           "ST_S",
      //                           "ST_T",
      //                           "ST_X" };
      // if (st_cur < 3) {
      //    printf("%s:\t(%d,%d)  \t%f\n", states[st_cur], i, j, ST_MX(st_MX, state_num[st_cur], i, j) );
      // } 
      // else if (st_cur >= 3 && st_cur < 9) {
      //    printf("%s:\t(%d,%d)  \t%f\n", states[st_cur], i, j, SP_MX(sp_MX, state_num[st_cur], i) );
      // }
      // else {
      //    printf("%s:\t(%d,%d)\n", states[st_cur], i, j );
      // }

      traceback_Append(tr, st_cur, i, j);

      /* For NCJ, we deferred i decrement. */
      if ( (st_cur == N_ST || st_cur == J_ST || st_cur == C_ST) && st_cur == st_prv) {
         i--;
      }

      /* Update previous state */
      st_prv = st_cur;
   }

   /* reverse order of traceback */
   traceback_Reverse(tr);

   /* find end and begin alignment points (first and last match state) */
   for (i = 0; i < tr->N; ++i) {
      if (tr->traces[i].st == M_ST) {
         tr->end = i;
         tr->first_m.i = tr->traces[i].i;
         tr->first_m.j = tr->traces[i].j;
         break;
      }
   }
   for (i = tr->N - 1; i >= 0; --i) {
      if (tr->traces[i].st == M_ST) {
         tr->start = i;
         tr->last_m.i = tr->traces[i].i;
         tr->last_m.j = tr->traces[i].j;
         break;
      }
   }

   return;
}


/*
 *  FUNCTION:  traceback_Reverse()
 *  SYNOPSIS:  Reverse Traceback from backwards to forwards.
 *
 *  PURPOSE:
 *
 *  ARGS:      <tr> traceback object
 *
 *  RETURN:
 */
void traceback_Reverse (TRACEBACK* tr)
{
   TRACE *tr1, *tr2, tmp;
   int N = tr->N;
   int st_cur, st_nxt;

   /* */
   for (int i = 0; i < tr->N; i++) 
   {
      tr1 = &tr->traces[i];
      tr2 = &tr->traces[i+1];

      if (tr1->st == tr2->st && (tr1->st == N_ST || tr1->st == C_ST || tr1->st == J_ST) ) 
      {
         if (tr1->i == 0 && tr1->st > 0) 
         {
            tr1->i = tr2->i;
            tr2->i = 0;
         }
      }
   }

   /* reverse traceback inplace */
   for (int i = 0; i < (tr->N / 2); ++i) 
   {
      tr1 = &tr->traces[i];
      tr2 = &tr->traces[N - 1 - i];

      tmp.st = tr1->st;
      tmp.i = tr1->i;
      tmp.j = tr1->j;

      tr1->st = tr2->st;
      tr1->i = tr2->i;
      tr1->j = tr2->j;

      tr2->st = tmp.st;
      tr2->i = tmp.i;
      tr2->j = tmp.j;
   }
}


/*
 *  FUNCTION:  traceback_Append()
 *  SYNOPSIS:  Append next state to Optimal Alignment.
 *
 *  PURPOSE:
 *
 *  ARGS:      <tr> traceback object
 *             <st> state
 *             <i>  index in sequence
 *             <j>  index in model
 *
 *  RETURN:
 */
void traceback_Append (TRACEBACK* tr,
                       int st,
                       int i,
                       int j)
{
   static char * states[] = {"ST_M",
                             "ST_I",
                             "ST_D",
                             "ST_E",
                             "ST_N",
                             "ST_J",
                             "ST_C",
                             "ST_B",
                             "ST_S",
                             "ST_T",
                             "ST_X" };
   int N = tr->N;
   int size = tr->size;
   // printf("[%d](%s,%d,%d) -> \n", N, states[st], i, j);

   /* grow trace length if needed */
   if (tr->N >= tr->size - 1) {
      tr->traces = (TRACE *)realloc(tr->traces, tr->size * 2 * sizeof(TRACE));
   }

   /* jump from current state to the prev state */
   switch (st)
   {
      /* Emit-on-Transition States: */
      case N_ST:
      case C_ST:
      case J_ST:
         tr->traces[N].i = ( (tr->traces[N - 1].st == st) ? i : 0 );
         tr->traces[N].j = 0;
         break;

      /* Non-Emitting States, not in Main Model: */
      case X_ST:
      case S_ST:
      case B_ST:
      case E_ST:
      case T_ST:
         tr->traces[N].i = 0;
         tr->traces[N].j = 0;
         break;

      /* Non-Emitting States, but in Main Model: */
      case D_ST:
         tr->traces[N].i = 0;
         tr->traces[N].j = j;
         break;

      /* Emitting States: */
      case M_ST:
      case I_ST:
         tr->traces[N].i = i;
         tr->traces[N].j = j;
         break;

      default:
         printf("ERROR: Hit Bogus State!!!\n");
         exit(0);
   }

   tr->traces[N].st = st;
   tr->N++;
   return;
}


/*
 *  FUNCTION:  traceback_Show()
 *  SYNOPSIS:  Build path into dynamic programming matrix according to traceback (1-index).
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix
 *             <tr>        Trace
 *
 *  RETURN:
 */
void traceback_Show (const int Q, const int T,
                     float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                     float sp_MX[ NUM_SPECIAL_STATES * (Q + 1) ],
                     TRACEBACK *tr)
{
   int          i, j, k, m;
   int          st;
   int          N = tr->N;

   /* Zero out Matrices */
   for (i = 0; i <= Q; i++)
   {
      for (j = 0; j <= T; j++)
      {
         MMX(i, j) = IMX(i, j) = DMX(i, j) = 0;
      }
      XMX(SP_N, i) = 0;
      XMX(SP_J, i) = 0;
      XMX(SP_E, i) = 0;
      XMX(SP_C, i) = 0;
      XMX(SP_B, i) = 0;
   }

   /*
   p7T_BOGUS =  0,
   p7T_M     =  1,
   p7T_D     =  2,
   p7T_I     =  3,
   p7T_S     =  4,
   p7T_N     =  5,
   p7T_B     =  6,
   p7T_E     =  7,
   p7T_C     =  8,
   p7T_T     =  9,
   p7T_J     = 10,
   p7T_X     = 11,
   */

   /* Input vals along the trace */
   for (k = 0; k < N; k++)
   {
      st = tr->traces[k].st;
      i = tr->traces[k].i;
      j = tr->traces[k].j;
      m = 1;
      // printf("%d/%d: %d(%d,%d) -> ", k, N, st, i, j);

      switch (st) {
      case M_ST:
         MMX(i, j) = m;
         break;
      case I_ST:
         MMX(i, j) = m;
         break;
      case D_ST:
         MMX(i, j) = m;
         break;
      case N_ST:
         XMX(SP_N, i) = m;
         break;
      case B_ST:
         XMX(SP_B, i) = m;
         break;
      case E_ST:
         XMX(SP_E, i) = m;
         break;
      case C_ST:
         XMX(SP_C, i) = m;
         break;
      case J_ST:
         XMX(SP_J, i) = m;
         break;
      default:
         break;
      }
   }
}




/*
 *  FUNCTION: traceback_Print()
 *  SYNOPSIS: Print Traceback object
 *
 *  PURPOSE:
 *
 *  ARGS:      <tr>      TRACEBACK Object
 *
 *  RETURN:
 */
void traceback_Print(TRACEBACK *tr)
{
   printf("traceback: \n");
   static char * states[] = {"ST_M",
                             "ST_I",
                             "ST_D",
                             "ST_E",
                             "ST_N",
                             "ST_J",
                             "ST_C",
                             "ST_B",
                             "ST_S",
                             "ST_T",
                             "ST_X" };

   // printf("LENGTH: %d / START: [%d](%d, %d) / END: [%d](%d, %d)\n", tr->N, tr->start, tr->first_m.i, tr->first_m.j, tr->end, tr->last_m.i, tr->last_m.j);
   
   for (unsigned int i = 0; i < tr->N; ++i)
   {
      printf("[%d](%s,%d,%d)-> \n", i, states[tr->traces[i].st], tr->traces[i].i, tr->traces[i].j);
   }
   printf("\n");
}

/*
 *  FUNCTION: traceback_Print()
 *  SYNOPSIS: Print Traceback object
 *
 *  PURPOSE:
 *
 *  ARGS:      <tr>      TRACEBACK Object
 *             <f>       Filename
 *
 *  RETURN:
 */
void traceback_Save(TRACEBACK *tr,
                    const char *_filename_)
{
   FILE *fp;
   fp = fopen(_filename_, "w");

   static char * states[] = {"ST_M",
                             "ST_I",
                             "ST_D",
                             "ST_E",
                             "ST_N",
                             "ST_J",
                             "ST_C",
                             "ST_B",
                             "ST_S",
                             "ST_T",
                             "ST_X" };

   // printf("LENGTH: %d / START: [%d](%d, %d) / END: [%d](%d, %d)\n", tr->N, tr->start, tr->first_m.i, tr->first_m.j, tr->end, tr->last_m.i, tr->last_m.j);
   
   for (unsigned int i = 0; i < tr->N; ++i)
   {
      int st = tr->traces[i].st;
      fprintf(fp, "[%d](%s,%d,%d)\n", i, states[st], tr->traces[i].i, tr->traces[i].j);

      // if (st == M_ST || st == I_ST || st == D_ST) {
      //    fprintf(fp, "[%d](%s,%d,%d) \n", i, states[st], tr->traces[i].i, tr->traces[i].j);
      // } else {
      //    fprintf(fp, "[%d](%s) \n", i, states[st]);
      // }
   }

   fprintf(fp, "\n");
   fclose(fp);
}
