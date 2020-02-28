/*******************************************************************************
 *
 *  FILE:    traceback.c
 *  PURPOSE: Traceback for Viterbi Algorithm (QUADRATIC SPACE).
 *
 *  AUTHOR:  Dave Rich
 *  BUG:     Lots.
 *
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* objects */
#include "objects/structs.h"
#include "objects/sequence.h"
#include "objects/hmm_profile.h"
#include "objects/alignment.h"

/* local imports */
#include "utility.h"
#include "hmm_parser.h"

/* header */
#include "traceback_quad.h"

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
 *             <st_MX>     Normal State (Match, Insert, Delete) Maalnix,
 *             <sp_MX>     Special State (J,N,B,C,E) Maalnix,
 *             <aln>       Traceback Alignment
 *
 *  RETURN:    No Return.
 */
void traceback_Build(const SEQUENCE*     query,
                     const HMM_PROFILE*  target,
                     const int           Q, 
                     const int           T,
                     float*              st_MX,
                     float*              sp_MX,
                     ALIGNMENT*          aln)
{
   int    i = Q;              /* position in query seq (row) (1...L) */
   int    j = 0;              /* position in target model (col) (1...M) */
   char   a;                  /* store current character in sequence */
   int    A;                  /* store int value of current character in sequence */
   int    st_cur, st_prv;     /* current, previous state in alnace */
   float  tol = 1e-5;         /* acceptable tolerance range for "equality tests" */
   char   *seq = query->seq;  /* alias for getting seq */

   /* local or global? */
   bool is_local = target->isLocal;

   /* allocate memory for alnace */
   const int min_size = 256;
   aln->N = 0;
   aln->Nalloc = min_size;
   aln->traces = malloc(min_size * sizeof(TRACE));

   /* Backalnacing, so C is end state */
   traceback_Append(aln, T_ST, i, j);
   traceback_Append(aln, C_ST, i, j);
   st_prv = C_ST;

   const char * states[] = {"ST_M",
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

   /* End of alnace is S state */
   while (st_prv != S_ST)
   {
      if (i == 0) 
      {
         traceback_Append(aln, S_ST, i, j);
         break;
      } 
      
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
               printf("ERROR: Failed to alnace from B_ST at (%d,%d)\n", i, j);
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
               st_cur = M_ST;    /* can't come from D, in a *local* Viterbi alnace. */
               for (j = T; j >= 1; j--) {
                  // printf("testing E at (%d, %d) => (%.2f v. %.2f) \n", i, j, XMX(SP_E, i), MMX(i, j) );
                  if ( CMP_TOL( XMX(SP_E, i), MMX(i, j) ) )
                     break;
               }
               if (j == 0) {
                  printf("ERROR: Failed to alnace from E_ST at (%d,%d)\n", i, j);
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
                  printf("ERROR: Failed to alnace from E_ST at (%d,%d)\n", i, j);
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
               printf("ERROR: Failed to alnace from M_ST at (%d,%d)\n", i, j);
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
               printf("ERROR: Failed to alnace from D_ST at (%d,%d)\n", i, j);
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
               printf("ERROR: Failed to alnace from I_ST at (%d,%d)\n", i, j);
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
               printf("ERROR: Failed to alnace from B_ST at (%d,%d)\n", i, j);
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
               printf("ERROR: Failed to alnace from J_ST at (%d,%d)\n", i, j);
               exit(1);
            }
            break;

         default:
            printf("ERROR: Hit Bogus State!!!\n");
            exit(1);
      }

      /* Add new state and (i,j) to alnace */
      int state_num[] = {MAT_ST, INS_ST, DEL_ST, SP_E, SP_N, SP_J, SP_C, SP_B, -1, -1};
      // if (st_cur < 3) {
      //    printf("%s:\t(%d,%d)  \t%f\n", states[st_cur], i, j, ST_MX(st_MX, state_num[st_cur], i, j) );
      // } 
      // else if (st_cur >= 3 && st_cur < 9) {
      //    printf("%s:\t(%d,%d)  \t%f\n", states[st_cur], i, j, SP_MX(sp_MX, state_num[st_cur], i) );
      // }
      // else {
      //    printf("%s:\t(%d,%d)\n", states[st_cur], i, j );
      // }

      traceback_Append(aln, st_cur, i, j);

      /* For NCJ, we deferred i decrement. */
      if ( (st_cur == N_ST || st_cur == J_ST || st_cur == C_ST) && st_cur == st_prv) {
         i--;
      }

      /* Update previous state */
      st_prv = st_cur;
   }

   /* reverse order of alnaceback */
   traceback_Reverse(aln);

   /* find end and begin alignment points (first and last match state) */
   for (i = 0; i < aln->N; ++i) {
      if (aln->traces[i].st == M_ST) {
         aln->beg = i;
         break;
      }
   }
   for (i = aln->N - 1; i >= 0; --i) {
      if (aln->traces[i].st == M_ST) {
         aln->end = i;
         break;
      }
   }

   return;
}

/*
 *  FUNCTION:  traceback_Append()
 *  SYNOPSIS:  Append next state to Optimal Alignment.
 *
 *  PURPOSE:
 *
 *  ARGS:      <aln>    Traceback Alignment
 *             <st>     state
 *             <i>      index in sequence
 *             <j>      index in model
 *
 *  RETURN:    No Return.
 */
void traceback_Append(ALIGNMENT*  aln,
                      int          st,
                      int          i,
                      int          j)
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
   int N = aln->N;
   int size = aln->Nalloc;
   // printf("[%d](%s,%d,%d) -> \n", N, states[st], i, j);

   /* grow alnace length if needed */
   if (aln->N >= aln->Nalloc - 1) {
      aln->traces = (TRACE *)realloc(aln->traces, aln->Nalloc * 2 * sizeof(TRACE));
   }

   /* jump from current state to the prev state */
   switch (st)
   {
      /* Emit-on-Transition States: */
      case N_ST:
      case C_ST:
      case J_ST:
         aln->traces[N].i = ( (aln->traces[N - 1].st == st) ? i : 0 );
         aln->traces[N].j = 0;
         break;

      /* Non-Emitting States, not in Main Model: */
      case X_ST:
      case S_ST:
      case B_ST:
      case E_ST:
      case T_ST:
         aln->traces[N].i = 0;
         aln->traces[N].j = 0;
         break;

      /* Non-Emitting States, but in Main Model: */
      case D_ST:
         aln->traces[N].i = 0;
         aln->traces[N].j = j;
         break;

      /* Emitting States: */
      case M_ST:
      case I_ST:
         aln->traces[N].i = i;
         aln->traces[N].j = j;
         break;

      default:
         printf("ERROR: Hit Bogus State!!!\n");
         exit(0);
   }

   aln->traces[N].st = st;
   aln->N++;
   return;
}


/*
 *  FUNCTION:  traceback_Reverse()
 *  SYNOPSIS:  Reverse Traceback from backwards to forwards.
 *
 *  PURPOSE:
 *
 *  ARGS:      <aln>    Traceback Alignment
 *
 *  RETURN:    No Return.
 */
void traceback_Reverse (ALIGNMENT* aln)
{
   TRACE *aln1, *aln2, tmp;
   int N = aln->N;
   int st_cur, st_nxt;

   /* */
   for (int i = 0; i < aln->N; i++) 
   {
      aln1 = &aln->traces[i];
      aln2 = &aln->traces[i+1];

      if (aln1->st == aln2->st && (aln1->st == N_ST || aln1->st == C_ST || aln1->st == J_ST) ) 
      {
         if (aln1->i == 0 && aln1->st > 0) 
         {
            aln1->i = aln2->i;
            aln2->i = 0;
         }
      }
   }

   /* reverse alnaceback inplace */
   for (int i = 0; i < (aln->N / 2); ++i) 
   {
      aln1 = &aln->traces[i];
      aln2 = &aln->traces[N - 1 - i];

      tmp.st = aln1->st;
      tmp.i = aln1->i;
      tmp.j = aln1->j;

      aln1->st = aln2->st;
      aln1->i = aln2->i;
      aln1->j = aln2->j;

      aln2->st = tmp.st;
      aln2->i = tmp.i;
      aln2->j = tmp.j;
   }
}