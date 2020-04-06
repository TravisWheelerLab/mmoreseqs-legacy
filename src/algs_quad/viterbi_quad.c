/*******************************************************************************
 *  @file viterbi.c
 *  @brief The Viterbi Algorithm and Traceback for Sequence Alignment Search.
 *
 *  @author Dave Rich
 *  @bug Lots.
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
#include "objects/hmm_profile.h"
#include "objects/sequence.h"

/* local imports */
#include "utilities/utility.h"
#include "hmm_parser.h"

/* header */
#include "viterbi_quad.h"

/*
 *  FUNCTION:  viterbi_Run()
 *  SYNOPSIS:  Run Viterbi Algorithm (Seq-to-Profile, general unoptimized)
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence,
 *             <target>    HMM model,
 *             <Q>         Query length,
 *             <T>         Target length,
 *             <st_MX>     Normal state matrix,
 *             <sp_MX>     Special state matrix
 *
 *  RETURN:    No Return.
 */
float viterbi_Quad(const SEQUENCE*    query,
                   const HMM_PROFILE* target,
                   const int          Q, 
                   const int          T,
                   float*             st_MX,
                   float*             sp_MX,
                   float*             sc_final)
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
   *sc_final = sc_best;
   return sc_best;
}