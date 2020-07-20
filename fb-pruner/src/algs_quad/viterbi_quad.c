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

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "objects.h"
#include "algs_quad.h"

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
int viterbi_Quad(    const SEQUENCE*    query,
                     const HMM_PROFILE* target,
                     const int          Q, 
                     const int          T, 
                     MATRIX_3D*         st_MX,
                     MATRIX_2D*         sp_MX,
                     float*             sc_final )
{
   char     a;          /* store current character in sequence */
   int      A;          /* store int value of character */
   int      i,j,k = 0;  /* row, column indices */
   char*    seq;        /* alias for getting seq */

   float    prev_mat, prev_del, prev_ins; 
   float    prev_beg, prev_end, prev_sum;
   float    prev_best;
   float    sc, sc_1, sc_2;
   float    sc_max, sc_best;
   float    sc1, sc2, sc3, sc4;
   float    sc_M, sc_I, sc_D;
   int      x_0, x_1;
   int      r_0, r_1;                 /* d (mod 3) for assigning prev array ptrs */
   int      j_0, j_1;

   /* local or global (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_3D*  test_MX;

   /* initialize debugging data tools */
   #if DEBUG
   {
      cloud_MX = debugger->cloud_MX;
      test_MX  = debugger->test_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
   }
   #endif   

   /* --------------------------------------------------------------------------------- */

   /* query sequence */
   seq = query->seq;

   x_0 = 0;
   r_0 = x_0;

   /* initialize special states (?) */
   XMX(SP_N, x_0) = 0;                                        /* S->N, p=1             */
   XMX(SP_B, x_0) = XSC(SP_N, SP_MOVE);                       /* S->N->B, no N-tail    */
   XMX(SP_E, x_0) = XMX(SP_C, x_0) = XMX(SP_J, x_0) = -INF;   /* need seq to get here  */

   /* initialize zero row (top-edge) */
   for (j = 0; j <= T; j++) { 
      j_0 = j;
      MMX(r_0, j_0) = IMX(r_0, j_0) = DMX(r_0, j_0) = -INF;          /* need seq to get here  */
   }

   /* FOR every position in QUERY seq */
   for (i = 1; i <= Q; i++)
   {
      x_0 = i;
      x_1 = i-1;
      r_0 = x_0;
      r_1 = x_1;

      j_0 = 0;

      /* Get next character in Query */
      a = seq[i-1];
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX(r_0, j_0) = IMX(r_0, j_0) = DMX(r_0, j_0) = -INF;
      XMX(SP_E, x_0) = -INF;

      /* FOR every position in TARGET profile */
      for (j = 1; j < T; j++)
      {
         j_0 = j;
         j_1 = j-1;

         /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         prev_mat = MMX(r_1, j_1) + TSC(j_1, M2M);
         prev_ins = IMX(r_1, j_1) + TSC(j_1, I2M);
         prev_del = DMX(r_1, j_1) + TSC(j_1, D2M);
         prev_beg = XMX(SP_B, x_1) + TSC(j_1, B2M); /* from begin match state (new alignment) */

         /* best-to-match */
         prev_best = calc_Max( 
                        calc_Max( prev_mat, prev_ins ), 
                        calc_Max( prev_del, prev_beg ) );
         MMX(r_0, j_0)  = prev_best + MSC(j, A);

         /* FIND BEST PATH TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX(r_1, j_0) + TSC(j_0, M2I);
         prev_ins = IMX(r_1, j_0) + TSC(j_0, I2I);
         /* best-to-insert */
         prev_best = calc_Max(prev_mat, prev_ins);
         IMX(r_0, j_0) = prev_best + ISC(j_0, A);

         /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX(r_0, j_1) + TSC(j_1, M2D);
         prev_del = DMX(r_0, j_1) + TSC(j_1, D2D);
         /* best-to-delete */
         prev_best = calc_Max(prev_mat, prev_del);
         DMX(r_0, j_0) = prev_best;

         /* UPDATE E STATE */
         sc1 = XMX(SP_E, x_0);
         sc2 = MMX(r_0, j_0) + sc_E;
         sc4 = XMX(SP_E, x_0) = calc_Max( XMX(SP_E, x_0), MMX(r_0, j_0) + sc_E );

         /* embed linear row into quadratic test matrix */
         #if DEBUG
         {
            MX_2D(cloud_MX, x_0, j) = 1.0;
            MX_3D(test_MX, MAT_ST, x_0, j_0) = MMX(r_0, j_0);
            MX_3D(test_MX, INS_ST, x_0, j_0) = IMX(r_0, j_0);
            MX_3D(test_MX, DEL_ST, x_0, j_0) = DMX(r_0, j_0);
         }
         #endif
      }

      /* UNROLLED FINAL LOOP ITERATION */
      j = T;
      j_0 = j;
      j_1 = j-1;

      /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
      /* best previous state transition (match takes the diag element of each prev state) */
      prev_mat = MMX(r_1, j_1)  + TSC(j_1, M2M);
      prev_ins = IMX(r_1, j_1)  + TSC(j_1, I2M);
      prev_del = DMX(r_1, j_1)  + TSC(j_1, D2M);
      prev_beg = XMX(SP_B, x_1)  + TSC(j_1, B2M); /* from begin match state (new alignment) */
      /* best-to-match */
      prev_best = calc_Max(
                     calc_Max( prev_mat, prev_ins ),
                     calc_Max( prev_del, prev_beg ) );
      MMX(r_0, j_0) = prev_best + MSC(j_0, A);

      /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
      /* previous states (match takes the left element of each state) */
      prev_mat = MMX(r_0, j_1) + TSC(j_1, M2D);
      prev_del = DMX(r_0, j_1) + TSC(j_1, D2D);
      /* best-to-delete */
      prev_best = calc_Max( prev_mat, prev_del );
      DMX(r_0, j_0) = prev_best;

      /* UPDATE E STATE */
      sc1 = XMX(SP_E, x_0);
      sc2 = MMX(r_0, j_0);
      sc3 = DMX(r_0, j_0);
      XMX(SP_E, x_0) = calc_Max( sc1, 
                           calc_Max( sc2, sc3 ) );

      /* SPECIAL STATES */
      /* J state */
      sc_1 = XMX(SP_J, x_1) + XSC(SP_J, SP_LOOP);     /* J->J */
      sc_2 = XMX(SP_E, x_0) + XSC(SP_E, SP_LOOP);     /* E->J is E's "loop" */
      XMX(SP_J, x_0) = calc_Max( sc_1, sc_2 );

      /* C state */
      sc_1 = XMX(SP_C, x_1) + XSC(SP_C, SP_LOOP);
      sc_2 = XMX(SP_E, x_0) + XSC(SP_E, SP_MOVE);
      XMX(SP_C, x_0) = calc_Max( sc_1, sc_2 );

      /* N state */
      XMX(SP_N, x_0) = XMX(SP_N, x_1) + XSC(SP_N, SP_LOOP);

      /* B state */
      sc_1 = XMX(SP_N, x_0) + XSC(SP_N, SP_MOVE);     /* N->B is N's move */
      sc_2 = XMX(SP_J, x_0) + XSC(SP_J, SP_MOVE);     /* J->B is J's move */
      XMX(SP_B, x_0) = calc_Max( sc_1, sc_2 );

      /* embed linear row into quadratic test matrix */
      #if DEBUG
      {
         MX_2D(cloud_MX, x_0, j_0) = 1.0;
         MX_3D(test_MX, MAT_ST, x_0, j) = MMX(r_0, j_0);
         MX_3D(test_MX, INS_ST, x_0, j) = IMX(r_0, j_0);
         MX_3D(test_MX, DEL_ST, x_0, j) = DMX(r_0, j_0);
      }
      #endif
   }

   /* T state (stores final state score) */
   sc_best = XMX(SP_C, Q) + XSC(SP_C, SP_MOVE);
   *sc_final = sc_best;

   return STATUS_SUCCESS;
}