/*******************************************************************************
 *  FILE:      fwdback_linear.c
 *  SYNOPSIS:  The Forward-Backward Algorithm for Sequence Alignment Search.
 *             (Linear Space Alg)
 *
 *  AUTHOR:    Dave Rich
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
#include "algs_linear.h"

/* header */
#include "fwdback_quad.h"


/* ****************************************************************************************** * 
 *  FUNCTION: forward_Run()
 *  SYNOPSIS: Perform Forward part of Forward-Backward Algorithm. (LINEAR ALG)
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX3>    Normal State (Match, Insert, Delete) Matrix,
 *                         [ NORMAL_STATES * 3 * (Q+T+1) ]
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *                         [ SPECIAL_STATES * (Q+1) ]
 *             <res>       Results Data
 *
 *  RETURN: 
 * ****************************************************************************************** */
float forward_Linear(const SEQUENCE*   query, 
                     const HMM_PROFILE* target, 
                     const int          Q, 
                     const int          T, 
                     MATRIX_3D*         st_MX3,
                     MATRIX_2D*         sp_MX,
                     float*             sc_final)
{
   logsum_Init();

   char     a;          /* store current character in sequence */
   int      A;          /* store int value of character */
   int      i,j,k = 0;  /* row, column indices */
   char*    seq;        /* alias for getting seq */

   float    prev_mat, prev_del, prev_ins; 
   float    prev_beg, prev_end, prev_sum;
   float    sc, sc_1, sc_2;
   float    sc_max, sc_best;
   float    sc1, sc2, sc3, sc4;
   float    sc_M, sc_I, sc_D;
   int      x_0, x_1;
   int      r_0, r_1;                 /* d (mod 3) for assigning prev array ptrs */

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
      dbfp     = fopen( debugger->dbfp_path, "w+" );
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

   /* initialize special states (?) */
   XMX(SP_N,0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E,0) = XMX(SP_C,0) = XMX(SP_J,0) = -INF;          /* need seq to get here (?)  */

   /* initialize 0 row (top-edge) */
   r_0 = 0 % 2;
   for (j = 0; j < T; j++)
      MMX3(0,j) = IMX3(0,j) = DMX3(0,j) = -INF;      /* need seq to get here (?)  */

   /* FOR every position in QUERY seq */
   for (i = 1; i <= Q; i++)
   {  
      x_1 = x_0;
      x_0 = i;
      r_1 = r_0;
      r_0 = i % 2;

      /* Get next sequence character */
      a = seq[i-1];
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX3(r_0,0) = IMX3(r_0,0) = DMX3(r_0,0) = -INF;
      XMX(SP_E,i) = -INF;

      /* MAIN RECURSION */
      /* FOR every position in TARGET profile */
      for (j = 1; j < T; j++)
      {
         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         sc1 = prev_mat = MMX3(r_1,j-1)  + TSC(j-1,M2M);
         sc2 = prev_ins = IMX3(r_1,j-1)  + TSC(j-1,I2M);
         sc3 = prev_del = DMX3(r_1,j-1)  + TSC(j-1,D2M);
         sc4 = prev_beg = XMX(SP_B,i-1)  + TSC(j-1,B2M); /* from begin match state (new alignment) */

         /* best-to-match */
         prev_sum = logsum( 
                        logsum( prev_mat, prev_ins ),
                        logsum( prev_beg, prev_del ) );
         MMX3(r_0,j) = prev_sum + MSC(j,A);

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX3(r_0,j) + TSC(j,M2I);
         prev_ins = IMX3(r_0,j) + TSC(j,I2I);
         /* best-to-insert */
         prev_sum = logsum( prev_mat, prev_ins );
         IMX3(r_0,j) = prev_sum + ISC(j,A);

         /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX3(r_0,j-1) + TSC(j-1,M2D);
         prev_del = DMX3(r_0,j-1) + TSC(j-1,D2D);
         /* best-to-delete */
         prev_sum = logsum( prev_mat, prev_del );
         DMX3(r_0,j) = prev_sum;

         /* UPDATE E STATE */
         sc1 = XMX(SP_E,i);
         sc2 = prev_mat = MMX3(r_0,j) + sc_E;
         sc4 = prev_del = DMX3(r_0,j) + sc_E;
         XMX(SP_E,i) = logsum( 
                           logsum( prev_mat, prev_del ),
                           XMX(SP_E,i) );

         /* embed linear row into quadratic test matrix */
         #if DEBUG
         {
            MX_2D(cloud_MX, x_0, j) = 1.0;
            MX_3D(test_MX, MAT_ST, x_0, j) = MMX3(r_0, j);
            MX_3D(test_MX, INS_ST, x_0, j) = IMX3(r_0, j);
            MX_3D(test_MX, DEL_ST, x_0, j) = DMX3(r_0, j);
         }
         #endif
      }

      /* UNROLLED FINAL LOOP ITERATION */
      j = T; 

      /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
      /* best previous state transition (match takes the diag element of each prev state) */
      sc1 = prev_mat = MMX3(r_1,j-1)  + TSC(j-1,M2M);
      sc2 = prev_ins = IMX3(r_1,j-1)  + TSC(j-1,I2M);
      sc3 = prev_del = DMX3(r_1,j-1)  + TSC(j-1,D2M);
      sc4 = prev_beg = XMX(SP_B,i-1)  + TSC(j-1,B2M);    /* from begin match state (new alignment) */
      /* sum-to-match */
      prev_sum = logsum( 
                        logsum( prev_mat, prev_ins ),
                        logsum( prev_del, prev_beg ) );
      MMX3(r_0,j) = prev_sum + MSC(j,A);

      /* FIND SUM OF PATHS TO INSERT STATE (unrolled) */
      IMX3(r_0,j) = -INF;

      /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) (unrolled) */
      /* previous states (match takes the left element of each state) */
      prev_mat = MMX3(r_0,j-1) + TSC(j-1,M2D);
      prev_del = DMX3(r_0,j-1) + TSC(j-1,D2D);
      /* sum-to-delete */
      prev_sum = logsum( prev_mat, prev_del );
      DMX3(r_0,j) = prev_sum;

      /* UPDATE E STATE (unrolled) */
      sc1 = XMX(SP_E,i);
      sc2 = MMX3(r_0,j);
      sc4 = DMX3(r_0,j);
      /* best-to-begin */
      XMX(SP_E,i) = logsum( 
                        logsum( DMX3(r_0,j), MMX3(r_0,j) ),
                        XMX(SP_E,i) );

      /* SPECIAL STATES */
      /* J state */
      sc_1 = XMX(SP_J,i-1) + XSC(SP_J,SP_LOOP);       /* J->J */
      sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_LOOP);       /* E->J is E's "loop" */
      XMX(SP_J,i) = logsum( sc_1, sc_2 );         

      /* C state */
      sc_1 = XMX(SP_C,i-1) + XSC(SP_C,SP_LOOP);
      sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_MOVE);
      XMX(SP_C,i) = logsum( sc_1, sc_2 );

      /* N state */
      XMX(SP_N,i) = XMX(SP_N,i-1) + XSC(SP_N,SP_LOOP);

      /* B state */
      sc_1 = XMX(SP_N,i) + XSC(SP_N,SP_MOVE);         /* N->B is N's move */
      sc_2 = XMX(SP_J,i) + XSC(SP_J,SP_MOVE);         /* J->B is J's move */
      XMX(SP_B,i) = logsum( sc_1, sc_2 );     

      /* embed linear row into quadratic test matrix */
      #if DEBUG
      {
         MX_2D(cloud_MX, x_0, j) = 1.0;
         MX_3D(test_MX, MAT_ST, x_0, j) = MMX3(r_0, j);
         MX_3D(test_MX, INS_ST, x_0, j) = IMX3(r_0, j);
         MX_3D(test_MX, DEL_ST, x_0, j) = DMX3(r_0, j);
      }
      #endif
   }

   #if DEBUG
      fclose(dbfp);
   #endif

   /* T state */
   sc_best = XMX(SP_C,Q) + XSC(SP_C,SP_MOVE);
   *sc_final = sc_best; 
   return sc_best;
}

/* ****************************************************************************************** *
 * FUNCTION: backward_Linear()
 * SYNOPSIS: Perform Backward part of Forward-Backward Algorithm. (LINEAR ALG)
 *
 * PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX3>    Normal State (Match, Insert, Delete) Matrix,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <res>       Results Data
 *
 * RETURN: 
 * ****************************************************************************************** */
float backward_Linear(const SEQUENCE*    query, 
                     const HMM_PROFILE* target, 
                     const int          Q, 
                     const int          T, 
                     MATRIX_3D*         st_MX3, 
                     MATRIX_2D*         sp_MX,
                     float*             sc_final)
{
   logsum_Init();

   char     a;                        /* store current character in sequence */
   int      A;                        /* store int value of character */
   int      i,j,k = 0;                /* row, column indices */
   char*    seq;        /* alias for getting seq */

   float    prev_mat, prev_del, prev_ins; 
   float    prev_beg, prev_end, prev_sum;
   float    sc, sc_1, sc_2;
   float    sc_max, sc_best;
   float    sc1, sc2, sc3, sc4;
   float    sc_M, sc_I, sc_D;
   int      x_0, x_1;
   int      r_0, r_1;                 /* d (mod 3) for assigning prev array ptrs */

   /* local or global (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* --------------------------------------------------------------------------------- */

   /* query sequence */
   seq = query->seq;

   /* Initialize the Q row. */
   r_0 = Q % 2;

   XMX(SP_J,Q) = XMX(SP_B,Q) = XMX(SP_N,Q) = -INF;
   XMX(SP_C,Q) = XSC(SP_C,SP_MOVE);
   XMX(SP_E,Q) = XMX(SP_C,Q) + XSC(SP_E,SP_MOVE);

   MMX3(r_0,T) = DMX3(r_0,T) = XMX(SP_E,Q);
   IMX3(r_0,T) = -INF;

   for (j = T-1; j >= 1; j--)
   {
      sc1 = XMX(SP_E,Q) + sc_E;
      sc2 = DMX3(r_0,j+1)  + TSC(j,M2D);
      MMX3(r_0,j) = logsum( XMX(SP_E,Q) + sc_E, 
                              DMX3(r_0,j+1)  + TSC(j,M2D) );

      sc1 = XMX(SP_E,Q) + sc_E;
      sc2 = DMX3(r_0,j+1)  + TSC(j,D2D);
      DMX3(r_0,j) = logsum( XMX(SP_E,Q) + sc_E,
                              DMX3(r_0,j+1)  + TSC(j,D2D) );

      IMX3(r_0,j) = -INF;
   }

   /* MAIN RECURSION */
   /* FOR every position in QUERY seq */
   for (i = Q-1; i >= 1; i--)
   {
      x_1 = x_0;
      x_0 = i;
      r_1 = r_0;
      r_0 = i % 2;

      /* Get next sequence character */
      a = seq[i];
      A = AA_REV[a];

      /* SPECIAL STATES */
      j = 1; int x_0 = i; 
      XMX(SP_B,i) = MMX3(r_1,j) + TSC(j-1,B2M) + MSC(j,A);

      /* B -> MATCH */
      for (j = 2; j <= T; j++) {
         XMX(SP_B,i) = logsum( XMX(SP_B,i),
                                    MMX3(r_1,j) + TSC(j-1,B2M) + MSC(j,A) );
      }
      // printf("\n");

      XMX(SP_J,i) = logsum( XMX(SP_J,i+1) + XSC(SP_J,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_J,SP_MOVE) );

      XMX(SP_C,i) = XMX(SP_C,i+1) + XSC(SP_C,SP_LOOP);

      XMX(SP_E,i) = logsum( XMX(SP_J,i) + XSC(SP_E,SP_LOOP),
                                 XMX(SP_C,i) + XSC(SP_E,SP_MOVE) );

      XMX(SP_N,i) = logsum( XMX(SP_N,i+1) + XSC(SP_N,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_N,SP_MOVE) );

      MMX3(r_0,T) = DMX3(r_0,T) = XMX(SP_E,i);
      IMX3(r_0,T) = -INF;

      x_0 = i;

      /* FOR every position in TARGET profile */
      for (j = T-1; j >= 1; j--)
      {
         sc_M = MSC(j+1,A);
         sc_I = ISC(j,A);

         /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
         sc1 = prev_mat = MMX3(r_1,j+1) + TSC(j,M2M) + sc_M;
         sc2 = prev_ins = IMX3(r_1,j)   + TSC(j,M2I) + sc_I;
         sc3 = prev_del = DMX3(r_0,j+1) + TSC(j,M2D);
         sc4 = prev_end = XMX(SP_E,i)   + sc_E;     /* from end match state (new alignment) */
         /* best-to-match */
         prev_sum = logsum( 
                           logsum( prev_mat, prev_ins ),
                           logsum( prev_end, prev_del ) );
         MMX3(r_0,j) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
         sc1 = prev_mat = MMX3(r_1,j+1) + TSC(j,I2M) + sc_M;
         sc2 = prev_ins = IMX3(r_1,j)   + TSC(j,I2I) + sc_I;
         /* best-to-insert */
         prev_sum = logsum( prev_mat, prev_ins );
         IMX3(r_0,j) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
         sc1 = prev_mat = MMX3(r_1,j+1) + TSC(j,D2M) + sc_M;
         sc2 = prev_del = DMX3(r_0,j+1) + TSC(j,D2D);
         sc3 = prev_end = XMX(SP_E,i)  + sc_E;
         /* best-to-delete */
         prev_sum = logsum( 
                           prev_mat, 
                           logsum( prev_del, prev_end ) );
         DMX3(r_0,j) = prev_sum;
      }
   }

   /* FINAL i = 0 row */
   x_1 = x_0;
   x_0 = 0;
   r_1 = r_0;
   r_0 = 0;
   /* At i=0, only N,B states are reachable. */
   a = seq[1];
   A = AA_REV[a];

   /* t_BM index is 0 because it's stored off-by-one. */
   XMX(SP_B,0) = MMX3(r_1,1) + TSC(0,B2M) + MSC(1,A);

   for (j = 2; j >= T; j++) {
      XMX(SP_B,0) = logsum( XMX(SP_B,0),
                                 MMX3(r_1,j) + TSC(j-1,B2M) + MSC(j,A) );
   }

   XMX(SP_J,i) = -INF;
   XMX(SP_C,i) = -INF;
   XMX(SP_E,i) = -INF;

   XMX(SP_N,i) = logsum( XMX(SP_N,1) + XSC(SP_N,SP_LOOP),
                              XMX(SP_B,0) + XSC(SP_N,SP_MOVE) );

   for (j = T; j >= 1; j--) {
      MMX3(r_0,j) = IMX3(r_0,j) = DMX3(r_0,j) = -INF;
   }

   sc_best = XMX(SP_N,0);
   *sc_final = sc_best;
   return sc_best;
}


