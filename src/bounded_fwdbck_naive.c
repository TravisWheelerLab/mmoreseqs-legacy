/*******************************************************************************
 *  @file bounded.c
 *  @brief Bounded Forward-Backward Algorithm for Cloud Search (NAIVE).
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

// imports
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

// local imports (after struct declarations)
#include "structs.h"
#include "misc.h"
#include "hmm_parser.h"
#include "viterbi.h"
#include "forward_backward.h"
#include "cloud_search_quad.h"
#include "bounded_fwdbck_naive.h"

/*  
 *  FUNCTION: forward_Bounded_Naive_Run()
 *  SYNOPSIS: Perform Forward part of Forward-Backward Algorithm.
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <res>       Results Data
 *
 *  RETURN: 
 */
float forward_Bounded_Naive_Run (const SEQ* query, 
                               const HMM_PROFILE* target, 
                               int Q, int T, 
                               float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                               float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                               float st_MX_cloud[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                               float *sc_final)
{
   init_Logsum();

   char   a;           /* store current character in sequence */
   int    A;           /* store int value of character */
   int    i,j,k = 0;   /* row, column indices */
   char   *seq = query->seq; /* alias for getting seq */
   int    x_0, r_0, x_1, r_1;

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_sum;
   float  sc, sc_1, sc_2;
   float  sc_max, sc_best;
   float  sc1, sc2, sc3, sc4;
   float  sc_cloud;
   COORDS tr_end; /* ending match state of optimal alignment (for traceback) */

   /* local or global (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* DEBUG OUTPUT */
   FILE *tfp = fopen("test_output/fwd_bounded.naive.txt", "w+");

   /* --------------------------------------------------------------------------------- */

   /* clear matrix */
   dp_matrix_Clear(Q, T, st_MX, sp_MX);

   /* initialize special states (?) */
   XMX(SP_N,0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E,0) = XMX(SP_C,0) = XMX(SP_J,0) = -INF;          /* need seq to get here (?)  */

   /* initialize 0 row (top-edge) */
   for (j = 0; j < T; j++)
      { MMX(0,j) = IMX(0,j) = DMX(0,j) = -INF; }       /* need seq to get here (?)  */

   /* FOR every position in QUERY seq */
   for (i = 1; i <= Q; i++)
   {  
      fprintf(tfp, "(i=%d): ", i);
      /* Get next sequence character */
      a = seq[i-1];
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX(i,0) = IMX(i,0) = DMX(i,0) = -INF;
      XMX(SP_E,i) = -INF;

      /* MAIN RECURSION */
      /* FOR every position in TARGET profile */
      for (j = 1; j < T; j++)
      {
         /* if matrix cell is in the cloud, then compute cell */
         sc_cloud = ST_MX( st_MX_cloud, I_ST, i, j);
         if ( sc_cloud > 0 ) 
         {
            fprintf(tfp, "%d...", j);
            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            sc1 = prev_mat = MMX(i-1,j-1)  + TSC(j-1,M2M);
            sc2 = prev_ins = IMX(i-1,j-1)  + TSC(j-1,I2M);
            sc3 = prev_del = DMX(i-1,j-1)  + TSC(j-1,D2M);
            sc4 = prev_beg = XMX(SP_B,i-1) + TSC(j-1,B2M); /* from begin match state (new alignment) */
            /* best-to-match */
            prev_sum = calc_Logsum( 
                           calc_Logsum( prev_mat, prev_ins ),
                           calc_Logsum( prev_beg, prev_del ) );
            MMX(i,j) = prev_sum + MSC(j,A);
            // fprintf(tfp, "MMX(%d,%d): m=%f, i=%f, d=%f, b=%f\n", i, j, sc1, sc2, sc3, sc4 );

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX(i-1,j) + TSC(j,M2I);
            prev_ins = IMX(i-1,j) + TSC(j,I2I);
            /* best-to-insert */
            prev_sum = calc_Logsum( prev_mat, prev_ins );
            IMX(i,j) = prev_sum + ISC(j,A);

            /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX(i,j-1) + TSC(j-1,M2D);
            prev_del = DMX(i,j-1) + TSC(j-1,D2D);
            /* best-to-delete */
            prev_sum = calc_Logsum( prev_mat, prev_del );
            DMX(i,j) = prev_sum;

            /* UPDATE E STATE */
            sc1 = XMX(SP_E,i);
            sc2 = prev_mat = MMX(i,j) + sc_E;
            sc4 = prev_del = DMX(i,j) + sc_E;

            XMX(SP_E,i) = calc_Logsum( 
                              calc_Logsum( prev_mat, prev_del ),
                              XMX(SP_E,i) );

            x_0 = r_0 = i;
            // fprintf(tfp, "PRE (%d,%d) -> E: %f, Esc: %f, MMX: %f, DMX: %f, MSC: %f \n", x_0, j, XMX(SP_E, x_0), sc_E, MMX(r_0, j), DMX(r_0, j), MSC(j,A) );
         }
      }

      sc_cloud = ST_MX( st_MX_cloud, M_ST, i, T);
      if ( sc_cloud > 0 )
      {
         /* UNROLLED FINAL LOOP ITERATION */
         j = T; 
         fprintf(tfp, "%d!...", j);

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         sc1 = prev_mat = MMX(i-1,j-1)  + TSC(j-1,M2M);
         sc2 = prev_ins = IMX(i-1,j-1)  + TSC(j-1,I2M);
         sc3 = prev_del = DMX(i-1,j-1)  + TSC(j-1,D2M);
         sc4 = prev_beg = XMX(SP_B,i-1) + TSC(j-1,B2M);    /* from begin match state (new alignment) */
         /* sum-to-match */
         prev_sum = calc_Logsum( 
                           calc_Logsum( prev_mat, prev_ins ),
                           calc_Logsum( prev_del, prev_beg )
                        );
         MMX(i,j) = prev_sum + MSC(j,A);

         /* FIND SUM OF PATHS TO INSERT STATE */
         IMX(i,j) = -INF;

         /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) (unrolled) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX(i,j-1) + TSC(j-1,M2D);
         prev_del = DMX(i,j-1) + TSC(j-1,D2D);
         /* sum-to-delete */
         prev_sum = calc_Logsum( prev_mat, prev_del );
         DMX(i,j) = prev_sum;

         /* UPDATE E STATE (unrolled) */
         sc1 = XMX(SP_E,i);
         sc2 = MMX(i,j);
         sc4 = DMX(i,j);
         /* best-to-begin */
         XMX(SP_E,i) = calc_Logsum( 
                           calc_Logsum( DMX(i,j), MMX(i,j) ),
                           XMX(SP_E,i) );
      }
      fprintf(tfp, "...done\n");
      
      /* SPECIAL STATES */
      /* J state */
      sc_1 = XMX(SP_J,i-1) + XSC(SP_J,SP_LOOP);       /* J->J */
      sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_LOOP);       /* E->J is E's "loop" */
      XMX(SP_J,i) = calc_Logsum( sc_1, sc_2 );         

      /* C state */
      sc_1 = XMX(SP_C,i-1) + XSC(SP_C,SP_LOOP);
      sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_MOVE);
      XMX(SP_C,i) = calc_Logsum( sc_1, sc_2 );

      /* N state */
      XMX(SP_N,i) = XMX(SP_N,i-1) + XSC(SP_N,SP_LOOP);

      /* B state */
      sc_1 = XMX(SP_N,i) + XSC(SP_N,SP_MOVE);         /* N->B is N's move */
      sc_2 = XMX(SP_J,i) + XSC(SP_J,SP_MOVE);         /* J->B is J's move */
      XMX(SP_B,i) = calc_Logsum( sc_1, sc_2 );      
   }

   /* T state */
   sc_best = XMX(SP_C,Q) + XSC(SP_C,SP_MOVE);
   sc_final[0] = sc_best;

   fclose(tfp);
   return sc_best;
}

/* FUNCTION: backward_Bounded_Naive_Run()
 * SYNOPSIS: Perform Backward part of Forward-Backward Algorithm.
 *
 * PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <res>       Results Data
 *
 * RETURN: 
*/
float backward_Bounded_Naive_Run (const SEQ* query, 
                                const HMM_PROFILE* target, 
                                int Q, int T, 
                                float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                                float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                                float st_MX_cloud[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                                float *sc_final)
{
   init_Logsum();

   char   a;           /* store current character in sequence */
   int    A;           /* store int value of character */
   int    i,j,k = 0;   /* row, column indices */
   char   *seq = query->seq; /* alias for getting seq */

   float  prev_mat, prev_del, prev_ins, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_M, sc_I;
   float  sc_max, sc_best;
   float  sc1, sc2, sc3, sc4;
   float  sc_cloud;

   /* local or global (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* DEBUG OUTPUT */
   FILE *tfp = fopen("test_output/bck_bounded.naive.txt", "w+");

   /* --------------------------------------------------------------------------------- */

   /* Initialize the Q row. */
   XMX(SP_J,Q) = XMX(SP_B,Q) = XMX(SP_N,Q) = -INF;
   XMX(SP_C,Q) = XSC(SP_C,SP_MOVE);
   XMX(SP_E,Q) = XMX(SP_C,Q) + XSC(SP_E,SP_MOVE);

   MMX(Q,T) = DMX(Q,T) = XMX(SP_E,Q);
   IMX(Q,T) = -INF;

   for (j = T-1; j >= 1; j--)
   {
      sc1 = XMX(SP_E,Q) + sc_E;
      sc2 = DMX(Q,j+1)  + TSC(j,M2D);
      MMX(Q,j) = calc_Logsum( XMX(SP_E,Q) + sc_E, 
                              DMX(Q,j+1)  + TSC(j,M2D) );

      sc1 = XMX(SP_E,Q) + sc_E;
      sc2 = DMX(Q,j+1)  + TSC(j,D2D);
      DMX(Q,j) = calc_Logsum( XMX(SP_E,Q) + sc_E,
                              DMX(Q,j+1)  + TSC(j,D2D) );

      IMX(Q,j) = -INF;
   }

   /* MAIN RECURSION */
   /* FOR every position in QUERY seq */
   for (i = Q-1; i >= 1; i--)
   {
      fprintf(tfp, "(i=%d): ", i);
      /* Get next sequence character */
      a = seq[i];
      A = AA_REV[a];

      /* SPECIAL STATES */
      XMX(SP_B,i) = MMX(i+1,1) + TSC(0,B2M) + MSC(1,A);

      /* B -> MATCH */
      for (j = 2; j <= T; j++)
      {
         XMX(SP_B,i) = calc_Logsum( XMX(SP_B,i),
                                    MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
      }

      XMX(SP_J,i) = calc_Logsum( XMX(SP_J,i+1) + XSC(SP_J,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_J,SP_MOVE) );

      XMX(SP_C,i) = XMX(SP_C,i+1) + XSC(SP_C,SP_LOOP);

      XMX(SP_E,i) = calc_Logsum( XMX(SP_J,i) + XSC(SP_E,SP_LOOP),
                                 XMX(SP_C,i) + XSC(SP_E,SP_MOVE) );

      XMX(SP_N,i) = calc_Logsum( XMX(SP_N,i+1) + XSC(SP_N,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_N,SP_MOVE) );

      MMX(i,T) = DMX(i,T) = XMX(SP_E,i);
      IMX(i,T) = -INF;

      /* FOR every position in TARGET profile */
      for (j = T-1; j >= 1; j--)
      {
         sc_cloud = ST_MX( st_MX_cloud, I_ST, i, j);
         if ( sc_cloud > 0 ) 
         {

            fprintf(tfp, "%d...", j);
            sc_M = MSC(j+1,A);
            sc_I = ISC(j,A);

            /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
            sc1 = prev_mat = MMX(i+1,j+1) + TSC(j,M2M) + sc_M;
            sc2 = prev_ins = IMX(i+1,j)   + TSC(j,M2I) + sc_I;
            sc3 = prev_del = DMX(i,j+1)   + TSC(j,M2D);
            sc4 = prev_end = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
            /* best-to-match */
            prev_sum = calc_Logsum( 
                              calc_Logsum( prev_mat, prev_ins ),
                              calc_Logsum( prev_end, prev_del )
                        );
            MMX(i,j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            sc1 = prev_mat = MMX(i+1,j+1) + TSC(j,I2M) + sc_M;
            sc2 = prev_ins = IMX(i+1,j)   + TSC(j,I2I) + sc_I;
            /* best-to-insert */
            prev_sum = calc_Logsum( prev_mat, prev_ins );
            IMX(i,j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            sc1 = prev_mat = MMX(i+1,j+1) + TSC(j,D2M) + sc_M;
            sc2 = prev_del = DMX(i,j+1)   + TSC(j,D2D);
            sc3 = prev_end = XMX(SP_E,i)  + sc_E;
            /* best-to-delete */
            prev_sum = calc_Logsum( 
                              prev_mat, 
                              calc_Logsum( prev_del, prev_end ) );
            DMX(i,j) = prev_sum;
         }
      }
      fprintf(tfp, "...done\n");
      fprintf(tfp, "N(0): %f\n", XMX(SP_N, 0));
   }

   /* FINAL i = 0 row */
   /* At i=0, only N,B states are reachable. */
   a = seq[1];
   A = AA_REV[a];

   /* t_BM index is 0 because it's stored off-by-one. */
   XMX(SP_B,0) = MMX(1,1) + TSC(0,B2M) + MSC(1,A);

   for (j = 2; j >= T; j++) {
      XMX(SP_B,0) = calc_Logsum( XMX(SP_B,0),
                                 MMX(1,j) + TSC(j-1,B2M) + MSC(j,A) );
   }

   XMX(SP_J,i) = -INF;
   XMX(SP_C,i) = -INF;
   XMX(SP_E,i) = -INF;

   XMX(SP_N,i) = calc_Logsum( XMX(SP_N,1) + XSC(SP_N,SP_LOOP),
                              XMX(SP_B,0) + XSC(SP_N,SP_MOVE) );

   for (j = T; j >= 1; j--) {
      MMX(i,j) = IMX(i,j) = DMX(i,j) = -INF;
   }

   fprintf(tfp, "a=%d, A=%d, MSC=%f\n", a, A, MSC(1,A) );

   sc_best = XMX(SP_N,0);
   sc_final[0] = sc_best;
   fclose(tfp);
   return sc_best;
}








