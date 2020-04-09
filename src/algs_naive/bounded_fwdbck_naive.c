/*******************************************************************************
 *  @file bounded.c
 *  @brief Bounded Forward-Backward Algorithm for Cloud Search (NAIVE).
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
#include "objects/edgebound.h"
#include "objects/hmm_profile.h"
#include "objects/sequence.h"
#include "objects/matrix/matrix_2d.h"
#include "objects/matrix/matrix_3d.h"

/* local imports */ 
#include "utilities/utility.h"
#include "testing.h"

/* header */
#include "bounded_fwdbck_naive.h"

/*  
 *  FUNCTION: bound_Forward_Naive()
 *  SYNOPSIS: Perform Forward part of Forward-Backward Algorithm.
 *
 *  ARGS:      <query>        query sequence, 
 *             <target>       HMM model,
 *             <Q>            query length, 
 *             <T>            target length,
 *             <st_MX>        Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>        Special State (J,N,B,C,E) Matrix,
 *             <st_MX_cloud>  Boolean Matrix of Cloud Bounds,
 *             <sc_final>     (OUTPUT) Final Score
 *
 *  RETURN:    Returns the Forward Score.
 */
float bound_Forward_Naive(const SEQUENCE*    query, 
                          const HMM_PROFILE* target, 
                          const int          Q, 
                          const int          T, 
                          MATRIX_3D*         st_MX, 
                          MATRIX_2D*         sp_MX,
                          MATRIX_3D*         st_MX_cloud, 
                          float*             sc_final)
{
   init_Logsum();

   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */
   int    i,j,k = 0;             /* row, column indices */
   char   *seq = query->seq;     /* alias for getting seq */
   int    x_0, r_0, x_1, r_1;    /* */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_sum;
   float  sc, sc_1, sc_2;
   float  sc_max, sc_best;
   float  sc1, sc2, sc3, sc4;
   float  sc_cloud;
   COORDS tr_end;       /* ending match state of optimal alignment (for traceback) */

   /* local or global (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* DEBUG OUTPUT */
   FILE *tfp = fopen("test_output/fwd_bounded.naive.txt", "w+");

   /* --------------------------------------------------------------------------------- */

   /* initialize special states (?) */
   XMX_M(SP_N,0) = 0;                                         /* S->N, p=1             */
   XMX_M(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX_M(SP_E,0) = XMX_M(SP_C,0) = XMX_M(SP_J,0) = -INF;          /* need seq to get here (?)  */

   /* initialize 0 row (top-edge) */
   for (j = 0; j < T; j++)
      { MMX_M(0,j) = IMX_M(0,j) = DMX_M(0,j) = -INF; }       /* need seq to get here (?)  */

   /* FOR every position in QUERY seq */
   for (i = 1; i <= Q; i++)
   {  
      fprintf(tfp, "(i=%d): ", i);
      /* Get next sequence character */
      a = seq[i-1];
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX_M(i,0) = IMX_M(i,0) = DMX_M(i,0) = -INF;
      XMX_M(SP_E,i) = -INF;

      /* MAIN RECURSION */
      /* FOR every position in TARGET profile */
      for (j = 1; j < T; j++)
      {
         /* if matrix cell is in the cloud, then compute cell */
         sc_cloud = ST_MX_M( st_MX_cloud, I_ST, i, j);
         if ( sc_cloud > 0 ) 
         {
            fprintf(tfp, "%d...", j);
            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            sc1 = prev_mat = MMX_M(i-1,j-1)  + TSC(j-1,M2M);
            sc2 = prev_ins = IMX_M(i-1,j-1)  + TSC(j-1,I2M);
            sc3 = prev_del = DMX_M(i-1,j-1)  + TSC(j-1,D2M);
            sc4 = prev_beg = XMX_M(SP_B,i-1) + TSC(j-1,B2M); /* from begin match state (new alignment) */
            /* best-to-match */
            prev_sum = calc_Logsum( 
                           calc_Logsum( prev_mat, prev_ins ),
                           calc_Logsum( prev_beg, prev_del ) );
            MMX_M(i,j) = prev_sum + MSC(j,A);
            // fprintf(tfp, "MMX_M(%d,%d): m=%f, i=%f, d=%f, b=%f\n", i, j, sc1, sc2, sc3, sc4 );

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX_M(i-1,j) + TSC(j,M2I);
            prev_ins = IMX_M(i-1,j) + TSC(j,I2I);
            /* best-to-insert */
            prev_sum = calc_Logsum( prev_mat, prev_ins );
            IMX_M(i,j) = prev_sum + ISC(j,A);

            /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX_M(i,j-1) + TSC(j-1,M2D);
            prev_del = DMX_M(i,j-1) + TSC(j-1,D2D);
            /* best-to-delete */
            prev_sum = calc_Logsum( prev_mat, prev_del );
            DMX_M(i,j) = prev_sum;

            /* UPDATE E STATE */
            sc1 = XMX_M(SP_E,i);
            sc2 = prev_mat = MMX_M(i,j) + sc_E;
            sc4 = prev_del = DMX_M(i,j) + sc_E;

            XMX_M(SP_E,i) = calc_Logsum( 
                              calc_Logsum( prev_mat, prev_del ),
                              XMX_M(SP_E,i) );

            x_0 = r_0 = i;
            // fprintf(tfp, "PRE (%d,%d) -> E: %f, Esc: %f, MMX: %f, DMX: %f, MSC: %f \n", x_0, j, XMX_M(SP_E, x_0), sc_E, MMX_M(r_0, j), DMX_M(r_0, j), MSC(j,A) );
         }
      }

      sc_cloud = ST_MX_M( st_MX_cloud, M_ST, i, T);
      if ( sc_cloud > 0 )
      {
         /* UNROLLED FINAL LOOP ITERATION */
         j = T; 
         fprintf(tfp, "%d!...", j);

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         sc1 = prev_mat = MMX_M(i-1,j-1)  + TSC(j-1,M2M);
         sc2 = prev_ins = IMX_M(i-1,j-1)  + TSC(j-1,I2M);
         sc3 = prev_del = DMX_M(i-1,j-1)  + TSC(j-1,D2M);
         sc4 = prev_beg = XMX_M(SP_B,i-1) + TSC(j-1,B2M);    /* from begin match state (new alignment) */
         /* sum-to-match */
         prev_sum = calc_Logsum( 
                           calc_Logsum( prev_mat, prev_ins ),
                           calc_Logsum( prev_del, prev_beg )
                        );
         MMX_M(i,j) = prev_sum + MSC(j,A);

         /* FIND SUM OF PATHS TO INSERT STATE */
         IMX_M(i,j) = -INF;

         /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) (unrolled) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX_M(i,j-1) + TSC(j-1,M2D);
         prev_del = DMX_M(i,j-1) + TSC(j-1,D2D);
         /* sum-to-delete */
         prev_sum = calc_Logsum( prev_mat, prev_del );
         DMX_M(i,j) = prev_sum;

         /* UPDATE E STATE (unrolled) */
         sc1 = XMX_M(SP_E,i);
         sc2 = MMX_M(i,j);
         sc4 = DMX_M(i,j);
         /* best-to-begin */
         XMX_M(SP_E,i) = calc_Logsum( 
                           calc_Logsum( DMX_M(i,j), MMX_M(i,j) ),
                           XMX_M(SP_E,i) );
      }
      fprintf(tfp, "...done\n");
      
      /* SPECIAL STATES */
      /* J state */
      sc_1 = XMX_M(SP_J,i-1) + XSC(SP_J,SP_LOOP);       /* J->J */
      sc_2 = XMX_M(SP_E,i)   + XSC(SP_E,SP_LOOP);       /* E->J is E's "loop" */
      XMX_M(SP_J,i) = calc_Logsum( sc_1, sc_2 );         

      /* C state */
      sc_1 = XMX_M(SP_C,i-1) + XSC(SP_C,SP_LOOP);
      sc_2 = XMX_M(SP_E,i)   + XSC(SP_E,SP_MOVE);
      XMX_M(SP_C,i) = calc_Logsum( sc_1, sc_2 );

      /* N state */
      XMX_M(SP_N,i) = XMX_M(SP_N,i-1) + XSC(SP_N,SP_LOOP);

      /* B state */
      sc_1 = XMX_M(SP_N,i) + XSC(SP_N,SP_MOVE);         /* N->B is N's move */
      sc_2 = XMX_M(SP_J,i) + XSC(SP_J,SP_MOVE);         /* J->B is J's move */
      XMX_M(SP_B,i) = calc_Logsum( sc_1, sc_2 );      
   }

   /* T state */
   sc_best = XMX_M(SP_C,Q) + XSC(SP_C,SP_MOVE);
   sc_final[0] = sc_best;

   fclose(tfp);
   return sc_best;
}

/* 
 *  FUNCTION: bound_Backward_Naive()
 *  SYNOPSIS: Perform Backward part of Forward-Backward Algorithm.
 *
 *  ARGS:      <query>        query sequence, 
 *             <target>       HMM model,
 *             <Q>            query length, 
 *             <T>            target length,
 *             <st_MX>        Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>        Special State (J,N,B,C,E) Matrix,
 *             <st_MX_cloud>  Boolean Matrix of Cloud Bounds,
 *             <sc_final>     (OUTPUT) Final Score
 *
 * RETURN:     Returns the Forward Score.
*/
float bound_Backward_Naive(const SEQUENCE*    query, 
                          const HMM_PROFILE* target, 
                          const int          Q, 
                          const int          T, 
                          MATRIX_3D*         st_MX, 
                          MATRIX_2D*         sp_MX,
                          MATRIX_3D*         st_MX_cloud, 
                          float*             sc_final)
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
   XMX_M(SP_J,Q) = XMX_M(SP_B,Q) = XMX_M(SP_N,Q) = -INF;
   XMX_M(SP_C,Q) = XSC(SP_C,SP_MOVE);
   XMX_M(SP_E,Q) = XMX_M(SP_C,Q) + XSC(SP_E,SP_MOVE);

   MMX_M(Q,T) = DMX_M(Q,T) = XMX_M(SP_E,Q);
   IMX_M(Q,T) = -INF;

   for (j = T-1; j >= 1; j--)
   {
      sc1 = XMX_M(SP_E,Q) + sc_E;
      sc2 = DMX_M(Q,j+1)  + TSC(j,M2D);
      MMX_M(Q,j) = calc_Logsum( XMX_M(SP_E,Q) + sc_E, 
                              DMX_M(Q,j+1)  + TSC(j,M2D) );

      sc1 = XMX_M(SP_E,Q) + sc_E;
      sc2 = DMX_M(Q,j+1)  + TSC(j,D2D);
      DMX_M(Q,j) = calc_Logsum( XMX_M(SP_E,Q) + sc_E,
                              DMX_M(Q,j+1)  + TSC(j,D2D) );

      IMX_M(Q,j) = -INF;
   }

   /* MAIN RECURSION */
   /* FOR every position in QUERY seq */
   for (i = Q-1; i >= 1; i--)
   {
      fprintf(tfp, "(i=%d): ", i);
      /* Get next sequence character */
      a = seq[i];
      A = AA_REV[a];

      /* UPDATE SPECIAL STATES at the start of EACH ROW */

      /* B STATE (COMPLETE) */
      XMX_M(SP_B,i) = -INF;
      for (j = 1; j <= T; j++)
      {
         XMX_M(SP_B,i) = calc_Logsum( XMX_M(SP_B,i),
                                    MMX_M(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
      }

      /* B STATE (SPARSE) */
      // XMX_M(SP_B,i) = -INF;
      // for (j = 1; j <= T; j++)
      // {
      //    sc_cloud = ST_MX_M( st_MX_cloud, I_ST, i, j);
      //    if ( sc_cloud > 0 ) {
      //       XMX_M(SP_B,i) = calc_Logsum( XMX_M(SP_B,i),
      //                                  MMX_M(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
      //    }
      // }

      XMX_M(SP_J,i) = calc_Logsum( XMX_M(SP_J,i+1) + XSC(SP_J,SP_LOOP),
                                 XMX_M(SP_B,  i) + XSC(SP_J,SP_MOVE) );

      XMX_M(SP_C,i) = XMX_M(SP_C,i+1) + XSC(SP_C,SP_LOOP);

      XMX_M(SP_E,i) = calc_Logsum( XMX_M(SP_J,i) + XSC(SP_E,SP_LOOP),
                                 XMX_M(SP_C,i) + XSC(SP_E,SP_MOVE) );

      XMX_M(SP_N,i) = calc_Logsum( XMX_M(SP_N,i+1) + XSC(SP_N,SP_LOOP),
                                 XMX_M(SP_B,  i) + XSC(SP_N,SP_MOVE) );

      MMX_M(i,T) = DMX_M(i,T) = XMX_M(SP_E,i);
      IMX_M(i,T) = -INF;

      /* FOR every position in TARGET profile */
      for (j = T-1; j >= 1; j--)
      {
         sc_cloud = ST_MX_M( st_MX_cloud, I_ST, i, j);
         if ( sc_cloud > 0 ) 
         {

            fprintf(tfp, "%d...", j);
            sc_M = MSC(j+1,A);
            sc_I = ISC(j,A);

            /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
            sc1 = prev_mat = MMX_M(i+1,j+1) + TSC(j,M2M) + sc_M;
            sc2 = prev_ins = IMX_M(i+1,j)   + TSC(j,M2I) + sc_I;
            sc3 = prev_del = DMX_M(i,j+1)   + TSC(j,M2D);
            sc4 = prev_end = XMX_M(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
            /* best-to-match */
            prev_sum = calc_Logsum( 
                              calc_Logsum( prev_mat, prev_ins ),
                              calc_Logsum( prev_end, prev_del )
                        );
            MMX_M(i,j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            sc1 = prev_mat = MMX_M(i+1,j+1) + TSC(j,I2M) + sc_M;
            sc2 = prev_ins = IMX_M(i+1,j)   + TSC(j,I2I) + sc_I;
            /* best-to-insert */
            prev_sum = calc_Logsum( prev_mat, prev_ins );
            IMX_M(i,j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            sc1 = prev_mat = MMX_M(i+1,j+1) + TSC(j,D2M) + sc_M;
            sc2 = prev_del = DMX_M(i,j+1)   + TSC(j,D2D);
            sc3 = prev_end = XMX_M(SP_E,i)  + sc_E;
            /* best-to-delete */
            prev_sum = calc_Logsum( 
                              prev_mat, 
                              calc_Logsum( prev_del, prev_end ) );
            DMX_M(i,j) = prev_sum;
         }
      }
      fprintf(tfp, "...done\n");
      fprintf(tfp, "N(0): %f\n", XMX_M(SP_N, 0));
   }

   /* FINAL i = 0 row */
   /* At i=0, only N,B states are reachable. */
   a = seq[1];
   A = AA_REV[a];

   /* t_BM index is 0 because it's stored off-by-one. */
   j = 1;
   XMX_M(SP_B,0) = MMX_M(1,j) + TSC(j-1,B2M) + MSC(j,A);
   for (j = 2; j >= T; j++) {
      XMX_M(SP_B,0) = calc_Logsum( XMX_M(SP_B,0),
                                 MMX_M(1,j) + TSC(j-1,B2M) + MSC(j,A) );
   }

   /* B STATE (SPARSE) */
   // XMX_M(SP_B,i) = -INF;
   // for (j = 1; j <= T; j++)
   // {
   //    sc_cloud = ST_MX_M( st_MX_cloud, I_ST, i, j);
   //    if ( sc_cloud > 0 ) {
   //       XMX_M(SP_B,i) = calc_Logsum( XMX_M(SP_B,i),
   //                                  MMX_M(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
   //    }
   // }

   XMX_M(SP_J,i) = -INF;
   XMX_M(SP_C,i) = -INF;
   XMX_M(SP_E,i) = -INF;

   XMX_M(SP_N,i) = calc_Logsum( XMX_M(SP_N,1) + XSC(SP_N,SP_LOOP),
                              XMX_M(SP_B,0) + XSC(SP_N,SP_MOVE) );

   for (j = T; j >= 1; j--) {
      MMX_M(i,j) = IMX_M(i,j) = DMX_M(i,j) = -INF;
   }

   fprintf(tfp, "a=%d, A=%d, MSC=%f\n", a, A, MSC(1,A) );

   sc_best = XMX_M(SP_N,0);
   sc_final[0] = sc_best;
   fclose(tfp);
   return sc_best;
}