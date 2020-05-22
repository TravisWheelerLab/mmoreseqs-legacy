/*******************************************************************************
 *  @file bound_fwdbck_naive.c
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

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "bound_fwdbck_naive.h"

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
                          MATRIX_2D*         st_MX_cloud, 
                          float*             sc_final)
{
   /* extract data from objects */
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   int      i, j, k;                         /* row, column indices */
   char*    seq;                             /* alias for getting seq */

   int      x, y1, y2;                       /* row, leftcol and rightcol bounds in row */
   int      x_0, row_cur, x_1, row_prv;      /* real index of current and previous rows */
   int      r_0, r_1;                        /* row offset -> r_0: row_cur % 2, r_1: row_prv % 2 */
   int      r_0b, r_0e, r_1b, r_1e;          /* begin and end indices for row in edgebound list */
   int      d_0, d_1, d_2;                   /* d (mod 3) for assigning diag array ptrs */
   bool     y2_re;                           /* checks if edge touches rightbound */

   float    prev_mat, prev_del, prev_ins;    /* temp placeholder sums */
   float    prev_beg, prev_end, prev_sum;    /* temp placeholder sums */
   float    sc, sc1, sc2, sc3, sc4;          /* temp placeholder sums (testing) */
   float    sc_M, sc_I, sc_D;                /* match, insert, delete scores */

   float    sc_cloud;

   /* local or global (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* --------------------------------------------------------------------------------- */

   /* clear leftover data */
   DP_MATRIX_Fill( Q, T, st_MX, sp_MX, -INF );

   /* query sequence */
   seq = query->seq;

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
         sc_cloud = *MATRIX_2D_Get( st_MX_cloud, i, j );
         if ( sc_cloud > 0 ) 
         {
            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            prev_mat = MMX(i-1,j-1)  + TSC(j-1,M2M);
            prev_ins = IMX(i-1,j-1)  + TSC(j-1,I2M);
            prev_del = DMX(i-1,j-1)  + TSC(j-1,D2M);
            prev_beg = XMX(SP_B,i-1) + TSC(j-1,B2M); /* from begin match state (new alignment) */
            /* best-to-match */
            prev_sum = logsum( 
                           logsum( prev_mat, prev_ins ),
                           logsum( prev_beg, prev_del ) );
            MMX(i,j) = prev_sum + MSC(j,A);

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX(i-1,j) + TSC(j,M2I);
            prev_ins = IMX(i-1,j) + TSC(j,I2I);
            /* best-to-insert */
            prev_sum = logsum( prev_mat, prev_ins );
            IMX(i,j) = prev_sum + ISC(j,A);

            /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX(i,j-1) + TSC(j-1,M2D);
            prev_del = DMX(i,j-1) + TSC(j-1,D2D);
            /* best-to-delete */
            prev_sum = logsum( prev_mat, prev_del );
            DMX(i,j) = prev_sum;

            /* UPDATE E STATE */
            prev_mat = MMX(i,j) + sc_E;
            prev_del = DMX(i,j) + sc_E;
            XMX(SP_E,i) = logsum( 
                              logsum( prev_mat, prev_del ),
                              XMX(SP_E,i) );

            x_0 = r_0 = i;
         }
      }

      /* if matrix cell is in the cloud, then compute cell */
      sc_cloud = *MATRIX_2D_Get( st_MX_cloud, i, T );
      if ( sc_cloud > 0 )
      {
         /* UNROLLED FINAL LOOP ITERATION */
         j = T; 

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         prev_mat = MMX(i-1,j-1)  + TSC(j-1,M2M);
         prev_ins = IMX(i-1,j-1)  + TSC(j-1,I2M);
         prev_del = DMX(i-1,j-1)  + TSC(j-1,D2M);
         prev_beg = XMX(SP_B,i-1) + TSC(j-1,B2M);    /* from begin match state (new alignment) */
         /* sum-to-match */
         prev_sum = logsum( 
                           logsum( prev_mat, prev_ins ),
                           logsum( prev_del, prev_beg )
                        );
         MMX(i,j) = prev_sum + MSC(j,A);

         /* FIND SUM OF PATHS TO INSERT STATE */
         IMX(i,j) = -INF;

         /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) (unrolled) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX(i,j-1) + TSC(j-1,M2D);
         prev_del = DMX(i,j-1) + TSC(j-1,D2D);
         /* sum-to-delete */
         prev_sum = logsum( prev_mat, prev_del );
         DMX(i,j) = prev_sum;

         /* UPDATE E STATE (unrolled) */
         /* best-to-begin */
         XMX(SP_E,i) = logsum( 
                           logsum( DMX(i,j), MMX(i,j) ),
                           XMX(SP_E,i) );
      }
      
      /* SPECIAL STATES */
      /* J state */
      sc1 = XMX(SP_J,i-1) + XSC(SP_J,SP_LOOP);       /* J->J */
      sc2 = XMX(SP_E,i)   + XSC(SP_E,SP_LOOP);       /* E->J is E's "loop" */
      XMX(SP_J,i) = logsum( sc1, sc2 );         

      /* C state */
      sc1 = XMX(SP_C,i-1) + XSC(SP_C,SP_LOOP);
      sc2 = XMX(SP_E,i)   + XSC(SP_E,SP_MOVE);
      XMX(SP_C,i) = logsum( sc1, sc2 );

      /* N state */
      XMX(SP_N,i) = XMX(SP_N,i-1) + XSC(SP_N,SP_LOOP);

      /* B state */
      sc1 = XMX(SP_N,i) + XSC(SP_N,SP_MOVE);         /* N->B is N's move */
      sc2 = XMX(SP_J,i) + XSC(SP_J,SP_MOVE);         /* J->B is J's move */
      XMX(SP_B,i) = logsum( sc1, sc2 );      
   }

   /* T state */
   *sc_final = XMX(SP_C,Q) + XSC(SP_C,SP_MOVE);
   return *sc_final;
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
float bound_Backward_Naive(const SEQUENCE*      query, 
                          const HMM_PROFILE*    target, 
                          const int             Q, 
                          const int             T, 
                          MATRIX_3D*            st_MX, 
                          MATRIX_2D*            sp_MX,
                          MATRIX_2D*            st_MX_cloud, 
                          float*                sc_final)
{
   /* extract data from objects */
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   int      i, j, k;                         /* row, column indices */
   char*    seq;                             /* alias for getting seq */

   int      x, y1, y2;                       /* row, leftcol and rightcol bounds in row */
   int      x_0, row_cur, x_1, row_prv;      /* real index of current and previous rows */
   int      r_0, r_1;                        /* row offset -> r_0: row_cur % 2, r_1: row_prv % 2 */
   int      r_0b, r_0e, r_1b, r_1e;          /* begin and end indices for row in edgebound list */
   int      d_0, d_1, d_2;                   /* d (mod 3) for assigning diag array ptrs */
   bool     y2_re;                           /* checks if edge touches rightbound */

   float    prev_mat, prev_del, prev_ins;    /* temp placeholder sums */
   float    prev_beg, prev_end, prev_sum;    /* temp placeholder sums */
   float    sc, sc1, sc2, sc3, sc4;          /* temp placeholder sums (testing) */
   float    sc_M, sc_I, sc_D;                /* match, insert, delete scores */

   float    sc_cloud;

   /* local or global (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* --------------------------------------------------------------------------------- */

   /* clear leftover data */
   DP_MATRIX_Fill( Q, T, st_MX, sp_MX, -INF );

   /* query sequence */
   seq = query->seq;

   /* Initialize the Q row. */
   XMX(SP_J,Q) = XMX(SP_B,Q) = XMX(SP_N,Q) = -INF;
   XMX(SP_C,Q) = XSC(SP_C,SP_MOVE);
   XMX(SP_E,Q) = XMX(SP_C,Q) + XSC(SP_E,SP_MOVE);

   /* if matrix cell is in the cloud, then compute cell */
   sc_cloud = *MATRIX_2D_Get( st_MX_cloud, Q, T );
   if ( sc_cloud > 0 )
   {
      MMX(Q,T) = DMX(Q,T) = XMX(SP_E,Q);
      IMX(Q,T) = -INF;
   }

   for (j = T-1; j >= 1; j--)
   {
      /* if matrix cell is in the cloud, then compute cell */
      sc_cloud = *MATRIX_2D_Get( st_MX_cloud, Q, j );
      if ( sc_cloud > 0 )
      {
         MMX(Q,j) = logsum( XMX(SP_E,Q) + sc_E, 
                           DMX(Q,j+1)  + TSC(j,M2D) );
         DMX(Q,j) = logsum( XMX(SP_E,Q) + sc_E,
                              DMX(Q,j+1)  + TSC(j,D2D) );
         IMX(Q,j) = -INF;
      }
   }

   /* MAIN RECURSION */
   /* FOR every position in QUERY seq */
   for (i = Q-1; i >= 1; i--)
   {
      /* Get next sequence character */
      a = seq[i];
      A = AA_REV[a];

      /* UPDATE SPECIAL STATES at the start of EACH ROW */

      /* B STATE (COMPLETE) */
      XMX(SP_B,i) = -INF;
      for (j = 1; j <= T; j++)
      {
         sc_cloud = *MATRIX_2D_Get( st_MX_cloud, i+1, j );
         if ( sc_cloud > 0 )
         {
            XMX(SP_B,i) = logsum( XMX(SP_B,i),
                                    MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
         }
      }

      XMX(SP_J,i) = logsum( XMX(SP_J,i+1) + XSC(SP_J,SP_LOOP),
                                 XMX(SP_B,  i) + XSC(SP_J,SP_MOVE) );

      XMX(SP_C,i) = XMX(SP_C,i+1) + XSC(SP_C,SP_LOOP);

      XMX(SP_E,i) = logsum( XMX(SP_J,i) + XSC(SP_E,SP_LOOP),
                                 XMX(SP_C,i) + XSC(SP_E,SP_MOVE) );

      XMX(SP_N,i) = logsum( XMX(SP_N,i+1) + XSC(SP_N,SP_LOOP),
                                 XMX(SP_B,  i) + XSC(SP_N,SP_MOVE) );

      /* if matrix cell is in the cloud, then compute cell */
      sc_cloud = *MATRIX_2D_Get( st_MX_cloud, i, T );
      if ( sc_cloud > 0 )
      {
         MMX(i,T) = DMX(i,T) = XMX(SP_E,i);
         IMX(i,T) = -INF;
      }

      /* FOR every position in TARGET profile */
      for (j = T-1; j >= 1; j--)
      {
         /* if matrix cell is in the cloud, then compute cell */
         sc_cloud = *MATRIX_2D_Get( st_MX_cloud, i, j );
         if ( sc_cloud > 0 ) 
         {
            sc_M = MSC(j+1,A);
            sc_I = ISC(j,A);

            /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
            sc1 = prev_mat = MMX(i+1,j+1) + TSC(j,M2M) + sc_M;
            sc2 = prev_ins = IMX(i+1,j)   + TSC(j,M2I) + sc_I;
            sc3 = prev_del = DMX(i,j+1)   + TSC(j,M2D);
            sc4 = prev_end = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
            /* best-to-match */
            prev_sum = logsum( 
                              logsum( prev_mat, prev_ins ),
                              logsum( prev_end, prev_del )
                        );
            MMX(i,j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            sc1 = prev_mat = MMX(i+1,j+1) + TSC(j,I2M) + sc_M;
            sc2 = prev_ins = IMX(i+1,j)   + TSC(j,I2I) + sc_I;
            /* best-to-insert */
            prev_sum = logsum( prev_mat, prev_ins );
            IMX(i,j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            sc1 = prev_mat = MMX(i+1,j+1) + TSC(j,D2M) + sc_M;
            sc2 = prev_del = DMX(i,j+1)   + TSC(j,D2D);
            sc3 = prev_end = XMX(SP_E,i)  + sc_E;
            /* best-to-delete */
            prev_sum = logsum( 
                              prev_mat, 
                              logsum( prev_del, prev_end ) );
            DMX(i,j) = prev_sum;
         }
      }
   }

   /* FINAL i = 0 row */
   /* At i=0, only N,B states are reachable. */
   a = seq[1];
   A = AA_REV[a];

   /* t_BM index is 0 because it's stored off-by-one. */
   j = 1;
   XMX(SP_B,0) = MMX(1,j) + TSC(j-1,B2M) + MSC(j,A);
   for (j = 2; j >= T; j++) {
      XMX(SP_B,0) = logsum( XMX(SP_B,0),
                                 MMX(1,j) + TSC(j-1,B2M) + MSC(j,A) );
   }

   /* B STATE (SPARSE) */
   // XMX(SP_B,i) = -INF;
   // for (j = 1; j <= T; j++)
   // {
   //    sc_cloud = ST_MX( st_MX_cloud, I_ST, i, j);
   //    if ( sc_cloud > 0 ) {
   //       XMX(SP_B,i) = logsum( XMX(SP_B,i),
   //                                  MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
   //    }
   // }

   XMX(SP_J,i) = -INF;
   XMX(SP_C,i) = -INF;
   XMX(SP_E,i) = -INF;

   XMX(SP_N,i) = logsum( XMX(SP_N,1) + XSC(SP_N,SP_LOOP),
                              XMX(SP_B,0) + XSC(SP_N,SP_MOVE) );

   for (j = T; j >= 1; j--) {
      MMX(i,j) = IMX(i,j) = DMX(i,j) = -INF;
   }

   *sc_final = XMX(SP_N,0);
   return *sc_final;
}
