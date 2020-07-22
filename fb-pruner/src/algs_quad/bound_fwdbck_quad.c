/*******************************************************************************
 *  FILE:      bounded_fwdbck_quad.c
 *  PURPOSE:   Bounded Forward-Backward Algorithm (QUADRATIC SPACE)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
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
#include "bound_fwdbck_quad.h"

/*  
 *  FUNCTION: forward_bounded_Run()
 *  SYNOPSIS: Perform Edge-Bounded Forward part of Cloud Search Algorithm.
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <edg>       Bounds Data (row-wise)
 *             <sc_final>  Final Score
 *
 *  RETURN: 
 */
float bound_Forward_Quad(  const SEQUENCE*      query, 
                           const HMM_PROFILE*   target,
                           const int            Q, 
                           const int            T, 
                           MATRIX_3D*           st_MX,
                           MATRIX_2D*           sp_MX, 
                           EDGEBOUNDS*          edg,
                           float*               sc_final )
{
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   int      i,j,k = 0;                       /* row, column indices */
   char     *seq = query->seq;               /* alias for getting seq */
   int      N = edg->N;                      /* length of edgebound list */

   int      x, y1, y2;                       /* row, left-edge and right-edge bounds in row */
   int      row_cur, row_prv;                /* real index of current and previous rows */
   int      x_0, x_1;
   int      r_0, r_1;                        /* row offset -> r_0: row_cur % 2, r_1: row_prv % 2 */
   int      c_0, c_1;
   int      r_0b, r_0e, r_1b, r_1e;          /* begin and end indices for row in edgebound list */
   bool     y2_re;

   float    prev_mat, prev_del, prev_ins;    /* temp placeholder sums */
   float    prev_beg, prev_end, prev_esc;    /* temp placeholder sums */
   float    prev_loop, prev_move, prev_sum;  /* temp placeholder sums */
   float    sc, sc1, sc2, sc3, sc4;          /* temp placeholder sums (testing) */
   float    sc_best;                         /* alignment score (return value) */
   float    sc_M, sc_I, sc_D;

   int      d_mod3, d_0, d_1, d_2;           /* d (mod 3) for assigning prev array ptrs */

   /* local or global? (multiple alignments) */
   bool     is_local = target->isLocal;
   float    sc_E = (is_local) ? 0 : -INF;

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

   /* clear leftover data */
   DP_MATRIX_Fill( Q, T, st_MX, sp_MX, -INF );

   /* query sequence */
   seq = query->seq;

   /* Clear top-row */
   row_cur = 0;
   x_0 = 0;             /* current row in matrix */
   r_0 = x_0;           /* for use in linear space alg (mod-mapping) */

   /* Initialize special states (?) */
   XMX(SP_N, x_0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B, x_0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E, x_0) = XMX(SP_C, x_0) = XMX(SP_J, x_0) = -INF;    /* need seq to get here (?)  */

   /* Initialize zero row (top-edge) */
   for (j = 0; j < T; j++) {
      MMX(r_0, j) = IMX(r_0, j) = DMX(r_0, j) = -INF;          /* need seq to get here (?)  */
   }

   /* pass over top-row (i=0) edgebounds from list */
   k    = 0;            /* current index in edgebounds */
   r_0b = 0;            /* beginning index for current row in list */
   while ( k < N && EDG_X(edg,k).id == row_cur ) {
      k++;
   }
   r_0e = k;            /* ending index for current row in list */

   /* init lookback 1 row */
   row_prv = row_cur;
   x_1 = x_0;
   r_1 = r_0;
   r_1b = r_0b;
   r_1e = r_0e;

   /* MAIN RECURSION */
   /* FOR every position in QUERY sequence (row in matrix) */
   for (x_0 = 1; x_0 <= Q; x_0++)
   {
      /* convert quadratic space row index to linear space row index (ex % 2) */
      row_cur = x_0;
      r_0     = x_0;        /* for use in linear space alg (mod-mapping) */

      /* add every edgebound from current row */
      r_0b = k;
      while ( k < N && EDG_X(edg,k).id == row_cur ) {
         k++;
      }
      r_0e = k;

      /* Get next sequence character */
      a = seq[x_1];  /* off-by-one */
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX(r_0, 0) = IMX(r_0, 0) = DMX(r_0, 0) = -INF;
      XMX(SP_E, x_0) = -INF;

      /* FOR every EDGEBOUND in current ROW */
      for (i = r_0b; i < r_0e; i++)
      {
         /* in this context, "diag" represents the "row" */
         x = EDG_X(edg,i).id;                 /* NOTE: this is always the same as cur_row, x_0 */
         y1 = MAX(1, EDG_X(edg,i).lb);        /* can't overflow the left edge */
         y2 = EDG_X(edg,i).rb;
         y2_re = (y2 > T);                      /* check if cloud touches right edge */
         y2 = MIN(y2, T);                       /* can't overflow the right edge */

         /* MAIN RECURSION */
         /* FOR every position in TARGET profile */
         for (j = y1; j < y2; j++)
         {
            c_0 = j;
            c_1 = j-1;

            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            prev_mat = MMX(r_1, c_1)  + TSC(c_1, M2M);
            prev_ins = IMX(r_1, c_1)  + TSC(c_1, I2M);
            prev_del = DMX(r_1, c_1)  + TSC(c_1, D2M);
            prev_beg = XMX(SP_B, x_1) + TSC(c_1, B2M); /* from begin match state (new alignment) */
            /* best-to-match */
            prev_sum = logsum( 
                           logsum( prev_mat, prev_ins ),
                           logsum( prev_del, prev_beg ) );
            MMX(r_0, c_0) = prev_sum + MSC(j,A);

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX(r_1, c_0) + TSC(j,M2I);
            prev_ins = IMX(r_1, c_0) + TSC(j,I2I);
            /* best-to-insert */
            prev_sum = logsum( prev_mat, prev_ins );
            IMX(r_0, c_0) = prev_sum + ISC(c_0, A);

            /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX(r_0, c_1) + TSC(c_1, M2D);
            prev_del = DMX(r_0, c_1) + TSC(c_1, D2D);
            /* best-to-delete */
            prev_sum = logsum(prev_mat, prev_del);
            DMX(r_0, c_0) = prev_sum;

            /* UPDATE E STATE */
            prev_esc = XMX(SP_E, x_0);
            prev_mat = MMX(r_0, j) + sc_E;
            prev_del = DMX(r_0, j) + sc_E;
            XMX(SP_E, x_0) = logsum( 
                                    logsum( prev_mat, prev_del ),
                                       prev_esc );

            #if DEBUG 
            {
               MX_2D(cloud_MX, x_0, c_0) = 1.0;
               MX_3D(test_MX, MAT_ST, x_0, c_0) = MMX(r_0, c_0);
               MX_3D(test_MX, INS_ST, x_0, c_0) = IMX(r_0, c_0);
               MX_3D(test_MX, DEL_ST, x_0, c_0) = DMX(r_0, c_0);
            }
            #endif
         }

         /* UNROLLED FINAL LOOP ITERATION */
         if ( y2_re ) 
         {
            j = T; 
            c_0 = j;
            c_1 = j-1;

            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            prev_mat = MMX(r_1, c_1)  + TSC(c_1, M2M);
            prev_ins = IMX(r_1, c_1)  + TSC(c_1, I2M);
            prev_del = DMX(r_1, c_1)  + TSC(c_1, D2M);
            prev_beg = XMX(SP_B, x_1) + TSC(c_1, B2M);    /* from begin match state (new alignment) */
            /* best-to-match */
            prev_sum = logsum( 
                              logsum( prev_mat, prev_ins ),
                              logsum( prev_del, prev_beg ) );
            MMX(r_0, c_0) = prev_sum + MSC(c_0, A);

            /* FIND SUM OF PATHS TO INSERT STATE */
            IMX(r_0, c_0) = -INF;

            /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX(r_0, c_1) + TSC(c_1, M2D);
            prev_del = DMX(r_0, c_1) + TSC(c_1, D2D);
            /* best-to-delete */
            prev_sum = logsum( prev_mat, prev_del );
            DMX(r_0, c_0) = prev_sum;

            /* UPDATE E STATE */
            prev_esc = XMX(SP_E, x_0);
            prev_mat = MMX(r_0, c_0);
            prev_del = DMX(r_0, c_0);
            XMX(SP_E, x_0) = logsum( 
                                 logsum( prev_mat, prev_del ),
                                 prev_esc );

            #if DEBUG 
            {
               MX_2D(cloud_MX, x_0, c_0) += 1.0;
               MX_3D(test_MX, MAT_ST, x_0, c_0) = MMX(r_0, c_0);
               MX_3D(test_MX, INS_ST, x_0, c_0) = IMX(r_0, c_0);
               MX_3D(test_MX, DEL_ST, x_0, c_0) = DMX(r_0, c_0);
            }
            #endif
         }
      }

      /* ONCE ROW IS COMPLETED, UPDATE SPECIAL STATES */
      /* SPECIAL STATES */
      /* J state */
      sc1 = XMX(SP_J, x_1) + XSC(SP_J, SP_LOOP);   /* J->J */
      sc2 = XMX(SP_E, x_0) + XSC(SP_E, SP_LOOP);   /* E->J is E's "loop" */
      XMX(SP_J, x_0) = logsum( sc1, sc2 );

      /* C state */
      sc1 = XMX(SP_C, x_1) + XSC(SP_C, SP_LOOP);
      sc2 = XMX(SP_E, x_0) + XSC(SP_E, SP_MOVE);
      XMX(SP_C, x_0) = logsum( sc1, sc2 );

      /* N state */
      XMX(SP_N, x_0) = XMX(SP_N, x_1) + XSC(SP_N, SP_LOOP);

      /* B state */
      sc1 = XMX(SP_N, x_0) + XSC(SP_N, SP_MOVE);       /* N->B is N's move */
      sc2 = XMX(SP_J, x_0) + XSC(SP_J, SP_MOVE);       /* J->B is J's move */
      XMX(SP_B, x_0) = logsum( sc1, sc2 );

      /* SET CURRENT ROW TO PREVIOUS ROW */
      row_prv = row_cur;
      x_1  = x_0;
      r_1  = r_0;
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* T state */
   sc_best = XMX(SP_C, Q) + XSC(SP_C, SP_MOVE);
   *sc_final = sc_best;

   return sc_best;
}


/*  
 *  FUNCTION: backward_bounded_Run()
 *  SYNOPSIS: Perform Edge-Bounded Backward part of Cloud Search Algorithm.
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <bnd>       Bounds Data (row-wise)
 *             <sc_final>  Final Score
 *
 *  RETURN: 
 */
float bound_Backward_Quad( const SEQUENCE*      query, 
                           const HMM_PROFILE*   target,
                           const int            Q, 
                           const int            T, 
                           MATRIX_3D*           st_MX,
                           MATRIX_2D*           sp_MX, 
                           EDGEBOUNDS*          edg,
                           float*               sc_final )
{
   FILE* fp = NULL;

   char   a;                              /* store current character in sequence */
   int    A;                              /* store int value of character */
   int    i,j,k = 0;                      /* row, column indices */
   char*  seq;                           /* alias for getting seq */
   int    N = edg->N;                     /* length of edgebound list */

   int    x, y1, y2;                      /* row, leftcol and rightcol bounds in row */
   int    x_0, row_cur, x_1, row_prv;     /* real index of current and previous rows */
   int    r_0, r_1;                       /* row offset -> r_0: row_cur % 2, r_1: row_prv % 2 */
   int    r_0b, r_0e, r_1b, r_1e;         /* begin and end indices for row in edgebound list */

   float  prev_mat, prev_del, prev_ins;   /* temp placeholder sums */
   float  prev_beg, prev_end, prev_sum;   /* temp placeholder sums */
   float  sc, sc1, sc2, sc3, sc4;         /* temp placeholder sums (testing) */
   float  sc_best;                        /* alignment score (return value) */
   float  sc_M, sc_I, sc_D;

   /* local or global? (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* --------------------------------------------------------------------------------- */

   /* clear leftover data */
   DP_MATRIX_Fill( Q, T, st_MX, sp_MX, -INF );

   /* query sequence */
   seq = query->seq;

   /* Initialize the Q row. */
   row_cur = Q;         
   x_0 = Q;             /* current row in matrix */
   r_0 = x_0;           /* for use in linear space alg (mod-mapping) */

   /* Initialize special states */
   XMX(SP_J, x_0) = XMX(SP_B, x_0) = XMX(SP_N, x_0) = -INF;
   XMX(SP_C, x_0) = XSC(SP_C, SP_MOVE);
   XMX(SP_E, x_0) = XMX(SP_C, x_0) + XSC(SP_E, SP_MOVE);

   /* pass over (Q) bottom-row edgebounds from list */
   k    = N-1;          /* current index in edgebounds */
   r_0b = N-1;          /* beginning index for current row in list */
   while ( k >= 0 && EDG_X(edg,k).id == row_cur ) {
      k--;
   }
   r_0e = k;            /* ending index for current row in list */

   /* if there is a bound on row and the right-most bound spans T */
   if ( (r_0b - r_0e > 0) && (EDG_X(edg,r_0b).rb > T) )
   {
      MMX(r_0, T) = DMX(r_0, T) = XMX(SP_E, x_0);
      IMX(r_0, T) = -INF;
   }

   // /* Initialize normal states (naive) */
   // for (j = T-1; j >= 1; j--)
   // {
   //    sc1 = XMX(SP_E, x_0) + sc_E;
   //    sc2 = DMX(r_0, j+1)  + TSC(j, M2D);
   //    MMX(r_0, j) = logsum( XMX(SP_E, x_0) + sc_E, 
   //                               DMX(r_0, j+1)  + TSC(j, M2D) );

   //    sc1 = XMX(SP_E, x_0) + sc_E;
   //    sc2 = DMX(r_0, j+1)  + TSC(j, D2D);
   //    DMX(r_0, j) = logsum( XMX(SP_E, x_0) + sc_E,
   //                               DMX(r_0, j+1)  + TSC(j, D2D) );

   //    // IMX(r_0,j) = -INF;
   // }

   /* Initialize normal states (sparse) */
   for (i = r_0b; i > r_0e; i--) 
   {
      y1 = MAX(0, EDG_X(edg,i).lb);     /* can't overflow the left edge */
      y2 = MIN(EDG_X(edg,i).rb, T);     /* can't overflow the right edge */

      for (j = y2-1; j >= y1; j--)
      {
         sc1 = XMX(SP_E, x_0) + sc_E;
         sc2 = DMX(r_0, j+1)  + TSC(j, M2D);
         MMX(r_0, j) = logsum( XMX(SP_E, x_0) + sc_E, 
                                    DMX(r_0, j+1)  + TSC(j, M2D) );

         sc1 = XMX(SP_E, x_0) + sc_E;
         sc2 = DMX(r_0, j+1)  + TSC(j, D2D);
         DMX(r_0, j) = logsum( XMX(SP_E, x_0) + sc_E,
                                    DMX(r_0, j+1)  + TSC(j, D2D) );

         IMX(r_0,j) = -INF;
      }
   }

   /* init lookback 1 row */
   row_prv = row_cur;
   x_1 = x_0;
   r_1 = r_0;
   r_1b = r_0b;
   r_1e = r_0e;

   /* MAIN RECURSION */
   /* FOR every bound in EDGEBOUND */
   for (x_0 = Q-1; x_0 > 0; x_0--)
   {
      // fprintf(tfp, "(i=%d): ", x_0);
      /* convert quadratic space row index to linear space row index */
      row_cur = x_0;             
      r_0 = x_0;                 /* for use in linear space alg (mod-mapping) */

      /* add every edgebound from current row */
      r_0b = k;
      // printf("k: %d => row_cur: %d, k.id: %d\n", k, row_cur, EDG_X(edg,k).id);
      while ( k >= 0 && EDG_X(edg,k).id >= row_cur ) {
         k--;
      }
      r_0e = k;

      /* Get next sequence character */
      a = seq[x_0];
      A = AA_REV[a];

      /* UPDATE SPECIAL STATES at the start of EACH ROW */

      /* B STATE (NAIVE) */
      // XMX(SP_B, x_0) = -INF;
      // for (j = 1; j <= T; j++) {
      //    XMX(SP_B, x_0) = logsum( XMX(SP_B, x_0),
      //                                  MMX(r_1, j) + TSC(j-1, B2M) + MSC(j, A) );
      // }

      /* B STATE (SPARSE) */
      XMX(SP_B, x_0) = -INF;
      for (i = r_1b; i > r_1e; i--) 
      {
         y1 = MAX(1, EDG_X(edg,i).lb);         /* can't overflow the left edge */
         y2 = MIN(EDG_X(edg,i).rb, T);         /* can't overflow the right edge */

         for (j = y2-1; j >= y1; j--)
         {
            XMX(SP_B, x_0) = logsum( XMX(SP_B, x_0),
                                          MMX(r_1, j) + TSC(j-1, B2M) + MSC(j, A) );
         }
      }

      XMX(SP_J, x_0) = logsum( XMX(SP_J, x_1) + XSC(SP_J, SP_LOOP),
                                    XMX(SP_B, x_0) + XSC(SP_J, SP_MOVE) );

      XMX(SP_C, x_0) = XMX(SP_C, x_1) + XSC(SP_C, SP_LOOP);

      XMX(SP_E, x_0) = logsum( XMX(SP_J, x_0) + XSC(SP_E, SP_LOOP),
                                    XMX(SP_C, x_0) + XSC(SP_E, SP_MOVE) );

      XMX(SP_N, x_0) = logsum( XMX(SP_N, x_1) + XSC(SP_N, SP_LOOP),
                                    XMX(SP_B, x_0) + XSC(SP_N, SP_MOVE) );

      /* if there is a bound on row and the right-most bound spans T */
      if ( (r_0b - r_0e > 0) && (EDG_X(edg,r_0b).rb > T) )
      {
         MMX(r_0, T) = DMX(r_0, T) = XMX(SP_E, x_0);
         IMX(r_0, T) = -INF;
      }

      /* FOR every EDGEBOUND in current ROW */
      for (i = r_0b; i > r_0e; i--)
      {
         /* in this context, "diag" represents the "row" */
         // x_0  = EDG_X(edg,i).id;
         y1 = MAX(0, EDG_X(edg,i).lb);     /* can't overflow the left edge */
         y2 = MIN(EDG_X(edg,i).rb, T);     /* can't overflow the right edge */

         /* FOR every position in TARGET profile */
         for (j = y2-1; j >= y1; j--)
         {
            sc_M = MSC(j+1,A);
            sc_I = ISC(j,A);

            /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
            prev_mat = MMX(r_1, j+1)  + TSC(j, M2M) + sc_M;
            prev_ins = IMX(r_1,   j)  + TSC(j, M2I) + sc_I;
            prev_del = DMX(r_0, j+1)  + TSC(j, M2D);
            prev_end = XMX(SP_E, x_0) + sc_E;     /* from end match state (new alignment) */
            /* best-to-match */
            prev_sum = logsum( 
                              logsum( prev_mat, prev_ins ),
                              logsum( prev_end, prev_del ) );
            MMX(r_0, j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            prev_mat = MMX(r_1, j+1) + TSC(j, I2M) + sc_M;
            prev_ins = IMX(r_1,   j) + TSC(j, I2I) + sc_I;
            /* best-to-insert */
            prev_sum = logsum( prev_mat, prev_ins );
            IMX(r_0, j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            prev_mat = MMX(r_1, j+1)  + TSC(j,D2M) + sc_M;
            prev_del = DMX(r_0, j+1)  + TSC(j,D2D);
            prev_end = XMX(SP_E,x_0)  + sc_E;
            /* best-to-delete */
            prev_sum = logsum( 
                              prev_mat,
                              logsum( prev_del, prev_end ) );
            DMX(r_0, j) = prev_sum;
         }
      }

      /* SET CURRENT ROW TO PREVIOUS ROW */
      row_prv = row_cur;
      x_1  = x_0;
      r_1  = r_0;
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* FINAL ROW */
   row_cur = x_0 = i = 0;
   r_0 = x_0;

   /* pass over (Q) bottom-row edgebounds from list */
   r_0b = k;            /* beginning index for current row in list */
   while ( k >= 0 && EDG_X(edg,k).id == row_cur ) {
      k--;
   }
   r_0e = k;            /* ending index for current row in list */

   /* FINAL i = 0 row */
   a = seq[1];
   A = AA_REV[a];

   /* B STATE (NAIVE) */
   // XMX(SP_B, x_0) = -INF;
   // for (j = 1; j >= T; j++) {
   //    XMX(SP_B,x_0) = logsum( XMX(SP_B,x_0),
   //                                 MMX(r_1,j) + TSC(j-1,B2M) + MSC(j,A) );
   // }

   /* B STATE (SPARSE) */
   XMX(SP_B, x_0) = -INF;
   for (i = r_1b; i > r_1e; i--) {
      y1 = MAX(1, EDG_X(edg,i).lb);         /* can't overflow the left edge */
      y2 = MIN(EDG_X(edg,i).rb, T);         /* can't overflow the right edge */

      for (j = y2-1; j >= y1; j--)
      {
         XMX(SP_B, x_0) = logsum( XMX(SP_B, x_0),
                                       MMX(r_1, j) + TSC(j-1, B2M) + MSC(j, A) );
      }
   }

   XMX(SP_J, r_0) = -INF;
   XMX(SP_C, r_0) = -INF;
   XMX(SP_E, r_0) = -INF;

   XMX(SP_N, x_0) = logsum( XMX(SP_N, x_1) + XSC(SP_N, SP_LOOP),
                                 XMX(SP_B, x_0) + XSC(SP_N, SP_MOVE) );

   // for (j = T; j >= 1; j--) {
   //    MMX(x_0, j) = IMX(x_0, j) = DMX(x_0, j) = -INF;
   // }

   sc_best = XMX(SP_N,0);
   *sc_final = sc_best;
   // fclose(tfp);
   return sc_best;
}
