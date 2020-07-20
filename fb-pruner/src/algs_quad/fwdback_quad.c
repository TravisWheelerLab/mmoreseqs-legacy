/*******************************************************************************
 *  FILE:      fwdback_quad.c
 *  PURPOSE:   The Forward-Backward Algorithm for Sequence Alignment Search.
 *             (Quadratic Space Alg)
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
#include "fwdback_quad.h"


/*  
 *  FUNCTION: forward_Run()
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
int forward_Quad(    const SEQUENCE*    query, 
                     const HMM_PROFILE* target, 
                     const int          Q, 
                     const int          T, 
                     MATRIX_3D*         st_MX, 
                     MATRIX_2D*         sp_MX,
                     float*             sc_final )
{
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   int      b, d, i, j, k;                   /* diagonal, row, column indices */
   char*    seq;                             /* alias for getting seq */
   int      N;                               /* length of edgebound list */
   bool     is_local;                        /* whether using local or global alignments */

   int      row_cur, row_prv;                /* current and previous rows */
   int      x, y1, y2;                       /* row, leftcol and rightcol bounds in row */
   int      x_0, x_1;                        /* real index of current and previous rows */
   int      r_0, r_1;                        /* row offset -> r_0: row_cur % 2, r_1: row_prv % 2 */
   int      c_0, c_1;                        /* real index of current and previous columns */
   int      r_0b, r_0e, r_1b, r_1e;          /* begin and end indices for row in edgebound list */
   int      d_0, d_1, d_2;                   /* d (mod 3) for assigning diag array ptrs */
   bool     y2_re;                           /* checks if edge touches rightbound */

   float    prev_mat, prev_del, prev_ins;    /* temp placeholder sums */
   float    prev_beg, prev_end, prev_esc;    /* temp placeholder sums */
   float    prev_loop, prev_sum;             /* temp placeholder sums */
   float    sc, sc_1, sc_2, sc_3, sc_4;      /* temp placeholder sums */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, e-state scores */
   float    sc_best;

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_2D*  cloud_MX3;
   MATRIX_3D*  test_MX;
   MATRIX_3D*  test_MX3;
   int         num_writes;
   int         num_clears;

   /* initialize debugging matrix */
   #if DEBUG
   {
      cloud_MX    = debugger->cloud_MX;
      cloud_MX3   = debugger->cloud_MX3;
      test_MX     = debugger->test_MX;
      test_MX3    = debugger->test_MX3;

      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
      MATRIX_2D_Fill( cloud_MX3, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
      MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
      MATRIX_3D_Fill( test_MX3, -INF );

      num_writes = 0;
      num_clears = 0;
   }
   #endif

   /* verify memory is clean (all cells set to -INF) */
   // #if MEMCHECK
   // {
   //    int cmp =  MATRIX_3D_Check_Clean( st_MX );
   //    printf("PRE-CHECK CLEAN -> FORWARD?\t%d\n", cmp);
   //    if ( cmp != 0 ) {
   //       MATRIX_3D_Clean( st_MX );
   //    }
   // }
   // #endif 

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   /* query sequence */
   seq = query->seq;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   /* clear leftover data */
   if ( st_MX->clean == false ) {
      DP_MATRIX_Clean( Q, T, st_MX, sp_MX );
      st_MX->clean = true;
   }
   st_MX->clean = true;

   /* initialize special states (?) */
   XMX(SP_N,0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E,0) = XMX(SP_C,0) = XMX(SP_J,0) = -INF;          /* need seq to get here (?)  */

   /* initialize 0 row (top-edge) */
   for (j = 0; j < T; j++)
      MMX(0,j) = IMX(0,j) = DMX(0,j) = -INF;             /* need seq to get here (?)  */

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
      }

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
                     logsum( prev_del, prev_beg ) );
      MMX(i,j) = prev_sum + MSC(j,A);

      /* FIND SUM OF PATHS TO INSERT STATE (unrolled) */
      IMX(i,j) = -INF;

      /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) (unrolled) */
      /* previous states (match takes the left element of each state) */
      prev_mat = MMX(i,j-1) + TSC(j-1,M2D);
      prev_del = DMX(i,j-1) + TSC(j-1,D2D);
      /* sum-to-delete */
      prev_sum = logsum( prev_mat, prev_del );
      DMX(i,j) = prev_sum;

      /* UPDATE E STATE (unrolled) */
      prev_esc = XMX(SP_E,i);
      prev_mat = MMX(i,j);
      prev_del = DMX(i,j);
      /* best-to-begin */
      XMX(SP_E,i) = logsum( 
                        logsum( prev_del, prev_mat ),
                        prev_esc );

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
   }

   /* T state */
   sc_best = XMX(SP_C,Q) + XSC(SP_C,SP_MOVE);
   *sc_final = sc_best; 

   return STATUS_SUCCESS;
}

/* FUNCTION: backward_Run()
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
int backward_Quad(   const SEQUENCE*    query, 
                     const HMM_PROFILE* target, 
                     const int          Q, 
                     const int          T, 
                     MATRIX_3D*         st_MX, 
                     MATRIX_2D*         sp_MX,
                     float*             sc_final)
{
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   int      b, d, i, j, k;                   /* diagonal, row, column indices */
   char*    seq;                             /* alias for getting seq */
   int      N;                               /* length of edgebound list */
   bool     is_local;                        /* whether using local or global alignments */

   int      row_cur, row_prv;                /* current and previous rows */
   int      x, y1, y2;                       /* row, leftcol and rightcol bounds in row */
   int      x_0, x_1;                        /* real index of current and previous rows */
   int      r_0, r_1;                        /* row offset -> r_0: row_cur % 2, r_1: row_prv % 2 */
   int      c_0, c_1;                        /* real index of current and previous columns */
   int      r_0b, r_0e, r_1b, r_1e;          /* begin and end indices for row in edgebound list */
   int      d_0, d_1, d_2;                   /* d (mod 3) for assigning diag array ptrs */
   bool     y2_re;                           /* checks if edge touches rightbound */

   float    prev_mat, prev_del, prev_ins;    /* temp placeholder sums */
   float    prev_beg, prev_end, prev_esc;    /* temp placeholder sums */
   float    prev_loop, prev_sum;             /* temp placeholder sums */
   float    sc, sc_1, sc_2, sc_3, sc_4;      /* temp placeholder sums */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, e-state scores */
   float    sc_best;

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_2D*  cloud_MX3;
   MATRIX_3D*  test_MX;
   MATRIX_3D*  test_MX3;
   int         num_writes;
   int         num_clears;

   /* initialize debugging matrix */
   #if DEBUG
   {
      cloud_MX    = debugger->cloud_MX;
      cloud_MX3   = debugger->cloud_MX3;
      test_MX     = debugger->test_MX;
      test_MX3    = debugger->test_MX3;

      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
      MATRIX_2D_Fill( cloud_MX3, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
      MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
      MATRIX_3D_Fill( test_MX3, -INF );

      num_writes = 0;
      num_clears = 0;
   }
   #endif

   /* verify memory is clean (all cells set to -INF) */
   #if MEMCHECK
   {
      int cmp =  MATRIX_3D_Check_Clean( st_MX );
      printf("PRE-CHECK CLEAN -> FORWARD?\t%d\n", cmp);
      if ( cmp != 0 ) {
         MATRIX_3D_Clean( st_MX );
      }
   }
   #endif 

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   /* query sequence */
   seq = query->seq;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   /* clear leftover data */
   if ( st_MX->clean == false ) {
      DP_MATRIX_Clean( Q, T, st_MX, sp_MX );
      st_MX->clean = true;
   }
   st_MX->clean = true;

   /* Initialize the Q row. */
   row_cur = Q;         
   x_0 = Q;                /* current row in matrix */
   r_0 = x_0 % 2;          /* for use in linear space alg (mod-mapping) */

   XMX(SP_J,Q) = XMX(SP_B,Q) = XMX(SP_N,Q) = -INF;
   XMX(SP_C,Q) = XSC(SP_C,SP_MOVE);
   XMX(SP_E,Q) = XMX(SP_C,Q) + XSC(SP_E,SP_MOVE);

   MMX(Q,T) = DMX(Q,T) = XMX(SP_E,Q);
   IMX(Q,T) = -INF;

   for (j = T-1; j >= 1; j--)
   {
      c_0 = j;
      c_1 = j+1;

      prev_esc = XMX(SP_E,Q) + sc_E;
      prev_del = DMX(Q,j+1)  + TSC(j,M2D);
      MMX(Q,j) = logsum( prev_esc, prev_del );

      prev_esc = XMX(SP_E,Q) + sc_E;
      prev_del = DMX(Q,j+1)  + TSC(j,D2D);
      DMX(Q,j) = logsum( prev_esc, prev_del );

      IMX(Q,j) = -INF;
   }

   /* MAIN RECURSION */
   /* FOR every position in QUERY seq */
   for (i = Q-1; i >= 1; i--)
   {
      x_0 = i;
      x_1 = i-1;
      r_0 = x_0 % 2;
      r_1 = x_1 % 2;

      j = 1;
      c_0 = j;
      c_1 = j-1;

      /* Get next sequence character */
      a = seq[x_0];
      A = AA_REV[a];

      /* SPECIAL STATES */
      XMX(SP_B,i) = MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A);

      /* B -> MATCH */
      for (j = 2; j <= T; j++) {
         XMX(SP_B,i) = logsum( XMX(SP_B,i),
                                    MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
      }

      XMX(SP_J,i) = logsum( XMX(SP_J,i+1) + XSC(SP_J,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_J,SP_MOVE) );

      XMX(SP_C,i) = XMX(SP_C,i+1) + XSC(SP_C,SP_LOOP);

      XMX(SP_E,i) = logsum( XMX(SP_J,i) + XSC(SP_E,SP_LOOP),
                                 XMX(SP_C,i) + XSC(SP_E,SP_MOVE) );

      XMX(SP_N,i) = logsum( XMX(SP_N,i+1) + XSC(SP_N,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_N,SP_MOVE) );

      MMX(i,T) = DMX(i,T) = XMX(SP_E,i);
      IMX(i,T) = -INF;

      x_0 = i;
      /* FOR every position in TARGET profile */
      for (j = T-1; j >= 1; j--)
      {
         sc_M = MSC(j+1,A);
         sc_I = ISC(j,A);

         /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
         prev_mat = MMX(i+1,j+1) + TSC(j,M2M) + sc_M;
         prev_ins = IMX(i+1,j)   + TSC(j,M2I) + sc_I;
         prev_del = DMX(i,j+1)   + TSC(j,M2D);
         prev_end = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
         /* best-to-match */
         prev_sum = logsum( 
                        logsum( prev_mat, prev_ins ),
                        logsum( prev_end, prev_del ) );
         MMX(i,j) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
         prev_mat = MMX(i+1,j+1) + TSC(j,I2M) + sc_M;
         prev_ins = IMX(i+1,j)   + TSC(j,I2I) + sc_I;
         /* best-to-insert */
         prev_sum = logsum( prev_mat, prev_ins );
         IMX(i,j) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
         prev_mat = MMX(i+1,j+1) + TSC(j,D2M) + sc_M;
         prev_del = DMX(i,  j+1) + TSC(j,D2D);
         prev_end = XMX(SP_E,i)  + sc_E;
         /* best-to-delete */
         prev_sum = logsum( prev_mat, 
                        logsum( prev_del, prev_end ) );
         DMX(i,j) = prev_sum;
      }
   }

   /* FINAL ROW (i = 0) */
   /* At i=0, only N,B states are reachable. */
   a = seq[1];
   A = AA_REV[a];

   /* t_BM index is 0 because it's stored off-by-one. */
   XMX(SP_B,0) = MMX(1,1) + TSC(0,B2M) + MSC(1,A);

   for (j = 2; j >= T; j++) {
      XMX(SP_B,0) = logsum( XMX(SP_B,0),
                                 MMX(1,j) + TSC(j-1,B2M) + MSC(j,A) );
   }

   XMX(SP_J,i) = -INF;
   XMX(SP_C,i) = -INF;
   XMX(SP_E,i) = -INF;

   XMX(SP_N,i) = logsum( XMX(SP_N,1) + XSC(SP_N,SP_LOOP),
                              XMX(SP_B,0) + XSC(SP_N,SP_MOVE) );

   for (j = T; j >= 1; j--) {
      MMX(i,j) = IMX(i,j) = DMX(i,j) = -INF;
   }

   #if DEBUG 
   {
      
   }
   #endif

   sc_best = XMX(SP_N,0);
   *sc_final = sc_best;
   return sc_best;
}








