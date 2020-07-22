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
 *  FUNCTION: forward_Quad()
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
   float    prev_loop, prev_move, prev_sum;  /* temp placeholder sums */
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

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   /* query sequence */
   seq         = query->seq;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   /* first row */
   x_0 = 0;
   r_0 = 0;

   /* initialize special states (?) */
   XMX(SP_N, x_0) = 0;                                            /* S->N, p=1             */
   XMX(SP_B, x_0) = XSC(SP_N, SP_MOVE);                           /* S->N->B, no N-tail    */
   XMX(SP_E, x_0) = XMX(SP_C, x_0) = XMX(SP_J, x_0) = -INF;       /* need seq to get here (?)  */

   /* initialize 0 row (top-edge) */
   for (j = 0; j < T; j++) { 
      c_0 = j;
      MMX(r_0, c_0) = IMX(r_0, c_0) = DMX(r_0, c_0) = -INF;    /* need seq to get here (?)  */
   }

   /* FOR every position in QUERY seq */
   for (i = 1; i <= Q; i++)
   {  
      x_0 = i;
      x_1 = i-1;
      r_0 = x_0;
      r_1 = x_1;

      c_0 = 0;

      /* Get next sequence character */
      a = seq[x_1];
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX(r_0, c_0) = IMX(r_0, c_0) = DMX(r_0, c_0) = -INF;
      XMX(SP_E, x_0) = -INF;

      /* MAIN RECURSION */
      /* FOR every position in TARGET profile */
      for (j = 1; j < T; j++)
      {
         c_0 = j;
         c_1 = j-1;

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         prev_mat = MMX(r_1, c_1)  + TSC(c_1, M2M);
         prev_ins = IMX(r_1, c_1)  + TSC(c_1, I2M);
         prev_del = DMX(r_1, c_1)  + TSC(c_1, D2M);
         prev_beg = XMX(SP_B, x_1)  + TSC(c_1, B2M); /* from begin match state (new alignment) */
         /* best-to-match */
         prev_sum = logsum( 
                        logsum( prev_mat, prev_ins ),
                        logsum( prev_beg, prev_del ) );
         MMX(r_0, c_0) = prev_sum + MSC(c_0, A);

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the previous row (upper) of each state) */
         prev_mat = MMX(r_1, c_0) + TSC(c_0, M2I);
         prev_ins = IMX(r_1, c_0) + TSC(c_0, I2I);
         /* best-to-insert */
         prev_sum = logsum( prev_mat, prev_ins );
         IMX(r_0, c_0) = prev_sum + ISC(c_0, A);

         /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
         /* previous states (match takes the previous column (left) of each state) */
         prev_mat = MMX(r_0, c_1) + TSC(c_1, M2D);
         prev_del = DMX(r_0, c_1) + TSC(c_1, D2D);
         /* best-to-delete */
         prev_sum = logsum( prev_mat, prev_del );
         DMX(r_0, c_0) = prev_sum;

         /* UPDATE E STATE */
         prev_mat = MMX(r_0, c_0) + sc_E;
         prev_del = DMX(r_0, c_0) + sc_E;
         /* best-to-e-state */
         prev_esc = XMX(SP_E, x_0);
         XMX(SP_E, x_0) = logsum( 
                              logsum( prev_mat, prev_del ),
                              prev_esc );

         /* embed linear row into quadratic test matrix */
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
      j = T; 
      c_0 = j;
      c_1 = j-1;

      /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
      /* best previous state transition (match takes the diag element of each prev state) */
      prev_mat = MMX(r_1, c_1)  + TSC(c_1, M2M);
      prev_ins = IMX(r_1, c_1)  + TSC(c_1, I2M);
      prev_del = DMX(r_1, c_1)  + TSC(c_1, D2M);
      prev_beg = XMX(SP_B, x_1)  + TSC(c_1, B2M);    /* from begin match state (new alignment) */
      /* sum-to-match */
      prev_sum = logsum( 
                     logsum( prev_mat, prev_ins ),
                     logsum( prev_del, prev_beg ) );
      MMX(r_0, c_0) = prev_sum + MSC(c_0, A);

      /* FIND SUM OF PATHS TO INSERT STATE (unrolled) */
      IMX(r_0, c_0) = -INF;

      /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) (unrolled) */
      /* previous states (match takes the left element of each state) */
      prev_mat = MMX(r_0, c_1) + TSC(c_1, M2D);
      prev_del = DMX(r_0, c_1) + TSC(c_1, D2D);
      /* sum-to-delete */
      prev_sum = logsum( prev_mat, prev_del );
      DMX(r_0, c_0) = prev_sum;

      /* UPDATE E STATE (unrolled) */
      prev_esc = XMX(SP_E, x_0);
      prev_mat = MMX(r_0, c_0);
      prev_del = DMX(r_0, c_0);
      /* best-to-begin */
      XMX(SP_E, x_0) = logsum( 
                           logsum( prev_del, prev_mat ),
                           prev_esc );

      /* SPECIAL STATES */
      /* J state */
      sc_1 = XMX(SP_J, x_1) + XSC(SP_J, SP_LOOP);       /* J->J */
      sc_2 = XMX(SP_E, x_0) + XSC(SP_E, SP_LOOP);       /* E->J is E's "loop" */
      XMX(SP_J, x_0) = logsum( sc_1, sc_2 );         

      /* C state */
      sc_1 = XMX(SP_C, x_1) + XSC(SP_C, SP_LOOP);
      sc_2 = XMX(SP_E, x_0) + XSC(SP_E, SP_MOVE);
      XMX(SP_C, x_0) = logsum( sc_1, sc_2 );

      /* N state */
      XMX(SP_N, x_0) = XMX(SP_N, x_1) + XSC(SP_N, SP_LOOP);

      /* B state */
      sc_1 = XMX(SP_N, x_0) + XSC(SP_N, SP_MOVE);         /* N->B is N's move */
      sc_2 = XMX(SP_J, x_0) + XSC(SP_J, SP_MOVE);         /* J->B is J's move */
      XMX(SP_B, x_0) = logsum( sc_1, sc_2 );     

      /* embed linear row into quadratic test matrix */
      #if DEBUG
      {
         MX_2D(cloud_MX, x_0, c_0) = 1.0;
         MX_3D(test_MX, MAT_ST, x_0, c_0) = MMX(r_0, c_0);
         MX_3D(test_MX, INS_ST, x_0, c_0) = IMX(r_0, c_0);
         MX_3D(test_MX, DEL_ST, x_0, c_0) = DMX(r_0, c_0);
      }
      #endif
   }

   /* T state */
   sc_best     = XMX(SP_C, Q) + XSC(SP_C, SP_MOVE);
   *sc_final   = sc_best; 

   /* flag matrices that they contain dirty values (not -INF) */
   st_MX->clean = false;
   sp_MX->clean = false;

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

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   /* query sequence */
   seq         = query->seq;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   /* Initialize the Q row. */
   x_0 = Q;
   r_0 = x_0;
   c_0 = T;

   XMX(SP_J, x_0) = XMX(SP_B, x_0) = XMX(SP_N, x_0) = -INF;
   XMX(SP_C, x_0) = XSC(SP_C, SP_MOVE);
   XMX(SP_E, x_0) = XMX(SP_C, x_0) + XSC(SP_E, SP_MOVE);

   MMX(r_0, c_0) = DMX(r_0, c_0) = XMX(SP_E, x_0);
   IMX(r_0, c_0) = -INF;

   for (j = T-1; j >= 1; j--)
   {
      c_0 = j;
      c_1 = j+1;

      prev_esc = XMX(SP_E, Q) + sc_E;
      prev_del = DMX(r_0, c_1)  + TSC(c_0, M2D);
      MMX(r_0, c_0) = logsum( prev_esc, prev_del );

      prev_esc = XMX(SP_E, Q) + sc_E;
      prev_del = DMX(r_0, c_1)  + TSC(c_0, D2D);
      DMX(r_0, c_0) = logsum( prev_esc, prev_del );

      IMX(r_0, c_0) = -INF;
   }

   /* MAIN RECURSION */
   /* FOR every position in QUERY seq */
   for (i = Q-1; i >= 1; i--)
   {
      x_0 = i;
      x_1 = i+1;
      r_0 = x_0;
      r_1 = x_1;

      j = 0;
      c_0 = j;
      c_1 = j+1;

      /* Get next sequence character */
      a = seq[x_0];
      A = AA_REV[a];

      /* SPECIAL STATES */

      /* B -> MATCH */
      /* NOTE: When j = 0, MMX and MSC do not match HMMER p7_GBackward() implementation.   */
      XMX(SP_B, x_0) = MMX(r_1, 1) + TSC(0, B2M) + MSC(1, A);
      
      for (j = 2; j <= T; j++) {
         c_0 = j;
         c_1 = j-1;
         prev_sum = XMX(SP_B, x_0);
         prev_mat = MMX(r_1, c_0) + TSC(c_1, B2M) + MSC(c_0, A);
         XMX(SP_B, x_0) = logsum( prev_sum, prev_mat);
      }

      XMX(SP_J, x_0) = logsum( XMX(SP_J, x_1) + XSC(SP_J, SP_LOOP),
                                 XMX(SP_B, x_0) + XSC(SP_J, SP_MOVE) );

      XMX(SP_C, x_0) = XMX(SP_C, x_1) + XSC(SP_C, SP_LOOP);

      XMX(SP_E, x_0) = logsum( XMX(SP_J, x_0) + XSC(SP_E, SP_LOOP),
                                 XMX(SP_C, x_0) + XSC(SP_E, SP_MOVE) );

      XMX(SP_N, x_0) = logsum( XMX(SP_N, x_1) + XSC(SP_N, SP_LOOP),
                                 XMX(SP_B, x_0)   + XSC(SP_N, SP_MOVE) );

      MMX(r_0, T) = DMX(r_0, T) = XMX(SP_E, x_0);
      IMX(r_0, T) = -INF;

      #if DEBUG 
      {
         MX_2D(cloud_MX, x_0, c_0) = 1.0;
         MX_3D(test_MX, MAT_ST, x_0, c_0) = MMX(r_0, c_0);
         MX_3D(test_MX, INS_ST, x_0, c_0) = IMX(r_0, c_0);
         MX_3D(test_MX, DEL_ST, x_0, c_0) = DMX(r_0, c_0);
      }
      #endif

      /* FOR every position in TARGET profile */
      for (j = T-1; j >= 1; j--)
      {
         c_0 = j;
         c_1 = j+1;

         /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
         prev_mat = MMX(r_1, c_1) + TSC(c_0, M2M) + MSC(c_1, A);
         prev_ins = IMX(r_1, c_0) + TSC(c_0, M2I) + ISC(c_1, A);
         prev_del = DMX(r_0, c_1) + TSC(c_0, M2D);
         prev_end = XMX(SP_E, x_0) + sc_E;     /* from end match state (new alignment) */
         /* best-to-match */
         prev_sum = logsum( 
                        logsum( prev_mat, prev_ins ),
                        logsum( prev_end, prev_del ) );
         MMX(r_0, c_0) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
         prev_mat = MMX(r_1, c_1) + TSC(c_0, I2M) + MSC(c_1, A);
         prev_ins = IMX(r_1, c_0) + TSC(c_0, I2I) + ISC(c_0, A);
         /* best-to-insert */
         prev_sum = logsum( prev_mat, prev_ins );
         IMX(r_0, c_0) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
         prev_mat = MMX(r_1, c_1) + TSC(c_0, D2M) + MSC(c_1, A);
         prev_del = DMX(r_0, c_1) + TSC(c_0, D2D);
         prev_end = XMX(SP_E, x_0) + sc_E;
         /* best-to-delete */
         prev_sum = logsum( prev_mat, 
                        logsum( prev_del, prev_end ) );
         DMX(r_0, c_0) = prev_sum;

         #if DEBUG 
         {
            MX_2D(cloud_MX, x_0, c_0) = 1.0;
            MX_3D(test_MX, MAT_ST, x_0, c_0) = MMX(r_0, c_0);
            MX_3D(test_MX, INS_ST, x_0, c_0) = IMX(r_0, c_0);
            MX_3D(test_MX, DEL_ST, x_0, c_0) = DMX(r_0, c_0);
         }
         #endif
      }
   }

   /* FINAL ROW (i = 0) */
   /* At i=0, only N,B states are reachable. */
   i = 0;
   x_0 = i;
   x_1 = i+1;
   r_0 = x_0;
   r_1 = x_1;

   j = 0;
   c_0 = j;
   c_1 = j+1;

   a = seq[x_0];
   A = AA_REV[a];

   /* t_BM index is 0 because it's stored off-by-one. */
   XMX(SP_B, 0) = MMX(1, 1) + TSC(0, B2M) + MSC(1, A);

   for (j = 2; j <= T; j++) {
      c_0 = j;
      c_1 = j-1;
      prev_sum = XMX(SP_B, x_0);
      prev_mat = MMX(x_1, c_0) + TSC(c_1, B2M) + MSC(c_0, A);
      XMX(SP_B, x_0) = logsum( prev_sum, prev_mat );
   }

   XMX(SP_J, x_0) = -INF;
   XMX(SP_C, x_0) = -INF;
   XMX(SP_E, x_0) = -INF;

   XMX(SP_N, x_0) = logsum( XMX(SP_N, x_1) + XSC(SP_N,SP_LOOP),
                              XMX(SP_B, x_0) + XSC(SP_N,SP_MOVE) );

   for (j = T; j >= 1; j--) {
      c_0 = j;
      MMX(r_0, c_0) = IMX(r_0, c_0) = DMX(r_0, c_0) = -INF;
   }

   #if DEBUG 
   {
      MX_2D(cloud_MX, x_0, c_0) = 1.0;
      MX_3D(test_MX, MAT_ST, x_0, c_0) = MMX(r_0, c_0);
      MX_3D(test_MX, INS_ST, x_0, c_0) = IMX(r_0, c_0);
      MX_3D(test_MX, DEL_ST, x_0, c_0) = DMX(r_0, c_0);
   }
   #endif

   sc_best = XMX(SP_N,0);
   *sc_final = sc_best;

   /* flag matrices that they contain dirty values (not -INF) */
   st_MX->clean = false;
   sp_MX->clean  = false;

   return STATUS_SUCCESS;
}

