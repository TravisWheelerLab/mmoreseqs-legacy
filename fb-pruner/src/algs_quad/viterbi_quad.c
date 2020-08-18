/*******************************************************************************
 *  FILE:      viterbi_quad.c
 *  PURPOSE:   The Viterbi Algorithm and Traceback for Sequence Alignment Search.
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
#include "viterbi_quad.h"

/*
 *  FUNCTION:  run_Viterbi_Quad()
 *  SYNOPSIS:  Run Viterbi Algorithm (Seq-to-Profile, general unoptimized)
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Viterbi_Quad(   const SEQUENCE*    query,        /* query sequence */
                        const HMM_PROFILE* target,       /* target hmm model */
                        const int          Q,            /* query length */
                        const int          T,            /* target length */
                        MATRIX_3D*         st_MX,        /* normal matrix */
                        MATRIX_2D*         sp_MX,        /* special matrix */
                        float*             sc_final )    /* final max score */
{
   /* vars for accessing query/target data structs */
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   char*    seq;                             /* alias for getting seq */
   int      N;                               /* length of edgebound list */
   bool     is_local;                        /* whether using local or global alignments */

   /* vars for indexing into data matrices by row-col */
   int      b, d, i, j, k;                   /* antidiagonal, row, column indices */
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* mod mapping of column index into data matrix (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */

   /* vars for indexing into data matrices by anti-diag */
   int      d_0, d_1, d_2;                   /* real index of current and previous antidiagonals */
   int      dx0, dx1, dx2;                   /* mod mapping of antidiagonal index into data matrix */
   int      k_0, k_1;                        /* offset into antidiagonal */
   int      d_st, d_end, d_cnt;              /* starting and ending diagonal indices */
   int      dim_T, dim_Q, dim_TOT;           /* dimensions of submatrix being searched */
   int      dim_min, dim_max;                /* diagonal index where num cells reaches highest point and diminishing point */ 
   int      num_cells;                       /* number of cells in current diagonal */

   /* vars for indexing into edgebound lists */
   BOUND*   bnd;                             /* current bound */
   int      id;                              /* id in edgebound list (row/diag) */
   int      r_0;                             /* current index for current row */
   int      r_0b, r_0e;                      /* begin and end indices for current row in edgebound list */
   int      r_1;                             /* current index for previous row */
   int      r_1b, r_1e;                      /* begin and end indices for current row in edgebound list */
   int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   bool     rb_T;                            /* checks if edge touches right bound of matrix */

   /* vars for recurrance scores */
   float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;            /* previous loop and move for special states */
   float    prv_sum, prv_best;             /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */
   
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

   /* query sequence */
   seq         = query->seq;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   q_0 = 0;
   qx0 = q_0;

   /* initialize special states (?) */
   XMX(SP_N, q_0) = 0;                                        /* S->N, p=1             */
   XMX(SP_B, q_0) = XSC(SP_N, SP_MOVE);                       /* S->N->B, no N-tail    */
   XMX(SP_E, q_0) = XMX(SP_C, q_0) = XMX(SP_J, q_0) = -INF;   /* need seq to get here  */

   /* initialize zero row (top-edge) */
   for (t_0 = 0; t_0 <= T; t_0++) { 
      MMX(qx0, t_0) = IMX(qx0, t_0) = DMX(qx0, t_0) = -INF;          /* need seq to get here  */
   }

   /* FOR every position in QUERY seq */
   for (q_0 = 1; q_0 <= Q; q_0++)
   {
      q_1 = q_0-1;
      qx0 = q_0;
      qx1 = q_1;

      t_0 = 0;

      /* Get next character in Query */
      a = seq[q_1];
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX(qx0, t_0) = IMX(qx0, t_0) = DMX(qx0, t_0) = -INF;
      XMX(SP_E, q_0) = -INF;

      /* FOR every position in TARGET profile */
      for (t_0 = 1; t_0 < T; t_0++)
      {
         t_1 = t_0-1;

         /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         prv_M = MMX(qx1, t_1) + TSC(t_1, M2M);
         prv_I = IMX(qx1, t_1) + TSC(t_1, I2M);
         prv_D = DMX(qx1, t_1) + TSC(t_1, D2M);
         prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M); /* from begin match state (new alignment) */
         /* best-to-match */
         prv_best = calc_Max( 
                        calc_Max( prv_M, prv_I ), 
                        calc_Max( prv_D, prv_B ) );
         MMX(qx0, t_0)  = prv_best + MSC(t_0, A);

         /* FIND BEST PATH TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         prv_M = MMX(qx1, t_0) + TSC(t_0, M2I);
         prv_I = IMX(qx1, t_0) + TSC(t_0, I2I);
         /* best-to-insert */
         prv_best = calc_Max(prv_M, prv_I);
         IMX(qx0, t_0) = prv_best + ISC(t_0, A);

         /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         prv_M = MMX(qx0, t_1) + TSC(t_1, M2D);
         prv_D = DMX(qx0, t_1) + TSC(t_1, D2D);
         /* best-to-delete */
         prv_best = calc_Max(prv_M, prv_D);
         DMX(qx0, t_0) = prv_best;

         /* UPDATE E STATE */
         prv_E = XMX(SP_E, q_0);
         prv_M = MMX(qx0, t_0) + sc_E;
         /* best-to-e-state */
         XMX(SP_E, q_0) = calc_Max( prv_E, prv_M );

         /* embed linear row into quadratic test matrix */
         #if DEBUG
         {
            MX_2D(cloud_MX, q_0, t_0) = 1.0;
            MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
            MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
            MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
         }
         #endif
      }

      /* UNROLLED FINAL LOOP ITERATION */
      t_0 = T;
      t_1 = t_0 - 1;

      /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
      /* best previous state transition (match takes the diag element of each prev state) */
      prv_M = MMX(qx1, t_1)  + TSC(t_1, M2M);
      prv_I = IMX(qx1, t_1)  + TSC(t_1, I2M);
      prv_D = DMX(qx1, t_1)  + TSC(t_1, D2M);
      prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M); /* from begin match state (new alignment) */
      /* best-to-match */
      prv_best = calc_Max(
                     calc_Max( prv_M, prv_I ),
                     calc_Max( prv_D, prv_B ) );
      MMX(qx0, t_0) = prv_best + MSC(t_0, A);

      /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
      /* previous states (match takes the left element of each state) */
      prv_M = MMX(qx0, t_1) + TSC(t_1, M2D);
      prv_D = DMX(qx0, t_1) + TSC(t_1, D2D);
      /* best-to-delete */
      prv_best = calc_Max( prv_M, prv_D );
      DMX(qx0, t_0) = prv_best;

      /* UPDATE E STATE */
      prv_E = XMX(SP_E, q_0);
      prv_M = MMX(qx0, t_0);
      prv_D = DMX(qx0, t_0);
      XMX(SP_E, q_0) = calc_Max( prv_E, 
                           calc_Max( prv_M, prv_D ) );

      /* SPECIAL STATES */
      /* J state */
      prv_J = XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP);     /* J->J */
      prv_E  = XMX(SP_E, q_0) + XSC(SP_E, SP_LOOP);     /* E->J is E's "loop" */
      XMX(SP_J, q_0) = calc_Max( prv_J, prv_E );

      /* C state */
      prv_C = XMX(SP_C, q_1) + XSC(SP_C, SP_LOOP);
      prv_E  = XMX(SP_E, q_0) + XSC(SP_E, SP_MOVE);
      XMX(SP_C, q_0) = calc_Max( prv_C, prv_E );

      /* N state */
      prv_N = XMX(SP_N, q_1) + XSC(SP_N, SP_LOOP);
      XMX(SP_N, q_0) = prv_N;

      /* B state */
      prv_N = XMX(SP_N, q_0) + XSC(SP_N, SP_MOVE);     /* N->B is N's move */
      prv_J = XMX(SP_J, q_0) + XSC(SP_J, SP_MOVE);     /* J->B is J's move */
      XMX(SP_B, q_0) = calc_Max( prv_N, prv_J );

      /* embed linear row into quadratic test matrix */
      #if DEBUG
      {
         MX_2D(cloud_MX, q_0, t_0) = 1.0;
         MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
         MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
         MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
      }
      #endif
   }

   /* T state (stores final state score) */
   sc_best = XMX(SP_C, Q) + XSC(SP_C, SP_MOVE);
   *sc_final = sc_best;

   return STATUS_SUCCESS;
}

/*
 *  FUNCTION:  run_Viterbi_Reverse_Quad()
 *  SYNOPSIS:  Run Backward version of Viterbi Algorithm (Seq-to-Profile, general unoptimized)
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Viterbi_Reverse_Quad(    const SEQUENCE*    query,        /* query sequence */
                                 const HMM_PROFILE* target,       /* target hmm model */
                                 const int          Q,            /* query length */
                                 const int          T,            /* target length */
                                 MATRIX_3D*         st_MX,        /* normal matrix */
                                 MATRIX_2D*         sp_MX,        /* special matrix */
                                 float*             sc_final )    /* final max score */
{
   /* vars for accessing query/target data structs */
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   char*    seq;                             /* alias for getting seq */
   int      N;                               /* length of edgebound list */
   bool     is_local;                        /* whether using local or global alignments */

   /* vars for indexing into data matrices by row-col */
   int      b, d, i, j, k;                   /* antidiagonal, row, column indices */
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* mod mapping of column index into data matrix (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */

   /* vars for indexing into data matrices by anti-diag */
   int      d_0, d_1, d_2;                   /* real index of current and previous antidiagonals */
   int      dx0, dx1, dx2;                   /* mod mapping of antidiagonal index into data matrix */
   int      k_0, k_1;                        /* offset into antidiagonal */
   int      d_st, d_end, d_cnt;              /* starting and ending diagonal indices */
   int      dim_T, dim_Q, dim_TOT;           /* dimensions of submatrix being searched */
   int      dim_min, dim_max;                /* diagonal index where num cells reaches highest point and diminishing point */ 
   int      num_cells;                       /* number of cells in current diagonal */

   /* vars for indexing into edgebound lists */
   BOUND*   bnd;                             /* current bound */
   BOUND    bnd_new;                         /* for adding new bound to edgebound list */
   int      id;                              /* id in edgebound list (row/diag) */
   int      r_0;                             /* current index in edgebound list */
   int      r_0b, r_0e;                      /* begin and end indices for current row in edgebound list */
   int      r_1b, r_1e;                      /* begin and end indices for current row in edgebound list */
   int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   bool     rb_T;                            /* checks if edge touches right bound of matrix */

   /* vars for recurrance scores */
   float    prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;              /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C; /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;            /* previous loop and move for special states */
   float    prv_sum, prv_best;             /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_2D*  cloud_MX3;
   MATRIX_3D*  test_MX;
   MATRIX_3D*  test_MX3;
   EDGEBOUNDS* test_edg;
   int         num_writes;
   int         num_clears;

   /* initialize debugging matrix */
   #if DEBUG
   {
      cloud_MX    = debugger->cloud_MX;
      cloud_MX3   = debugger->cloud_MX3;
      test_MX     = debugger->test_MX;
      test_MX3    = debugger->test_MX3;
      test_edg    = debugger->test_edg;

      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
      MATRIX_2D_Fill( cloud_MX3, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
      MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
      MATRIX_3D_Fill( test_MX3, -INF );
      EDGEBOUNDS_Reuse( test_edg, Q, T );

      num_writes = 0;
      num_clears = 0;
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* query sequence */
   seq         = query->seq;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   /* Initialize the Q row. */
   q_0 = Q;
   qx0 = q_0;
   t_0 = T;

   XMX(SP_J, q_0) = XMX(SP_B, q_0) = XMX(SP_N, q_0) = -INF;
   XMX(SP_C, q_0) = XSC(SP_C, SP_MOVE);
   XMX(SP_E, q_0) = XMX(SP_C, q_0) + XSC(SP_E, SP_MOVE);

   MMX(qx0, t_0) = DMX(qx0, t_0) = XMX(SP_E, q_0);
   IMX(qx0, t_0) = -INF;

   #if DEBUG 
   {
      MX_2D(cloud_MX, q_0, t_0) = 1.0;
      MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
      MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
      MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
   }
   #endif

   for (t_0 = T-1; t_0 >= 1; t_0--)
   {
      t_1 = t_0 + 1;

      prv_E = XMX(SP_E, Q) + sc_E;
      prv_D = DMX(qx0, t_1)  + TSC(t_0, M2D);
      MMX(qx0, t_0) = calc_Max( prv_E, prv_D );

      prv_E = XMX(SP_E, Q) + sc_E;
      prv_D = DMX(qx0, t_1)  + TSC(t_0, D2D);
      DMX(qx0, t_0) = calc_Max( prv_E, prv_D );

      IMX(qx0, t_0) = -INF;

      #if DEBUG 
      {
         MX_2D(cloud_MX, q_0, t_0) = 1.0;
         MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
         MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
         MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
      }
      #endif
   }

   /* MAIN RECURSION */
   /* FOR every position in QUERY seq */
   for (q_0 = Q-1; q_0 >= 1; q_0--)
   {
      q_1 = q_0 + 1;
      qx0 = q_0;
      qx1 = q_1;

      /* Get next sequence character */
      a = seq[q_0];
      A = AA_REV[a];

      /* SPECIAL STATES */

      /* B STATE -> MATCH */
      /* NOTE: When j = 0, MMX and MSC do not match HMMER p7_GBackward() implementation.   */
      XMX(SP_B, q_0) = MMX(qx1, 1) + TSC(0, B2M) + MSC(1, A);
      for (t_0 = 2; t_0 <= T; t_0++) {
         t_1 = t_0 - 1;
         prv_sum = XMX(SP_B, q_0);
         prv_M = MMX(qx1, t_0) + TSC(t_1, B2M) + MSC(t_0, A);
         XMX(SP_B, q_0) = calc_Max( prv_sum, prv_M);
      }

      prv_J = XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP);
      prv_B  = XMX(SP_B, q_0) + XSC(SP_J, SP_MOVE);
      XMX(SP_J, q_0) = calc_Max( prv_J, prv_B );

      prv_C = XMX(SP_C, q_1) + XSC(SP_C, SP_LOOP);
      XMX(SP_C, q_0) = prv_C;

      prv_J = XMX(SP_J, q_0) + XSC(SP_E, SP_LOOP);
      prv_C = XMX(SP_C, q_0) + XSC(SP_E, SP_MOVE);
      XMX(SP_E, q_0) = calc_Max( prv_J, prv_C );

      prv_N = XMX(SP_N, q_1) + XSC(SP_N, SP_LOOP);
      prv_B  = XMX(SP_B, q_0) + XSC(SP_N, SP_MOVE);
      XMX(SP_N, q_0) = calc_Max( prv_N, prv_B );

      t_0 = T;
      MMX(qx0, T) = DMX(qx0, T) = XMX(SP_E, q_0);
      IMX(qx0, T) = -INF;

      #if DEBUG 
      {
         MX_2D(cloud_MX, q_0, t_0) = 1.0;
         MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
         MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
         MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
      }
      #endif

      /* FOR every position in TARGET profile */
      for (t_0 = T-1; t_0 >= 1; t_0--)
      {
         t_1 = t_0 + 1;

         /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
         prv_M = MMX(qx1, t_1) + TSC(t_0, M2M) + MSC(t_1, A);
         prv_I = IMX(qx1, t_0) + TSC(t_0, M2I) + ISC(t_1, A);
         prv_D = DMX(qx0, t_1) + TSC(t_0, M2D);
         prv_E = XMX(SP_E, q_0) + sc_E;     /* from end match state (new alignment) */
         /* best-to-match */
         prv_sum = calc_Max( 
                        calc_Max( prv_M, prv_I ),
                        calc_Max( prv_E, prv_D ) );
         MMX(qx0, t_0) = prv_sum;

         /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
         prv_M = MMX(qx1, t_1) + TSC(t_0, I2M) + MSC(t_1, A);
         prv_I = IMX(qx1, t_0) + TSC(t_0, I2I) + ISC(t_0, A);
         /* best-to-insert */
         prv_sum = calc_Max( prv_M, prv_I );
         IMX(qx0, t_0) = prv_sum;

         /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
         prv_M = MMX(qx1, t_1) + TSC(t_0, D2M) + MSC(t_1, A);
         prv_D = DMX(qx0, t_1) + TSC(t_0, D2D);
         prv_E = XMX(SP_E, q_0) + sc_E;
         /* best-to-delete */
         prv_sum = calc_Max( prv_M, 
                        calc_Max( prv_D, prv_E ) );
         DMX(qx0, t_0) = prv_sum;

         #if DEBUG 
         {
            MX_2D(cloud_MX, q_0, t_0) = 1.0;
            MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
            MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
            MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
         }
         #endif
      }
   }

   /* FINAL ROW (i = 0) */
   /* At q_0 = 0, only N,B states are reachable. */
   q_0 = 0;
   q_1 = q_0 + 1;
   qx0 = q_0;
   qx1 = q_1;

   t_0 = 0;
   t_1 = t_0 + 1;

   a = seq[q_0];
   A = AA_REV[a];

   /* SPECIAL STATES */

   /* B STATE -> MATCH */
   prv_M = MMX(q_1, t_1) + TSC(0, B2M) + MSC(1, A);
   XMX(SP_B, q_0) = prv_M;
   for (t_0 = 2; t_0 <= T; t_0++) {
      t_1 = t_0-1;
      prv_sum = XMX(SP_B, q_0);
      prv_M = MMX(q_1, t_0) + TSC(t_1, B2M) + MSC(t_0, A);
      XMX(SP_B, q_0) = calc_Max( prv_sum, prv_M );
   }

   XMX(SP_J, q_0) = -INF;
   XMX(SP_C, q_0) = -INF;
   XMX(SP_E, q_0) = -INF;

   prv_N = XMX(SP_N, q_1) + XSC(SP_N,SP_LOOP);
   prv_B  = XMX(SP_B, q_0) + XSC(SP_N,SP_MOVE);
   XMX(SP_N, q_0) = calc_Max( prv_N, prv_B );

   for (t_0 = T; t_0 >= 1; t_0--) {
      MMX(qx0, t_0) = IMX(qx0, t_0) = DMX(qx0, t_0) = -INF;
   }

   #if DEBUG 
   {
      MX_2D(cloud_MX, q_0, t_0) = 1.0;
      MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(qx0, t_0);
      MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(qx0, t_0);
      MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(qx0, t_0);
   }
   #endif

   sc_best = XMX(SP_N, 0);
   *sc_final = sc_best;

   /* flag matrices that they contain dirty values (not -INF) */
   st_MX->clean = false;
   sp_MX->clean  = false;

   return STATUS_SUCCESS;
}
