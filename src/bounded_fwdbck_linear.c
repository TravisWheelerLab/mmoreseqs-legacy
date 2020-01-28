/*******************************************************************************
 *  @file bounded_fwdbck_linear.c
 *  @brief Cloud Forward-Backward Algorithm (Linear Space Implementation)
 *
 *  @synopsis
 *       NOTE: HOW TO CONVERT row-coords to diag-coords
 *       MMX3(i-1,j-1) => MMX3(, d_2)
 *       MMX3(i,  j-1) => MMX3(, d_1)
 *       MMX3(i,  j  ) => MMX3(, d_1)
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
#include "forward_backward.h"
#include "bounded_fwdbck_linear.h"

// macros
// #define getName(var) #var
// #define SCALE_FACTOR 1000

// macro functions
// NOTE: wrap all macro vars in parens!!
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))


/*  
 *  FUNCTION: forward_bounded_Run()
 *  SYNOPSIS: Perform Edge-Bounded Forward part of Cloud Search Algorithm.
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
float forward_bounded_Run3(const SEQ* query, 
                           const HMM_PROFILE* target,
                           int Q, int T, 
                           float* st_MX3,
                           float* st_MX,
                           float* sp_MX, 
                           EDGEBOUNDS* edg,
                           float *sc_final )
{
   char   a;                              /* store current character in sequence */
   int    A;                              /* store int value of character */
   int    i,j,k = 0;                      /* row, column indices */
   char   *seq = query->seq;              /* alias for getting seq */
   int    N = edg->N;                     /* length of edgebound list */

   int    x, y1, y2;                      /* row, leftcol and rightcol bounds in row */
   int    x_0, row_cur, x_1, row_prv;     /* real index of current and previous rows */
   int    r_0, r_1;                       /* row offset -> r_0: row_cur % 2, r_1: row_prv % 2 */
   int    r_0b, r_0e, r_1b, r_1e;         /* begin and end indices for row in edgebound list */
   bool   y2_re;

   float  prev_mat, prev_del, prev_ins;   /* temp placeholder sums */
   float  prev_beg, prev_end, prev_sum;   /* temp placeholder sums */
   float  sc, sc1, sc2, sc3, sc4;         /* temp placeholder sums (testing) */
   float  sc_best;                        /* alignment score (return value) */
   float  sc_M, sc_I, sc_D;

   int    d_mod3, d_0, d_1, d_2; /* d (mod 3) for assigning prev array ptrs */

   /* local or global? (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* DEBUG OUTPUT */
   FILE *tfp = fopen("test_output/fwd_bounded.linear.txt", "w+");

   /* --------------------------------------------------------------------------------- */

   /* clear all pre-existing data from matrix */
   dp_matrix_Clear3(Q, T, st_MX3, sp_MX);

   /* Clear top-row (0th row) */
   x_0 = 0;
   r_0 = x_0 % 2;

   /* Initialize special states (?) */
   XMX(SP_N, x_0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B, x_0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E, x_0) = XMX(SP_C, x_0) = XMX(SP_J, x_0) = -INF;          /* need seq to get here (?)  */

   /* Initialize zero row (top-edge) */
   for (j = 0; j < T; j++)
      MMX3(r_0, j) = IMX3(r_0, j) = DMX3(r_0, j) = -INF;              /* need seq to get here (?)  */

   /* pass over top-row (i=0) edgebounds from list */
   row_cur = r_0 = 0;   /* current row in matrix */
   k = 0;               /* current index in edgebounds */
   r_0b = 0;
   while ( k < N && edg->bounds[k].diag == row_cur ) {
      k++;
   }
   r_0e = k;

   /* init look back 1 (r_1) */
   x_1 = 0;
   r_1 = r_0;
   r_1b = r_0b;
   r_1e = r_0e;

   /* MAIN RECURSION */
   /* FOR every position in QUERY sequence (row in matrix) */
   for (x_0 = 1; x_0 <= Q; x_0++)
   {
      fprintf(tfp, "(i=%d): ", x_0);
      /* convert quadratic space row index to linear space row index (ex % 2) */
      row_cur = x_0;
      x_1 = x_0 - 1;
      r_0 = x_0 % 2;        /* for use in linear space alg (mod-mapping) */
      r_1 = x_1 % 2;        /* for use in linear space alg (mod-mapping) */

      /* add every edgebound from current row */
      r_0b = k;
      // fprintf("k: %d => row_cur: %d, k.diag: %d\n", k, row_cur, edg->bounds[k].diag);
      while ( k < N && edg->bounds[k].diag == row_cur ) {
         k++;
      }
      r_0e = k;

      /* Get next sequence character */
      a = seq[x_1];  /* off-by-one */
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX3(r_0, 0) = IMX3(r_0, 0) = DMX3(r_0, 0) = -INF;
      XMX(SP_E, x_0) = -INF;

      /* FOR every EDGEBOUND in current ROW */
      for (i = r_0b; i < r_0e; i++)
      {
         /* in this context, "diag" represents the "row" */
         x = edg->bounds[i].diag;               /* NOTE: this is always the same as cur_row, x_0 */
         y1 = max(1, edg->bounds[i].lb);        /* can't overflow the left edge */
         y2 = edg->bounds[i].rb;
         y2_re = (y2 > T);                      /* check if cloud touches right edge */
         y2 = min(y2, T);                       /* can't overflow the right edge */
         
         /* MAIN RECURSION */
         /* FOR every position in TARGET profile */
         for (j = y1; j < y2; j++)
         {
            fprintf(tfp, "%d...", j);
            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            sc1 = prev_mat = MMX3(r_1,j-1)  + TSC(j-1,M2M);
            sc2 = prev_ins = IMX3(r_1,j-1)  + TSC(j-1,I2M);
            sc3 = prev_del = DMX3(r_1,j-1)  + TSC(j-1,D2M);
            sc4 = prev_beg = XMX(SP_B, x_1) + TSC(j-1,B2M); /* from begin match state (new alignment) */
            /* best-to-match */
            prev_sum = calc_Logsum( 
                           calc_Logsum( prev_mat, prev_ins ),
                           calc_Logsum( prev_del, prev_beg ) );
            MMX3(r_0,j) = prev_sum + MSC(j,A);
            // fprintf(tfp, "MMX3(%d,%d): m=%f, i=%f, d=%f, b=%f\n", x_0, j, sc1, sc2, sc3, sc4 );

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX3(r_1,j) + TSC(j,M2I);
            prev_ins = IMX3(r_1,j) + TSC(j,I2I);
            /* best-to-insert */
            prev_sum = calc_Logsum( prev_mat, prev_ins );
            IMX3(r_0,j) = prev_sum + ISC(j,A);

            /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX3(r_0,j-1) + TSC(j-1,M2D);
            prev_del = DMX3(r_0,j-1) + TSC(j-1,D2D);
            /* best-to-delete */
            prev_sum = calc_Logsum(prev_mat, prev_del);
            DMX3(r_0,j) = prev_sum;

            /* UPDATE E STATE */
            prev_mat = MMX3(r_0, j) + sc_E;
            prev_del = DMX3(r_0, j) + sc_E;
            XMX(SP_E, x_0) = calc_Logsum( 
                                    calc_Logsum( prev_mat, prev_del ),
                                    XMX(SP_E, x_0) );

            fprintf(tfp, "PRE (%d,%d) -> E: %f, Esc: %f, MMX: %f, DMX: %f, MSC: %f \n", x_0, j, XMX(SP_E, x_0), sc_E, MMX3(r_0, j), DMX3(r_0, j), MSC(j,A) );
         }

         /* UNROLLED FINAL LOOP ITERATION */
         if ( y2_re ) 
         {
            fprintf(tfp, "%d!...", j);
            j = T; 

            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            prev_mat = MMX3(r_1,j-1)  + TSC(j-1,M2M);
            prev_ins = IMX3(r_1,j-1)  + TSC(j-1,I2M);
            prev_del = DMX3(r_1,j-1)  + TSC(j-1,D2M);
            prev_beg = XMX(SP_B,r_1) + TSC(j-1,B2M);    /* from begin match state (new alignment) */
            /* best-to-match */
            prev_sum = calc_Logsum( 
                              calc_Logsum( prev_mat, prev_ins ),
                              calc_Logsum( prev_del, prev_beg )
                        );
            MMX3(r_0,j) = prev_sum + MSC(j,A);

            /* FIND SUM OF PATHS TO INSERT STATE */
            IMX3(r_0,j) = -INF;

            /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX3(r_0,j-1) + TSC(j-1,M2D);
            prev_del = DMX3(r_0,j-1) + TSC(j-1,D2D);
            /* best-to-delete */
            prev_sum = calc_Logsum( prev_mat, prev_del );
            DMX3(r_0,j) = prev_sum;

            /* UPDATE E STATE */
            prev_mat = MMX3(r_0, j);
            prev_del = DMX3(r_0, j);
            XMX(SP_E, x_0) = calc_Logsum( 
                                 calc_Logsum( prev_mat, prev_del ),
                                 XMX(SP_E, x_0) );
         }
      }
      fprintf(tfp, "...done\n");

      /* ONCE ROW IS COMPLETED, UPDATE SPECIAL STATES */
      {
         /* SPECIAL STATES */
         /* J state */
         sc1 = XMX(SP_J, x_1) + XSC(SP_J, SP_LOOP);   /* J->J */
         sc2 = XMX(SP_E, x_0) + XSC(SP_E, SP_LOOP);   /* E->J is E's "loop" */
         XMX(SP_J, x_0) = calc_Logsum( sc1, sc2 );

         /* C state */
         sc1 = XMX(SP_C, x_1) + XSC(SP_C, SP_LOOP);
         sc2 = XMX(SP_E, x_0) + XSC(SP_E, SP_MOVE);
         XMX(SP_C, x_0) = calc_Logsum( sc1, sc2 );
         // printf("x_0: %d -> C(-1): %f, Cloop: %f, E(0): %f, Emove: %f\n", x_0, XMX(SP_C, x_1), XSC(SP_C, SP_LOOP), XMX(SP_E, x_0), XSC(SP_E, SP_MOVE) );

         /* N state */
         XMX(SP_N, x_0) = XMX(SP_N, x_1) + XSC(SP_N, SP_LOOP);

         /* B state */
         sc1 = XMX(SP_N, x_0) + XSC(SP_N, SP_MOVE);       /* N->B is N's move */
         sc2 = XMX(SP_J, x_0) + XSC(SP_J, SP_MOVE);       /* J->B is J's move */
         XMX(SP_B, x_0) = calc_Logsum( sc1, sc2 );

         // printf("x_0: %d -> J: %f, C: %f, N: %f, B: %f\n", x_0, XMX(SP_J, x_0), XMX(SP_C, x_0), XMX(SP_N, x_0), XMX(SP_N, x_0));
      }

      /* SCRUB PREVIOUS (LOOK-BACK-2) ROW VALUES */
      for (i = r_1b; i < r_1e; i++) 
      {
         /* in this context, "diag" represents the "row" */
         // x  = edg->bounds[i].diag;
         y1 = max(0, edg->bounds[i].lb);     /* can't overflow the left edge */
         y2 = min(edg->bounds[i].rb, T);     /* can't overflow the right edge */

         for (j = y1; j <= y2; j++) {
            MMX3(r_1, j) = IMX3(r_1, j) = DMX3(r_1, j) = -INF;
         }
      }

      // /* NAIVE SCRUB (removes all values) - FOR DEBUGGING */
      // for (j = 0; j <= T; j++) {
      //    MMX3(r_1, j) = IMX3(r_1, j) = DMX3(r_1, j) = -INF;
      // }

      /* EMBED LINEAR MX into QUAD MX - FOR DEBUGGING */
      for (j = 0; j <= T; j++) {
         MMX(x_0,j) = MMX3(r_0,j);
         IMX(x_0,j) = IMX3(r_0,j);
         DMX(x_0,j) = DMX3(r_0,j);
      }

      /* set curr rows to prv rows */
      row_prv = x_1 = row_cur;
      r_1 = r_0;
      r_1b = r_0b;
      r_1e = r_0e;
   }

   x_0 = Q;
   x_1 = Q - 1;
   r_0 = x_0 % 2;
   r_1 = x_1 % 2;

   /* T state */
   sc_best = XMX(SP_C, x_0) + XSC(SP_C,SP_MOVE);
   sc_final[0] = sc_best;

   fclose(tfp);
   return sc_best;
}


/*  
 *  FUNCTION: backward_bounded_Run3()
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
float backward_bounded_Run3(const SEQ* query, 
                           const HMM_PROFILE* target,
                           int Q, int T, 
                           float* st_MX3,
                           float* st_MX,
                           float* sp_MX, 
                           EDGEBOUNDS* edg,
                           float *sc_final )
{
   char   a;                              /* store current character in sequence */
   int    A;                              /* store int value of character */
   int    i,j,k = 0;                      /* row, column indices */
   char   *seq = query->seq;              /* alias for getting seq */
   int    N = edg->N;                     /* length of edgebound list */

   int    x, y1, y2;                    /* row, leftcol and rightcol bounds in row */
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

   /* DEBUG OUTPUT */
   FILE *tfp = fopen("test_output/bck_bounded.linear.txt", "w+");

   /* --------------------------------------------------------------------------------- */

   /* clear all pre-existing data from matrix */
   dp_matrix_Clear3(Q, T, st_MX3, sp_MX);

   /* Initialize the Q row. */
   row_cur = x_0 = Q;
   r_1 = r_0 + 1;
   r_0 = x_0 % 2;     /* for use in linear space alg (mod-mapping) */
   r_1 = x_1 % 2;     /* for use in linear space alg (mod-mapping) */

   XMX(SP_J, x_0) = XMX(SP_B, x_0) = XMX(SP_N, x_0) = -INF;
   XMX(SP_C, x_0) = XSC(SP_C,SP_MOVE);
   XMX(SP_E, x_0) = XMX(SP_C, x_0) + XSC(SP_E,SP_MOVE);

   MMX3(r_0, T) = DMX3(r_0, T) = XMX(SP_E, r_0);
   IMX3(r_0, T) = -INF;

   for (j = T-1; j >= 1; j--)
   {
      sc1 = XMX(SP_E, x_0) + sc_E;
      sc2 = DMX3(r_0, j+1)  + TSC(j, M2D);
      MMX3(r_0, j) = calc_Logsum( XMX(SP_E, x_0) + sc_E, 
                              DMX3(r_0, j+1)  + TSC(j, M2D) );

      sc1 = XMX(SP_E, x_0) + sc_E;
      sc2 = DMX3(r_0, j+1)  + TSC(j, D2D);
      DMX3(r_0, j) = calc_Logsum( XMX(SP_E, x_0) + sc_E,
                              DMX3(r_0, j+1)  + TSC(j, D2D) );

      IMX3(r_0,j) = -INF;
   }

   /* pass over (Q) bottom-row edgebounds from list */
   k = N-1;
   r_0b = N-1;
   while ( k >= 0 && edg->bounds[k].diag == row_cur ) {
      k--;
   }
   r_0e = k;

   /* MAIN RECURSION */
   /* FOR every bound in EDGEBOUND */
   for (x_0 = Q-1; x_0 > 0; --x_0)
   {
      fprintf(tfp, "(i=%d): ", x_0);
      /* convert quadratic space row index to linear space row index (ex % 2) */
      row_cur = x_0;
      x_1 = x_0 + 1;
      r_0 = x_0 % 2;    /* for use in linear space alg (mod-mapping) */
      r_1 = x_1 % 2;    /* for use in linear space alg (mod-mapping) */

      /* add every edgebound from current row */
      r_0b = k;
      // printf("k: %d => row_cur: %d, k.diag: %d\n", k, row_cur, edg->bounds[k].diag);
      while ( k >= 0 && edg->bounds[k].diag >= row_cur ) {
         k--;
      }
      r_0e = k;

      /* Get next sequence character */
      a = seq[x_0];
      A = AA_REV[a];

      /* UPDATE SPECIAL STATES at the start of EACH ROW */
      {
         /* SPECIAL STATES */
         j = 1;
         XMX(SP_B, x_0) = MMX3(r_1, j) + TSC(j-1, B2M) + MSC(j, A);
         // printf("B(%d)\t (%d)%f\t", x_0, j, XMX(SP_B, x_0) );

         /* B -> MATCH */
         for (j = 2; j <= T; j++) {
            XMX(SP_B, x_0) = calc_Logsum( XMX(SP_B, x_0),
                                       MMX3(r_1, j) + TSC(j-1, B2M) + MSC(j, A) );
            // printf("B(%d,%d): B=%f, M=%f, TSC=%f, MSC=%f\n", x_0, j, XMX(SP_B, x_0), MMX3(r_1, j), TSC(j-1, B2M), MSC(j, A) );
         }
         // printf("\n");

         XMX(SP_J, x_0) = calc_Logsum( XMX(SP_J, x_1) + XSC(SP_J, SP_LOOP),
                                       XMX(SP_B, x_0) + XSC(SP_J, SP_MOVE) );

         XMX(SP_C, x_0) = XMX(SP_C, x_1) + XSC(SP_C, SP_LOOP);

         XMX(SP_E, x_0) = calc_Logsum( XMX(SP_J, x_0) + XSC(SP_E,SP_LOOP),
                                       XMX(SP_C, x_0) + XSC(SP_E,SP_MOVE) );

         XMX(SP_N, x_0) = calc_Logsum( XMX(SP_N, x_1) + XSC(SP_N,SP_LOOP),
                                       XMX(SP_B, x_0) + XSC(SP_N,SP_MOVE) );

         MMX3(r_0, T) = DMX3(r_0, T) = XMX(SP_E, x_0);
         IMX3(r_0, T) = -INF;

         // printf("x_0: %d -> B: %f, J: %f, C: %f, E: %f, N: %f, MSC(%d): %f \n", x_0, XMX(SP_B, x_0), XMX(SP_J, x_0), XMX(SP_C, x_0),  XMX(SP_E, x_0), XMX(SP_N, x_0), A, MSC(j, A) );
      }

      // fprintf(tfp, "x_0: %d (%d, %d) -> (%d, %d)\n", x_0, r_0b, r_0e, edg->bounds[r_0b], edg->bounds[r_0e-1]);

      /* FOR every EDGEBOUND in current ROW */
      for (i = r_0b; i > r_0e; i--)
      {
         /* in this context, "diag" represents the "row" */
         // x_0  = edg->bounds[i].diag;
         y1 = max(1, edg->bounds[i].lb);     /* can't overflow the left edge */
         y2 = min(edg->bounds[i].rb, T);     /* can't overflow the right edge */

         /* FOR every position in TARGET profile */
         for (j = y2-1; j >= y1; --j)
         {
            fprintf(tfp, "%d...", j);
            sc_M = MSC(j+1,A);
            sc_I = ISC(j+1,A);
            // printf("(%d,%d): A=%d, MSC=%f, ISC=%f\n", x_0, j, A, sc_M, sc_I);

            /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
            prev_mat = MMX3(r_1, j+1)  + TSC(j, M2M) + sc_M;
            prev_ins = IMX3(r_1, j)    + TSC(j, M2I) + sc_I;
            prev_del = DMX3(r_0, j+1)  + TSC(j, M2D);
            prev_end = XMX(SP_E, x_0) + sc_E;     /* from end match state (new alignment) */
            /* best-to-match */
            prev_sum = calc_Logsum( 
                              calc_Logsum( prev_mat, prev_ins ),
                              calc_Logsum( prev_end, prev_del ) );
            MMX3(r_0, j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            prev_mat = MMX3(r_1, j+1) + TSC(j, I2M) + sc_M;
            prev_ins = IMX3(r_1, j)   + TSC(j, I2I) + sc_I;
            /* best-to-insert */
            prev_sum = calc_Logsum( prev_mat, prev_ins );
            IMX3(r_0,j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            prev_mat = MMX3(r_1, j+1)  + TSC(j,D2M) + sc_M;
            prev_del = DMX3(r_0, j+1)  + TSC(j,D2D);
            prev_end = XMX(SP_E, x_0) + sc_E;
            /* best-to-delete */
            prev_sum = calc_Logsum( 
                              prev_mat,
                              calc_Logsum( prev_del, prev_end ) );
            DMX3(r_0,j) = prev_sum;

            // printf("(%d,%d): M=%f, I=%f, D=%f\n", r_0, j, MMX3(r_0, j), IMX3(r_0,j), DMX3(r_0,j) );
         }
      }

      /* SCRUB PREVIOUS ROW VALUES (backward) */
      for (i = r_1b; i > r_1e; i--) 
      {
         /* in this context, "diag" represents the "row" */
         // x  = edg->bounds[i].diag;
         y1 = max(0, edg->bounds[i].lb);     /* can't overflow the left edge */
         y2 = min(edg->bounds[i].rb, T);     /* can't overflow the right edge */

         for (j = y1-1; j >= y2; --j) {
            MMX3(r_1, j) = IMX3(r_1, j) = DMX3(r_1, j) = -INF;
         }
      }

      // /* NAIVE SCRUB (removes all values) - FOR DEBUGGING */
      // for (i = 0; i <= T; i++) {
      //    MMX3(r_1, j) = IMX3(r_1, j) = DMX3(r_1, j) = -INF;
      // }

      /* EMBED LINEAR MX into QUAD MX - FOR DEBUGGING */
      for (i = 0; i <= T; i++) {
         MMX(x_0,j) = MMX3(r_0,j);
         IMX(x_0,j) = IMX3(r_0,j);
         DMX(x_0,j) = DMX3(r_0,j);
      }

      /* SET CURRENT ROW TO PREVIOUS ROW */
      row_prv = x_1 = row_cur;
      r_1 = r_0;
      r_1b = r_0b;
      r_1e = r_0e;

      fprintf(tfp, "...done\n");
      fprintf(tfp, "N(0): %f\n", XMX(SP_N, 0));
   }

   /* FINAL ROW */
   row_cur = x_0 = i = 0;
   x_1 = x_0 + 1;
   r_0 = x_0 % 2;
   r_1 = x_1 % 2;

   /* FINAL i = 0 row */
   a = seq[1];
   A = AA_REV[a];

   j = 1;
   XMX(SP_B,x_0) = MMX3(r_1,j) + TSC(x_0,B2M) + MSC(j,A);

   for (j = 2; j >= T; j++) {
      XMX(SP_B,x_0) = calc_Logsum( 
                           XMX(SP_B,x_0),
                           MMX3(r_1,j) + TSC(j-1,B2M) + MSC(j,A) );
   }

   XMX(SP_J,x_0) = -INF;
   XMX(SP_C,x_0) = -INF;
   XMX(SP_E,x_0) = -INF;

   XMX(SP_N,x_0) = calc_Logsum( XMX(SP_N,x_1) + XSC(SP_N,SP_LOOP),
                                XMX(SP_B,x_0) + XSC(SP_N,SP_MOVE) );

   for (j = T; j >= 1; j--) {
      MMX3(x_0,j) = IMX3(x_0,j) = DMX3(x_0,j) = -INF;
   }

   fprintf(tfp, "a=%d, A=%d, MSC=%f\n", a, A, MSC(1,A) );

   /* EMBED LINEAR MX into QUAD MX - FOR DEBUGGING */
   for (j = 0; j <= T; j++) {
      MMX(x_0,j) = MMX3(r_0,j);
      IMX(x_0,j) = IMX3(r_0,j);
      DMX(x_0,j) = DMX3(r_0,j);
   }

   sc_best = XMX(SP_N,0);
   sc_final[0] = sc_best;
   fclose(tfp);
   return sc_best;
}

