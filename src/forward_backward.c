/*******************************************************************************
 *  @file forward_backward.c
 *  @brief The Forward-Backward Algorithm for Sequence Alignment Search.
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

/*  
 *  FUNCTION: forward_Run()
 *  SYNOPSIS: Perform Forward-Backward Algorithm (general, unoptimized)
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <res>       Results Data
 *
 *  RETURN: 
 */
void fwdbck_Run (const SEQ* query, 
                  const HMM_PROFILE* target,
                  int Q, int T, 
                  float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                  float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                  RESULTS* res)
{
   init_Logsum();

   /* allocate DP matrix */

   /* NORMAL STATE (MATCH, INSERT, DELETE) MATRIX */
   float st_MX2[ NUM_NORMAL_STATES * (Q+1) * (T+1) ]; 
   /* SPECIAL STATES MATRIX */
   float sp_MX2[ NUM_SPECIAL_STATES * (Q+1) ];

   /* NORMAL STATE (MATCH, INSERT, DELETE) MATRIX */
   float st_MX3[ NUM_NORMAL_STATES * (Q+1) * (T+1) ]; 
   /* SPECIAL STATES MATRIX */
   float sp_MX3[ NUM_SPECIAL_STATES * (Q+1) ];

   for (int i = 0; i <= Q; i++)
   {
      for (int j = 0; j <= T; j++)
      {
         ST_MX(st_MX,MAT_ST,i,j)  = ST_MX(st_MX,INS_ST,i,j)  = ST_MX(st_MX,DEL_ST,i,j)  = -INF;
         ST_MX(st_MX2,MAT_ST,i,j) = ST_MX(st_MX2,INS_ST,i,j) = ST_MX(st_MX2,DEL_ST,i,j) = -INF;
         ST_MX(st_MX3,MAT_ST,i,j) = ST_MX(st_MX3,INS_ST,i,j) = ST_MX(st_MX3,DEL_ST,i,j) = -INF;
      }

      for (int j = 0; j < NUM_SPECIAL_STATES; j++)
      {
         SP_MX(sp_MX,j,i) = -INF;
         SP_MX(sp_MX2,j,i) = -INF;
         SP_MX(sp_MX3,j,i) = -INF;
      }
   }

	forward_Run(query, target, Q, T, st_MX, sp_MX, res);
   backward_Run(query, target, Q, T, st_MX2, sp_MX2, res);

   int i,j;
   for (i = 0; i < Q; i++)
   {
      for (j = 0; j < T; j++)
      {
         ST_MX(st_MX3,MAT_ST,i,j) = ST_MX(st_MX,MAT_ST,i,j) + ST_MX(st_MX2,MAT_ST,i,j);
         ST_MX(st_MX3,INS_ST,i,j) = ST_MX(st_MX,INS_ST,i,j) + ST_MX(st_MX2,INS_ST,i,j);
         ST_MX(st_MX3,DEL_ST,i,j) = ST_MX(st_MX,DEL_ST,i,j) + ST_MX(st_MX2,DEL_ST,i,j);
      }
   }
}

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
float forward_Run (const SEQ* query, 
                   const HMM_PROFILE* target, 
                   int Q, int T, 
                   float st_MX[ NUM_NORMAL_STATES * (Q+1) * 3 ], 
                   float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                   RESULTS* res)
{
   init_Logsum();

   char   a;           /* store current character in sequence */
   int    A;           /* store int value of character */
   int    i,j,k = 0;   /* row, column indices */
   char   *seq = query->seq; /* alias for getting seq */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_sum;
   float  sc, sc_1, sc_2;
   float  sc_max, sc_best;
   float  sc1, sc2, sc3, sc4;
   COORDS tr_end; /* ending match state of optimal alignment (for traceback) */

   /* local or global (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* --------------------------------------------------------------------------------- */

   /* clear all pre-existing data from matrix */
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
         sc1 = prev_mat = MMX(i-1,j-1)  + TSC(j-1,M2M);
         sc2 = prev_ins = IMX(i-1,j-1)  + TSC(j-1,I2M);
         sc3 = prev_del = DMX(i-1,j-1)  + TSC(j-1,D2M);
         sc4 = prev_beg = XMX(SP_B,i-1) + TSC(j-1,B2M); /* from begin match state (new alignment) */

         /* best-to-match */
         prev_sum = calc_Logsum( 
                        calc_Logsum( prev_mat, prev_ins ),
                        calc_Logsum( prev_beg, prev_del )
                     );
         MMX(i,j) = prev_sum + MSC(j,A);

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
      }

      /* UNROLLED FINAL LOOP ITERATION */
      j = T; 

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

      /* FIND SUM OF PATHS TO INSERT STATE (unrolled) */
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

   return sc_best;
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
float backward_Run (const SEQ* query, 
                    const HMM_PROFILE* target, 
                    int Q, int T, 
                    float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                    float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                    RESULTS* res)
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

   /* local or global (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

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
      /* Get next sequence character */
      a = seq[i];
      A = AA_REV[a];

      /* SPECIAL STATES */
      j = 1; int x_0 = i; 
      XMX(SP_B,i) = MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A);
      // printf("B(%d)\t (%d)%f\t", x_0, j, XMX(SP_B, x_0) );

      /* B -> MATCH */
      for (j = 2; j <= T; j++) {
         XMX(SP_B,i) = calc_Logsum( XMX(SP_B,i),
                                    MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
         // printf("B(%d,%d): B=%f, M=%f, TSC=%f, MSC=%f\n", x_0, j, XMX(SP_B, x_0), MMX(i+1, j), TSC(j-1, B2M), MSC(j, A) );
      }
      // printf("\n");

      XMX(SP_J,i) = calc_Logsum( XMX(SP_J,i+1) + XSC(SP_J,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_J,SP_MOVE) );

      XMX(SP_C,i) = XMX(SP_C,i+1) + XSC(SP_C,SP_LOOP);

      XMX(SP_E,i) = calc_Logsum( XMX(SP_J,i) + XSC(SP_E,SP_LOOP),
                                 XMX(SP_C,i) + XSC(SP_E,SP_MOVE) );

      XMX(SP_N,i) = calc_Logsum( XMX(SP_N,i+1) + XSC(SP_N,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_N,SP_MOVE) );

      MMX(i,T) = DMX(i,T) = XMX(SP_E,i);
      IMX(i,T) = -INF;

      x_0 = i;
      // printf("x_0: %d -> B: %f, J: %f, C: %f, E: %f, N: %f, MSC(%d): %f \n", x_0, XMX(SP_B, x_0), XMX(SP_J, x_0), XMX(SP_C, x_0),  XMX(SP_E, x_0), XMX(SP_N, x_0), A, MSC(j, A));

      /* FOR every position in TARGET profile */
      for (j = T-1; j >= 1; j--)
      {
         sc_M = MSC(j+1,A);
         sc_I = ISC(j,A);
         // printf("(%d,%d): A=%d, MSC=%f, ISC=%f\n", x_0, j, A, sc_M, sc_I);

         /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
         sc1 = prev_mat = MMX(i+1,j+1) + TSC(j,M2M) + sc_M;
         sc2 = prev_ins = IMX(i+1,j)   + TSC(j,M2I) + sc_I;
         sc3 = prev_del = DMX(i,j+1)   + TSC(j,M2D);
         sc4 = prev_end = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
         /* best-to-match */
         prev_sum = calc_Logsum( 
                           calc_Logsum( prev_mat, prev_ins ),
                           calc_Logsum( prev_end, prev_del ) );
         MMX(i,j) = prev_sum;

         // printf("M (%d,%d)  \tMM: %.5f\t IM: %.5f\t DM: %.5f\t BM: %.5f\t MSC: %.5f\t ISC: %.5f\t M2M: %.5f\t M2I: %.5f\t M2D: %.5f\t DMX: %.5f\n", 
         //    i, j, sc1, sc2, sc3, sc4, sc_M, sc_I, TSC(j,M2M), TSC(j,M2I), TSC(j,M2D), DMX(i,j+1) );

         /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
         sc1 = prev_mat = MMX(i+1,j+1) + TSC(j,I2M) + sc_M;
         sc2 = prev_ins = IMX(i+1,j)   + TSC(j,I2I) + sc_I;
         /* best-to-insert */
         prev_sum = calc_Logsum( prev_mat, prev_ins );
         IMX(i,j) = prev_sum;

         // printf("I (%d,%d)  \tIM: %.5f\t II: %.5f\t I2M: %.5f\t I2I: %.5f\n", 
            // i, j, sc1, sc2, TSC(j,I2M), TSC(j,I2I) );

         /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
         sc1 = prev_mat = MMX(i+1,j+1) + TSC(j,D2M) + sc_M;
         sc2 = prev_del = DMX(i,  j+1) + TSC(j,D2D);
         sc3 = prev_end = XMX(SP_E,i)  + sc_E;
         /* best-to-delete */
         prev_sum = calc_Logsum( 
                           prev_mat, 
                           calc_Logsum( prev_del, prev_end ) );
         DMX(i,j) = prev_sum;

         int r_0 = i;
         // printf("(%d,%d): M=%f, I=%f, D=%f\n", r_0, j, MMX(r_0, j), IMX(r_0,j), DMX(r_0,j) );
      }
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

   // printf("i\tN\tJ\tC\tE\tB\n");
   // for (i = 0; i < Q; ++i) {
   //    printf("%d\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n", i, XMX(SP_N,i), XMX(SP_J,i), XMX(SP_C,i), XMX(SP_E,i), XMX(SP_B,i));
   // }

   sc_best = XMX(SP_N,0);

   return sc_best;
}








