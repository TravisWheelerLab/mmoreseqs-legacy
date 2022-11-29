/*******************************************************************************
 *  - FILE:  viterbi_linear.c
 *  - DESC:  The Viterbi Algorithm and Traceback for Sequence Alignment Search.
 *******************************************************************************/

/* imports */
#include <stdbool.h>
#include <math.h>

/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"

/* header */
#include "viterbi_linear.h"

int run_Viterbi_Linear(const SEQUENCE* query,
                       const HMM_PROFILE* target,
                       const int Q,
                       const int T,
                       MATRIX_3D* st_MX3,
                       MATRIX_2D* sp_MX,
                       float* sc_final) {
  /* vars for accessing query/target data structs */
  char a;        /* store current character in sequence */
  int A;         /* store int value of character */
  char* seq;     /* alias for getting seq */
  bool is_local; /* whether using local or global alignments */

  /* vars for indexing into data matrices by row-col */
  int q_0, q_1;      /* real index of current and previous rows (query) */
  int qx0, qx1;      /* mod mapping of column index into data matrix (query) */
  int t_0, t_1;      /* real index of current and previous columns (target) */


  /* vars for recurrance scores */
  float prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
  float prv_B, prv_E;           /* previous (B) begin and (E) end states */
  float prv_J, prv_N, prv_C;    /* previous (J) jump, (N) initial, and (C) terminal states */
  float prv_best;      /* temp subtotaling vars */
  float sc_best;                /* final best scores */
  float sc_E; /* match, insert, delete, end scores */


/* initialize debugging matrix */
#if DEBUG
  {
    /* debugger tools */
    FILE* dbfp;
    MATRIX_2D* cloud_MX;
    MATRIX_2D* cloud_MX3;
    MATRIX_3D* test_MX;
    MATRIX_3D* test_MX3;
    int num_writes;
    int num_clears;

    cloud_MX = debugger->cloud_MX;
    cloud_MX3 = debugger->cloud_MX3;
    test_MX = debugger->test_MX;
    test_MX3 = debugger->test_MX3;

    MATRIX_2D_Reuse(cloud_MX, Q + 1, T + 1);
    MATRIX_2D_Fill(cloud_MX, 0);
    MATRIX_2D_Reuse(cloud_MX3, 3, (Q + 1) + (T + 1));
    MATRIX_2D_Fill(cloud_MX3, 0);
    MATRIX_3D_Reuse(test_MX, NUM_NORMAL_STATES, Q + 1, T + 1);
    MATRIX_3D_Fill(test_MX, -INF);
    MATRIX_3D_Reuse(test_MX3, NUM_NORMAL_STATES, 3, (Q + 1) + (T + 1));
    MATRIX_3D_Fill(test_MX3, -INF);

    num_writes = 0;
    num_clears = 0;
  }
#endif

  /* --------------------------------------------------------------------------------- */

  /* query sequence */
  seq = query->seq;
  /* local or global alignments? */
  is_local = target->isLocal;
  sc_E = (is_local) ? 0 : -INF;

  q_0 = 0;
  qx0 = q_0 % 2;

  /* initialize special states (?) */
  XMX(SP_N, q_0) = 0;                                      /* S->N, p=1             */
  XMX(SP_B, q_0) = XSC(SP_N, SP_MOVE);                     /* S->N->B, no N-tail    */
  XMX(SP_E, q_0) = XMX(SP_C, q_0) = XMX(SP_J, q_0) = -INF; /* need seq to get here  */

  /* initialize zero row (top-edge) */
  for (t_0 = 0; t_0 <= T; t_0++) {
    MMX3(qx0, t_0) = IMX3(qx0, t_0) = DMX3(qx0, t_0) = -INF; /* need seq to get here  */
  }

  /* FOR every position in QUERY seq */
  for (q_0 = 1; q_0 <= Q; q_0++) {
    q_1 = q_0 - 1;
    qx0 = q_0 % 2;
    qx1 = q_1 % 2;

    t_0 = 0;

    /* Get next character in Query */
    a = seq[q_1];
    A = AA_REV[a];

    /* Initialize zero column (left-edge) */
    MMX3(qx0, t_0) = IMX3(qx0, t_0) = DMX3(qx0, t_0) = -INF;
    XMX(SP_E, q_0) = -INF;

    /* FOR every position in TARGET profile */
    for (t_0 = 1; t_0 < T; t_0++) {
      t_1 = t_0 - 1;

      // printf("q_0,t_0 = (%d,%d)\n", q_0, t_0 );

      /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
      /* best previous state transition (match takes the diag element of each prev state) */
      prv_M = MMX3(qx1, t_1) + TSC(t_1, M2M);
      prv_I = IMX3(qx1, t_1) + TSC(t_1, I2M);
      prv_D = DMX3(qx1, t_1) + TSC(t_1, D2M);
      prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M); /* from begin match state (new alignment) */
      /* best-to-match */
      prv_best = MATH_Max(
          MATH_Max(prv_M, prv_I),
          MATH_Max(prv_D, prv_B));
      MMX3(qx0, t_0) = prv_best + MSC(t_0, A);

      /* FIND BEST PATH TO INSERT STATE (FROM MATCH OR INSERT) */
      /* previous states (match takes the left element of each state) */
      prv_M = MMX3(qx1, t_0) + TSC(t_0, M2I);
      prv_I = IMX3(qx1, t_0) + TSC(t_0, I2I);
      /* best-to-insert */
      prv_best = MATH_Max(prv_M, prv_I);
      IMX3(qx0, t_0) = prv_best + ISC(t_0, A);

      /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
      /* previous states (match takes the left element of each state) */
      prv_M = MMX3(qx0, t_1) + TSC(t_1, M2D);
      prv_D = DMX3(qx0, t_1) + TSC(t_1, D2D);
      /* best-to-delete */
      prv_best = MATH_Max(prv_M, prv_D);
      DMX3(qx0, t_0) = prv_best;

      /* UPDATE E STATE */
      prv_E = XMX(SP_E, q_0);
      prv_M = MMX3(qx0, t_0) + sc_E;
      /* best-to-e-state */
      XMX(SP_E, q_0) = MATH_Max(prv_E, prv_M);

/* embed linear row into quadratic test matrix */
#if DEBUG
      {
        MX_2D(cloud_MX, q_0, t_0) = 1.0;
        MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
        MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
        MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);
      }
#endif
    }

    /* UNROLLED FINAL LOOP ITERATION */
    t_0 = T;
    t_1 = t_0 - 1;

    /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
    /* best previous state transition (match takes the diag element of each prev state) */
    prv_M = MMX3(qx1, t_1) + TSC(t_1, M2M);
    prv_I = IMX3(qx1, t_1) + TSC(t_1, I2M);
    prv_D = DMX3(qx1, t_1) + TSC(t_1, D2M);
    prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M); /* from begin match state (new alignment) */
    /* best-to-match */
    prv_best = MATH_Max(
        MATH_Max(prv_M, prv_I),
        MATH_Max(prv_D, prv_B));
    MMX3(qx0, t_0) = prv_best + MSC(t_0, A);

    /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
    /* previous states (match takes the left element of each state) */
    prv_M = MMX3(qx0, t_1) + TSC(t_1, M2D);
    prv_D = DMX3(qx0, t_1) + TSC(t_1, D2D);
    /* best-to-delete */
    prv_best = MATH_Max(prv_M, prv_D);
    DMX3(qx0, t_0) = prv_best;

    /* UPDATE E STATE */
    prv_E = XMX(SP_E, q_0);
    prv_M = MMX3(qx0, t_0);
    prv_D = DMX3(qx0, t_0);
    XMX(SP_E, q_0) = MATH_Max(prv_E,
                              MATH_Max(prv_M, prv_D));

    /* SPECIAL STATES */
    /* J state */
    prv_J = XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP); /* J->J */
    prv_E = XMX(SP_E, q_0) + XSC(SP_E, SP_LOOP); /* E->J is E's "loop" */
    XMX(SP_J, q_0) = MATH_Max(prv_J, prv_E);

    /* C state */
    prv_C = XMX(SP_C, q_1) + XSC(SP_C, SP_LOOP);
    prv_E = XMX(SP_E, q_0) + XSC(SP_E, SP_MOVE);
    XMX(SP_C, q_0) = MATH_Max(prv_C, prv_E);

    /* N state */
    prv_N = XMX(SP_N, q_1) + XSC(SP_N, SP_LOOP);
    XMX(SP_N, q_0) = prv_N;

    /* B state */
    prv_N = XMX(SP_N, q_0) + XSC(SP_N, SP_MOVE); /* N->B is N's move */
    prv_J = XMX(SP_J, q_0) + XSC(SP_J, SP_MOVE); /* J->B is J's move */
    XMX(SP_B, q_0) = MATH_Max(prv_N, prv_J);

/* embed linear row into quadratic test matrix */
#if DEBUG
    {
      MX_2D(cloud_MX, q_0, t_0) = 1.0;
      MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
      MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
      MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);

      // MATRIX_3D_Dump( test_MX, stdout );
    }
#endif
  }

  /* T state (stores final state score) */
  sc_best = XMX(SP_C, Q) + XSC(SP_C, SP_MOVE);
  *sc_final = sc_best;

  return STATUS_SUCCESS;
}
