/*******************************************************************************
 *  - FILE:      fwdback_vectorized.c
 *  SYNOPSIS:  The Forward-Backward Algorithm for Sequence Alignment Search.
 *             ( Linear O(Q), SIMD Vectorized )
 *  NOTES:
 *    - WIP. Not implementated.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* vectorization */
#include <xmmintrin.h> /* SSE  */
#include <emmintrin.h> /* SSE2 */

/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"

/* header */
#include "_algs_vectorized.h"
#include "fwdback_vectorized.h"

STATUS_FLAG
run_Forward_Vectorized(const SEQUENCE* query,
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
  int N;         /* length of edgebound list */
  bool is_local; /* whether using local or global alignments */

  /* vars for indexing into data matrices by row-col */
  int b, d, i, j, k; /* antidiagonal, row, column indices */
  int q_0, q_1;      /* real index of current and previous rows (query sequence) */
  int qx0, qx1;      /* mod-mapped current and previous row index into data matrix (query sequence) (for linear space) */
  int t_0, t_1;      /* real index of current and previous columns (target profile) */
  int tx0, tx1;      /* vector-mapped current and previous columns (target profile) (for vectorization) */
  int Qv, Tv;        /* quad-mapped number of vectors to span each */

  /* vars for indexing into data matrices by anti-diag */
  int d_0, d_1, d_2;         /* real index of current and previous antidiagonals */
  int dx0, dx1, dx2;         /* mod mapping of antidiagonal index into data matrix */
  int k_0, k_1;              /* offset into antidiagonal */
  int d_st, d_end, d_cnt;    /* starting and ending diagonal indices */
  int dim_T, dim_Q, dim_TOT; /* dimensions of submatrix being searched */
  int dim_min, dim_max;      /* diagonal index where num cells reaches highest point and diminishing point */
  int num_cells;             /* number of cells in current diagonal */

  /* vars for indexing into edgebound lists */
  BOUND* bnd;     /* current bound */
  BOUND bnd_new;  /* for adding new bound to edgebound list */
  int id;         /* id in edgebound list (row/diag) */
  int r_0;        /* current index in edgebound list */
  int r_0b, r_0e; /* begin and end indices for current row in edgebound list */
  int r_1b, r_1e; /* begin and end indices for current row in edgebound list */
  int le_0, re_0; /* right/left matrix bounds of current diag */
  int lb_0, rb_0; /* bounds of current search space on current diag */
  int lb_1, rb_1; /* bounds of current search space on previous diag */
  int lb_2, rb_2; /* bounds of current search space on 2-back diag */
  bool rb_T;      /* checks if edge touches right bound of matrix */

  /* vars for recurrance scores */
  float prv_M, prv_I, prv_D; /* previous (M) match, (I) insert, (D) delete states */
  float prv_B, prv_E;        /* previous (B) begin and (E) end states */
  float prv_N, prv_C, prv_J; /* previous (N) initial, (C) terminal, and (J) jump states */
  float prv_loop, prv_move;  /* previous loop and move for special states */
  float prv_sum, prv_best;   /* temp subtotaling vars */
  float sc_best;             /* final best scores */
  float sc_M, sc_I, sc_D;    /* normal state scores: match, insert, delete */
  float sc_B, sc_E;          /* special states scores: begin, end */
  float sc_N, sc_C, sc_J;    /* special states scores: initial, terminal, jump */

  /* vars for vectorization */
  register __VECTOR_FLT prv_M_vec, prv_I_vec, prv_D_vec; /* previous row, "upper values" */
  register __VECTOR_FLT curr_vec;                        /* current row */
  register __VECTOR_FLT prv_col_D_vec;                   /* shifted column for "left values" for updating delete state */
  register __VECTOR_FLT prv_E_vec;                       /* previous E state => conituously updates max */
  register __VECTOR_FLT prv_B_vec;                       /* previous B state => for BEGIN->MATCH */
  __VECTOR_FLT zero_vec;                                 /* vector of all zeros:  0.0 */

  /* debugger tools */
  FILE* dbfp;
  MATRIX_2D* cloud_MX;
  MATRIX_2D* cloud_MX3;
  MATRIX_3D* test_MX;
  MATRIX_3D* test_MX3;
  EDGEBOUNDS* test_edg;
  int num_writes;
  int num_clears;

/* initialize debugging matrix */
#if DEBUG
  {
    cloud_MX = debugger->cloud_MX;
    cloud_MX3 = debugger->cloud_MX3;
    test_MX = debugger->test_MX;
    test_MX3 = debugger->test_MX3;
    test_edg = debugger->test_edg;

    MATRIX_2D_Reuse(cloud_MX, Q + 1, T + 1);
    MATRIX_2D_Fill(cloud_MX, 0);
    MATRIX_2D_Reuse(cloud_MX3, 3, (Q + 1) + (T + 1));
    MATRIX_2D_Fill(cloud_MX3, 0);
    MATRIX_3D_Reuse(test_MX, NUM_NORMAL_STATES, Q + 1, T + 1);
    MATRIX_3D_Fill(test_MX, -INF);
    MATRIX_3D_Reuse(test_MX3, NUM_NORMAL_STATES, 3, (Q + 1) + (T + 1));
    MATRIX_3D_Fill(test_MX3, -INF);
    EDGEBOUNDS_Reuse(test_edg, Q, T);

    num_writes = 0;
    num_clears = 0;
  }
#endif

  /* --------------------------------------------------------------------------------- */

  /* initialize logsum lookup table if it has not already been */
  MATH_Logsum_Init();

  /* query sequence */
  seq = query->seq;
  Qv = NUM_VEC_FOR_SEQ(Q, FLOATS_PER_VEC);
  Tv = NUM_VEC_FOR_SEQ(T, FLOATS_PER_VEC);
  /* local or global alignments? */
  is_local = target->isLocal;
  sc_E = (is_local) ? 0 : -INF;

  /* Clear top-row */
  q_0 = 0;
  qx0 = 0;

  /* initialize special states (?) */
  XMX(SP_N, q_0) = 0;                  /* S->N, p=1             */
  XMX(SP_B, q_0) = XSC(SP_N, SP_MOVE); /* S->N->B, no N-tail    */
  XMX(SP_E, q_0) = XMX(SP_C, q_0) = XMX(SP_J, q_0) = -INF;

  /* initialize 0 row (top-edge) */
  for (t_0 = 0; t_0 < T; t_0++) {
    MMX3(qx0, t_0) = IMX3(qx0, t_0) = DMX3(qx0, t_0) = -INF;
  }

  /* FOR every position in QUERY seq */
  for (q_0 = 1; q_0 <= Q; q_0++) {
    /* MAIN RECURSION */
    /* FOR every position in TARGET profile */
    for (t_0 = 1; t_0 < T; t_0++) {
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

  /* flag matrices that they contain dirty values (not -INF) */
  st_MX3->clean = false;
  sp_MX->clean = false;

  return STATUS_SUCCESS;
}

STATUS_FLAG
run_Backward_Vectorized(const SEQUENCE* query,
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
  int N;         /* length of edgebound list */
  bool is_local; /* whether using local or global alignments */

  /* vars for indexing into data matrices by row-col */
  int b, d, i, j, k; /* antidiagonal, row, column indices */
  int q_0, q_1;      /* real index of current and previous rows (query) */
  int qx0, qx1;      /* mod mapping of column index into data matrix (query) */
  int t_0, t_1;      /* real index of current and previous columns (target) */

  /* vars for indexing into data matrices by anti-diag */
  int d_0, d_1, d_2;         /* real index of current and previous antidiagonals */
  int dx0, dx1, dx2;         /* mod mapping of antidiagonal index into data matrix */
  int k_0, k_1;              /* offset into antidiagonal */
  int d_st, d_end, d_cnt;    /* starting and ending diagonal indices */
  int dim_T, dim_Q, dim_TOT; /* dimensions of submatrix being searched */
  int dim_min, dim_max;      /* diagonal index where num cells reaches highest point and diminishing point */
  int num_cells;             /* number of cells in current diagonal */

  /* v   /* vars for indexing into edgebound lists */
  BOUND* bnd;     /* current bound */
  BOUND bnd_new;  /* for adding new bound to edgebound list */
  int id;         /* id in edgebound list (row/diag) */
  int r_0;        /* current index in edgebound list */
  int r_0b, r_0e; /* begin and end indices for current row in edgebound list */
  int r_1b, r_1e; /* begin and end indices for current row in edgebound list */
  int le_0, re_0; /* right/left matrix bounds of current diag */
  int lb_0, rb_0; /* bounds of current search space on current diag */
  int lb_1, rb_1; /* bounds of current search space on previous diag */
  int lb_2, rb_2; /* bounds of current search space on 2-back diag */
  bool rb_T;      /* checks if edge touches right bound of matrix */

  /* vars for recurrance scores */
  float prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
  float prv_B, prv_E;           /* previous (B) begin and (E) end states */
  float prv_J, prv_N, prv_C;    /* previous (J) jump, (N) initial, and (C) terminal states */
  float prv_loop, prv_move;     /* previous loop and move for special states */
  float prv_sum, prv_best;      /* temp subtotaling vars */
  float sc_best;                /* final best scores */
  float sc_M, sc_I, sc_D, sc_E; /* match, insert, delete, end scores */

  /* debugger tools */
  FILE* dbfp;
  MATRIX_2D* cloud_MX;
  MATRIX_2D* cloud_MX3;
  MATRIX_3D* test_MX;
  MATRIX_3D* test_MX3;
  EDGEBOUNDS* test_edg;
  int num_writes;
  int num_clears;

/* initialize debugging matrix */
#if DEBUG
  {
    cloud_MX = debugger->cloud_MX;
    cloud_MX3 = debugger->cloud_MX3;
    test_MX = debugger->test_MX;
    test_MX3 = debugger->test_MX3;
    test_edg = debugger->test_edg;

    MATRIX_2D_Reuse(cloud_MX, Q + 1, T + 1);
    MATRIX_2D_Fill(cloud_MX, 0);
    MATRIX_2D_Reuse(cloud_MX3, 3, (Q + 1) + (T + 1));
    MATRIX_2D_Fill(cloud_MX3, 0);
    MATRIX_3D_Reuse(test_MX, NUM_NORMAL_STATES, Q + 1, T + 1);
    MATRIX_3D_Fill(test_MX, -INF);
    MATRIX_3D_Reuse(test_MX3, NUM_NORMAL_STATES, 3, (Q + 1) + (T + 1));
    MATRIX_3D_Fill(test_MX3, -INF);
    EDGEBOUNDS_Reuse(test_edg, Q, T);

    num_writes = 0;
    num_clears = 0;
  }
#endif

  /* --------------------------------------------------------------------------------- */

  /* initialize logsum lookup table if it has not already been */
  MATH_Logsum_Init();

  /* query sequence */
  seq = query->seq;
  /* local or global alignments? */
  is_local = target->isLocal;
  sc_E = (is_local) ? 0 : -INF;

  /* Initialize the Q row. */
  q_0 = Q;
  qx0 = q_0 % 2;
  t_0 = T;

  XMX(SP_J, q_0) = XMX(SP_B, q_0) = XMX(SP_N, q_0) = -INF;
  XMX(SP_C, q_0) = XSC(SP_C, SP_MOVE);
  XMX(SP_E, q_0) = XMX(SP_C, q_0) + XSC(SP_E, SP_MOVE);

  MMX3(qx0, t_0) = DMX3(qx0, t_0) = XMX(SP_E, q_0);
  IMX3(qx0, t_0) = -INF;

#if DEBUG
  {
    MX_2D(cloud_MX, q_0, t_0) = 1.0;
    MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
    MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
    MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);
  }
#endif

  for (t_0 = T - 1; t_0 >= 1; t_0--) {
    t_1 = t_0 + 1;

    prv_E = XMX(SP_E, Q) + sc_E;
    prv_D = DMX3(qx0, t_1) + TSC(t_0, M2D);
    MMX3(qx0, t_0) = MATH_Sum(prv_E, prv_D);

    prv_E = XMX(SP_E, Q) + sc_E;
    prv_D = DMX3(qx0, t_1) + TSC(t_0, D2D);
    DMX3(qx0, t_0) = MATH_Sum(prv_E, prv_D);

    IMX3(qx0, t_0) = -INF;

#if DEBUG
    {
      MX_2D(cloud_MX, q_0, t_0) = 1.0;
      MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
      MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
      MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);
    }
#endif
  }

  /* MAIN RECURSION */
  /* FOR every position in QUERY seq */
  for (q_0 = Q - 1; q_0 >= 1; q_0--) {
    q_1 = q_0 + 1;
    qx0 = q_0 % 2;
    qx1 = q_1 % 2;

    /* Get next sequence character */
    a = seq[q_0];
    A = AA_REV[a];

    /* SPECIAL STATES */

    /* B STATE -> MATCH */
    /* NOTE: When j = 0, MMX and MSC do not match HMMER p7_GBackward() implementation.   */
    XMX(SP_B, q_0) = MMX3(qx1, 1) + TSC(0, B2M) + MSC(1, A);
    for (t_0 = 2; t_0 <= T; t_0++) {
      t_1 = t_0 - 1;
      prv_sum = XMX(SP_B, q_0);
      prv_M = MMX3(qx1, t_0) + TSC(t_1, B2M) + MSC(t_0, A);
      XMX(SP_B, q_0) = MATH_Sum(prv_sum, prv_M);
    }

    prv_J = XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP);
    prv_B = XMX(SP_B, q_0) + XSC(SP_J, SP_MOVE);
    XMX(SP_J, q_0) = MATH_Sum(prv_J, prv_B);

    prv_C = XMX(SP_C, q_1) + XSC(SP_C, SP_LOOP);
    XMX(SP_C, q_0) = prv_C;

    prv_J = XMX(SP_J, q_0) + XSC(SP_E, SP_LOOP);
    prv_C = XMX(SP_C, q_0) + XSC(SP_E, SP_MOVE);
    XMX(SP_E, q_0) = MATH_Sum(prv_J, prv_C);

    prv_N = XMX(SP_N, q_1) + XSC(SP_N, SP_LOOP);
    prv_B = XMX(SP_B, q_0) + XSC(SP_N, SP_MOVE);
    XMX(SP_N, q_0) = MATH_Sum(prv_N, prv_B);

    t_0 = T;
    MMX3(qx0, T) = DMX3(qx0, T) = XMX(SP_E, q_0);
    IMX3(qx0, T) = -INF;

#if DEBUG
    {
      MX_2D(cloud_MX, q_0, t_0) = 1.0;
      MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
      MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
      MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);
    }
#endif

    /* FOR every position in TARGET profile */
    for (t_0 = T - 1; t_0 >= 1; t_0--) {
      t_1 = t_0 + 1;

      /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
      prv_M = MMX3(qx1, t_1) + TSC(t_0, M2M) + MSC(t_1, A);
      prv_I = IMX3(qx1, t_0) + TSC(t_0, M2I) + ISC(t_1, A);
      prv_D = DMX3(qx0, t_1) + TSC(t_0, M2D);
      prv_E = XMX(SP_E, q_0) + sc_E; /* from end match state (new alignment) */
      /* best-to-match */
      prv_sum = MATH_Sum(
          MATH_Sum(prv_M, prv_I),
          MATH_Sum(prv_E, prv_D));
      MMX3(qx0, t_0) = prv_sum;

      /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
      prv_M = MMX3(qx1, t_1) + TSC(t_0, I2M) + MSC(t_1, A);
      prv_I = IMX3(qx1, t_0) + TSC(t_0, I2I) + ISC(t_0, A);
      /* best-to-insert */
      prv_sum = MATH_Sum(prv_M, prv_I);
      IMX3(qx0, t_0) = prv_sum;

      /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
      prv_M = MMX3(qx1, t_1) + TSC(t_0, D2M) + MSC(t_1, A);
      prv_D = DMX3(qx0, t_1) + TSC(t_0, D2D);
      prv_E = XMX(SP_E, q_0) + sc_E;
      /* best-to-delete */
      prv_sum = MATH_Sum(prv_M,
                         MATH_Sum(prv_D, prv_E));
      DMX3(qx0, t_0) = prv_sum;

#if DEBUG
      {
        MX_2D(cloud_MX, q_0, t_0) = 1.0;
        MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
        MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
        MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);
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
  prv_M = MMX3(q_1, t_1) + TSC(0, B2M) + MSC(1, A);
  XMX(SP_B, q_0) = prv_M;
  for (t_0 = 2; t_0 <= T; t_0++) {
    t_1 = t_0 - 1;
    prv_sum = XMX(SP_B, q_0);
    prv_M = MMX3(q_1, t_0) + TSC(t_1, B2M) + MSC(t_0, A);
    XMX(SP_B, q_0) = MATH_Sum(prv_sum, prv_M);
  }

  XMX(SP_J, q_0) = -INF;
  XMX(SP_C, q_0) = -INF;
  XMX(SP_E, q_0) = -INF;

  prv_N = XMX(SP_N, q_1) + XSC(SP_N, SP_LOOP);
  prv_B = XMX(SP_B, q_0) + XSC(SP_N, SP_MOVE);
  XMX(SP_N, q_0) = MATH_Sum(prv_N, prv_B);

  for (t_0 = T; t_0 >= 1; t_0--) {
    MMX3(qx0, t_0) = IMX3(qx0, t_0) = DMX3(qx0, t_0) = -INF;
  }

#if DEBUG
  {
    MX_2D(cloud_MX, q_0, t_0) = 1.0;
    MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
    MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
    MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);
  }
#endif

  sc_best = XMX(SP_N, 0);
  *sc_final = sc_best;

  /* flag matrices that they contain dirty values (not -INF) */
  st_MX3->clean = false;
  sp_MX->clean = false;

  return STATUS_SUCCESS;
}
