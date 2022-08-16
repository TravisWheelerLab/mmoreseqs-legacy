/*******************************************************************************
 *  - FILE:  traceback_sparse.c
 *  - DESC:  Traceback for Viterbi Algorithm for Sparse Matrices.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../algs_quad/_algs_quad.h"

/* header */
#include "_algs_sparse.h"
#include "viterbi_traceback_sparse.h"

/* private functions */
static inline float
MY_Sum(const float x, const float y);

static inline float
MY_Prod(const float x, const float y);

static inline float
MY_Zero();

static inline float
MY_One();

/*! FUNCTION:  run_Viterbi_Traceback_Sparse()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *
 *    RETURN:  Return <STATUS_SUCCESS> if no errors.
 */
STATUS_FLAG
run_Viterbi_Traceback_Sparse(const SEQUENCE* query,                 /* query sequence */
                             const HMM_PROFILE* target,             /* HMM model */
                             const int Q,                           /* query/seq length */
                             const int T,                           /* target/model length */
                             EDGEBOUNDS* edg,                       /* edgebounds of sparse matrix */
                             const RANGE* dom_range,                /* (OPTIONAL) domain range for computing fwd/bck on specific domain. If NULL, computes complete fwd/bck. */
                             MATRIX_3D_SPARSE* restrict st_SMX_vit, /* Normal State (Match, Insert, Delete) Matrix */
                             MATRIX_2D* restrict sp_MX_vit,         /* Special State (J,N,B,C,E) Matrix */
                             ALIGNMENT* aln)                        /* OUTPUT: Traceback Alignment */
{
  printf("=== run_Viterbi_Traceback_Sparse() [BEGIN] ===\n");
  // exit(0);

  /* output file pointer */
  FILE* fp;

  /* generic dp matrix pointers for macros */
  MATRIX_3D_SPARSE* st_SMX = st_SMX_vit;
  MATRIX_2D* sp_MX = sp_MX_vit;

  /* vars for accessing query/target data structs */
  char a;        /* store current character in sequence */
  int A;         /* store int value of character */
  char* seq;     /* alias for getting seq */
  bool is_local; /* whether */

  /* vars for indexing into data matrices by row-col */
  int b, d, i, j, k; /* antidiagonal, row, column indices */
  int q_0, q_1;      /* real index of current and previous rows (query) */
  int qx0, qx1;      /* maps column index into data index (query) */
  int t_0, t_1;      /* real index of current and previous columns (target) */
  int tx0, tx1;      /* maps target index into data index (target)  */
  int t_range;       /* range of targets on current row */

  /* vars for indexing into edgebound lists */
  BOUND bnd;       /* current bound */
  int id;          /* id in edgebound list (row/diag) */
  int r_0;         /* current index for current row */
  int r_0b, r_0e;  /* begin and end indices for current row in edgebound list */
  int r_1;         /* current index for previous row */
  int r_1b, r_1e;  /* begin and end indices for current row in edgebound list */
  int le_0, re_0;  /* right/left matrix bounds of current diag */
  int lb_0, rb_0;  /* bounds of current search space on current diag */
  int lb_1, rb_1;  /* bounds of current search space on previous diag */
  int lb_2, rb_2;  /* bounds of current search space on 2-back diag */
  bool lb_T, rb_T; /* checks if edge touches right bound of matrix */

  /* vars for recurrance scores */
  float sc_prv, sc_cur, sc_nxt; /* current score */
  float prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
  float prv_B, prv_E;           /* previous (B) begin and (E) end states */
  float prv_J, prv_N, prv_C;    /* previous (J) jump, (N) initial, and (C) terminal states */
  float prv_loop, prv_move;     /* previous loop and move for special states */
  float prv_sum, prv_best;      /* temp subtotaling vars */
  float sc_best, sc_max;        /* final best scores */
  float sc_M, sc_I, sc_D, sc_E; /* match, insert, delete, end scores */

  /* vars for sparse matrix */
  MATRIX_3D_SPARSE* mx;
  EDGEBOUNDS* edg_inner;    /* edgebounds for search space of backward/forward */
  EDGEBOUNDS* edg_outer;    /* edgebounds for sparse matrix shape */
  int N;                    /* size of edgebounds */
  RANGE T_range;            /* target range */
  RANGE Q_range;            /* query range */
  bool is_q_0_in_dom_range; /* checks if current query position is inside the domain range */
  bool is_q_1_in_dom_range; /* checks if previous query position is inside the domain range */

  /* vars for alignment traceback */
  int st_cur, st_prv, st_nxt; /* current, previous state in traceback */
  int q_prv, t_prv;           /* previous  */
  TRACE* tr;                  /* trace object for appending */
  float min_diff;             /* minimum difference between previous and current score, accounting for transition */
  // const float    tol = 1e-4;                /* acceptable tolerance range for "equality tests" */
  const float tol = 1e-4; /* acceptable tolerance range for "equality tests" */
  bool is_found;          /* determines whether cell is found */
  bool is_multi_aln;      /* allow jumps / multiple alignments? */
  int aln_len;            /* length of alignment */
  int t_max, t_found;     /* best scoring option */
  float* ref;             /* reference cell */
  float sc_post;          /* posterior score */
  float max_error;        /* maximum error during entire traceback */

  /* --------------------------------------------------------------------------- */

  /* get sequence and trace */
  seq = query->seq;
  /* local or global? */
  is_local = target->isLocal;
  /* allow jumps / multiple alignments ? */
  is_multi_aln = false;

  /* domain range (query sequence) */
  if (dom_range == NULL) {
    Q_range.beg = 0;
    Q_range.end = Q;
  } else {
    Q_range = *dom_range;
  }
  /* target range */
  T_range.beg = 0;
  T_range.end = T + 1;

  /* Initialize traceback */
  {
    /* initial indexes */
    q_0 = q_prv = Q_range.end;
    t_0 = t_prv = 0;

    /* get edgebound range */
    r_0b = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0 + 1);
    r_0e = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0);

    /* clear memory for trace */
    ALIGNMENT_Reuse(aln, Q, T);

    /* Backtracing, so begins at the C (exit) state */
    ALIGNMENT_AppendTrace(aln, T_ST, q_0, t_0);
    ALIGNMENT_AppendTrace(aln, C_ST, q_0, t_0);
    st_prv = C_ST;
  }

  /* Run traceback until S (entry) state */
  while (st_prv != S_ST) {
    // printf("{%s,(%d,%d)}\n", STATE_NAMES[st_prv], q_0, t_0);
    /* if we have reached the end of the query sequence */
    if (q_0 == Q_range.beg) {
      ALIGNMENT_AppendTrace(aln, S_ST, q_0, t_0);
      break;
    }

    /* if q_0 has been decremented, then get edgebound range of next row  */
    r_0b = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0 + 1);
    r_0e = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0);

    /* if in the core model */
    if (st_prv == M_ST || st_prv == I_ST || st_prv == D_ST) {
      is_found = STATUS_SUCCESS;
      // fprintf(stderr, "Seeking: {%s,(%d,%d)}, r_0:(%d,%d)\n", STATE_NAMES[st_prv], q_0, t_0, r_0b, r_0e);
      /* if either t_0 or q_0 have been decremented, then lookup index offset */
      if (t_prv != t_0 || q_prv != q_0) {
        is_found = MATRIX_3D_SPARSE_GetOffset_ByCoords(st_SMX_vit, q_0, t_0, &qx1, &qx0, NULL, &tx0);
      }
      /* if not found, throw error */
      if (is_found == STATUS_FAILURE) {
        fprintf(stderr, "ERROR: Impossible position reached in posterior backtrace at { %s, (%d, %d) }\n",
                STATE_NAMES[st_prv], q_0, t_0);
        fprintf(stderr, "BOUNDS: (%d,%d)\n", r_0b, r_0e);
        for (r_0 = r_0b; r_0 > r_0e; r_0--) {
          /* find bound that contains current index */
          bnd = MATRIX_3D_SPARSE_GetBound_byIndex(st_SMX, r_0);
          fprintf(stderr, "BOUND[%d]:{%d,(%d,%d)} ", r_0, bnd.id, bnd.lb, bnd.rb);
        }
        fprintf(stderr, "\n");
        ERRORCHECK_exit(EXIT_FAILURE);
      }
    }

    /* for finding next maximal state, init max to -inf and init state to BOGUS STATE */
    st_cur = X_ST;
    sc_max = -INF;
    min_diff = INF;

    /* update previous */
    q_prv = q_0;
    t_prv = t_0;

    /* previous target and query sequence */
    q_1 = q_0 - 1;
    t_1 = t_0 - 1;
    tx1 = tx0 - 1;

    /* get next sequence character */
    a = seq[q_1];
    A = AA_REV[a];

    /* jump from current state to the prev state */
    switch (st_prv) {
      /* C STATE to {C,E} */
      case C_ST: /* C(q_0) comes from C(q_1) or E(q_0) */
      {
        /* current state */
        sc_prv = XMX(SP_C, q_0);

        if (sc_prv == -INF) {
          fprintf(stderr, "ERROR: Impossible C_ST reached at {%s,%d,%d}\n",
                  STATE_NAMES[st_prv], q_0, t_0);
          ERRORCHECK_exit(EXIT_FAILURE);
        }

        /* possible previous states */
        prv_C = MY_Prod(XMX(SP_C, q_1),
                        XSC(SP_C, SP_LOOP));
        prv_E = MY_Prod(XMX(SP_E, q_0),
                        XSC(SP_E, SP_MOVE));

        /* find transition with minimum difference */
        if (min_diff > ABS(sc_prv - prv_C)) {
          min_diff = ABS(sc_prv - prv_C);
          st_cur = C_ST;
          sc_cur = prv_C;
        }
        if (min_diff > ABS(sc_prv - prv_E)) {
          min_diff = ABS(sc_prv - prv_E);
          st_cur = E_ST;
          sc_cur = prv_E;
        }

        if (min_diff > tol) {
          fprintf(stderr, "ERROR: Failed to trace from C_ST at (%d,%d)\n", q_0, t_0);
          fprintf(stderr, "CLOSEST NEXT STATE: (%s -> %s), DIFF: %.7f, TOL: %.7f\n",
                  STATE_NAMES[st_prv], STATE_NAMES[st_cur], min_diff, tol);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      } break;

      /* E STATE to {M, D} */
      case E_ST: /* E connects from any M state. t_0 is set here. */
      {
        /* current state */
        sc_prv = XMX(SP_E, q_0);
        t_max = -1;

        if (sc_prv == -INF) {
          fprintf(stderr, "ERROR: Impossible E_ST reached at (%d,%d)\n", q_0, t_0);
          ERRORCHECK_exit(EXIT_FAILURE);
        }

        if (is_local) /* local mode: ends in M */
        {
          /* can't come from D, in a *local* Viterbi alignment. */
          st_prv = M_ST;

          /* possible previous states (any M state) */
          /* FOR every BOUND in current ROW */
          for (r_0 = r_0b; r_0 > r_0e; r_0--) {
            /* get bound data */
            bnd = MATRIX_3D_SPARSE_GetBound_byIndex(st_SMX, r_0);
            lb_0 = MAX(bnd.lb, T_range.beg); /* can't overflow left edge */
            rb_0 = MIN(bnd.rb, T_range.end); /* can't overflow right edge */

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(st_SMX, r_0);
            qx1 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Nxt(st_SMX, r_0);

            /* location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd.lb; /* total_offset = offset_location - starting_location */

            for (t_0 = lb_0; t_0 < rb_0; t_0++) {
              t_1 = t_0 + 1;
              tx0 = t_0 - bnd.lb;
              tx1 = tx0 + 1;

              /* possible previous state */
              prv_M = MSMX(qx0, tx0);
              prv_D = DSMX(qx0, tx0);

              /* verifies if scores agree with true previous state in alignment */
              if (min_diff > ABS(sc_prv - prv_M)) {
                min_diff = ABS(sc_prv - prv_M);
                st_cur = M_ST;
                sc_cur = prv_M;
                t_max = t_0;
              }
            }
          }

          /* if no entry point into M found */
          if (min_diff > tol) {
            fprintf(stderr, "ERROR: Failed to trace from E_ST at (%d,%d)\n", q_0, t_0);
            fprintf(stderr, "CLOSEST NEXT STATE: (%s -> %s), DIFF: %.7f, TOL: %.7f\n",
                    STATE_NAMES[st_prv], STATE_NAMES[st_cur], min_diff, tol);
            fprintf(stderr, "FOUND: %d\n", t_max);
            ERRORCHECK_exit(EXIT_FAILURE);
          } else {
            t_0 = t_max;
          }
        } else /*! TODO: (?) glocal mode: we either come from D_M or M_M */
        {
          fprintf(stderr, "ERROR: Glocal mode not supported for sparse alignment.");
          ERRORCHECK_exit(STATUS_FAILURE);
        }
      } break;

      /* M STATE to {B,M,I,D} */
      case M_ST: /* M connects from (q_1, t_1), or B */
      {
        /* current state */
        sc_prv = MSMX(qx0, tx0);

        /* No valid alignment goes to -INF */
        if (sc_prv == -INF) {
          fprintf(stderr, "ERROR: Impossible M_ST reached at (%d,%d)\n", q_0, t_0);
          fprintf(stderr, "CLOSEST NEXT STATE: (%s -> %s), DIFF: %.7f, TOL: %.7f\n",
                  STATE_NAMES[st_prv], STATE_NAMES[st_cur], min_diff, tol);
          ERRORCHECK_exit(EXIT_FAILURE);
        }

        /* possible previous states */
        prv_M = MY_Prod(MSMX(qx1, tx1),
                        MY_Prod(TSC(t_1, M2M), MSC(t_0, A)));
        prv_I = MY_Prod(ISMX(qx1, tx1),
                        MY_Prod(TSC(t_1, I2M), MSC(t_0, A)));
        prv_D = MY_Prod(DSMX(qx1, tx1),
                        MY_Prod(TSC(t_1, D2M), MSC(t_0, A)));
        prv_B = MY_Prod(XMX(SP_B, q_1),
                        MY_Prod(TSC(t_1, B2M), MSC(t_0, A)));

#if DEBUG
        {
          // printf("(t_0,q_0)=%d,%d | msc(t_0=%d,A=%d): %.3f\n",
          //    t_0, q_0, t_0, A, MSC(t_0, A));
          // printf("CUR_M: %.3f, B: %.3f, M: %.3f, I: %.3f, D: %.3f\n",
          //    sc_cur, XMX(SP_B, q_1), MSMX(qx1, tx1), ISMX(qx1, tx1), DSMX(qx1, tx1) );
          // printf("CUR_M: %.3f, B: %.3f, M: %.3f, I: %.3f, D: %.3f\n",
          //    sc_cur, prv_B, prv_M, prv_I, prv_D );
          // printf("CUR_M: %.3f, B: %.3f, M: %.3f, I: %.3f, D: %.3f\n",
          //    ABS( sc_cur - sc_cur ), ABS( sc_cur - prv_B), ABS( sc_cur - prv_M ), ABS( sc_cur - prv_I ), ABS( sc_cur - prv_D ) );
        }
#endif

        /* verifies if scores agree with true previous state in alignment */
        if (min_diff > ABS(sc_prv - prv_M)) {
          min_diff = ABS(sc_prv - prv_M);
          sc_cur = prv_M;
          st_cur = M_ST;
        }
        if (min_diff > ABS(sc_prv - prv_I)) {
          min_diff = ABS(sc_prv - prv_I);
          sc_cur = prv_I;
          st_cur = I_ST;
        }
        if (min_diff > ABS(sc_prv - prv_D)) {
          min_diff = ABS(sc_prv - prv_D);
          sc_cur = prv_D;
          st_cur = D_ST;
        }
        if (min_diff > ABS(sc_prv - prv_B)) {
          min_diff = ABS(sc_prv - prv_B);
          sc_cur = prv_B;
          st_cur = B_ST;
        }

        if (min_diff > tol) {
          fprintf(stderr, "ERROR: Failed to trace from M_ST at (%d,%d)\n", t_0, q_0);
          fprintf(stderr, "CLOSEST NEXT STATE: (%s -> %s), DIFF: %.7f, TOL: %.7f\n",
                  STATE_NAMES[st_prv], STATE_NAMES[st_cur], min_diff, tol);
          ERRORCHECK_exit(EXIT_FAILURE);
        }

        /* update index to previous state */
        t_0--;
        q_0--;
      } break;

      /* D STATE to {M,D} */
      case D_ST: /* D connects from M,D at (q_0, t_1) */
      {
        /* current state */
        sc_prv = DSMX(qx0, tx0);

        /* No valid alignment goes to -INF */
        if (sc_prv == -INF) {
          fprintf(stderr, "ERROR: Impossible D_ST reached at (%d,%d)\n", q_0, t_0);
          ERRORCHECK_exit(EXIT_FAILURE);
        }

        /* possible previous states */
        prv_M = MY_Prod(MSMX(qx0, tx1),
                        TSC(t_1, M2D));
        prv_D = MY_Prod(DSMX(qx0, tx1),
                        TSC(t_1, D2D));

        /* verifies if scores agree with true previous state in alignment */
        if (min_diff > ABS(sc_prv - prv_M)) {
          min_diff = ABS(sc_prv - prv_M);
          st_cur = M_ST;
          sc_cur = prv_M;
        }
        if (min_diff > ABS(sc_prv - prv_D)) {
          min_diff = ABS(sc_prv - prv_D);
          st_cur = D_ST;
          sc_cur = prv_D;
        }

        if (min_diff > tol) {
          fprintf(stderr, "ERROR: Failed to trace from D_ST at (%d,%d)\n", q_0, t_0);
          fprintf(stderr, "CLOSEST NEXT STATE: (%s -> %s), DIFF: %.7f, TOL: %.7f\n",
                  STATE_NAMES[st_prv], STATE_NAMES[st_cur], min_diff, tol);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
        /* update index to previous state */
        t_0--;
      } break;

      /* I STATE to {M,I} */
      case I_ST: /* I connects from M,I at (q_1, t_0) */
      {
        /* current state */
        sc_prv = ISMX(qx0, tx0);

        /* No valid alignment goes to -INF */
        if (sc_prv == -INF) {
          fprintf(stderr, "ERROR: Impossible I_ST reached at (%d,%d)\n", q_0, t_0);
          ERRORCHECK_exit(EXIT_FAILURE);
        }

        /* possible previous states */
        prv_M = MY_Prod(MSMX(qx1, tx0),
                        MY_Prod(TSC(t_0, M2I), ISC(t_0, A)));
        prv_I = MY_Prod(ISMX(qx1, tx0),
                        MY_Prod(TSC(t_0, I2I), ISC(t_0, A)));

        /* verifies if scores agree with true previous state in alignment */
        if (min_diff > ABS(sc_prv - prv_M)) {
          min_diff = ABS(sc_prv - prv_M);
          st_cur = M_ST;
          sc_cur = prv_M;
        }
        if (min_diff > ABS(sc_prv - prv_I)) {
          min_diff = ABS(sc_prv - prv_I);
          st_cur = I_ST;
          sc_cur = prv_I;
        }

        if (min_diff > tol) {
          fprintf(stderr, "ERROR: Failed to trace from I_ST at (%d,%d)\n", q_0, t_0);
          fprintf(stderr, "CLOSEST NEXT STATE: (%s -> %s), DIFF: %.7f, TOL: %.7f\n",
                  STATE_NAMES[st_prv], STATE_NAMES[st_cur], min_diff, tol);
          ERRORCHECK_exit(EXIT_FAILURE);
        }

        /* update index to previous state */
        q_0--;
      } break;

      /* N STATE to {N,S} */
      case N_ST: /* N connects from S,N */
      {
        /* current state */
        sc_prv = XMX(SP_N, q_0);

        /* No valid alignment goes to -INF */
        if (sc_prv == -INF) {
          fprintf(stderr, "ERROR: Impossible N_ST reached at (%d,%d)\n", q_0, t_0);
          ERRORCHECK_exit(EXIT_FAILURE);
        }

        /* if at beginning of query sequence, then alignment completes at S state, else N state */
        if (q_0 <= 0) {
          st_cur = S_ST;
        } else {
          st_cur = N_ST;
        }
      } break;

      /* B STATE to {N,J} */
      case B_ST: /* B connects from N,J */
      {
        /* current state */
        sc_prv = XMX(SP_B, q_0);

        /* possible previous states */
        prv_N = MY_Prod(XMX(SP_N, q_0),
                        XSC(SP_N, SP_MOVE));
        prv_J = MY_Prod(XMX(SP_J, q_0),
                        XSC(SP_J, SP_MOVE));

        /* verifies that scores agree with true previous state in alignment */
        if (min_diff > ABS(sc_prv - prv_N)) {
          min_diff = ABS(sc_prv - prv_N);
          st_cur = N_ST;
          sc_cur = prv_N;
        }
        if (min_diff > ABS(sc_prv - prv_J)) {
          min_diff = ABS(sc_prv - prv_J);
          st_cur = J_ST;
          sc_cur = prv_J;
        }

        if (min_diff > tol) {
          fprintf(stderr, "ERROR: Failed to trace from B_ST at (%d,%d)\n", q_0, t_0);
          fprintf(stderr, "CLOSEST NEXT STATE: (%s -> %s), DIFF: %.7f, TOL: %.7f\n",
                  STATE_NAMES[st_prv], STATE_NAMES[st_cur], min_diff, tol);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      } break;

      /* J STATE to {J,E} */
      case J_ST: /* J connects from E(q_0) or J(q_1) */
      {
        /* current state */
        sc_prv = XMX(SP_J, q_0);

        if (sc_prv == -INF) {
          fprintf(stderr, "ERROR: Impossible J_ST reached at (%d,%d)\n", q_0, t_0);
          ERRORCHECK_exit(EXIT_FAILURE);
        }

        /* possible previous states */
        prv_J = MY_Prod(XMX(SP_J, q_1), XSC(SP_J, SP_LOOP));
        prv_E = MY_Prod(XMX(SP_E, q_0), XSC(SP_E, SP_LOOP));

        /* verifies that scores agree with true previous state in alignment */
        if (min_diff > ABS(sc_prv - prv_J)) {
          min_diff = ABS(sc_prv - prv_J);
          st_cur = J_ST;
          sc_cur = prv_J;
        }
        if (min_diff > ABS(sc_prv - prv_E)) {
          min_diff = ABS(sc_prv - prv_E);
          st_cur = E_ST;
          sc_cur = prv_E;
        }
        if (min_diff > tol) {
          fprintf(stderr, "ERROR: Failed to trace from J_ST at (%d,%d)\n", q_0, t_0);
          fprintf(stderr, "CLOSEST NEXT STATE: (%s -> %s), DIFF: %.7f, TOL: %.7f\n",
                  STATE_NAMES[st_prv], STATE_NAMES[st_cur], min_diff, tol);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      } break;

      default: {
        fprintf(stderr, "ERROR: Hit Bogus State!!! => %d\n", st_prv);
        ERRORCHECK_exit(EXIT_FAILURE);
      }
    }

#if DEBUG
    {
      // printf("ADDED: {%s, (%d, %d)}\n", STATE_NAMES[st_cur], q_0, t_0);
    }
#endif

    /* push new trace onto the alignment */
    ALIGNMENT_AppendScore(aln, sc_cur);
    ALIGNMENT_AppendTrace(aln, st_cur, q_0, t_0);

    /* For {N,C,J}, we deferred q_0 decrement. */
    if ((st_cur == N_ST || st_cur == J_ST || st_cur == C_ST) && (st_cur == st_prv)) {
      q_0--;
    }

    fprintf(stderr, "{%s,(%d,%d)} -> {%s,(%d,%d)}\n",
            STATE_NAMES[st_prv], q_prv, t_prv, STATE_NAMES[st_cur], q_0, t_0);

    /* Update previous state */
    st_prv = st_cur;
    sc_prv = sc_cur;
  }

  /* reverse order of traceback */
  ALIGNMENT_Reverse(aln);

  /* Set main model regions of alignment */
  ALIGNMENT_FindRegions(aln);
  ALIGNMENT_ScoreRegions(aln);

  /* find end and begin alignment points (first and last match state) */
  N = aln->traces->N;
  tr = aln->traces->data;
  for (i = 0; i < N; i++) {
    if (tr[i].st == B_ST) {
      VECTOR_INT_Pushback(aln->tr_beg, i + 1);
    }
    if (tr[i].st == E_ST) {
      VECTOR_INT_Pushback(aln->tr_end, i - 1);
    }
  }
  aln->beg = aln->tr_beg->data[0];
  aln->end = aln->tr_end->data[0];

#if DEBUG
  {
    MATRIX_2D* cloud_MX = debugger->cloud_MX;
    MATRIX_2D_Reuse(cloud_MX, Q + 1, T + 1);
    MATRIX_2D_Fill(cloud_MX, 0);
    for (int i = 0; i < N; i++) {
      if (tr[i].st == M_ST || tr[i].st == I_ST || tr[i].st == D_ST)
        MX_2D(cloud_MX, tr[i].q_0, tr[i].t_0) = -1.0;
    }
  }
#endif

  return STATUS_SUCCESS;
}

/* MATH RULES: These determine how probilities are summed, multiplied, and certain identities */

static inline float
MY_Sum(const float x,
       const float y) {
  return MATH_Sum(x, y);
}

static inline float
MY_Prod(const float x,
        const float y) {
  return MATH_Prod(x, y);
}

static inline float
MY_Max(const float x,
       const float y) {
  return MATH_Max(x, y);
}

static inline float
MY_Zero() {
  return MATH_Zero();
}

static inline float
MY_One() {
  return MATH_One();
}

static inline float
MY_Compare_Tol(const float x,
               const float y) {
  return MATH_Max(x, y);
}
