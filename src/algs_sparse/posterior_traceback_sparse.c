/*******************************************************************************
 *  - FILE:  posterior_traceback_sparse.h
 *  - DESC:  The Maximum Posterior Probability and Optimal Alignment.
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
#include "../algs_linear/_algs_linear.h"
#include "../parsers/_parsers.h"

/* header */
#include "_algs_sparse.h"
#include "posterior_traceback_sparse.h"

/* This determines whether optimal alignment will be computed for cell, prevent NaN errors */
#define TSC_DELTA(t_0, tr) (KRON_DELTA(TSC((t_0), (tr))))
#define XSC_DELTA(sp, tr) (KRON_DELTA(XSC((sp), (tr))))
/* Kronecker Delta: 1 if finite, else FLT_MIN instead of -INF to prevent NaN errors. */
#define KRON_DELTA(val) (((val) == -INF) ? FLT_MIN : 1.0)

/* Math functions */
static inline float
MY_Sum(const float x, const float y);
static inline float
MY_Prod(const float x, const float y);
static inline float
MY_Zero();
static inline float
MY_One();

/* Select Patch from Previous State */
static inline STATUS_FLAG
Select_Path_from_M(const SEQUENCE* query,
                   const HMM_PROFILE* target,
                   const MATRIX_3D_SPARSE* st_SMX,
                   const MATRIX_2D* sp_MX,
                   int q_0,
                   int t_0,
                   int qx0,
                   int tx0,
                   int st_cur,
                   float sc_cur,
                   int* st_prv,
                   float* sc_prv,
                   float min_diff) {
  float prv_M, prv_I, prv_D, prv_B;
  int q_1, t_1;
  int qx1, tx1;

  int st_paths[] = {M_ST, I_ST, D_ST, B_ST};
  float sc_paths[NUM_ALL_STATES];

  q_1 = q_0 - 1;
  t_1 = t_0 - 1;
  qx1 = qx0 - 1;
  tx1 = tx1 - 1;
}

static inline STATUS_FLAG
Select_Path_from_I(MATRIX_3D_SPARSE* st_SMX,
                   MATRIX_2D* sp_MX,
                   int q_0,
                   int t_0,
                   int qx0,
                   int tx0,
                   int st_cur,
                   float sc_cur,
                   int* st_prv,
                   float* sc_prv,
                   float min_diff) {
  int paths[] = {M_ST, I_ST};
}

static inline STATUS_FLAG
Select_Path_from_D(MATRIX_3D_SPARSE* st_SMX,
                   MATRIX_2D* sp_MX,
                   int q_0,
                   int t_0,
                   int qx0,
                   int tx0,
                   int st_cur,
                   float sc_cur,
                   int* st_prv,
                   float* sc_prv,
                   float min_diff) {
  int paths[] = {M_ST, D_ST};
}

static inline STATUS_FLAG
Select_Path_from_N(MATRIX_3D_SPARSE* st_SMX,
                   MATRIX_2D* sp_MX,
                   int q_0,
                   int t_0,
                   int qx0,
                   int tx0,
                   int st_cur,
                   float sc_cur,
                   int* st_prv,
                   float* sc_prv,
                   float min_diff) {
  int paths[] = {N_ST, S_ST};
}

static inline STATUS_FLAG
Select_Path_from_C(MATRIX_3D_SPARSE* st_SMX,
                   MATRIX_2D* sp_MX,
                   int q_0,
                   int t_0,
                   int qx0,
                   int tx0,
                   int st_cur,
                   float sc_cur,
                   int* st_prv,
                   float* sc_prv,
                   float min_diff) {
  int paths[] = {C_ST, E_ST};
}

static inline STATUS_FLAG
Select_Path_from_J(MATRIX_3D_SPARSE* st_SMX,
                   MATRIX_2D* sp_MX,
                   int q_0,
                   int t_0,
                   int qx0,
                   int tx0,
                   int st_cur,
                   float sc_cur,
                   int* st_prv,
                   float* sc_prv,
                   float min_diff) {
  int paths[] = {J_ST, E_ST};
}

static inline STATUS_FLAG
Select_Path_from_E(MATRIX_3D_SPARSE* st_SMX,
                   MATRIX_2D* sp_MX,
                   int q_0,
                   int t_0,
                   int qx0,
                   int tx0,
                   int st_cur,
                   float sc_cur,
                   int* st_prv,
                   float* sc_prv,
                   float min_diff) {
  int paths[] = {M_ST, D_ST};
}

static inline STATUS_FLAG
Select_Path_from_B(MATRIX_3D_SPARSE* st_SMX,
                   MATRIX_2D* sp_MX,
                   int q_0,
                   int t_0,
                   int qx0,
                   int tx0,
                   int st_cur,
                   float sc_cur,
                   int* st_prv,
                   float* sc_prv,
                   float min_diff) {
  int paths[] = {N_ST, J_ST};
}

static inline STATUS_FLAG
MinDiff_Update(float sc_prv,
               float sc_new,
               int st_prv,
               int st_new,
               float* diff_min,
               float* sc_min,
               int* st_min) {
  float diff_new = ABS(sc_prv - sc_new);
  if (*diff_min > diff_new) {
    *diff_min = diff_new;
    *sc_min = sc_new;
    *st_min = st_new;
  }
}

static inline STATUS_FLAG
MinDiff_Position_Update(float sc_prv,
                        float sc_new,
                        int st_prv,
                        int st_new,
                        float* diff_min,
                        float* sc_min,
                        int* st_min) {
  float diff_new = ABS(sc_prv - sc_new);
  if (*diff_min > diff_new) {
    *diff_min = diff_new;
    *sc_min = sc_new;
    *st_min = st_new;
  }
}

static inline STATUS_FLAG
MinDiff_CheckError(int q_0,
                   int t_0,
                   float sc_prv,
                   int st_prv,
                   float sc_cur,
                   int st_cur,
                   float diff_min,
                   float diff_tol) {
  /* verifies that scores agree with true previous state in alignment */
  if (diff_min > diff_tol) {
    fprintf(stderr, "ERROR: Failed to trace from M_ST at (%d,%d)\n", t_0, q_0);
    fprintf(stderr, "CLOSEST NEXT STATE: (%s -> %s), DIFF: %.7f, TOL: %.7f\n",
            STATE_NAMES[st_prv], STATE_NAMES[st_cur], diff_min, diff_tol);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
}

static inline STATUS_FLAG
MaxScore_Update(float sc_new,
                int st_new,
                float* sc_max,
                int* st_max) {
  if (*sc_max < sc_new) {
    *sc_max = sc_new;
    *st_max = st_new;
  }
}

static inline STATUS_FLAG
MaxScore_Position_Update(float sc_new,
                         int st_new,
                         int pos_new,
                         float* sc_max,
                         int* st_max,
                         int* pos_max) {
  if (*sc_max < sc_new) {
    *sc_max = sc_new;
    *st_max = st_new;
    *pos_max = pos_new;
  }
}

STATUS_FLAG
run_OptimalAccuracy_Traceback_Sparse(const SEQUENCE* query,         /* query sequence */
                                     const HMM_PROFILE* target,     /* target hmm model */
                                     const int Q,                   /* query length */
                                     const int T,                   /* target length */
                                     EDGEBOUNDS* edg,               /* edgebounds */
                                     RANGE* dom_range,              /* query span of bounds */
                                     MATRIX_3D_SPARSE* st_SMX_post, /* posterior normal matrix */
                                     MATRIX_2D* sp_MX_post,         /* posterior special matrix */
                                     MATRIX_3D_SPARSE* st_SMX_opt,  /* optimal accuracy normal matrix */
                                     MATRIX_2D* sp_MX_opt,          /* optimal accuracy special matrix */
                                     ALIGNMENT* aln)                /* OUTPUT: optimal alignment */
{
  printf("=== run_Posterior_Traceback_Sparse() [BEGIN] ===\n");
  FILE* fp;

  /* generic dp matrix pointers for macros */
  MATRIX_3D_SPARSE* st_SMX = st_SMX_opt;
  MATRIX_2D* sp_MX = sp_MX_opt;

  /* vars for accessing query/target data structs */
  char a;        /* store current character in sequence */
  int A;         /* store int value of character */
  char* seq;     /* alias for getting seq */
  bool is_local; /* whether local or global alignment */

  /* vars for indexing into data matrices by row-col */
  int b, d, i, j, k; /* antidiagonal, row, column indices */
  int q_0, q_1;      /* real index of current and previous rows (query) */
  int qx0, qx1;      /* maps column index into data index (query) */
  int t_0, t_1;      /* real index of current and previous columns (target) */
  int tx0, tx1;      /* maps target index into data index (target)  */
  int t_range;       /* range of targets on current row */

  /* vars for indexing into edgebound lists */
  BOUND bnd;       /* current bound */
  int id_0;        /* id in edgebound list (row/diag) */
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
  float sc_prv, sc_cur, sc_max; /* current and maximum score */
  float prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
  float prv_B, prv_E;           /* previous (B) begin and (E) end states */
  float prv_J, prv_N, prv_C;    /* previous (J) jump, (N) initial, and (C) terminal states */
  float prv_loop, prv_move;     /* previous loop and move for special states */
  float prv_max, prv_best;      /* temp subtotaling vars */
  float sc_best;                /* final best scores */
  float sc_M, sc_I, sc_D, sc_E; /* match, insert, delete, end scores */

  /* vars for sparse matrix */
  MATRIX_3D_SPARSE* mx;     /* */
  EDGEBOUNDS* edg_inner;    /* edgebounds for search space of backward/forward */
  EDGEBOUNDS* edg_outer;    /* edgebounds for sparse matrix shape */
  RANGE T_range;            /* target range */
  RANGE Q_range;            /* query range */
  bool is_q_0_in_dom_range; /* checks if current query position is inside the domain range */
  bool is_q_1_in_dom_range; /* checks if previous query position is inside the domain range */

  /* vars for alignment traceback */
  int st_cur, st_prv, st_max; /* current, previous state in traceback */
  int q_prv, t_prv;           /* previous and maximum query and target positions */
  TRACE* tr;                  /* trace object for appending */
  float min_diff;             /* minimum difference between previous and current score, accounting for transition */
  // const float    tol = 1e-4;                /* acceptable tolerance range for "equality tests" */
  const float tol = 1; /* acceptable tolerance range for "equality tests" */
  bool is_found;       /* determines whether cell is found */
  bool is_multi_aln;   /* allow jumps / multiple alignments? */
  int aln_len;         /* length of alignment */
  int t_max, t_found;  /* best scoring option */
  float* ref;          /* reference cell */
  float sc_post;       /* posterior score */
  float max_error;     /* maximum error during entire traceback */

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

  /* intial indexes */
  q_0 = q_prv = Q_range.end;
  t_0 = t_prv = 0;

  /* get edgebound range */
  r_0b = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0 + 1);
  r_0e = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0);

  /* clear memory for trace */
  ALIGNMENT_Reuse(aln, Q, T);
  sc_cur = 0.0f;

  /* Backtracing, so begins at T, the model exit state */
  ALIGNMENT_AppendTrace(aln, T_ST, q_0, t_0);
  ALIGNMENT_AppendScore(aln, 0.0f);

  /* Immediate after, it transitions to C, the post main model / post END state */
  ALIGNMENT_AppendTrace(aln, C_ST, q_0, t_0);
  ALIGNMENT_AppendScore(aln, 0.0f);
  st_prv = C_ST;

  int cnt = 0;

  /* End of trace is S, the model entry state */
  while (st_prv != S_ST) {
    cnt++;
    printf("[%d]: (q_0,t_0) = (%d,%d) | st = %d = %s | sc = %.9f | diff = %.9f, %.9f\n",
           cnt, q_0, t_0, st_prv, STATE_NAMES[st_prv], sc_cur, min_diff, max_error);

    /* if we have reached the end of the query sequence */
    if (q_0 == Q_range.beg) {
      ALIGNMENT_AppendTrace(aln, S_ST, q_0, t_0);
      ALIGNMENT_AppendScore(aln, 0.0f);
      break;
    }

    /* if q_0 has been decremented, then get edgebound range of next row  */
    r_0b = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0 + 1);
    r_0e = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0);

    /* if in the core model, find location */
    if (st_prv == M_ST || st_prv == I_ST || st_prv == D_ST) {
      is_found = STATUS_SUCCESS;
      // fprintf(stderr, "Seeking: {%s,(%d,%d)}, r_0:(%d,%d)\n", STATE_NAMES[st_prv], q_0, t_0, r_0b, r_0e);
      /* if either t_0 or q_0 have been decremented, then lookup index offset */
      if (t_prv != t_0 || q_prv != q_0) {
        is_found = MATRIX_3D_SPARSE_GetOffset_ByCoords(st_SMX_post, q_0, t_0, &qx1, &qx0, NULL, &tx0);
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
    st_max = X_ST;
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

        /* possible previous states */
        prv_C = XSC_DELTA(SP_C, SP_LOOP) *
                XMX(SP_C, q_1);
        prv_E = XSC_DELTA(SP_E, SP_MOVE) *
                XMX(SP_E, q_0);

        // MinDiff_Update(sc_prv, prv_C, C_ST, C_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_Update(sc_prv, prv_E, C_ST, E_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_CheckError(q_0, t_0, sc_prv, st_prv, sc_cur, st_cur, min_diff, tol);

        MaxScore_Update(prv_C, C_ST, &sc_max, &st_max);
        MaxScore_Update(prv_E, E_ST, &sc_max, &st_max);
        printf("st_max: %d\n", st_max);
      } break;

      /* E STATE to {M, D} */
      case E_ST: /* E connects from any M state. t_0 is set here. */
      {
        /* current state */
        sc_prv = XMX(SP_E, q_0);
        t_max = -1;

        if (is_local) /* local mode: ends in M */
        {
          /* best possible previous M or D state */
          /* FOR every BOUND in current ROW */
          for (r_0 = r_0b; r_0 > r_0e; r_0--) {
            /* get bound data */
            bnd = MATRIX_3D_SPARSE_GetBound_byIndex(st_SMX, r_0);
            lb_0 = MAX(bnd.lb, T_range.beg); /* can't overflow left edge */
            rb_0 = MIN(bnd.rb, T_range.end); /* can't overflow right edge */

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(st_SMX, r_0); /* (q_0, t_0) location offset */
            qx1 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Prv(st_SMX, r_0); /* (q_1, t_0) location offset */

            /* initial location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd.lb; /* total_offset = offset_location - starting_location */

            for (t_0 = lb_0; t_0 < rb_0; t_0++) {
              t_1 = t_0 - 1;
              tx0 = t_0 - bnd.lb;
              tx1 = tx0 - 1;

              /* possible previous state */
              prv_M = MSMX(qx0, tx0);
              prv_D = DSMX(qx0, tx0);

              // printf("prv_M, prv_D [%d,%d]: { %.7f, %.7f }\n", q_0, t_0, prv_M, prv_D);

              // if (min_diff > ABS(sc_prv - prv_M))
              // {
              //    min_diff = ABS(sc_prv - prv_M);
              //    st_cur = M_ST;
              //    sc_cur = prv_M;
              //    t_max = t_0;
              // }

              if (sc_max < prv_M) {
                sc_max = prv_M;
                st_cur = M_ST;
                sc_cur = prv_M;
                t_max = t_0;
              }
            }
          }

          /* verifies if scores agree with true previous state in alignment */
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
        /* current score */
        sc_prv = MSMX(qx0, tx0);

        /* possible previous states */
        prv_M = TSC_DELTA(t_1, M2M) *
                MSMX(qx1, tx1);
        prv_I = TSC_DELTA(t_1, I2M) *
                ISMX(qx1, tx1);
        prv_D = TSC_DELTA(t_1, D2M) *
                DSMX(qx1, tx1);
        prv_B = TSC_DELTA(t_1, B2M) *
                XMX(SP_B, q_1);

        // fprintf(stdout, "CUR: %.3f, M: %.3f, I: %.3f, D: %.3f, B: %.3f\n",
        //         MSMX(qx0, tx0), MSMX(qx1, tx1), ISMX(qx1, tx1), DSMX(qx1, tx1), XMX(SP_B, q_1));

        // MinDiff_Update(sc_prv, prv_M, M_ST, M_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_Update(sc_prv, prv_I, M_ST, I_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_Update(sc_prv, prv_D, M_ST, D_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_Update(sc_prv, prv_B, M_ST, B_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_CheckError(q_0, t_0, sc_prv, st_prv, sc_cur, st_cur, min_diff, tol);

        MaxScore_Update(prv_M, M_ST, &sc_max, &st_max);
        MaxScore_Update(prv_I, I_ST, &sc_max, &st_max);
        MaxScore_Update(prv_D, D_ST, &sc_max, &st_max);
        MaxScore_Update(prv_B, B_ST, &sc_max, &st_max);

        /* update index to previous state */
        t_0--;
        q_0--;
      } break;

      /* D STATE to {M,D} */
      case D_ST: /* D connects from M,D at (q_0, t_1) */
      {
        /* current state */
        sc_prv = DSMX(qx0, tx0);

        /* possible previous states */
        prv_M = TSC_DELTA(t_1, M2D) *
                MSMX(qx0, tx1);
        prv_D = TSC_DELTA(t_1, D2D) *
                DSMX(qx0, tx1);

        // MinDiff_Update(sc_prv, prv_M, D_ST, M_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_Update(sc_prv, prv_D, D_ST, D_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_CheckError(q_0, t_0, sc_prv, st_prv, sc_cur, st_cur, min_diff, tol);

        MaxScore_Update(prv_M, M_ST, &sc_max, &st_max);
        MaxScore_Update(prv_D, D_ST, &sc_max, &st_max);

        /* update index to previous state */
        t_0--;
      } break;

      /* I STATE to {M,I} */
      case I_ST: /* I connects from M,I at (q_1, t_0) */
      {
        /* current state */
        sc_prv = ISMX(qx0, tx0);

        /* possible previous states */
        prv_M = TSC_DELTA(t_0, M2I) *
                MSMX(qx1, tx0);
        prv_I = TSC_DELTA(t_0, I2I) *
                ISMX(qx1, tx0);

        // MinDiff_Update(sc_prv, prv_M, I_ST, M_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_Update(sc_prv, prv_I, I_ST, I_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_CheckError(q_0, t_0, sc_prv, st_prv, sc_cur, st_cur, min_diff, tol);

        MaxScore_Update(prv_M, M_ST, &sc_max, &st_max);
        MaxScore_Update(prv_I, I_ST, &sc_max, &st_max);

        /* update index to previous state */
        q_0--;
      } break;

      /* N STATE to {N,S} */
      case N_ST: /* N connects from S,N */
      {
        /* current state */
        sc_prv = XMX(SP_N, q_0);

        /* if at beginning of query sequence, then alignment completes at S state, else loops on N state */
        if (q_0 <= 0) {
          st_cur = S_ST;
        } else {
          st_cur = N_ST;
        }
      } break;

      /* B STATE to {N,J} */
      case B_ST: /* B connects from N, J */
      {
        /* current state */
        sc_prv = XMX(SP_B, q_0);

        /* possible previous states */
        prv_N = XSC_DELTA(SP_N, SP_MOVE) *
                XMX(SP_N, q_0);
        prv_J = XSC_DELTA(SP_J, SP_MOVE) *
                XMX(SP_J, q_0);

        // MinDiff_Update(sc_prv, prv_N, B_ST, N_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_Update(sc_prv, prv_J, B_ST, J_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_CheckError(q_0, t_0, sc_prv, st_prv, sc_cur, st_cur, min_diff, tol);

        MaxScore_Update(prv_N, N_ST, &sc_max, &st_max);
        MaxScore_Update(prv_J, J_ST, &sc_max, &st_max);
      } break;

      /* J STATE to {J,E} */
      case J_ST: /* J connects from E(q_0) or J(q_1) */
      {
        /* current state */
        sc_prv = XMX(SP_J, q_0);

        /* possible previous states */
        prv_J = XSC_DELTA(SP_J, SP_LOOP) *
                XMX(SP_J, q_1);
        prv_E = XSC_DELTA(SP_E, SP_LOOP) *
                XMX(SP_E, q_0);

        // MinDiff_Update(sc_prv, prv_J, J_ST, J_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_Update(sc_prv, prv_E, J_ST, E_ST, &min_diff, &sc_cur, &st_cur);
        // MinDiff_CheckError(q_0, t_0, sc_prv, st_prv, sc_cur, st_cur, min_diff, tol);

        MaxScore_Update(prv_J, J_ST, &sc_max, &st_max);
        MaxScore_Update(prv_E, E_ST, &sc_max, &st_max);
      } break;

      default: {
        fprintf(stderr, "ERROR: Hit Bogus State!!! => %d\n", st_prv);
        ERRORCHECK_exit(EXIT_FAILURE);
      }
    }

    // /* error check: if no state was selected for next in trace. */
    st_cur = st_max;
    if (st_cur == X_ST) {
      fprintf(stderr, "ERROR: Traceback failed from state %s==%d:(%d,%d).\n",
              STATE_NAMES[st_prv], st_prv, q_0, t_0);
      fprintf(stderr, "ERROR: Traceback failed to state %s==%d:(%d,%d).\n",
              STATE_NAMES[st_cur], st_cur, q_0, t_0);
      exit(EXIT_FAILURE);
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

    // fprintf(stderr, "{%s,(%d,%d)} -> {%s,(%d,%d)} => ERROR: %.9f\n",
    //         STATE_NAMES[st_prv], q_prv, t_prv, STATE_NAMES[st_cur], q_0, t_0, min_diff);

    /* Update previous state */
    st_prv = st_cur;
    sc_prv = sc_cur;

    max_error = MAX(max_error, min_diff);
  }

  /* reverse order of traceback */
  ALIGNMENT_Reverse(aln);

  /* Set main model regions of alignment */
  ALIGNMENT_FindRegions(aln);
  ALIGNMENT_ScoreRegions(aln);

  int N = VECTOR_INT_GetSize(aln->tr_beg);
  VECTOR_FLT_Dump_byOpt(aln->tr_score, ",", "TR_SCORE", stdout);
  VECTOR_INT_Dump_byOpt(aln->tr_beg, ",", "TR_BEG", stdout);
  VECTOR_INT_Dump_byOpt(aln->tr_end, ",", "TR_END", stdout);
  for (int i = 0; i < N; i++) {
    TRACE tr = ALIGNMENT_GetTrace(aln, i);
    fprintf(stderr, "-%s-", STATE_CHARS[tr.st]);
    if (((i + 1) % 25) == 0) {
      fprintf(stderr, "\n");
    }
  }
  printf("\n");
  exit(0);

#if DEBUG
  {
    MATRIX_2D* cloud_MX = debugger->cloud_MX;
    MATRIX_2D_Reuse(cloud_MX, Q + 1, T + 1);
    MATRIX_2D_Fill(cloud_MX, 0);
    for (int i = 0; i < N; i++) {
      TRACE tr = ALIGNMENT_GetTrace(aln, i);
      if (tr.st == M_ST || tr.st == I_ST || tr.st == D_ST)
        MX_2D(cloud_MX, tr.q_0, tr.t_0) = -1.0;
    }
  }
#endif

  return STATUS_SUCCESS;
}

STATUS_FLAG
run_OptimalAccuracy_Sparse(const SEQUENCE* query,         /* query sequence */
                           const HMM_PROFILE* target,     /* target HMM model */
                           const int Q,                   /* query length */
                           const int T,                   /* target length */
                           EDGEBOUNDS* edg,               /* edgebounds */
                           RANGE* dom_range,              /* OPTIONAL: domain range for computing fwd/bck on specific domain. If NULL, computes complete fwd/bck. */
                           MATRIX_3D_SPARSE* st_SMX_post, /* posterior normal state matrix */
                           MATRIX_2D* sp_MX_post,         /* posterior special state matrix */
                           MATRIX_3D_SPARSE* st_SMX_opt,  /* OUTPUT: optimal normal state matrix */
                           MATRIX_2D* sp_MX_opt,          /* OUTPUT: optimal special state matrix */
                           float* sc_final)               /* OUTPUT: final score */
{
  printf("=== run_Posterior_Optimal_Accuracy_Sparse() [BEGIN] ===\n");

  /* output file pointer */
  FILE* fp;

  /* vars for matrix access for macros */
  MATRIX_3D_SPARSE* st_SMX = st_SMX_opt; /* normal state matrix */
  MATRIX_2D* sp_MX = sp_MX_opt;          /* special state matrix */

  /* vars for accessing query/target data structs */
  char a;        /* store current character in sequence */
  int A;         /* store int value of character */
  char* seq;     /* alias for getting seq */
  int N;         /* length of edgebound list */
  bool is_local; /* whether using local or global alignments */

  /* vars for indexing into data matrices by row-col */
  int b, d, i, j, k; /* antidiagonal, row, column indices */
  int q_0, q_1;      /* real index of current and previous rows (query) */
  int qx0, qx1;      /* maps column index into data index (query) */
  int t_0, t_1;      /* real index of current and previous columns (target) */
  int tx0, tx1;      /* maps target index into data index (target)  */
  int t_range;       /* range of targets on current row */

  /* vars for indexing into data matrices by anti-diag */
  int d_0, d_1, d_2;         /* real index of current and previous antidiagonals */
  int dx0, dx1, dx2;         /* mod mapping of antidiagonal index into data matrix */
  int k_0, k_1;              /* offset into antidiagonal */
  int d_st, d_end, d_cnt;    /* starting and ending diagonal indices */
  int dim_T, dim_Q, dim_TOT; /* dimensions of submatrix being searched */
  int dim_min, dim_max;      /* diagonal index where num cells reaches highest point and diminishing point */
  int num_cells;             /* number of cells in current diagonal */

  /* vars for indexing into edgebound lists */
  BOUND bnd;       /* current bound */
  int id_0;        /* id in edgebound list (row/diag) */
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
  float prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
  float prv_B, prv_E;           /* previous (B) begin and (E) end states */
  float prv_J, prv_N, prv_C;    /* previous (J) jump, (N) initial, and (C) terminal states */
  float prv_loop, prv_move;     /* previous loop and move for special states */
  float prv_max, prv_best;      /* temp subtotaling vars */
  float sc_best;                /* final best scores */
  float sc_M, sc_I, sc_D, sc_E; /* match, insert, delete, end scores */

  /* vars for sparse matrix */
  MATRIX_3D_SPARSE* mx;
  EDGEBOUNDS* edg_inner; /* edgebounds for search space of backward/forward */
  EDGEBOUNDS* edg_outer; /* edgebounds for sparse matrix shape */
  RANGE T_range;
  RANGE Q_range;
  bool is_q_0_in_dom_range;
  bool is_q_1_in_dom_range;

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
    MATRIX_2D_Fill(sp_MX, -INF);

    num_writes = 0;
    num_clears = 0;
  }
#endif

  /* --------------------------------------------------------------------------------- */

  // MATRIX_3D_SPARSE_Fill(st_SMX_opt, 0.0f);
  // MATRIX_2D_Fill(sp_MX_opt, 0.0f);

  /* initialize logsum lookup table if it has not already been */
  MATH_Logsum_Init();

  /* query sequence */
  seq = query->seq;
  N = EDGEBOUNDS_GetSize(edg);
  /* local or global alignments? */
  is_local = target->isLocal;
  sc_E = (is_local) ? 1.0 : 0.0;

  /* domain range (query sequence) */
  if (dom_range == NULL) {
    Q_range.beg = 0;
    Q_range.end = Q + 1;
  } else {
    Q_range = *dom_range;
  }
  /* target range */
  T_range.beg = 0;
  T_range.end = T + 1;

  /* init index */
  q_0 = 0;
  r_0b = 0;
  r_0e = 0;

  /* UNROLLED INITIAL ROW */
  q_0 = Q_range.beg;
  {
    /* check if query position is in domain */
    is_q_0_in_dom_range = IS_IN_RANGE(Q_range.beg, Q_range.end, q_0);
    /* get edgebound range */
    r_0b = EDGEBOUNDS_GetIndex_byRow_Fwd(edg, q_0);
    r_0e = EDGEBOUNDS_GetIndex_byRow_Fwd(edg, q_0 + 1);

    /* initialize special states */
    XMX_X(sp_MX_opt, SP_E, q_0) = MY_Zero();
    XMX_X(sp_MX_opt, SP_J, q_0) = MY_Zero();
    XMX_X(sp_MX_opt, SP_C, q_0) = MY_Zero();
    /* S->N, p=1 */
    XMX_X(sp_MX_opt, SP_N, q_0) = MY_One();
    /* S->N->B, no N-tail */
    XMX_X(sp_MX_opt, SP_B, q_0) = MY_One();

    /* only compute if in domain range */
    if (is_q_0_in_dom_range == true) {
      /* FOR every BOUND in zero row (-INF values set during initialization, so unneccessary) */
      for (r_0 = r_0b; r_0 < r_0e; r_0++) {
        /* get bound data */
        bnd = EDG_X(edg, r_0);
        id_0 = bnd.id;
        lb_0 = MAX(bnd.lb - 1, T_range.beg); /* can't overflow the left edge */
        rb_0 = MIN(bnd.rb, T_range.end);     /* can't overflow the right edge */

        /* fetch data mapping bound start location to data block in sparse matrix */
        qx0 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(st_SMX, r_0);
        /* initial location for square matrix and mapping to sparse matrix */
        t_0 = lb_0;
        tx0 = t_0 - bnd.lb; /* total_offset = offset_location - starting_location */

        /* FOR every position in TARGET profile */
        for (t_0 = lb_0; t_0 < rb_0; t_0++, tx0++) {
          tx0 = t_0 - bnd.lb;
          /* zero column is -inf in logspace.  We can skip this step and convert to normal space now. */
          MSMX_X(st_SMX_opt, qx0, tx0) = -INF;
          ISMX_X(st_SMX_opt, qx0, tx0) = -INF;
          DSMX_X(st_SMX_opt, qx0, tx0) = -INF;

/* embed linear row into quadratic test matrix */
#if DEBUG
          {
            MX_2D(cloud_MX, q_0, t_0) += 2.0;
            MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
            MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
            MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
          }
#endif
        }
      }

      /* init lookback 1 row */
      r_1b = r_0b;
      r_1e = r_0e;
    }
  }

  /* MAIN RECURSION */
  /* FOR every position in QUERY sequence (row in matrix) */
  for (q_0 = Q_range.beg + 1; q_0 < Q_range.end; q_0++) {
    q_1 = q_0 - 1;
    t_0 = 0;

    /* check if query position is in domain */
    is_q_0_in_dom_range = IS_IN_RANGE(Q_range.beg, Q_range.end, q_0);
    /* get edgebound range */
    r_0b = EDGEBOUNDS_GetIndex_byRow_Fwd(edg, q_0);
    r_0e = EDGEBOUNDS_GetIndex_byRow_Fwd(edg, q_0 + 1);

    /* Get next sequence character */
    a = seq[q_1]; /* off-by-one */
    A = AA_REV[a];

    XMX(SP_E, q_0) = -INF;

    /* only compute if in domain range */
    if (is_q_0_in_dom_range == true) {
      /* FOR every BOUND in current ROW */
      for (r_0 = r_0b; r_0 < r_0e; r_0++) {
        /* get bound data */
        bnd = EDG_X(edg, r_0);
        id_0 = bnd.id;
        lb_0 = MAX(bnd.lb - 1, T_range.beg); /* can't overflow left edge */
        rb_0 = MIN(bnd.rb, T_range.end);     /* can't overflow right edge */

        /* fetch data mapping bound start location to data block in sparse matrix */
        qx0 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(st_SMX, r_0);
        qx1 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Prv(st_SMX, r_0);

        /* initial location for square matrix and mapping to sparse matrix */
        t_0 = lb_0;
        tx0 = t_0 - bnd.lb; /* total_offset = offset_location - starting_location */
        tx1 = tx0 - 1;

        /* unrolled first loop: special case for left edge of range */
        t_0 = lb_0;
        {
          tx0 = t_0 - bnd.lb;

          /* zero column is -inf in logspace.  We can skip this step and convert to normal space now. */
          MSMX_X(st_SMX_opt, qx0, tx0) = -INF;
          ISMX_X(st_SMX_opt, qx0, tx0) = -INF;
          DSMX_X(st_SMX_opt, qx0, tx0) = -INF;

/* embed linear row into quadratic test matrix */
#if DEBUG
          {
            if (debugger->is_viz == true) {
              MX_2D(cloud_MX, q_0, t_0) = 1.0;
              MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
              MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
              MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
            }
          }
#endif
        }

        /* MAIN RECURSION */
        /* FOR every position in TARGET profile */
        for (t_0 = lb_0 + 1; t_0 < rb_0 - 1; t_0++) {
          t_1 = t_0 - 1;
          tx0 = t_0 - bnd.lb;
          tx1 = tx0 - 1;

          /* FIND MAX OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
          /* best previous state transition (match takes the diag element of each prev state) */
          prv_M = TSC_DELTA(t_1, M2M) *
                  (MSMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
          prv_I = TSC_DELTA(t_1, I2M) *
                  (ISMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
          prv_D = TSC_DELTA(t_1, D2M) *
                  (DSMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
          prv_B = TSC_DELTA(t_1, B2M) *
                  (XMX_X(sp_MX_opt, SP_B, q_1) + MSMX_X(st_SMX_post, qx0, tx0));
          /* best-to-match */
          prv_max = MAX(MAX(prv_M, prv_I),
                        MAX(prv_B, prv_D));
          MSMX_X(st_SMX_opt, qx0, tx0) = prv_max;

          /* FIND MAX OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
          /* previous states (match takes the previous row (upper) of each state) */
          prv_M = TSC_DELTA(t_0, M2I) *
                  (MSMX_X(st_SMX_opt, qx1, tx0) + ISMX_X(st_SMX_post, qx0, tx0));
          prv_I = TSC_DELTA(t_0, I2I) *
                  (ISMX_X(st_SMX_opt, qx1, tx0) + ISMX_X(st_SMX_post, qx0, tx0));
          /* best-to-insert */
          prv_max = MAX(prv_M, prv_I);
          ISMX_X(st_SMX_opt, qx0, tx0) = prv_max;

          /* FIND MAX OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
          /* previous states (match takes the previous column (left) of each state) */
          prv_M = TSC_DELTA(t_1, M2D) *
                  (MSMX_X(st_SMX_opt, qx0, tx1));
          prv_D = TSC_DELTA(t_1, D2D) *
                  (DSMX_X(st_SMX_opt, qx0, tx1));
          /* best-to-delete */
          prv_max = MAX(prv_M, prv_D);
          DSMX_X(st_SMX_opt, qx0, tx0) = prv_max;

          /* UPDATE E STATE */
          prv_M = MSMX(qx0, tx0) * sc_E;
          prv_E = XMX_X(sp_MX_opt, SP_E, q_0);
          /* best-to-end */
          prv_max = MAX(prv_E, prv_M);
          XMX_X(sp_MX_opt, SP_E, q_0) = prv_max;

/* embed linear row into quadratic test matrix */
#if DEBUG
          {
            if (debugger->is_viz == true) {
              MX_2D(cloud_MX, q_0, t_0) = 1.0;
              MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
              MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
              MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
            }
          }
#endif
        }

        /* UNROLLED FINAL LOOP ITERATION */
        if (rb_0 > 1) {
          t_0 = rb_0 - 1;
          t_1 = t_0 - 1;
          tx0 = t_0 - bnd.lb;
          tx1 = tx0 - 1;

          /* FIND MAX OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
          /* best previous state transition (match takes the diag element of each prev state) */
          prv_M = TSC_DELTA(t_1, M2M) *
                  (MSMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
          prv_I = TSC_DELTA(t_1, I2M) *
                  (ISMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
          prv_D = TSC_DELTA(t_1, D2M) *
                  (DSMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
          prv_B = TSC_DELTA(t_1, B2M) *
                  (XMX_X(sp_MX_opt, SP_B, q_1) + MSMX_X(st_SMX_post, qx0, tx0));
          /* best-to-match */
          prv_max = MAX(MAX(prv_M, prv_I),
                        MAX(prv_B, prv_D));
          MSMX_X(st_SMX_opt, qx0, tx0) = prv_max;

          /* FIND SUM OF PATHS TO INSERT STATE (unrolled) */
          ISMX(qx0, tx0) = -INF;

          /* FIND MAX OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
          /* previous states (match takes the previous column (left) of each state) */
          prv_M = TSC_DELTA(t_1, M2D) *
                  (MSMX_X(st_SMX_opt, qx0, tx1));
          prv_D = TSC_DELTA(t_1, D2D) *
                  (DSMX_X(st_SMX_opt, qx0, tx1));
          /* best-to-delete */
          prv_max = MAX(prv_M, prv_D);
          DSMX_X(st_SMX_opt, qx0, tx0) = prv_max;

          /* UPDATE E STATE */
          prv_M = MSMX(qx0, tx0) * sc_E;
          prv_E = XMX_X(sp_MX_opt, SP_E, q_0);
          /* best-to-end */
          prv_max = MAX(prv_E, prv_M);
          XMX_X(sp_MX_opt, SP_E, q_0) = prv_max;

/* embed linear row into quadratic test matrix */
#if DEBUG
          {
            if (debugger->is_viz == true) {
              MX_2D(cloud_MX, q_0, t_0) = 1.0;
              MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
              MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
              MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
            }
          }
#endif
        }
      }
    }

    /* SPECIAL STATES */
    t_0 = T;

    prv_E = XMX_X(sp_MX_opt, SP_E, q_0);
    prv_M = MSMX_X(st_SMX_post, qx0, tx0);
    prv_D = DSMX_X(st_SMX_post, qx0, tx0);
    prv_max = MAX(prv_E,
                  MAX(prv_M, prv_D));

    /* J state ( J->J or E->J (E's loop) ) */
    prv_J = XSC_DELTA(SP_J, SP_LOOP) *
            (XMX_X(sp_MX_opt, SP_J, q_1) + XMX_X(sp_MX_post, SP_J, q_0));
    prv_E = XSC_DELTA(SP_E, SP_LOOP) *
            (XMX_X(sp_MX_opt, SP_E, q_0));
    prv_max = MAX(prv_J, prv_E);
    XMX_X(sp_MX_opt, SP_J, q_0) = prv_max;

    /* C state ( C->C or E->C ) */
    prv_C = XSC_DELTA(SP_C, SP_LOOP) *
            (XMX_X(sp_MX_opt, SP_C, q_1) + XMX_X(sp_MX_post, SP_C, q_0));
    prv_E = XSC_DELTA(SP_E, SP_MOVE) *
            (XMX_X(sp_MX_opt, SP_E, q_0));
    prv_max = MAX(prv_C, prv_E);
    XMX(SP_C, q_0) = prv_max;

    /* N state ( N->N (tail) ) */
    prv_N = XSC_DELTA(SP_N, SP_LOOP) *
            (XMX_X(sp_MX_opt, SP_N, q_1) + XMX_X(sp_MX_post, SP_N, q_0));
    XMX(SP_N, q_0) = prv_N;

    /* B state ( N->B (N's move) or J->B (J's move) ) */
    prv_N = XSC_DELTA(SP_N, SP_MOVE) *
            (XMX_X(sp_MX_opt, SP_N, q_0));
    prv_J = XSC_DELTA(SP_J, SP_MOVE) *
            (XMX_X(sp_MX_opt, SP_J, q_0));
    prv_max = MAX(prv_N, prv_J);
    XMX(SP_B, q_0) = prv_max;

    /* SET CURRENT ROW TO PREVIOUS ROW */
    r_1b = r_0b;
    r_1e = r_0e;
  }

  /* T state */
  sc_best = XMX_X(sp_MX_post, SP_C, Q);
  *sc_final = sc_best;

  return STATUS_SUCCESS;
}

/* MATH RULES: These determine how probilities are summed, multiplied, and certain identities */

static inline float
MY_Sum(const float x,
       const float y) {
  return MATH_LogSum(x, y);
}

static inline float
MY_Prod(const float x,
        const float y) {
  return MATH_LogProd(x, y);
}

static inline float
MY_Max(const float x,
       const float y) {
  return MATH_Max(x, y);
}

static inline float
MY_Zero() {
  return MATH_LogZero();
}

static inline float
MY_One() {
  return MATH_LogOne();
}
