/*******************************************************************************
 *  - FILE:  bound_fwdbck_sparse.c
 *  - DESC:  Bounded Forward/Backward Algorithm
 *           (Sparse Space)
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

/* self header */
#include "_algs_sparse.h"
#include "bound_fwdbck_sparse.h"

/*!  NOTE:  HOW TO CONVERT row-coords to diag-coords
 *          MMX3(i-1,j-1) => MMX3(, d_2)
 *          MMX3(i,  j-1) => MMX3(, d_1)
 *          MMX3(i,  j  ) => MMX3(, d_1)
 */

/* private functions */
static inline float
MY_Sum(const float x, const float y);

static inline float
MY_Prod(const float x, const float y);

static inline float
MY_Zero();

static inline float
MY_One();

STATUS_FLAG
run_Bound_Forward_Sparse(const SEQUENCE* query,                 /* query sequence */
                         const HMM_PROFILE* target,             /* target HMM model */
                         const int Q,                           /* query length */
                         const int T,                           /* target length */
                         MATRIX_3D_SPARSE* restrict st_SMX_fwd, /* normal state matrix */
                         MATRIX_2D* restrict sp_MX_fwd,         /* special state matrix */
                         const EDGEBOUNDS* edg,                 /* edgebounds */
                         // const DP_OPTS*                options,       /* dynamic programming options */
                         const RANGE* dom_range, /* (OPTIONAL) domain range for computing fwd/bck on specific domain. If NULL, computes complete fwd/bck. */
                         float* sc_final)        /* (OUTPUT) final score */
{
  /* vars for matrix access for macros */
  MATRIX_3D_SPARSE* st_SMX = st_SMX_fwd; /* normal state matrix */
  MATRIX_2D* sp_MX = sp_MX_fwd;          /* special state matrix */

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
  BOUND bnd;      /* current bound */
  int id_0;       /* id in edgebound list (row/diag) */
  int r_0;        /* current index for current row */
  int r_0b, r_0e; /* begin and end indices for current row in edgebound list */
  int r_1;        /* current index for previous row */
  int r_1b, r_1e; /* begin and end indices for current row in edgebound list */
  int le_0, re_0; /* right/left matrix bounds of current diag */
  int lb_0, rb_0; /* bounds of current search space on current diag */
  int lb_1, rb_1; /* bounds of current search space on previous diag */
  int lb_2, rb_2; /* bounds of current search space on 2-back diag */
  int lb_T, rb_T; /* truncated bounds within domain */

  /* vars for recurrance scores */
  float prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
  float prv_B, prv_E;           /* previous (B) begin and (E) end states */
  float prv_J, prv_N, prv_C;    /* previous (J) jump, (N) initial, and (C) terminal states */
  float prv_loop, prv_move;     /* previous loop and move for special states */
  float prv_sum, prv_best;      /* temp subtotaling vars */
  float sc_best;                /* final best scores */
  float sc_M, sc_I, sc_D, sc_E; /* match, insert, delete, end scores */

  /* vars for sparse matrix */
  MATRIX_3D_SPARSE* mx;
  EDGEBOUNDS* edg_inner;    /* edgebounds for search space of backward/forward */
  EDGEBOUNDS* edg_outer;    /* edgebounds for sparse matrix shape */
  RANGE T_range;            /* target range */
  RANGE Q_range;            /* query range */
  bool is_q_0_in_dom_range; /* checks if current query position is inside the domain range */
  bool is_q_1_in_dom_range; /* checks if previous query position is inside the domain range */

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
    MATRIX_2D_Fill(cloud_MX, 0.0);
    int num_cells = MATRIX_2D_Cloud_Fill(cloud_MX, edg, 1.0);
    int num_cells2 = MATRIX_2D_Cloud_Count(cloud_MX);
    float pretotal = MATRIX_2D_Total(cloud_MX);
    // printf("Pre-Forward Cell Totals = %d %d %f\n", num_cells, num_cells2, pretotal);

    // MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
    // MATRIX_2D_Fill( cloud_MX3, 0 );
    // MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
    // MATRIX_3D_Fill( test_MX, -INF );
    // MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
    // MATRIX_3D_Fill( test_MX3, -INF );
    // EDGEBOUNDS_Reuse( test_edg, Q, T );
    // MATRIX_2D_Fill( sp_MX, -INF );

    num_writes = 0;
    num_clears = 0;
  }
#endif

  /* --------------------------------------------------------------------------------- */

  /* initialize logsum lookup table if it has not already been */
  MATH_Logsum_Init();

  /* query sequence */
  mx = st_SMX;
  seq = query->seq;
  N = EDGEBOUNDS_GetSize(edg);
  /* local or global alignments? */
  is_local = target->isLocal;
  sc_E = (is_local) ? MY_One() : MY_Zero();
  if (is_local == true) {
    MY_One();
  } else {
    MY_Zero();
  }

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

  /* UNROLLED INITIAL QUERY LOOP */
  q_0 = Q_range.beg;
  {
    /* check if query position is in domain */
    is_q_0_in_dom_range = IS_IN_RANGE(Q_range.beg, Q_range.end, q_0);
    /* get edgebound range */
    r_0b = EDGEBOUNDS_GetIndex_byRow_Fwd(edg, q_0);
    r_0e = EDGEBOUNDS_GetIndex_byRow_Fwd(edg, q_0 + 1);

    /* initialize special states */
    XMX(SP_E, q_0) = MY_Zero();
    XMX(SP_J, q_0) = MY_Zero();
    XMX(SP_C, q_0) = MY_Zero();
    /* S->N, p=1 */
    XMX(SP_N, q_0) = MY_One();
    /* S->N->B, no N-tail */
    XMX(SP_B, q_0) = XSC(SP_N, SP_MOVE);

    /* only compute if in domain range */
    if (is_q_0_in_dom_range == true) {
      /* FOR every BOUND in zero row */
      for (r_0 = r_0b; r_0 < r_0e; r_0++) {
        /* get bound data */
        bnd = MATRIX_3D_SPARSE_GetBound_byIndex(st_SMX, r_0);
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
          MSMX(qx0, tx0) = MY_Zero();
          ISMX(qx0, tx0) = MY_Zero();
          DSMX(qx0, tx0) = MY_Zero();

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
    }

    /* lookback one row */
    r_1b = r_0b;
    r_1e = r_0e;
  }

  /* MAIN QUERY LOOP */
  /* FOR every position in QUERY sequence (row in matrix) */
  for (q_0 = Q_range.beg + 1; q_0 < Q_range.end; q_0++) {
    q_1 = q_0 - 1;
    t_0 = 0;

    /* check if query position is in domain */
    is_q_0_in_dom_range = IS_IN_RANGE(Q_range.beg + 1, Q_range.end, q_0);
    // is_q_0_in_dom_range = (q_0 >= Q_range.beg && q_0 < Q_range.end);
    /* get edgebound range */
    r_0b = EDGEBOUNDS_GetIndex_byRow_Fwd(edg, q_0);
    r_0e = EDGEBOUNDS_GetIndex_byRow_Fwd(edg, q_0 + 1);

    /* Get next sequence character */
    a = seq[q_1]; /* off-by-one */
    A = AA_REV[a];

    /* Init E state for current row */
    XMX(SP_E, q_0) = MY_Zero();

    /* only compute if in domain range */
    if (is_q_0_in_dom_range == true) {
      /* FOR every BOUND in current ROW */
      for (r_0 = r_0b; r_0 < r_0e; r_0++) {
        /* get bound data */
        bnd = MATRIX_3D_SPARSE_GetBound_byIndex(st_SMX, r_0);
        id_0 = bnd.id;
        lb_0 = MAX(bnd.lb - 1, T_range.beg); /* can't overflow left edge. the leftmost cell will be set to zero, so (-1) adds left padding cell.  */
        rb_0 = MIN(bnd.rb, T_range.end);     /* can't overflow right edge */

        /* fetch data mapping bound start location to data block in sparse matrix */
        qx0 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(st_SMX, r_0);
        qx1 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Prv(st_SMX, r_0);
        /* initial location for square matrix and mapping to sparse matrix */
        t_0 = lb_0;
        tx0 = t_0 - bnd.lb; /* total_offset = offset_location - starting_location */
        tx1 = tx0 - 1;

        /* UNROLLED INITIAL TARGET LOOP: special case for left edge of range */
        t_0 = lb_0;
        {
          tx0 = t_0 - bnd.lb;

          MSMX(qx0, tx0) = MY_Zero();
          ISMX(qx0, tx0) = MY_Zero();
          DSMX(qx0, tx0) = MY_Zero();

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

        /* MAIN TARGET LOOP */
        /* FOR every position in TARGET profile */
        for (t_0 = lb_0 + 1; t_0 < rb_0 - 1; t_0++) {
          t_1 = t_0 - 1;
          tx0 = t_0 - bnd.lb;
          tx1 = tx0 - 1;

          /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
          /* best previous state transition (match takes the diag element of each prev state) */
          prv_M = MY_Prod(MSMX(qx1, tx1), TSC(t_1, M2M));
          prv_I = MY_Prod(ISMX(qx1, tx1), TSC(t_1, I2M));
          prv_D = MY_Prod(DSMX(qx1, tx1), TSC(t_1, D2M));
          prv_B = MY_Prod(XMX(SP_B, q_1), TSC(t_1, B2M)); /* from begin match state (new alignment) */
          /* best-to-match */
          prv_sum = MY_Sum(MY_Sum(prv_M, prv_I),
                           MY_Sum(prv_B, prv_D));
          MSMX(qx0, tx0) = MY_Prod(prv_sum, MSC(t_0, A));

          /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
          /* previous states (match takes the previous row (upper) of each state) */
          prv_M = MY_Prod(MSMX(qx1, tx0), TSC(t_0, M2I));
          prv_I = MY_Prod(ISMX(qx1, tx0), TSC(t_0, I2I));
          /* best-to-insert */
          prv_sum = MY_Sum(prv_M, prv_I);
          ISMX(qx0, tx0) = MY_Prod(prv_sum, ISC(t_0, A));

          /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
          /* previous states (match takes the previous column (left) of each state) */
          prv_M = MY_Prod(MSMX(qx0, tx1), TSC(t_1, M2D));
          prv_D = MY_Prod(DSMX(qx0, tx1), TSC(t_1, D2D));
          /* best-to-delete */
          prv_sum = MY_Sum(prv_M, prv_D);
          DSMX(qx0, tx0) = prv_sum;

          /* UPDATE E STATE */
          prv_M = MY_Prod(MSMX(qx0, tx0), sc_E);
          prv_D = MY_Prod(DSMX(qx0, tx0), sc_E);
          prv_E = XMX(SP_E, q_0);
          /* best-to-e-state */
          prv_sum = MY_Sum(MY_Sum(prv_M, prv_D),
                           prv_E);
          XMX(SP_E, q_0) = prv_sum;

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

        /* UNROLLED FINAL TARGET LOOP: special case for right edge of range (only when range is greater than one cell) */
        if (t_0 = rb_0 - 1, (rb_0 - lb_0) > 1) {
          t_1 = t_0 - 1;
          tx0 = t_0 - bnd.lb;
          tx1 = tx0 - 1;

          /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
          /* best previous state transition (match takes the diag element of each prev state) */
          prv_M = MY_Prod(MSMX(qx1, tx1), TSC(t_1, M2M));
          prv_I = MY_Prod(ISMX(qx1, tx1), TSC(t_1, I2M));
          prv_D = MY_Prod(DSMX(qx1, tx1), TSC(t_1, D2M));
          prv_B = MY_Prod(XMX(SP_B, q_1), TSC(t_1, B2M)); /* from begin match state (new alignment) */
          /* sum-to-match */
          prv_sum = MY_Sum(MY_Sum(prv_M, prv_I),
                           MY_Sum(prv_D, prv_B));
          MSMX(qx0, tx0) = MY_Prod(prv_sum, MSC(t_0, A));

          /* FIND SUM OF PATHS TO INSERT STATE (unrolled) */
          ISMX(qx0, tx0) = MY_Zero();

          /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) (unrolled) */
          /* previous states (match takes the left element of each state) */
          prv_M = MY_Prod(MSMX(qx0, tx1), TSC(t_1, M2D));
          prv_D = MY_Prod(DSMX(qx0, tx1), TSC(t_1, D2D));
          /* sum-to-delete */
          prv_sum = MY_Sum(prv_M, prv_D);
          DSMX(qx0, tx0) = prv_sum;

          /* UPDATE E STATE (unrolled) */
          prv_E = XMX(SP_E, q_0);
          prv_M = MSMX(qx0, tx0);
          prv_D = DSMX(qx0, tx0);
          /* best-to-begin */
          prv_sum = MY_Sum(MY_Sum(prv_D, prv_M),
                           prv_E);
          XMX(SP_E, q_0) = prv_sum;

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
    }

    /* SPECIAL STATES */
    /* J state */
    prv_J = MY_Prod(XMX(SP_J, q_1), XSC(SP_J, SP_LOOP)); /* J->J */
    prv_E = MY_Prod(XMX(SP_E, q_0), XSC(SP_E, SP_LOOP)); /* E->J is E's "loop" */
    prv_sum = MY_Sum(prv_J, prv_E);
    XMX(SP_J, q_0) = prv_sum;

    /* C state */
    prv_C = MY_Prod(XMX(SP_C, q_1), XSC(SP_C, SP_LOOP));
    prv_E = MY_Prod(XMX(SP_E, q_0), XSC(SP_E, SP_MOVE));
    prv_sum = MY_Sum(prv_C, prv_E);
    XMX(SP_C, q_0) = prv_sum;

    /* N state */
    prv_N = MY_Prod(XMX(SP_N, q_1), XSC(SP_N, SP_LOOP));
    XMX(SP_N, q_0) = prv_N;

    /* B state */
    prv_N = MY_Prod(XMX(SP_N, q_0), XSC(SP_N, SP_MOVE)); /* N->B is N's move */
    prv_J = MY_Prod(XMX(SP_J, q_0), XSC(SP_J, SP_MOVE)); /* J->B is J's move */
    prv_sum = MY_Sum(prv_N, prv_J);
    XMX(SP_B, q_0) = prv_sum;

    /* SET CURRENT ROW TO PREVIOUS ROW */
    r_1b = r_0b;
    r_1e = r_0e;
  }

  /* T state */
  sc_best = MY_Prod(XMX(SP_C, Q_range.end - 1), XSC(SP_C, SP_MOVE));
  *sc_final = sc_best;

/* output test matrix */
#if DEBUG
  {
  }
#endif

  return STATUS_SUCCESS;
}

/** FUNCTION:  run_Bound_Backward_Sparse()
 *  SYNOPSIS:  Perform Edge-Bounded Backward step of Cloud Search Algorithm.
 *             Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of
 *             (left-bound, right-bound) pairs sorted by row.
 *             Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_().
 *             Final score produced by Forward is stored in <sc_final>.
 *
 *    RETURN:  Returns the final score of the Forward Algorithm.
 */
STATUS_FLAG
run_Bound_Backward_Sparse(const SEQUENCE* query,                 /* query sequence */
                          const HMM_PROFILE* target,             /* target HMM model */
                          const int Q,                           /* query length */
                          const int T,                           /* target length */
                          MATRIX_3D_SPARSE* restrict st_SMX_bck, /* normal state matrix */
                          MATRIX_2D* restrict sp_MX_bck,         /* special state matrix */
                          const EDGEBOUNDS* edg,                 /* edgebounds */
                          const RANGE* dom_range,                /* (OPTIONAL) domain range for computing fwd/bck on specific domain. If NULL, computes complete fwd/bck. */
                          float* sc_final)                       /* (OUTPUT) final score */
{
  /* vars for matrix access for macros */
  MATRIX_3D_SPARSE* st_SMX = st_SMX_bck; /* normal state matrix */
  MATRIX_2D* sp_MX = sp_MX_bck;          /* special state matrix */

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

  /* vars for indexing into data matrices by anti-diag */
  int d_0, d_1, d_2;         /* real index of current and previous antidiagonals */
  int dx0, dx1, dx2;         /* mod mapping of antidiagonal index into data matrix */
  int k_0, k_1;              /* offset into antidiagonal */
  int d_st, d_end, d_cnt;    /* starting and ending diagonal indices */
  int dim_T, dim_Q, dim_TOT; /* dimensions of submatrix being searched */
  int dim_min, dim_max;      /* diagonal index where num cells reaches highest point and diminishing point */
  int num_cells;             /* number of cells in current diagonal */

  /* vars for indexing into edgebound lists */
  BOUND bnd;      /* current bound */
  int id;         /* id in edgebound list (row/diag) */
  int r_0;        /* current index for current row */
  int r_0b, r_0e; /* begin and end indices for current row in edgebound list */
  int r_1;        /* current index for previous row */
  int r_1b, r_1e; /* begin and end indices for current row in edgebound list */
  int le_0, re_0; /* right/left matrix bounds of current diag */
  int lb_0, rb_0; /* bounds of current search space on current diag */
  int lb_1, rb_1; /* bounds of current search space on previous diag */
  int lb_2, rb_2; /* bounds of current search space on 2-back diag */
  int lb_T, rb_T; /* truncated bounds within domain */

  /* vars for recurrance scores */
  float prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
  float prv_B, prv_E;           /* previous (B) begin and (E) end states */
  float prv_J, prv_N, prv_C;    /* previous (J) jump, (N) initial, and (C) terminal states */
  float prv_loop, prv_move;     /* previous loop and move for special states */
  float prv_sum, prv_best;      /* temp subtotaling vars */
  float sc_best;                /* final best scores */
  float sc_M, sc_I, sc_D, sc_E; /* match, insert, delete, end scores */

  /* vars for sparse matrix */
  MATRIX_3D_SPARSE* mx;
  EDGEBOUNDS* edg_inner;    /* edgebounds for search space of backward/forward */
  EDGEBOUNDS* edg_outer;    /* edgebounds for sparse matrix shape */
  RANGE T_range;            /* target range */
  RANGE Q_range;            /* query range */
  bool is_q_0_in_dom_range; /* checks if current query position is inside the domain range */
  bool is_q_1_in_dom_range; /* checks if previous query position is inside the domain range */

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
  mx = st_SMX;
  seq = query->seq;
  N = EDGEBOUNDS_GetSize(edg);
  /* local or global alignments? */
  is_local = target->isLocal;
  sc_E = (is_local) ? MY_One() : MY_Zero();
  if (is_local == true) {
    MY_One();
  } else {
    MY_Zero();
  }

  /* domain range (query sequence) */
  if (dom_range == NULL) {
    Q_range.beg = 0;
    Q_range.end = Q;
  } else {
    Q_range = *dom_range;
  }
  /* valid target range */
  T_range.beg = 1;
  T_range.end = T;

  /* UNROLLED INITIAL QUERY LOOP */
  q_0 = Q_range.end;
  {
    /* if inside domain */
    is_q_0_in_dom_range = IS_IN_RANGE(Q_range.beg, Q_range.end, q_0);
    /* get edgebound range */
    r_0b = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0 + 1);
    r_0e = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0);

    /* INIT SPECIAL STATES */
    XMX(SP_J, q_0) = MY_Zero();
    XMX(SP_B, q_0) = MY_Zero();
    XMX(SP_N, q_0) = MY_Zero();
    XMX(SP_C, q_0) = XSC(SP_C, SP_MOVE);
    XMX(SP_E, q_0) = MY_Prod(XMX(SP_C, q_0), XSC(SP_E, SP_MOVE));

    /* if sequence position is in domain range */
    if (is_q_0_in_dom_range == true) {
      bnd = MATRIX_3D_SPARSE_GetBound_byIndex(st_SMX, r_0b);

      /* FOR every SPAN in current ROW */
      for (r_0 = r_0b; r_0 > r_0e; r_0--) {
        /* get bound data */
        // bnd   = EDG_X(edg, r_0);               /* bounds for current bound */
        bnd = MATRIX_3D_SPARSE_GetBound_byIndex(st_SMX, r_0);
        lb_0 = MAX(bnd.lb, T_range.beg);     /* can't overflow left edge */
        rb_0 = MIN(bnd.rb, T_range.end + 1); /* can't overflow right edge */

        /* fetch data mapping bound start location to data block in sparse matrix */
        qx0 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(st_SMX, r_0);
        qx1 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Nxt(st_SMX, r_0);
        /* location for square matrix and mapping to sparse matrix */
        t_0 = lb_0;
        tx0 = t_0 - bnd.lb;

        /* UNROLLED INITIAL TARGET LOOP */
        t_0 = rb_0 - 1;
        {
          t_1 = t_0 + 1;
          tx0 = t_0 - bnd.lb;
          tx1 = tx0 + 1;

          MSMX(qx0, tx0) = XMX(SP_E, q_0);
          ISMX(qx0, tx0) = MY_Zero();
          DSMX(qx0, tx0) = XMX(SP_E, q_0);

#if DEBUG
          {
            MX_2D(cloud_MX, q_0, t_0) += 2.0;
            MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
            MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
            MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
          }
#endif
        }

        /* MAIN TARGET LOOP */
        /* FOR every position of TARGET in SPAN */
        for (t_0 = rb_0 - 2; t_0 >= lb_0; t_0--) {
          t_1 = t_0 + 1;
          tx0 = t_0 - bnd.lb;
          tx1 = tx0 + 1;

          prv_E = MY_Prod(XMX(SP_E, Q), sc_E);
          prv_D = MY_Prod(DSMX(qx0, tx1), TSC(t_0, M2D));
          MSMX(qx0, tx0) = MY_Sum(prv_E, prv_D);

          ISMX(qx0, tx0) = MY_Zero();

          prv_E = MY_Prod(XMX(SP_E, Q), sc_E);
          prv_D = MY_Prod(DSMX(qx0, tx1), TSC(t_0, D2D));
          DSMX(qx0, tx0) = MY_Sum(prv_E, prv_D);

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
    }

    /* init lookback 1 row */
    r_1b = r_0b;
    r_1e = r_0e;
  }

  /* MAIN QUERY LOOP */
  /* FOR every position in QUERY */
  for (q_0 = Q_range.end - 1; q_0 > Q_range.beg; q_0--) {
    q_1 = q_0 + 1;

    /* if inside domain */
    is_q_0_in_dom_range = IS_IN_RANGE(Q_range.beg, Q_range.end, q_0);
    is_q_1_in_dom_range = IS_IN_RANGE(Q_range.beg, Q_range.end, q_1);
    /* get edgebound range */
    r_0b = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0 + 1);
    r_0e = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0);

    /* Get next sequence character */
    a = seq[q_0];
    A = AA_REV[a];

    /* UPDATE B STATE */
    XMX(SP_B, q_0) = MY_Zero();
    /* if previous q is in domain range, update B state */
    if (is_q_1_in_dom_range == true) {
      for (r_1 = r_1b; r_1 > r_1e; r_1--) {
        /* get bound data */
        // bnd   = EDG_X(edg, r_1);               /* bounds for current bound */
        bnd = MATRIX_3D_SPARSE_GetBound_byIndex(st_SMX, r_1);
        lb_0 = MAX(bnd.lb, T_range.beg);     /* can't overflow left edge */
        rb_0 = MIN(bnd.rb, T_range.end + 1); /* can't overflow right edge */
        // printf("q_0=%d, r_0=%d, b_0=(%d,%d)\n", q_0, r_1, lb_0, rb_0);

        /* fetch data location to bound start location (in offset) */
        // qx1 = VECTOR_INT_Get( st_SMX->imap_cur, r_1 );    /* (q_0, t_0) location offset */
        qx1 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(st_SMX, r_1);
        /* location for square matrix and mapping to sparse matrix */
        t_0 = lb_0;
        tx0 = t_0 - bnd.lb; /* total_offset = offset_location - starting_location */

        for (t_0 = rb_0 - 1; t_0 >= lb_0; t_0--) {
          t_1 = t_0 - 1;
          tx0 = t_0 - bnd.lb;
          tx1 = tx0 - 1;

          prv_sum = XMX(SP_B, q_0);
          prv_M = MY_Prod(MSMX(qx1, tx0),
                          MY_Prod(TSC(t_1, B2M), MSC(t_0, A)));
          XMX(SP_B, q_0) = MY_Sum(prv_sum, prv_M);
          // printf("[%d]%f %f %f\n", t_0, XMX(SP_B, q_0), prv_M, prv_sum );
        }
        // printf("bsc = %f\n", XMX(SP_B, q_0) );
      }
    }

    /* UPDATE SPECIAL STATES */
    prv_J = MY_Prod(XMX(SP_J, q_1), XSC(SP_J, SP_LOOP));
    prv_B = MY_Prod(XMX(SP_B, q_0), XSC(SP_J, SP_MOVE));
    XMX(SP_J, q_0) = MY_Sum(prv_J, prv_B);

    prv_C = MY_Prod(XMX(SP_C, q_1), XSC(SP_C, SP_LOOP));
    XMX(SP_C, q_0) = prv_C;

    prv_J = MY_Prod(XMX(SP_J, q_0), XSC(SP_E, SP_LOOP));
    prv_C = MY_Prod(XMX(SP_C, q_0), XSC(SP_E, SP_MOVE));
    XMX(SP_E, q_0) = MY_Sum(prv_J, prv_C);

    prv_N = MY_Prod(XMX(SP_N, q_1), XSC(SP_N, SP_LOOP));
    prv_B = MY_Prod(XMX(SP_B, q_0), XSC(SP_N, SP_MOVE));
    XMX(SP_N, q_0) = MY_Sum(prv_N, prv_B);

    if (is_q_0_in_dom_range == true) {
      /* FOR every SPAN in current ROW */
      for (r_0 = r_0b; r_0 > r_0e; r_0--) {
        /* get bound data */
        // bnd   = EDG_X(edg, r_0);       /* bounds for current bound */
        bnd = MATRIX_3D_SPARSE_GetBound_byIndex(st_SMX, r_0);
        lb_0 = MAX(bnd.lb, T_range.beg);     /* can't overflow left edge */
        rb_0 = MIN(bnd.rb, T_range.end + 1); /* can't overflow right edge */
        // printf("q_0=%d, r_0=%d, b_0=(%d,%d)\n", q_0, r_0, lb_0, rb_0);

        /* fetch data location to bound start location (in offset) */
        // qx0   = VECTOR_INT_Get( st_SMX->imap_cur, r_0 );    /* (q_0, t_0) location offset */
        // qx1   = VECTOR_INT_Get( st_SMX->imap_nxt, r_0 );    /* (q_0, t_0) location offset */
        qx0 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(st_SMX, r_0);
        qx1 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Nxt(st_SMX, r_0);
        /* location for square matrix and mapping to sparse matrix */
        t_0 = T;
        tx0 = t_0 - bnd.lb; /* total_offset = offset_location - starting_location */

        /* UNROLLED INITIAL TARGET LOOP */
        t_0 = rb_0 - 1;
        {
          t_1 = t_0 + 1;
          tx0 = t_0 - bnd.lb;
          tx1 = tx0 + 1;

          MSMX(qx0, tx0) = XMX(SP_E, q_0);
          ISMX(qx0, tx0) = MY_Zero();
          DSMX(qx0, tx0) = XMX(SP_E, q_0);

#if DEBUG
          {
            MX_2D(cloud_MX, q_0, t_0) += 2.0;
            MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
            MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
            MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
          }
#endif
        }

        /* MAIN TARGET LOOP */
        /* FOR every position of TARGET in SPAN */
        for (t_0 = rb_0 - 2; t_0 >= lb_0; t_0--) {
          t_1 = t_0 + 1;
          tx0 = t_0 - bnd.lb;
          tx1 = tx0 + 1;

          /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
          prv_M = MY_Prod(MSC(t_1, A),
                          MY_Prod(MSMX(qx1, tx1), TSC(t_0, M2M)));
          prv_I = MY_Prod(ISC(t_1, A),
                          MY_Prod(ISMX(qx1, tx0), TSC(t_0, M2I)));
          prv_D = MY_Prod(DSMX(qx0, tx1), TSC(t_0, M2D));
          prv_E = MY_Prod(XMX(SP_E, q_0), sc_E); /* from end match state (new alignment) */
          /* best-to-match */
          prv_sum = MY_Sum(MY_Sum(prv_M, prv_I),
                           MY_Sum(prv_E, prv_D));
          MSMX(qx0, tx0) = prv_sum;

          // printf("(q_0,t_0)=(%d,%d) => a=(%c,%d) tsc=%f %f %f msc=%f mmx=%f %f imx=%f dmx=%f %f emx=%f\n",
          //    q_0, t_0, a, A, TSC(t_1, M2M), TSC(t_0, M2I), TSC(t_0, M2D), MSC(t_0, A),
          //    MSMX(qx0, tx0), prv_M, prv_I, prv_D, DSMX(qx0, tx1), prv_E );

          /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
          prv_M = MY_Prod(MSMX(qx1, tx1),
                          MY_Prod(TSC(t_0, I2M), MSC(t_1, A)));
          prv_I = MY_Prod(ISMX(qx1, tx0),
                          MY_Prod(TSC(t_0, I2I), ISC(t_0, A)));
          /* best-to-insert */
          prv_sum = MY_Sum(prv_M, prv_I);
          ISMX(qx0, tx0) = prv_sum;

          /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
          prv_M = MY_Prod(MSMX(qx1, tx1),
                          MY_Prod(TSC(t_0, D2M), MSC(t_1, A)));
          prv_D = MY_Prod(DSMX(qx0, tx1), TSC(t_0, D2D));
          prv_E = MY_Prod(XMX(SP_E, q_0), sc_E);
          /* best-to-delete */
          prv_sum = MY_Sum(prv_M,
                           MY_Sum(prv_D, prv_E));
          DSMX(qx0, tx0) = prv_sum;
          // printf("dmx=%f => %f %f %f\n", prv_sum, prv_M, prv_D, prv_E);

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
    }

    /* SET CURRENT ROW TO PREVIOUS ROW */
    r_1b = r_0b;
    r_1e = r_0e;
  }

  /* UNROLLED FINAL ROW */
  q_0 = Q_range.beg;
  {
    /* At q_0 = 0, only N,B states are reachable. */
    q_1 = q_0 + 1;

    /* get edgebound range */
    r_0b = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0 + 1);
    r_0e = EDGEBOUNDS_GetIndex_byRow_Bck(edg, q_0);

    /* get sequence */
    a = seq[q_0];
    A = AA_REV[a];

    /* UPDATE B STATE */
    XMX(SP_B, q_0) = MY_Zero();

    /* if previous q is in domain, update B state */
    if (is_q_1_in_dom_range == true) {
      for (r_1 = r_1b; r_1 > r_1e; r_1--) {
        // bnd   = EDG_X(edg, r_1);       /* bounds for current bound */
        bnd = MATRIX_3D_SPARSE_GetBound_byIndex(st_SMX, r_1);
        lb_0 = MAX(bnd.lb, T_range.beg); /* can't overflow left edge */
        rb_0 = MIN(bnd.rb, T_range.end); /* can't overflow right edge */

        /* fetch data location to bound start location (in offset) */
        // qx1 = VECTOR_INT_Get( st_SMX->imap_cur, r_1 );    /* (q_0, t_0) location offset */
        qx1 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur(st_SMX, r_1);

        /* location for square matrix and mapping to sparse matrix */
        t_0 = lb_0;
        tx0 = t_0 - bnd.lb; /* total_offset = offset_location - starting_location */

        for (t_0 = rb_0 - 1; t_0 >= lb_0; t_0--) {
          t_1 = t_0 - 1;
          /* calculate offset from beginning of sparse data block */
          tx0 = t_0 - bnd.lb;
          tx1 = tx0 - 1;

          prv_sum = XMX(SP_B, q_0);
          prv_M = MY_Prod(MSMX(qx1, tx0),
                          MY_Prod(TSC(t_1, B2M), MSC(t_0, A)));
          XMX(SP_B, q_0) = MY_Sum(prv_sum, prv_M);
        }
      }
    }

    /* UPDATE SPECIAL STATES */
    XMX(SP_J, q_0) = MY_Zero();
    XMX(SP_C, q_0) = MY_Zero();
    XMX(SP_E, q_0) = MY_Zero();

    prv_N = MY_Prod(XMX(SP_N, q_1), XSC(SP_N, SP_LOOP));
    prv_B = MY_Prod(XMX(SP_B, q_0), XSC(SP_N, SP_MOVE));
    XMX(SP_N, q_0) = MY_Sum(prv_N, prv_B);
  }

  sc_best = XMX(SP_N, Q_range.beg);
  *sc_final = sc_best;

#if DEBUG
  {
  }
#endif

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
