/*******************************************************************************
 *    - FILE:       cloud_search_linear.c
 *    - DESC:     Cloud Search for Forward-Backward Pruning Algorithm
 *                (Linear Space Alg)
 *    NOTES:
 *       - Intermittant Bug
==1889360==ERROR: AddressSanitizer: heap-buffer-overflow on address 0x60c000003fe0 at pc 0x56500ca90ba9 bp 0x7fff9d366fd0 sp 0x7fff9d366fc0
WRITE of size 4 at 0x60c000003fe0 thread T0
    #0 0x56500ca90ba8 in run_Cloud_Forward_Linear src/algs_linear/cloud_search_linear.c:373
    #1 0x56500c9d13a6 in WORK_cloud_search_linear src/work/work_cloud_search.c:86
    #2 0x56500ca497ae in mmore_searchmmore_pipeline src/pipelines/pipeline_mmore_main.c:107
    #3 0x56500c9e43d4 in main src/application.c:74
    #4 0x7fa683d950b2 in __libc_start_main (/lib/x86_64-linux-gnu/libc.so.6+0x270b2)
    #5 0x56500c97a08d in _start (/home/devreckas/Google-Drive/Wheeler-Labs/personal-work/fbpruner-project/mmore/build/mmoreseqs-DEBUG+0x2208d)

Address 0x60c000003fe0 is a wild pointer.
SUMMARY: AddressSanitizer: heap-buffer-overflow src/algs_linear/cloud_search_linear.c:373 in run_Cloud_Forward_Linear

 *    TODO:
 *       - lb_vec and rb_vec should be created in main routine so we don't need to create/destroy every routine.
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

/* header */
#include "_algs_linear.h"
#include "cloud_search_linear.h"

/**
 *      NOTE: HOW TO CONVERT row-coords to diag-coords
 *       MMX3(q_1, t_1) => MMX3(d_2, k_1)
 *       MMX3(q_0, t_1) => MMX3(d_1, k_0)
 *       MMX3(q_1, t_0) => MMX3(d_1, k_1)
 */

/**
 *       NOTE: CONVERSION - row-coords => diag-coords
 *       MX_M(i-1, j-1) => MX3_M(d_2, k-1)
 *       MX_M(i  , j-1) => MX3_M(d_1, k  )
 *       MX_M(i-1, j  ) => MX3_M(d_1, k-1)
 *
 *       MX_M(i+1, j+1) => MX3_M(d_2, k+1)
 *       MX_M(i  , j+1) => MX3_M(d_1, k  )
 *       MX_M(i+1, j  ) => MX3_M(d_1, k+1)
 */

/* private functions */
/* Math functions */
static inline float
MY_Sum(float x, float y);

static inline float
MY_Prod(float x, float y);

static inline float
MY_Zero();

static inline float
MY_One();

STATUS_FLAG run_Cloud_Forward_Linear(
  const SEQUENCE* query,     /* query sequence */
  const HMM_PROFILE* target, /* target hmm model */
  const int Q,               /* query length */
  const int T,               /* target length */
  MATRIX_3D* st_MX3,         /* normal state matrix */
  MATRIX_2D* sp_MX,          /* special state matrix */
  const ALIGNMENT* tr,       /* viterbi traceback */
  EDGEBOUND_ROWS* rows,      /* temporary edgebounds by-row vector */
  EDGEBOUNDS* edg,           /* OUTPUT: edgebounds of cloud search space */
  CLOUD_PARAMS* params,      /* pruning parameters */
  float* inner_sc,           /* OUTPUT: maximum score inside viterbi bounds */
  float* max_sc              /* OUTPUT: highest score found during search */
) {

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
  int d_0, d_1, d_2;              /* real index of current and previous antidiagonals */
  int dx0, dx1, dx2;              /* mod mapping of antidiagonal index into data matrix */
  int k_0, k_1;                   /* offset into antidiagonal */
  int d_st, d_end, d_cnt, d_last; /* starting and ending diagonal indices */
  int dim_T, dim_Q, dim_TOT;      /* dimensions of submatrix being searched */
  int dim_min, dim_max;           /* diagonal index where num cells reaches highest point and diminishing point */
  int num_cells;                  /* number of cells in current diagonal */

  /* vars for indexing into edgebound lists */
  BOUND* bnd;      /* current bound */
  BOUND bnd_new;   /* for adding new bound to edgebound list */
  int id;          /* id in edgebound list (row/diag) */
  int r_0;         /* current index in edgebound list */
  int r_0b, r_0e;  /* begin and end indices for current row in edgebound list */
  int r_1;         /* current index for previous row */
  int r_1b, r_1e;  /* begin and end indices for current row in edgebound list */
  int le_0, re_0;  /* right/left matrix bounds of current diag */
  int lb_0, rb_0;  /* bounds of current search space on current diag */
  int lb_1, rb_1;  /* bounds of current search space on previous diag */
  int lb_2, rb_2;  /* bounds of current search space on 2-back diag */
  bool lb_T, rb_T; /* checks if edge touches bounds of matrix */

  /* vars for recurrance scores */
  float prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
  float prv_B, prv_E;           /* previous (B) begin and (E) end states */
  float prv_J, prv_N, prv_C;    /* previous (J) jump, (N) initial, and (C) terminal states */
  float prv_loop, prv_move;     /* previous loop and move for special states */
  float prv_sum, prv_best;      /* temp subtotaling vars */
  float sc_best;                /* final best scores */
  float sc_M, sc_I, sc_D, sc_E; /* match, insert, delete, end scores */

  /* vars for traceback */
  TRACE* beg; /* beginning of the alignment */
  TRACE* end; /* end of the alignment */

  /* vars for pruning */
  bool is_term_flag;             /* termination flag for end of search */
  float cell_max, diag_max;      /* temp maximum scores */
  float inner_max, total_max;    /* maximum score found in matrix */
  float total_limit, diag_limit; /* threshold determined by max_scores - alpha */
  VECTOR_INT* lb_vec[3];         /* left bound list for previous 3 antdiags */
  VECTOR_INT* rb_vec[3];         /* right bound list for previous 3 antidiags */
  VECTOR_INT* lb_vec_tmp;        /* left swap pointer */
  VECTOR_INT* rb_vec_tmp;        /* right swap pointer */

  /* pruning parameters */
  float alpha;
  float beta;
  int gamma;
  float hard_limit;
  /* antidiag range for the start/end points in the input viterbi alignment */
  RANGE vit_range;

  /* in order to approximate the score accurately, we need to know the high score position */
  COORDS coords_max;
  COORDS coords_innermax;

  /* debugger tools */
  FILE* dbfp;
  MATRIX_2D* cloud_MX;
  MATRIX_2D* cloud_MX3;
  MATRIX_3D* test_MX;
  MATRIX_3D* test_MX3;
  int num_writes;
  int num_clears;

/* initialize debugging matrix */
#if DEBUG
  {
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

  /* initialize logsum lookup table if it has not already been */
  MATH_Logsum_Init();
  int test_max_re = 0;
  int test_max_width = 0;

/* check if data is cleaned */
#if DEBUG
  {
    int cmp = MATRIX_3D_Check_Clean(st_MX3);
    printf("PRE-CHECK CLEAN -> CLOUD FWD?\t%d\n", cmp);
  }
#endif

  /* clear all old data from data matrix if necessary */
  if (st_MX3->clean = false) {
    MATRIX_3D_Clean(st_MX3);
  }

  /* query sequence */
  seq = query->seq;
  /* local or global alignments? */
  is_local = target->isLocal;
  sc_E = (is_local) ? MY_One() : MY_Zero();

  /* get start and end points of viterbi alignment */
  beg = &(tr->traces->data[tr->beg]);
  end = &(tr->traces->data[tr->end]);

  /* get pruning parameters */
  alpha = params->alpha;
  beta = params->beta;
  gamma = params->gamma;
  hard_limit = params->hard_limit;
  /* start and end points of input viterbi alignment */
  vit_range = (RANGE){beg->q_0 + beg->t_0, end->q_0 + end->t_0};
  /* position of maximum score in cloud */
  coords_innermax = (COORDS){-1, -1};
  coords_max = (COORDS){-1, -1};

  /* initialize edges and bounds */
  le_0 = 0;
  re_0 = 0;
  lb_0 = lb_1 = lb_2 = 0;
  rb_0 = rb_1 = rb_2 = 0;

  /* set edgebound dimensions and orientation */
  EDGEBOUNDS_Reuse(edg, Q, T);
#if (CLOUD_METHOD == CLOUD_DIAGS)
  {
    edg->edg_mode = EDG_DIAG;
  }
#elif (CLOUD_METHOD == CLOUD_ROWS)
  {
    EDGEBOUND_ROWS_Reuse(rows, Q, T, (RANGE){0, Q});
    edg->edg_mode = EDG_ROW;
  }
#endif

  /* malloc dynamic memory */
  for (i = 0; i < 3; i++) {
    lb_vec[i] = VECTOR_INT_Create();
    rb_vec[i] = VECTOR_INT_Create();
  }

  /* verify that starting points are valid */
  if (beg->q_0 < 0 || beg->q_0 > Q || beg->t_0 < 0 || beg->t_0 > T) {
    fprintf(stderr, "# ERROR: Invalid start points for Cloud Forward Search: BEG(%d,%d) -> END(%d,%d)\n", beg->q_0, beg->t_0, end->q_0, end->t_0);
    fprintf(stderr, "# Query Length: %d, Target Length: %d\n", Q, T);
#if DEBUG
    {
      /* temporary override: this should be handled in pipeline */
      if (end->q_0 > Q) {
        end->q_0 = Q;
      }
      if (end->t_0 > T) {
        end->t_0 = T;
      }
    }
#else
    {
      ERRORCHECK_exit(EXIT_FAILURE);
    }
#endif
  }

  /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
  if (beg->q_0 == 0 || beg->t_0 == 0) {
    beg->q_0 += 1;
    beg->t_0 += 1;
  }

  /* TODO: initialize previous start state */
  q_0 = beg->q_0;
  q_1 = q_0 - 1;
  t_0 = beg->t_0;
  t_1 = t_0 - 1;

  /* dimension of submatrix */
  dim_TOT = Q + T; /* antidiag dimensions */
  dim_Q = Q - beg->q_0;
  dim_T = T - beg->t_0;

  /* diag index at corners of dp matrix */
  d_st = 0;
  d_end = dim_TOT;

  /* diag index of different start points, creating submatrix */
  d_st = beg->q_0 + beg->t_0;
  // d_end = end->q_0 + end->t_0;

  /* diag index where num cells reaches highest point and begins diminishing */
  dim_min = MIN(d_st + dim_Q, d_st + dim_T);
  dim_max = MAX(d_st + dim_Q, d_st + dim_T);

  /* set bounds of starting cell */
  lb_0 = beg->q_0;
  rb_0 = beg->q_0;
  VECTOR_INT_Pushback(lb_vec[1], lb_0);
  VECTOR_INT_Pushback(rb_vec[1], rb_0);
  num_cells = 0;

  /* keeps largest number seen on current diagonal */
  is_term_flag = false;
  diag_max = MY_Zero();
  total_max = MY_Zero();
  /* number of passes through antidiags */
  d_cnt = 0;

  /* begin state probability begins at zero (free to start alignment) */
  prv_B = MY_One();
  // prv_E = MY_One();

  /* ITERATE THROUGH ANTI-DIAGONALS */
  /* Note: added num_cells condition */
  for (d_0 = d_st; d_0 <= d_end + 1; d_0++, d_cnt++) {
    d_1 = d_0 - 1; /* look back 1 antidiagonal */
    d_2 = d_0 - 2; /* look back 2 antidiagonal */
    /* mod-mapping of antidiagonals into linear space */
    dx0 = d_0 % 3;
    dx1 = d_1 % 3;
    dx2 = d_2 % 3;

    /* is dp matrix diagonal growing or shrinking? */
    if (d_0 <= dim_min) {
      num_cells++;
    }
    if (d_0 > dim_max) {
      num_cells--;
    }

    /* Edgecheck updates: determine antidiag indices within matrix bounds */
    le_0 = MAX(beg->q_0, d_0 - T);
    re_0 = le_0 + num_cells;

// test_max_re = MAX(test_max_re, re_0);
// test_max_width = MAX(test_max_width, re_0 - le_0 + 1);

/* Macro-controlled - Bounds pruning method */
#if (PRUNER == PRUNER_XDROP_EDGETRIM)
    {
      /* prune bounds using x-drop, no bifurcating */
      PRUNER_via_xdrop_edgetrim_Linear(
          st_MX3, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec);
    }
#elif (PRUNER == PRUNER_XDROP_BIFURCATE)
    {
      /* prune bounds using x-drop, bifurcating */
      PRUNER_via_xdrop_bifurcate_Linear(
          st_MX3, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec);
    }
#elif (PRUNER == PRUNER_DBL_XDROP_EDGETRIM_OR_DIE)
    {
      /* prune bounds using local and global x-drop, edgetrimming or terminating search */
      PRUNER_edgetrim_by_global_and_diag_Linear(
          st_MX3, sp_MX, alpha, beta, gamma, hard_limit,
          vit_range, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0,
          &total_max, &coords_max, &is_term_flag, lb_vec, rb_vec);
    }
#endif

    /* if currently inside viterbi range, update max inner_sc */
    if (d_0 < vit_range.end) {
      inner_max = total_max;
      coords_innermax = coords_max;
    }

    /* Add pruned bounds to edgebound list */
    for (i = 0; i < lb_vec[0]->N; i++) {
      lb_0 = VEC_X(lb_vec[0], i);
      rb_0 = VEC_X(rb_vec[0], i);

      /* Update bounds (spans all cells adjacent to previous antidiagonals cells that were not pruned) */
      lb_0 = lb_0;
      rb_0 = rb_0 + 1;

      /* Update bounds to account for dp matrix bounds */
      lb_0 = MAX(lb_0, le_0);
      rb_0 = MIN(rb_0, re_0);

      /* Update changes to list */
      VEC_X(lb_vec[0], i) = lb_0;
      VEC_X(rb_vec[0], i) = rb_0;

      /* Bound to be added */
      bnd_new = (BOUND){d_0, lb_0, rb_0};

/* macro-controlled: how to store cloud -> antidiag or row-wise */
#if (CLOUD_METHOD == CLOUD_DIAGS)
      {
        /* add new bounds to edgebounds as antidiag-wise */
        EDGEBOUNDS_Pushback(edg, bnd_new);
      }
/* antidiag is still WIP */
#elif (CLOUD_METHOD == CLOUD_ROWS)
      {
        /* reorient new bounds from antidiag-wise to row-wise and integrate it into row-wise edgebound list */
        EDGEBOUND_ROWS_IntegrateDiag_Fwd(rows, &bnd_new);

/* add new bounds to edgebounds as antidiag-wise (for comparative testing) */
#if DEBUG
        {
          EDGEBOUNDS_Pushback(test_edg, bnd_new);
        }
#endif
      }
#endif
    }

    /* If diagonal set is empty, then all branches have been pruned, so we're done */
    if (lb_vec[0]->N <= 0) {
      break;
    }

    /* MAIN RECURSION */
    /* Iterate through list of antidiagonal ranges */
    for (i = 0; i < lb_vec[0]->N; i++) {
      lb_0 = VEC_X(lb_vec[0], i);
      rb_0 = VEC_X(rb_vec[0], i);

      /* Iterate through cells in range */
      for (k_0 = lb_0; k_0 < rb_0; k_0++) {
        k_1 = k_0 - 1;

        /* row-col coords */
        q_0 = k_0;
        q_1 = q_0 - 1;
        t_0 = d_0 - k_0;
        t_1 = t_0 - 1;

        a = seq[q_0];
        A = AA_REV[a];

        /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
        /* best previous state transition (match takes the diag element of each prev state) */
        /* NOTE: Convert (q-1,t-1) <=> (d-2,k-1) */
        prv_M = MY_Prod(MMX3(dx2, k_1), TSC(t_1, M2M));
        prv_I = MY_Prod(IMX3(dx2, k_1), TSC(t_1, I2M));
        prv_D = MY_Prod(DMX3(dx2, k_1), TSC(t_1, D2M));
        /* Free to begin match state (new alignment) */
        /* NOTE: only allow begin transition at start of viterbi alignment */
        // prv_B only assigned once at start */
        /* best-to-match */
        prv_sum = MY_Sum(MY_Sum(prv_M, prv_I),
                         MY_Sum(prv_D, prv_B));
        MMX3(dx0, k_0) = prv_sum + MSC(t_0, A);

        /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
        /* previous states (match takes the left element of each state) */
        /* NOTE: Convert (q-1,t) <=> (d-1,k-1) */
        prv_M = MY_Prod(MMX3(dx1, k_1), TSC(t_0, M2I));
        prv_I = MY_Prod(IMX3(dx1, k_1), TSC(t_0, I2I));
        /* best-to-insert */
        prv_sum = MY_Sum(prv_M, prv_I);
        IMX3(dx0, k_0) = MY_Prod(prv_sum, ISC(t_0, A));

        /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
        /* previous states (match takes the left element of each state) */
        /* NOTE: Convert (q,t-1) <=> (d-1, k) */
        prv_M = MY_Prod(MMX3(dx1, k_0), TSC(t_1, M2D));
        prv_D = MY_Prod(DMX3(dx1, k_0), TSC(t_1, D2D));
        /* best-to-delete */
        prv_sum = MY_Sum(prv_M, prv_D);
        DMX3(dx0, k_0) = prv_sum;

/* embed cell data in quadratic matrix */
#if DEBUG
        {
          MX_2D(cloud_MX, q_0, t_0) += 1.0;
          MX_2D(cloud_MX3, dx0, k_0) += 1.0;

          MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(dx0, k_0);
          MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(dx0, k_0);
          MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(dx0, k_0);
        }
#endif
      }
    }

    /* Scrub values from 2-back bound data */
    for (i = 0; i < lb_vec[2]->N; i++) {
      lb_2 = lb_vec[2]->data[i];
      rb_2 = rb_vec[2]->data[i];

      for (k_0 = lb_2; k_0 < rb_2; k_0++) {
        q_0 = k_0;
        t_0 = d_2 - k_0;
        MMX3(dx2, k_0) = MY_Zero();
        IMX3(dx2, k_0) = MY_Zero();
        DMX3(dx2, k_0) = MY_Zero();
#if DEBUG
        {
          MX_2D(cloud_MX, q_0, t_0) += 2.0;
          MX_2D(cloud_MX3, dx2, k_0) -= 1.0;
        }
#endif
      }
    }

/* check that all necessary cells have been cleared */
#if MEMCHECK
    {
      bool is_clean = false;

      for (int k_0 = 0; k_0 < (Q + 1) + (T + 1); k_0++) {
        is_clean = false;
        is_clean += ((MMX3(dx2, k_0) == -INF) == false);
        is_clean += ((MMX3(dx2, k_0) == -INF) == false);
        is_clean += ((DMX3(dx2, k_0) == -INF) == false);
        if (is_clean != 0) {
          ERRORCHECK_memcheck(d_2, k_0, MMX3(dx2, k_0), IMX3(dx2, k_0), DMX3(dx2, k_0));
          MMX3(dx2, k_0) = IMX3(dx2, k_0) = DMX3(dx2, k_0) = -INF;
        }
      }
    }
#endif

    /* Shift bounds */
    lb_vec_tmp = lb_vec[2];
    rb_vec_tmp = rb_vec[2];
    lb_vec[2] = lb_vec[1];
    rb_vec[2] = rb_vec[1];
    lb_vec[1] = lb_vec[0];
    rb_vec[1] = rb_vec[0];
    lb_vec[0] = lb_vec_tmp;
    rb_vec[0] = rb_vec_tmp;

    /* scrub 2-back bound data (linear only) */
    VECTOR_INT_Reuse(lb_vec[0]);
    VECTOR_INT_Reuse(rb_vec[0]);

    /* disallow starting new alignments after first pass */
    prv_B = -INF;
    // prv_E = -INF;

    /* if termination condition has been triggered, then break out of loop */
    if (is_term_flag == true) {
      break;
    }
  }

  /* Scrub last two rows */
  d_last = d_0;
  for (d_0 = d_last; d_0 < d_last + 2; d_0++) {
    d_1 = d_0 - 1; /* look back 1 antidiagonal */
    d_2 = d_0 - 2; /* look back 2 antidiagonal */
    /* map antidiagonals to data matrix */
    dx0 = d_0 % 3;
    dx1 = d_1 % 3;
    dx2 = d_2 % 3;

    /* Scrub 2-back bound data */
    for (i = 0; i < lb_vec[2]->N; i++) {
      lb_2 = lb_vec[2]->data[i];
      rb_2 = rb_vec[2]->data[i];

      for (k_0 = lb_2; k_0 < rb_2; k_0++) {
        q_0 = k_0;
        t_0 = d_2 - k_0;
        MMX3(dx2, k_0) = MY_Zero();
        IMX3(dx2, k_0) = MY_Zero();
        DMX3(dx2, k_0) = MY_Zero();

#if DEBUG
        {
          MX_2D(cloud_MX, q_0, t_0) += 2.0;
          MX_2D(cloud_MX3, dx2, k_0) -= 1.0;
        }
#endif
      }
    }

    /* Shift bounds */
    lb_vec_tmp = lb_vec[2];
    rb_vec_tmp = rb_vec[2];
    lb_vec[2] = lb_vec[1];
    rb_vec[2] = rb_vec[1];
    lb_vec[1] = lb_vec[0];
    rb_vec[1] = rb_vec[0];
    lb_vec[0] = lb_vec_tmp;
    rb_vec[0] = rb_vec_tmp;
  }

/* check that all necessary cells have been cleared */
#if MEMCHECK
  {
    bool is_clean = false;

    for (dx0 = 0; dx0 < 3; dx0++) {
      for (int k_0 = 0; k_0 < (Q + 1) + (T + 1); k_0++) {
        is_clean = false;
        is_clean += ((MMX3(dx0, k_0) == -INF) == false);
        is_clean += ((MMX3(dx0, k_0) == -INF) == false);
        is_clean += ((DMX3(dx0, k_0) == -INF) == false);
        if (is_clean != 0) {
          ERRORCHECK_memcheck(dx0, k_0, MMX3(dx0, k_0), IMX3(dx0, k_0), DMX3(dx0, k_0));
          MMX3(dx0, k_0) = IMX3(dx0, k_0) = DMX3(dx0, k_0) = -INF;
        }
      }
    }
  }
#endif

/* if storing row-wise, compare to diag-wise */
#if (CLOUD_METHOD == CLOUD_ROWS)
  {
    /* output rows to edgebounds */
    EDGEBOUND_ROWS_Convert(rows, edg);

#if DEBUG
    {
      /* compare cloud rows method to antidiagonal method */
      int cmp = EDGEBOUNDS_Compare_by_Cloud_Single(cloud_MX, edg, test_edg);
      printf("COMPARE (rows vs antidiag):\t%s\n", (cmp == 0) ? "PASS" : "FAIL");
      if (cmp != 0) {
        EDGEBOUNDS_Dump(edg, stdout);
        EDGEBOUNDS_Dump(test_edg, stdout);

        printf("=== ROW-WISE ===\n");
        MATRIX_2D_Fill(cloud_MX, 0);
        MATRIX_2D_Cloud_Fill(cloud_MX, edg, 1);
        DP_MATRIX_VIZ_Dump(cloud_MX, stdout);

        printf("=== ROW-WISE ===\n");
        MATRIX_2D_Fill(cloud_MX, 0);
        MATRIX_2D_Cloud_Fill(cloud_MX, test_edg, 1);
        DP_MATRIX_VIZ_Dump(cloud_MX, stdout);

        printf("=== OVERLAY ===\n");
        DP_MATRIX_VIZ_Compare(cloud_MX, edg, test_edg);
        DP_MATRIX_VIZ_Dump(cloud_MX, stdout);
      }
    }
#endif
  }
#endif

/* check if data is cleaned */
#if DEBUG
  {
    int cmp = MATRIX_3D_Check_Clean(st_MX3);
    printf("POST-CHECK CLEAN -> CLOUD FWD?\t%d\n", cmp);
    printf("PRESCORES:: inner_max = %f, outer_max = %f\n", inner_max, total_max);
  }
#endif

  /* free dynamic memory ( TODO: create/destroy vector in main routine. ) */
  for (i = 0; i < 3; i++) {
    VECTOR_INT_Destroy(lb_vec[i]);
    VECTOR_INT_Destroy(rb_vec[i]);
  }

  /* after search, all cells are set to -INF */
  st_MX3->clean = true;

  bool is_score_correction = true;
  /* score correction: we need B to simulate proper model states */
  if (is_score_correction == true) {
    float presc, postsc;

    /* pre-core model: S->N->...->N->B->(M */
    presc = 0.0f;
    for (q_0 = 1; q_0 < beg->q_0; q_0++) {
      /* N loop */
      presc += XSC(SP_N, SP_LOOP);
    }
    t_1 = beg->t_0 - 1;
    presc = MY_Prod(presc, TSC(t_1, B2M));
    total_max = MY_Prod(total_max, presc);
    inner_max = MY_Prod(inner_max, presc);

    /* since total and inner ERRORCHECK_exit core model at different points, we compute their post-core corrections separately */
    /* post-core model: M)->E->C->...->C->T */
    postsc = 0.0f;
    postsc = MY_Prod(postsc, XSC(SP_E, SP_MOVE));
    for (q_0 = coords_max.q_0; q_0 <= Q; q_0++) {
      postsc = MY_Prod(postsc, XSC(SP_C, SP_LOOP));
    }
    postsc = MY_Prod(postsc, XSC(SP_C, SP_MOVE));
    total_max = MY_Prod(total_max, postsc);

    /* post-core model: M)->E->C->...->C->T */
    postsc = 0.0f;
    postsc = MY_Prod(postsc, XSC(SP_E, SP_MOVE));
    for (q_0 = coords_innermax.q_0; q_0 <= Q; q_0++) {
      postsc = MY_Prod(postsc, XSC(SP_C, SP_LOOP));
    }
    postsc = MY_Prod(postsc, XSC(SP_C, SP_MOVE));
    inner_max = MY_Prod(inner_max, postsc);
  }

#if DEBUG
  {
    printf("POSTSCORES:: inner_max = %f, outer_max = %f\n", inner_max, total_max);
  }
#endif

  /* highest score found in cloud search */
  *max_sc = total_max;
  *inner_sc = inner_max;

  // printf("## TEST RE => Q,T: %d,%d max_diag: %d...%d/%d, max_re: %d, max_width: %d\n", Q, T, d_st, d_last, dim_TOT, test_max_re, test_max_width);

  return STATUS_SUCCESS;
}

STATUS_FLAG
run_Cloud_Backward_Linear(const SEQUENCE* query,     /* query sequence */
                          const HMM_PROFILE* target, /* target hmm model */
                          const int Q,               /* query length */
                          const int T,               /* target length */
                          MATRIX_3D* st_MX3,         /* normal state matrix */
                          MATRIX_2D* sp_MX,          /* special state matrix */
                          const ALIGNMENT* tr,       /* viterbi traceback */
                          EDGEBOUND_ROWS* rows,      /* temporary edgebounds by-row */
                          EDGEBOUNDS* edg,           /* (OUTPUT) */
                          CLOUD_PARAMS* params,      /* pruning parameters */
                          float* inner_sc,           /* OUTPUT: maximum score inside viterbi bounds */
                          float* max_sc)             /* OUTPUT: highest score found during search */
{
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
  int d_0, d_1, d_2;              /* real index of current and previous antidiagonals */
  int dx0, dx1, dx2;              /* mod mapping of antidiagonal index into data matrix */
  int k_0, k_1;                   /* real index offset into diagonals */
  int d_st, d_end, d_cnt, d_last; /* starting and ending diagonal indices */
  int dim_T, dim_Q, dim_TOT;      /* dimensions of submatrix being searched */
  int dim_min, dim_max;           /* diagonal index where num cells reaches highest point and diminishing point */
  int num_cells;                  /* number of cells in current diagonal */

  /* vars for indexing into edgebound lists */
  BOUND* bnd;      /* current bound */
  BOUND bnd_new;   /* for adding new bound to edgebound list */
  int id;          /* id in edgebound list (row/diag) */
  int r_0;         /* current index in edgebound list */
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
  float prv_sum, prv_best;      /* temp subtotaling vars */
  float sc_best;                /* final best scores */
  float sc_M, sc_I, sc_D, sc_E; /* match, insert, delete, end scores */

  /* vars for traceback */
  TRACE* beg; /* beginning of the alignment */
  TRACE* end; /* end of the alignment */

  /* vars for pruning */
  bool is_term_flag;             /* termination flag set by pruner */
  float cell_max, diag_max;      /* temp maximum scores */
  float inner_max, total_max;    /* maximum score found in matrix */
  float total_limit, diag_limit; /* threshold determined by max_scores - alpha */
  VECTOR_INT* lb_vec[3];         /* left bound list for previous 3 antdiags */
  VECTOR_INT* rb_vec[3];         /* right bound list for previous 3 antidiags */
  VECTOR_INT* lb_vec_tmp;        /* left swap pointer */
  VECTOR_INT* rb_vec_tmp;        /* right swap pointer */

  /* pruning parameters */
  float alpha;
  float beta;
  int gamma;
  float hard_limit;
  RANGE vit_range;

  /* in order to approximate the score accurately, we need to know the high score position */
  COORDS coords_max;
  COORDS coords_innermax;

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

/* check if data is cleaned */
#if DEBUG
  {
    int cmp = MATRIX_3D_Check_Clean(st_MX3);
    printf("PRE-CHECK CLEAN -> CLOUD BCK?\t%d\n", cmp);
  }
#endif

  /* clear all old data from data matrix if necessary */
  if (st_MX3->clean = false) {
    MATRIX_3D_Clean(st_MX3);
  }

  /* query sequence */
  seq = query->seq;
  /* local or global alignments? */
  is_local = target->isLocal;
  sc_E = (is_local) ? 0 : -INF;

  /* get start and end points of viterbi alignment */
  beg = &(tr->traces->data[tr->beg]);
  end = &(tr->traces->data[tr->end]);

  /* get pruning parameters */
  alpha = params->alpha;
  beta = params->beta;
  gamma = params->gamma;
  hard_limit = params->hard_limit;
  /* antidiag range for the start/end points in the input viterbi alignment */
  vit_range = (RANGE){beg->q_0 + beg->t_0, end->q_0 + end->t_0};

  /* initialize edges and bounds */
  le_0 = 0;
  re_0 = 0;
  lb_0 = lb_1 = lb_2 = 0;
  rb_0 = rb_1 = rb_2 = 0;

  /* set edgebound dimensions and orientation */
  EDGEBOUNDS_Reuse(edg, Q, T);
#if (CLOUD_METHOD == CLOUD_DIAGS)
  {
    edg->edg_mode = EDG_DIAG;
  }
#elif (CLOUD_METHOD == CLOUD_ROWS)
  {
    EDGEBOUND_ROWS_Reuse(rows, Q, T, (RANGE){0, Q});
    edg->edg_mode = EDG_ROW;
  }
#endif

  /* malloc dynamic memory */
  for (i = 0; i < 3; i++) {
    lb_vec[i] = VECTOR_INT_Create();
    rb_vec[i] = VECTOR_INT_Create();
  }

  /* verify that starting points are valid */
  if (end->q_0 < 0 || end->q_0 > Q || end->t_0 < 0 || end->t_0 > T) {
    fprintf(stderr, "# ERROR: Invalid start points for Cloud Backward Search: BEG(%d,%d) -> END(%d,%d)\n", beg->q_0, beg->t_0, end->q_0, end->t_0);
    fprintf(stderr, "# Query Length: %d, Target Length: %d\n", Q, T);
#if DEBUG
    {
      /* temporary override: this should be handled in pipeline */
      if (end->q_0 > Q) {
        end->q_0 = Q;
      }
      if (end->t_0 > T) {
        end->t_0 = T;
      }
    }
#else
    {
      ERRORCHECK_exit(EXIT_FAILURE);
    }
#endif
  }

  /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
  if (end->q_0 == Q || end->t_0 == T) {
    end->q_0 -= 1;
    end->t_0 -= 1;
  }

  /* TODO: initialize previous start state */
  q_0 = end->q_0;
  q_1 = q_0 + 1;
  t_0 = end->t_0;
  t_1 = t_0 + 1;

  /* dimension of submatrix */
  dim_TOT = Q + T;
  dim_Q = end->q_0;
  dim_T = end->t_0;

  /* diag index at corners of dp matrix */
  d_st = 0;
  d_end = dim_TOT;

  /* diag index of different start points, creating submatrix */
  // d_st = beg->q_0 + beg->t_0;
  d_end = end->q_0 + end->t_0;

  /* diag index where num cells reaches highest point and begins diminishing */
  dim_min = MIN(dim_T, dim_Q);
  dim_max = MAX(dim_T, dim_Q);

  /* set bounds of starting cell */
  lb_0 = end->q_0;
  rb_0 = end->q_0 + 1;
  VECTOR_INT_Pushback(lb_vec[1], lb_0);
  VECTOR_INT_Pushback(rb_vec[1], rb_0);
  num_cells = 0;

  /* keeps largest number seen on current diagonal */
  is_term_flag = false;
  diag_max = -INF;
  total_max = -INF;
  /* number of antidiags passed through */
  d_cnt = 0;

  /* begin state probability begins at zero (free to start alignment) */
  // prv_B = MY_One();
  prv_E = MY_One();

  /* ITERATE THROUGHT ANTI-DIAGONALS */
  for (d_0 = d_end; d_0 >= d_st; d_0--, d_cnt++) {
    d_1 = d_0 + 1; /* look back 1 diagonal */
    d_2 = d_0 + 2; /* look back 2 diagonals */
    /* mod-mapping of antidiagonals into linear space */
    dx0 = d_0 % 3;
    dx1 = d_1 % 3;
    dx2 = d_2 % 3;

    /* Is dp matrix diagonal growing or shrinking? */
    if (d_0 >= dim_max) {
      num_cells++;
    }
    if (d_0 < dim_min) {
      num_cells--;
    }

    /* Edgecheck updates: determine antidiag indices within matrix bounds */
    le_0 = MAX(end->q_0 - (d_end - d_0), 0);
    re_0 = le_0 + num_cells;

/* Prune bounds */
#if (PRUNER == PRUNER_XDROP_EDGETRIM)
    {
      /* prune bounds using x-drop, no bifurcating */
      PRUNER_via_xdrop_edgetrim_Linear(
          st_MX3, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec);
    }
#elif (PRUNER == PRUNER_XDROP_BIFURCATE)
    {
      /* prune bounds using x-drop, bifurcating */
      PRUNER_via_xdrop_bifurcate_Linear(
          st_MX3, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec);
    }
#elif (PRUNER == PRUNER_DBL_XDROP_EDGETRIM_OR_DIE)
    {
      /* prune bounds using both local and global x-drop, edgetrimming or terminating search */
      PRUNER_edgetrim_by_global_and_diag_Linear(
          st_MX3, sp_MX, alpha, beta, gamma, hard_limit,
          vit_range, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0,
          &total_max, &coords_max, &is_term_flag, lb_vec, rb_vec);
    }
#endif

    /* if currently inside viterbi range, update max inner_sc */
    if (d_0 >= vit_range.beg) {
      inner_max = total_max;
      coords_innermax = coords_max;
    }

    /* Add pruned bounds to edgebound list */
    for (i = 0; i < lb_vec[0]->N; i++) {
      /* pull bounds from list */
      lb_0 = VEC_X(lb_vec[0], i);
      rb_0 = VEC_X(rb_vec[0], i);

      /* Update bounds (spans all cells adjacent to previous antidiagonals cells that were not pruned) */
      lb_0 = lb_0 - 1;
      rb_0 = rb_0;

      /* Update bounds to account for dp matrix bounds */
      lb_0 = MAX(lb_0, le_0);
      rb_0 = MIN(rb_0, re_0);
      /* zeroth row (when d_0 == k_0) is not a valid state in backward. */
      rb_0 = MIN(rb_0, d_0);

      /* Update changes to list */
      VEC_X(lb_vec[0], i) = lb_0;
      VEC_X(rb_vec[0], i) = rb_0;

      bnd_new = (BOUND){d_0, lb_0, rb_0};

#if (CLOUD_METHOD == CLOUD_DIAGS)
      {
        /* add new bounds to edgebounds as antidiag-wise */
        EDGEBOUNDS_Pushback(edg, bnd_new);
      }
#endif

#if (CLOUD_METHOD == CLOUD_ROWS)
      {
        /* reorient new bounds from antidiag-wise to row-wise and integrate it into row-wise edgebound list */
        EDGEBOUND_ROWS_IntegrateDiag_Bck(rows, &bnd_new);

#if DEBUG
        {
          /* add new bounds to edgebounds as antidiag-wise (for comparative testing) */
          EDGEBOUNDS_Pushback(test_edg, bnd_new);
        }
#endif
      }
#endif
    }

    /* If diagonal set is empty, then all branches have been pruned, so we're done */
    if (lb_vec[0]->N <= 0) {
      break;
    }

    /* MAIN RECURSION */
    for (i = 0; i < lb_vec[0]->N; i++) {
      lb_0 = VEC_X(lb_vec[0], i);
      rb_0 = VEC_X(rb_vec[0], i);
      /* if we are in an edgecase */
      // lb_T = ( lb_0 <= 0 );
      // rb_T = ( rb_0 >= d_0 );
      /* zeroth row (when d_0 == k_0) is not a valid state in backward. */
      // rb_0 = MIN( rb_0, d_0 );

      /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
      for (k_0 = lb_0; k_0 < rb_0; k_0++) {
        k_1 = k_0 + 1;

        /* get x-y coords */
        q_0 = k_0;
        q_1 = k_0 + 1;
        t_0 = d_0 - k_0;
        t_1 = t_0 + 1;

        /*    === ROW-WISE to DIAG_WISE ===
         *    MX_M(i+1, j+1) => MX3_M(d_2, k+1)
         *    MX_M(i  , j+1) => MX3_M(d_1, k  )
         *    MX_M(i+1, j  ) => MX3_M(d_1, k+1)
         */

        /* next sequence character */
        a = seq[q_0];
        A = AA_REV[a];

        /* match and insertion scores */
        sc_M = MSC(t_1, A);
        sc_I = ISC(t_1, A);

        /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
        prv_M = MY_Prod(MMX3(dx2, k_1),
                        MY_Prod(TSC(t_0, M2M), sc_M));
        prv_I = MY_Prod(IMX3(dx1, k_1),
                        MY_Prod(TSC(t_0, M2I), sc_I));
        prv_D = MY_Prod(DMX3(dx1, k_0), TSC(t_0, M2D));
        // prv_E = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
        // prv_E = sc_E;
        /* best-to-match */
        prv_sum = MY_Sum(MY_Sum(prv_M, prv_I), MY_Sum(prv_D, prv_E));
        MMX3(dx0, k_0) = prv_sum;

        // printf("(q_0,t_0)=%d,%d => a=%c, A=%d, tsc=%f, msc=%f mmx=%f\n",
        //    q_0, t_0, a, A, TSC(t_1, M2M), MSC(t_0, A), MMX3(dx0, k_0) );

        /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
        prv_M = MY_Prod(MMX3(dx2, k_1), MY_Prod(TSC(t_0, I2M), sc_M));
        prv_I = MY_Prod(IMX3(dx1, k_1),
                        MY_Prod(TSC(t_0, I2I), sc_I));
        /* best-to-insert */
        prv_sum = MY_Sum(prv_M, prv_I);
        IMX3(dx0, k_0) = prv_sum;

        /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
        prv_M = MY_Prod(MMX3(dx2, k_1),
                        MY_Prod(TSC(t_0, D2M), sc_M));
        prv_D = MY_Prod(DMX3(dx1, k_0), TSC(t_0, D2D));
        /* best-to-delete */
        prv_sum = MY_Sum(prv_M, prv_D);
        prv_sum = MY_Sum(prv_sum, prv_E);
        DMX3(dx0, k_0) = prv_sum;

/* embed cell data in quadratic matrix */
#if DEBUG
        {
          MX_2D(cloud_MX, q_0, t_0) += 1.0;
          MX_2D(cloud_MX3, dx0, k_0) += 1.0;

          MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(dx0, k_0);
          MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(dx0, k_0);
          MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(dx0, k_0);
        }
#endif
      }
    }

    /* Scrub 2-back bound data */
    for (i = 0; i < lb_vec[2]->N; i++) {
      lb_2 = lb_vec[2]->data[i];
      rb_2 = rb_vec[2]->data[i];

      for (k_0 = lb_2; k_0 < rb_2; k_0++) {
        q_0 = k_0;
        t_0 = d_2 - q_0;

        MMX3(dx2, k_0) = IMX3(dx2, k_0) = DMX3(dx2, k_0) = -INF;

#if DEBUG
        {
          MX_2D(cloud_MX, q_0, t_0) += 2.0;
          MX_2D(cloud_MX3, dx2, k_0) -= 1.0;
        }
#endif
      }
    }

/* check that all necessary cells have been cleared */
#if MEMCHECK
    {
      bool is_clean = false;

      for (int k_0 = 0; k_0 < (Q + 1) + (T + 1); k_0++) {
        is_clean = false;
        is_clean += ((MMX3(dx2, k_0) == -INF) == false);
        is_clean += ((MMX3(dx2, k_0) == -INF) == false);
        is_clean += ((DMX3(dx2, k_0) == -INF) == false);

        if (is_clean != 0) {
          ERRORCHECK_memcheck(d_2, k_0, MMX3(dx2, k_0), IMX3(dx2, k_0), DMX3(dx2, k_0));
          MMX3(dx2, k_0) = IMX3(dx2, k_0) = DMX3(dx2, k_0) = -INF;
          is_clean = 0;
        }
      }
    }
#endif

    /* Shift bounds */
    lb_vec_tmp = lb_vec[2];
    rb_vec_tmp = rb_vec[2];
    lb_vec[2] = lb_vec[1];
    rb_vec[2] = rb_vec[1];
    lb_vec[1] = lb_vec[0];
    rb_vec[1] = rb_vec[0];
    lb_vec[0] = lb_vec_tmp;
    rb_vec[0] = rb_vec_tmp;

    /* scrub 2-back bound data (linear only) */
    VECTOR_INT_Reuse(lb_vec[0]);
    VECTOR_INT_Reuse(rb_vec[0]);

    /* disallow starting new alignments after first pass */
    // prv_B = MY_Zero();
    prv_E = MY_Zero();

    /* if termination flag is set, break out of loop */
    if (is_term_flag == true) {
      break;
    }
  }

  /* scrub last two rows */
  d_last = d_0;
  for (d_0 = d_last; d_0 > d_last - 2; d_0--) {
    d_1 = d_0 + 1; /* look back 1 diagonal */
    d_2 = d_0 + 2; /* look back 2 diagonals */
    /* mod-mapping of antidiagonals into linear space */
    dx0 = d_0 % 3;
    dx1 = d_1 % 3;
    dx2 = d_2 % 3;

    /* Scrub 2-back bound data */
    for (i = 0; i < lb_vec[2]->N; i++) {
      lb_2 = lb_vec[2]->data[i];
      rb_2 = rb_vec[2]->data[i];

      for (k_0 = lb_2; k_0 < rb_2; k_0++) {
        q_0 = k_0;
        t_0 = d_2 - q_0;

        MMX3(dx2, k_0) = MY_Zero();
        IMX3(dx2, k_0) = MY_Zero();
        DMX3(dx2, k_0) = MY_Zero();

#if DEBUG
        {
          MX_2D(cloud_MX, q_0, t_0) += 2.0;
          MX_2D(cloud_MX3, dx2, k_0) -= 1.0;
        }
#endif
      }
    }

    /* Shift bounds */
    lb_vec_tmp = lb_vec[2];
    rb_vec_tmp = rb_vec[2];
    lb_vec[2] = lb_vec[1];
    rb_vec[2] = rb_vec[1];
    lb_vec[1] = lb_vec[0];
    rb_vec[1] = rb_vec[0];
    lb_vec[0] = lb_vec_tmp;
    rb_vec[0] = rb_vec_tmp;
  }

/* check that all necessary cells have been cleared */
#if MEMCHECK
  {
    bool is_clean = false;

    for (dx0 = 0; dx0 < 3; dx0++) {
      for (int k_0 = 0; k_0 < (Q + 1) + (T + 1); k_0++) {
        is_clean = false;
        is_clean += ((MMX3(dx0, k_0) == -INF) == false);
        is_clean += ((MMX3(dx0, k_0) == -INF) == false);
        is_clean += ((DMX3(dx0, k_0) == -INF) == false);
        if (is_clean != 0) {
          ERRORCHECK_memcheck(dx0, k_0, MMX3(dx0, k_0), IMX3(dx0, k_0), DMX3(dx0, k_0));
          MMX3(dx0, k_0) = IMX3(dx0, k_0) = DMX3(dx0, k_0) = -INF;
        }
      }
    }
  }
#endif

  /* reverse order of diagonals */
  EDGEBOUNDS_Reverse(edg);

  /* free dynamic memory */
  for (i = 0; i < 3; i++) {
    VECTOR_INT_Destroy(lb_vec[i]);
    VECTOR_INT_Destroy(rb_vec[i]);
  }

#if (CLOUD_METHOD == CLOUD_ROWS)
  {
    /* output rows to edgebounds */
    EDGEBOUND_ROWS_Convert(rows, edg);

#if DEBUG
    {
      /* compare cloud rows method to antidiagonal method */
      int cmp = EDGEBOUNDS_Compare_by_Cloud_Single(cloud_MX, edg, test_edg);
      printf("COMPARE (rows vs antidiag):\t%s\n", (cmp == 0) ? "PASS" : "FAIL");
      if (cmp != 0) {
        EDGEBOUNDS_Dump(edg, stdout);
        EDGEBOUNDS_Dump(test_edg, stdout);

        printf("=== DIAG-WISE ===\n");
        MATRIX_2D_Fill(cloud_MX, 0);
        MATRIX_2D_Cloud_Fill(cloud_MX, edg, 1);
        DP_MATRIX_VIZ_Dump(cloud_MX, stdout);

        printf("=== ROW-WISE ===\n");
        MATRIX_2D_Fill(cloud_MX, 0);
        MATRIX_2D_Cloud_Fill(cloud_MX, test_edg, 1);
        DP_MATRIX_VIZ_Dump(cloud_MX, stdout);

        printf("=== OVERLAY ===\n");
        DP_MATRIX_VIZ_Compare(cloud_MX, edg, test_edg);
        DP_MATRIX_VIZ_Dump(cloud_MX, stdout);
      }
    }
#endif
  }
#endif

/* check if data is cleaned */
#if DEBUG
  {
    int cmp = MATRIX_3D_Check_Clean(st_MX3);
    printf("POST-CHECK CLEAN -> CLOUD BCK?\t%d\n", cmp);
    printf("PRESCORES:: inner_max = %f, outer_max = %f\n", inner_max, total_max);
  }
#endif

  st_MX3->clean = true;

  bool is_score_correction = true;
  /* score correction: we need B to simulate proper model states */
  if (is_score_correction == true) {
    float presc, postsc;

    /* since total and inner enter the core model at different points, compute them separately */
    /* pre-core model: S->N->...->N->B->(M */
    presc = 0.0f;
    for (q_0 = 1; q_0 < coords_max.q_0; q_0++) {
      /* N loop */
      presc = MY_Prod(presc, XSC(SP_N, SP_LOOP));
    }
    t_1 = beg->t_0 - 1;
    /* N->B->M */
    presc = MY_Prod(presc, TSC(t_1, B2M));
    total_max = MY_Prod(total_max, presc);

    /* pre-core model: S->N->...->N->B->(M */
    presc = 0.0f;
    for (q_0 = 1; q_0 < coords_innermax.q_0; q_0++) {
      /* N loop */
      presc = MY_Prod(presc, XSC(SP_N, SP_LOOP));
    }
    t_1 = beg->t_0 - 1;
    /* N->B->M */
    presc = MY_Prod(presc, TSC(t_1, B2M));
    inner_max = MY_Prod(inner_max, presc);

    /* since total and inner exit at the same point, compute them together */
    /* post-core model: M)->E->C->...->C->T */
    postsc = 0.0f;
    /* M->E->C */
    postsc = MY_Prod(postsc, XSC(SP_E, SP_MOVE));
    for (q_0 = end->q_0; q_0 <= Q; q_0++) {
      /* C loop */
      postsc = MY_Prod(postsc, XSC(SP_C, SP_LOOP));
    }
    /* C->T */
    postsc = MY_Prod(postsc, XSC(SP_C, SP_MOVE));
    total_max = MY_Prod(total_max, postsc);
    inner_max = MY_Prod(inner_max, postsc);
  }

#if DEBUG
  {
    printf("POSTSCORES:: inner_max = %f, outer_max = %f\n", inner_max, total_max);
  }
#endif

  /* highest score found in cloud search */
  *max_sc = total_max;
  *inner_sc = inner_max;

  return STATUS_SUCCESS;
}

/* MATH RULES: These determine how probilities are summed and certain identities */

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
MY_Zero() {
  return MATH_LogZero();
}

static inline float
MY_One() {
  return MATH_LogOne();
}
