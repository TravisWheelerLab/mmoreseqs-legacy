/*******************************************************************************
 *  - FILE:  pruning_linear.c
 *  - DESC:  Pruning methods for Cloud Search.
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
#include "pruning_linear.h"

inline STATUS_FLAG
PRUNER_via_xdrop_edgetrim_Linear(MATRIX_3D* st_MX3,     /* normal state matrix */
                                 MATRIX_2D* sp_MX,      /* special state matrix */
                                 const float alpha,     /* x-drop value */
                                 const int gamma,       /* number of antidiagonals before pruning */
                                 const int d_1,         /* previous antidiagonal */
                                 const int d_0,         /* current antidiagonal */
                                 const int dx1,         /* previous antidiag (mod-mapped) */
                                 const int dx0,         /* current antidiag (mod-mapped) */
                                 const int d_cnt,       /* number of antidiags traversed */
                                 const int le,          /* right edge of dp matrix on current antidiag */
                                 const int re,          /* left edge of dp matrix on current antidiag */
                                 float* total_max,      /* (UPDATED) current maximum score */
                                 VECTOR_INT* lb_vec[3], /* OUTPUT: current list of left-bounds */
                                 VECTOR_INT* rb_vec[3]) /* OUTPUT: current list of right-bounds */
{
  int i, j, k;                   /* indexes */
  int q_0;                       /* row index (query position) */
  int t_0;                       /* column index (target position) */
  int k_0;                       /* offset into antidiagonal */
  int lb_0, rb_0;                /* left/right bounds of current antidiagonal */
  int lb_1, rb_1;                /* left/right bounds of previous antidiagonal */
  float diag_max, cell_max;      /* max score for all normal states in a given cell/antidiagonal */
  float diag_limit, total_limit; /* pruning thresholds for current limit */

  /* reset int vectors */
  VECTOR_INT_Reuse(lb_vec[0]);
  VECTOR_INT_Reuse(rb_vec[0]);

  /* Update maximum score using antidiagonal */
  for (i = 0; i < lb_vec[1]->N; i++) {
    lb_1 = VEC_X(lb_vec[1], i);
    rb_1 = VEC_X(rb_vec[1], i);

    diag_max = -INF;
    for (k_0 = lb_1; k_0 < rb_1; k_0++) {
      q_0 = k_0;
      t_0 = d_1 - k_0; /* looking back one diag */

      diag_max = MATH_Max(
          MATH_Max(diag_max, MMX3(dx1, k_0)),
          MATH_Max(IMX3(dx1, k_0), DMX3(dx1, k_0)));
    }

    /* Total max records largest cell score seen so far */
    *total_max = MAX(*total_max, diag_max);
  }

  /* Set score threshold for pruning */
  total_limit = *total_max - alpha;

  /* All edgebounds in previous antidiagonal */
  for (i = 0; i < lb_vec[1]->N; i++) {
    lb_1 = VEC_X(lb_vec[1], i);
    rb_1 = VEC_X(rb_vec[1], i);

    /* If free passes are not complete, no pruning */
    if (gamma >= d_cnt) {
      VECTOR_INT_Pushback(lb_vec[0], lb_1);
      VECTOR_INT_Pushback(rb_vec[0], rb_1);
    } else /* If free passes are complete (gamma < d), prune and set new edgebounds */
    {
      /* impossible state */
      lb_0 = INT_MIN;
      rb_0 = INT_MIN;

      /* Find the first cell from the left which passes above threshold */
      for (q_0 = lb_1; q_0 < rb_1; q_0++) {
        q_0 = k_0;
        t_0 = d_1 - k_0; /* looking back one diag */

        cell_max = MATH_Max(MMX3(dx1, k_0),
                            MATH_Max(IMX3(dx1, k_0), DMX3(dx1, k_0)));

        /* prune in left edgebound */
        if (cell_max >= total_limit) {
          lb_0 = k_0;
          VECTOR_INT_Pushback(lb_vec[0], lb_0);
          break;
        }
      }

      /* If no boundary edges are found on diag, then branch is pruned entirely */
      if (lb_0 == INT_MIN)
        continue;

      /* Find the first cell from the right which passes above threshold */
      for (k = rb_1 - 1; k >= lb_1; k--) {
        q_0 = k_0;
        t_0 = d_1 - k_0; /* looking back one diag */

        cell_max = MATH_Max(MMX3(dx1, k_0),
                            MATH_Max(IMX3(dx1, k_0), DMX3(dx1, k_0)));

        /* prune in right edgebound */
        if (cell_max >= total_limit) {
          rb_0 = k_0 + 1;
          VECTOR_INT_Pushback(rb_vec[0], rb_0);
          break;
        }
      }
    }
  }

  return STATUS_SUCCESS;
}

STATUS_FLAG
PRUNER_via_xdrop_bifurcate_Linear(MATRIX_3D* st_MX3,     /* normal state matrix */
                                  MATRIX_2D* sp_MX,      /* special state matrix */
                                  const float alpha,     /* x-drop value */
                                  const int gamma,       /* number of antidiagonals before pruning */
                                  const int d_1,         /* previous antidiagonal */
                                  const int d_0,         /* current antidiagonal */
                                  const int dx1,         /* previous antidiag (mod-mapped) */
                                  const int dx0,         /* current antidiag (mod-mapped) */
                                  const int d_cnt,       /* number of antidiags traversed */
                                  const int le,          /* right edge of dp matrix on current antidiag */
                                  const int re,          /* left edge of dp matrix on current antidiag */
                                  float* total_max,      /* (UPDATED) current maximum score */
                                  VECTOR_INT* lb_vec[3], /* OUTPUT: current list of left-bounds */
                                  VECTOR_INT* rb_vec[3]) /* OUTPUT: current list of right-bounds */
{
  int i, j, k;                   /* indexes */
  int q_0;                       /* row index (query position) */
  int t_0;                       /* column index (target position) */
  int k_0;                       /* offset into antidiagonal */
  int lb_0, rb_0;                /* left/right bounds of current antidiagonal */
  int lb_1, rb_1;                /* left/right bounds of previous antidiagonal */
  float diag_max, cell_max;      /* max score for all normal states in a given cell/antidiagonal */
  float diag_limit, total_limit; /* pruning thresholds for current limit */

  /* reset int vectors */
  VECTOR_INT_Reuse(lb_vec[0]);
  VECTOR_INT_Reuse(rb_vec[0]);

  /* Update maximum score using antidiagonal */
  for (i = 0; i < lb_vec[1]->N; i++) {
    lb_1 = lb_vec[1]->data[i];
    rb_1 = rb_vec[1]->data[i];

    diag_max = -INF;
    for (k_0 = lb_1; k_0 < rb_1; k_0++) {
      q_0 = k_0;
      t_0 = d_1 - k_0; /* looking back one diag */

      diag_max = MATH_Max(MATH_Max(diag_max, MMX3(dx1, k)),
                          MATH_Max(IMX3(dx1, k), DMX3(dx1, k)));
    }

    /* Total max records largest cell score seen so far */
    *total_max = MAX(*total_max, diag_max);
  }

  /* Set score threshold for pruning */
  total_limit = *total_max - alpha;

  /* All edgebounds in previous antidiagonal */
  for (i = 0; i < lb_vec[1]->N; i++) {
    lb_1 = lb_vec[1]->data[i];
    rb_1 = rb_vec[1]->data[i];

    /* If free passes are not complete, do no pruning */
    if (gamma >= d_cnt) {
      VECTOR_INT_Pushback(lb_vec[0], lb_1);
      VECTOR_INT_Pushback(rb_vec[0], rb_1);
    } else /* If free passes are complete (gamma < d), prune and set new edgebounds */
    {
      /* impossible state */
      lb_0 = INT_MIN;
      rb_0 = INT_MIN;

      /* Find the first cell from the left which passes above threshold */
      for (k = lb_1; k < rb_1; k++) {
        q_0 = k_0;
        t_0 = d_1 - k_0; /* looking back one diag */

        cell_max = MATH_Max(MMX3(dx1, k_0),
                            MATH_Max(IMX3(dx1, k_0), DMX3(dx1, k_0)));

        /* prune in left edgebound */
        if (cell_max >= total_limit) {
          lb_0 = k_0;
          VECTOR_INT_Pushback(lb_vec[0], lb_0);
          break;
        }
      }

      /* If no boundary edges are found on diag, then branch is pruned entirely */
      if (lb_0 == INT_MIN)
        continue;

      /* Find the first cell from the right which passes above threshold */
      for (k_0 = rb_1 - 1; k_0 >= lb_1; k_0--) {
        q_0 = k_0;
        t_0 = d_1 - k_0; /* looking back one diag */

        cell_max = MATH_Max(MMX3(dx1, k_0),
                            MATH_Max(IMX3(dx1, k_0), DMX3(dx1, k_0)));

        /* prune in right edgebound */
        if (cell_max >= total_limit) {
          rb_0 = k_0 + 1;
          VECTOR_INT_Pushback(rb_vec[0], rb_0);
          break;
        }
      }
    }
  }
  return STATUS_SUCCESS;
}

/*! FUNCTION: 	PRUNER_via_dbl_xdrop_edgetrim_by_global_and_diag_Linear()
 *  SYNOPSIS: 	Prunes antidiagonal of Cloud Search.
 * 				Uses x-drop and only trims in from left and right ends of search space. No bifurcation.
 *					Alpha: 		x-drop value determining whether cells are pruned in antidiagonal
 * 				Alpha-Max: 	x-drop value determining whether search is terminated
 *					Beta:  		number of free passes before pruning
 * 				(1) Updates the total_max, which stores the highest scoring cell in the matrix thus far.
 * 				(1b) If diag_max falls below score global threshold, terminate entire search.
 * 		      (2) Performs pruning from left-edge, moving right until a cell is found which all states fall below limit = (total_max - alpha).
 * 				(3) Performs pruning from right-edge, moving left.
 * 				(4) Stores edgebounds in left and right bound list.
 */
inline STATUS_FLAG
PRUNER_edgetrim_by_global_and_diag_Linear(MATRIX_3D* st_MX3,      /* normal state matrix */
                                          MATRIX_2D* sp_MX,       /* special state matrix */
                                          const float alpha,      /* x-drop value for by-diag prune */
                                          const float beta,       /* x-drop value for global prune */
                                          const int gamma,        /* number of antidiagonals before pruning */
                                          const float hard_limit, /* hard floor value for global prune */
                                          const RANGE vit_range,  /* antidiagonal locations for the start-end of the input viterbi alignment */
                                          const int d_1,          /* previous antidiagonal */
                                          const int d_0,          /* current antidiagonal */
                                          const int dx1,          /* previous antidiag (mod-mapped) */
                                          const int dx0,          /* current antidiag (mod-mapped) */
                                          const int d_cnt,        /* number of antidiags traversed */
                                          const int le,           /* right edge of dp matrix on current antidiag */
                                          const int re,           /* left edge of dp matrix on current antidiag */
                                          float* total_max,       /* UPDATED: current maximum score */
                                          COORDS* coords_max,     /* UPDATED: location of maximum score */
                                          bool* is_term_flag,     /* UPDATED: if termination trigger has been reached */
                                          VECTOR_INT* lb_vec[3],  /* OUTPUT: current list of left-bounds */
                                          VECTOR_INT* rb_vec[3])  /* OUTPUT: current list of right-bounds */
{
  int i, j, k;              /* indexes */
  int q_0;                  /* row index (query position) */
  int t_0;                  /* column index (target position) */
  int k_0;                  /* offset into antidiagonal */
  int lb_0, rb_0;           /* left/right bounds of current antidiagonal */
  int lb_1, rb_1;           /* left/right bounds of previous antidiagonal */
  float diag_max, cell_max; /* max score for all normal states in a given cell/antidiagonal */
  float diag_limit = -INF;  /* pruning threshold based on global max */
  float total_limit = -INF; /* pruning threshold based on antidiag max */
  float cell_limit = -INF;  /* pruning threshold based on all thresholds */
  bool is_d_0_in_viterbi;   /* checks if we have gone past the seed viterbi alignment */
  bool is_past_gamma;       /* checks if search has passed all non-pruned starting diagonals */
  bool is_past_limit;       /* checks if any cells score exceeds score pruning threshold */
  bool is_term;             /* checks if search is ended */
  float prv_max;            /* previous maximum to check for updates */
  COORDS coords_diagmax;    /* coordinates of max scoring cell on current antidiagonal */

  /* clear data int vectors (which will be used to create edgebounds) */
  VECTOR_INT_Reuse(lb_vec[0]);
  VECTOR_INT_Reuse(rb_vec[0]);

  /* update maximum score using antidiagonal */
  diag_max = -INF;
  prv_max = -INF;
  coords_diagmax = (COORDS){0, 0};

  for (i = 0; i < lb_vec[1]->N; i++) {
    lb_1 = VEC_X(lb_vec[1], i);
    rb_1 = VEC_X(rb_vec[1], i);

    for (k_0 = lb_1; k_0 < rb_1; k_0++) {
      q_0 = k_0;
      t_0 = d_1 - k_0; /* looking back one diag */
      /* NOTE: Can we presume the only max we care about is in the MATCH state? */
      diag_max = MATH_Max(MATH_Max(diag_max, MMX3(dx1, k_0)),
                          MATH_Max(IMX3(dx1, k_0), DMX3(dx1, k_0)));
      // diag_max = MATH_Max( diag_max, MMX3(dx1, k_0) );

      /* if maximum has increased, then update max cell (we really only care about q_0) */
      if (diag_max > prv_max) {
        coords_diagmax.q_0 = q_0;
        // coords_diagmax.t_0 = t_0;
        prv_max = diag_max;
      }
    }
  }

  // if (d_cnt > 10) exit(0);

  /* Update global_max if new maximum found */
  if (*total_max < diag_max) {
    *total_max = diag_max;
    *coords_max = coords_diagmax;
  }

  /* Set pruning threshold limit based on global maximum */
  total_limit = *total_max - beta;
  /* Set pruning threshold limit based on antidiag maximum */
  diag_limit = diag_max - alpha;
  /* Set score limit based on all pruning parameters */
  cell_limit = MATH_Max(hard_limit,
                        MATH_Max(total_limit, diag_limit));
  /* Check if search is past non-pruning antidiagonals */
  is_past_gamma = (d_cnt > gamma);
  /* Check if any cells exceed maximum pruning threshold */
  is_past_limit = (diag_max >= total_limit);

  /* All edgebounds in previous antidiagonal */
  for (i = 0; i < lb_vec[1]->N; i++) {
    lb_1 = lb_vec[1]->data[i];
    rb_1 = rb_vec[1]->data[i];

    /* If free passes are not complete, skip pruning */
    if (is_past_gamma == false) {
      VECTOR_INT_Pushback(lb_vec[0], lb_1);
      VECTOR_INT_Pushback(rb_vec[0], rb_1);
    }
    elif (is_past_limit == false) /* if no cells pass pruning threshold */
    {
      /* check if diagonal is inside seed viterbi alignment */
      is_d_0_in_viterbi = IS_IN_RANGE(vit_range.beg, vit_range.end, d_0);
      if (is_d_0_in_viterbi == true) /* inside Viterbi seed */
      {
        // DBG_PRINTF("In Viterbi, Below Threshold: (%.3f > %.3f) %d\n", cell_limit, diag_max, coords_diagmax.q_0 );
        /* NOTE: Option - Force search to continue and only select highest scoring cell */
        // VECTOR_INT_Pushback( lb_vec[0], coords_diagmax.q_0 );
        // VECTOR_INT_Pushback( rb_vec[0], coords_diagmax.q_0 );
        // return STATUS_SUCCESS;
        /* NOTE: Option - prune all search space and just terminate? */
        *is_term_flag = true;
        return STATUS_SUCCESS;
      } else /* outside Viterbi seed */
      {
        // fprintf(stdout, "# FLAG TRIGGERED!!! d_cnt=%d, d_0=%d, vit_range={%d,%d}\n", d_cnt, d_0, vit_range.beg, vit_range.end);
        *is_term_flag = true;
        return STATUS_SUCCESS;
      }
    }
    else /* If free passes are complete (gamma < d), prune and set new edgebounds */
    {
      /* impossible state */
      lb_0 = INT_MIN;
      rb_0 = INT_MIN;

      /* Find the first cell from the left which passes above threshold */
      for (k_0 = lb_1; k_0 < rb_1; k_0++) {
        q_0 = k_0;
        t_0 = d_1 - k_0; /* looking back one diag */

        cell_max = MATH_Max(MMX3(dx1, k_0),
                            MATH_Max(IMX3(dx1, k_0), DMX3(dx1, k_0)));
        // cell_max = 	MMX3(dx1, k_0);

        /* prune in left edgebound */
        if (cell_max >= cell_limit) {
          /* when first cell found above limit if found, set boundary and end search */
          lb_0 = k_0;
          VECTOR_INT_Pushback(lb_vec[0], lb_0);
          break;
        }
      }

      /* If no boundary edges are found on diag, then branch is pruned entirely */
      if (lb_0 == INT_MIN) {
        continue;
      }

      /* Find the first cell from the right which passes above threshold */
      for (k_0 = rb_1 - 1; k_0 >= lb_1; k_0--) {
        q_0 = k_0;
        t_0 = d_1 - k_0; /* looking back one diag */

        // cell_max = 	MATH_Max( MMX3(dx1, k_0),
        //             MATH_Max( IMX3(dx1, k_0), DMX3(dx1, k_0) ) );
        cell_max = MMX3(dx1, k_0);

        /* prune in right edgebound */
        if (cell_max >= cell_limit) {
          /* when first cell found above limit if found, set boundary and end search */
          rb_0 = k_0 + 1;
          VECTOR_INT_Pushback(rb_vec[0], rb_0);
          break;
        }
      }
    }
  }
  return STATUS_SUCCESS;
}
