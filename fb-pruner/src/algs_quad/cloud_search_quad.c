/*******************************************************************************
 *  FILE:       cloud_search_quad.c
 *  PURPOSE:    Cloud Search for Forward-Backward Pruning Alg. (QUADRATIC SPACE)
 *
 *  AUTHOR:     Dave Rich
 *  BUG:        Lots.
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
#include "cloud_search_quad.h"

/*
 *      NOTE: HOW TO CONVERT row-coords to diag-coords
 *       MMX3(q_1, t_1) => MMX3(, d_2)
 *       MMX3(q_0, t_1) => MMX3(, d_1)
 *       MMX3(q_1, t_0) => MMX3(, d_1)
 */

/* 
 *       NOTE: CONVERSION - row-coords => diag-coords
 *       MX_M(i-1, j-1) => MX3(d_2, k-1)
 *       MX_M(i  , j-1) => MX3(d_1, k  ) 
 *       MX_M(i-1, j  ) => MX3(d_1, k-1)
 *
 *       MX_M(i+1, j+1) => MX3(d_2, k+1)
 *       MX_M(i  , j+1) => MX3(d_1, k  )
 *       MX_M(i+1, j  ) => MX3(d_1, k+1)
 */

/*  
 *  FUNCTION:  run_Cloud_Forward_Run()
 *  SYNOPSIS:  Perform Forward part of Cloud Search Algorithm.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Cloud_Forward_Quad(   const SEQUENCE*      query,         /* query sequence */
                              const HMM_PROFILE*   target,        /* target hmm model */
                              const int            Q,             /* query length */
                              const int            T,             /* target length */
                              MATRIX_3D*           st_MX,         /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
                              MATRIX_2D*           sp_MX,         /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
                              const ALIGNMENT*     tr,            /* viterbi traceback */ 
                              EDGEBOUND_ROWS*      rows,          /* temporary edgebounds by-row vector */
                              VECTOR_INT*          lb_vec[3],     /* temporary left-bound vectors for pruning */
                              VECTOR_INT*          rb_vec[3],     /* temporary right-bound vectors for pruning */
                              EDGEBOUNDS*          edg,           /* OUTPUT: edgebounds of cloud search space */
                              CLOUD_PARAMS*        params )       /* pruning parameters */
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
   float    prv_loop, prv_move;              /* previous loop and move for special states */
   float    prv_sum, prv_best;               /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* vars for traceback */
   TRACE*         beg;                       /* beginning of the alignment */
   TRACE*         end;                       /* end of the alignment */

   /* vars for pruning */
   float          cell_max, diag_max, total_max;   /* maximum score found in matrix */
   float          total_limit, diag_limit;         /* threshold determined by max_scores - alpha */
   BOUND*         dp_bound;                        /* bounds for dp matrix of current antidiagonal */
   VECTOR_INT*    lb_vec_tmp;                      /* left swap pointer */
   VECTOR_INT*    rb_vec_tmp;                      /* right swap pointer */

   /* pruning parameters */
   float alpha;
   float beta;
   int   gamma;

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

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   /* reuse left and right-bound vectors */
   for ( int i = 0; i < 3; i++ ) {
      VECTOR_INT_Reuse( lb_vec[i] );
      VECTOR_INT_Reuse( rb_vec[i] );
   }

   /* clear all old data from data matrix if necessary */
   if ( st_MX->clean = false ) {
      MATRIX_3D_Clean( st_MX );
   }

   /* query sequence */
   seq         = query->seq;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   /* get pruning parameters */
   alpha = params->alpha;
   beta  = params->beta;
   gamma = params->gamma;

   /* initialize edges and bounds */
   le_0 = 0;
   re_0 = 0;
   lb_0 = lb_1 = lb_2 = 0;
   rb_0 = rb_1 = rb_2 = 0;

   /* set edgebound dimensions and orientation */
   EDGEBOUND_ROWS_Reuse( rows, Q, T );
   EDGEBOUNDS_Reuse( edg, Q, T );
   #if ( CLOUD_METHOD == CLOUD_DIAGS )
   {
      edg->edg_mode  = EDG_DIAG;
   }
   #elif ( CLOUD_METHOD == CLOUD_ROWS )
   {
      edg->edg_mode  = EDG_ROW;
   }
   #endif

   /* malloc dynamic memory */
   for ( i = 0; i < 3; i++ ) {
      lb_vec[i] = VECTOR_INT_Create();
      rb_vec[i] = VECTOR_INT_Create();
   }

   /* query sequence */
   seq = query->seq;

   /* get start and end points of viterbi alignment */
   beg = &(tr->traces->data[tr->beg]);
   end = &(tr->traces->data[tr->end]);

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if ( beg->i == 0 || beg->j == 0 ) {
      beg->i += 1;
      beg->j += 1;
   }

   /* TODO: initialize previous start state */
   q_0 = beg->i;
   q_1 = q_0 - 1;
   t_0 = beg->j;
   t_1 = t_0 - 1;

   /* dimension of submatrix */
   dim_TOT  = Q + T;
   dim_Q    = Q - beg->i;
   dim_T    = T - beg->j;

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = dim_TOT;

   /* diag index of different start points, creating submatrix */
   d_st = beg->i + beg->j;
   // d_end = end->i + end->j;

   /* diag index where num cells reaches highest point and begins diminishing */
   dim_min = MIN( d_st + dim_Q, d_st + dim_T );
   dim_max = MAX( d_st + dim_Q, d_st + dim_T );

   /* set bounds of starting cell */
   lb_0 = beg->i;
   rb_0 = beg->i + 1;
   VECTOR_INT_Pushback( lb_vec[1], lb_0 );
   VECTOR_INT_Pushback( rb_vec[1], rb_0 );
   num_cells = 0;

   /* keeps largest number seen on current diagonal */
   diag_max    = -INF;
   total_max   = -INF;
   /* number of passes through antidiags */
   d_cnt = 0;

   /* begin state probability begins at zero (free to start alignment) */
   prv_B = 0;
   // prv_E = 0;

   /* ITERATE THROUGH ANTI-DIAGONALS */
   for ( d_0 = d_st; d_0 <= d_end; d_0++, d_cnt++ )
   {
      d_1   = d_0 - 1;      /* look back 1 antidiagonal */
      d_2   = d_0 - 2;      /* look back 2 antidiagonal */
      /* map antidiagonals to data matrix */
      dx0   = d_0; 
      dx1   = d_1;
      dx2   = d_2;

      /* Is dp matrix diagonal growing or shrinking? */
      if ( d_0 <= dim_min )
         num_cells++;
      if ( d_0 > dim_max )
         num_cells--;

      /* Edgecheck updates: determine antidiag indices within matrix bounds */
      le_0 = MAX( beg->i, d_0 - T );
      re_0 = le_0 + num_cells;

      /* Prune bounds */
      #if ( PRUNER == PRUNER_XDROP_EDGETRIM )
      {
         /* prune bounds using x-drop, no bifurcating */
         prune_via_xdrop_edgetrim_Quad( 
            st_MX, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_XDROP_BIFURCATE )
      {
         /* prune bounds using x-drop, bifurcating */
         prune_via_xdrop_bifurcate_Quad( 
            st_MX, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_DBL_XDROP_EDGETRIM_OR_DIE )
      {
         /* prune bounds using local and global x-drop, edgetrimming or terminating search */
         prune_diag_by_xdrop_edgetrim_or_die_Quad( 
            st_MX, sp_MX, alpha, beta, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #endif

      /* Add pruned bounds to edgebound list */
      for ( i = 0; i < lb_vec[0]->N; i++ )
      {
         lb_0 = lb_vec[0]->data[i];
         rb_0 = rb_vec[0]->data[i];

         /* Update bounds (spans all cells adjacent to previous antidiagonals cells that were not pruned) */
         lb_0 = lb_0;
         rb_0 = rb_0 + 1;

         /* Update bounds to account for dp matrix bounds */
         lb_0 = MAX(lb_0, le_0);
         rb_0 = MIN(rb_0, re_0);

         /* Update changes to list */
         lb_vec[0]->data[b] = lb_0;
         rb_vec[0]->data[b] = rb_0;

         bnd_new = (BOUND){d_0, lb_0, rb_0};

         #if ( CLOUD_METHOD == CLOUD_DIAGS )
         {
            /* add new bounds to edgebounds as antidiag-wise */
            EDGEBOUNDS_Pushback( edg, &bnd_new );
         }
         #endif

         #if ( CLOUD_METHOD == CLOUD_ROWS )
         {
            /* reorient new bounds from antidiag-wise to row-wise and integrate it into row-wise edgebound list */
            EDGEBOUND_ROWS_Integrate_Antidiag_Fwd( rows, &bnd_new );

            /* add new bounds to edgebounds as antidiag-wise (for comparative testing) */
            #if DEBUG
            {
               EDGEBOUNDS_Pushback( test_edg, &bnd_new );
            }
            #endif
         }
         #endif
      }

      /* If diagonal set is empty, then all branches have been pruned, so we're done */
      if ( lb_vec[0]->N <= 0 ) break;

      /* MAIN RECURSION */
      /* Iterate the ranges of the antidiagonal */
      for ( i = 0; i < lb_vec[0]->N; i++ ) 
      {
         lb_0 = lb_vec[0]->data[i];
         rb_0 = rb_vec[0]->data[i];

         /* Iterate through cells in range */
         for ( k_0 = lb_0; k_0 < rb_0; k_0++ )
         {
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
            prv_M = MMX(q_1, t_1) + TSC(t_1, M2M);
            prv_I = IMX(q_1, t_1) + TSC(t_1, I2M);
            prv_D = DMX(q_1, t_1) + TSC(t_1, D2M);
            /* free to begin match state (new alignment) */
            // prv_B = 0; /* assigned once at start */
            /* best-to-match */
            prv_sum = logsum( 
                           logsum( prv_M, prv_I ),
                           logsum( prv_D, prv_B ) );
            MMX(q_0, t_0) = prv_sum + MSC(t_0, A);

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            prv_M = MMX(q_1, t_0) + TSC(t_0, M2I);
            prv_I = IMX(q_1, t_0) + TSC(t_0, I2I);
            /* best-to-insert */
            prv_sum = logsum( prv_M, prv_I );
            IMX(q_0, t_0) = prv_sum + ISC(t_0, A);

            /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            prv_M = MMX(q_0, t_1) + TSC(t_1, M2D);
            prv_D = DMX(q_0, t_1) + TSC(t_1, D2D);
            /* best-to-delete */
            prv_sum = logsum( prv_M, prv_D );
            DMX(q_0, t_0) = prv_sum;

            /* embed in test matrix */
            #if DEBUG
               MX_2D(cloud_MX, q_0, t_0) = 1.0;
               MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(q_0, t_0);
               MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(q_0, t_0);
               MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(q_0, t_0);
            #endif 
         }
      }

      /* Shift bounds */
      lb_vec_tmp  = lb_vec[2];
      rb_vec_tmp  = rb_vec[2];
      lb_vec[2]   = lb_vec[1];
      rb_vec[2]   = rb_vec[1];
      lb_vec[1]   = lb_vec[0];
      rb_vec[1]   = rb_vec[0];
      lb_vec[0]   = lb_vec_tmp;
      rb_vec[0]   = rb_vec_tmp;

      /* scrub 2-back bound data (linear only) */
      VECTOR_INT_Reuse( lb_vec[0] );
      VECTOR_INT_Reuse( rb_vec[0] );

      /* disallow starting new alignments after first pass */
      prv_B = -INF;
      // prv_E = -INF;
   }
   
   /* show visualization of search cloud */
   #if DEBUG
   {
      DP_MATRIX_VIZ_Trace( cloud_MX, tr );
   }
   #endif

   return STATUS_SUCCESS;
}

/*  
 *  FUNCTION:  run_Cloud_Backward_Run()
 *  SYNOPSIS:  Perform Backward part of Cloud Search Algorithm.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Cloud_Backward_Quad(     const SEQUENCE*      query,         /* query sequence */
                                 const HMM_PROFILE*   target,        /* target hmm model */
                                 const int            Q,             /* query length */
                                 const int            T,             /* target length */
                                 MATRIX_3D*           st_MX,         /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
                                 MATRIX_2D*           sp_MX,         /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
                                 const ALIGNMENT*     tr,            /* viterbi traceback */ 
                                 EDGEBOUND_ROWS*      rows,          /* temporary edgebounds by-row vector */
                                 VECTOR_INT*          lb_vec[3],     /* temporary left-bound vectors for pruning */
                                 VECTOR_INT*          rb_vec[3],     /* temporary right-bound vectors for pruning */
                                 EDGEBOUNDS*          edg,           /* OUTPUT: edgebounds of cloud search space */
                                 CLOUD_PARAMS*        params )       /* pruning parameters */
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
   int      k_0, k_1;                        /* real index offset into diagonals */
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
   float    prv_loop, prv_move;              /* previous loop and move for special states */
   float    prv_sum, prv_best;               /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* vars for traceback */
   TRACE*         beg;                       /* beginning of the alignment */
   TRACE*         end;                       /* end of the alignment */

   /* vars for pruning */
   float          cell_max, diag_max, total_max;   /* maximum score found in matrix */
   float          total_limit, diag_limit;         /* threshold determined by max_scores - alpha */
   BOUND*         dp_bound;                        /* bounds for dp matrix of current antidiagonal */
   VECTOR_INT*    lb_vec_tmp;                      /* left swap pointer */
   VECTOR_INT*    rb_vec_tmp;                      /* right swap pointer */

   /* pruning parameters */
   float alpha;
   float beta;
   int   gamma;

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

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   /* reuse left and right-bound vectors */
   for ( int i = 0; i < 3; i++ ) {
      VECTOR_INT_Reuse( lb_vec[i] );
      VECTOR_INT_Reuse( rb_vec[i] );
   }

   /* clear all old data from data matrix if necessary */
   if ( st_MX->clean = false ) {
      MATRIX_3D_Clean( st_MX );
   }

   /* query sequence */
   seq         = query->seq;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   /* get pruning parameters */
   alpha = params->alpha;
   beta  = params->beta;
   gamma = params->gamma;

   /* initialize edges and bounds */
   le_0 = 0;
   re_0 = 0;
   lb_0 = lb_1 = lb_2 = 0;
   rb_0 = rb_1 = rb_2 = 0;

   /* set edgebound dimensions and orientation */
   EDGEBOUND_ROWS_Reuse( rows, Q, T );
   EDGEBOUNDS_Reuse( edg, Q, T );
   #if ( CLOUD_METHOD == CLOUD_DIAGS )
   {
      edg->edg_mode  = EDG_DIAG;
   }
   #elif ( CLOUD_METHOD == CLOUD_ROWS )
   {
      edg->edg_mode  = EDG_ROW;
   }
   #endif

   /* query sequence */
   seq = query->seq;

   /* get start and end points of viterbi alignment */
   beg = &(tr->traces->data[tr->beg]);
   end = &(tr->traces->data[tr->end]);

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (end->i == Q || end->j == T) {
      end->i -= 1;
      end->j -= 1;
   }

   /* TODO: initialize previous start state */
   q_0 = end->i;
   q_1 = q_0 + 1;
   t_0 = end->j;
   t_1 = t_0 + 1;

   /* dimension of submatrix */
   dim_TOT  = Q + T;
   dim_Q    = end->i;
   dim_T    = end->j;

   /* diag index at corners of dp matrix */
   d_st     = 0;
   d_end    = dim_TOT;

   /* diag index of different start points, creating submatrix */
   // d_st = beg->i + beg->j;
   d_end = end->i + end->j;

   /* diag index where num cells reaches highest point and begins diminishing */
   dim_min  = MIN(dim_T, dim_Q);
   dim_max  = MAX(dim_T, dim_Q);

   /* set bounds of starting cell */
   lb_0 = end->i;
   rb_0 = end->i + 1;
   VECTOR_INT_Pushback( lb_vec[1], lb_0 );
   VECTOR_INT_Pushback( rb_vec[1], rb_0 );
   num_cells = 0;

   /* keeps largest number seen on current diagonal */
   diag_max    = -INF;
   total_max   = -INF;
   /* number of antidiags passed through */
   d_cnt    = 0;

   /* begin state probability begins at zero (free to start alignment) */
   // prv_B = 0;
   prv_E = 0;

   /* ITERATE THROUGHT ANTI-DIAGONALS */
   for ( d_0 = d_end; d_0 >= d_st; d_0--, d_cnt++ )
   {
      d_1 = d_0 + 1;    /* look back 1 diagonal */
      d_2 = d_0 + 2;    /* look back 2 diagonals */
      /* mod-mapping of antidiagonals into linear space */
      dx0 = d_0 % 3; 
      dx1 = d_1 % 3;
      dx2 = d_2 % 3;

      /* Is dp matrix diagonal growing or shrinking? */
      if ( d_0 >= dim_max )
         num_cells++;
      if ( d_0 < dim_min )
         num_cells--;

      /* TODO: is there a closed form for edges, aka not using num_cells? */
      /* Edgechecks: find antidiag cell range that is inside matrix bounds */
      le_0 = MAX(end->i - (d_end - d_0), 0);
      re_0 = le_0 + num_cells;

      /* Prune bounds */
      #if ( PRUNER == PRUNER_XDROP_EDGETRIM )
      {
         /* prune bounds using x-drop, no bifurcating */
         prune_via_xdrop_edgetrim_Quad( 
            st_MX, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_XDROP_BIFURCATE )
      {
         /* prune bounds using x-drop, bifurcating */
         prune_via_xdrop_bifurcate_Quad( 
            st_MX, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_DBL_XDROP_EDGETRIM_OR_DIE )
      {
         /* prune bounds using local and global x-drop, edgetrimming or terminating search */
         prune_diag_by_xdrop_edgetrim_or_die_Quad( 
            st_MX, sp_MX, alpha, beta, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #endif

      /* Add pruned bounds to edgebound list */
      for ( i = 0; i < lb_vec[0]->N; i++ )
      {
         lb_0 = lb_vec[0]->data[i];
         rb_0 = rb_vec[0]->data[i];

         /* Update bounds (spans all cells adjacent to previous antidiagonals cells that were not pruned) */
         lb_0 = lb_0 - 1;
         rb_0 = rb_0;

         /* Update bounds to account for dp matrix bounds */
         lb_0 = MAX(lb_0, le_0);
         rb_0 = MIN(rb_0, re_0);

         /* Update changes to list */
         lb_vec[0]->data[b] = lb_0;
         rb_vec[0]->data[b] = rb_0;

         bnd_new = (BOUND){d_0, lb_0, rb_0};

         #if ( CLOUD_METHOD == CLOUD_DIAGS )
         {
            /* add new bounds to edgebounds as antidiag-wise */
            EDGEBOUNDS_Pushback( edg, &bnd_new );
         }
         #endif

         #if ( CLOUD_METHOD == CLOUD_ROWS )
         {
            /* reorient new bounds from antidiag-wise to row-wise and integrate it into row-wise edgebound list */
            EDGEBOUND_ROWS_Integrate_Antidiag_Fwd( rows, &bnd_new );

            /* add new bounds to edgebounds as antidiag-wise (for comparative testing) */
            #if DEBUG
            {
               EDGEBOUNDS_Pushback( test_edg, &bnd_new );
            }
            #endif
         }
         #endif
      }

      /* If diagonal set is empty, then all branches have been pruned, so we're done */
      if ( lb_vec[0]->N <= 0 ) break;

      /* MAIN RECURSION */
      for ( i = 0; i < lb_vec[0]->N; i++ )
      {
         lb_0 = lb_vec[0]->data[i];
         rb_0 = rb_vec[0]->data[i];

         /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
         for ( k_0 = lb_0; k_0 < rb_0; k_0++ )
         {
            /* get row-col coords */
            q_0 = k_0;
            q_1 = q_0 + 1;
            t_0 = d_0 - k_0;
            t_1 = t_0 + 1;

            /* next sequence character */
            a = seq[q_0];
            A = AA_REV[a];
            
            /* match and insertion scores */
            sc_M = MSC(t_1, A);
            sc_I = ISC(t_1, A);

            /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
            prv_M = MMX(q_1, t_1) + TSC(t_0, M2M) + sc_M;
            prv_I = IMX(q_1, t_0) + TSC(t_0, M2I) + sc_I;
            prv_D = DMX(q_0, t_1) + TSC(t_0, M2D);
            // prv_E = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
            // prv_E = sc_E;
            /* best-to-match */
            prv_sum = logsum( 
                           logsum( prv_M, prv_I ),
                           logsum( prv_D, prv_E ) );
            MMX(q_0, t_0) = prv_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            prv_M = MMX(q_1, t_1) + TSC(t_0, I2M) + sc_M;
            prv_I = IMX(q_1, t_0) + TSC(t_0, I2I) + sc_I;
            /* best-to-insert */
            prv_sum = logsum( prv_M, prv_I );
            IMX(q_0, t_0) = prv_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            prv_M = MMX(q_1, t_1) + TSC(t_0, D2M) + sc_M;
            prv_D = DMX(q_0, t_1) + TSC(t_0, D2D);
            /* best-to-delete */
            prv_sum = logsum( prv_M, prv_D );
            prv_sum = logsum( prv_sum, prv_E );
            DMX(q_0, t_0) = prv_sum;

            /* embed cell data in quadratic matrix */
            #if DEBUG
            {
               MX_2D( cloud_MX, q_0, t_0 ) = 1.0;
               MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX(q_0, t_0);
               MX_3D(test_MX, INS_ST, q_0, t_0) = IMX(q_0, t_0);
               MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX(q_0, t_0);
            }
            #endif 
         }
      }

      /* Shift bounds */
      lb_vec_tmp  = lb_vec[2];
      rb_vec_tmp  = rb_vec[2];
      lb_vec[2]   = lb_vec[1];
      rb_vec[2]   = rb_vec[1];
      lb_vec[1]   = lb_vec[0];
      rb_vec[1]   = rb_vec[0];
      lb_vec[0]   = lb_vec_tmp;
      rb_vec[0]   = rb_vec_tmp;

      /* scrub 2-back bound data (linear only) */
      VECTOR_INT_Reuse( lb_vec[0] );
      VECTOR_INT_Reuse( rb_vec[0] );

      /* disallow starting new alignments after first pass */
      // prv_B = -INF;
      prv_E = -INF;
   }

   /* reverse order of diagonals */
   EDGEBOUNDS_Reverse(edg);

   /* show visualization of search cloud */
   #if DEBUG
   {
      DP_MATRIX_VIZ_Trace( cloud_MX, tr );
   }
   #endif 

   return total_max;
} 