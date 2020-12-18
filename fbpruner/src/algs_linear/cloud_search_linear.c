/*******************************************************************************
 *    FILE:       cloud_search_linear.c
 *    PURPOSE:    Cloud Search for Forward-Backward Pruning Alg.
 *                (Linear Space Alg)
 *
 *    AUTHOR:     Dave Rich
 *    BUG:        No known.
 *    TODO:       lb_vec and rb_vec should be created in main routine so we don't need to create/destroy every routine.       
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
#include "../utilities/utilities.h"
#include "../objects/objects.h"

#include "algs_linear.h"

/* for debugging */
#include "../algs_quad/algs_quad.h"

/* header */
#include "cloud_search_linear.h"

/*
 *      NOTE: HOW TO CONVERT row-coords to diag-coords
 *       MMX3(q_1, t_1) => MMX3(d_2, k_1)
 *       MMX3(q_0, t_1) => MMX3(d_1, k_0)
 *       MMX3(q_1, t_0) => MMX3(d_1, k_1)
 */

/* 
 *       NOTE: CONVERSION - row-coords => diag-coords
 *       MX_M(i-1, j-1) => MX3_M(d_2, k-1)
 *       MX_M(i  , j-1) => MX3_M(d_1, k  ) 
 *       MX_M(i-1, j  ) => MX3_M(d_1, k-1)
 *
 *       MX_M(i+1, j+1) => MX3_M(d_2, k+1)
 *       MX_M(i  , j+1) => MX3_M(d_1, k  )
 *       MX_M(i+1, j  ) => MX3_M(d_1, k+1)
 */

/*
 *  FUNCTION: run_Cloud_Forward_Linear()
 *  SYNOPSIS: Perform Forward part of Cloud Search Algorithm.
 *            Traverses the dynamic programming matrix antidiagonally, running the
 *            Forward algorithm, starting at the Viterbi alignment beginning.  
 *            At the end of each antidiagonal, compares each cell against the current maximum
 *            scoring cell.  If cell falls below (MAX_SCORE * alpha), then cell is removed 
 *            from search space.  Terminates when reaches the final cell in dp matrix or 
 *            all cells in current antidiag have been pruned.  
 *            Stores final edgebound data in <edg>.
 *  RETURN:   Returns <STATUS_SUCCESS> if no errors.
 */
int run_Cloud_Forward_Linear(    const SEQUENCE*      query,        /* query sequence */
                                 const HMM_PROFILE*   target,       /* target hmm model */
                                 const int            Q,            /* query length */
                                 const int            T,            /* target length */
                                 MATRIX_3D*           st_MX3,       /* normal state matrix */
                                 MATRIX_2D*           sp_MX,        /* special state matrix */
                                 const ALIGNMENT*     tr,           /* viterbi traceback */
                                 EDGEBOUND_ROWS*      rows,         /* temporary edgebounds by-row vector */
                                 EDGEBOUNDS*          edg,          /* OUTPUT: edgebounds of cloud search space */
                                 CLOUD_PARAMS*        params )      /* pruning parameters */
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
   int      d_st, d_end, d_cnt, d_last;      /* starting and ending diagonal indices */
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
   TRACE*         beg;                             /* beginning of the alignment */
   TRACE*         end;                             /* end of the alignment */

   /* vars for pruning */
   bool           is_term_flag;                    /* termination flag for end of search */
   float          cell_max, diag_max, total_max;   /* maximum score found in matrix */
   float          total_limit, diag_limit;         /* threshold determined by max_scores - alpha */
   BOUND*         dp_bound;                        /* bounds for dp matrix of current antidiagonal */
   VECTOR_INT*    lb_vec[3];                       /* left bound list for previous 3 antdiags */
   VECTOR_INT*    rb_vec[3];                       /* right bound list for previous 3 antidiags */
   VECTOR_INT*    lb_vec_tmp;                      /* left swap pointer */
   VECTOR_INT*    rb_vec_tmp;                      /* right swap pointer */

   /* pruning parameters */
   float       alpha;
   float       beta;
   int         gamma;

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

   /* check if data is cleaned */
   #if DEBUG
   {
      int cmp = MATRIX_3D_Check_Clean( st_MX3 );
      printf("PRE-CHECK CLEAN -> CLOUD FWD?\t%d\n", cmp );
   }
   #endif 

   /* clear all old data from data matrix if necessary */
   if ( st_MX3->clean = false ) {
      MATRIX_3D_Clean( st_MX3 );
   }

   /* query sequence */
   seq         = query->seq;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   /* get pruning parameters */
   alpha       = params->alpha;
   beta        = params->beta;
   gamma       = params->gamma;

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

   /* TODO: add this edgecheck to all  */
   /* verify that starting points are valid */
   if ( beg->q_0 < 0 || beg->q_0 > Q || beg->t_0 < 0 || beg->t_0 > T ) {
      fprintf(stderr, "# ERROR: Invalid start points for Cloud Forward Search: (%d,%d)\n", beg->q_0, beg->t_0 );
      fprintf(stderr, "# Query Length: %d, Target Length: %d\n", Q, T );
      // exit(EXIT_FAILURE);
   }
   if ( beg->q_0 > Q ) {
      beg->q_0 = Q;
   }
   if ( beg->t_0 > T ) {
      beg->t_0 = T;
   }

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if ( beg->q_0 == 0 || beg->t_0 == 0 ) {
      beg->q_0 += 1;
      beg->t_0 += 1;
   }

   /* TODO: initialize previous start state */
   q_0 = beg->q_0;
   q_1 = q_0 - 1;
   t_0 = beg->t_0;
   t_1 = t_0 - 1;

   /* dimension of submatrix */
   dim_TOT  = Q + T;
   dim_Q    = Q - beg->q_0;
   dim_T    = T - beg->t_0;

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = dim_TOT;

   /* diag index of different start points, creating submatrix */
   d_st = beg->q_0 + beg->t_0;
   // d_end = end->q_0 + end->t_0;

   /* diag index where num cells reaches highest point and begins diminishing */
   dim_min = MIN( d_st + dim_Q, d_st + dim_T );
   dim_max = MAX( d_st + dim_Q, d_st + dim_T );

   /* set bounds of starting cell */
   lb_0 = beg->q_0;
   rb_0 = beg->q_0;
   VECTOR_INT_Pushback( lb_vec[1], lb_0 );
   VECTOR_INT_Pushback( rb_vec[1], rb_0 );
   num_cells = 0;

   /* keeps largest number seen on current diagonal */
   is_term_flag = false;
   diag_max    = -INF;
   total_max   = -INF;
   /* number of passes through antidiags */
   d_cnt = 0;

   /* begin state probability begins at zero (free to start alignment) */
   prv_B = 0;
   // prv_E = 0;

   /* ITERATE THROUGH ANTI-DIAGONALS */
   for (d_0 = d_st; d_0 <= d_end; d_0++, d_cnt++)
   {
      d_1   = d_0 - 1;      /* look back 1 antidiagonal */
      d_2   = d_0 - 2;      /* look back 2 antidiagonal */
      /* mod-mapping of antidiagonals into linear space */
      dx0   = d_0 % 3; 
      dx1   = d_1 % 3;
      dx2   = d_2 % 3;

      /* is dp matrix diagonal growing or shrinking? */
      if ( d_0 <= dim_min )
         num_cells++;
      if ( d_0 > dim_max )
         num_cells--;

      /* Edgecheck updates: determine antidiag indices within matrix bounds */
      le_0 = MAX( beg->q_0, d_0 - T );
      re_0 = le_0 + num_cells;

      /* Macro-controlled - Bounds pruning method */
      #if ( PRUNER == PRUNER_XDROP_EDGETRIM )
      {
         /* prune bounds using x-drop, no bifurcating */
         prune_via_xdrop_edgetrim_Linear( 
            st_MX3, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_XDROP_BIFURCATE )
      {
         /* prune bounds using x-drop, bifurcating */
         prune_via_xdrop_bifurcate_Linear( 
            st_MX3, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_DBL_XDROP_EDGETRIM_OR_DIE )
      {
         /* prune bounds using local and global x-drop, edgetrimming or terminating search */
         prune_via_dbl_xdrop_edgetrim_or_die_Linear( 
            st_MX3, sp_MX, alpha, beta, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, &is_term_flag, lb_vec, rb_vec );
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

         /* macro-controlled: how to store cloud -> antidiag or row-wise */
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
            /* NOTE: Convert (q-1,t-1) <=> (d-2,k-1) */ 
            prv_M = MMX3(dx2, k_1)  + TSC(t_1, M2M);
            prv_I = IMX3(dx2, k_1)  + TSC(t_1, I2M);
            prv_D = DMX3(dx2, k_1)  + TSC(t_1, D2M);
            /* Free to begin match state (new alignment) */
            /* NOTE: only allow begin transition at start of viterbi alignment */
            // prv_B = 0; /* assigned once at start */
            /* best-to-match */
            prv_sum = logsum( 
                           logsum( prv_M, prv_I ),
                           logsum( prv_D, prv_B ) );
            MMX3(dx0, k_0) = prv_sum + MSC(t_0, A);

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            /* NOTE: Convert (q-1,t) <=> (d-1,k-1) */
            prv_M = MMX3(dx1, k_1) + TSC(t_0, M2I);
            prv_I = IMX3(dx1, k_1) + TSC(t_0, I2I);
            /* best-to-insert */
            prv_sum = logsum( prv_M, prv_I );
            IMX3(dx0, k_0) = prv_sum + ISC(t_0, A);

            /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            /* NOTE: Convert (q,t-1) <=> (d-1, k) */
            prv_M = MMX3(dx1, k_0) + TSC(t_1, M2D);
            prv_D = DMX3(dx1, k_0) + TSC(t_1, D2D);
            /* best-to-delete */
            prv_sum = logsum( prv_M, prv_D );
            DMX3(dx0, k_0) = prv_sum;

            /* embed cell data in quadratic matrix */
            #if DEBUG
            {
               MX_2D( cloud_MX, q_0, t_0 ) += 1.0;
               MX_2D( cloud_MX3, dx0, k_0 ) += 1.0;

               MX_3D( test_MX, MAT_ST, q_0, t_0 ) = MMX3(dx0, k_0);
               MX_3D( test_MX, INS_ST, q_0, t_0 ) = IMX3(dx0, k_0);
               MX_3D( test_MX, DEL_ST, q_0, t_0 ) = DMX3(dx0, k_0);
            }
            #endif 
         }
      }

      /* Scrub 2-back bound data */
      for ( i = 0; i < lb_vec[2]->N; i++ )
      {
         lb_2 = lb_vec[2]->data[i];
         rb_2 = rb_vec[2]->data[i];

         for ( k_0 = lb_2; k_0 < rb_2; k_0++ ) 
         {
            q_0 = k_0;
            t_0 = d_2 - k_0;
            MMX3(dx2, k_0) = IMX3(dx2, k_0) = DMX3(dx2, k_0) = -INF;
            #if DEBUG 
            {
               MX_2D( cloud_MX, q_0, t_0 ) += 2.0;
               MX_2D( cloud_MX3, dx2, k_0 ) -= 1.0;
            }
            #endif
         }
      }

      /* check that all necessary cells have been cleared */
      #if MEMCHECK
      {
         bool is_clean = false;
         
         for (int k_0 = 0; k_0 < (Q+1)+(T+1); k_0++) 
         {
            is_clean = false;
            is_clean += (( MMX3(dx2, k_0) == -INF ) == false);
            is_clean += (( MMX3(dx2, k_0) == -INF ) == false);
            is_clean += (( DMX3(dx2, k_0) == -INF ) == false);
            if ( is_clean != 0 ) {
               ERRORCHECK_memcheck( d_2, k_0, MMX3(dx2, k_0), IMX3(dx2, k_0), DMX3(dx2, k_0) );
               MMX3(dx2, k_0) = IMX3(dx2, k_0) = DMX3(dx2, k_0) = -INF;
            }
         }
      }
      #endif 

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

      /* if termination flag is set, break out of loop */
      if (is_term_flag == true) {
         break;
      }
   }

   /* Scrub last two rows */
   d_last = d_0;
   for (d_0 = d_last; d_0 < d_last + 2; d_0++)
   {
      d_1 = d_0 - 1;      /* look back 1 antidiagonal */
      d_2 = d_0 - 2;      /* look back 2 antidiagonal */
      /* map antidiagonals to data matrix */
      dx0  = d_0 % 3; 
      dx1  = d_1 % 3;
      dx2  = d_2 % 3;

      /* Scrub 2-back bound data */
      for ( i = 0; i < lb_vec[2]->N; i++ )
      {
         lb_2 = lb_vec[2]->data[i];
         rb_2 = rb_vec[2]->data[i];

         for ( k_0 = lb_2; k_0 < rb_2; k_0++ ) 
         {
            q_0 = k_0;
            t_0 = d_2 - k_0;
            MMX3(dx2, k_0) = IMX3(dx2, k_0) = DMX3(dx2, k_0) = -INF;
            #if DEBUG
            {
               MX_2D( cloud_MX, q_0, t_0 ) += 2.0;
               MX_2D( cloud_MX3, dx2, k_0 ) -= 1.0;
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
   }

   /* check that all necessary cells have been cleared */
   #if MEMCHECK
   {
      bool is_clean = false;

      for (dx0 = 0; dx0 < 3; dx0++)
      {
         for (int k_0 = 0; k_0 < (Q+1)+(T+1); k_0++) 
         {
            is_clean = false;
            is_clean += (( MMX3(dx0, k_0) == -INF ) == false);
            is_clean += (( MMX3(dx0, k_0) == -INF ) == false);
            is_clean += (( DMX3(dx0, k_0) == -INF ) == false);
            if ( is_clean != 0 ) {
               ERRORCHECK_memcheck( dx0, k_0, MMX3(dx0, k_0), IMX3(dx0, k_0), DMX3(dx0, k_0) );
               MMX3(dx0, k_0) = IMX3(dx0, k_0) = DMX3(dx0, k_0) = -INF;
            }
         }
      }
   }
   #endif 

   /* if storing row-wise, compare to diag-wise */
   #if ( CLOUD_METHOD == CLOUD_ROWS ) 
   {
      /* output rows to edgebounds */
      EDGEBOUND_ROWS_Convert( rows, edg );

      #if DEBUG
      {
         /* compare cloud rows method to antidiagonal method */
         int cmp = EDGEBOUNDS_Compare_by_Cloud_Single( cloud_MX, edg, test_edg );
         printf("COMPARE (rows vs antidiag):\t%s\n", (cmp == 0) ? "PASS" : "FAIL");
         if ( cmp != 0 ) {
            EDGEBOUNDS_Dump( edg, stdout );
            EDGEBOUNDS_Dump( test_edg, stdout );

            printf("=== ROW-WISE ===\n");
            MATRIX_2D_Fill( cloud_MX, 0 );
            MATRIX_2D_Cloud_Fill( cloud_MX, edg, 1 );
            DP_MATRIX_VIZ_Dump( cloud_MX, stdout );

            printf("=== ROW-WISE ===\n");
            MATRIX_2D_Fill( cloud_MX, 0 );
            MATRIX_2D_Cloud_Fill( cloud_MX, test_edg, 1 );
            DP_MATRIX_VIZ_Dump( cloud_MX, stdout );

            printf("=== OVERLAY ===\n");
            DP_MATRIX_VIZ_Compare( cloud_MX, edg, test_edg );
            DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
         }
      }
      #endif
   }
   #endif
   
   /* check if data is cleaned */
   #if DEBUG
   {
      DP_MATRIX_Save(Q, T, test_MX, sp_MX, "test_output/cloud_search_fwd.mx");

      int cmp = MATRIX_3D_Check_Clean( st_MX3 );
      printf("POST-CHECK CLEAN -> CLOUD FWD?\t%d\n", cmp );
      printf("MAX CLOUD_FWD SCORE: %f, LIMIT: %f\n", total_max, total_limit);
   }
   #endif 

   /* free dynamic memory ( TODO: create/destroy vector in main routine. ) */
   for (i = 0; i < 3; i++) {
      VECTOR_INT_Destroy( lb_vec[i] );
      VECTOR_INT_Destroy( rb_vec[i] );
   }

   /* after search, all cells are set to -INF */
   st_MX3->clean = true;

   return STATUS_SUCCESS;
}


/*
 *  FUNCTION: run_Cloud_Backward_Linear()
 *  SYNOPSIS: Perform Backward part of Cloud Search Algorithm.
 *            Traverses the dynamic programming matrix antidiagonally, running the
 *            Forward algorithm, starting at the Viterbi alignment ending.  
 *            At the end of each antidiagonal, compares each cell against the current maximum
 *            scoring cell.  If cell falls below (MAX_SCORE * alpha), then cell is removed 
 *            from search space.  Terminates when reaches the final cell in dp matrix or 
 *            all cells in current antidiag have been pruned.  
 *            Stores final edgebound data in <edg>.
 *  RETURN:   Maximum score.
 */
int run_Cloud_Backward_Linear(   const SEQUENCE*      query,        /* query sequence */
                                 const HMM_PROFILE*   target,       /* target hmm model */
                                 const int            Q,            /* query length */
                                 const int            T,            /* target length */
                                 MATRIX_3D*           st_MX3,       /* normal state matrix */
                                 MATRIX_2D*           sp_MX,        /* special state matrix */
                                 const ALIGNMENT*     tr,           /* viterbi traceback */
                                 EDGEBOUND_ROWS*      rows,         /* temporary edgebounds by-row */
                                 EDGEBOUNDS*          edg,          /* (OUTPUT) */
                                 CLOUD_PARAMS*        params )      /* pruning parameters */
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
   int      d_st, d_end, d_cnt, d_last;      /* starting and ending diagonal indices */
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
   bool           is_term_flag;                    /* termination flag set by  */
   float          cell_max, diag_max, total_max;   /* maximum score found in matrix */
   float          total_limit, diag_limit;         /* threshold determined by max_scores - alpha */
   BOUND*         dp_bound;                        /* bounds for dp matrix of current antidiagonal */
   VECTOR_INT*    lb_vec[3];                       /* left bound list for previous 3 antdiags */
   VECTOR_INT*    rb_vec[3];                       /* right bound list for previous 3 antidiags */
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
   EDGEBOUNDS* test_edg;

   int         num_writes;
   int         num_clears;

   /* initialize debugging matrix */
   #if DEBUG
   {
      cloud_MX    = debugger->cloud_MX;
      cloud_MX3   = debugger->cloud_MX3;
      test_MX     = debugger->test_MX;
      test_MX3    = debugger->test_MX3;
      test_edg    = debugger->test_edg;

      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
      MATRIX_2D_Fill( cloud_MX3, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
      MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
      MATRIX_3D_Fill( test_MX3, -INF );
      EDGEBOUNDS_Reuse( test_edg, Q, T );

      num_writes = 0;
      num_clears = 0;
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   /* check if data is cleaned */
   #if DEBUG
   {
      int cmp = MATRIX_3D_Check_Clean( st_MX3 );
      printf("PRE-CHECK CLEAN -> CLOUD BCK?\t%d\n", cmp );
   }
   #endif 

   /* clear all old data from data matrix if necessary */
   if ( st_MX3->clean = false ) {
      MATRIX_3D_Clean( st_MX3 );
   }

   /* query sequence */
   seq         = query->seq;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   /* get pruning parameters */
   alpha       = params->alpha;
   beta        = params->beta;
   gamma       = params->gamma;

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
   dp_bound = (BOUND*) malloc( sizeof(BOUND) );
   for ( i = 0; i < 3; i++ ) {
      lb_vec[i] = VECTOR_INT_Create();
      rb_vec[i] = VECTOR_INT_Create();
   }

   /* query sequence */
   seq = query->seq;

   /* get start and end points of viterbi alignment */
   beg = &(tr->traces->data[tr->beg]);
   end = &(tr->traces->data[tr->end]);

   /* TODO: add this edgecheck to all  */
   /* verify that starting points are valid */
   // printf("q_range = {%d,%d}/%d || t_range = {%d,%d}/%d\n");
   if ( end->q_0 < 0 || end->q_0 > Q || end->t_0 < 0 || end->t_0 > T ) {
      fprintf(stderr, "# ERROR: Invalid start points for Cloud Forward Search: (%d,%d)\n", beg->q_0, beg->t_0 );
      fprintf(stderr, "# Query Length: %d, Target Length: %d\n", Q, T );
      // ERRORCHECK_exit(EXIT_FAILURE);
   }
   /* temporary override: this should be handled in pipeline */
   if ( end->q_0 > Q ) {
      end->q_0 = Q;
   }
   if ( end->t_0 > T ) {
      end->t_0 = T;
   }

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if ( end->q_0 == Q || end->t_0 == T ) {
      end->q_0 -= 1;
      end->t_0 -= 1;
   }

   /* TODO: initialize previous start state */
   q_0 = end->q_0;
   q_1 = q_0 + 1;
   t_0 = end->t_0;
   t_1 = t_0 + 1;

   /* dimension of submatrix */
   dim_TOT  = Q + T;
   dim_Q    = end->q_0;
   dim_T    = end->t_0;

   /* diag index at corners of dp matrix */
   d_st     = 0;
   d_end    = dim_TOT;

   /* diag index of different start points, creating submatrix */
   // d_st = beg->q_0 + beg->t_0;
   d_end = end->q_0 + end->t_0;

   /* diag index where num cells reaches highest point and begins diminishing */
   dim_min  = MIN(dim_T, dim_Q);
   dim_max  = MAX(dim_T, dim_Q);

   /* set bounds of starting cell */
   lb_0 = end->q_0;
   rb_0 = end->q_0 + 1;
   VECTOR_INT_Pushback( lb_vec[1], lb_0 );
   VECTOR_INT_Pushback( rb_vec[1], rb_0 );
   num_cells = 0;

   /* keeps largest number seen on current diagonal */
   is_term_flag = false;
   diag_max    = -INF;
   total_max   = -INF;
   /* number of antidiags passed through */
   d_cnt    = 0;

   /* begin state probability begins at zero (free to start alignment) */
   // prv_B = 0;
   prv_E = 0;

   /* ITERATE THROUGHT ANTI-DIAGONALS */
   for (d_0 = d_end; d_0 >= d_st; d_0--, d_cnt++)
   {
      d_1 = d_0 + 1;   /* look back 1 diagonal */
      d_2 = d_0 + 2;   /* look back 2 diagonals */
      /* mod-mapping of antidiagonals into linear space */
      dx0  = d_0 % 3; 
      dx1  = d_1 % 3;
      dx2  = d_2 % 3;

      /* Is dp matrix diagonal growing or shrinking? */
      if (d_0 >= dim_max)
         num_cells++;
      if (d_0 < dim_min)
         num_cells--;

      /* Edgecheck updates: determine antidiag indices within matrix bounds */
      le_0 = MAX( end->q_0 - (d_end - d_0), 0 );
      re_0 = le_0 + num_cells;

      /* Prune bounds */
      #if ( PRUNER == PRUNER_XDROP_EDGETRIM )
      {
         /* prune bounds using x-drop, no bifurcating */
         prune_via_xdrop_edgetrim_Linear( 
            st_MX3, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_XDROP_BIFURCATE )
      {
         /* prune bounds using x-drop, bifurcating */
         prune_via_xdrop_bifurcate_Linear( 
            st_MX3, sp_MX, alpha, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_DBL_XDROP_EDGETRIM_OR_DIE )
      {
         /* prune bounds using local and global x-drop, edgetrimming or terminating search */
         prune_via_dbl_xdrop_edgetrim_or_die_Linear( 
            st_MX3, sp_MX, alpha, beta, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, &is_term_flag, lb_vec, rb_vec );
      }
      #endif

      /* Add pruned bounds to edgebound list */
      for ( i = 0; i < lb_vec[0]->N; i++ )
      {
         /* pull bounds from list */
         lb_0 = lb_vec[0]->data[i];
         rb_0 = rb_vec[0]->data[i];

         /* Update bounds (spans all cells adjacent to previous antidiagonals cells that were not pruned) */
         lb_0 = lb_0 - 1;
         rb_0 = rb_0;

         /* Update bounds to account for dp matrix bounds */
         lb_0 = MAX(lb_0, le_0);
         rb_0 = MIN(rb_0, re_0);

         /* Update changes to list */
         lb_vec[0]->data[i] = lb_0;
         rb_vec[0]->data[i] = rb_0;

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
            EDGEBOUND_ROWS_Integrate_Antidiag_Bck( rows, &bnd_new );

            #if DEBUG
            {
               /* add new bounds to edgebounds as antidiag-wise (for comparative testing) */
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
            k_1 = k_0 + 1;

            /* get x-y coords */
            q_0 = k_0;
            q_1 = k_0 + 1;
            t_0 = d_0 - k_0;
            t_1 = t_0 + 1;

            /*    
             *    === ROW-WISE to DIAG_WISE ===
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
            prv_M = MMX3(dx2, k_1) + TSC(t_0, M2M) + sc_M;
            prv_I = IMX3(dx1, k_1) + TSC(t_0, M2I) + sc_I;
            prv_D = DMX3(dx1, k_0) + TSC(t_0, M2D);
            // prv_E = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
            // prv_E = sc_E;
            /* best-to-match */
            prv_sum = logsum( 
                           logsum( prv_M, prv_I ),
                           logsum( prv_D, prv_E ) );
            MMX3(dx0, k_0) = prv_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            prv_M = MMX3(dx2, k_1) + TSC(t_0, I2M) + sc_M;
            prv_I = IMX3(dx1, k_1) + TSC(t_0, I2I) + sc_I;
            /* best-to-insert */
            prv_sum = logsum( prv_M, prv_I );
            IMX3(dx0, k_0) = prv_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            prv_M = MMX3(dx2, k_1) + TSC(t_0, D2M) + sc_M;
            prv_D = DMX3(dx1, k_0) + TSC(t_0, D2D);
            /* best-to-delete */
            prv_sum = logsum( prv_M, prv_D );
            prv_sum = logsum( prv_sum, prv_E );
            DMX3(dx0, k_0) = prv_sum;

            /* embed cell data in quadratic matrix */
            #if DEBUG
            {
               MX_2D( cloud_MX, q_0, t_0 )  += 1.0;
               MX_2D( cloud_MX3, dx0, k_0 ) += 1.0;

               MX_3D( test_MX, MAT_ST, q_0, t_0 ) = MMX3(dx0, k_0);
               MX_3D( test_MX, INS_ST, q_0, t_0 ) = IMX3(dx0, k_0);
               MX_3D( test_MX, DEL_ST, q_0, t_0 ) = DMX3(dx0, k_0);
            }
            #endif 
         }
      }

      /* Scrub 2-back bound data */
      for ( i = 0; i < lb_vec[2]->N; i++ )
      {
         lb_2 = lb_vec[2]->data[i];
         rb_2 = rb_vec[2]->data[i];

         for ( k_0 = lb_2; k_0 < rb_2; k_0++ ) 
         {
            q_0 = k_0;
            t_0 = d_2 - q_0;

            MMX3(dx2, k_0) = IMX3(dx2, k_0) = DMX3(dx2, k_0) = -INF;

            #if DEBUG 
            {
               MX_2D( cloud_MX, q_0, t_0 ) += 2.0;
               MX_2D( cloud_MX3, dx2, k_0 ) -= 1.0;
            }
            #endif
         }
      }

      /* check that all necessary cells have been cleared */
      #if MEMCHECK
      {
         bool is_clean = false;

         for (int k_0 = 0; k_0 < (Q+1)+(T+1); k_0++) 
         {
            is_clean = false;
            is_clean += (( MMX3(dx2, k_0) == -INF ) == false);
            is_clean += (( MMX3(dx2, k_0) == -INF ) == false);
            is_clean += (( DMX3(dx2, k_0) == -INF ) == false);

            if ( is_clean != 0 ) {
               ERRORCHECK_memcheck( d_2, k_0, MMX3(dx2, k_0), IMX3(dx2, k_0), DMX3(dx2, k_0) );
               MMX3(dx2, k_0) = IMX3(dx2, k_0) = DMX3(dx2, k_0) = -INF;
               is_clean = 0;
            }
         }
      }
      #endif 

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

      /* if termination flag is set, break out of loop */
      if (is_term_flag == false) {
         break;
      }
   }

   /* scrub last two rows */
   d_last = d_0;
   for (d_0 = d_last; d_0 > d_last - 2; d_0--)
   {
      d_1 = d_0 + 1;   /* look back 1 diagonal */
      d_2 = d_0 + 2;   /* look back 2 diagonals */
      /* mod-mapping of antidiagonals into linear space */
      dx0  = d_0 % 3; 
      dx1  = d_1 % 3;
      dx2  = d_2 % 3;

      /* Scrub 2-back bound data */
      for ( i = 0; i < lb_vec[2]->N; i++ )
      {
         lb_2 = lb_vec[2]->data[i];
         rb_2 = rb_vec[2]->data[i];

         for ( k_0 = lb_2; k_0 < rb_2; k_0++ ) 
         {
            q_0 = k_0;
            t_0 = d_2 - q_0;

            MMX3(dx2, k_0) = IMX3(dx2, k_0) = DMX3(dx2, k_0) = -INF;

            #if DEBUG 
            {
               MX_2D( cloud_MX, q_0, t_0 ) += 2.0;
               MX_2D( cloud_MX3, dx2, k_0 ) -= 1.0;
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
   }

   /* check that all necessary cells have been cleared */
   #if MEMCHECK
   {
      bool is_clean = false;

      for (dx0 = 0; dx0 < 3; dx0++)
      {
         for (int k_0 = 0; k_0 < (Q+1)+(T+1); k_0++) 
         {
            is_clean = false;
            is_clean += (( MMX3(dx0, k_0) == -INF ) == false);
            is_clean += (( MMX3(dx0, k_0) == -INF ) == false);
            is_clean += (( DMX3(dx0, k_0) == -INF ) == false);
            if ( is_clean != 0 ) {
               ERRORCHECK_memcheck( dx0, k_0, MMX3(dx0, k_0), IMX3(dx0, k_0), DMX3(dx0, k_0) );
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
      VECTOR_INT_Destroy( lb_vec[i] );
      VECTOR_INT_Destroy( rb_vec[i] );
   }

   #if ( CLOUD_METHOD == CLOUD_ROWS ) 
   {
      /* output rows to edgebounds */
      EDGEBOUND_ROWS_Convert( rows, edg );

      #if DEBUG
      {
         /* compare cloud rows method to antidiagonal method */
         int cmp = EDGEBOUNDS_Compare_by_Cloud_Single( cloud_MX, edg, test_edg );
         printf("COMPARE (rows vs antidiag):\t%s\n", (cmp == 0) ? "PASS" : "FAIL");
         if ( cmp != 0 ) {
            EDGEBOUNDS_Dump( edg, stdout );
            EDGEBOUNDS_Dump( test_edg, stdout );

            printf("=== DIAG-WISE ===\n");
            MATRIX_2D_Fill( cloud_MX, 0 );
            MATRIX_2D_Cloud_Fill( cloud_MX, edg, 1 );
            DP_MATRIX_VIZ_Dump( cloud_MX, stdout );

            printf("=== ROW-WISE ===\n");
            MATRIX_2D_Fill( cloud_MX, 0 );
            MATRIX_2D_Cloud_Fill( cloud_MX, test_edg, 1 );
            DP_MATRIX_VIZ_Dump( cloud_MX, stdout );

            printf("=== OVERLAY ===\n");
            DP_MATRIX_VIZ_Compare( cloud_MX, edg, test_edg );
            DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
         }
      }
      #endif
   }
   #endif

   /* check if data is cleaned */
   #if DEBUG
   {
      DP_MATRIX_Save(Q, T, test_MX, sp_MX, "test_output/cloud_search_bck.mx");

      int cmp = MATRIX_3D_Check_Clean( st_MX3 );
      printf("POST-CHECK CLEAN -> CLOUD BCK?\t%d\n", cmp );
      printf("MAX CLOUD_BCK SCORE: %f, LIMIT: %f\n", total_max, total_limit);
   }
   #endif 

   st_MX3->clean = true;

   return STATUS_SUCCESS;
}

