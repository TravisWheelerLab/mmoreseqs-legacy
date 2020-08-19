/*******************************************************************************
 *  FILE:      cloud_search_linear.c
 *  PURPOSE:   Cloud Search for Forward-Backward Pruning Alg.
 *             (Linear Space Alg)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
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
#include "algs_linear.h"
#include "algs_quad.h"

/* header */
#include "cloud_search_linear.h"

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
 *  RETURN:   Maximum score.
 */
int  run_Cloud_Forward_Linear_Rows(    const SEQUENCE*    query,        /* query sequence */
                                       const HMM_PROFILE* target,       /* target hmm model */
                                       const int          Q,            /* query length */
                                       const int          T,            /* target length */
                                       MATRIX_3D*         st_MX3,       /* normal state matrix */
                                       MATRIX_2D*         sp_MX,        /* special state matrix */
                                       const ALIGNMENT*   tr,           /* viterbi traceback */
                                       EDGEBOUND_ROWS*    rows,         /* temporary edgebounds by-row */
                                       EDGEBOUNDS*        edg,          /* (OUTPUT) */
                                       CLOUD_PARAMS*      params )      /* pruning parameters */
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
   int      d_st, d_end, d_cnt, d_last;      /* starting and ending diagonal indices */
   int      dim_T, dim_Q, dim_TOT;           /* dimensions of submatrix being searched */
   int      dim_min, dim_max;                /* diagonal index where num cells reaches highest point and diminishing point */ 
   int      num_cells;                       /* number of cells in current diagonal */

   /* vars for indexing into edgebound lists */
   int      x, y1, y2;                       /* row, leftcol and rightcol bounds in row (edgebounds) */
   int      e_0b, e_0e;                      /* begin and end indices for current row in edgebound list */
   int      e_1b, e_1e;                      /* begin and end indices for current row in edgebound list */
   int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   bool     y2_re;                           /* checks if edge touches right bound of matrix */

   /* vars for recurrance scores */
   float    prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;              /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C; /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;            /* previous loop and move for special states */
   float    prv_sum, prv_best;             /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* vars for traceback */
   TRACE*         beg;                       /* beginning of the alignment */
   TRACE*         end;                       /* end of the alignment */

   /* vars for pruning */
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
      EDGEBOUNDS_Reuse( test_edg, Q ,T );

      num_writes = 0;
      num_clears = 0;
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

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
   
   /* set edgebound dimensions and orientation */
   edg->Q         = Q;
   edg->T         = T;
   #if ( CLOUD_METHOD == CLOUD_DIAGS )
   {
      EDGEBOUNDS_Reuse( edg, Q, T );
      edg->edg_mode  = EDG_DIAG;
   }
   #elif ( CLOUD_METHOD == CLOUD_ROWS )
   {
      /* clear rows data */
      EDGEBOUND_ROWS_Reuse( rows, Q, T );
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
   beg = &(tr->traces[tr->beg]);
   end = &(tr->traces[tr->end]);

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if ( beg->i == 0 || beg->j == 0 ) {
      beg->i += 1;
      beg->j += 1;
   }

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

   /* TODO: set initial match to zero ( single free entry point ) */
   d = d_st-1;

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
   prv_E = 0;

   /* ITERATE THROUGH ANTI-DIAGONALS */
   for (d = d_st; d <= d_end; d++, d_cnt++)
   {
      d_0 = d;          /* current antidiagonal */
      d_1 = (d-1);      /* look back 1 antidiagonal */
      d_2 = (d-2);      /* look back 2 antidiagonal */
      /* mod-mapping of antidiagonals into linear space */
      dx0  = d_0 % 3; 
      dx1  = d_1 % 3;
      dx2  = d_2 % 3;

      /* is dp matrix diagonal growing or shrinking? */
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
         prune_via_xdrop_edgetrim_Linear( st_MX3, sp_MX, alpha, beta, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_XDROP_BIFURCATE )
      {
         /* prune bounds using x-drop, bifurcating */
         prune_via_xdrop_bifurcate_Linear( st_MX3, sp_MX, alpha, beta, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_DBL_XDROP_EDGETRIM_OR_DIE )
      {
         /* prune bounds using local and global x-drop, edgetrimming or terminating search */
         prune_via_dbl_xdrop_edgetrim_or_die_Linear( st_MX3, sp_MX, alpha, beta, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #endif

      /* Add pruned bounds to edgebound list */
      for ( b = 0; b < lb_vec[0]->N; b++ )
      {
         lb_0 = lb_vec[0]->data[b];
         rb_0 = rb_vec[0]->data[b];

         /* Update bounds (spans all cells adjacent to previous antidiagonals cells that were not pruned) */
         lb_0 = lb_0;
         rb_0 = rb_0 + 1;

         /* test: no pruning */
         #if DEBUG && PRUNER_NONE
         {
            lb_0 = lb_vec[1]->data[b];
            rb_0 = rb_vec[1]->data[b];
            lb_0 = lb_0;
            rb_0 = rb_0 + 1;
         }
         #endif

         /* Update bounds to account for dp matrix bounds */
         lb_0 = MAX(lb_0, le_0);
         rb_0 = MIN(rb_0, re_0);

         /* Update changes to list */
         lb_vec[0]->data[b] = lb_0;
         rb_vec[0]->data[b] = rb_0;

         BOUND bnd = (BOUND){d_0,lb_0,rb_0};

         #if ( CLOUD_METHOD == CLOUD_DIAGS )
         {
            /* Add new bounds to edgebounds (save for testing) */
            EDGEBOUNDS_Pushback( edg, &bnd );
         }
         #elif ( CLOUD_METHOD == CLOUD_ROWS )
         {
            /* Reorient new diagonal bounds and integrate it into row-wise edgebound list */
            EDGEBOUND_ROWS_Integrate_Antidiag_Fwd( rows, &bnd );
         }
         #endif
         #if DEBUG && ( CLOUD_METHOD == CLOUD_ROWS )
         {
            /* Add new bounds to edgebounds (save for testing) */
            EDGEBOUNDS_Pushback( test_edg, &bnd );
         }
         #endif
      }

      /* If diagonal set is empty, then all branches have been pruned, so we're done */
      // printf("lb_vec_length: %d\n", lb_vec[0]->N );
      if ( lb_vec[0]->N <= 0 ) break;

      /* MAIN RECURSION */
      /* Iterate the bound ranges of the current antidiagonal */
      for ( b = 0; b < lb_vec[0]->N; b++ ) 
      {
         lb_0 = lb_vec[0]->data[b];
         rb_0 = rb_vec[0]->data[b];

         /* Iterate through cells in range */
         for ( k = lb_0; k < rb_0; k++ )
         {
            /* quadratic coords */
            i = k;
            j = d_0 - i;

            a = seq[i];
            A = AA_REV[a];

            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            /* NOTE: Convert (i-1,j-1) <=> (d-2,k-1) */ 
            prv_M = MMX3(dx2,k-1)  + TSC(j-1,M2M);
            prv_I = IMX3(dx2,k-1)  + TSC(j-1,I2M);
            prv_D = DMX3(dx2,k-1)  + TSC(j-1,D2M);
            /* free to begin match state (new alignment) */
            // prv_B = 0; /* assigned once at start */
            /* best-to-match */
            prv_sum = logsum( 
                           logsum( prv_M, prv_I ),
                           logsum( prv_D, prv_B ) );
            MMX3(dx0,k) = prv_sum + MSC(j,A);

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            /* NOTE: Convert (i-1,j) <=> (d-1,k-1) */
            prv_M = MMX3(dx1,k-1) + TSC(j,M2I);
            prv_I = IMX3(dx1,k-1) + TSC(j,I2I);
            /* best-to-insert */
            prv_sum = logsum( prv_M, prv_I );
            IMX3(dx0,k) = prv_sum + ISC(j,A);

            /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            /* NOTE: Convert (i,j-1) <=> (d-1, k) */
            prv_M = MMX3(dx1,k) + TSC(j-1,M2D);
            prv_D = DMX3(dx1,k) + TSC(j-1,D2D);
            /* best-to-delete */
            prv_sum = logsum(prv_M, prv_D);
            DMX3(dx0,k) = prv_sum;

            /* embed cell data in quadratic matrix */
            #if DEBUG
            {
               MX_2D( cloud_MX, i, j ) += DIRTY_VAL;
               MX_3D( test_MX, MAT_ST, i, j ) = MMX3(dx0, k);
               MX_3D( test_MX, INS_ST, i, j ) = IMX3(dx0, k);
               MX_3D( test_MX, DEL_ST, i, j ) = DMX3(dx0, k);
            }
            #endif 
         }
      }

      // /* Naive Scrub */
      // for (k = 0; k < (T+1)+(Q+1); k++) {
      //    MMX3(dx2, k) = IMX3(dx2, k) = DMX3(dx2, k) = -INF;
      // }

      /* Scrub 2-back bound data */
      for ( b = 0; b < lb_vec[2]->N; b++ )
      {
         lb_2 = lb_vec[2]->data[b];
         rb_2 = rb_vec[2]->data[b];

         for ( k = lb_2; k < rb_2; k++ ) 
         {
            i = k;
            j = d_2 - i;
            MMX3(dx2,k) = IMX3(dx2,k) = DMX3(dx2,k) = -INF;
            #if DEBUG
            {
               MX_2D( cloud_MX, i, j ) += SCRUB_VAL;
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
      prv_B = -INF;
      // prv_E = -INF;
   }

   /* TODO: Scrub last two rows */
   for (d = (d_end + 1); d < (d_end + 3); d++)
   {
      d_0 = d;          /* current antidiagonal */
      d_1 = (d-1);      /* look back 1 antidiagonal */
      d_2 = (d-2);      /* look back 2 antidiagonal */
      /* mod-mapping of antidiagonals into linear space */
      dx0  = d_0 % 3; 
      dx1  = d_1 % 3;
      dx2  = d_2 % 3;

      /* Scrub 2-back bound data */
      for ( b = 0; b < lb_vec[2]->N; b++ )
      {
         lb_2 = lb_vec[2]->data[b];
         rb_2 = rb_vec[2]->data[b];

         for ( k = lb_2; k < rb_2; k++ ) 
         {
            i = k;
            j = d_2 - i;
            MMX3(dx2,k) = IMX3(dx2,k) = DMX3(dx2,k) = -INF;
            #if DEBUG
            {
               MX_2D( cloud_MX, i, j ) += SCRUB_VAL;
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

   /* free dynamic memory */
   for (i = 0; i < 3; i++) {
      VECTOR_INT_Destroy( lb_vec[i] );
      VECTOR_INT_Destroy( rb_vec[i] );
   }

   #if ( CLOUD_METHOD == CLOUD_ROWS ) 
   {
      /* output rows to edgebounds */
      EDGEBOUND_ROWS_Convert( rows, edg );
   }
   #endif
   #if DEBUG && ( CLOUD_METHOD == CLOUD_ROWS )
   {
      /* compare cloud rows method to antidiagonal method */
      int cmp = EDGEBOUNDS_Compare_by_Cloud_Single( cloud_MX, edg, test_edg );
      printf("COMPARE (rows vs antidiag):\t%s\n", (cmp == 0) ? "PASS" : "FAIL");
      if ( cmp != 0 ) {
         DP_MATRIX_VIZ_Compare( cloud_MX, edg, test_edg );
         // DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
      }
   }
   #endif
   /* show visualization of search cloud */
   #if DEBUG
   {
      // DP_MATRIX_VIZ_Trace( cloud_MX, tr );
      // DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
      // DP_MATRIX_Trace_Dump( Q, T, test_MX, sp_MX, tr, stdout );

      /* final test that all cells are cleared */
      int cmp = MATRIX_3D_Check_Clean( st_MX3 );
      printf("POST-CHECK CLEAN -> CLOUD FWD?\t %d\n", cmp );
      if ( cmp != 0 )  {
         // MATRIX_3D_Dump( st_MX3, stdout );
         MATRIX_3D_Clean( st_MX3 );
      }

      /* close debugger tools */
      printf("# EDGEBOUNDS:\n");
      EDGEBOUNDS_Dump( edg, stdout );
      test_edg = EDGEBOUNDS_Destroy( test_edg );
      fclose( dbfp );
   }
   #endif

   return total_max;
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
int run_Cloud_Backward_Linear_Rows(     const SEQUENCE*   query,         /* query sequence */
                                    const HMM_PROFILE* target,       /* target hmm model */
                                    const int          Q,            /* query length */
                                    const int          T,            /* target length */
                                    MATRIX_3D*         st_MX3,       /* normal state matrix */
                                    MATRIX_2D*         sp_MX,        /* special state matrix */
                                    const ALIGNMENT*   tr,           /* viterbi traceback */
                                    EDGEBOUND_ROWS*    rows,         /* temporary edgebounds by-row */
                                    EDGEBOUNDS*        edg,          /* (OUTPUT) edgebounds */
                                    CLOUD_PARAMS*      params )      /* pruning parameters */
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
   int      d_st, d_end, d_cnt, d_last;      /* starting and ending diagonal indices */
   int      dim_T, dim_Q, dim_TOT;           /* dimensions of submatrix being searched */
   int      dim_min, dim_max;                /* diagonal index where num cells reaches highest point and diminishing point */ 
   int      num_cells;                       /* number of cells in current diagonal */

   /* vars for indexing into edgebound lists */
   int      x, y1, y2;                       /* row, leftcol and rightcol bounds in row (edgebounds) */
   int      e_0b, e_0e;                      /* begin and end indices for current row in edgebound list */
   int      e_1b, e_1e;                      /* begin and end indices for current row in edgebound list */
   int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   bool     y2_re;                           /* checks if edge touches right bound of matrix */

   /* vars for recurrance scores */
   float    prv_M, prv_I, prv_D;    /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;              /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C; /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;            /* previous loop and move for special states */
   float    prv_sum, prv_best;             /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* vars for traceback */
   TRACE*         beg;                       /* beginning of the alignment */
   TRACE*         end;                       /* end of the alignment */

   /* vars for pruning */
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
   
   /* set edgebound dimensions and orientation */
   edg->Q         = Q;
   edg->T         = T;

   #if ( CLOUD_METHOD == CLOUD_DIAGS )
   {
      EDGEBOUNDS_Reuse( edg, Q, T );
      edg->edg_mode  = EDG_DIAG;
   }
   #elif ( CLOUD_METHOD == CLOUD_ROWS )
   {
      /* clear rows data */
      EDGEBOUND_ROWS_Reuse( rows, Q, T );
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
   beg = &(tr->traces[tr->beg]);
   end = &(tr->traces[tr->end]);

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (end->i == Q || end->j == T) {
      end->i -= 1;
      end->j -= 1;
   }

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

   /* TODO: set initial match to zero ( single free entry point ) */
   d = d_end+1;

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

   /* begin state probability begins at zero (free to start first alignment position) */
   prv_B = 0;
   prv_E = 0; 

   /* ITERATE THROUGHT ANTI-DIAGONALS */
   for (d = d_end; d >= d_st; d--, d_cnt++)
   {
      d_0 = d;       /* current diagonal */
      d_1 = (d+1);   /* look back 1 diagonal */
      d_2 = (d+2);   /* look back 2 diagonals */
      /* mod-mapping of antidiagonals into linear space */
      dx0  = d_0 % 3; 
      dx1  = d_1 % 3;
      dx2  = d_2 % 3;

      /* Is dp matrix diagonal growing or shrinking? */
      if (d >= dim_max)
         num_cells++;
      if (d < dim_min)
         num_cells--;

      /* Edgecheck updates: determine antidiag indices within matrix bounds */
      le_0 = MAX( end->i - (d_end - d_0), 0 );
      re_0 = le_0 + num_cells;

      /* Prune bounds */
      #if ( PRUNER == PRUNER_XDROP_EDGETRIM )
      {
         /* prune bounds using x-drop, no bifurcating */
         prune_via_xdrop_edgetrim_Linear( st_MX3, sp_MX, alpha, beta, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_XDROP_BIFURCATE )
      {
         /* prune bounds using x-drop, bifurcating */
         prune_via_xdrop_bifurcate_Linear( st_MX3, sp_MX, alpha, beta, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_DBL_XDROP_EDGETRIM_OR_DIE )
      {
         /* prune bounds using local and global x-drop, edgetrimming or terminating search */
         prune_via_dbl_xdrop_edgetrim_or_die_Linear( st_MX3, sp_MX, alpha, beta, gamma, d_1, d_0, dx1, dx0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #endif

      /* Add pruned bounds to edgebound list */
      for ( b = 0; b < lb_vec[0]->N; b++ )
      {
         lb_0 = lb_vec[0]->data[b];
         rb_0 = rb_vec[0]->data[b];

         /* Update bounds (spans all cells adjacent to previous antidiagonals cells that were not pruned) */
         lb_0 = lb_0 - 1;
         rb_0 = rb_0;

         /* test: no pruning */
         #if DEBUG
         {
            // lb_0 = lb_vec[1]->data[b];
            // rb_0 = rb_vec[1]->data[b];
            // lb_0 = lb_0 - 1;
            // rb_0 = rb_0;
         }
         #endif

         /* Update bounds to account for dp matrix bounds */
         lb_0 = MAX(lb_0, le_0);
         rb_0 = MIN(rb_0, re_0);

         /* Update changes to list */
         lb_vec[0]->data[b] = lb_0;
         rb_vec[0]->data[b] = rb_0;

         BOUND bnd = (BOUND){d_0,lb_0,rb_0};

         #if ( CLOUD_METHOD == CLOUD_DIAGS )
         {
            /* add new bounds to edgebounds as antidiag-wise */
            EDGEBOUNDS_Pushback( edg, &bnd );
         }
         #endif
         #if ( CLOUD_METHOD == CLOUD_ROWS )
         {
            /* reorient new bounds from antidiag-wise to row-wise and integrate it into row-wise edgebound list */
            EDGEBOUND_ROWS_Integrate_Antidiag_Bck( rows, &bnd );
         }
         #endif
         #if DEBUG && ( CLOUD_METHOD == CLOUD_ROWS )
         {
            /* add new bounds to edgebounds as antidiag-wise (for comparative testing) */
            EDGEBOUNDS_Pushback( test_edg, &bnd );
         }
         #endif
      }

      /* MAIN RECURSION */
      /* Iterate through the range bounds of the current antidiagonal */
      for ( b = 0; b < lb_vec[0]->N; b++ )
      {
         lb_0 = lb_vec[0]->data[b];
         rb_0 = rb_vec[0]->data[b];

         /* Iterate through cells in current bound */
         for ( k = lb_0; k < rb_0; k++ )
         {
            /* get x-y coords */
            i = k;
            j = d_0 - i;

            /*    
             *    === ROW-WISE to DIAG_WISE ===
             *    MX_M(i+1, j+1) => MX3_M(d_2, k+1)
             *    MX_M(i  , j+1) => MX3_M(d_1, k  )
             *    MX_M(i+1, j  ) => MX3_M(d_1, k+1)
             */

            /* next sequence character */
            a = seq[i];
            A = AA_REV[a];

            /* match and insertion scores */
            sc_M = MSC(j+1, A);
            sc_I = ISC(j+1, A);

            /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
            prv_M = MMX3(dx2, k+1) + TSC(j, M2M) + sc_M;
            prv_I = IMX3(dx1, k+1) + TSC(j, M2I) + sc_I;
            prv_D = DMX3(dx1, k  ) + TSC(j, M2D);
            // prv_E = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
            // prv_E = sc_E;
            /* best-to-match */
            prv_sum = logsum( 
                           logsum( prv_M, prv_I ),
                           logsum( prv_D, prv_E ) );
            MMX3(dx0,k) = prv_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            sc_I = ISC(j, A);

            prv_M = MMX3(dx2, k+1) + TSC(j, I2M) + sc_M;
            prv_I = IMX3(dx1, k+1) + TSC(j, I2I) + sc_I;
            /* best-to-insert */
            prv_sum = logsum( prv_M, prv_I );
            IMX3(dx0,k) = prv_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            prv_M = MMX3(dx2, k+1) + TSC(j, D2M) + sc_M;
            prv_D = DMX3(dx1, k  ) + TSC(j, D2D);
            /* best-to-delete */
            prv_sum = logsum( prv_M, prv_D );
            prv_sum = logsum( prv_sum, prv_E );
            DMX3(dx0,k) = prv_sum;

            /* embed cell data in quadratic matrix */
            #if DEBUG
            {
               MX_2D( cloud_MX, i, j ) += DIRTY_VAL;
               MX_3D( test_MX, MAT_ST, i, j ) = MMX3(dx0, k);
               MX_3D( test_MX, INS_ST, i, j ) = IMX3(dx0, k);
               MX_3D( test_MX, DEL_ST, i, j ) = DMX3(dx0, k);
            }
            #endif 
         }  
      }

      // /* Naive Scrub */
      // for (k = 0; k < (T+1)+(Q+1); k++) {
      //    MMX3(dx2, k) = IMX3(dx2, k) = DMX3(dx2, k) = -INF;
      // }

      /* Scrub 2-back bound data */
      for ( b = 0; b < lb_vec[2]->N; b++ )
      {
         lb_2 = lb_vec[2]->data[b];
         rb_2 = rb_vec[2]->data[b];

         for ( k = lb_2; k < rb_2; k++ ) 
         {
            i = k;
            j = d_2 - i;
            MMX3(dx2,k) = IMX3(dx2,k) = DMX3(dx2,k) = -INF;
            #if DEBUG
            {
               MX_2D( cloud_MX, i, j ) += SCRUB_VAL;
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

   /* TODO: scrub last two rows */
   for (d = (d_st - 1); d > (d_st - 3); d--)
   {
      d_0 = d;       /* current diagonal */
      d_1 = (d+1);   /* look back 1 diagonal */
      d_2 = (d+2);   /* look back 2 diagonals */
      /* mod-mapping of antidiagonals into linear space */
      dx0  = d_0 % 3; 
      dx1  = d_1 % 3;
      dx2  = d_2 % 3;

      /* Scrub 2-back bound data */
      for ( b = 0; b < lb_vec[2]->N; b++ )
      {
         lb_2 = lb_vec[2]->data[b];
         rb_2 = rb_vec[2]->data[b];

         for ( k = lb_2; k < rb_2; k++ ) 
         {
            i = k;
            j = d_2 - i;
            MMX3(dx2,k) = IMX3(dx2,k) = DMX3(dx2,k) = -INF;
            #if DEBUG
            {
               MX_2D( cloud_MX, i, j ) += SCRUB_VAL;
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

   /* free dynamic memory */
   for (i = 0; i < 3; i++) {
      VECTOR_INT_Destroy( lb_vec[i] );
      VECTOR_INT_Destroy( rb_vec[i] );
   }

   #if ( CLOUD_METHOD == CLOUD_ROWS ) 
   {
      /* output rows to edgebounds */
      EDGEBOUND_ROWS_Convert( rows, edg );
   }
   #endif
   #if DEBUG && ( CLOUD_METHOD == CLOUD_ROWS )
   {
      /* compare cloud rows method to antidiagonal method */
      int cmp = EDGEBOUNDS_Compare_by_Cloud_Single( cloud_MX, edg, test_edg );
      printf("COMPARE (rows vs antidiag):\t%s\n", (cmp == 0) ? "PASS" : "FAIL");
      if ( cmp != 0 ) {
         DP_MATRIX_VIZ_Compare( cloud_MX, edg, test_edg );
         DP_MATRIX_VIZ_Dump( cloud_MX, dbfp );
      }
   }
   #endif

   #if DEBUG
   {
      /* show visualization of search cloud */
      DP_MATRIX_VIZ_Trace( cloud_MX, tr );
      DP_MATRIX_VIZ_Dump( cloud_MX, dbfp );
      // DP_MATRIX_Trace_Dump( Q, T, test_MX, sp_MX, tr, stdout );

      /* final test that all cells are cleared */
      int cmp = MATRIX_3D_Check_Clean( st_MX3 );
      printf("POST-CHECK CLEAN -> CLOUD BCK?\t %d\n", cmp );
      if ( cmp != 0 ) {
         // MATRIX_3D_Dump( st_MX3, stdout );
         MATRIX_3D_Clean( st_MX3 );
      }
   }
   #endif 

   return total_max;
}
