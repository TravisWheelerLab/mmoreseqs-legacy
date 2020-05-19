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
 *  FUNCTION: cloud_Forward_Linear()
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
float cloud_Forward_Linear_Rows( const SEQUENCE*    query,        /* query sequence */
                                 const HMM_PROFILE* target,       /* target hmm model */
                                 const int          Q,            /* query length */
                                 const int          T,            /* target length */
                                 MATRIX_3D*         st_MX3,       /* normal state matrix */
                                 MATRIX_2D*         sp_MX,        /* special state matrix */
                                 const ALIGNMENT*   tr,           /* viterbi traceback */
                                 EDGEBOUND_ROWS*    rows,         /* temporary helper edgebounds by-row */
                                 EDGEBOUNDS*        edg,          /* (OUTPUT) */
                                 const float        alpha,        /* PARAM: pruning drop */
                                 const int          beta )        /* PARAM: free passes before pruning */
{
   /* vars for navigating matrix */
   int            b, d, i, j, k;             /* diagonal, row, column indices */
   int            le_0, re_0;                /* right/left matrix bounds of current diag */
   int            lb_0, rb_0;                /* bounds of current diag */
   int            lb_1, rb_1;                /* bounds of previous diag */
   int            lb_2, rb_2;                /* bounds of 2-back diag */
   int            num_cells;                 /* number of cells in diagonal */
   int            d_st, d_end, d_cnt;        /* starting and ending diagonal indices */
   int            dim_min, dim_max;          /* diagonal index where num cells reaches highest point and diminishing point */ 
   int            dim_T, dim_Q, dim_TOT;     /* dimensions of submatrix being searched */
   TRACE*         beg;                       /* beginning of the alignment */
   TRACE*         end;                       /* end of the alignment */

   /* vars for computing cells */
   char           a;                         /* store current character in sequence */
   int            A;                         /* store int value of character */
   char*          seq;                       /* alias for getting seq */

   /* vars for recurrance */
   int            d_0, d_1, d_2;             /* for assigning prev array ptrs */
   int            d0, d1, d2;                /* for assigning prev array ptrs (in mod3 for linear space) */
   float          prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float          sc, sc_1, sc_2, sc_best, sc_max;
   float          sc_M, sc_I, sc_D;

   /* vars for pruning */
   float          cell_max, diag_max, total_max;   /* maximum score found in matrix */
   float          total_limit, diag_limit;         /* threshold determined by max_scores - alpha */
   BOUND*         dp_bound;                        /* bounds for dp matrix of current antidiagonal */
   VECTOR_INT*    lb_vec[3];                       /* left bound list for previous 3 antdiags */
   VECTOR_INT*    rb_vec[3];                       /* right bound list for previous 3 antidiags */
   VECTOR_INT*    lb_vec_tmp;                      /* left swap pointer */
   VECTOR_INT*    rb_vec_tmp;                      /* right swap pointer */

   /* local or global? (multiple alignments) */
   bool   is_local   = target->isLocal;
   float  sc_E       = (is_local) ? 0 : -INF;

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_3D*  test_MX;
   EDGEBOUNDS* test_edg;

   /* initialize debugging matrix */
   #if DEBUG
   {
      printf("CLOUD_METHOD: %d\n", CLOUD_METHOD);

      dbfp     = fopen( debugger->dbfp_path, "w+" );
      cloud_MX = debugger->cloud_MX;
      test_MX  = debugger->test_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );

      /* second edgebound for comparisons */
      test_edg = EDGEBOUNDS_Create();
      test_edg->edg_mode = EDG_DIAG;
      EDGEBOUNDS_Reuse( test_edg, Q, T );
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   if ( st_MX3->clean == false ) {
      /* clear leftover data */
      DP_MATRIX_Fill( Q, T, st_MX3, sp_MX, -INF );
      st_MX3->clean = false;
   } 
   st_MX3->clean = true;

   #if DEBUG
   {
      int cmp =  MATRIX_3D_Check_Clean( st_MX3 );
      printf("PRE-CHECK CLEAN  -> CLOUD FWD?\t %d\n", cmp);
      if ( cmp != 0 ) {
         MATRIX_3D_Clean( st_MX3 );
      }
   }
   #endif 
   
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
   prev_beg = 0;
   prev_end = 0;

   /* ITERATE THROUGH ANTI-DIAGONALS */
   for (d = d_st; d <= d_end; d++, d_cnt++)
   {
      d_0 = d;          /* current antidiagonal */
      d_1 = (d-1);      /* look back 1 antidiagonal */
      d_2 = (d-2);      /* look back 2 antidiagonal */
      /* mod-mapping of antidiagonals into linear space */
      d0  = d_0 % 3; 
      d1  = d_1 % 3;
      d2  = d_2 % 3;

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
         prune_via_xdrop_edgetrim_Linear( st_MX3, sp_MX, alpha, beta, d_1, d_0, d1, d0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_XDROP_SPLIT )
      {
         /* prune bounds using x-drop, bifurcating */
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
            prev_mat = MMX3(d2,k-1)  + TSC(j-1,M2M);
            prev_ins = IMX3(d2,k-1)  + TSC(j-1,I2M);
            prev_del = DMX3(d2,k-1)  + TSC(j-1,D2M);
            /* free to begin match state (new alignment) */
            // prev_beg = 0; /* assigned once at start */
            /* best-to-match */
            prev_sum = logsum( 
                           logsum( prev_mat, prev_ins ),
                           logsum( prev_del, prev_beg ) );
            MMX3(d0,k) = prev_sum + MSC(j,A);

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            /* NOTE: Convert (i-1,j) <=> (d-1,k-1) */
            prev_mat = MMX3(d1,k-1) + TSC(j,M2I);
            prev_ins = IMX3(d1,k-1) + TSC(j,I2I);
            /* best-to-insert */
            prev_sum = logsum( prev_mat, prev_ins );
            IMX3(d0,k) = prev_sum + ISC(j,A);

            /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            /* NOTE: Convert (i,j-1) <=> (d-1, k) */
            prev_mat = MMX3(d1,k) + TSC(j-1,M2D);
            prev_del = DMX3(d1,k) + TSC(j-1,D2D);
            /* best-to-delete */
            prev_sum = logsum(prev_mat, prev_del);
            DMX3(d0,k) = prev_sum;

            /* embed cell data in quadratic matrix */
            #if DEBUG
            {
               MX_2D( cloud_MX, i, j ) += DIRTY_VAL;
               MX_3D( test_MX, MAT_ST, i, j ) = MMX3(d0, k);
               MX_3D( test_MX, INS_ST, i, j ) = IMX3(d0, k);
               MX_3D( test_MX, DEL_ST, i, j ) = DMX3(d0, k);
            }
            #endif 
         }
      }

      // /* Naive Scrub */
      // for (k = 0; k < (T+1)+(Q+1); k++) {
      //    MMX3(d2, k) = IMX3(d2, k) = DMX3(d2, k) = -INF;
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
            MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = -INF;
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
      prev_beg = -INF;
      // prev_end = -INF;
   }

   /* TODO: Scrub last two rows */
   for (d = (d_end + 1); d < (d_end + 3); d++)
   {
      d_0 = d;          /* current antidiagonal */
      d_1 = (d-1);      /* look back 1 antidiagonal */
      d_2 = (d-2);      /* look back 2 antidiagonal */
      /* mod-mapping of antidiagonals into linear space */
      d0  = d_0 % 3; 
      d1  = d_1 % 3;
      d2  = d_2 % 3;

      /* Scrub 2-back bound data */
      for ( b = 0; b < lb_vec[2]->N; b++ )
      {
         lb_2 = lb_vec[2]->data[b];
         rb_2 = rb_vec[2]->data[b];

         for ( k = lb_2; k < rb_2; k++ ) 
         {
            i = k;
            j = d_2 - i;
            MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = -INF;
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
         DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
      }
   }
   #endif
   /* show visualization of search cloud */
   #if DEBUG
   {
      DP_MATRIX_VIZ_Trace( cloud_MX, tr );
      DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
      // DP_MATRIX_Trace_Dump( Q, T, test_MX, sp_MX, tr, stdout );

      /* final test that all cells are cleared */
      int cmp = MATRIX_3D_Check_Clean( st_MX3 );
      printf("POST-CHECK CLEAN -> CLOUD FWD?\t %d\n", cmp );
      if ( cmp != 0 )  {
         // MATRIX_3D_Dump( st_MX3, stdout );
         MATRIX_3D_Clean( st_MX3 );
      }

      /* close debugger tools */
      fclose( dbfp );
   }
   #endif

   return total_max;
}


/*
 *  FUNCTION: cloud_Backward_Linear()
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
float cloud_Backward_Linear_Rows(   const SEQUENCE*   query,         /* query sequence */
                                    const HMM_PROFILE* target,       /* target hmm model */
                                    const int          Q,            /* query length */
                                    const int          T,            /* target length */
                                    MATRIX_3D*         st_MX3,       /* normal state matrix */
                                    MATRIX_2D*         sp_MX,        /* special state matrix */
                                    const ALIGNMENT*   tr,           /* viterbi traceback */
                                    EDGEBOUND_ROWS*    rows,         /* temporary edgebounds by-row */
                                    EDGEBOUNDS*        edg,          /* (OUTPUT) edgebounds */
                                    const float        alpha,        /* PARAM: pruning drop */
                                    const int          beta )        /* PARAM: free passes before pruning */
{
   /* vars for navigating matrix */
   int            b, d, i, j, k;             /* diagonal, row, column indices */
   int            le_0, re_0;                /* right/left matrix bounds of current diag */
   int            lb_0, rb_0;                /* bounds of current diag */
   int            lb_1, rb_1;                /* bounds of previous diag */
   int            lb_2, rb_2;                /* bounds of 2-back diag */
   int            num_cells;                 /* number of cells in diagonal */
   int            d_st, d_end, d_cnt;        /* starting and ending diagonal indices */
   int            dim_min, dim_max;          /* diagonal index where num cells reaches highest point and diminishing point */ 
   int            dim_T, dim_Q, dim_TOT;     /* dimensions of submatrix being searched */
   TRACE*         beg;                       /* beginning of the alignment */
   TRACE*         end;                       /* end of the alignment */

   /* vars for computing cells */
   char           a;                         /* store current character in sequence */
   int            A;                         /* store int value of character */
   char*          seq;                       /* alias for getting seq */

   /* vars for recurrance */
   int            d_0, d_1, d_2;             /* for assigning prev array ptrs */
   int            d0, d1, d2;                /* for assigning prev array ptrs (in mod3 for linear space) */
   float          prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float          sc, sc_1, sc_2, sc_best, sc_max;
   float          sc_M, sc_I, sc_D;

   /* vars for pruning */
   float          cell_max, diag_max, total_max;   /* maximum score found in matrix */
   float          total_limit, diag_limit;         /* threshold determined by max_scores - alpha */
   BOUND*         dp_bound;                        /* bounds for dp matrix of current antidiagonal */
   VECTOR_INT*    lb_vec[3];                       /* left bound list for previous 3 antdiags */
   VECTOR_INT*    rb_vec[3];                       /* right bound list for previous 3 antidiags */
   VECTOR_INT*    lb_vec_tmp;                      /* left swap pointer */
   VECTOR_INT*    rb_vec_tmp;                      /* right swap pointer */

   /* local or global? (multiple alignments) */
   bool   is_local   = target->isLocal;
   float  sc_E       = (is_local) ? 0 : -INF;

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_3D*  test_MX;
   EDGEBOUNDS* test_edg;

   /* initialize debugging matrix */
   #if DEBUG
   {
      printf("CLOUD_METHOD: %d\n", CLOUD_METHOD);

      dbfp     = fopen( debugger->dbfp_path, "w+" );
      cloud_MX = debugger->cloud_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      test_MX  = debugger->test_MX;
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );

      test_edg = EDGEBOUNDS_Create();
      test_edg->edg_mode = EDG_DIAG;
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* clear leftover data if necessary */
   if ( st_MX3->clean == false ) {
      DP_MATRIX_Fill( Q, T, st_MX3, sp_MX, -INF );
      st_MX3->clean = true;
   }
   st_MX3->clean = true;

   /* verify data is clean */
   #if DEBUG
   {
      int cmp =  MATRIX_3D_Check_Clean( st_MX3 );
      printf("PRE-CHECK CLEAN  -> CLOUD BCK?\t %d\n", cmp);
      if ( cmp != 0 ) {
         MATRIX_3D_Clean( st_MX3 );
      }
   }
   #endif 
   
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

   /* set edgebound dimensions and orientation */
   edg->Q         = Q;
   edg->T         = T;
   edg->edg_mode  = EDG_DIAG;

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
   prev_beg = 0;
   prev_end = 0; 

   /* ITERATE THROUGHT ANTI-DIAGONALS */
   for (d = d_end; d >= d_st; d--, d_cnt++)
   {
      d_0 = d;       /* current diagonal */
      d_1 = (d+1);   /* look back 1 diagonal */
      d_2 = (d+2);   /* look back 2 diagonals */
      /* mod-mapping of antidiagonals into linear space */
      d0  = d_0 % 3; 
      d1  = d_1 % 3;
      d2  = d_2 % 3;

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
         prune_via_xdrop_edgetrim_Linear( st_MX3, sp_MX, alpha, beta, d_1, d_0, d1, d0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );
      }
      #elif ( PRUNER == PRUNER_XDROP_SPLIT )
      {
         /* prune bounds using x-drop, bifurcating */
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
            prev_mat = MMX3(d2, k+1) + TSC(j, M2M) + sc_M;
            prev_ins = IMX3(d1, k+1) + TSC(j, M2I) + sc_I;
            prev_del = DMX3(d1, k  ) + TSC(j, M2D);
            // prev_end = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
            // prev_end = sc_E;
            /* best-to-match */
            prev_sum = logsum( 
                           logsum( prev_mat, prev_ins ),
                           logsum( prev_del, prev_end ) );
            MMX3(d0,k) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            sc_I = ISC(j, A);

            prev_mat = MMX3(d2, k+1) + TSC(j, I2M) + sc_M;
            prev_ins = IMX3(d1, k+1) + TSC(j, I2I) + sc_I;
            /* best-to-insert */
            prev_sum = logsum( prev_mat, prev_ins );
            IMX3(d0,k) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            prev_mat = MMX3(d2, k+1) + TSC(j, D2M) + sc_M;
            prev_del = DMX3(d1, k  ) + TSC(j, D2D);
            /* best-to-delete */
            prev_sum = logsum( prev_mat, prev_del );
            prev_sum = logsum( prev_sum, prev_end );
            DMX3(d0,k) = prev_sum;

            /* embed cell data in quadratic matrix */
            #if DEBUG
            {
               MX_2D( cloud_MX, i, j ) += DIRTY_VAL;
               MX_3D( test_MX, MAT_ST, i, j ) = MMX3(d0, k);
               MX_3D( test_MX, INS_ST, i, j ) = IMX3(d0, k);
               MX_3D( test_MX, DEL_ST, i, j ) = DMX3(d0, k);
            }
            #endif 
         }  
      }

      // /* Naive Scrub */
      // for (k = 0; k < (T+1)+(Q+1); k++) {
      //    MMX3(d2, k) = IMX3(d2, k) = DMX3(d2, k) = -INF;
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
            MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = -INF;
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
      // prev_beg = -INF;
      prev_end = -INF;
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
      d0  = d_0 % 3; 
      d1  = d_1 % 3;
      d2  = d_2 % 3;

      /* Scrub 2-back bound data */
      for ( b = 0; b < lb_vec[2]->N; b++ )
      {
         lb_2 = lb_vec[2]->data[b];
         rb_2 = rb_vec[2]->data[b];

         for ( k = lb_2; k < rb_2; k++ ) 
         {
            i = k;
            j = d_2 - i;
            MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = -INF;
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
         DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
      }
   }
   #endif
   #if DEBUG
   {
      /* show visualization of search cloud */
      DP_MATRIX_VIZ_Trace( cloud_MX, tr );
      DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
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
   #if DEBUG
   {
      /* close necessary debugger tools */
      fclose(dbfp);
   }
   #endif

   return total_max;
}

