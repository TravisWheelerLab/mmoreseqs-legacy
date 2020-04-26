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
 *  FUNCTION:  cloud_Forward_Run()
 *  SYNOPSIS:  Perform Forward part of Cloud Search Algorithm.
 */
float cloud_Forward_Quad(  const SEQUENCE*      query,         /* query sequence */
                           const HMM_PROFILE*   target,        /* target model */
                           const int            Q,             /* query length */
                           const int            T,             /* target length */
                           MATRIX_3D*           st_MX,         /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
                           MATRIX_2D*           sp_MX,         /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
                           const ALIGNMENT*     tr,            /* viterbi traceback */ 
                           EDGEBOUNDS*          edg,           /* OUTPUT: edgebounds of cloud search space */
                           float                alpha,         /* PARAM: x-drop threshold */
                           int                  beta )         /* PARAM: free passes */
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

   /* debugging matrix */
   #if DEBUG
   {
      dbfp     = fopen( debugger->dbfp_path, "w+" );
      cloud_MX = debugger->cloud_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      test_MX  = debugger->test_MX;
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* clear leftover data */
   DP_MATRIX_Fill( Q, T, st_MX, sp_MX, -INF );

   /* set edgebound dimensions and orientation */
   edg->Q         = Q;
   edg->T         = T;
   edg->edg_mode  = EDG_DIAG;

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

   /* ITERATE THROUGH ANTI-DIAGONALS */
   for ( d = d_st; d <= d_end; d++, d_cnt++ )
   {
      d_0 = d;          /* current antidiagonal */
      d_1 = (d-1);      /* look back 1 antidiagonal */
      d_2 = (d-2);      /* look back 2 antidiagonal */
      /* mod-mapping of antidiagonals */
      d0  = d_0; 
      d1  = d_1;
      d2  = d_2;

      /* Is dp matrix diagonal growing or shrinking? */
      if ( d_0 <= dim_min )
         num_cells++;
      if ( d_0 > dim_max )
         num_cells--;

      /* Edgecheck updates: determine antidiag indices within matrix bounds */
      le_0 = MAX( beg->i, d_0 - T );
      re_0 = le_0 + num_cells;

      /* Prune bounds */
      prune_via_xdrop_edgetrim_Quad( st_MX, sp_MX, alpha, beta, d_1, d_0, d1, d0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );

      /* Add pruned bounds to edgebound list */
      for ( b = 0; b < lb_vec[0]->N; b++ )
      {
         lb_0 = lb_vec[0]->data[b];
         rb_0 = rb_vec[0]->data[b];

         /* Update bounds (spans all cells adjacent to previous antidiagonals cells that were not pruned) */
         lb_0 = lb_0;
         rb_0 = rb_0 + 1;

         /* test: no pruning */
         #if DEBUG
         {
            // lb_0 = lb_vec[1]->data[b];
            // rb_0 = rb_vec[1]->data[b];
            // lb_0 = lb_0;
            // rb_0 = rb_0 + 1;
         }
         #endif

         /* Update bounds to account for dp matrix bounds */
         lb_0 = MAX(lb_0, le_0);
         rb_0 = MIN(rb_0, re_0);

         /* Update changes to list */
         lb_vec[0]->data[b] = lb_0;
         rb_vec[0]->data[b] = rb_0;

         /* Add new bounds to edgebounds */
         EDGEBOUNDS_Pushback( edg, &( (BOUND){d_0,lb_0,rb_0} ) );

         /* Reorient new bounds and integrate it into edgebound list */
         // EDGEBOUNDS_Integrate( edg, (BOUND){d_0, lb_0, rb_0} );
      }

      /* MAIN RECURSION */
      /* Iterate the ranges of the antidiagonal */
      for ( b = 0; b < lb_vec[0]->N; b++ ) 
      {
         lb_0 = lb_vec[0]->data[b];
         rb_0 = rb_vec[0]->data[b];

         /* Iterate through cells in range */
         for ( k = lb_0; k < rb_0; k++ )
         {
            /* get x-y coords */
            i = k;
            j = d_0 - i;

            a = seq[i];
            A = AA_REV[a];

            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            prev_mat = MMX(i-1,j-1)  + TSC(j-1,M2M);
            prev_ins = IMX(i-1,j-1)  + TSC(j-1,I2M);
            prev_del = DMX(i-1,j-1)  + TSC(j-1,D2M);
            /* free to begin match state (new alignment) */
            // prev_beg = 0; /* assigned once at start */
            /* best-to-match */
            prev_sum = logsum( 
                           logsum( prev_mat, prev_ins ),
                           logsum( prev_del, prev_beg ) );
            MMX(i,j) = prev_sum + MSC(j,A);

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX(i-1,j) + TSC(j,M2I);
            prev_ins = IMX(i-1,j) + TSC(j,I2I);
            /* best-to-insert */
            prev_sum = logsum( prev_mat, prev_ins );
            IMX(i,j) = prev_sum + ISC(j,A);

            /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX(i,j-1) + TSC(j-1,M2D);
            prev_del = DMX(i,j-1) + TSC(j-1,D2D);
            /* best-to-delete */
            prev_sum = logsum(prev_mat, prev_del);
            DMX(i,j) = prev_sum;

            /* set each cell accessed */
            #if DEBUG
               MX_2D( cloud_MX, i, j ) = 1.0;
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
   }

   /* free dynamic memory */
   for (i = 0; i < 3; i++) {
      VECTOR_INT_Destroy( lb_vec[i] );
      VECTOR_INT_Destroy( rb_vec[i] );
   }
   
   /* show visualization of search cloud */
   #if DEBUG
   {
      DP_MATRIX_VIZ_Trace( cloud_MX, tr );
   }
   #endif 

   /* close necessary debugger tools */
   #if DEBUG
   {
      fclose(dbfp);
      fprintf(stdout, "==> ...Cloud Forward Quadratic Completed.\n");
   }
   #endif

   return total_max;
}

/*  
 *  FUNCTION: cloud_backward_Run()
 *  SYNOPSIS: Perform Backward part of Cloud Search Algorithm.
 */
float cloud_Backward_Quad( const SEQUENCE*     query,          /* query sequence */
                           const HMM_PROFILE*   target,        /* target model */
                           const int            Q,             /* query length */
                           const int            T,             /* target length */
                           MATRIX_3D*           st_MX,         /* normal state matrix, dim: ( NUM_NORMAL_STATES, Q+1, T+1 ) */
                           MATRIX_2D*           sp_MX,         /* special state matrix, dim: ( NUM_SPECIAL_STATES, Q+1 ) */
                           const ALIGNMENT*     tr,            /* viterbi traceback */ 
                           EDGEBOUNDS*          edg,           /* OUTPUT: edgebounds of cloud search space */
                           float                alpha,         /* PARAM: x-drop threshold */
                           int                  beta )         /* PARAM: free passes */
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

   /* debugging matrix */
   #if DEBUG
   {
      dbfp     = fopen( debugger->dbfp_path, "w+" );
      cloud_MX = debugger->cloud_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      test_MX  = debugger->test_MX;
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* clear leftover data */
   DP_MATRIX_Fill( Q, T, st_MX, sp_MX, -INF );

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
   prev_beg = 0;

   /* ITERATE THROUGHT ANTI-DIAGONALS */
   for ( d = d_end; d >= d_st; d--, d_cnt++ )
   {
      d_0 = d;       /* current diagonal */
      d_1 = (d+1);   /* look back 1 diagonal */
      d_2 = (d+2);   /* look back 2 diagonals */
      /* mod-mapping of antidiagonals into linear space */
      d0  = d_0 % 3; 
      d1  = d_1 % 3;
      d2  = d_2 % 3;

      /* Is dp matrix diagonal growing or shrinking? */
      if ( d >= dim_max )
         num_cells++;
      if ( d < dim_min )
         num_cells--;

      /* Edgechecks: find antidiag cell range that is inside matrix bounds */
      le_0 = MAX(end->i - (d_end - d_0), 0);
      re_0 = le_0 + num_cells;
      /* TODO: closed form? */

      /* Prune bounds */
      prune_via_xdrop_edgetrim_Quad( st_MX, sp_MX, alpha, beta, d_1, d_0, d1, d0, d_cnt, le_0, re_0, &total_max, lb_vec, rb_vec );

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

         /* Add new bounds to edgebounds */
         EDGEBOUNDS_Pushback( edg, &( (BOUND){d_0,lb_0,rb_0} ) );

         /* Reorient new bounds and integrate it into edgebound list */
         // EDGEBOUNDS_Integrate( edg, (BOUND){d_0, lb_0, rb_0} );
      }

      /* MAIN RECURSION */
      for ( b = 0; b < lb_vec[0]->N; b++ )
      {
         lb_0 = lb_vec[0]->data[b];
         rb_0 = rb_vec[0]->data[b];

         /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
         for ( k = lb_0; k < rb_0; k++ )
         {
            /* get x-y coords */
            i = k;
            j = d_0 - i;

            /* next sequence character */
            a = seq[i];
            A = AA_REV[a];
            
            /* match and insertion scores */
            sc_M = MSC(j+1, A);
            sc_I = ISC(j+1, A);

            /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
            prev_mat = MMX(i+1, j+1) + TSC(j, M2M) + sc_M;
            prev_ins = IMX(i+1, j  ) + TSC(j, M2I) + sc_I;
            prev_del = DMX(  i, j+1) + TSC(j, M2D);
            // prev_end = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
            // prev_end = sc_E;
            /* best-to-match */
            prev_sum = logsum( 
                           logsum( prev_mat, prev_ins ),
                           logsum( prev_del, prev_end ) );
            MMX(i,j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            prev_mat = MMX(i+1, j+1) + TSC(j, I2M) + sc_M;
            prev_ins = IMX(i+1, j)   + TSC(j, I2I) + sc_I;
            /* best-to-insert */
            prev_sum = logsum( prev_mat, prev_ins );
            IMX(i,j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            prev_mat = MMX(i+1, j+1) + TSC(j, D2M) + sc_M;
            prev_del = DMX(  i, j+1) + TSC(j, D2D);
            /* best-to-delete */
            prev_sum = logsum( prev_mat, prev_del );
            prev_sum = logsum( prev_sum, prev_end );
            DMX(i,j) = prev_sum;

            /* embed cell data in quadratic matrix */
            #if DEBUG
            {
               MX_2D( cloud_MX, i, j ) = 1.0;
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
   }

   /* reverse order of diagonals */
   EDGEBOUNDS_Reverse(edg);

   /* free dynamic memory */
   for (i = 0; i < 3; i++) {
      VECTOR_INT_Destroy( lb_vec[i] );
      VECTOR_INT_Destroy( rb_vec[i] );
   }

   /* show visualization of search cloud */
   #if DEBUG
   {
      DP_MATRIX_VIZ_Trace( cloud_MX, tr );
   }
   #endif 

   /* close necessary debugger tools */
   #if DEBUG
   {
      fclose(dbfp);
      fprintf(stdout, "==> ...Cloud Backward Quadratic Completed.\n" );
   }
   #endif

   return total_max;
} 