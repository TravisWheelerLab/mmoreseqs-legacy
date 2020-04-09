/*******************************************************************************
 *  FILE:      cloud_search_linear.h
 *  PURPOSE:   Cloud Search for Forward-Backward Pruning Algorithm. (LINEAR SPACE)
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

/* data stuctures and utility functions */
#include "objects/structs.h"
#include "utilities/utility.h"
#include "testing.h"

/* objects */
#include "objects/edgebound.h"
#include "objects/hmm_profile.h"
#include "objects/sequence.h"
#include "objects/alignment.h"
#include "objects/vectors/vector_range_2d.h"

/* file parsers */
#include "hmm_parser.h"

/* header */
#include "cloud_search_linear_rows.h"

/* 
 *       NOTE: CONVERSION - row-coords => diag-coords
 *       MX(i-1, j-1) => MX3(d_2, k-1)
 *       MX(i  , j-1) => MX3(d_1, k  ) 
 *       MX(i-1, j  ) => MX3(d_1, k-1)
 *
 *       MX(i+1, j+1) => MX3(d_2, k+1)
 *       MX(i  , j+1) => MX3(d_1, k  )
 *       MX(i+1, j  ) => MX3(d_1, k+1)
 */


/* ****************************************************************************************** *
 *  
 *  FUNCTION: cloud_Forward_Linear_Rows()
 *  SYNOPSIS: Perform Forward part of Cloud Search Algorithm.
 *            Traverses the dynamic programming matrix antidiagonally, running the
 *            Forward algorithm, starting at the Viterbi alignment beginning.  
 *            At the end of each antidiagonal, compares each cell against the current maximum
 *            scoring cell.  If cell falls below (MAX_SCORE * alpha), then cell is removed 
 *            from search space.  Terminates when reaches the final cell in dp matrix or 
 *            all cells in current antidiag have been pruned.  
 *            Stores final edgebound data in <edg>.
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix (quadratic space),
 *             <st_MX3>    Normal State (Match, Insert, Delete) Matrix (linear space),
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <tr>        Traceback Data,
 *             <edg>       (OUTPUT) Edgebounds Tracker Data,
 *             <alpha>     Pruning Drop,
 *             <beta>      Number of Passes before Pruning Begins
 *
 *  RETURN:    No Return.
 *
/* ****************************************************************************************** */
void cloud_Forward_Linear_Rows(  const SEQUENCE*      query, 
                                 const HMM_PROFILE*   target,
                                 const int            Q, 
                                 const int            T, 
                                 float*               st_MX, 
                                 float*               st_MX3,
                                 float*               sp_MX, 
                                 const ALIGNMENT*     tr,
                                 EDGEBOUNDS*          edg,
                                 const float          alpha, 
                                 const int            beta,
                                 const bool           test )
{
   /* vars for navigating matrix */
   int d,i,j,k;                        /* diagonal, row, column indices */
   /* NOTE: all left inclusive, right exclusive indexing!! */
   int lb, rb, le, re;                 /* right/left bounds and edges */
   int lb_1, rb_1, lb_2, rb_2;         /* for storing lookback bounds */
   int lb_new, rb_new;                 /* tmp for new right/left bounds */
   int *lb_list, rb_list;              /* tmp for forking right/left bounds */
   int num_cells;                      /* number of cells in diagonal */
   int d_st, d_end, d_cnt;             /* starting and ending diagonal indices */
   int dim_min, dim_max;               /* diagonal index where num cells reaches highest point and diminishing point */ 
   int dim_T, dim_Q, dim_TOT;          /* dimensions of submatrix being searched */
   int total_cnt;                      /* number of cells computed */
   TRACE beg, end;                    /* start and end point of alignment */

   /* storing cloud results */
   EDGEBOUNDS*       edg_diags;         /* stores the cloud space in diagonals */
   VECTOR_RANGE_2D*  vec_rows;          /* stores the cloud space in rows */

   /* vars for computing cells */
   char   a;                           /* store current character in sequence */
   int    A;                           /* store int value of character */
   char   *seq = query->seq;           /* alias for getting seq */

   /* vars for recurrance */
   int    d_0, d_1, d_2;               /* for assigning prev array ptrs */
   int    d0,  d1,  d2;                /* for assigning prev array ptrs in mod for linear space */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;

   /* vars for pruning */
   float  cell_max, total_max, total_limit, diag_max, diag_limit;

   /* local or global? (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* track number of cells computed */
   int cpu_num = 0;

   /* --------------------------------------------------------------------------------- */

   /* clear residual data */
   if (test) dp_matrix_Clear(Q, T, st_MX, sp_MX);
   dp_matrix_Clear3(Q, T, st_MX3, sp_MX);

   edg_diags = EDGEBOUNDS_Create();
   vec_rows  = VECTOR_RANGE_2D_Create_Size(T + 1);

   /* get start and end points of viterbi alignment */
   beg = tr->traces[tr->beg];
   end = tr->traces[tr->end];

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (beg.i == 0 || beg.j == 0) {
      beg.i += 1;
      beg.j += 1;
   }

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   /* diag index of different start points, creating submatrix */
   d_st = beg.i + beg.j;
   // d_end = end.i + end.j;

   /* dimension of submatrix */
   dim_TOT = Q + T;
   dim_Q = Q - beg.i;
   dim_T = T - beg.j;

   /* diag index where num cells reaches highest point and begins diminishing */
   dim_min = MIN(d_st + dim_Q, d_st + dim_T);
   dim_max = MAX(d_st + dim_Q, d_st + dim_T);

   /* set bounds using starting cell */
   lb = beg.i;
   rb = beg.i + 1;
   num_cells = 0;
   lb_1 = rb_1 = lb_2 = rb_2 = 0;

   /* keeps largest number seen on current diagonal */
   diag_max = -INF;
   total_max = -INF;

   /* begin state probability begins at 0 */
   prev_beg = 0;

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
      if (d <= dim_min)
         num_cells++;
      if (d > dim_max)
         num_cells--;

      /* UPDATE/PRUNE BOUNDS */
      /* if free passes are complete (beta < d), prune and set new edgebounds */
      /* beta must be > 1 */
      if (beta < d_cnt)
      {
         /* impossible state */
         lb_new = INT_MIN;    
         rb_new = INT_MIN;

         /* FIND MAX SCORE ON CURRENT DIAGONAL */
         diag_max = -INF;
         for (k = lb_1; k < rb_1; k++)
         {
            /* coords for quadratic matrix */
            i = k;
            j = d_1 - i;    /* back one diag */

            diag_max = calc_Max( 
                           calc_Max( diag_max, MMX3(d1,k) ),
                           calc_Max( IMX3(d1,k), DMX3(d1,k) ) );
         }

         /* Total max records largest cell score seen so far */
         total_max = MAX(total_max, diag_max);

         /* Set score threshold for pruning */
         diag_limit = diag_max - alpha;
         total_limit = total_max - alpha;
         // printf("total_max: %.2f\t total_limit: %.2f\t diag_max: %.2f\t diag_limit: %.2f\n", total_max, total_limit, diag_max, diag_limit);

         /* FIND FIRST SCORE TO EXCEED THRESHOLD FROM THE LEFT */
         for (k = lb_1; k < rb_1; k++)
         {
            /* coords for quadratic matrix */
            i = k;
            j = d_1 - i; /* looking back one diag */

            cell_max = calc_Max( MMX3(d1,k),
                           calc_Max( IMX3(d1,k), DMX3(d1,k) ) );

            /* prune in left edgebound */
            if( cell_max >= total_limit ) {
               lb_new = i;
               break;
            }
         }

         /* If no boundary edges are found on diag, then branch is pruned entirely and we are done */
         if (lb_new == INT_MIN)
            break;

         /* FIND FIRST SCORE TO EXCEED THRESHOLD FROM THE RIGHT */
         for (k = rb_1 - 1; k >= lb_1; k--)
         {
            /* coords for quadratic matrix */
            i = k;
            j = d_1 - i;

            cell_max = calc_Max( MMX3(d1,k),
                           calc_Max( IMX3(d1,k), DMX3(d1,k) ) );

            /* prune in right edgebound */
            if( cell_max >= total_limit )
            {
               rb_new = (i + 1);
               break;
            }
         }
         // printf("pruned -> (%d,%d)...\n", lb_new, rb_new);
      }
      else /* else edges expand in square pattern */
      {
         lb_new = lb;
         rb_new = rb;
      }

      /* Update bounds */
      lb = lb_new;
      rb = rb_new + 1;
      /* NOTE: FOR TESTING - THiS REMOVES ALL PRUNING */
      // lb = lb;
      // rb = rb + 1;

      /* Edge-checks: find if diag cells that are inside matrix bounds */
      le = MAX(beg.i, d - T);
      re = le + num_cells;

      /* Check that they dont exceed edges of matrix */
      lb = MAX(lb, le);
      rb = MIN(rb, re);

      /* ADD NEW ANTI-DIAG TO EDGEBOUNDS */
      EDGEBOUNDS_Pushback(edg, (BOUND){d,lb,rb});

      /* Add new edgebounds to rows and diags */
      EDGEBOUNDS_Pushback(edg_diags, (BOUND){d,lb,rb});
      VECTOR_RANGE_2D_MergeFwd(vec_rows, (BOUND){d,lb,rb});

      /* MAIN RECURSION */
      /* ITERATE THROUGH CELLS OF NEXT ANTI-DIAGONAL */
      for (k = lb; k < rb; k++, total_cnt++)
      {
         /* quadratic coords */
         i = k;
         j = d - i;

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
         prev_sum = calc_Logsum( 
                        calc_Logsum( prev_mat, prev_ins ),
                        calc_Logsum( prev_del, prev_beg ) );
         MMX3(d0,k) = prev_sum + MSC(j,A);

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: Convert (i-1,j) <=> (d-1,k-1) */
         prev_mat = MMX3(d1,k-1) + TSC(j,M2I);
         prev_ins = IMX3(d1,k-1) + TSC(j,I2I);
         /* best-to-insert */
         prev_sum = calc_Logsum( prev_mat, prev_ins );
         IMX3(d0,k) = prev_sum + ISC(j,A);

         /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: Convert (i,j-1) <=> (d-1, k) */
         prev_mat = MMX3(d1,k) + TSC(j-1,M2D);
         prev_del = DMX3(d1,k) + TSC(j-1,D2D);
         /* best-to-delete */
         prev_sum = calc_Logsum(prev_mat, prev_del);
         DMX3(d0,k) = prev_sum;
      }

      /* SCRUB 2-BACK ANTIDIAGONAL */
      for (k = lb_2; k < rb_2; k++) {
         MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = -INF;
      }
      rb_2 = rb_1;
      lb_2 = lb_1;
      rb_1 = rb;
      lb_1 = lb;

      // /* NAIVE SCRUB - FOR DEBUGGING */
      // for (k = 0; k < (T+1)+(Q+1); k++) {
      //    MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = -INF;
      // }

      /* Embed Current Row into Quadratic Array - FOR DEBUGGING */
      if (test) {
         for (k = le; k <= re; k++) {
            i = k;
            j = d_0 - i;

            MMX(i,j) = MMX3(d0,k);
            IMX(i,j) = IMX3(d0,k);
            DMX(i,j) = DMX3(d0,k);
         }
      }
   }
}


/* ****************************************************************************************** *
 *  
 *  FUNCTION: cloud_Backward_Linear()
 *  SYNOPSIS: Perform Backward part of Cloud Search Algorithm.
 *            Traverses the dynamic programming matrix antidiagonally, running the
 *            Forward algorithm, starting at the Viterbi alignment ending.  
 *            At the end of each antidiagonal, compares each cell against the current maximum
 *            scoring cell.  If cell falls below (MAX_SCORE * alpha), then cell is removed 
 *            from search space.  Terminates when reaches the final cell in dp matrix or 
 *            all cells in current antidiag have been pruned.  
 *            Stores final edgebound data in <edg>.
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix (quadratic space),
 *             <st_MX3>    Normal State (Match, Insert, Delete) Matrix (linear space),
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <tr>        Traceback Data,
 *             <edg>       (OUTPUT) Edgebounds Tracker Data,
 *             <alpha>     Pruning Drop,
 *             <beta>      Number of Passes before Pruning Begins
 *
 *  RETURN:    No Return.
 *
/* ****************************************************************************************** */
void cloud_Backward_Linear_Rows( const SEQUENCE*      query, 
                                 const HMM_PROFILE*   target,
                                 const int            Q, 
                                 const int            T, 
                                 float*               st_MX, 
                                 float*               st_MX3,
                                 float*               sp_MX, 
                                 const ALIGNMENT*     tr,
                                 EDGEBOUNDS*          edg,
                                 const float          alpha, 
                                 const int            beta,
                                 const bool           test )
{
   /* vars for navigating matrix */
   int d,i,j,k;                  /* diagonal, row, column, ... indices */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int lb_1, rb_1, lb_2, rb_2;   /* for storing lookback bounds */
   int lb_new, rb_new;           /* tmp vars for new right/left bounds */
   int num_cells;                /* number of cells in diagonal */
   int d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int dim_T, dim_Q;             /* dimensions of submatrix being searched */
   int total_cnt;                /* number of cells computed */
   TRACE beg, end;               /* start and end point of alignment */

   /* storing cloud results */
   EDGEBOUNDS*       edg_diags;         /* stores the cloud space in diagonals */
   VECTOR_RANGE_2D*  vec_rows;          /* stores the cloud space in rows */

   /* vars for computing cells */
   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */
   char   *seq = query->seq;     /* alias for getting seq */

   /* vars for recurrance */
   int    d_0, d_1, d_2;         /* for assigning prev array ptrs */
   int    d0,  d1,  d2;          /* for assigning prev array ptrs in mod for linear space */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;

   /* vars for pruning */
   float  cell_max, total_max, total_limit, diag_max, diag_limit;

   /* local or global? (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* track number of cells computed */
   int cpu_num = 0;

   /* --------------------------------------------------------------------------------- */

   /* clear residual data */
   if (test) dp_matrix_Clear(Q, T, st_MX, sp_MX);
   dp_matrix_Clear3(Q, T, st_MX3, sp_MX);

   edg_diags = EDGEBOUNDS_Create();
   vec_rows  = VECTOR_RANGE_2D_Create_Size(T + 1);

   /* get start and end points of viterbi alignment */
   beg = tr->traces[tr->beg];
   end = tr->traces[tr->end];

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (beg.i == Q+1 || beg.j == T+1) {
      beg.i -= 1;
      beg.j -= 1;
   }

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   /* diag index of different start points, creating submatrix */
   // d_st = beg.i + beg.j;
   d_end = end.i + end.j;

   /* dimension of submatrix */
   dim_Q = end.i;
   dim_T = end.j;

   /* diag index where num cells reaches highest point and begins diminishing */
   dim_min = MIN(end.i, end.j);
   dim_max = MAX(end.i, end.j);

   /* set bounds using starting cell */
   lb = end.i;
   rb = end.i + 1;
   num_cells = 0;
   lb_1 = rb_1 = lb_2 = rb_2 = 0;

   /* end state starts at 0 */
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

      /* UPDATE/PRUNE BOUNDS */
      /* if free passes are complete (beta < d), prune and set new edgebounds */
      /* beta must be > 1 */
      if (beta < d_cnt)
      {
         /* impossible state */
         lb_new = INT_MIN;    
         rb_new = INT_MIN;

         /* FIND MAX SCORE ON CURRENT DIAGONAL */
         diag_max = -INF;
         for (k = lb_1; k < rb_1; k++)
         {
            /* cartesian coords */
            i = k;
            j = d_1 - i;    /* back one diag */

            diag_max = calc_Max( 
                           calc_Max( diag_max, MMX3(d1,k) ),
                           calc_Max( IMX3(d1,k), DMX3(d1,k) ) );
         }

         /* Total max records largest cell score seen so far */
         total_max = MAX(total_max, diag_max);

         /* Set score threshold for pruning */
         diag_limit = diag_max - alpha;
         total_limit = total_max - alpha;

         /* FIND FIRST SCORE TO EXCEED THRESHOLD FROM THE LEFT */
         for (k = lb_1; k < rb_1; k++)
         {
            /* coords for quadratic matrix */
            i = k;
            j = d_1 - i; /* looking back one diag */

            cell_max = calc_Max( MMX3(d1,k),
                           calc_Max( IMX3(d1,k), DMX3(d1,k) ) );

            /* prune in left edgebound */
            if( cell_max >= total_limit ) {
               lb_new = i;
               break;
            }
         }

         /* If no boundary edges are found on diag, then branch is pruned entirely and we are done */
         if (lb_new == INT_MIN)
            break;

         /* FIND FIRST SCORE TO EXCEED THRESHOLD FROM THE RIGHT */
         for (k = rb_1 - 1; k >= lb_1; k--)
         {
            /* coords for quadratic matrix */
            i = k;
            j = d_1 - i;

            cell_max = calc_Max( MMX3(d1,k),
                           calc_Max( IMX3(d1,k), DMX3(d1,k) ) );

            /* prune in right edgebound */
            if( cell_max >= total_limit )
            {
               rb_new = (i + 1);
               break;
            }
         }
      }
      else /* else edges expand in square pattern */
      {
         lb_new = lb;
         rb_new = rb;
      }

      /* Update Bounds */
      lb = lb_new - 1;
      rb = rb_new;
      // /* NOTE: TESTING - THiS REMOVES ALL PRUNING */
      // lb = lb - 1;
      // rb = rb;

      /* Edge-check: find diag cells that are inside matrix bounds */
      le = MAX(end.i - (d_end - d), 0);
      re = le + num_cells;

      lb = MAX(lb, le);
      rb = MIN(rb, re);

      /* ADD NEW ANTI-DIAG TO EDGEBOUNDS */
      EDGEBOUNDS_Pushback(edg, (BOUND){d,lb,rb});

      /* Add new edgebounds to rows and diags */
      EDGEBOUNDS_Pushback(edg_diags, (BOUND){d,lb,rb});
      VECTOR_RANGE_2D_MergeFwd(vec_rows, (BOUND){d,lb,rb});

      /* MAIN RECURSION */
      /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
      for (k = lb; k < rb; k++, total_cnt++)
      {
         i = k;
         j = d - i;

         /*
          *    MX(i+1, j+1) => MX3(d_2, k+1)
          *    MX(i  , j+1) => MX3(d_1, k  )
          *    MX(i+1, j  ) => MX3(d_1, k+1)
          */

         /* Get next sequence character */
         a = seq[i];
         A = AA_REV[a];
         
         sc_M = MSC(j+1,A);
         sc_I = ISC(j+1,A);

         /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
         prev_mat = MMX3(d2,k+1) + TSC(j,M2M) + sc_M;
         prev_ins = IMX3(d1,k+1) + TSC(j,M2I) + sc_I;
         prev_del = DMX3(d1,k  ) + TSC(j,M2D);
         // prev_end = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
         // prev_end = sc_E;
         /* best-to-match */
         prev_sum = calc_Logsum( 
                           calc_Logsum( prev_mat, prev_ins ),
                           calc_Logsum( prev_del, prev_end ) );
         MMX3(d0,k) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
         prev_mat = MMX3(d2,k+1) + TSC(j,I2M) + sc_M;
         prev_ins = IMX3(d1,k+1) + TSC(j,I2I) + sc_I;
         /* best-to-insert */
         prev_sum = calc_Logsum( prev_mat, prev_ins );
         IMX3(d0,k) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
         prev_mat = MMX3(d2,k+1) + TSC(j,D2M) + sc_M;
         prev_del = DMX3(d1,k  ) + TSC(j,D2D);
         /* best-to-delete */
         prev_sum = calc_Logsum( prev_mat, prev_del );
         prev_sum = calc_Logsum( prev_sum, prev_end );
         DMX3(d0,k) = prev_sum;
      }

      /* SCRUB 2-BACK ANTIDIAGONAL */
      for (k = lb_2; k < rb_2; k++) {
         MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = -INF;
      }
      rb_2 = rb_1;
      lb_2 = lb_1;
      rb_1 = rb;
      lb_1 = lb;

      // /* NAIVE SCRUB - FOR DEBUGGING */
      // for (k = 0; k < (T+1)+(Q+1); k++) {
      //    MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = -INF;
      // }

      /* Embed Current Row into Quadratic Array - FOR DEBUGGING */
      if (test) {
         for (k = le; k <= re; k++) {
            i = k;
            j = d_0 - i;

            MMX(i,j) = MMX3(d0,k);
            IMX(i,j) = IMX3(d0,k);
            DMX(i,j) = DMX3(d0,k);
         }
      }
   }
   /* reverse order of diagonals */
   EDGEBOUNDS_Reverse(edg);
}
