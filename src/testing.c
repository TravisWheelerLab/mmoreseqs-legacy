/*******************************************************************************
 *  @file cloud_search.h
 *  @brief Testing for navigating through the matrices.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

// imports
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

// local imports (after struct declarations)
#include "structs.h"
#include "misc.h"
#include "hmm_parser.h"

/* test to cycle through all diags */
void fwd_test_cycle(int Q, int T,
                   float* st_MX,
                   float* sp_MX,
                   TRACEBACK* tr)
{
   /* vars for navigating matrix */
   int d,i,j,k,x;                  /* diagonal, row, column indices */
   /* NOTE: all left inclusive, right exclusive indexing!! */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int lb_new, rb_new;           /* tmp for new right/left bounds */\
   int num_cells;                /* number of cells in diagonal */
   int d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int dim_T, dim_Q, dim_TOT;    /* dimensions of submatrix being searched */
   int total_cnt;                /* number of cells computed */
   COORDS start, end;            /* start and end point of alignment */

   /* vars for computing cells */
   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */

   /* vars for recurrance */
   int    d_0, d_1, d_2;         /* for assigning prev array ptrs */
   int    d0,  d1,  d2;          /* for assigning prev array ptrs in mod for linear space */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;

   /* vars for pruning */
   float  cell_max, total_max, total_limit, diag_max, diag_limit;

   /* track number of cells computed */
   int cpu_num = 0;

   /* --------------------------------------------------------------------------------- */

   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);

   /* get start and end points of viterbi alignment */
   start = tr->first_m;
   end = tr->last_m;

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (start.i == 0 || start.j == 0) {
      start.i += 1;
      start.j += 1;
   }

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   /* diag index of different start points, creating submatrix */
   d_st = start.i + start.j;
   // d_end = end.i + end.j;

   /* dimension of submatrix */
   dim_TOT = Q + T;
   dim_Q = Q - start.i;
   dim_T = T - start.j;

   /* diag index where num cells reaches highest point and begins diminishing */
   // dim_min = MIN(Q + start.i, T + start.j);
   // dim_max = MAX(Q + start.i, T + start.j);
   dim_min = MIN(d_st + dim_Q, d_st + dim_T);
   dim_max = MAX(d_st + dim_Q, d_st + dim_T);

   // set bounds using starting cell
   lb = start.i;
   rb = start.i + 1;
   num_cells = 0;

   /* iterate through diags */
   for (d = d_st; d <= d_end; d++, rb++)
   {
      d_0 = d;          /* current antidiagonal */
      d_1 = (d-1);      /* look back 1 antidiagonal */
      d_2 = (d-2);      /* look back 2 antidiagonal */
      /* mod-mapping of antidiagonals into linear space */
      d0  = d_0; 
      d1  = d_1;
      d2  = d_2;

      /* is dp matrix diagonal growing or shrinking? */
      if (d <= dim_min)
         num_cells++;
      if (d > dim_max)
         num_cells--;

      /* find diag cells that are inside matrix bounds */
      le = MAX(start.i, d - T);
      re = le + num_cells;

      /* NOTE: TESTING - THiS REMOVES ALL PRUNING */
      lb = lb;
      rb = rb + 1;

      /* for testing, bounds = edges */
      lb = MAX(lb, le);
      rb = MIN(rb, re);

      /* iterate through cells of diag */
      for (k = lb; k < rb; k++)
      {
         /* quadratic coords */
         i = k;
         j = d - i;
         /* linear/antidiag coords */
         /* NOTE: For embedding in (T+1) length array, need to offset to keep in bounds. */
         /* NOTE: But this requires changing the coords of lookback antidiags at dim_min point */
         /* ex. x = k - le */
         x = k;

         if ( MMX(i,j) != 0 || DMX(i,j) != 0 || IMX(i,j) != 0)
            printf("OVERWRITE AT (%d,%d)!\n", i, j);

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         /* NOTE: (i-1,j-1) => (d-2,k-1) */ 
         MMX(i,j) = MMX(i-1,j-1) + 1;

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i-1,j) => (d-1,k-1) */
         IMX(i,j) = IMX(i-1,j) + 1;

         /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i,j-1) => (d-1, k) */
         DMX(i,j) = DMX(i,j-1) + 1;
      }
   }
}

void fwd_test_cycle3(int Q, int T, 
                     float* st_MX, 
                     float* st_MX3,
                     float* sp_MX, 
                     TRACEBACK* tr )
{
   /* vars for navigating matrix */
   int d,i,j,k,x;                  /* diagonal, row, column indices */
   /* NOTE: all left inclusive, right exclusive indexing!! */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int lb_new, rb_new;           /* tmp for new right/left bounds */\
   int num_cells;                /* number of cells in diagonal */
   int d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int dim_T, dim_Q, dim_TOT;    /* dimensions of submatrix being searched */
   int total_cnt;                /* number of cells computed */
   COORDS start, end;            /* start and end point of alignment */
   bool overwrite;

   /* vars for computing cells */
   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */

   /* vars for recurrance */
   int    d_0, d_1, d_2;         /* for assigning prev array ptrs */
   int    d0,  d1,  d2;          /* for assigning prev array ptrs in mod for linear space */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;

   /* vars for pruning */
   float  cell_max, total_max, total_limit, diag_max, diag_limit;

   /* track number of cells computed */
   int cpu_num = 0;

   /* --------------------------------------------------------------------------------- */

   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   dp_matrix_Clear_X3(Q, T, st_MX3, sp_MX, 0);

   /* get start and end points of viterbi alignment */
   start = tr->first_m;
   end = tr->last_m;

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (start.i == 0 || start.j == 0) {
      start.i += 1;
      start.j += 1;
   }

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   /* diag index of different start points, creating submatrix */
   d_st = start.i + start.j;
   // d_end = end.i + end.j;

   /* dimension of submatrix */
   dim_TOT = Q + T;
   dim_Q = Q - start.i;
   dim_T = T - start.j;

   /* diag index where num cells reaches highest point and begins diminishing */
   dim_min = MIN(d_st + dim_Q, d_st + dim_T);
   dim_max = MAX(d_st + dim_Q, d_st + dim_T);

   /* set bounds using starting cell */
   lb = start.i;
   rb = start.i + 1;
   num_cells = 0;

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

      /* NOTE: TESTING - THiS REMOVES ALL PRUNING */
      lb = lb;
      rb = rb + 1;

      /* Edge-checks: find if diag cells that are inside matrix bounds */
      le = MAX(start.i, d - T);
      re = le + num_cells;

      /* Check that they dont exceed edges of matrix */
      lb = MAX(lb, le);
      rb = MIN(rb, re);

      /* MAIN RECURSION */
      /* ITERATE THROUGH CELLS OF NEXT ANTI-DIAGONAL */
      for (k = lb; k < rb; k++, total_cnt++)
      {
         /* quadratic coords */
         i = k;
         j = d - i;
         /* linear/antidiag coords */
         /* NOTE: For embedding in (T+1) length array, need to offset to keep in bounds. */
         /* NOTE: But this requires changing the coords of lookback antidiags at dim_min point */
         /* ex. x = k - le */
         x = k;

         if ( MMX3(d0,x) != 0 || DMX3(d0,x) != 0 || IMX3(d0,x) != 0) {
            printf("OVERWRITE AT %d:(%d,%d)=%d!\n", d, i, j, DMX(d0,x) );
            overwrite = true;
         }

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         /* NOTE: (i-1,j-1) => (d-2,k-1) */ 
         MMX3(d0,x) = MMX3(d2,k-1) + 1;

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i-1,j) => (d-1,k-1) */
         IMX3(d0,x) = IMX3(d1,k-1) + 1;

         /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i,j-1) => (d-1, k) */
         DMX3(d0,x) = DMX3(d1,k) + 1;

         if (overwrite) {
            overwrite = false;
            MMX3(d0,x) = 30;
         }
      }

      /* NAIVE SCRUB */
      for (k = 0; k < ((T+1)+(Q+1)); k++) {
         MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = 0;
      }

      /* Embed Current Row into Quadratic Array - FOR DEBUGGING */
      for (k = le; k < re; k++) {
         /* quadratic coords */
         i = k;
         j = d_0 - i;
         /* linear coords */
         x = k;

         MMX(i,j) = MMX3(d0,x);
         IMX(i,j) = IMX3(d0,x);
         DMX(i,j) = DMX3(d0,x);
      }
   }
}

/* test to cycle through all diags in reverse */
void bck_test_cycle(int Q, int T,
                    float* st_MX,
                    float* sp_MX,
                    TRACEBACK* tr)
{
   /* vars for navigating matrix */
   int d,i,j,k,x;                  /* diagonal, row, column indices */
   /* NOTE: all left inclusive, right exclusive indexing!! */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int lb_new, rb_new;           /* tmp for new right/left bounds */\
   int num_cells;                /* number of cells in diagonal */
   int d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int dim_T, dim_Q, dim_TOT;    /* dimensions of submatrix being searched */
   int total_cnt;                /* number of cells computed */
   COORDS start, end;            /* start and end point of alignment */

   /* vars for computing cells */
   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */

   /* vars for recurrance */
   int    d_0, d_1, d_2;         /* for assigning prev array ptrs */
   int    d0,  d1,  d2;          /* for assigning prev array ptrs in mod for linear space */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;

   /* vars for pruning */
   float  cell_max, total_max, total_limit, diag_max, diag_limit;

   /* track number of cells computed */
   int cpu_num = 0;

   /* --------------------------------------------------------------------------------- */

   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);

   /* INIT ANTI-DIAGS */
   /* test coords */
   start = tr->first_m;
   end = tr->last_m;

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (start.i == Q+1 || start.j == T+1) {
      start.i -= 1;
      start.j -= 1;
   }

   // diag index at corners of dp matrix
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   // diag index of different start points, creating submatrix
   // d_st = start.i + start.j;
   d_end = end.i + end.j;

   // dimension of submatrix
   dim_TOT = Q + T;
   dim_Q = end.i;
   dim_T = end.j;

   // diag index where num cells reaches highest point and begins diminishing
   dim_min = MIN(end.i, end.j);
   dim_max = MAX(end.i, end.j);

   // set bounds using starting cell
   lb = end.i;
   rb = end.i + 1;
   num_cells = 0;

   /* ITERATE THROUGHT ANTI-DIAGONALS */
   for (d = d_end; d >= d_st; d--)
   {
      d_0 = d;       /* current diagonal */
      d_1 = (d+1);   /* look back 1 diagonal */
      d_2 = (d+2);   /* look back 2 diagonals */
      /* mod-mapping of antidiagonals into linear space */
      d0  = d_0; 
      d1  = d_1;
      d2  = d_2;

      /* is dp matrix diagonal growing or shrinking? */
      if (d >= dim_max)
         num_cells++;
      if (d < dim_min)
         num_cells--;

      /* find diag cells that are inside matrix bounds */
      le = MAX(end.i - (d_end - d), 0);
      re = le + num_cells;

      /* NOTE: TESTING - THiS REMOVES ALL PRUNING */
      lb = lb - 1;
      rb = rb;

      /* for testing, bounds = edges */
      lb = MAX(lb, le);
      rb = MIN(rb, re);

      /* MAIN RECURSION */
      /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
      for (k = lb; k < rb; k++)
      {
         /* cartesian coords */
         i = k;
         j = d - i;
         /* antidiag coords */
         /* NOTE: For embedding in (T+1) length array, need to offset to keep in bounds. */
         /* NOTE: But this requires changing the coords of lookback antidiags at dim_min point */
         /* ex. x = k - le */
         x = k;

         if ( MMX(i,j) != 0 || DMX(i,j) != 0 || IMX(i,j) != 0)
            printf("OVERWRITE AT (%d,%d)!\n", i, j);

         /*
          *   MAT: MX(i+1, j+1) => MX3(d_2, k+1)
          *   DEL: MX(i  , j+1) => MX3(d_1, k  )
          *   INS: MX(i+1, j  ) => MX3(d_1, k+1)
          */

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         /* NOTE: (i+1,j+1) => (d+2, k+1) */ 
         MMX(i,j) = MMX(i+1,j+1) + 1;

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i+1,j) => (d+1, k  ) */
         IMX(i,j) = IMX(i+1,j) + 1;

         /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i,j+1) => (d+1, k+1) */
         // DMX(i,j) = DMX(i,j+1) + 1;
         DMX(i,j) += d;
      }
   }
}

/* test to cycle through all diags in reverse */
void bck_test_cycle3(int Q, int T,
                    float* st_MX,
                    float* st_MX3,
                    float* sp_MX,
                    TRACEBACK* tr)
{
   /* vars for navigating matrix */
   int d,i,j,k,x;                /* diagonal, row, column, ... indices */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int lb_1, rb_1, lb_2, rb_2;   /* for storing lookback bounds */
   int lb_new, rb_new;           /* tmp vars for new right/left bounds */
   int num_cells;                /* number of cells in diagonal */
   int d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int dim_T, dim_Q, dim_TOT;    /* dimensions of submatrix being searched */
   int total_cnt;                /* number of cells computed */
   COORDS start, end;            /* start and end point of alignment */
   bool overwrite;

   /* vars for computing cells */
   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */

   /* vars for recurrance */
   int    d_0, d_1, d_2;         /* for assigning prev array ptrs */
   int    d0,  d1,  d2;          /* for assigning prev array ptrs in mod for linear space */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;

   /* vars for pruning */
   float  cell_max, total_max, total_limit, diag_max, diag_limit;

   /* track number of cells computed */
   int cpu_num = 0;

   /* --------------------------------------------------------------------------------- */

   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   dp_matrix_Clear_X3(Q, T, st_MX3, sp_MX, 0);

   /* get start and end points of viterbi alignment */
   start = tr->first_m;
   end = tr->last_m;

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (start.i == Q+1 || start.j == T+1) {
      start.i -= 1;
      start.j -= 1;
   }

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   /* diag index of different start points, creating submatrix */
   // d_st = start.i + start.j;
   d_end = end.i + end.j;

   /* dimension of submatrix */
   dim_TOT = Q + T;
   dim_Q = end.i;
   dim_T = end.j;

   /* diag index where num cells reaches highest point and begins diminishing */
   dim_min = MIN(end.i, end.j);
   dim_max = MAX(end.i, end.j);

   /* set bounds using starting cell */
   lb = end.i;
   rb = end.i + 1;
   num_cells = 0;

   /* end state starts at 0 */
   prev_end = 0;

   /* ITERATE THROUGHT ANTI-DIAGONALS */
   for (d = d_end; d >= d_st; d--)
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

      /* NOTE: TESTING - THiS REMOVES ALL PRUNING */
      lb = lb - 1;
      rb = rb;

      /* Edge-check: find diag cells that are inside matrix bounds */
      le = MAX(end.i - (d_end - d), 0);
      re = le + num_cells;

      lb = MAX(lb, le);
      rb = MIN(rb, re);

      /* MAIN RECURSION */
      /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
      for (k = lb; k < rb; k++)
      {
         i = k;
         j = d - i;
         /* linear/antidiag coords */
         /* NOTE: For embedding in (T+1) length array, need to offset to keep in bounds. */
         /* NOTE: But this requires changing the coords of lookback antidiags at dim_min point */
         /* ex. x = k - le */
         x = k;

         if ( MMX3(d0,x) != 0 || DMX3(d0,x) != 0 || IMX3(d0,x) != 0) {
            printf("OVERWRITE AT %d:(%d,%d)=%d!\n", d, i, j, DMX(d0,x) );
            overwrite = true;
         }

         /*
          *   MAT: MX(i+1, j+1) => MX3(d_2, k+1)
          *   DEL: MX(i  , j+1) => MX3(d_1, k  )
          *   INS: MX(i+1, j  ) => MX3(d_1, k+1)
          */

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         /* NOTE: (i+1,j+1) => (d+2,k+1) */ 
         MMX3(d0,x) = MMX3(d2,k+1) + 1;

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i+1,j) => (d+1,k+1) */
         IMX3(d0,x) = IMX3(d1,k+1) + 1;

         /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i,j+1) => (d+1, k) */
         // DMX3(d0,x) = DMX3(d1,k) + 1;
         DMX3(d0,x) += d;

         if (overwrite) {
            overwrite = false;
            // MMX3(d0,x) = 30;
         }
      }

      // /* SCRUB 2-BACK ANTIDIAGONAL */
      // for (k = lb_2; k < rb_2; k++) {
      //    MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = -INF;
      // }
      // rb_2 = rb_1;
      // lb_2 = lb_1;
      // rb_1 = rb;
      // lb_1 = lb;

      /* NAIVE SCRUB - FOR DEBUGGING */
      for (k = 0; k < (T+1)+(Q+1); k++) {
         MMX3(d2,k) = IMX3(d2,k) = DMX3(d2,k) = 0;
      }

      /* Embed Current Row into Quadratic Array - FOR DEBUGGING */
      for (k = le; k <= re; k++) {
         i = k;
         j = d_0 - i;

         MMX(i,j) = MMX3(d0,k);
         IMX(i,j) = IMX3(d0,k);
         DMX(i,j) = DMX3(d0,k);
      }
   }
}  



/* test to show the cloud area, fill with value, return number of used cells  */
int cloud_Fill(int Q, int T,
                float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                EDGEBOUNDS* edg,
                float val, 
                int mode)
{
   int x, d, i, j, k;
   int rb, lb;
   int num_cells = 0;
   int N = edg->N;
   int MAX = (DEL_ST*(Q+1)*(T+1)) + ((Q)*(T+1)) + (T);
   int IDX;

   if (mode == MODE_DIAG)
   {
      /* iterate over all bounds in edgebound list */
      for (x = 0; x < N; x++)
      {
         d  = edg->bounds[x].diag;
         lb = edg->bounds[x].lb;
         rb = edg->bounds[x].rb;

         /* insert value across diag in range */
         for (k = lb; k < rb; k++)
         {
            i = k;
            j = d - i;

            // if (i > Q || j > T ) 
            //    printf("OVER MAX: MAX=(Q=%d,T=%d) -> (%d,%d)\n", Q, T, i, j);

            MMX(i,j) += 1;
            // if (MMX(i,j) > 1)
            //    printf("DOUBLE HIT!  MAX=(Q=%d,T=%d) -> d:(%d,%d) = r:(%d,%d)\n", Q, T, d, k, i, j);

            IMX(i,j) = 1;
            DMX(i,j) += val;
            num_cells++;
         }
      }
   }
   else if (mode == MODE_ROW)
   {
      /* iterate over all bounds in edgebound list */
      for (x = 0; x < N; x++)
      {
         i  = edg->bounds[x].diag;
         lb = edg->bounds[x].lb;
         rb = edg->bounds[x].rb;

         /* insert value across row in range */
         for (j = lb; j < rb; j++)
         {
            IDX = (DEL_ST*(Q+1)*(T+1)) + ((i)*(T+1)) + (j);

            if (IDX > MAX)
               printf("OVER MAX: MAX=(%d,%d) -> (%d,%d)\n", Q, T, i, j);

            MMX(i,j) += 1;
            IMX(i,j) = 1;
            DMX(i,j) += val;
            num_cells++;
         }
      }
   }
   return num_cells;
}


/* test to show the cloud area, fill with value, return number of used cells  */
int cloud_Solid_Fill(int Q, int T,
                   float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                   float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                   EDGEBOUNDS* edg,
                   float val, 
                   int mode)
{
   int x, d, i, j, k;
   int rb, lb;
   int num_cells = 0;
   int N = edg->N;
   int MAX = (DEL_ST*(Q+1)*(T+1)) + ((Q)*(T+1)) + (T);
   int IDX;

   if (mode == MODE_DIAG)
   {
      /* iterate over all bounds in edgebound list */
      for (x = 0; x < N; x++)
      {
         d  = edg->bounds[x].diag;
         lb = edg->bounds[x].lb;
         rb = edg->bounds[x].rb;

         /* insert value across diag in range */
         for (k = lb; k < rb; k++)
         {
            i = k;
            j = d - i;

            if (i > Q || j > T ) 
               printf("OVER MAX: MAX=(Q=%d,T=%d) -> (%d,%d)\n", Q, T, i, j);

            MMX(i,j) = 1;
            IMX(i,j) = 1;
            DMX(i,j) = val;
            num_cells++;
         }
      }
   }
   else if (mode == MODE_ROW)
   {
      /* iterate over all bounds in edgebound list */
      for (x = 0; x < N; x++)
      {
         i  = edg->bounds[x].diag;
         lb = edg->bounds[x].lb;
         rb = edg->bounds[x].rb;

         /* insert value across row in range */
         for (j = lb; j < rb; j++)
         {
            IDX = (DEL_ST*(Q+1)*(T+1)) + ((i)*(T+1)) + (j);

            if (IDX > MAX)
               printf("OVER MAX: MAX=(%d,%d) -> (%d,%d)\n", Q, T, i, j);

            MMX(i,j) = 1;
            IMX(i,j) = 1;
            DMX(i,j) = val;
            num_cells++;
         }
      }
   }
   return num_cells;
}

/* test to show the cloud area, fill with value, return number of used cells  */
int cloud_Cell_Count(int Q, int T,
                   float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                   float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ])
{
   int i, j, num_cells;
   num_cells = 0;

   for (i = 0; i <= Q; i++)
   {
      for (j = 0; j <= T; j++)
      {
         if ( MMX(i,j) > 0 ) {
            num_cells++;
         }
      }
   }
   return num_cells;
}

/* print test and flush buffer */
void print_test()
{
   printf("test\n");
   fflush(stdout);
}