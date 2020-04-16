/*******************************************************************************
 *  @file cloud_search.h
 *  @brief Testing for navigating through the matrices.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* objects */
#include "objects/structs.h"
#include "objects/edgebound.h"
#include "objects/hmm_profile.h"
#include "objects/matrix/matrix_2d.h"
#include "objects/matrix/matrix_3d.h"

/* local imports */
#include "utilities/utility.h"
#include "hmm_parser.h"

/* header */
#include "testing.h"

/* cycle through all indices in quadratic matrix, diag-by-diag */
void fwd_test_cycle( const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX,
                     MATRIX_2D*  sp_MX,
                     ALIGNMENT*  tr )
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
   TRACE beg, end;            /* start and end point of alignment */

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
   // dim_min = MIN(Q + beg.i, T + beg.j);
   // dim_max = MAX(Q + beg.i, T + beg.j);
   dim_min = MIN(d_st + dim_Q, d_st + dim_T);
   dim_max = MAX(d_st + dim_Q, d_st + dim_T);

   // set bounds using starting cell
   lb = beg.i;
   rb = beg.i + 1;
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
      le = MAX(beg.i, d - T);
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

         if ( MMX_M(i,j) != 0 || DMX_M(i,j) != 0 || IMX_M(i,j) != 0)
            printf("OVERWRITE AT (%d,%d)!\n", i, j);

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         /* NOTE: (i-1,j-1) => (d-2,k-1) */ 
         MMX_M(i,j) = MMX_M(i-1,j-1) + 1;

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i-1,j) => (d-1,k-1) */
         IMX_M(i,j) = IMX_M(i-1,j) + 1;

         /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i,j-1) => (d-1, k) */
         DMX_M(i,j) = DMX_M(i,j-1) + 1;
      }
   }
}

/* cycle through all indices in linear matrix, diag-by-diag */
void fwd_test_cycle3(const int   Q, 
                     const int   T, 
                     MATRIX_3D*  st_MX, 
                     MATRIX_3D*  st_MX3,
                     MATRIX_2D*  sp_MX, 
                     ALIGNMENT*  tr )
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
   TRACE beg, end;            /* start and end point of alignment */
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
      le = MAX(beg.i, d - T);
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

         if ( MMX3_M(d0,x) != 0 || DMX3_M(d0,x) != 0 || IMX3_M(d0,x) != 0) {
            printf("OVERWRITE AT %d:(%d,%d)=%f!\n", d, i, j, DMX_M(d0,x) );
            overwrite = true;
         }

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         /* NOTE: (i-1,j-1) => (d-2,k-1) */ 
         MMX3_M(d0,x) = MMX3_M(d2,k-1) + 1;

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i-1,j) => (d-1,k-1) */
         IMX3_M(d0,x) = IMX3_M(d1,k-1) + 1;

         /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i,j-1) => (d-1, k) */
         DMX3_M(d0,x) = DMX3_M(d1,k) + 1;

         if (overwrite) {
            overwrite = false;
            MMX3_M(d0,x) = 30;
         }
      }

      /* NAIVE SCRUB */
      for (k = 0; k < ((T+1)+(Q+1)); k++) {
         MMX3_M(d2,k) = IMX3_M(d2,k) = DMX3_M(d2,k) = 0;
      }

      /* Embed Current Row into Quadratic Array - FOR DEBUGGING */
      for (k = le; k < re; k++) {
         /* quadratic coords */
         i = k;
         j = d_0 - i;
         /* linear coords */
         x = k;

         MMX_M(i,j) = MMX3_M(d0,x);
         IMX_M(i,j) = IMX3_M(d0,x);
         DMX_M(i,j) = DMX3_M(d0,x);
      }
   }
}

/* cycle through all indices in quadratic matrix in reverse, diag-by-diag */
void bck_test_cycle( const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX,
                     MATRIX_2D*  sp_MX,
                     ALIGNMENT*  tr )
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
   TRACE beg, end;            /* start and end point of alignment */

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
   beg = tr->traces[tr->beg];
   end = tr->traces[tr->end];

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (beg.i == Q+1 || beg.j == T+1) {
      beg.i -= 1;
      beg.j -= 1;
   }

   // diag index at corners of dp matrix
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   // diag index of different start points, creating submatrix
   // d_st = beg.i + beg.j;
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

         if ( MMX_M(i,j) != 0 || DMX_M(i,j) != 0 || IMX_M(i,j) != 0)
            printf("OVERWRITE AT (%d,%d)!\n", i, j);

         /*
          *   MAT: MX_M(i+1, j+1) => MX3_M(d_2, k+1)
          *   DEL: MX_M(i  , j+1) => MX3_M(d_1, k  )
          *   INS: MX_M(i+1, j  ) => MX3_M(d_1, k+1)
          */

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         /* NOTE: (i+1,j+1) => (d+2, k+1) */ 
         MMX_M(i,j) = MMX_M(i+1,j+1) + 1;

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i+1,j) => (d+1, k  ) */
         IMX_M(i,j) = IMX_M(i+1,j) + 1;

         /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i,j+1) => (d+1, k+1) */
         // DMX_M(i,j) = DMX_M(i,j+1) + 1;
         DMX_M(i,j) += d;
      }
   }
}

/* cycle through all indices in linear matrix in reverse, diag-by-diag */
void bck_test_cycle3(const int   Q, 
                     const int   T, 
                     MATRIX_3D*  st_MX, 
                     MATRIX_3D*  st_MX3,
                     MATRIX_2D*  sp_MX, 
                     ALIGNMENT*  tr )
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
   TRACE beg, end;            /* start and end point of alignment */
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

         if ( MMX3_M(d0,x) != 0 || DMX3_M(d0,x) != 0 || IMX3_M(d0,x) != 0) {
            printf("OVERWRITE AT %d:(%d,%d)=%f!\n", d, i, j, DMX_M(d0,x) );
            overwrite = true;
         }

         /*
          *   MAT: MX_M(i+1, j+1) => MX3_M(d_2, k+1)
          *   DEL: MX_M(i  , j+1) => MX3_M(d_1, k  )
          *   INS: MX_M(i+1, j  ) => MX3_M(d_1, k+1)
          */

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         /* NOTE: (i+1,j+1) => (d+2,k+1) */ 
         MMX3_M(d0,x) = MMX3_M(d2,k+1) + 1;

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i+1,j) => (d+1,k+1) */
         IMX3_M(d0,x) = IMX3_M(d1,k+1) + 1;

         /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         /* NOTE: (i,j+1) => (d+1, k) */
         // DMX3_M(d0,x) = DMX3_M(d1,k) + 1;
         DMX3_M(d0,x) += d;

         if (overwrite) {
            overwrite = false;
            // MMX3_M(d0,x) = 30;
         }
      }

      // /* SCRUB 2-BACK ANTIDIAGONAL */
      // for (k = lb_2; k < rb_2; k++) {
      //    MMX3_M(d2,k) = IMX3_M(d2,k) = DMX3_M(d2,k) = -INF;
      // }
      // rb_2 = rb_1;
      // lb_2 = lb_1;
      // rb_1 = rb;
      // lb_1 = lb;

      /* NAIVE SCRUB - FOR DEBUGGING */
      for (k = 0; k < (T+1)+(Q+1); k++) {
         MMX3_M(d2,k) = IMX3_M(d2,k) = DMX3_M(d2,k) = 0;
      }

      /* Embed Current Row into Quadratic Array - FOR DEBUGGING */
      for (k = le; k <= re; k++) {
         i = k;
         j = d_0 - i;

         MMX_M(i,j) = MMX3_M(d0,k);
         IMX_M(i,j) = IMX3_M(d0,k);
         DMX_M(i,j) = DMX3_M(d0,k);
      }
   }
}  



/* test to show the cloud area, fill with value, return number of used cells  */
int cloud_Fill(const int   Q, 
               const int   T,
               MATRIX_3D*  st_MX,
               MATRIX_2D*  sp_MX,
               EDGEBOUNDS* edg,
               float       val, 
               int         mode )
{
   int x, d, i, j, k;
   int rb, lb;
   int num_cells = 0;
   int N = edg->N;

   if (mode == MODE_DIAG)
   {
      /* iterate over all bounds in edgebound list */
      for (x = 0; x < N; x++)
      {
         d  = edg->bounds[x].id;
         lb = edg->bounds[x].lb;
         rb = edg->bounds[x].rb;

         /* insert value across diag in range */
         for (k = lb; k < rb; k++)
         {
            i = k;
            j = d - i;
            // printf("(i,j)=(%d,%d)\n", i, j);

            MMX_M(i,j) += 1;
            // if (MMX_M(i,j) > 1)
            //    printf("DOUBLE HIT!  MAX=(Q=%d,T=%d) -> d:(%d,%d) = r:(%d,%d)\n", Q, T, d, k, i, j);

            IMX_M(i,j) = 1;
            DMX_M(i,j) += val;
            num_cells++;
         }
      }
   }
   else if (mode == MODE_ROW)
   {
      /* iterate over all bounds in edgebound list */
      for (x = 0; x < N; x++)
      {
         i  = edg->bounds[x].id;
         lb = edg->bounds[x].lb;
         rb = edg->bounds[x].rb;

         /* insert value across row in range */
         for (j = lb; j < rb; j++)
         {
            MMX_M(i,j) += 1;
            IMX_M(i,j) =  1;
            DMX_M(i,j) += val;
            num_cells++;
         }
      }
   }
   return num_cells;
}


/* test to show the cloud area, fill with value, return number of used cells  */
int cloud_Solid_Fill(const int      Q, 
                     const int      T,
                     MATRIX_3D*     st_MX,
                     MATRIX_2D*     sp_MX,
                     EDGEBOUNDS*    edg,
                     float          val, 
                     int            mode )
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
         d  = edg->bounds[x].id;
         lb = edg->bounds[x].lb;
         rb = edg->bounds[x].rb;

         /* insert value across diag in range */
         for (k = lb; k < rb; k++)
         {
            i = k;
            j = d - i;

            if (i > Q || j > T ) 
               printf("OVER MAX: MAX=(Q=%d,T=%d) -> (%d,%d)\n", Q, T, i, j);

            MMX_M(i,j) = 1;
            IMX_M(i,j) = 1;
            DMX_M(i,j) = val;
            num_cells++;
         }
      }
   }
   else if (mode == MODE_ROW)
   {
      /* iterate over all bounds in edgebound list */
      for (x = 0; x < N; x++)
      {
         i  = edg->bounds[x].id;
         lb = edg->bounds[x].lb;
         rb = edg->bounds[x].rb;

         /* insert value across row in range */
         for (j = lb; j < rb; j++)
         {
            IDX = (DEL_ST*(Q+1)*(T+1)) + ((i)*(T+1)) + (j);

            if (IDX > MAX)
               printf("OVER MAX: MAX=(%d,%d) -> (%d,%d)\n", Q, T, i, j);

            MMX_M(i,j) = 1;
            IMX_M(i,j) = 1;
            DMX_M(i,j) = val;
            num_cells++;
         }
      }
   }
   return num_cells;
}

/* test to show the cloud area, fill with value, return number of used cells  */
int cloud_Cell_Count(const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX,
                     MATRIX_2D*  sp_MX )
{
   int i, j, num_cells;
   num_cells = 0;

   for (i = 0; i <= Q; i++)
   {
      for (j = 0; j <= T; j++)
      {
         if ( MMX_M(i,j) > 0 ) {
            num_cells++;
         }
      }
   }
   return num_cells;
}

/*
 *  FUNCTION:  dp_matrix_Print()
 *  SYNOPSIS:  Print out dynamic programming matrix to screen.
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix
 *
 *  RETURN:
 */
void dp_matrix_Print(const int         Q, 
                     const int         T,
                     MATRIX_3D*        st_MX,
                     MATRIX_2D*        sp_MX )
{
   /* PRINT resulting dp matrix */
   printf("\n\n==== DP MATRIX BEGIN ==== \n");
   /* Header */
   printf("\t");
   for (int i = 0; i <= T; i++)
   {
      printf("%d\t", i);
   }
   printf("\n");

   /* Row-by-Row */
   for (int i = 0; i < Q + 1; i++)
   {
      printf( "%d M\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%.3f\t", MMX_M(i, j) );
      }
      printf("\n");

      printf( "%d I\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%.3f\t", IMX_M(i, j) );
      }
      printf("\n");

      printf( "%d D\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%.3f\t", DMX_M(i, j) );
      }
      printf("\n\n");
   }

   printf("=== SPECIAL STATES ===\n");
   printf("N:\t");
   for (int i = 0; i <= Q; i++)
   { printf( "%.3f\t", XMX_M(SP_N, i) ); }
   printf("\n");
   printf("J:\t");
   for (int i = 0; i <= Q; i++)
   { printf( "%.3f\t", XMX_M(SP_J, i) ); }
   printf("\n");
   printf("E:\t");
   for (int i = 0; i <= Q; i++)
   { printf( "%.3f\t", XMX_M(SP_E, i) ); }
   printf("\n");
   printf("C:\t");
   for (int i = 0; i <= Q; i++)
   { printf( "%.3f\t", XMX_M(SP_C, i) ); }
   printf("\n");
   printf("B:\t");
   for (int i = 0; i <= Q; i++)
   { printf( "%.3f\t", XMX_M(SP_B, i) ); }
   printf("\n");

   printf("==== DP MATRIX END ==== \n\n");
}


/*
 *  FUNCTION:  dp_matrix_Print()
 *  SYNOPSIS:  Print out dynamic programming matrix to screen.
 *
 *  PURPOSE:
 *
 *  ARGS:      <T>         target length,
 *             <st_MX3>    Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix
 *
 *  RETURN:
 */
void dp_matrix_Print3(  const int         T, 
                        const int         Q,
                        const MATRIX_3D*  st_MX3 )
{
   /* PRINT resulting dp matrix */
   printf("\n\n==== DP MATRIX BEGIN ==== \n");
   /* Header */
   printf("\t");
   for (int i = 0; i <= T+1; i++)
   {
      printf("%d\t", i);
   }
   printf("\n");

   /* Row-by-Row */
   for (int i = 0; i < 3; i++)
   {
      printf( "%d M\t", i );
      for (int j = 0; j <= T+1; j++)
      {
         printf( "%.3f\t", MMX3_M(i, j) );
      }
      printf("\n");

      printf( "%d I\t", i );
      for (int j = 0; j <= T+1; j++)
      {
         printf( "%.3f\t", IMX3_M(i, j) );
      }
      printf("\n");

      printf( "%d D\t", i );
      for (int j = 0; j <= T+1; j++)
      {
         printf( "%.3f\t", DMX3_M(i, j) );
      }
      printf("\n\n");
   }

   printf("==== DP MATRIX END ==== \n\n");
}


/*
 *  FUNCTION:  test_matrix_Print()
 *  SYNOPSIS:  Print out dynamic programming matrix to screen.
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <test_MX>   matrix
 *
 *  RETURN:
 */
void test_matrix_Print( const int         Q, 
                        const int         T,
                        const MATRIX_3D*  test_MX )
{
   /* PRINT resulting dp matrix */
   printf("\n\n==== TEST MATRIX BEGIN ====\n");
   /* Header */
   printf("\t");
   for (int i = 0; i <= T; i++)
   {
      printf("%d\t", i);
   }
   printf("\n");

   /* Row-by-Row */
   for (int i = 0; i < Q + 1; i++)
   {
      printf( "%d\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%3.0f\t", TMX_M(i, j) );
      }
      printf("\n");
   }
   printf("==== TEST MATRIX END ====\n\n");
}

/* Clear all matrix values to -INF. (for testing) */
void dp_matrix_Clear(const int   Q, 
                     const int   T, 
                     MATRIX_3D*  st_MX, 
                     MATRIX_2D*  sp_MX )
{
   for (int i = 0; i <= Q; i++)
   {
      for (int j = 0; j <= T; j++) {
         MMX_M(i, j) = IMX_M(i, j) = DMX_M(i, j) = -INF;
      }

      for (int j = 0; j < NUM_SPECIAL_STATES; j++) {
         XMX_M(j, i) = -INF;
      }
   }
}

/* Clear all matrix values to -INF. (for testing) */
void dp_matrix_Clear3(  const int   Q, 
                        const int   T,
                        MATRIX_3D*  st_MX3,
                        MATRIX_2D*  sp_MX )
{
   for (int i = 0; i <= Q; i++)
   {
      for (int j = 0; j < NUM_SPECIAL_STATES; j++) {
         XMX_M(j, i) = -INF;
      }
   }

   for (int i = 0; i < 3; i++) {
      for (int j = 0; j < (T+1)+(Q+1); j++) {
         MMX3_M(i, j) = IMX3_M(i, j) = DMX3_M(i, j) = -INF;
      }
   }
}

/* Set all matrix values to val */
void dp_matrix_Clear_X( const int   Q, 
                        const int   T, 
                        MATRIX_3D*  st_MX, 
                        MATRIX_2D*  sp_MX,
                        float       val )
{
   for (int i = 0; i <= Q; i++)
   {
      for (int j = 0; j <= T; j++) {
         MMX_M(i, j) = IMX_M(i, j) = DMX_M(i, j) = val;
      }

      for (int j = 0; j < NUM_SPECIAL_STATES; j++) {
         XMX_M(j, i) = val;
      }
   }
}


/* Clear all matrix values to -INF. (for testing) */
void dp_matrix_Clear_X3(const int   Q, 
                        const int   T,
                        MATRIX_3D*  st_MX3,
                        MATRIX_2D*  sp_MX,
                        int         val)
{
   for (int i = 0; i <= Q; i++)
   {
      for (int j = 0; j < NUM_SPECIAL_STATES; j++) {
         XMX_M(j, i) = val;
      }
   }

   for (int i = 0; i < 3; i++) {
      for (int j = 0; j < (T+1)+(Q+1); j++) {
         MMX3_M(i, j) = IMX3_M(i, j) = DMX3_M(i, j) = val;
      }
   }
}


/* Set all matrix values to val */
int dp_matrix_Compare ( const int   Q, 
                        const int   T,
                        MATRIX_3D*  st_MX_1,
                        MATRIX_2D*  sp_MX_1,
                        MATRIX_3D*  st_MX_2,
                        MATRIX_2D*  sp_MX_2 )
{
   int i, j, st;

   for (i = 0; i <= Q; i++)
   {
      for (j = 0; j <= T; j++) 
      {
         for (st = 0; st < NUM_NORMAL_STATES; st++) 
         {
            if ( ST_MX_M(st_MX_1, st, i, j) != ST_MX_M(st_MX_2, st, i, j) ) 
            {
               float val1 = ST_MX_M(st_MX_1, st, i, j);
               float val2 = ST_MX_M(st_MX_2, st, i, j);
               printf("MATRIX NOT EQUAL at (%d,%d): %f vs %f\n", i, j, val1, val2);
               return false;
            } 
         }
      }

      for (st = 0; st < NUM_SPECIAL_STATES; st++) 
      {
         if ( SP_MX_M(sp_MX_1, st, i) != SP_MX_M(sp_MX_2, st, i) ) 
         {
            printf("MATRIX NOT EQUAL IN SPECIAL at (%d)\n", i);
            return false;
         }
      }
   }
   return true;
}

/* Copy source matrix into destination matrix */
void dp_matrix_Copy (const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX_src,
                     MATRIX_2D*  sp_MX_src,
                     MATRIX_3D*  st_MX_dst,
                     MATRIX_2D*  sp_MX_dst )
{
   int i, j, st;

   for (i = 0; i <= Q; i++)
   {
      for (j = 0; j <= T; j++) 
      {
         for (st = 0; st < NUM_NORMAL_STATES; st++) 
         {
            ST_MX_M(st_MX_dst, st, i, j) = ST_MX_M(st_MX_src, st, i, j);
         }
      }

      for (st = 0; st < NUM_SPECIAL_STATES; st++) 
      {
         SP_MX_M(sp_MX_dst, st, i) = SP_MX_M(sp_MX_src, st, i);
      }
   }
}


/*
 *  FUNCTION:  dp_matrix_Save()
 *  SYNOPSIS:  Save dynamic programming matrix to file.
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State
 *             <f>         Filename
 *
 *  RETURN:
 */
void dp_matrix_Save( const int         Q, 
                     const int         T, 
                     const MATRIX_3D*  st_MX, 
                     const MATRIX_2D*  sp_MX,
                     const char*       _filename_ )
{
   printf("Saving matrix...\n");
   FILE *fp;
   fp = fopen(_filename_, "w");
   dp_matrix_Dump(Q, T, st_MX, sp_MX, fp);
   fclose(fp);
   printf("Saved matrix to: '%s'\n", _filename_);
}

/*
 *  FUNCTION:  dp_matrix_Dump()
 *  SYNOPSIS:  Save dynamic programming matrix to file.
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State
 *             <fp>        File Pointer
 *
 *  RETURN:
 */
void dp_matrix_Dump( const int         Q, 
                     const int         T, 
                     const MATRIX_3D*  st_MX, 
                     const MATRIX_2D*  sp_MX,
                     FILE*             fp )
{
   /* PRINT resulting dp matrix */
   fprintf(fp, "##### DP MATRIX ##### \n");
   fprintf(fp, "XDIM\t%d\t%d\n\n", Q, T);

   /* Header */
   fprintf(fp, "##### NORMAL STATES #####\n");
   fprintf(fp, "XMATRIX\n");
   /* Header Indices */
   fprintf(fp, "#\t");
   for (int i = 0; i <= T; i++)
   {
      fprintf(fp, "%d\t", i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row Values */
   for (int i = 0; i < Q + 1; i++)
   {
      fprintf(fp, "%d\tM\t%d", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%9.4f ", MMX_M(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "\tI\t" );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%9.4f ", IMX_M(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "\tD\t" );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%9.4f ", DMX_M(i, j) );
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "/\n\n");

   fprintf(fp, "###### SPECIAL STATES #####\n");
   fprintf(fp, "N\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%9.4f ", XMX_M(SP_N, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "J\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%9.4f ", XMX_M(SP_J, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "E\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%9.4f ", XMX_M(SP_E, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "C\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%9.4f ", XMX_M(SP_C, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "B\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%9.4f ", XMX_M(SP_B, i) ); }
   fprintf(fp, "\n");
}

/*
 *  FUNCTION:  dp_matrix_trace_Save()
 *  SYNOPSIS:  Save dynamic programming matrix to file.
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State
 *             <tr>        TRACE object
 *             <f>         Filename
 *
 *  RETURN:
 */
void dp_matrix_trace_Save( const int         Q, 
                           const int         T, 
                           const MATRIX_3D*  st_MX, 
                           const MATRIX_2D*  sp_MX,
                           const ALIGNMENT*  tr,
                           const char*       _filename_ )
{
   printf("Saving matrix...\n");
   FILE *fp;
   fp = fopen(_filename_, "w");
   dp_matrix_trace_Dump( Q, T, st_MX, sp_MX, tr, fp );
   fclose( fp );
   printf("Saved Matrix to: '%s'\n", _filename_);
}


/*
 *  FUNCTION:  dp_matrix_trace_Save()
 *  SYNOPSIS:  Save dynamic programming matrix to file.
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State
 *             <tr>        TRACE object
 *             <f>         Filename
 *
 *  RETURN:
 */
void dp_matrix_trace_Dump( const int         Q, 
                           const int         T, 
                           const MATRIX_3D*  st_MX, 
                           const MATRIX_2D*  sp_MX,
                           const ALIGNMENT*  tr,
                           const FILE*       fp )
{
   /* PRINT resulting dp matrix */
   fprintf(fp, "##### DP MATRIX ##### \n");
   fprintf(fp, "XDIM\t%d\t%d\n\n", Q, T);

   /* Traceback */
   fprintf(fp, "XTRACE\n");
   /* Traceback */
   for (int i = 0; i < tr->N; i++) 
   {
      int st = tr->traces[i].st;
      if ( st == M_ST ) {
         fprintf(fp, "[%d]%s\t%d\t%d\n", i, STATE_NAMES[st], tr->traces[i].i, tr->traces[i].j);
      }
   }
   fprintf(fp, "/\n\n");

   /* Header */
   fprintf(fp, "##### NORMAL STATES #####\n");
   fprintf(fp, "XMATRIX\n");
   /* Header Indices */
   fprintf(fp, "#\t");
   for (int i = 0; i <= T; i++)
   {
      fprintf(fp, "%d\t", i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row Values */
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "M %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", MMX_M(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "I %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", IMX_M(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "D %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", DMX_M(i, j) );
      }
      fprintf(fp, "\n\n");
   }
   fprintf(fp, "/\n\n");

   fprintf(fp, "###### SPECIAL STATES #####\n");
   fprintf(fp, "XSPECIAL\n");
   fprintf(fp, "N\t");
   for (int i = 0; i <= Q; i++)
   { 
      fprintf(fp, "%.3f\t", XMX_M(SP_N, i) ); 
   }
   fprintf(fp, "\n");
   fprintf(fp, "J\t");
   for (int i = 0; i <= Q; i++)
   { 
      fprintf(fp, "%.3f\t", XMX_M(SP_J, i) ); 
   }
   fprintf(fp, "\n");
   fprintf(fp, "E\t");
   for (int i = 0; i <= Q; i++)
   { 
      fprintf(fp, "%.3f\t", XMX_M(SP_E, i) ); 
   }
   fprintf(fp, "\n");
   fprintf(fp, "C\t");
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "%.3f\t", XMX_M(SP_C, i) ); 
   }
   fprintf(fp, "\n");
   fprintf(fp, "B\t");
   for (int i = 0; i <= Q; i++)
   { 
      fprintf(fp, "%.3f\t", XMX_M(SP_B, i) ); 
   }
   fprintf(fp, "\n/\n");
}