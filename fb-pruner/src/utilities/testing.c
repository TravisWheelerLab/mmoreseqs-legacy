/*******************************************************************************
 *  FILE:      cloud_search.h
 *  PURPOSE:   Testing for navigating through the matrices.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
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
#include "objects.h"

/* header */
#include "utilities.h"

/*
 *  FUNCTION:  TEST_set_color()
 *  SYNOPSIS:  Test output of all colors.
 */
void TEST_cycle_colors()
{
   int num_colors = 7;
   char* colors[7]   = { "default", "red", "green", "yellow", "blue", "magenta", "cyan" };
   /* cycle through all colors */
   for ( int i = 0; i < num_colors; i++ ) {
      for ( int j = 0; j < 2; j++ ) {
         TEST_set_color( colors[i], j );
         printf("Hello World!\n");
      }
   }
   /* return to default */
   TEST_set_color( colors[0], 0 );
   printf("Color test complete. Color reset to default.\n\n");
}

/*
 *  FUNCTION:  TEST_set_color()
 *  SYNOPSIS:  Set console text color using color name string <color> and boolean <bold>.
 */
void TEST_set_color( char*    color,
                     bool     bold )
{
   int   tbl_size    = 7;
   /* code that needs to printed to console to change color */
   char* pcodes[7]   = { "", ";31", ";32", ";33", ";34", ";35", ";36" };
   /* color names */
   char* colors[7]   = { "default", "red", "green", "yellow", "blue", "magenta", "cyan" };

   for ( int i = 0; i < tbl_size; i++ ) 
   {
      if ( strcmp( color, colors[i] ) == 0 ) 
      {
         TEST_set_color_num(i, bold);
         return;
      }
   }

   printf("ERROR: color '%s' is not supported.\n", color);
   return;
}

/*
 *  FUNCTION:  TEST_set_color_num()
 *  SYNOPSIS:  Set console text color by index number <color_id> and boolean <bold>.
 */
void TEST_set_color_num(  int    color_id,
                         bool     bold )
{
   int   tbl_size    = 7;
   /* code that needs to printed to console to change color */
   char* pcodes[7]   = { "", ";31", ";32", ";33", ";34", ";35", ";36" };
   /* color names */
   char* colors[7]   = { "default", "red", "green", "yellow", "blue", "magenta", "cyan" };

   if ( color_id < tbl_size ) {
      printf("\033[%d%sm", bold, pcodes[color_id] );
      return;
   }

   printf("ERROR: color id '%d' is not supported.\n", color_id);
   return;
}

/*
 *  FUNCTION:  TEST_fwd_cycle()
 *  SYNOPSIS:  Cycle through all indices in quadratic matrix in forward direction, antidiag-by-antidiag 
 */
void TEST_fwd_cycle( const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX,
                     MATRIX_2D*  sp_MX,
                     ALIGNMENT*  tr )
{
   /* vars for navigating matrix */
   int         d,i,j,k,x;                /* diagonal, row, column indices */
   int         lb, rb, le, re;           /* right/left bounds and edges */
   int         lb_new, rb_new;           /* tmp for new right/left bounds */
   int         num_cells;                /* number of cells in diagonal */
   int         d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int         dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int         dim_T, dim_Q, dim_TOT;    /* dimensions of submatrix being searched */
   int         total_cnt;                /* number of cells computed */
   TRACE*      beg;                      /* start point of alignment */
   TRACE*      end;                       /* end point of alignment */
   /* vars for recurrance */
   int         d_0, d_1, d_2;             /* for assigning prev array ptrs */
   int         d0,  d1,  d2;              /* for assigning prev array ptrs in mod for linear space */

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_3D*  test_MX;
   bool        overwrite;

   /* initialize debugging matrix */
   #if DEBUG
   {
      dbfp     = fopen( debugger->dbfp_path, "w+" );
      cloud_MX = debugger->cloud_MX;
      test_MX  = debugger->test_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   DP_MATRIX_Fill( Q, T, st_MX, sp_MX, 0 );

   /* get start and end points of viterbi alignment */
   beg = &(tr->traces->data[tr->beg]);
   end = &(tr->traces->data[tr->end]);

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (beg->i == 0 || beg->j == 0) {
      beg->i += 1;
      beg->j += 1;
   }

   /* TEST: overwrite alignment points */
   beg->i = 100;
   beg->j = 50;
   end->i = Q;
   end->j = T;

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   /* diag index of different start points, creating submatrix */
   d_st = beg->i + beg->j;
   // d_end = end->i + end->j;

   /* dimension of submatrix */
   dim_TOT = Q + T;
   dim_Q = Q - beg->i;
   dim_T = T - beg->j;

   /* diag index where num cells reaches highest point and begins diminishing */
   // dim_min = MIN(Q + beg->i, T + beg->j);
   // dim_max = MAX(Q + beg->i, T + beg->j);
   dim_min = MIN(d_st + dim_Q, d_st + dim_T);
   dim_max = MAX(d_st + dim_Q, d_st + dim_T);

   // set bounds using starting cell
   lb = beg->i;
   rb = beg->i + 1;
   num_cells = 0;

   /* iterate through diags */
   for (d = d_st; d <= d_end; d++)
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
      le = MAX( beg->i, d - T );
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
         /* cartesian coords */
         i = k;
         j = d - i;
         x = k;

         if ( MMX(i,j) != 0 || DMX(i,j) != 0 || IMX(i,j) != 0)
            printf("OVERWRITE AT (%d,%d)!\n", i, j);

         /*
          *   MAT: MX_M(i+1, j+1) => MX3_M(d_2, k+1)
          *   DEL: MX_M(i  , j+1) => MX3_M(d_1, k  )
          *   INS: MX_M(i+1, j  ) => MX3_M(d_1, k+1)
          */

         // MMX(i,j) = MMX(i+1,j+1) + 1;
         // IMX(i,j) = IMX(i+1,j) + 1;
         DMX(i,j) += d;

         #if DEBUG 
         {
            MX_2D( cloud_MX, i, j ) += 1.0;
         }
         #endif
      }
   }

   /* set each cell accessed to 1.0 */
   #if DEBUG
      // MX_2D( cloud_MX, beg->i, beg->j ) = -1.0;
      // MX_2D( cloud_MX, end->i, end->j ) = -1.0;
      DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
      fclose( dbfp );
   #endif
}

/*
 *  FUNCTION:  TEST_fwd_cycle3()
 *  SYNOPSIS:  Cycle through all indices in linear matrix, antidiag-by-antidiag 
 */
void TEST_fwd_cycle3(const int   Q, 
                     const int   T, 
                     MATRIX_3D*  st_MX, 
                     MATRIX_3D*  st_MX3,
                     MATRIX_2D*  sp_MX, 
                     ALIGNMENT*  tr )
{
   /* vars for navigating matrix */
   int         d,i,j,k,x;                /* diagonal, row, column indices */
   int         lb, rb, le, re;           /* right/left bounds and edges */
   int         lb_new, rb_new;           /* tmp for new right/left bounds */
   int         num_cells;                /* number of cells in diagonal */
   int         d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int         dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int         dim_T, dim_Q, dim_TOT;    /* dimensions of submatrix being searched */
   int         total_cnt;                /* number of cells computed */
   TRACE*      beg;                      /* start point of alignment */
   TRACE*      end;                       /* end point of alignment */
   /* vars for recurrance */
   int         d_0, d_1, d_2;             /* for assigning prev array ptrs */
   int         d0,  d1,  d2;              /* for assigning prev array ptrs in mod for linear space */

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_3D*  test_MX;
   bool        overwrite;

   /* initialize debugging matrix */
   #if DEBUG
   {
      dbfp     = fopen( debugger->dbfp_path, "w+" );
      cloud_MX = debugger->cloud_MX;
      test_MX  = debugger->test_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   DP_MATRIX_Fill( Q, T, st_MX, sp_MX, 0 );

   /* get start and end points of viterbi alignment */
   beg = &(tr->traces->data[tr->beg]);
   end = &(tr->traces->data[tr->end]);

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (beg->i == 0 || beg->j == 0) {
      beg->i += 1;
      beg->j += 1;
   }

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   /* diag index of different start points, creating submatrix */
   d_st = beg->i + beg->j;
   // d_end = end->i + end->j;

   /* dimension of submatrix */
   dim_TOT = Q + T;
   dim_Q = Q - beg->i;
   dim_T = T - beg->j;

   /* diag index where num cells reaches highest point and begins diminishing */
   // dim_min = MIN(Q + beg->i, T + beg->j);
   // dim_max = MAX(Q + beg->i, T + beg->j);
   dim_min = MIN(d_st + dim_Q, d_st + dim_T);
   dim_max = MAX(d_st + dim_Q, d_st + dim_T);

   // set bounds using starting cell
   lb = beg->i;
   rb = beg->i + 1;
   num_cells = 0;

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
      le = MAX(beg->i, d - T);
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
         x = k;

         if ( MMX3(d0,x) != 0 || DMX3(d0,x) != 0 || IMX3(d0,x) != 0) {
            printf("OVERWRITE AT %d:(%d,%d)=%f!\n", d, i, j, DMX(d0,x) );
            overwrite = true;
         }
 
         MMX3(d0,x) = MMX3(d2,k-1) + 1;
         IMX3(d0,x) = IMX3(d1,k-1) + 1;
         DMX3(d0,x) = DMX3(d1,k) + 1;
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

   /* set each cell accessed to 1.0 */
   #if DEBUG
      MX_2D( cloud_MX, beg->i, beg->j ) = 55.5;
      MX_2D( cloud_MX, end->i, end->j ) = 77.7;
      DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
      fclose( dbfp );
   #endif
}

/*
 *  FUNCTION:  TEST_bck_cycle()
 *  SYNOPSIS:  Cycle through all indices in quadratic matrix in backward direction, antidiag-by-antidiag 
 */
void TEST_bck_cycle( const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX,
                     MATRIX_2D*  sp_MX,
                     ALIGNMENT*  tr )
{
   /* vars for navigating matrix */
   int         d,i,j,k,x;                /* diagonal, row, column indices */
   int         lb, rb, le, re;           /* right/left bounds and edges */
   int         lb_new, rb_new;           /* tmp for new right/left bounds */
   int         num_cells;                /* number of cells in diagonal */
   int         d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int         dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int         dim_T, dim_Q, dim_TOT;    /* dimensions of submatrix being searched */
   int         total_cnt;                /* number of cells computed */
   TRACE*      beg;                      /* start point of alignment */
   TRACE*      end;                       /* end point of alignment */
   /* vars for recurrance */
   int         d_0, d_1, d_2;             /* for assigning prev array ptrs */
   int         d0,  d1,  d2;              /* for assigning prev array ptrs in mod for linear space */

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_3D*  test_MX;
   bool        overwrite;

   /* initialize debugging matrix */
   #if DEBUG
   {
      dbfp     = fopen( debugger->dbfp_path, "w+" );
      cloud_MX = debugger->cloud_MX;
      test_MX  = debugger->test_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   DP_MATRIX_Fill( Q, T, st_MX, sp_MX, 0 );

   /* get start and end points of viterbi alignment */
   beg = &(tr->traces->data[tr->beg]);
   end = &(tr->traces->data[tr->end]);

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (end->i == Q || end->j == T) {
      end->i -= 1;
      end->j -= 1;
   }

   /* diag index of different start points, creating submatrix */
   // d_st = beg->i + beg->j;
   d_end = end->i + end->j;

   /* dimension of submatrix */
   dim_TOT  = Q + T;
   dim_Q    = end->i;
   dim_T    = end->j;

   /* diag index at corners of dp matrix */
   d_st     = 0;
   d_end    = dim_TOT;
   d_cnt    = 0;

   /* diag index where num cells reaches highest point and begins diminishing */
   // dim_min = MIN(Q + beg->i, T + beg->j);
   // dim_max = MAX(Q + beg->i, T + beg->j);
   dim_min = MIN(d_st + dim_Q, d_st + dim_T);
   dim_max = MAX(d_st + dim_Q, d_st + dim_T);

   // set bounds using starting cell
   lb = end->i;
   rb = end->i + 1;
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
      le = MAX(end->i - (d_end - d), 0);
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
         x = k;

         if ( MMX(i,j) != 0 || DMX(i,j) != 0 || IMX(i,j) != 0)
            printf("OVERWRITE AT (%d,%d)!\n", i, j);

         /*
          *   MAT: MX_M(i+1, j+1) => MX3_M(d_2, k+1)
          *   DEL: MX_M(i  , j+1) => MX3_M(d_1, k  )
          *   INS: MX_M(i+1, j  ) => MX3_M(d_1, k+1)
          */

         MMX(i,j) = MMX(i+1,j+1) + 1;
         IMX(i,j) = IMX(i+1,j) + 1;
         DMX(i,j) += d;
      }
   }

   /* set each cell accessed to 1.0 */
   #if DEBUG
      MX_2D( cloud_MX, beg->i, beg->j ) = 55.5;
      MX_2D( cloud_MX, end->i, end->j ) = 77.7;
      DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
      fclose( dbfp );
   #endif
}

/*
 *  FUNCTION:  TEST_bck_cycle3()
 *  SYNOPSIS:  Cycle through all indices in linear matrix in backward direction, antidiag-by-antidiag 
 */
void TEST_bck_cycle3(const int   Q, 
                     const int   T, 
                     MATRIX_3D*  st_MX, 
                     MATRIX_3D*  st_MX3,
                     MATRIX_2D*  sp_MX, 
                     ALIGNMENT*  tr )
{
   /* vars for navigating matrix */
   int         d,i,j,k,x;                /* diagonal, row, column indices */
   int         lb, rb, le, re;           /* right/left bounds and edges */
   int         lb_new, rb_new;           /* tmp for new right/left bounds */
   int         num_cells;                /* number of cells in diagonal */
   int         d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int         dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int         dim_T, dim_Q, dim_TOT;    /* dimensions of submatrix being searched */
   int         total_cnt;                /* number of cells computed */
   TRACE*      beg;                      /* start point of alignment */
   TRACE*      end;                       /* end point of alignment */
   /* vars for recurrance */
   int         d_0, d_1, d_2;             /* for assigning prev array ptrs */
   int         d0,  d1,  d2;              /* for assigning prev array ptrs in mod for linear space */

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_3D*  test_MX;
   bool        overwrite;

   /* initialize debugging matrix */
   #if DEBUG
   {
      dbfp     = fopen( debugger->dbfp_path, "w+" );
      cloud_MX = debugger->cloud_MX;
      test_MX  = debugger->test_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   DP_MATRIX_Fill( Q, T, st_MX, sp_MX, 0 );

   /* get start and end points of viterbi alignment */
   beg = &(tr->traces->data[tr->beg]);
   end = &(tr->traces->data[tr->end]);

   /* We don't want to start on the edge and risk out-of-bounds (go to next match state) */
   if (end->i == Q || end->j == T) {
      end->i -= 1;
      end->j -= 1;
   }

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   /* diag index of different start points, creating submatrix */
   // d_st = beg->i + beg->j;
   d_end = end->i + end->j;

   /* dimension of submatrix */
   dim_TOT = Q + T;
   dim_Q = Q - beg->i;
   dim_T = T - beg->j;

   /* diag index where num cells reaches highest point and begins diminishing */
   // dim_min = MIN(Q + beg->i, T + beg->j);
   // dim_max = MAX(Q + beg->i, T + beg->j);
   dim_min = MIN(d_st + dim_Q, d_st + dim_T);
   dim_max = MAX(d_st + dim_Q, d_st + dim_T);

   // set bounds using starting cell
   lb = end->i;
   rb = end->i + 1;
   num_cells = 0;

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
      le = MAX(end->i - (d_end - d), 0);
      re = le + num_cells;

      lb = MAX(lb, le);
      rb = MIN(rb, re);

      /* MAIN RECURSION */
      /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
      for (k = lb; k < rb; k++)
      {
         i = k;
         j = d - i;
         x = k;

         if ( MMX3(d0,x) != 0 || DMX3(d0,x) != 0 || IMX3(d0,x) != 0) {
            printf("OVERWRITE AT %d:(%d,%d)=%f!\n", d, i, j, DMX(d0,x) );
            overwrite = true;
         }

         /*
          *   MAT: MX_M(i+1, j+1) => MX3_M(d_2, k+1)
          *   DEL: MX_M(i  , j+1) => MX3_M(d_1, k  )
          *   INS: MX_M(i+1, j  ) => MX3_M(d_1, k+1)
          */

         MMX3(d0,x) = MMX3(d2,k+1) + 1;
         IMX3(d0,x) = IMX3(d1,k+1) + 1;
         DMX3(d0,x) += d;
      }

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

   /* set each cell accessed to 1.0 */
   #if DEBUG
      MX_2D( cloud_MX, beg->i, beg->j ) = 55.5;
      MX_2D( cloud_MX, end->i, end->j ) = 77.7;
      DP_MATRIX_VIZ_Dump( cloud_MX, stdout );
      fclose( dbfp );
   #endif
} 

/*
 *  FUNCTION:  MATRIX_2D_Cloud_Fill()
 *  SYNOPSIS:  Increment MATRIX_2D with value according to EDGEBOUNDS, returns number of cells covered by EDGEBOUNDS 
 */
int MATRIX_2D_Cloud_Fill(  MATRIX_2D*     cloud_MX,
                           EDGEBOUNDS*    edg,
                           float          val )
{
   int x, d, i, j, k;
   int rb, lb;
   int num_cells = 0;
   int N = edg->N;
   int mode = edg->edg_mode;

   if (mode == EDG_DIAG)
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

            MX_2D( cloud_MX, i, j ) += val;
            num_cells++;
         }
      }
   }
   else if (mode == EDG_ROW)
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
            MX_2D( cloud_MX, i, j ) += val;
            num_cells++;
         }
      }
   }
   return num_cells;
}

/*
 *  FUNCTION:  MATRIX_2D_Cloud_Compare()
 *  SYNOPSIS:  Compare <cloud_MX_a> and <cloud_MX_b> ( equality is if cells are both positive, negative, or zero ).
 *             When debugging, stores heatmap of differences into <debugger->cloud_MX>.
 *  RETURN:    If equal, returns 0.
 *             Otherwize, return 1. (or number of inequal cells in DEBUG mode). 
 */
int MATRIX_2D_Cloud_Compare(  MATRIX_2D*     cloud_MX_a,
                              MATRIX_2D*     cloud_MX_b )
{
   MATRIX_2D* cloud_MX_diff = debugger->cloud_MX;
   #if DEBUG
   {
      MATRIX_2D_Fill( cloud_MX_diff, 0 );
   }
   #endif

   int cmp = 0;
   for ( int i = 0; i < cloud_MX_a->R; i++ ) 
   {
      for ( int j = 0; j < cloud_MX_b->C; j++ ) 
      {
         float mx_a  = MX_2D( cloud_MX_a, i, j );
         float mx_b  = MX_2D( cloud_MX_b, i, j );
         bool  zero  = ( mx_a == 0 && mx_b == 0 );    /* both cells are zero */
         bool  pos   = ( mx_a > 0 && mx_b > 0 );      /* both cells are positive values */
         bool  neg   = ( mx_a < 0 && mx_b < 0 );      /* both cells are negative values */

         /* if both numbers are neither both zero, both positive, or both negative */
         if ( !( zero || pos || neg ) ) {
            #if DEBUG 
            {
               MX_2D( cloud_MX_diff, i, j ) = -1.0;
               if ( cmp == 0 )  /* only reports first inequality */
                  fprintf(stdout, "CLOUD INEQUALITY found at (%d,%d) => MX_A: %f vs MX_B: %f\n", i, j, mx_a, mx_b);
            }
            #endif
            cmp += 1;
            #if !DEBUG 
            {
               return cmp;
            }
            #endif
         }
      }
   }
   return cmp;
}

/*
 *  FUNCTION:  MATRIX_2D_Cloud_Count()
 *  SYNOPSIS:  Count number of cells in MATRIX_2D with positive values. 
 */
int MATRIX_2D_Cloud_Count(  MATRIX_2D*  cloud_MX )
{
   int i, j, num_cells;
   num_cells = 0;

   for (i = 0; i < cloud_MX->R; i++) {
      for (j = 0; j < cloud_MX->C; j++) {
         if ( MX_2D(cloud_MX, i, j) > 0 ) {
            num_cells++;
         }
      }
   }
   return num_cells;
}

/*
 *  FUNCTION: EDGEBOUNDS_Compare_by_Cloud()
 *  SYNOPSIS: Compare two EDGEBOUNDS by filling cloud matrices.
 */
int EDGEBOUNDS_Compare_by_Cloud( EDGEBOUNDS*    edg_a,
                                 MATRIX_2D*     mx_a,
                                 EDGEBOUNDS*    edg_b,
                                 MATRIX_2D*     mx_b )
{
   MATRIX_2D_Fill( mx_a, 0 );
   MATRIX_2D_Cloud_Fill( mx_a, edg_a, 1 );
   MATRIX_2D_Fill( mx_b, 0 );
   MATRIX_2D_Cloud_Fill( mx_b, edg_b, 1 );
   return MATRIX_2D_Cloud_Compare( mx_a, mx_b );
}


/*
 *  FUNCTION:  EDGEBOUNDS_Compare_by_Cloud_Single()
 *  SYNOPSIS:  Compares two EDGEBOUNDS by filling single cloud matrix.
 *             If equal, returns 0.  Else number of inequal cells.
 */
int EDGEBOUNDS_Compare_by_Cloud_Single(   MATRIX_2D*     mx,
                                          EDGEBOUNDS*    edg_a,
                                          EDGEBOUNDS*    edg_b )
{
   MATRIX_2D_Fill( mx, 0 );
   MATRIX_2D_Cloud_Fill( mx, edg_a, 1 );
   MATRIX_2D_Cloud_Fill( mx, edg_b, -1 );
   return MATRIX_2D_Check_Value( mx, 0 );
}

