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

// macro functions
// NOTE: wrap all macro vars in parens!!
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

/* test to cycle through all diags */
void fwd_test_cycle(int Q, int T,
                   float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                   float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                   TRACEBACK* tr)
{
   int d, i, j, k;               /* diag, row, col indices */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int num_cells;                /* number of cells in current diagonal */
   int d_cnt, total_cnt;         /* count of cells in current diag, count of total cells */
   int d_st, d_end;              /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and begins diminishing */
   int dim_T, dim_Q, dim_TOT;    /* dimensions of submatrix being searched */
   COORDS start, end;            /* start and end point of alignment */

   total_cnt = 0;

   /* INIT ANTI-DIAGS */
   /* test coords */
   start = tr->first_m;
   end = tr->last_m;

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;

   /* diag index of different start points, creating submatrix */
   d_st = start.i + start.j;
   // d_end = end.i + end.j;

   /* dimension of submatrix */
   dim_TOT = Q + T;
   dim_Q = Q - start.i;
   dim_T = T - start.j;

   /* diag index where num cells reaches highest point and begins diminishing */
   // dim_min = min(Q + start.i, T + start.j);
   // dim_max = max(Q + start.i, T + start.j);
   dim_min = min(d_st + dim_Q, d_st + dim_T);
   dim_max = max(d_st + dim_Q, d_st + dim_T);

   printf("TRACEBACK START: (%d, %d)\n", start.i, start.j);
   printf("MATRIX DIMS: (Q=%d, T=%d)\n", Q, T);
   printf("SUBMATRIX DIMS: (Q=%d, T=%d)\n", dim_Q, dim_T);

   // set bounds using starting cell
   lb = start.i;
   rb = start.i + 1;
   num_cells = 0;

   // printf("test cycle. (Q=%d,T=%d) (dim_Q=%d, dim_T=%d)...\n", Q, T, dim_Q, dim_T);

   /* TESTING VARS */
   bool growing = true;
   bool not_shrinking = true;

   /* iterate through diags */
   for (d = d_st; d <= d_end; d++, rb++)
   {
      // printf("d=%d, dim_min=%d, dim_max=%d, num_cells=%d \n", d, dim_min, dim_max, num_cells);

      d_cnt = 0;

      /* is dp matrix diagonal growing or shrinking? */
      if (d <= dim_min)
         num_cells++;
      if (d > dim_max)
         num_cells--;

      // printf("d: %d, num_cells: %d\n", d, num_cells);

      if (d > dim_min && growing) {
         growing = false;
         printf("END growing => dim_min: %d, d: %d, num_cells: %d\n", dim_min, d, num_cells);
         fflush(stdout);
      }
      if (d > dim_max && not_shrinking) {
         not_shrinking = false;
         printf("START shrinking => dim_max: %d, d: %d, num_cells: %d\n", dim_max, d, num_cells);
         fflush(stdout);
      }

      // printf("[%d,%d]\n", lb, rb);

      /* find diag cells that are inside matrix bounds */
      le = max(start.i, d - T);
      re = le + num_cells;

      // printf("[%d,%d]\n", le, re);

      /* for testing, bounds = edges */
      if (lb < le)
         lb = le;
      if (rb > re)
         rb = re;

      // printf("%d: -> edg[%d, %d] -> bounds[%d, %d]\n", d, le, re, lb, rb );

      // printf("[%d,%d]\n", lb, rb);

      /* iterate through cells of diag */
      for (k = lb; k < rb; k++)
      {
         i = k;
         j = d - i;
         // printf("(%d,%d)->", i, j);

         MMX(i,j) += 2;
         // MMX(i, j) = i;
         // IMX(i, j) = j;
         // DMX(i, j) = total_cnt;

         if ( (i - j) % 5 == 0 ) {
            MMX(i, j) += 2;
         }

         d_cnt++;
         total_cnt++;
      }
      // printf("\n");
   }
}

/* test to cycle through all diags in reverse */
void bck_test_cycle(int Q, int T,
                    float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                    float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                    TRACEBACK* tr)
{
   int d, i, j, k;
   int lb, rb, le, re;                    /* right/left bounds and edges */
   int num_cells;                         /* number of cells in current diagonal */
   int d_cnt, d_total_cnt, total_cnt;     /* count of diags, count of cells in current diag, count of total cells */
   int d_st, d_end;                       /* starting and ending diagonal indices */
   int dim_min, dim_max;                  /* diagonal index where num cells reaches highest point and begins diminishing */
   int dim_T, dim_Q, dim_TOT;             /* dimensions of submatrix being searched */
   COORDS start, end;                     /* start and end point of alignment */

   total_cnt = 0;

   /* INIT ANTI-DIAGS */
   /* test coords */
   start = tr->first_m;
   end = tr->last_m;

   // diag index at corners of dp matrix
   d_st = 0;
   d_end = Q + T;

   // diag index of different start points, creating submatrix
   // d_st = start.i + start.j;
   d_end = end.i + end.j;

   // dimension of submatrix
   dim_TOT = Q + T;
   dim_Q = end.i;
   dim_T = end.j;

   // diag index where num cells reaches highest point and begins diminishing
   dim_min = min(end.i, end.j);
   dim_max = max(end.i, end.j);

   // set bounds using starting cell
   lb = end.i;
   rb = end.i + 1;
   num_cells = 0;

   // printf("test cycle. (Q=%d,T=%d) (dim_Q=%d, dim_T=%d)...\n", Q, T, dim_Q, dim_T);

   /* iterate through diags */
   for (d = d_end; d >= d_st; d--, lb--)
   {
      // printf("d=%d, dim_min=%d, dim_max=%d, num_cells=%d \n", d, dim_min, dim_max, num_cells);

      d_cnt = 0;

      /* is dp matrix diagonal growing or shrinking? */
      if (d >= dim_max)
         num_cells++;
      if (d < dim_min)
         num_cells--;

      // printf("[%d,%d]\n", lb, rb);

      /* find diag cells that are inside matrix bounds */
      le = max(end.i - (d_end - d), 0);
      re = le + num_cells;

      // printf("[%d,%d]\n", le, re);

      /* for testing, bounds = edges */
      if (lb < le)
         lb = le;
      if (rb > re)
         rb = re;

      // printf("[%d,%d]\n", lb, rb);

      /* iterate through cells of diag */
      for (k = lb; k < rb; k++)
      {
         i = k;
         j = d - i;
         // printf("(%d,%d)->", i, j);

         MMX(i, j) += 1;
         // MMX(i, j) = i;
         // IMX(i, j) = j;
         // DMX(i, j) = total_cnt;

         if ( (i - j) % 7 == 0 ) {
            MMX(i, j) += 1;
         }

         d_cnt++;
         total_cnt++;
      }
      // printf("\n");
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

            MMX(i, j) += 1;
            IMX(i, j) = 1;
            DMX(i, j) += val;
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
               printf("OVER MAX: (%d,%d) -> (%d,%d)\n", Q, T, i, j);

            MMX(i, j) += 1;
            IMX(i, j) = 1;
            DMX(i, j) += val;
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