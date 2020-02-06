/*******************************************************************************
 *  @file merge_reorient_quad.c
 *  @brief  Functions for merging multiple EDGEBOUND objects and reorienting from antidiagonal-wise to row-wise.
 *
 *  @synopsis
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
#include "testing.h"
#include "objects/edgebound.h"
#include "merge_reorient_quad.h"

/*
 *  FUNCTION: edgebounds_Merge_Reorient_Quad()
 *  SYNOPSIS: Merge and Reorient edgebounds from Matrix Cloud (NAIVE).
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg_fwd>      Forward Edgebounds,
 *             <edg_bck>      Backward Edgebounds,
 *             <edg_res>      Result Edgebounds,
 *             <st_MX>        State Matrix
 *
 *  RETURN:
 */
int edgebounds_Merge_Reorient_Naive(EDGEBOUNDS*edg_fwd,
                                  EDGEBOUNDS*edg_bck,
                                  EDGEBOUNDS*edg_diag,
                                  EDGEBOUNDS*edg_row,
                                  int Q, int T,
                                  float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                                  float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ])
{
   /* merge edgebounds */
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_fwd, 1, MODE_DIAG);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_bck, 1, MODE_DIAG);
   
   /* reorient from diag-wise to row-wise */
   edgebounds_Build_From_Cloud(edg_row, Q, T, st_MX, MODE_ROW);
   edgebounds_Build_From_Cloud(edg_diag, Q, T, st_MX, MODE_DIAG);

   int num_cells = cloud_Cell_Count(Q, T, st_MX, sp_MX);

   return num_cells;
}


/*
 *  FUNCTION: edgebounds_Build_From_Cloud()
 *  SYNOPSIS: Create edgebounds from a Matrix Cloud.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>        Edgebound,
 *             <st_MX>      State Matrix
 *             <mode>       Diagonal or Row-wise Edgebound
 *
 *  RETURN:
 */
void edgebounds_Build_From_Cloud( EDGEBOUNDS* edg,
                                 int Q, int T,
                                 float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                                 int mode)
{
   int d, i, j, k;
   int x, y1, y2;
   int d_st, d_end;
   int dim_min, dim_max; 
   int le, re;                         
   bool in_cloud;
   int num_cells;
   BOUND bnd;
   
   /* create edgebound in antidiagonal-wise order */
   if (mode == MODE_DIAG) 
   {
      d_st = 0;
      d_end = Q + T;

      dim_min = MIN(Q,T);
      dim_max = MAX(Q,T);

      le = 0;
      re = 1;
      num_cells = 0;

      /* iterate through diags */
      for (d = d_st; d <= d_end; d++)
      {
         /* is dp matrix diagonal growing or shrinking? */
         if (d <= dim_min)
            num_cells++;
         if (d > dim_max)
            num_cells--;

         /* find diag cells that are inside matrix bounds */
         le = MAX(0, d - T);
         re = le + num_cells;

         /* iterate through cells of diag */
         for (k = le; k < re; k++)
         {
            i = k;
            j = d - i;

            if (in_cloud)
            {
               /* end of bound, add bound to edgebound list */
               if ( IMX(i,j) <= 0 ) 
               {
                  in_cloud = false;
                  y2 = k;
                  edgebounds_Add(edg, (BOUND) {x,y1,y2});
               }
            }
            else
            {
               /* start of new bound */
               if ( IMX(i,j) > 0 ) {
                  in_cloud = true;
                  x = d;
                  y1 = k;
               }
            }
         }

         /* if still in bound at end of diag, then at end of of bound, so add to edgebound list */
         if (in_cloud) 
         {
            in_cloud = false;
            y2 = re;
            edgebounds_Add(edg, (BOUND) {x,y1,y2});
         }
      }
   }
   /* create edgebound in row-wise order */
   else 
   if (mode == MODE_ROW) 
   {
      printf("");
      /* iterate through rows */
      for (i = 0; i <= Q; i++)
      {
         /* iterate through cells in row */
         for (j = 0; j <= T; j++)
         {
            if (in_cloud)
            {
               /* end of bound, add bound to edgebound list */
               if ( IMX(i,j) <= 0 ) 
               {
                  in_cloud = false;
                  y2 = j;

                  edgebounds_Add(edg, (BOUND) {x,y1,y2});
               }
            }
            else
            {
               /* start of new bound */
               if ( IMX(i,j) > 0 ) {
                  in_cloud = true;
                  x = i;
                  y1 = j;
                  // printf("Starting [%d]:%d,(%d,?)...\n", edg->N, x, y1);
               }
            }
         }

         /* if still in bound at end of row, then at end of of bound, so add to edgebound list */
         if (in_cloud) 
         {
            in_cloud = false;
            y2 = T+1;

            edgebounds_Add(edg, (BOUND) {x,y1,y2});
         }
      }
   }
}