/*******************************************************************************
 *  @file merge_reorient_linear.c
 *  @brief  Functions for merging multiple EDGEBOUND objects and reorienting from antidiagonal-wise to row-wise.
 *
 *  @synopsis
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

/* local imports */
#include "structs.h"
#include "misc.h"
#include "hmm_parser.h"
#include "objects/edgebound.h"
#include "merge_reorient_linear.h"

/*
 *  FUNCTION: edgebounds_Reflect()
 *  SYNOPSIS: Reflect antidiagonal bounds.
 *
 *  PURPOSE:
 *
 *  ARGS:      <bnd>      Bounds Object
 *
 *  RETURN:
 */
void edgebounds_Reflect(EDGEBOUNDS *edg)
{
   int d,lb,rb,rb_new,lb_new;
   for (int i = 0; i < edg->N; i++) {
      d = edg->bounds[i].diag;
      lb = edg->bounds[i].lb;
      rb = edg->bounds[i].rb;

      rb = d - lb + 1;
      lb = d - rb;
      edg->bounds[i] = (BOUND){d,lb,rb};
   }
}


/*
 *  FUNCTION: edgebounds_Merge()
 *  SYNOPSIS: Combine two edgebound lists into one.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg_fwd>      Forward Edgebounds  (by-diag, sorted ascending),
 *             <edg_bck>      Backward Edgebounds (by-diag, sorted ascending),
 *             <edg_res>      Result Edgebounds
 *
 *  RETURN:
 */
void edgebounds_Merge(int Q, int T,
                      EDGEBOUNDS *edg_1,
                      EDGEBOUNDS *edg_2,
                      EDGEBOUNDS *edg_out)
{
   int i, j, x, y, d, lb, rb, d_max;
   int st_x, end_x, st_y, end_y;
   bool has_merged;                  /* checks whether to merge or add bounds */
   const int tol = 0;               /* if clouds are within tolerance range, clouds are merged */
   const int num_input = 2;
   EDGEBOUNDS *edg;
   EDGEBOUNDS *edg_in[] = {edg_1, edg_2};
   EDGEBOUNDS *edg_cur = edgebounds_Create();
   int starts[] = {0,0};
   BOUND bnd_1, bnd_2;

   /* iterate over all diags */
   st_y = 0;
   d_max = (Q+1)+(T+1);
   for (d = 0; d < d_max; d++) 
   {
      end_y = st_y;

      /* iterate over input edgebound lists and find all bounds on current diagonal */
      for (i = 0; i < num_input; i++) 
      {
         edg = edg_in[i];
         st_x = starts[i];

         /* find all in edg in current diag */
         for (x = st_x; x < edg->N; x++) 
         {
            if (edg->bounds[x].diag != d) {
               break;
            } else {
               edgebounds_Add(edg_cur, edg->bounds[x]);
               st_x++;
            }
         }
         starts[i] = st_x;
      }

      /* merge down current list as much as possible (loop until no merges occur) */
      has_merged = true;
      while (has_merged)
      {
         has_merged = false;
         /* check if merge possible for all pairs */
         for (i = 0; i < edg_cur->N; i++) 
         {
            bnd_1 = edg_cur->bounds[i];
            for (j = i+1; j < edg_cur->N; j++) 
            {
               bnd_2 = edg_cur->bounds[j];
                /* if bounds overlap, then merge them */
               if ( (bnd_1.lb <= bnd_2.rb) && (bnd_2.lb <= bnd_1.rb) ) 
               {
                  lb = MIN(bnd_1.lb, bnd_2.lb);
                  rb = MAX(bnd_1.rb, bnd_2.rb);
                  edgebounds_Insert(edg_cur, (BOUND){d,lb,rb}, i);
                  edgebounds_Delete(edg_cur, edg_cur->bounds[edg->N-1], j);
                  has_merged = true;
                  continue;
               }
            }
         }
      }

      /* TODO: sort list */
      // edgbounds_Sort(edg_cur);

      /* insert all in current list into the output list */
      for (int i = 0; i < edg_cur->N; i++)
      {
         edgebounds_Add(edg_out, edg_cur->bounds[i]);
      }

      /* empty current list */
      edg_cur->N = 0;
   }
   edgebounds_Destroy(edg_cur);
}


/*
 *  FUNCTION: edgebounds_Reorient()
 *  SYNOPSIS: Change edgebounds coords from by-diagonal to by-row.
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>       
 *             <edg_in>      Edgebounds (by-diag, sorted ascending) (INPUT)
 *             <edg_out>     Edgebounds (by-row)                    (OUTPUT)
 *
 *  RETURN:
 */
void edgebounds_Reorient(int Q, int T,
                         EDGEBOUNDS *edg_in,
                         EDGEBOUNDS *edg_out)
{
   int x, y, y_st, y_end;
   int i,j,k;
   int d,lb,rb;
   BOUND bnd_in = (BOUND){0,0,0}; 
   BOUND bnd_out = (BOUND){0,0,0};
   bool in_cloud = false;  
   bool is_covered = false;   
   const int tol = 1;         /* max distance between two forks before merging */

   /* range for possible antidiags */
   y_st = 0;
   y_end = edg_in->N;

   /* for each row */
   for (x = 0; x < Q+1; x++)
   {
      // printf("== ROW %d == y=(%d,%d)\n", x, y_st, y_end);
      in_cloud = false;

      /* compare against the antidiags */
      for (y = y_st; y < y_end; y++)
      {
         bnd_in = edg_in->bounds[y];
         d = bnd_in.diag;

         /* cartesian coords */
         i = x;
         j = d - i;

         /* TODO: update if at impossible coords (unreachable antidiags) */
         // if (j < 0) {
         //    y_st = y;
         //    continue;
         // }
         // if (j > T+1) {
         //    break;
         // }

         /* check if cell is covered by anti-diag... */
         rb = d - bnd_in.lb;
         lb = d - bnd_in.rb + 1;
         /* TODO: Look at off-by-one bridging and possibly solve without reflection? */
         is_covered = (j >= lb && j <= rb);
         // printf("y:%d, d:%d => (%d,%d) => is %d in (%d,%d) : %s \n", y, d, i, j, j, bnd_in.lb, bnd_in.rb, is_covered ? "YES" : "NO");
         if (is_covered) {
            /* if in cloud, update bounds */
            if (in_cloud) {
               bnd_out.rb = j+1;
            }
            /* if not, create new cloud and bound */
            else 
            {
               bnd_out.diag = i;
               bnd_out.lb = j;
               bnd_out.rb = j+1;
               // printf("CLOUD BEGIN: ");
               // bound_Print(bnd_out);
               in_cloud = true;
            }
         } 
         /* if cell not covered by anti-diag */
         else 
         {
            /* if current antidiag doesn't contain next cell, we're at the end of current cloud */
            d = bnd_out.diag + bnd_out.rb;   /* antidiag containing current cell */
            if (in_cloud && bnd_in.diag > d ) 
            {
               // printf("CLOUD END: ");
               // bound_Print(bnd_out);

               edgebounds_Add(edg_out, bnd_out);
               in_cloud = false;
            }
         }
      }

      /* if the end of the row is reached and still in cloud */
      if (in_cloud) {
         // printf("CLOUD END: ");
         // bound_Print(bnd_out);

         edgebounds_Add(edg_out, bnd_out);
         in_cloud = false;
      }
   }
}