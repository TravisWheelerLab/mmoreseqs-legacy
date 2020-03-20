/*******************************************************************************
 *  FILE:      merge_reorient_linear.c
 *  SYNOPSIS:  Functions for merging multiple EDGEBOUND objects and 
 *             reorienting from diagonal-wise to row-wise.
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
#include "utility.h"
#include "testing.h"

/* objects */
#include "objects/edgebound.h"
#include "objects/hmm_profile.h"
#include "objects/sequence.h"
#include "objects/alignment.h"

/* file parsers */
#include "hmm_parser.h"

/* header */
#include "merge_reorient_linear.h"

/*
 *  FUNCTION: EDGEBOUNDS_Reflect()
 *  SYNOPSIS: Reflect antidiagonal bounds.
 *
 *  ARGS:      <bnd>      Bounds Object
 *
 *  RETURN:    No Return.
 */
void EDGEBOUNDS_Reflect(EDGEBOUNDS *edg) 
{
   int d,lb,rb,rb_new,lb_new;
   for (int i = 0; i < edg->N; i++) {
      d =  edg->bounds[i].id;
      lb = edg->bounds[i].lb;
      rb = edg->bounds[i].rb;

      rb = d - lb + 1;
      lb = d - rb;
      edg->bounds[i] = (BOUND){d,lb,rb};
   }
}


/*
 *  FUNCTION: EDGEBOUNDS_Merge()
 *  SYNOPSIS: Combine two edgebound lists into one.
 *
 *  ARGS:      <Q>            Query length,
 *             <T>            Target length,
 *             <edg_1>        Forward Edgebounds  (by-diag, sorted ascending),
 *             <edg_2>        Backward Edgebounds (by-diag, sorted ascending),
 *
 *  RETURN:
 */
EDGEBOUNDS* EDGEBOUNDS_Merge(const int         Q, 
                             const int         T,
                             const EDGEBOUNDS* edg_1,
                             const EDGEBOUNDS* edg_2)
{
   int         i, j;
   int         x, y; 
   int         d, d_max;
   int         lb, rb;
   int         st_x, end_x; 
   int         st_y, end_y;

   bool        has_merged = false;                 /* checks whether to merge or add bounds */
   const int   tol         = 0;                    /* if clouds are within tolerance range, clouds are merged */
   const int   num_input   = 2;                    /* number of input edgebounds */
   int         starts[]    = {0,0};                /* head pointers for current edgebound ranges */
   BOUND       bnd_1;
   BOUND       bnd_2;
   const 
   EDGEBOUNDS* edg;                                /* pointer for looping through all edgebounds */
   const 
   EDGEBOUNDS* edg_in[]  = {edg_1, edg_2};         /* list of edgebound sets */

   EDGEBOUNDS* edg_out   = EDGEBOUNDS_Create();    /* final output edgebound */
   EDGEBOUNDS* edg_cur   = EDGEBOUNDS_Create();    /* new edgebounds which the two old edgebounds will be merged into */

   /* iterate over all diags */
   st_y = 0;
   d_max = (Q+1)+(T+1);
   for (d = 0; d < d_max; d++) 
   {
      end_y = st_y;

      /* iterate over input edgebound lists and find all bounds on current diagonal */
      for (int i = 0; i < num_input; i++) 
      {
         edg = edg_in[i];
         st_x = starts[i];

         /* find all in edg in current diag */
         end_x  = edg->N;
         for (x = st_x; x < edg->N; x++) 
         {
            if (edg->bounds[x].id != d) {
               break;
            } else {
               EDGEBOUNDS_Pushback(edg_cur, edg->bounds[x]);
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
         for (int i = 0; i < edg_cur->N; i++) 
         {
            bnd_1 = edg_cur->bounds[i];
            for (int j = i+1; j < edg_cur->N; j++) 
            {
               bnd_2 = edg_cur->bounds[j];
                /* if bounds overlap, then merge them */
               if ( (bnd_1.lb <= bnd_2.rb) && (bnd_2.lb <= bnd_1.rb) ) 
               {
                  lb = MIN(bnd_1.lb, bnd_2.lb);
                  rb = MAX(bnd_1.rb, bnd_2.rb);
                  EDGEBOUNDS_Insert(edg_cur, (BOUND){d,lb,rb}, i);
                  EDGEBOUNDS_Delete(edg_cur, edg_cur->bounds[edg->N-1], j);
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
         EDGEBOUNDS_Pushback(edg_out, edg_cur->bounds[i]);
      }

      /* empty current list */
      edg_cur->N = 0;
   }
   EDGEBOUNDS_Destroy(edg_cur);

   return edg_out;
}


/*
 *  FUNCTION: EDGEBOUNDS_Reorient()
 *  SYNOPSIS: Change edgebounds coords from by-diagonal to by-row.
 *
 *  ARGS:      <Q>           Query length
 *             <T>           Target length
 *             <edg_in>      Edgebounds (by-diag, sorted ascending)           
 *
 *  RETURN:    <edg_out>     Edgebounds (by-row)          
 */
EDGEBOUNDS* EDGEBOUNDS_Reorient(const int         Q, 
                                const int         T,
                                const EDGEBOUNDS* edg_in)
{
   int       x,y;
   int       y_st, y_end;
   int       i,j,k;
   int       d,lb,rb;
   BOUND     bnd_in     = (BOUND){0,0,0}; 
   BOUND     bnd_out    = (BOUND){0,0,0};
   bool      in_cloud   = false;  
   bool      is_covered = false;   
   const int tol        = 1;         /* max distance between two forks before merging */

   EDGEBOUNDS* edg_out = EDGEBOUNDS_Create();

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
         d = bnd_in.id;

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
               bnd_out.id = i;
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
            d = bnd_out.id + bnd_out.rb;   /* antidiag containing current cell */
            if (in_cloud && bnd_in.id > d ) 
            {
               // printf("CLOUD END: ");
               // bound_Print(bnd_out);

               EDGEBOUNDS_Pushback(edg_out, bnd_out);
               in_cloud = false;
            }
         }
      }

      /* if the end of the row is reached and still in cloud */
      if (in_cloud) {
         // printf("CLOUD END: ");
         // bound_Print(bnd_out);

         EDGEBOUNDS_Pushback(edg_out, bnd_out);
         in_cloud = false;
      }
   }

   return edg_out;
}