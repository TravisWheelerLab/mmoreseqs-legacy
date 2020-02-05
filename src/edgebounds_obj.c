/*******************************************************************************
 *  @file edgebounds_obj.c
 *  @brief Functions for EDGEBOUNDS object.
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
#include "viterbi.h"
#include "forward_backward.h"
#include "testing.h"
#include "edgebounds_obj.h"

/*
 *  FUNCTION: edgebounds_Create()
 *  SYNOPSIS: Create new edgebounds object and return pointer.
 *
 *  PURPOSE:   
 *
 *  ARGS:      
 *
 *  RETURN:    <edg>      Edgebounds Object
 */
EDGEBOUNDS *edgebounds_Create()
{
   EDGEBOUNDS *edg;
   const int min_size = 16;
   edg = (EDGEBOUNDS *)malloc(sizeof(EDGEBOUNDS));
   edg->N = 0;
   edg->size = min_size;
   edg->bounds = (BOUND *)malloc(min_size * sizeof(BOUND));
   return edg;
}

/*
 *  FUNCTION: edgebounds_Init()
 *  SYNOPSIS: Initialize new edgebounds object pointer.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Init(EDGEBOUNDS **edg)
{
   const int min_size = 16;
   (*edg) = (EDGEBOUNDS *)malloc(sizeof(EDGEBOUNDS));
   (*edg)->N = 0;
   (*edg)->size = min_size;
   (*edg)->bounds = (BOUND *)malloc(min_size * sizeof(BOUND));
}

/*
 *  FUNCTION: edgebounds_Destroy()
 *  SYNOPSIS: Destroy edgebounds object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Destroy(EDGEBOUNDS *edg)
{
   free(edg->bounds);
   free(edg);
}


/*
 *  FUNCTION: edgebounds_Add()
 *  SYNOPSIS: Add bound to Edgebound list.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>       Edgebounds,
 *             <bnd>       Bound
 *
 *  RETURN:
 */
void edgebounds_Add(EDGEBOUNDS *edg,
                    BOUND bnd)
{
   /* resize if necessary */
   edg->N += 1;
   if (edg->N >= edg->size) edgebounds_Resize(edg);

   edg->bounds[edg->N-1].diag = bnd.diag;
   edg->bounds[edg->N-1].lb = bnd.lb;
   edg->bounds[edg->N-1].rb = bnd.rb;
}


/*
 *  FUNCTION: edgebounds_Resize()
 *  SYNOPSIS: Resize number of bounds in edgebound object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Resize(EDGEBOUNDS *edg)
{
   const int growth_factor = 2;
   edg->size *= growth_factor;
   edg->bounds = (BOUND *)realloc(edg->bounds, edg->size * sizeof(BOUND));
}


/*
 *  FUNCTION: edgebounds_Print()
 *  SYNOPSIS: Print EDGEBOUND object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Print(EDGEBOUNDS *edg)
{
   printf("printing edgebounds...\n");
   printf("N: %d, Nalloc: %d\n", edg->N, edg->size);
   for (unsigned int i = 0; i < edg->N; ++i)
   {
      printf("[%d] ", i);
      bound_Print(edg->bounds[i]);
   }
}


/*
 *  FUNCTION: bound_Print()
 *  SYNOPSIS: Print BOUND object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <bnd>      Bounds Object
 *
 *  RETURN:
 */
void bound_Print(BOUND bnd)
{
   printf("d: %d, lb: %d, rb: %d\n", bnd.diag, bnd.lb, bnd.rb);
}


/*
 *  FUNCTION: edgebounds_Save()
 *  SYNOPSIS: Save edgebound printout to file.
 *
 *  PURPOSE:
 *
 *  ARGS:      <bnd>      Bounds Object
 *             <f>        Filename
 *
 *  RETURN:
 */
void edgebounds_Save(EDGEBOUNDS *edg,
                      const char *_filename_)
{
   FILE *fp;
   fp = fopen(_filename_, "w");

   for (unsigned int i = 0; i < edg->N; ++i)
   {
      fprintf(fp, "[%d] x: %d\t y: (%d, %d)\n", i, edg->bounds[i].diag, edg->bounds[i].lb, edg->bounds[i].rb);
   }
   fclose(fp);
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
                      EDGEBOUNDS *edg_new)
{
   int i, x, y, d, lb, rb, d_max;
   int st_x, end_x, st_y, end_y;
   bool is_merged;                  /* checks whether to merge or add bounds */
   const int tol = 0;               /* if clouds are within tolerance range, clouds are merged */
   const int num_input = 2;
   EDGEBOUNDS *edg;
   EDGEBOUNDS *edg_in[] = {edg_1, edg_2};
   int starts[] = {0,0};
   BOUND bnd_in = (BOUND){0,0,0};
   BOUND bnd_out = (BOUND){0,0,0};

   /* iterate over all diags */
   st_y = 0;
   d_max = (Q+1)+(T+1);
   for (d = 0; d < d_max; d++) 
   {
      end_y = st_y;

      /* iterate over input edgebounds */
      for (i = 0; i < num_input; i++) 
      {
         edg = edg_in[i];
         st_x = starts[i];

         /* find all in edg in current diag */
         for (x = st_x; x < edg->N; x++) 
         {
            if (edg->bounds[x].diag != d) {
               break;
            }
         }
         end_x = x;

         /* integrate bounds from edg into edg_new */
         for (x = st_x; x < end_x; x++) 
         {
            bnd_in = edg->bounds[x];

            is_merged = false;
            for (y = st_y; y < edg_new->N; y++) 
            {
               bnd_out = edg_new->bounds[y];
               /* if bounds intersect, then merge */
               if ( (bnd_in.lb <= bnd_out.rb) && (bnd_out.lb <= bnd_in.rb) ) 
               {
                  lb = MIN(bnd_in.lb, bnd_out.lb);
                  rb = MAX(bnd_in.rb, bnd_out.rb);
                  edg_new->bounds[y] = (BOUND){d,lb,rb};
                  is_merged = true;
               }
            } 

            /* if there were no bounds to merge into, add it separately */
            if (!is_merged) {
               edgebounds_Add(edg_new, bnd_in);
               end_y++;
            }
         }
         starts[i] = end_x;
      }
      st_y = end_y;
   }
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
   int x,y;
   int i,j,k;
   int d,lb_1,rb_1,lb_2,rb_2;
   BOUND bnd_in = (BOUND){0,0,0}; 
   BOUND bnd_out = (BOUND){0,0,0};
   bool merged = false;
   bool in_cloud = false;     
   const int tol = 1;         /* max distance between two forks before merging */

   /* for each row */
   for (x = 0; x < (Q+1); x++)
   {
      in_cloud = false;
      for (y = 0; y < edg_in->N; y++)
      {
         bnd_in = edg_in->bounds[y];
         /* cartesian coords */
         i = x;
         j = bnd_in.diag - i;

         /* check if cell is covered by anti-diag... */
         if (j >= bnd_in.lb && j < bnd_in.rb) {
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
               edgebounds_Add(edg_out, bnd_out);
               in_cloud = false;
            }
         }
      }

      /* if the end of the row is reached and still in cloud */
      if (in_cloud) {
         edgebounds_Add(edg_out, bnd_out);
         in_cloud = false;
      }
   }
}


/*
 *  FUNCTION: edgebounds_Merge_Reorient()
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
int edgebounds_Merge_Reorient_Cloud(EDGEBOUNDS*edg_fwd,
                                  EDGEBOUNDS*edg_bck,
                                  EDGEBOUNDS*edg_new,
                                  int Q, int T,
                                  float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                                  float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ])
{
   /* merge edgebounds */
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_fwd, 1, MODE_DIAG);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_bck, 1, MODE_DIAG);
   
   /* reorient from diag to row-wise */
   edgebounds_Build_From_Cloud(edg_new, Q, T, st_MX, MODE_ROW);
   int num_cells = cloud_Cell_Count(Q, T, st_MX, sp_MX);

   return num_cells;
}


/*
 *  FUNCTION: edgebounds_Reorient_Cloud()
 *  SYNOPSIS: Reorient edgebounds from Matrix Cloud.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Forward Edgebounds,
 *             <edg_res>      Result Edgebounds,
 *             <st_MX>        State Matrix
 *
 *  RETURN:
 */
void edgebounds_Reorient_Cloud( EDGEBOUNDS*edg_old,
                               EDGEBOUNDS*edg_new,
                               int Q, int T,
                               float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                               float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ], 
                               int old_mode, int new_mode)
{
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_old, 1, old_mode);
   edgebounds_Build_From_Cloud(edg_new, Q, T, st_MX, new_mode);
}

/*
 *  FUNCTION: edgebounds_Merge_Cloud()
 *  SYNOPSIS: Reorient edgebounds from Matrix Cloud.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg_1>        Edgebound Set 1,
 *             <edg_2>        Edgebound Set 2,
 *             <edg_res>      Result Edgebounds,
 *             <st_MX>        State Matrix
 *
 *  RETURN:
 */
void edgebounds_Merge_Cloud( EDGEBOUNDS*edg_1,
                            EDGEBOUNDS*edg_2,
                            EDGEBOUNDS*edg_res,
                            int Q, int T,
                            float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                            float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                            int mode)
{
   /* initialize new edgebound */
   edgebounds_Create(edg_res);
   // static int min_size = 128;
   // edg_res->N = 0;
   // edg_res->size = min_size;
   // edg_res->bounds = (BOUND *)malloc(min_size * sizeof(BOUND));

   /* merge edgebounds */
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_1, 1, mode); 
   cloud_Fill(Q, T, st_MX, sp_MX, edg_2, 1, mode); 
   
   edgebounds_Build_From_Cloud(edg_res, Q, T, st_MX, mode);
}

/*
 *  FUNCTION: edgebounds_Build_From_Cloud()
 *  SYNOPSIS: Create edgebounds from a Matrix Cloud.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>        Edgebound,
 *             <st_MX>        State Matrix
 *             <mode>         Diagonal or Row-wise Edgebound
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
   
   /* create edgebound in antidiagonal-wise order */
   if (mode == MODE_DIAG) 
   {
      d_st = 0;
      d_end = (Q+1) + (T+1) - 1;
      dim_min = MIN(Q,T);
      dim_max = MAX(Q,T);

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
                  edg->bounds[edg->N] = (BOUND) {x, y1, y2};
                  edg->N++;
                  if (edg->N >= edg->size-1) {
                     edg->size *= 2;
                     edg->bounds = (BOUND *)realloc(edg->bounds, edg->size * sizeof(BOUND));
                  }
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
            edg->bounds[edg->N] = (BOUND) {x, y1, y2};
            edg->N++;
            if (edg->N >= edg->size-1) {
               edg->size *= 2;
               edg->bounds = (BOUND *)realloc(edg->bounds, edg->size * sizeof(BOUND));
            }
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
                  edg->bounds[edg->N].diag = x;
                  edg->bounds[edg->N].lb = y1;
                  edg->bounds[edg->N].rb = y2;
                  edg->N += 1;
                  // printf("N=%d X=%d\n", x, edg->bounds[0].diag);
                  if (edg->N >= edg->size-1) {
                     edg->size *= 2;
                     edg->bounds = (BOUND *)realloc(edg->bounds, edg->size * sizeof(BOUND));
                  }
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
            edg->bounds[edg->N].diag = x;
            edg->bounds[edg->N].lb = y1;
            edg->bounds[edg->N].rb = y2;
            edg->N += 1;
            if (edg->N >= edg->size-1) {
               edg->size *= 2;
               edg->bounds = (BOUND *)realloc(edg->bounds, edg->size * sizeof(BOUND));
            }
         }
      }
   }
}

/*
 *  FUNCTION: bounds_Compare()
 *  SYNOPSIS: Compare two Bounds, first by diag, then by lb, then by rb
 *
 *  PURPOSE:
 *
 *  ARGS:      <a>        Bound,
 *             <b>        Bound
 *
 *  RETURN:    1 if (a > b), 0 if equal, -1 if (a < b)
 */
int bounds_Compare(BOUND a, BOUND b)
{
   if (a.diag > b.diag) {
      return 1;
   } else 
   if (a.diag < b.diag) {
      return -1;
   }

   if (a.lb > b.lb) {
      return 1;
   } else
   if (a.lb < b.lb) {
      return -1;
   }

   if (a.rb > b.rb) {
      return 1;
   } else
   if (a.rb < b.rb) {
      return -1;
   }
   
   return 0;
}


