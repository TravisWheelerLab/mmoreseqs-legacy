/*******************************************************************************
 *  @file forward_backward.c
 *  @brief Functions for EDGEBOUNDS object.
 *
 *  @synopsis
 *
 *  @author Dave Rich (devrek)
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

// macros
#define getName(var) #var
#define SCALE_FACTOR 1000

// macro functions
// NOTE: wrap all macro vars in parens!!
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))


/*
 *  FUNCTION: edgebounds_Create()
 *  SYNOPSIS: Create new edgebounds object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Create(EDGEBOUNDS *edg)
{
   static int min_size = 128;
   edg = (EDGEBOUNDS *)malloc(sizeof(EDGEBOUNDS));
   edg->N = 0;
   edg->size = min_size;
   edg->bounds = (BOUND *)malloc(min_size * sizeof(BOUND));
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
 *  SYNOPSIS: Add bound to Edgebound
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>       Edgebounds,
 *             <bnd>       Bound
 *
 *  RETURN:
 */
void edgebounds_Add(EDGEBOUNDS *edg,
                    BOUND *bnd)
{
   edg->N++;

   if (edg->N == edg->size) {
      edg->bounds = (BOUND *)realloc(edg->bounds, edg->size * 2);
   }

   edg->bounds[edg->N] = *bnd;
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
   int size = edg->size * 2;
   edg->bounds = (BOUND *)realloc(edg->bounds, size * sizeof(BOUND));
}


/*
 *  FUNCTION: edgebounds_Print()
 *  SYNOPSIS: Print edgebound object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <bnd>      Bounds Object
 *
 *  RETURN:
 */
void edgebounds_Print(EDGEBOUNDS *edg)
{
   printf("printing edgebounds...\n");
   for (unsigned int i = 0; i < edg->N; ++i)
   {
      printf("[%d] d: %d, lb: %d, rb: %d\n", i, edg->bounds[i].diag, edg->bounds[i].lb, edg->bounds[i].rb);
   }
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
 *  ARGS:      <edg_fwd>      Forward Edgebounds (sorted by diag ascending),
 *             <edg_bck>      Backward Edgebounds (sorted by diag ascending),
 *             <edg_res>      Result Edgebounds
 *
 *  RETURN:
 */
void edgebounds_Merge(EDGEBOUNDS *edg_fwd,
                      EDGEBOUNDS *edg_bck,
                      EDGEBOUNDS *edg_new)
{
   int i, j, k, x;
   bool not_merged;          /* checks whether to merge or add bounds */
   int diag_cur;             /* currently examined diagonal */
   static int tol = 0;       /* if clouds are within tolerance range, clouds are merged */

   /* allocate new edgebounds */
   int min_size = 128;
   edg_new = (EDGEBOUNDS *)malloc( sizeof(EDGEBOUNDS) );
   edg_new->size = min_size;
   edg_new->N = 0;
   edg_new->bounds = (BOUND *)malloc( min_size * sizeof(BOUND) );
   // edgebounds_Create(edg_new);

   /* temp edgebound for merging current diagonal */
   EDGEBOUNDS *edg_tmp;
   edg_tmp = (EDGEBOUNDS *)malloc( sizeof(EDGEBOUNDS) );
   edg_tmp->size = min_size;
   edg_tmp->N = 0;
   edg_tmp->bounds = (BOUND *)malloc( min_size * sizeof(BOUND) );
   // edgebounds_Create(edg_tmp);

   printf("FORWARD:\n");
   edgebounds_Print(edg_fwd);
   printf("BACKWARD:\n");
   edgebounds_Print(edg_bck);

   BOUND bound_cur;

   /* begin with minimum diagonal */
   diag_cur = min(edg_fwd->bounds[0].diag, edg_bck->bounds[0].diag);

   printf("TOTALS: fwd: %d, bck: %d\n", edg_fwd->N, edg_bck->N);

   /* iterate over all bounds */
   i = 0; j = 0; k = 0;
   while (i < edg_fwd->N || j < edg_bck->N)
   {
      printf("diag_cur: %d, diag_fwd: %d, diag_bck: %d\n", diag_cur, edg_fwd->bounds[i].diag, edg_bck->bounds[j].diag);

      printf("edg_tmp(k): %d, edg_fwd(i): %d, edg_bck(j): %d\n", k, i, j);

      /* merge all forward bounds from current diagonal */
      while (i < edg_fwd->N && edg_fwd->bounds[i].diag == diag_cur)
      {
         /* merge bounds if applicable */
         bound_cur = edg_fwd->bounds[i];
         not_merged = true;
         for (x = 0; x < k; ++x)
         {
            /* if bounds overlap, merge it into and break from loop */
            if ( !(bound_cur.rb - edg_tmp->bounds[x].lb <= tol || bound_cur.lb - edg_tmp->bounds[x].rb >= tol ) )
            {
               printf("merge fwd (%d)...\n", i);
               not_merged = false;
               edg_tmp->bounds[x].lb = min(bound_cur.lb, edg_tmp->bounds[x].lb);
               edg_tmp->bounds[x].rb = max(bound_cur.rb, edg_tmp->bounds[x].rb);
               break;
            }
         }
         /* if bounds don't overlap with any previous, add new bounds  */
         if (not_merged)
         {
            printf("add fwd (%d)...\n", i);
            edg_tmp->bounds[k] = edg_fwd->bounds[i];
            edg_tmp->N++;
            ++k;
         }
         ++i;
      }

      printf("edg_tmp(k): %d, edg_fwd(i): %d, edg_bck(j): %d\n", k, i, j);
      /* merge all backward bounds from current diagonal */
      while (j < edg_bck->N && edg_bck->bounds[i].diag == diag_cur)
      {
         /* get next bound */
         bound_cur = edg_bck->bounds[i];
         not_merged = true;
         /* compare against every other  */
         for (x = 0; x < k; ++x)
         {
            /* if bounds overlap, merge it into */
            if ( !(bound_cur.rb - edg_tmp->bounds[x].lb <= tol || bound_cur.lb - edg_tmp->bounds[x].rb >= tol ) )
            {
               printf("merge bck (%d)...\n", j);
               not_merged = false;
               edg_tmp->bounds[x].lb = min(bound_cur.lb, edg_tmp->bounds[x].lb);
               edg_tmp->bounds[x].rb = max(bound_cur.rb, edg_tmp->bounds[x].rb);
               break;
            }
         }
         /* if bounds don't overlap with any previous, add new bounds  */
         if (not_merged)
         {
            printf("add bck (%d)...\n", j);
            edg_tmp->bounds[k] = edg_bck->bounds[j];
            edg_tmp->N++;
            ++k;
         }
         ++j;
      }

      int curr_diag = min(edg_fwd->bounds[i].diag, edg_bck->bounds[j].diag);
      exit(0);
   }
}


/*
 *  FUNCTION: edgebounds_Reorient()
 *  SYNOPSIS: Change edgebounds from by-diagonal to by-row.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg_diag>      Edgebounds (by-diag)
 *             <edg_row>       Edgebounds (by-row)
 *
 *  RETURN:
 */
void edgebounds_Reorient(EDGEBOUNDS *edg_src,
                         EDGEBOUNDS *edg_dest)
{
   int i,j;
   edgebounds_Create(edg_dest);

   /* convert edgebounds from (diag, leftbound, rightbound) to {(x1,y1),(x2,y2)} coords */
   for (i = 0; i < edg_src->N; ++i)
   {

   }
}


/*
 *  FUNCTION: edgebounds_Merge_Reorient()
 *  SYNOPSIS: Merge and Reorient edgebounds from Matrix Cloud.
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
   /* initialize new edgebound */
   static int min_size = 128;
   edg_new->N = 0;
   edg_new->size = min_size;
   edg_new->bounds = (BOUND *)malloc(min_size * sizeof(BOUND));

   /* merge edgebounds */
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
 *  ARGS:      <edg_1>        Edgebound,
 *             <st_MX>        State Matrix
 *             <mode>         Diagonal or Row-wise Edgebound
 *
 *  RETURN:
 */
void edgebounds_Build_From_Cloud( EDGEBOUNDS*edg,
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

   /* initialize new edgebound */
   static int min_size = 128;
   edg->N = 0;
   edg->size = min_size;
   edg->mode = mode;
   edg->bounds = (BOUND *)malloc(min_size * sizeof(BOUND));
   
   /* create edgebound in antidiagonal-wise order */
   if (mode == MODE_DIAG) 
   {
      d_st = 0;
      d_end = Q + T + 1;
      dim_min = min(Q,T);
      dim_max = max(Q,T);

      /* iterate through diags */
      for (d = d_st; d <= d_end; d++)
      {
         /* is dp matrix diagonal growing or shrinking? */
         if (d <= dim_min)
            num_cells++;
         if (d > dim_max)
            num_cells--;

         /* find diag cells that are inside matrix bounds */
         le = max(0, d - T);
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
   // printf("TEST: [0].diag = %d\n", edg->bounds[0].diag);
}





