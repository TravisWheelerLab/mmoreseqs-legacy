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

/* local imports */
#include "../objects/structs.h"
#include "../utilities/utilities.h"
#include "../objects/objects.h"

#include "algs_linear.h"

/* header */
#include "merge_reorient_linear.h"

/*
 *  FUNCTION:  EDGEBOUNDS_Reflect()
 *  SYNOPSIS:  Reflect antidiagonal bounds.
 */
void EDGEBOUNDS_Reflect(EDGEBOUNDS *edg) 
{
   int d, lb, rb, rb_new, lb_new;

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
 *  FUNCTION:  EDGEBOUNDS_Merge_Together()
 *  SYNOPSIS:  Combine two edgebound lists into one. 
 *             Assumes input lists are sorted and both oriented by-antidiagonal.
 */
void EDGEBOUNDS_Merge_Together(     const int           Q,             /* query length */
                                    const int           T,             /* target length */
                                    const EDGEBOUNDS*   edg_in_1,      /* edgebounds (fwd, sorted ascending) */
                                    const EDGEBOUNDS*   edg_in_2,      /* edgebounds (bck, sorted ascending) */
                                    EDGEBOUNDS*         edg_out )      /* OUTPUT: merged edgebounds (sorted ascending) */
{
   int         i, j;                /* indexes */
   int         x, y;                /* indexes */
   int         d, d_max;            /* diagonals */
   int         lb, rb;              /* left and right bound */
   int         st_x, end_x;         /* range of x indexes */
   int         st_y, end_y;         /* range of y indexes */

   bool        has_merged  = false;                /* checks whether to merge or add bounds */
   const int   tol         = 0;                    /* if clouds are within tolerance range, clouds are merged */
   const int   num_input   = 2;                    /* number of input edgebounds (current locked to two) */
   int         starts[]    = {0,0};                /* head pointers for current edgebound ranges */

   /* edgebound meta data */

   /* pointer for looping through all edgebounds */   
   const EDGEBOUNDS* edg;
   BOUND             bnd_1;
   BOUND             bnd_2;

   /* list of edgebound sets */
   const EDGEBOUNDS* edg_in[]  = {edg_in_1, edg_in_2};

   /* reset output edgebounds */
   EDGEBOUNDS_Reuse( edg_out, Q, T );
   
   /* verify that all input edgebounds are the same mode */
   #if DEBUG 
   {
      for ( i = 0; i < num_input - 1; i++ ) {
         if ( edg_in[i]->edg_mode != edg_in[i+1]->edg_mode ) 
         {
            fprintf( stderr, "ERROR: Not all edgebounds being merged have same orientation!");
            exit(EXIT_FAILURE);
         }
      }
   }
   #endif
   edg_out->edg_mode = edg_in_1->edg_mode;

   /* new edgebounds which the two old edgebounds will be merged into */
   EDGEBOUNDS* edg_cur  = EDGEBOUNDS_Create();

   /* assert edgebounds are sorted */
   #if DEBUG
   {
      for ( int i = 0; i < num_input; i++ ) 
      {
         edg = edg_in[i];

         if ( edg->N > 1)
         {
            for ( int j = 1; j < edg->N; j++ ) 
            {
               if ( edg->bounds[j].id < edg->bounds[j-1].id ) {
                  fprintf(stderr, "ERROR: edgebounds are not sorted.\n");
                  fprintf(stderr, "EDG->N = %d\n", edg->N );
                  fprintf(stderr, "EDG[%d]: [%d].id = %d VS [%d].id = %d\n", 
                     i, j-1, edg->bounds[j-1].id, j, edg->bounds[j].id );

                  EDGEBOUNDS_Dump( edg, stderr );
                  
                  exit(EXIT_FAILURE);
               }
            }
         }
      }
   }
   #endif 

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
            /* if edge has surpassed the current diagonal, then break */
            if (edg->bounds[x].id > d) {
               break;
            } 
            // /* if edge preceeds the current, skip to next diagonal (should never occur) */
            // else if (edg->bounds[x].id < d) {
            //    continue;
            // }
            /* otherwise, add edgebound to list and increment starting location */
            else {
               EDGEBOUNDS_Pushback(edg_cur, &(edg->bounds[x]) );
               st_x++;
            }
         }
         /* update next starting location to the end location of the last */
         starts[i] = st_x;
      }

      /* merge down current list as much as possible (loop until no merges occur) */
      has_merged = true;
      while ( has_merged )
      {
         has_merged = false;
         /* check if merge possible for all pairs (n^2 in number of bounds in ) */
         for (int i = 0; i < edg_cur->N; i++) 
         {
            bnd_1 = edg_cur->bounds[i];
            for (int j = i+1; j < edg_cur->N; j++) 
            {
               bnd_2 = edg_cur->bounds[j];
                /* if bounds overlap, then merge them */
               if ( (bnd_1.lb <= bnd_2.rb) && (bnd_2.lb <= bnd_1.rb) ) 
               {
                  /* take mins and maxes to build new range */
                  lb = MIN(bnd_1.lb, bnd_2.lb);
                  rb = MAX(bnd_1.rb, bnd_2.rb);
                  /* replace the edgebound at first index */
                  EDGEBOUNDS_Insert(edg_cur, i, &( (BOUND){d,lb,rb} ) );
                  /* delete the edgebound at the second index (backfills from end of list) */
                  EDGEBOUNDS_Delete(edg_cur, j);
                  /* flag that a merge has occurred */
                  has_merged = true;
                  /* go back to start of list */
                  continue;
               }
            }
         }
      }

      /* TODO: sort list (not necessary, since inputs assumed to be sorted) */
      // edgbounds_Sort(edg_cur);

      /* insert all in current list into the output list */
      for (int i = 0; i < edg_cur->N; i++)
      {
         EDGEBOUNDS_Pushback( edg_out, &(edg_cur->bounds[i]) );
      }

      /* empty current list */
      edg_cur->N = 0;
   }
   EDGEBOUNDS_Destroy( edg_cur );
}


/*
 *  FUNCTION: EDGEBOUNDS_Reorient_to_Row()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.       
 */
void EDGEBOUNDS_Reorient_to_Row( const int           Q,          /* query length */
                                 const int           T,          /* target length */
                                 EDGEBOUNDS*         edg_in,     /* edgebounds (antidiag-wise, sorted ascending) */
                                 EDGEBOUNDS*         edg_out )   /* OUPUT: edgebounds (row-wise, sorted ascending) */
{
   int       x, y;                              /* indexes into edgebounds */
   int       y_st, y_end;                       /* index range */
   int       i,j,k;                             /* indexes for diagonal and x-y coords */
   int       d,lb,rb;                           /* diag/row, left-bound, right-bound */
   BOUND     bnd_in     = (BOUND){0,0,0};       /* */
   BOUND     bnd_out    = (BOUND){0,0,0};       /* */
   bool      in_cloud   = false;                /* flags whether currently extending bounds along a row */
   bool      is_covered = false;                /* */
   const int tol        = 0;                    /* max distance between two row indexes to merge */

   /* first, verify that input edgebounds are stored by-diagonal */
   if ( edg_in->edg_mode != EDG_DIAG ) {
      // /* if not, then swap pointers and return */
      EDGEBOUNDS* swap = edg_in;
      edg_in = edg_out;
      edg_out = swap;
      /* if not, make a copy of input edgebounds for output */
      // EDGEBOUNDS_Copy( edg_out, edg_in );
      return;
   }

   /* reuse edgebounds */
   EDGEBOUNDS_Reuse( edg_out, Q, T );
   edg_out->edg_mode = EDG_ROW;

   /* range for possible antidiag indexes */
   y_st = 0;
   y_end = edg_in->N;

   /* for each row in dp_matrix */
   for (x = 0; x < Q+1; x++)
   {
      in_cloud = false;

      /* compare against the antidiags at the current row */
      for (y = y_st; y < y_end; y++)
      {
         bnd_in   = edg_in->bounds[y];
         d        = bnd_in.id;

         /* x-y coords */
         i = x;         /* row index */
         j = d - i;     /* column index == offset in antidiagonal */

         /* if j < 0, it is outside the left-bounds of the dp_matrix */
         if (j < 0) {
            /* if this antidiag is outside the dp_matrix now, it will be outside for all lower rows as well; begin future row's search at this antidiag */
            y_st = y;
            continue;
         }
         /* if j > T+1, it is output the right-bounds of the dp_matrix */
         if (j > T+1) {
            /* if this antidiag is outside the dp_matrix now, all subsequent antidiags will be as well; so go to next row */
            break;
         }

         /* check if cell is covered by anti-diag... */
         rb = d - bnd_in.lb;
         lb = d - bnd_in.rb + 1;
         is_covered = (j >= lb && j <= rb);
         // printf("y:%d, d:%d => (%d,%d) => is %d in (%d,%d) : %s \n", y, d, i, j, j, bnd_in.lb, bnd_in.rb, is_covered ? "YES" : "NO");
         
         /* if antidiag range covers cell in row... */
         if (is_covered) {
            /* if in cloud, update bounds */
            if (in_cloud) {
               bnd_out.rb = j + 1;
            }
            /* if not, create new cloud and bound */
            else 
            {
               bnd_out.id = i;
               bnd_out.lb = j;
               bnd_out.rb = j + 1;
               in_cloud = true;
            }
         } 
         /* if cell not covered by anti-diag */
         else 
         {
            /* if current antidiag doesn't contain next cell, we're at the end of current cloud */
            d = bnd_out.id + bnd_out.rb;   /* antidiag containing current cell (unused) */
            /* if currently building a new cloud in row */
            if ( in_cloud ) /* removed condition => bnd_in.id > d */
            {
               EDGEBOUNDS_Pushback(edg_out, &(bnd_out) );
               in_cloud = false;
            }
         }
      }

      /* if the end of the row is reached and still in cloud */
      if ( in_cloud ) 
      {
         EDGEBOUNDS_Pushback(edg_out, &(bnd_out) );
         in_cloud = false;
      }
   }
}

/*
 *  FUNCTION: EDGEBOUNDS_Reorient_Pushback()
 *  SYNOPSIS: Add antidiag-wise BOUND to row-wise EDGEBOUNDS.
 */
void EDGEBOUNDS_Reorient_Pushback(EDGEBOUNDS*         edg,     /* edgebounds (row-wise, sorted ascending) */
                                  BOUND*              bnd )    /* bound to be inserted (antdiag-wise) */
{
}