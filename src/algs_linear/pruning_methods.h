/*******************************************************************************
 *  FILE:      pruning_methods.h
 *  PURPOSE:   Assorted Pruning Methods for Cloud Search
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

/* ****************************************************************************************** *
 *  
 *  FUNCTION:  Prune_Xdrop_Trim_Edges
 *  SYNOPSIS:        
 *
 *  ARGS:      <edg>          EDGEBOUNDS to be appended to,
 *             <beta>         Current antidiagonal being pruned,
 *             <>
 *
 *  RETURN:    No return, but <edg> is updated.
 *
/* ****************************************************************************************** */
void Prune_Xdrop_Trim_Edges( EDGEBOUNDS*  edg,
                             int          beta,
                             float        alpha )
{
   // /* UPDATE/PRUNE BOUNDS */
   // /* if free passes are complete (beta < d), prune and set new edgebounds */
   // /* beta must be > 1 */
   // if (beta < d_cnt)
   // {
   //    /* impossible state */
   //    lb_new = INT_MIN;    
   //    rb_new = INT_MIN;

   //    /* FIND MAX SCORE ON CURRENT DIAGONAL */
   //    diag_max = -INF;
   //    for (k = lb_1; k < rb_1; k++)
   //    {
   //       /* coords for quadratic matrix */
   //       i = k;
   //       j = d_1 - i;    /* back one diag */

   //       diag_max = calc_Max( 
   //                      calc_Max( diag_max, MMX3(d1,k) ),
   //                      calc_Max( IMX3(d1,k), DMX3(d1,k) ) );
   //    }

   //    /* Total max records largest cell score seen so far */
   //    total_max = MAX(total_max, diag_max);

   //    /* Set score threshold for pruning */
   //    diag_limit = diag_max - alpha;
   //    total_limit = total_max - alpha;
   //    // printf("total_max: %.2f\t total_limit: %.2f\t diag_max: %.2f\t diag_limit: %.2f\n", total_max, total_limit, diag_max, diag_limit);

   //    /* FIND FIRST SCORE TO EXCEED THRESHOLD FROM THE LEFT */
   //    for (k = lb_1; k < rb_1; k++)
   //    {
   //       /* coords for quadratic matrix */
   //       i = k;
   //       j = d_1 - i; /* looking back one diag */

   //       cell_max = calc_Max( MMX3(d1,k),
   //                      calc_Max( IMX3(d1,k), DMX3(d1,k) ) );

   //       /* prune in left edgebound */
   //       if( cell_max >= total_limit ) {
   //          lb_new = i;
   //          break;
   //       }
   //    }

   //    /* If no boundary edges are found on diag, then branch is pruned entirely and we are done */
   //    if (lb_new == INT_MIN)
   //       break;

   //    /* FIND FIRST SCORE TO EXCEED THRESHOLD FROM THE RIGHT */
   //    for (k = rb_1 - 1; k >= lb_1; k--)
   //    {
   //       /* coords for quadratic matrix */
   //       i = k;
   //       j = d_1 - i;

   //       cell_max = calc_Max( MMX3(d1,k),
   //                      calc_Max( IMX3(d1,k), DMX3(d1,k) ) );

   //       /* prune in right edgebound */
   //       if( cell_max >= total_limit )
   //       {
   //          rb_new = (i + 1);
   //          break;
   //       }
   //    }
   //    // printf("pruned -> (%d,%d)...\n", lb_new, rb_new);
   // }
   // else /* else edges expand in square pattern */
   // {
   //    lb_new = lb;
   //    rb_new = rb;
   // }

   // /* Update bounds */
   // lb = lb_new;
   // rb = rb_new + 1;
   // /* NOTE: FOR TESTING - THiS REMOVES ALL PRUNING */
   // // lb = lb;
   // // rb = rb + 1;

   // /* Edge-checks: find if diag cells that are inside matrix bounds */
   // le = MAX(beg.i, d - T);
   // re = le + num_cells;

   // /* Check that they dont exceed edges of matrix */
   // lb = MAX(lb, le);
   // rb = MIN(rb, re);

   // EDGEBOUNDS_Pushback(edg, (BOUND){d,lb,rb});
}