/*******************************************************************************
 *  FILE:      merge_reorient_linear.c
 *  SYNOPSIS:  Functions for merging multiple EDGEBOUND objects and 
 *             reorienting from diagonal-wise to row-wise.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *  TODO:      -  reorient_via_diff implementation: Finds the disjunctive union between antidiag(i) and antidiag(i-1).  
 *                Updates only cells on those rows. Can simplify via "bridging" aka, if each row only has one span, and only need to track the first and last index.
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
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "_algs_linear.h"

/* header */
#include "merge_reorient_linear.h"

/*! FUNCTION:  EDGEBOUNDS_Reflect()
 *  SYNOPSIS:  Reflect antidiagonal bounds.
 *             Antidiag <d_0> is indexed from query to target (or vice versa).
 *             <edg> must be oriented by-antidiagonal.
 */
STATUS_FLAG 
EDGEBOUNDS_Reflect( EDGEBOUNDS*   edg ) 
{
   int d_0, lb, rb, rb_new, lb_new;

   for (int i = 0; i < edg->N; i++) {
      d_0   = edg->bounds[i].id;
      lb    = edg->bounds[i].lb;
      rb    = edg->bounds[i].rb;

      rb = d_0 - lb + 1;
      lb = d_0 - rb;
      edg->bounds[i] = (BOUND){d_0, lb, rb};
   }

   return STATUS_SUCCESS;
}

/*! FUNCTION:  EDGEBOUNDS_Union()
 *  SYNOPSIS:  Combine two edgebound lists into one to cover the union.
 *             Assumes input lists are sorted ascending.
 *             This is the method selector.
 */
STATUS_FLAG 
EDGEBOUNDS_Union(    const int           Q,             /* query length */
                     const int           T,             /* target length */
                     EDGEBOUNDS*         edg_in_1,      /* edgebounds (fwd, sorted ascending) */
                     EDGEBOUNDS*         edg_in_2,      /* edgebounds (bck, sorted ascending) */
                     EDGEBOUNDS*         edg_out )      /* OUTPUT: merged edgebounds (sorted ascending) */
{
   STATUS_FLAG status;
   status = EDGEBOUNDS_Union_byRow( Q, T, edg_in_1, edg_in_2, edg_out );
   return status;
}


/*! FUNCTION:  EDGEBOUNDS_Union_byRow()
 *  SYNOPSIS:  Combine two edgebound lists into one to cover the union.
 *             Assumes input lists are sorted ascending.
 *  METHOD:    Going by row/diag <id_0>, by accumulating a sublist of all bounds from <edg_in_1> and <edg_in_2> in <id_0>.
 *             Then loops through all bounds in sublist pairwise, merging down overlapping edges. Continues until no more merges can be performed.
 *             Pushes the sublist onto <edg_out> and continues as such through each row/diag.
 */
STATUS_FLAG 
EDGEBOUNDS_Union_byRow(    const int           Q,             /* query length */
                           const int           T,             /* target length */
                           EDGEBOUNDS*         edg_in_1,      /* edgebounds (fwd, sorted ascending) */
                           EDGEBOUNDS*         edg_in_2,      /* edgebounds (bck, sorted ascending) */
                           EDGEBOUNDS*         edg_out )      /* OUTPUT: merged edgebounds (sorted ascending) */
{
   EDGEBOUNDS*    edg_in[2];              /* list of input edgebounds */
   int            edg_head[2];            /* head pointers for current edgebound ranges */
   int            edg_minmax[2];          /* min/max antidiag for edg_1 and edg_2 */
   EDGEBOUNDS*    edg;                    /* tmp pointer for edgebounds */
   EDGEBOUNDS*    edg_cur;                /* working space */
   int            i, j;                   /* indexes */
   RANGE          d_range;                /* antidiag range of full matrix */
   int            r_0, r_0b, r_0e;        /* antidiag range of edgebounds */
   int            id_0, lb, rb;           /* edgebound data */
   /* constants */
   const int      tol         = 0;        /* if clouds are within tolerance range, clouds are merged */
   const int      num_input   = 2;        /* number of input edgebounds (current locked to two) */

   /* list of edgebound sets */
   edg_head[0] = 0;
   edg_in[0]   = edg_in_1; 
   edg_head[1] = 0;
   edg_in[1]   = edg_in_2;
   /* find ranges */
   edg_minmax[0] = EDG_X( edg_in[0], 0 ).id;
   edg_minmax[1] = EDG_X( edg_in[1], 0 ).id;
   d_range.beg = MIN( edg_minmax[0], edg_minmax[1] );
   edg_minmax[0] = EDG_X( edg_in[0], edg_in[0]->N - 1 ).id;
   edg_minmax[1] = EDG_X( edg_in[1], edg_in[1]->N - 1 ).id;
   d_range.end = MAX( edg_minmax[0], edg_minmax[1] );

   /* reset output edgebounds */
   EDGEBOUNDS_Reuse( edg_out, Q, T );
   
   /* verify that all input edgebounds are the same mode */
   if ( edg_in[0]->edg_mode != edg_in[1]->edg_mode ) 
   {
      fprintf( stderr, "ERROR: Not all edgebounds being merged have same orientation!");
      exit(EXIT_FAILURE);
   }
   edg_out->edg_mode = edg_in_1->edg_mode;

   /* working space: edgebounds which the two old edgebounds will be merged into (by row) */
   edg_cur  = EDGEBOUNDS_Create();

   /* init */
   r_0b = r_0e = 0;

   /* iterate over all antidiag/row id */
   for (id_0 = d_range.beg; id_0 < d_range.end; id_0++) 
   {
      /* iterate over input edgebound lists and find all bounds on current diagonal */
      for (int i = 0; i < num_input; i++) 
      {
         /* load edgebounds */
         edg   = edg_in[i];
         r_0b  = edg_head[i];
         r_0e  = edg_head[i];

         /* find range of edgebounds in current antidiag */
         EDGEBOUNDS_NxtRow( edg, &r_0b, &r_0e, id_0 );

         /* add all edgebounds to merge list */
         for (r_0 = r_0b; r_0 < r_0e; r_0++) {
            BOUND* bnd = &EDG_X( edg, r_0 );
            EDGEBOUNDS_Pushback( edg_cur, bnd );
         }
         /* update next starting location to the end location of the last row */
         edg_head[i] = r_0e;
      }

      /* sort output edgebounds, since inter-row ordering could be lost */
      EDGEBOUNDS_Sort( edg_cur );

      /* merge down current list as much as possible (loop until no merges occur) */
      bool has_merged = true;
      while ( has_merged == true )
      {
         has_merged = false;
         BOUND bnd_1, bnd_2;
         /* check if merge possible for all pairs (n^2 for number of bounds in antidiag) */
         for (int i = 0; i < edg_cur->N; i++) 
         {
            bnd_1 = EDG_X( edg_cur, i );
            
            for (int j = i+1; j < edg_cur->N; j++) 
            {
               bnd_2 = EDG_X( edg_cur, j );

               /* if bounds overlap, then merge them (NOTE: logic can probably be optimized) */
               bool overlap = IS_IN_RANGE(bnd_1.lb, bnd_1.rb, bnd_2.lb) ||
                              IS_IN_RANGE(bnd_1.lb, bnd_1.rb, bnd_2.rb);
               if ( overlap == true ) 
               {
                  /* take mins and maxes to build new range */
                  lb = MIN(bnd_1.lb, bnd_2.lb);
                  rb = MAX(bnd_1.rb, bnd_2.rb);
                  BOUND bnd = (BOUND){ id_0, lb, rb };
                  
                  /* replace the edgebound at first index */
                  EDGEBOUNDS_Insert(edg_cur, i, &bnd );
                  /* delete the edgebound at the second index (backfills from end of list) */
                  EDGEBOUNDS_Delete(edg_cur, j);
                  
                  /* flag that a merge has occurred */
                  has_merged = true;
               }
            }
         }
      }

      /* insert all in current list into the output list */
      for (int i = 0; i < edg_cur->N; i++) {
         BOUND* bnd = &(edg_cur->bounds[i]);
         EDGEBOUNDS_Pushback( edg_out, bnd );
      }

      /* empty current list */
      EDGEBOUNDS_Clear( edg_cur );
   }

   /* free working space */
   EDGEBOUNDS_Destroy( edg_cur );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  EDGEBOUNDS_Union_via_Bridge()
 *  SYNOPSIS:  Combine two edgebound lists into one. Bridges all bounds on each antidiagonal into single bound.
 *             Assumes input lists are sorted and both oriented by-antidiagonal.
 */
STATUS_FLAG 
EDGEBOUNDS_Union_Abridged(  const int           Q,             /* query length */
                           const int           T,             /* target length */
                           EDGEBOUNDS*         edg_in_1,      /* edgebounds (fwd, sorted ascending) */
                           EDGEBOUNDS*         edg_in_2,      /* edgebounds (bck, sorted ascending) */
                           EDGEBOUNDS*         edg_out )      /* OUTPUT: merged edgebounds (sorted ascending) */
{
   EDGEBOUNDS*    edg_in[2];              /* list of input edgebounds */
   int            edg_head[2];            /* head pointers for current edgebound ranges */
   int            edg_minmax[2];          /* min/max antidiag for edg_1 and edg_2 */
   EDGEBOUNDS*    edg;                    /* tmp pointer for edgebounds */
   BOUND          bnd_abridged;            /* bound for storing bridged bound */
   int            i, j;                   /* indexes */
   RANGE          d_range;                /* antidiag range of full matrix */
   int            r_0, r_0b, r_0e;        /* antidiag range of edgebounds */
   int            id_0, lb, rb;           /* edgebound data */
   int            bnd_cnt;                /* number of bounds of current antidiag */
   /* constants */
   const int      tol         = 0;        /* if clouds are within tolerance range, clouds are merged */
   const int      num_input   = 2;        /* number of input edgebounds (current locked to two) */

   /* list of edgebound sets */
   edg_head[0] = 0;
   edg_in[0]   = edg_in_1; 
   edg_head[1] = 0;
   edg_in[1]   = edg_in_2;
   /* find ranges */
   edg_minmax[0] = EDG_X( edg_in[0], 0 ).id;
   edg_minmax[1] = EDG_X( edg_in[1], 0 ).id;
   d_range.beg = MIN( edg_minmax[0], edg_minmax[1] );
   edg_minmax[0] = EDG_X( edg_in[0], edg_in[0]->N - 1 ).id;
   edg_minmax[1] = EDG_X( edg_in[1], edg_in[1]->N - 1 ).id;
   d_range.end = MAX( edg_minmax[0], edg_minmax[1] );

   /* reset output edgebounds */
   EDGEBOUNDS_Reuse( edg_out, Q, T );
   
   /* verify that all input edgebounds are the same mode */
   if ( edg_in[0]->edg_mode != edg_in[1]->edg_mode ) 
   {
      fprintf( stderr, "ERROR: Not all edgebounds being merged have same orientation!");
      exit(EXIT_FAILURE);
   }
   edg_out->edg_mode = edg_in_1->edg_mode;

   /* init */
   r_0b = r_0e = 0;
   /* iterate over all antidiag/row id */
   for (id_0 = d_range.beg; id_0 < d_range.end; id_0++) 
   {
      /* start bound with impossible bounds */
      bnd_abridged = (BOUND) { id_0, INT_MAX, INT_MIN };
      bnd_cnt = 0;
      /* iterate over input edgebound lists and find all bounds on current diagonal */
      for (int i = 0; i < num_input; i++) 
      {
         /* load edgebounds */
         edg   = edg_in[i];
         r_0b  = edg_head[i];
         r_0e  = edg_head[i];

         /* find range of edgebounds in current antidiag */
         EDGEBOUNDS_NxtRow( edg, &r_0b, &r_0e, id_0 );
         bnd_cnt += (r_0e - r_0b);

         /* account for all edgebounds to bridged bound */
         for (r_0 = r_0b; r_0 < r_0e; r_0++) {
            BOUND bnd   = EDG_X( edg, r_0 );
            /* take the minimum and maximum of all bounds in antidiag */
            bnd_abridged.lb  = MIN( bnd_abridged.lb, bnd.lb );
            bnd_abridged.rb  = MAX( bnd_abridged.rb, bnd.rb );
         }
         /* update next starting location to the end location of the last row */
         edg_head[i] = r_0e;
      }

      /* add bound to new list */
      EDGEBOUNDS_Pushback( edg_out, &bnd_abridged );
   }

   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_ReorientToRow()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.
 *            This is a method selector. 
 */
STATUS_FLAG 
EDGEBOUNDS_ReorientToRow(     const int           Q,              /* query length */
                              const int           T,              /* target length */
                              EDGEBOUNDS*         edg_in,         /* edgebounds (antidiag-wise, sorted ascending) */
                              EDGEROWS*     edg_builder,    /* temporary working space */
                              EDGEBOUNDS*         edg_out )       /* OUPUT: edgebounds (row-wise, sorted ascending) */
{
   STATUS_FLAG status;
   status = EDGEBOUNDS_ReorientToRow_byDiag( Q, T, edg_in, edg_builder, edg_out );
   return status;
}

/*! FUNCTION: EDGEBOUNDS_ReorientToRow_byRow()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.
 *            Approaches the problem row-wise.
 *            For each row, looks at each viable antidiag and checks if that antidiag intersects the row.
 *            Adds those to growing bound until break is found, then adds it to the row-wise list.   
 */
STATUS_FLAG 
EDGEBOUNDS_ReorientToRow_byRow(     const int           Q,              /* query length */
                                    const int           T,              /* target length */
                                    EDGEBOUNDS*         edg_in,         /* edgebounds (antidiag-wise, sorted ascending) */
                                    EDGEROWS*     edg_rows,       /* temporary working space */
                                    EDGEBOUNDS*         edg_out )       /* OUPUT: edgebounds (row-wise, sorted ascending) */
{
   int         x, y;                              /* indexes into edgebounds */
   RANGE       q_range;                           /* index range for query row */
   RANGE       y_beg, y_end;                      /* index range for edgebound list */
   int         y_beg_alt, y_end_alt;              /* alt method for getting index range for edgebound list */
   bool        y_end_set;                         /* verifies <y_end> only updated once per loop */
   int         y_min, y_max;                      /* minimum and maximum possible values for y */
   RANGE       d_range;                           /* index range for valid antidiag ids */
   int         i, j, k;                           /* indexes for diagonal and x-y coords */
   int         id, lb, rb;                        /* diag/row, left-bound, right-bound */
   int         t_0, q_0;                          /* row-wise indexes */
   int         d_0, k_0;                          /* antidiag-wise indexes */
   BOUND       bnd_in, bnd_out;                   /* bounds for */
   bool        in_cloud;                          /* flags whether currently extending bounds along a row */
   bool        is_covered;                        /* whether antidiagonal range intersects row range */
   bool        is_first_row;                      /* whether a edgebound cloud has been found on row yet */
   int         gap_count;                         /* number of gaps between ranges in antidiags */
   int         cell_count;                        /* counts number of times cell was touched */ 
   const int   gap_tolerance  = 0;                /* max distance between two row indexes to merge */
   const bool  outside_loop_method = true;        /* determines method to compute y_range */

   cell_count  = 0;

   /* first, verify that input edgebounds are stored by-diagonal */
   if ( edg_in->edg_mode == EDG_ROW ) {
      /* if already stored by row, just make a copy of input edgebounds for output */
      EDGEBOUNDS_Copy( edg_out, edg_in );
      return STATUS_SUCCESS;
   }

   /* reuse edgebounds */
   EDGEBOUNDS_Reuse( edg_out, Q, T );
   edg_out->edg_mode = EDG_ROW;

   /* range for possible antidiag indexes */
   q_range.beg = 0;
   q_range.end = Q+1;
   /* range in edgebound list covering <q_range> */
   y_min       = 0;
   y_max       = edg_in->N;
   y_beg       = (RANGE){0, 0};
   y_end       = (RANGE){0, 0};
   y_beg_alt   = 0;
   y_end_alt   = edg_in->N;

   /* for each row in matrix (position in query) */
   for (q_0 = q_range.beg; q_0 < q_range.end; q_0++)
   {
      y_end_set = false;

      /* initialize bound to invalid range */
      bnd_out  = (BOUND){ q_0, -1, -1 };

      /* this outside loop method finds the valid start and end points of loop before it starts, so no innerloop if statements */
      if ( outside_loop_method == true )
      {
         /* find lowest value <d_0> such that its <t_0> intersection of <q_0>th row of matrix is inside matrix bounds */
         d_range.beg = 0 + q_0;
         /* find greatest value <d_0> such that its <t_0> intersection of <q_0>th row of matrix is inside matrix bounds */
         d_range.end = ((T+1) + q_0);
         /* given start and end points of row, we can increment <y> to cover start of valid edgebounds in list */
         y_beg.beg   = y_beg.end;
         EDGEBOUNDS_NxtRow( edg_in, &y_beg.beg, &y_beg.end, d_range.beg );
         /* given start and end points of row, we can increment <y> to cover start of valid edgebounds in list */
         y_end.beg   = y_end.end;
         EDGEBOUNDS_NxtRow( edg_in, &y_end.beg, &y_end.end, d_range.end );
         y_min = y_beg.beg;
         y_max = y_end.end;
      }

      /* this inner loop method finds the valid start and end points by checking at every increment of y */
      if ( outside_loop_method == false )
      {
         y_min = y_beg_alt;
      }
      
      /* compare against all antidiags which could intersect row */
      for (y = y_min; y < y_max; y++)
      {
         /* get edgebound */
         bnd_in   = EDG_X( edg_in, y );
         d_0      = bnd_in.id;

         /* row-wise coords */
         k_0   = q_0;         /* row index = offset index */
         t_0   = d_0 - k_0;   /* column index where antidiagonal intersects row */
         cell_count++;

         /* edge-checks (NOTE: removed b/c checks now handled outside the loop) */
         if ( outside_loop_method == false )
         {
            /* if t_0 < 0, it is outside the left-bounds of the dp_matrix */
            if ( t_0 < 0 ) {
               /* if this antidiag is outside the dp_matrix now, it will be outside for all lower rows as well; begin future row's search at this antidiag */
               y_beg_alt = y;
               continue;
            }
            /* if t_0 > T+1, it is outside the right-bounds of the dp_matrix */
            if ( t_0 > T+1 ) {
               /* if this antidiag is outside the dp_matrix now, all subsequent antidiags will be as well; so go to next row */
               if ( y_end_set == false ) {
                  y_end_set = true;
                  y_end_alt = y;
               }
               break;
            }
         }
         
         /* check if cell is covered by anti-diag... */
         rb = d_0 - bnd_in.lb;
         lb = d_0 - bnd_in.rb + 1;

         /* if antidiag's intersection with point t_0 is inside span */
         is_covered = IS_IN_RANGE( lb, rb, t_0 );
         /* if antidiag range covers cell... */
         if ( is_covered == true ) 
         {
            /* if is adjacent (or with gap tolerance) to current cloud, extend bounds */
            if ( t_0 <= bnd_out.rb + gap_tolerance ) {
               bnd_out.rb = t_0 + 1;
            }
            /* if this is not the first bound on row, push previous bound onto list, and start new */
            elif ( is_first_row = (bnd_out.lb < 0), is_first_row == true ) {
               bnd_out.lb = t_0;
               bnd_out.rb = t_0 + 1;
            }
            /* if this is not the first bound on row, close current bound */ 
            else {
               EDGEBOUNDS_Pushback(edg_out, &bnd_out);
               bnd_out.lb = t_0;
               bnd_out.rb = t_0 + 1;
            }
         } 
      }
      /* if any bound has been found on row, it will need to be closed at the end. */
      if ( is_first_row = (bnd_out.lb < 0), is_first_row == false ) 
      {
         EDGEBOUNDS_Pushback(edg_out, &bnd_out);
      }
   }

   printf("CELLS TOUCHED: %d %d %f\n", cell_count, Q*T, (float)cell_count/(float)(Q*T) );
   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_ReorientToRow_byDiag()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.
 *            For each antidiag, look at each cell and integrates it into a growing <edg_rows>.
 */
STATUS_FLAG 
EDGEBOUNDS_ReorientToRow_byDiag(    const int           Q,          /* query length */
                                    const int           T,          /* target length */
                                    EDGEBOUNDS*         edg_in,     /* edgebounds (antidiag-wise, sorted ascending) */
                                    EDGEROWS*     edg_rows,   /* temporary working space */
                                    EDGEBOUNDS*         edg_out )   /* OUPUT: edgebounds (row-wise, sorted ascending) */
{
   int         x, y;                              /* indexes into edgebounds */
   RANGE       q_range;                           /* index range for query row */
   RANGE       y_beg, y_end;                      /* index range for edgebound list */
   int         y_beg_alt, y_end_alt;              /* alt method for getting index range for edgebound list */
   bool        y_end_set;                         /* verifies <y_end> only updated once per loop */
   int         y_min, y_max;                      /* minimum and maximum possible values for y */
   RANGE       d_range;                           /* index range for valid antidiag ids */
   int         i, j, k;                           /* indexes for diagonal and x-y coords */
   int         id, lb, rb;                        /* diag/row, left-bound, right-bound */
   int         t_0, q_0;                          /* row-wise indexes */
   int         d_0, k_0;                          /* antidiag-wise indexes */
   int         k_min, k_max;                      /* minimum and maximum k_0 indexes reached by antidiag */
   int         lb_0, rb_0;                        /* left and right bounds */
   BOUND*      bnd_in;                            /* bound from <edg_in> */
   BOUND*      bnd_rows;                          /* bound from <edg_rows> */
   BOUND*      bnd_out;                           /* bound from <edg_rows> */
   BOUND       bnd_new;                           /* space for creating new bound */
   bool        in_cloud;                          /* flags whether currently extending bounds along a row */
   bool        is_covered;                        /* whether antidiagonal range intersects row range */
   bool        is_first_row;                      /* whether a edgebound cloud has been found on row yet */
   int         gap_count;                         /* number of gaps between ranges in antidiags */
   int         cell_count;                        /* counts number of times cell was touched */ 
   const int   gap_tolerance  = 0;                /* max distance between two row indexes to merge */
   const bool  outside_loop_method = true;        /* determines method to compute y_range */
   const bool  abridged_method = true;            /* determines whether to allow gaps on row bounds */

   cell_count  = 0;

   /* first, verify that input edgebounds are stored by-diagonal */
   if ( edg_in->edg_mode == EDG_ROW ) {
      /* if already stored by row, just make a copy of input edgebounds for output */
      EDGEBOUNDS_Copy( edg_out, edg_in );
      return STATUS_SUCCESS;
   }

   /* range of <edg_in> edgebound list */
   y_min = 0;
   y_max = edg_in->N;

   /* find the minimum and maximum rows reached by edgebounds */
   q_range.beg = Q+1;
   q_range.end = 0;
   for (y = y_min; y < y_max; y++) 
   {
      bnd_in   = EDGEBOUNDS_Get( edg_in, y );
      d_0      = bnd_in->id;
      k_min    = bnd_in->lb;
      k_max    = bnd_in->rb - 1;
      /* check if antidiag bound has the lowest row rank */
      q_0 = k_min;
      t_0 = d_0 - q_0;
      q_range.beg = MIN( q_range.beg, q_0 );
      /* check if antidiag bound has the lowest row rank */
      q_0 = k_max;
      t_0 = d_0 - q_0;
      q_range.end = MAX( q_range.end, q_0 + 1 );
   }

   /* set edgebound builder only to cover the rows that cloud touches */
   printf("Q,T=(%d,%d), Q_range=(%d,%d)\n", Q, T, q_range.beg, q_range.end );
   EDGEROWS_Reuse( edg_rows, Q, T, q_range );

   /* reuse edgebounds */
   EDGEBOUNDS_Reuse( edg_out, Q, T );
   edg_out->edg_mode = EDG_ROW;
   
   /* for every bound in edgebounds list */
   for (y = y_min; y < y_max; y++)
   {
      bnd_in   = EDGEBOUNDS_Get( edg_in, y );
      d_0      = bnd_in->id;
      lb_0     = bnd_in->lb;
      rb_0     = bnd_in->rb;

      /* for every cell that is spanned by antidiag, integrate into row-wise list */
      for (k_0 = lb_0; k_0 < rb_0; k_0++)
      {
         /* row-wise coords */
         q_0 = k_0;
         t_0 = d_0 - q_0;
         cell_count++;

         /* get latest bound in requested row */
         bnd_rows = EDGEROWS_GetLast_byRow( edg_rows, q_0 );

         /* determine whether to expand current row bounds or create new */
         bool is_expand_row;

         /* non-abrided method will create a new edgebound range for each continuous span on row */
         if ( abridged_method == false ) {
            /* if row is not empty AND row bounds is adjacent (within tolerance) to new cell, merge them */
            is_expand_row = ( (bnd_rows != NULL) && (t_0 <= bnd_rows->rb + gap_tolerance) );
         }
         /* abridged method spans the entire row from the minimum to maximum of cloud. Same as infinite gap tolerance */
         if ( abridged_method == true ) {
            /* if row is not empty AND row bounds is adjacent (within tolerance) to new cell, merge them */
            is_expand_row = ( bnd_rows != NULL );
         }
         
         if ( is_expand_row == true ) 
         {
            bnd_rows->rb = t_0 + 1;
         }
         /* otherwise, create new bound and to row */
         else
         {
            bnd_new = (BOUND){ q_0, t_0, t_0+1};
            EDGEROWS_Pushback( edg_rows, q_0, &bnd_new );
         }
      }
   }

   printf("CELLS TOUCHED: %d %d %f\n", cell_count, Q*T, (float)cell_count/(float)(Q*T) );
   EDGEROWS_Convert( edg_rows, edg_out );
}

/*! FUNCTION: EDGEBOUNDS_ReorientToRow_byDiff()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.
 *            For each antidiagonal look at the disjunctive union to its previous antidiagonal, 
 *            and only updates those spans at each pass.
 * TODO: WIP!!!
 */
STATUS_FLAG 
EDGEBOUNDS_ReorientToRow_byDiff(    const int           Q,          /* query length */
                                    const int           T,          /* target length */
                                    EDGEBOUNDS*         edg_in,     /* edgebounds (antidiag-wise, sorted ascending) */
                                    EDGEROWS*     edg_rows,   /* temporary working space */
                                    EDGEBOUNDS*         edg_out )   /* OUPUT: edgebounds (row-wise, sorted ascending) */
{
   
}

/*! FUNCTION: EDGEBOUNDS_Reorient_Pushback()
 *  SYNOPSIS: Add antidiag-wise BOUND to row-wise EDGEBOUNDS.
 */
STATUS_FLAG 
EDGEBOUNDS_Reorient_Pushback(    EDGEBOUNDS*         edg,     /* edgebounds (row-wise, sorted ascending) */
                                 BOUND*              bnd )    /* bound to be inserted (antdiag-wise) */
{

   return STATUS_SUCCESS;
}