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
 */
STATUS_FLAG 
EDGEBOUNDS_Union(    const int           Q,             /* query length */
                     const int           T,             /* target length */
                     EDGEBOUNDS*         edg_in_1,      /* edgebounds (fwd, sorted ascending) */
                     EDGEBOUNDS*         edg_in_2,      /* edgebounds (bck, sorted ascending) */
                     EDGEBOUNDS*         edg_out )      /* OUTPUT: merged edgebounds (sorted ascending) */
{
   EDGEBOUNDS*    edg_in[2];              /* list of input edgebounds */
   int            edg_head[2];            /* head pointers for current edgebound ranges */
   int            edg_minmax[2];          /* min/max antidiag for edg_1 and edg_2 */
   EDGEBOUNDS*    edg;                    /* tmp pointer for edgebounds */
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
   EDGEBOUNDS* edg_cur  = EDGEBOUNDS_Create();

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


/*! FUNCTION: EDGEBOUNDS_Reorient_to_Row()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.
 *            For each row, looks at each viable antidiag and checks if that antidiag intersects the row.
 *            Adds those to growing bound until break is found, then adds it to the row-wise list.   
 */
STATUS_FLAG 
EDGEBOUNDS_Reorient_to_Row(   const int           Q,          /* query length */
                              const int           T,          /* target length */
                              EDGEBOUNDS*         edg_in,     /* edgebounds (antidiag-wise, sorted ascending) */
                              EDGEBOUNDS*         edg_out )   /* OUPUT: edgebounds (row-wise, sorted ascending) */
{
   int         x, y;                              /* indexes into edgebounds */
   RANGE       q_range;                           /* index range for query row */
   RANGE       y_beg, y_end;                      /* index range for edgebound list */
   RANGE       d_range;                           /* index range for valid antidiag ids */
   int         i, j, k;                           /* indexes for diagonal and x-y coords */
   int         id, lb, rb;                        /* diag/row, left-bound, right-bound */
   int         t_0, q_0;                          /* row-wise indexes */
   int         d_0, k_0;                          /* antidiag-wise indexes */
   BOUND       bnd_in, bnd_out;                   /* bounds for */
   bool        in_cloud;                          /* flags whether currently extending bounds along a row */
   bool        is_covered;                        /* whether antidiagonal range intersects row range */
   bool        is_first_row;                      /* whether a edgebound cloud has been found on row yet */
   const int   gap_tolerance  = 0;                /* max distance between two row indexes to merge */

   /* init */

   /* first, verify that input edgebounds are stored by-diagonal */
   if ( edg_in->edg_mode != EDG_DIAG ) {
      /* if not, make a copy of input edgebounds for output */
      EDGEBOUNDS_Copy( edg_out, edg_in );
      return STATUS_SUCCESS;
   }

   /* reuse edgebounds */
   EDGEBOUNDS_Reuse( edg_out, Q, T );
   edg_out->edg_mode = EDG_ROW;

   /* range for possible antidiag indexes */
   q_range.beg = 0;
   q_range.end = Q+1;
   y_beg       = (RANGE){0,0};
   y_end       = (RANGE){0, edg_in->N};

   /* for each row in matrix (position in query) */
   for (q_0 = q_range.beg; q_0 < q_range.end; q_0++)
   {
      /* initialize bound to invalid range */
      bnd_out  = (BOUND) { q_0, -1, -1 };

      /* constrain bounds to valid ranges */
      if (false)
      {
         d_range.beg = 0 + q_0;
         y_beg.beg   = y_beg.end;
         EDGEBOUNDS_NxtRow( edg_in, &y_beg.beg, &y_beg.end, d_0 );
         d_range.end = ((T+1) + q_0) +  1;
         y_end.beg = y_end.end;
         EDGEBOUNDS_NxtRow( edg_in, &y_end.beg, &y_end.end, d_0 );
      }
      
      /* compare against all antidiags which could intersect row */
      for (y = y_beg.beg; y < y_end.end; y++)
      {
         /* get edgebound */
         bnd_in   = EDG_X( edg_in, y );
         d_0      = bnd_in.id;

         /* row-wise coords */
         k_0   = q_0;         /* row index = offset index */
         t_0   = d_0 - k_0;   /* column index where antidiagonal intersects row */

         /* edge-checks (NOTE: these checks could be handled outside the loop if the edgebound ids were indexed) */
         /* if t_0 < 0, it is outside the left-bounds of the dp_matrix */
         if ( t_0 < 0 ) {
            /* if this antidiag is outside the dp_matrix now, it will be outside for all lower rows as well; begin future row's search at this antidiag */
            y_beg.beg = y;
            continue;
         }
         /* if t_0 > T+1, it is outside the right-bounds of the dp_matrix */
         if ( t_0 > T+1 ) {
            /* if this antidiag is outside the dp_matrix now, all subsequent antidiags will be as well; so go to next row */
            break;
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

   return STATUS_SUCCESS;
}

/*! FUNCTION: EDGEBOUNDS_Reorient_to_Row()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.
 *            For each antidiagonal look at the disjunctive union to its previous antidiagonal, 
 *            and only updates those spans at each pass.
 * TODO: WIP!!!
 */
STATUS_FLAG 
EDGEBOUNDS_Reorient_to_Row_via_Diff(   const int           Q,          /* query length */
                                       const int           T,          /* target length */
                                       EDGEBOUNDS*         edg_in,     /* edgebounds (antidiag-wise, sorted ascending) */
                                       EDGEBOUND_ROWS*     edg_rows,   /* OPTIONSAL: temporary working space */
                                       EDGEBOUNDS*         edg_out )   /* OUPUT: edgebounds (row-wise, sorted ascending) */
{
   int         x, y;                               /* indexes into edgebounds */
   RANGE       q_range;                            /* index range for query row */
   RANGE       y_beg, y_end;                       /* index range for edgebound list */
   RANGE       d_range;                            /* index range for valid antidiag ids */
   int         i, j, k;                            /* indexes for diagonal and x-y coords */
   int         id, lb, rb;                         /* diag/row, left-bound, right-bound */
   int         t_0, q_0;                           /* row-wise indexes */
   int         d_0, k_0;                           /* antidiag-wise indexes */
   BOUND       bnd_in, bnd_out;                    /* bounds for */
   bool        in_cloud;                           /* flags whether currently extending bounds along a row */
   bool        is_covered;                         /* whether antidiagonal range intersects row range */
   bool        is_first_row;                       /* whether a edgebound cloud has been found on row yet */
   const int   gap_tolerance  = 0;                 /* max distance between two row indexes to merge */
   RANGE       gaps[(MAX_BOUNDS_PER_ROW * 2) + 1]; /* gaps which indicate the bounds to be updated */ 
   int         gap_cnt;                            /* current number of gaps currently */

   /* init */
   gap_cnt = 0;
   /* if working space is not supplied, create it */
   if ( edg_rows == NULL ) {

   }

   /* first, verify that input edgebounds are stored by-diagonal */
   if ( edg_in->edg_mode != EDG_DIAG ) {
      /* if not, make a copy of input edgebounds for output */
      EDGEBOUNDS_Copy( edg_out, edg_in );
      return STATUS_SUCCESS;
   }

   /* reuse edgebounds */
   EDGEBOUNDS_Reuse( edg_out, Q, T );
   edg_out->edg_mode = EDG_ROW;

   /* range for possible antidiag indexes */
   q_range.beg = 0;
   q_range.end = Q+1;
   y_beg       = (RANGE){0,0};
   y_end       = (RANGE){0, edg_in->N};

   /* for each row in matrix (position in query) */
   for (q_0 = q_range.beg; q_0 < q_range.end; q_0++)
   {
      /* initialize bound to invalid range */
      bnd_out  = (BOUND) { q_0, -1, -1 };

      /* constrain bounds to valid ranges */
      if (false)
      {
         d_range.beg = 0 + q_0;
         y_beg.beg   = y_beg.end;
         EDGEBOUNDS_NxtRow( edg_in, &y_beg.beg, &y_beg.end, d_0 );
         d_range.end = ((T+1) + q_0) +  1;
         y_end.beg = y_end.end;
         EDGEBOUNDS_NxtRow( edg_in, &y_end.beg, &y_end.end, d_0 );
      }
      
      /* compare against all antidiags which could intersect row */
      for (y = y_beg.beg; y < y_end.end; y++)
      {
         /* get edgebound */
         bnd_in   = EDG_X( edg_in, y );
         d_0      = bnd_in.id;

         /* row-wise coords */
         k_0   = q_0;         /* row index = offset index */
         t_0   = d_0 - k_0;   /* column index where antidiagonal intersects row */

         /* edge-checks (NOTE: these checks could be handled outside the loop if the edgebound ids were indexed) */
         /* if t_0 < 0, it is outside the left-bounds of the dp_matrix */
         if ( t_0 < 0 ) {
            /* if this antidiag is outside the dp_matrix now, it will be outside for all lower rows as well; begin future row's search at this antidiag */
            y_beg.beg = y;
            continue;
         }
         /* if t_0 > T+1, it is outside the right-bounds of the dp_matrix */
         if ( t_0 > T+1 ) {
            /* if this antidiag is outside the dp_matrix now, all subsequent antidiags will be as well; so go to next row */
            break;
         }

         /* check if cell is covered by anti-diag... */
         rb = d_0 - bnd_in.lb;
         lb = d_0 - bnd_in.rb + 1;
         /* if intersection point t_0 is inside span */
         is_covered = IS_IN_RANGE( lb, rb, t_0 );
         // printf("y:%d, d:%d => (%d,%d) => is %d in (%d,%d) : %s \n", y, d, i, j, j, bnd_in.lb, bnd_in.rb, is_covered ? "YES" : "NO");
         
         /* if antidiag range covers cell in row... */
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
               EDGEBOUNDS_Pushback(edg_out, &(bnd_out) );
               bnd_out.lb = t_0;
               bnd_out.rb = t_0 + 1;
            }
         } 
      }
      /* if any bound has been found on row, it will need to be closed at the end. */
      if ( is_first_row = (bnd_out.lb < 0), is_first_row == true ) 
      {
         EDGEBOUNDS_Pushback(edg_out, &(bnd_out) );
      }
   }

   return STATUS_SUCCESS;
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