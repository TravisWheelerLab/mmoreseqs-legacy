/*******************************************************************************
 *  FILE:      merge_reorient_linear.c
 *  PURPOSE:   Functions for merging multiple EDGEBOUND objects and 
 *             reorienting from diagonal-wise to row-wise.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/


#ifndef _MERGE_REORIENT_LINEAR_H
#define _MERGE_REORIENT_LINEAR_H

/*! FUNCTION:  EDGEBOUNDS_Reflect()
 *  SYNOPSIS:  Reflect antidiagonal bounds.
 *             Antidiag <d_0> is indexed from query to target (or vice versa).
 *             <edg> must be oriented by-antidiagonal.
 */
STATUS_FLAG 
EDGEBOUNDS_Reflect( EDGEBOUNDS*   edg );

/*! FUNCTION:  EDGEBOUNDS_Merge_Together()
 *  SYNOPSIS:  Combine two edgebound lists into one. 
 *             Assumes input lists are sorted and both oriented by-antidiagonal.
 */
STATUS_FLAG 
EDGEBOUNDS_Union(    const int           Q,             /* query length */
                     const int           T,             /* target length */
                     EDGEBOUNDS*         edg_in_1,      /* edgebounds (fwd, sorted ascending) */
                     EDGEBOUNDS*         edg_in_2,      /* edgebounds (bck, sorted ascending) */
                     EDGEBOUNDS*         edg_out );     /* OUTPUT: union edgebounds (sorted ascending) */

/*! FUNCTION:  EDGEBOUNDS_Union_via_Bridge()
 *  SYNOPSIS:  Combine two edgebound lists into one. Bridges all bounds on each antidiagonal into single bound.
 *             Assumes input lists are sorted and both oriented by-antidiagonal.
 */
STATUS_FLAG 
EDGEBOUNDS_Union_via_Bridge(  const int           Q,             /* query length */
                              const int           T,             /* target length */
                              EDGEBOUNDS*         edg_in_1,      /* edgebounds (fwd, sorted ascending) */
                              EDGEBOUNDS*         edg_in_2,      /* edgebounds (bck, sorted ascending) */
                              EDGEBOUNDS*         edg_out );     /* OUTPUT: merged edgebounds (sorted ascending) */

/*! FUNCTION: EDGEBOUNDS_Reorient_to_Row()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.       
 */
STATUS_FLAG 
EDGEBOUNDS_Reorient_to_Row(   const int           Q,          /* query length */
                              const int           T,          /* target length */
                              EDGEBOUNDS*         edg_in,     /* edgebounds (antidiag-wise, sorted ascending) */
                              EDGEBOUNDS*         edg_out );  /* OUPUT: edgebounds (row-wise, sorted ascending) */

/*! FUNCTION: EDGEBOUNDS_Reorient_to_Row()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.
 *            For each antidiagonal look at the disjunctive union to its previous antidiagonal, 
 *            and only updates those spans at each pass.
 */
STATUS_FLAG 
EDGEBOUNDS_Reorient_to_Row_via_Diff(   const int           Q,          /* query length */
                                       const int           T,          /* target length */
                                       EDGEBOUNDS*         edg_in,     /* edgebounds (antidiag-wise, sorted ascending) */
                                       EDGEBOUND_ROWS*     edg_rows,   /* temporary working space */
                                       EDGEBOUNDS*         edg_out );  /* OUPUT: edgebounds (row-wise, sorted ascending) */


/*! FUNCTION: EDGEBOUNDS_Reorient_Pushback()
 *  SYNOPSIS: Add antidiag-wise BOUND to row-wise EDGEBOUNDS.
 */
STATUS_FLAG 
EDGEBOUNDS_Reorient_Pushback(   EDGEBOUNDS*         edg,     /* edgebounds (row-wise, sorted ascending) */
                                BOUND*              bnd );   /* bound to be inserted (antdiag-wise) */

#endif /* _MERGE_REORIENT_LINEAR_H */