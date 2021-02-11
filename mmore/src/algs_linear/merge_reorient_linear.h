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
 *             This is the method selector.
 */
STATUS_FLAG 
EDGEBOUNDS_Union(    const int           Q,             /* query length */
                     const int           T,             /* target length */
                     EDGEBOUNDS*         edg_in_1,      /* edgebounds (fwd, sorted ascending) */
                     EDGEBOUNDS*         edg_in_2,      /* edgebounds (bck, sorted ascending) */
                     EDGEBOUNDS*         edg_out );     /* OUTPUT: union edgebounds (sorted ascending) */

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
                           EDGEBOUNDS*         edg_out );     /* OUTPUT: merged edgebounds (sorted ascending) */

/*! FUNCTION:  EDGEBOUNDS_Union_via_Bridge()
 *  SYNOPSIS:  Combine two edgebound lists into one. 
 *             Bridges all bounds on each row/diag into single bound, spanning from the min to max index on that row/diag.
 *             Assumes input lists are sorted and both have common orientation (byrow/bydiag).
 *             WARNING: Sort is not checked.
 */
STATUS_FLAG 
EDGEBOUNDS_Union_Abridged( const int           Q,             /* query length */
                           const int           T,             /* target length */
                           EDGEBOUNDS*         edg_in_1,      /* edgebounds (fwd, sorted ascending) */
                           EDGEBOUNDS*         edg_in_2,      /* edgebounds (bck, sorted ascending) */
                           EDGEBOUNDS*         edg_out );     /* OUTPUT: merged edgebounds (sorted ascending) */

/*! FUNCTION: EDGEBOUNDS_ReorientToRow()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.       
 */
STATUS_FLAG 
EDGEBOUNDS_ReorientToRow(  const int            Q,             /* query length */
                           const int            T,             /* target length */
                           EDGEBOUNDS*          edg_in,        /* edgebounds (antidiag-wise, sorted ascending) */
                           EDGEBOUND_ROWS*      edg_builder,   /* edgebound working space */
                           EDGEBOUNDS*          edg_out );     /* OUPUT: edgebounds (row-wise, sorted ascending) */

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
                                    EDGEBOUND_ROWS*     edg_rows,       /* temporary working space */
                                    EDGEBOUNDS*         edg_out );      /* OUPUT: edgebounds (row-wise, sorted ascending) */

/*! FUNCTION: EDGEBOUNDS_ReorientToRow_byDiag()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.
 *            For each antidiag, look at each cell and integrates it into a growing <edg_rows>.
 */
STATUS_FLAG 
EDGEBOUNDS_ReorientToRow_byDiag(    const int           Q,          /* query length */
                                    const int           T,          /* target length */
                                    EDGEBOUNDS*         edg_in,     /* edgebounds (antidiag-wise, sorted ascending) */
                                    EDGEBOUND_ROWS*     edg_rows,   /* temporary working space */
                                    EDGEBOUNDS*         edg_out );  /* OUPUT: edgebounds (row-wise, sorted ascending) */

/*! FUNCTION: EDGEBOUNDS_ReorientToRow_byDiff()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.
 *            For each antidiagonal look at the disjunctive union to its previous antidiagonal, 
 *            and only updates those spans at each pass.
 */
STATUS_FLAG 
EDGEBOUNDS_ReorientToRow_byDiff(    const int           Q,          /* query length */
                                    const int           T,          /* target length */
                                    EDGEBOUNDS*         edg_in,     /* edgebounds (antidiag-wise, sorted ascending) */
                                    EDGEBOUND_ROWS*     edg_rows,   /* edgebound working space */
                                    EDGEBOUNDS*         edg_out );  /* OUPUT: edgebounds (row-wise, sorted ascending) */


/*! FUNCTION: EDGEBOUNDS_Reorient_Pushback()
 *  SYNOPSIS: Add antidiag-wise BOUND to row-wise EDGEBOUNDS.
 */
STATUS_FLAG 
EDGEBOUNDS_Reorient_Pushback(   EDGEBOUNDS*         edg,     /* edgebounds (row-wise, sorted ascending) */
                                BOUND*              bnd );   /* bound to be inserted (antdiag-wise) */

#endif /* _MERGE_REORIENT_LINEAR_H */