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

/*
 *  FUNCTION: EDGEBOUNDS_Reflect()
 *  SYNOPSIS: Reflect antidiagonal bounds.
 */
void EDGEBOUNDS_Reflect(	EDGEBOUNDS* 	edg );

/*
 *  FUNCTION: EDGEBOUNDS_Merge()
 *  SYNOPSIS: Combine two edgebound lists into one. Assumes both lists are sorted.
 */
void EDGEBOUNDS_Merge(	const int           Q,             /* query length */
                      	const int           T,             /* target length */
                      	const EDGEBOUNDS*   edg_in_1,      /* edgebounds (fwd, sorted ascending) */
                      	const EDGEBOUNDS*   edg_in_2,      /* edgebounds (bck, sored ascending) */
                      	EDGEBOUNDS*         edg_out );     /* OUTPUT: merged edgebounds */

/*
 *  FUNCTION: EDGEBOUNDS_Reorient_to_Row()
 *  SYNOPSIS: Reorient EDGEBOUNDS from by-diagonal to by-row.       
 */
void EDGEBOUNDS_Reorient_to_Row( const int           Q,          /* query length */
                                 const int           T,          /* target length */
                                 EDGEBOUNDS*         edg_in,     /* edgebounds (antidiag-wise, sorted ascending) */
                                 EDGEBOUNDS*         edg_out );  /* OUPUT: edgebounds (row-wise, sorted ascending) */

/*
 *  FUNCTION: EDGEBOUNDS_Reorient_Pushback()
 *  SYNOPSIS: Add antidiag-wise BOUND to row-wise EDGEBOUNDS.
 */
void EDGEBOUNDS_Reorient_Pushback(	EDGEBOUNDS*  	edg,    		/* edgebounds (row-wise, sorted ascending) */
                                  	BOUND*          bnd );    	/* bound to be inserted (antdiag-wise) */

#endif /* _MERGE_REORIENT_LINEAR_H */