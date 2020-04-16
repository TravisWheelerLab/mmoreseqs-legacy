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

/* */
void edgebounds_Reflect(EDGEBOUNDS *edg);

/* merge two edgebounds into a single edgebounds */
void EDGEBOUNDS_Merge(const int         	Q, 
                      const int         	T,
                      const EDGEBOUNDS* 	edg_in_1,
                      const EDGEBOUNDS* 	edg_in_2, 
                      EDGEBOUNDS* 			edg_out );

/* reorient edgebounds from by-diag to by-row */
void EDGEBOUNDS_Reorient(const int         	Q, 
                         const int         	T,
                         const EDGEBOUNDS* 	edg_in,
                         EDGEBOUNDS* 		edg_out );

#endif /* _MERGE_REORIENT_LINEAR_H */