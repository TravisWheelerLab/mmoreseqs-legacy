/*******************************************************************************
 *  FILE:      pruning_quad.c
 *  PURPOSE:   Pruning methods for Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _PRUNING_QUAD_H
#define _PRUNING_QUAD_H

/*
 *  FUNCTION: 	prune_diag_by_xdrop_edgetrim()
 *  SYNOPSIS: 	Prunes antidiagonal of Cloud Search.  
 * 				Uses x-drop and only trims in from left and right ends of search space. No bifurcation.
 * 				(1) Updates the total_max, which stores the highest scoring cell in the matrix thus far.
 * 		      (2) Performs pruning from left-edge, moving right until a cell is found which all states fall below limit = (total_max - alpha).
 * 				(3) Performs pruning from right-edge, moving left.
 * 				(4) Stores edgebounds in left and right bound list.
 */
void prune_via_xdrop_edgetrim_Quad( MATRIX_3D* 		st_MX,			/* normal state matrix */
												MATRIX_2D* 		sp_MX,			/* special state matrix */
	                               	const float     alpha,			/* x-drop value */
	                               	const int       beta,			/* number of antidiagonals before pruning */
												const int 		d_1,			/* previous antidiagonal */
												const int 		d_0,			/* current antidiagonal */
												const int 		d1,				/* previous antidiag (mod-mapped) */
												const int 		d0, 			/* current antidiag (mod-mapped) */
												const int 		d_cnt, 			/* number of antidiags traversed */
												const int 		le, 			/* right edge of dp matrix on current antidiag */
												const int 		re,				/* left edge of dp matrix on current antidiag */
		                              float*    		total_max,		/* (UPDATED) current maximum score */
												VECTOR_INT* 	lb_vec[3], 		/* OUTPUT: current list of left-bounds */
												VECTOR_INT* 	rb_vec[3] );	/* OUTPUT: current list of right-bounds */

#endif /* _PRUNING_QUAD_H */