/*******************************************************************************
 *  FILE:      pruning_linear.c
 *  PURPOSE:   Pruning methods for Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _PRUNING_LINEAR_H
#define _PRUNING_LINEAR_H

/*! FUNCTION: 	PRUNER_diag_by_xdrop_edgetrim()
 *  SYNOPSIS: 	Prunes antidiagonal of Cloud Search.  
 * 				Uses x-drop and only trims in from left and right ends of search space. No bifurcation.
 * 				(1) Updates the total_max, which stores the highest scoring cell in the matrix thus far.
 * 		      	(2) Performs pruning from left-edge, moving right until a cell is found which all states fall below limit = (total_max - alpha).
 * 				(3) Performs pruning from right-edge, moving left.
 * 				(4) Stores edgebounds in left and right bound list.
 */
STATUS_FLAG 
PRUNER_via_xdrop_edgetrim_Linear( 	MATRIX_3D* 		st_MX3,			/* normal state matrix */
									MATRIX_2D* 		sp_MX,			/* special state matrix */
									const float     alpha,			/* x-drop value */
									const int       gamma,			/* number of antidiagonals before pruning */
									const int 		d_1,			/* previous antidiagonal */
									const int 		d_0,			/* current antidiagonal */
									const int 		dx1,			/* previous antidiag (mod-mapped) */
									const int 		dx0, 			/* current antidiag (mod-mapped) */
									const int 		d_cnt, 			/* number of antidiags traversed */
									const int 		le, 			/* right edge of dp matrix on current antidiag */
									const int 		re,				/* left edge of dp matrix on current antidiag */
									float*    		total_max,		/* (UPDATED) current maximum score */
									VECTOR_INT* 	lb_vec[3], 		/* OUTPUT: current list of left-bounds */
									VECTOR_INT* 	rb_vec[3] );	/* OUTPUT: current list of right-bounds */


/*! FUNCTION: 	PRUNER_diag_by_xdrop_bifurcate_Linear()
 *  SYNOPSIS: 	Prunes antidiagonal of Cloud Search.
 * 				Uses x-drop and trims all cells which fall below pruning threshold
 * 				(1) Computes the diag_max on the given antidiagonal.
 * 				(2) Updates the total_max, which stores the highest scoring cell in the matrix thus far.
 * 		      	(3) Performs pruning from left-edge.
 * 				(3a) When cells go from below to above threshold, position added to left-bound list.
 * 				(3b) When a cell go from above to below threshold, position added to right-bound list.
 */
STATUS_FLAG 
PRUNER_via_xdrop_bifurcate_Linear( 	MATRIX_3D* 		st_MX3,			/* normal state matrix */
									MATRIX_2D* 		sp_MX,			/* special state matrix */
									const float     alpha,			/* x-drop value */
									const int       gamma,			/* number of antidiagonals before pruning */
									const int 		d_1,			/* previous antidiagonal */
									const int 		d_0,			/* current antidiagonal */
									const int 		dx1,			/* previous antidiag (mod-mapped) */
									const int 		dx0, 			/* current antidiag (mod-mapped) */
									const int 		d_cnt, 			/* number of antidiags traversed */
									const int 		le, 			/* right edge of dp matrix on current antidiag */
									const int 		re,				/* left edge of dp matrix on current antidiag */
									float*    		total_max,		/* (UPDATED) current maximum score */
									VECTOR_INT* 	lb_vec[3], 		/* OUTPUT: current list of left-bounds */
									VECTOR_INT* 	rb_vec[3] );	/* OUTPUT: current list of right-bounds */

/*! FUNCTION: 	PRUNER_diag_by_xdrop_edgetrim_or_die_Linear()
 *  SYNOPSIS: 	Prunes antidiagonal of Cloud Search.
 * 				Uses x-drop and only trims in from left and right ends of search space. No bifurcation.
 *				Alpha: 		value determining whether cells are pruned in antidiagonal
 * 				Alpha-Max: 	value determining whether search is terminated
 *				Beta:  		number of free passes before pruning
 * 				(1) Updates the total_max, which stores the highest scoring cell in the matrix thus far.
 * 				(1b) If diag_max falls below score global threshold, terminate entire search.
 * 		      	(2) Performs pruning from left-edge, moving right until a cell is found which all states fall below limit = (total_max - alpha).
 * 				(3) Performs pruning from right-edge, moving left.
 * 				(4) Stores edgebounds in left and right bound list.
 */
STATUS_FLAG
PRUNER_via_dbl_xdrop_edgetrim_or_die_Linear( 	MATRIX_3D* 		st_MX3,			/* normal state matrix */
												MATRIX_2D* 		sp_MX,			/* special state matrix */
												const float     alpha,			/* x-drop value for by-diag prune */
												const float 	beta, 			/* x-drop value for global prune */
												const int       gamma,			/* number of antidiagonals before pruning */
												const RANGE 	vit_range, 		/* antidiagonal locations for the start-end of the input viterbi alignment */ 
												const int 		d_1,			/* previous antidiagonal */
												const int 		d_0,			/* current antidiagonal */
												const int 		dx1,			/* previous antidiag (mod-mapped) */
												const int 		dx0, 			/* current antidiag (mod-mapped) */
												const int 		d_cnt, 			/* number of antidiags traversed */
												const int 		le, 			/* right edge of dp matrix on current antidiag */
												const int 		re,				/* left edge of dp matrix on current antidiag */
												float*    		total_max,		/* (UPDATED) current maximum score */
												bool* 			is_term_flag, 	/* (UPDATED) if termination trigger has been reached */ 
												VECTOR_INT* 	lb_vec[3], 		/* OUTPUT: current list of left-bounds */
												VECTOR_INT* 	rb_vec[3] );	/* OUTPUT: current list of right-bounds */
#endif /* _PRUNING_LINEAR_H */