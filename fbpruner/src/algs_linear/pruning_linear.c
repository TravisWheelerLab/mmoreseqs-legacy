/*******************************************************************************
 *  FILE:      pruning_linear.c
 *  PURPOSE:   Pruning methods for Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
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
#include "pruning_linear.h"

/*
 *  FUNCTION: 	PRUNER_diag_by_xdrop_edgetrim_Linear()
 *  SYNOPSIS: 	Prunes antidiagonal of Cloud Search.
 * 				Uses x-drop and only trims in from left and right ends of search space. No bifurcation.
 * 				(1) Updates the total_max, which stores the highest scoring cell in the matrix thus far.
 * 		      	(2) Performs pruning from left-edge, moving right until a cell is found which all states fall below limit = (total_max - alpha).
 * 				(3) Performs pruning from right-edge, moving left.
 * 				(4) Stores edgebounds in left and right bound list.
 */
inline
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
	                                	VECTOR_INT* 	rb_vec[3] )		/* OUTPUT: current list of right-bounds */
{
	int 		i, j, k; 					/* indexes */
	int 		q_0;						/* row index (query position) */
	int 		t_0;						/* column index (target position) */
	int 		k_0;						/* offset into antidiagonal */
	int 		lb_0, rb_0; 				/* left/right bounds of current antidiagonal */
	int 		lb_1, rb_1; 				/* left/right bounds of previous antidiagonal */
	float 		diag_max, cell_max;         /* max score for all normal states in a given cell/antidiagonal */
	float 		diag_limit, total_limit; 	/* pruning thresholds for current limit */

	/* reset int vectors */
	VECTOR_INT_Reuse( lb_vec[0] );
	VECTOR_INT_Reuse( rb_vec[0] );

	/* Update maximum score using antidiagonal */
	for ( i = 0; i < lb_vec[1]->N; i++ ) {
		lb_1 = VEC_X( lb_vec[1], i );
		rb_1 = VEC_X( rb_vec[1], i );

		diag_max = -INF;
		for ( k_0 = lb_1; k_0 < rb_1; k_0++ )
		{
			q_0 = k_0;
			t_0 = d_1 - k_0;    /* looking back one diag */
			
			diag_max = calc_Max(
			               calc_Max( diag_max, MMX3(dx1, k_0) ),
			               calc_Max( IMX3(dx1, k_0), DMX3(dx1, k_0) ) );
		}

		/* Total max records largest cell score seen so far */
		*total_max = MAX( *total_max, diag_max );
	}

	/* Set score threshold for pruning */
	total_limit = *total_max - alpha;

	/* All edgebounds in previous antidiagonal */
	for ( i = 0; i < lb_vec[1]->N; i++ )
	{
		lb_1 = VEC_X( lb_vec[1], i );
		rb_1 = VEC_X( rb_vec[1], i );

		/* If free passes are not complete, no pruning */
		if ( gamma >= d_cnt )
		{
			VECTOR_INT_Pushback( lb_vec[0], lb_1 );
			VECTOR_INT_Pushback( rb_vec[0], rb_1 );
		}
		else /* If free passes are complete (gamma < d), prune and set new edgebounds */
		{
			/* impossible state */
			lb_0 = INT_MIN;
			rb_0 = INT_MIN;

			/* Find the first cell from the left which passes above threshold */
			for ( q_0 = lb_1; q_0 < rb_1; q_0++ )
			{
				q_0 = k_0;
				t_0 = d_1 - k_0; 	/* looking back one diag */

				cell_max = 	calc_Max( MMX3(dx1, k_0),
				            calc_Max( IMX3(dx1, k_0), DMX3(dx1, k_0) ) );

				/* prune in left edgebound */
				if ( cell_max >= total_limit )
				{
					lb_0 = k_0;
					VECTOR_INT_Pushback( lb_vec[0], lb_0 );
					break;
				}
			}

			/* If no boundary edges are found on diag, then branch is pruned entirely */
			if ( lb_0 == INT_MIN )
				continue;

			/* Find the first cell from the right which passes above threshold */
			for ( k = rb_1 - 1; k >= lb_1; k-- )
			{
				q_0 = k_0;
				t_0 = d_1 - k_0; 	/* looking back one diag */

				cell_max = 	calc_Max( 	MMX3(dx1, k_0),
				            calc_Max(   IMX3(dx1, k_0), DMX3(dx1, k_0 ) ) );

				/* prune in right edgebound */
				if ( cell_max >= total_limit )
				{
					rb_0 = k_0 + 1;
					VECTOR_INT_Pushback( rb_vec[0], rb_0 );
					break;
				}
			}
		}
	}
}

/*
 *  FUNCTION: 	prune_diag_by_xdrop_bifurcate_Linear()
 *  SYNOPSIS: 	Prunes antidiagonal of Cloud Search.
 * 				Uses x-drop and trims all cells which fall below pruning threshold
 * 				(1) Computes the diag_max on the given antidiagonal.
 * 				(2) Updates the total_max, which stores the highest scoring cell in the matrix thus far.
 * 		      	(3) Performs pruning from left-edge.
 * 				(3a) When cells go from below to above threshold, position added to left-bound list.
 * 				(3b) When a cell go from above to below threshold, position added to right-bound list.
 */
inline
STATUS_FLAG 
PRUNER_via_xdrop_bifurcate_Linear( 	MATRIX_3D* 		st_MX3,			/* normal state matrix */
                                    MATRIX_2D* 		sp_MX,			/* special state matrix */
                                    const float    alpha,			/* x-drop value */
                                    const int      gamma,			/* number of antidiagonals before pruning */
                                    const int 		d_1,			   /* previous antidiagonal */
                                    const int 		d_0,			   /* current antidiagonal */
                                    const int 		dx1,			   /* previous antidiag (mod-mapped) */
                                    const int 		dx0, 			   /* current antidiag (mod-mapped) */
                                    const int 		d_cnt, 			/* number of antidiags traversed */
                                    const int 		le, 			   /* right edge of dp matrix on current antidiag */
                                    const int 		re,				/* left edge of dp matrix on current antidiag */
                                    float*    		total_max,		/* (UPDATED) current maximum score */
                                    VECTOR_INT* 	lb_vec[3], 		/* OUTPUT: current list of left-bounds */
                                    VECTOR_INT* 	rb_vec[3] )		/* OUTPUT: current list of right-bounds */
{
	int 		i, j, k; 					   /* indexes */
	int 		q_0;						      /* row index (query position) */
	int 		t_0;						      /* column index (target position) */
	int 		k_0;						      /* offset into antidiagonal */
	int 		lb_0, rb_0; 				   /* left/right bounds of current antidiagonal */
	int 		lb_1, rb_1; 				   /* left/right bounds of previous antidiagonal */
	float 	diag_max, cell_max;        /* max score for all normal states in a given cell/antidiagonal */
	float 	diag_limit, total_limit; 	/* pruning thresholds for current limit */

	/* reset int vectors */
	VECTOR_INT_Reuse( lb_vec[0] );
	VECTOR_INT_Reuse( rb_vec[0] );

	/* Update maximum score using antidiagonal */
	for ( i = 0; i < lb_vec[1]->N; i++ ) {
		lb_1 = lb_vec[1]->data[i];
		rb_1 = rb_vec[1]->data[i];

		diag_max = -INF;
		for ( k_0 = lb_1; k_0 < rb_1; k_0++ )
		{
			q_0 = k_0;
			t_0 = d_1 - k_0;    /* looking back one diag */
			
			diag_max = calc_Max( calc_Max( diag_max,     MMX3(dx1, k) ),
			                     calc_Max( IMX3(dx1, k), DMX3(dx1, k) ) );
		}

		/* Total max records largest cell score seen so far */
		*total_max = MAX( *total_max, diag_max );
	}

	/* Set score threshold for pruning */
	total_limit = *total_max - alpha;

	/* All edgebounds in previous antidiagonal */
	for ( i = 0; i < lb_vec[1]->N; i++ )
	{
		lb_1 = lb_vec[1]->data[i];
		rb_1 = rb_vec[1]->data[i];

		/* If free passes are not complete, do no pruning */
		if ( gamma >= d_cnt )
		{
			VECTOR_INT_Pushback( lb_vec[0], lb_1 );
			VECTOR_INT_Pushback( rb_vec[0], rb_1 );
		}
		else /* If free passes are complete (gamma < d), prune and set new edgebounds */
		{
			/* impossible state */
			lb_0 = INT_MIN;
			rb_0 = INT_MIN;

			/* Find the first cell from the left which passes above threshold */
			for ( k = lb_1; k < rb_1; k++ )
			{
				q_0 = k_0;
				t_0 = d_1 - k_0; 	/* looking back one diag */

				cell_max = 	calc_Max( MMX3(dx1, k_0),
				            calc_Max( IMX3(dx1, k_0), DMX3(dx1, k_0) ) );

				/* prune in left edgebound */
				if ( cell_max >= total_limit )
				{
					lb_0 = k_0;
					VECTOR_INT_Pushback( lb_vec[0], lb_0 );
					break;
				}
			}

			/* If no boundary edges are found on diag, then branch is pruned entirely */
			if ( lb_0 == INT_MIN )
				continue;

			/* Find the first cell from the right which passes above threshold */
			for ( k_0 = rb_1 - 1; k_0 >= lb_1; k_0-- )
			{
				q_0 = k_0;
				t_0 = d_1 - k_0; 	/* looking back one diag */

				cell_max = 	calc_Max( MMX3(dx1, k_0),
				            calc_Max( IMX3(dx1, k_0), DMX3(dx1, k_0) ) );

				/* prune in right edgebound */
				if ( cell_max >= total_limit )
				{
					rb_0 = k_0 + 1;
					VECTOR_INT_Pushback( rb_vec[0], rb_0 );
					break;
				}
			}
		}
	}
}

/*
 *  FUNCTION: 	prune_diag_by_xdrop_edgetrim_or_die_Linear()
 *  SYNOPSIS: 	Prunes antidiagonal of Cloud Search.
 * 				Uses x-drop and only trims in from left and right ends of search space. No bifurcation.
 *				Alpha: 		x-drop value determining whether cells are pruned in antidiagonal
 * 				Alpha-Max: 	x-drop value determining whether search is terminated
 *				Beta:  		number of free passes before pruning
 * 				(1) Updates the total_max, which stores the highest scoring cell in the matrix thus far.
 * 				(1b) If diag_max falls below score global threshold, terminate entire search.
 * 		      (2) Performs pruning from left-edge, moving right until a cell is found which all states fall below limit = (total_max - alpha).
 * 				(3) Performs pruning from right-edge, moving left.
 * 				(4) Stores edgebounds in left and right bound list.
 */
inline
STATUS_FLAG 
PRUNER_via_dbl_xdrop_edgetrim_or_die_Linear( 	MATRIX_3D* 		st_MX3,			/* normal state matrix */
                                                MATRIX_2D* 		sp_MX,			/* special state matrix */
                                                const float    alpha,			/* x-drop value for by-diag prune */
                                                const float 	beta, 			/* x-drop value for global prune */
                                                const int      gamma,			/* number of antidiagonals before pruning */
                                                const RANGE 	vit_range, 		/* antidiagonal locations for the start-end of the input viterbi alignment */ 
                                                const int 		d_1,			   /* previous antidiagonal */
                                                const int 		d_0,			   /* current antidiagonal */
                                                const int 		dx1,			   /* previous antidiag (mod-mapped) */
                                                const int 		dx0, 			   /* current antidiag (mod-mapped) */
                                                const int 		d_cnt, 			/* number of antidiags traversed */
                                                const int 		le, 			   /* right edge of dp matrix on current antidiag */
                                                const int 		re,				/* left edge of dp matrix on current antidiag */
                                                float*    		total_max,		/* (UPDATED) current maximum score */
                                                bool* 			is_term_flag, 	/* (UPDATED) if termination trigger has been reached */ 
                                                VECTOR_INT* 	lb_vec[3], 		/* OUTPUT: current list of left-bounds */
                                                VECTOR_INT* 	rb_vec[3] )		/* OUTPUT: current list of right-bounds */
{
	int 		i, j, k; 					   /* indexes */
	int 		q_0;						      /* row index (query position) */
	int 		t_0;						      /* column index (target position) */
	int 		k_0;						      /* offset into antidiagonal */
	int 		lb_0, rb_0; 				   /* left/right bounds of current antidiagonal */
	int 		lb_1, rb_1; 				   /* left/right bounds of previous antidiagonal */
	float 	diag_max, cell_max;        /* max score for all normal states in a given cell/antidiagonal */
	float 	diag_limit 	= -INF;			/* pruning threshold based on global max */
	float 	total_limit = -INF; 		   /* termination threshold based on antidiag max */
	bool 		is_d_0_in_viterbi; 			/* checks if we have gone past the seed viterbi alignment */

	/* clear data int vectors (which will be used to create edgebounds) */
	VECTOR_INT_Reuse( lb_vec[0] );
	VECTOR_INT_Reuse( rb_vec[0] );

	/* (1) update maximum score using antidiagonal */
	for ( i = 0; i < lb_vec[1]->N; i++ ) {
		lb_1 = VEC_X( lb_vec[1], i );
		rb_1 = VEC_X( rb_vec[1], i );

		diag_max = -INF;
		for ( k_0 = lb_1; k_0 < rb_1; k_0++ )
		{
			q_0 = k_0;
			t_0 = d_1 - k_0;    /* looking back one diag */
			
			diag_max = calc_Max( calc_Max( diag_max,        MMX3(dx1, k_0) ),
			                     calc_Max( IMX3(dx1, k_0),  DMX3(dx1, k_0) ) );
		}

		/* Total max records largest cell score seen so far */
		*total_max = MAX( *total_max, diag_max );
	}

	/* Set score limit for terminating search */
	total_limit = *total_max - beta;
	/* Set score limit threshold for pruning */
	diag_limit 	= diag_max - alpha;

	/* check if outside seed viterbi alignment */
	is_d_0_in_viterbi = IS_IN_RANGE(vit_range.beg, vit_range.end, d_0);
	/* if we have searched the length of the input viterbi alignment, check termination condition */
	if ( is_d_0_in_viterbi == false ) 
	{
		// fprintf(stdout, "# OUT OF RANGE: d_cnt=%d, d_0=%d, vit_range={%d,%d}\n", d_cnt, d_0, vit_range.beg, vit_range.end);
		/* if entire antidiagonal falls below termination threshold (total_limit), then remove all branches and terminate search */
		*is_term_flag = (diag_max < total_limit);
		if ( *is_term_flag == true ) 
		{
			fprintf(stdout, "# FLAG TRIGGERED!!! d_cnt=%d, d_0=%d, vit_range={%d,%d}\n", d_cnt, d_0, vit_range.beg, vit_range.end);
			return;
		}
	}
	
	/* All edgebounds in previous antidiagonal */
	for ( i = 0; i < lb_vec[1]->N; i++ )
	{
		lb_1 = lb_vec[1]->data[i];
		rb_1 = rb_vec[1]->data[i];

		/* If free passes are not complete, skip pruning */
		if ( gamma >= d_cnt )
		{
			VECTOR_INT_Pushback( lb_vec[0], lb_1 );
			VECTOR_INT_Pushback( rb_vec[0], rb_1 );
		}
		else /* If free passes are complete (gamma < d), prune and set new edgebounds */
		{
			/* impossible state */
			lb_0 = INT_MIN;
			rb_0 = INT_MIN;

			/* Find the first cell from the left which passes above threshold */
			for ( k_0 = lb_1; k_0 < rb_1; k_0++ )
			{
				q_0 = k_0;
				t_0 = d_1 - k_0; 	/* looking back one diag */

				cell_max = 	calc_Max( MMX3(dx1, k_0),
				            calc_Max( IMX3(dx1, k_0), DMX3(dx1, k_0) ) );

				/* prune in left edgebound */
				if ( cell_max >= diag_limit )
				{
					/* when first cell found above limit if found, set boundary and end search */
					lb_0 = k_0;
					VECTOR_INT_Pushback( lb_vec[0], lb_0 );
					break;
				}
			}

			/* If no boundary edges are found on diag, then branch is pruned entirely */
			if ( lb_0 == INT_MIN ) {
				continue;
			}

			/* Find the first cell from the right which passes above threshold */
			for ( k_0 = rb_1 - 1; k_0 >= lb_1; k_0-- )
			{
				q_0 = k_0;
				t_0 = d_1 - k_0; 	/* looking back one diag */

				cell_max = 	calc_Max( MMX3(dx1, k_0),
				            calc_Max( IMX3(dx1, k_0), DMX3(dx1, k_0) ) );

				/* prune in right edgebound */
				if ( cell_max >= diag_limit )
				{
					/* when first cell found above limit if found, set boundary and end search */
					rb_0 = k_0 + 1;
					VECTOR_INT_Pushback( rb_vec[0], rb_0 );
					break;
				}
			}
		}
	}
}