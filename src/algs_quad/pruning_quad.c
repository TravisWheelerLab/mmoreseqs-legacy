/*******************************************************************************
 *  FILE:      prune_methods_quad.c
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
#include "structs.h"
#include "utilities.h"
#include "objects.h"
#include "algs_quad.h"

/* header */
#include "pruning_quad.h"

/*
 *  FUNCTION: 	prune_diag_by_xdrop_edgetrim()
 *  SYNOPSIS: 	Prunes antidiagonal of Cloud Search.
 * 				Uses x-drop and only trims in from left and right ends of search space. No bifurcation.
 * 				(1) Updates the total_max, which stores the highest scoring cell in the matrix thus far.
 * 		      	(2) Performs pruning from left-edge, moving right until a cell is found which all states fall below limit = (total_max - alpha).
 * 				(3) Performs pruning from right-edge, moving left.
 * 				(4) Stores edgebounds in left and right bound list.
 */
inline
void prune_via_xdrop_edgetrim_Quad( MATRIX_3D* 		st_MX,			/* normal state matrix */
									MATRIX_2D* 		sp_MX,			/* special state matrix */
									const float    	alpha,			/* x-drop value */
									const int      	beta,			/* number of antidiagonals before pruning */
									const int 		d_1,			/* previous antidiagonal */
									const int 		d_0,			/* current antidiagonal */
									const int 		d1,				/* previous antidiag (mod-mapped) */
									const int 		d0, 			/* current antidiag (mod-mapped) */
									const int 		d_cnt, 			/* number of antidiags traversed */
									const int 		le, 			/* right edge of dp matrix on current antidiag */
									const int 		re,				/* left edge of dp matrix on current antidiag */
									float*    		total_max,		/* (UPDATED) current maximum score */
									VECTOR_INT* 	lb_vec[3], 		/* OUTPUT: current list of left-bounds */
									VECTOR_INT* 	rb_vec[3] )		/* OUTPUT: current list of right-bounds */
{
	int 		b, i, j, k; 				/* indexes */
	int 		lb_0, rb_0; 				/* left/right bounds of current antidiagonal */
	int 		lb_1, rb_1; 				/* left/right bounds of previous antidiagonal */
	float 		diag_max, cell_max;         /* max score for all normal states in a given cell/antidiagonal */
	float 		diag_limit, total_limit; 	/* pruning thresholds for current limit */

	/* reset int vectors */
	VECTOR_INT_Reuse( lb_vec[0] );
	VECTOR_INT_Reuse( rb_vec[0] );

	/* Update maximum score using antidiagonal */
	for ( b = 0; b < lb_vec[1]->N; b++ ) {
		lb_1 = lb_vec[1]->data[b];
		rb_1 = rb_vec[1]->data[b];

		diag_max = -INF;
		for ( k = lb_1; k < rb_1; k++ )
		{
			i = k;
			j = d_1 - i;    /* looking back one diag */
			diag_max = calc_Max(
			               calc_Max( diag_max, MMX(i, j) ),
			               calc_Max( IMX(i, j), DMX(i, j) ) );
		}

		/* Total max records largest cell score seen so far */
		*total_max = MAX( *total_max, diag_max );
	}

	/* Set score threshold for pruning */
	total_limit = *total_max - alpha;

	/* Prune each branch in previous antidiagonal */
	for ( b = 0; b < lb_vec[1]->N; b++ )
	{
		lb_1 = lb_vec[1]->data[b];
		rb_1 = rb_vec[1]->data[b];

		/* If free passes are not complete, do no pruning */
		if ( beta >= d_cnt )
		{
			VECTOR_INT_Pushback( lb_vec[0], lb_1 );
			VECTOR_INT_Pushback( rb_vec[0], rb_1 );
		}
		else /* If free passes are complete (beta < d), prune and set new edgebounds */
		{
			/* impossible state */
			lb_0 = INT_MIN;
			rb_0 = INT_MIN;

			/* Find the first cell from the left which passes above threshold */
			for ( k = lb_1; k < rb_1; k++ )
			{
				i = k;
				j = d_1 - i; 	/* looking back one diag */

				cell_max = 	calc_Max( MMX(i, j),
				                      calc_Max( IMX(i, j), DMX(i, j) ) );

				/* prune in left edgebound */
				if ( cell_max >= total_limit )
				{
					lb_0 = i;
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
				i = k;
				j = d_1 - i; 	/* looking back one diag */

				cell_max = 	calc_Max( MMX(i, j),
				                      calc_Max( IMX(i, j),   DMX(i, j) ) );

				/* prune in right edgebound */
				if ( cell_max >= total_limit )
				{
					rb_0 = i + 1;
					VECTOR_INT_Pushback( rb_vec[0], rb_0 );
					break;
				}
			}
		}
	}
}