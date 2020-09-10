/*******************************************************************************
 *  FILE:      pipeline_main.c
 *  PURPOSE:   Main Cloud Search Pipeline.
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
#include <ctype.h>
#include <time.h>

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "objects.h"
#include "parsers.h"
#include "algs_linear.h"
#include "algs_quad.h"
#include "algs_naive.h"

/* header */
#include "pipelines.h"

/* standard pipeline */
void main_pipeline( WORKER* worker ) 
{
	/* initialize data structures needed for tasks */
	WORK_init( worker );
	WORK_open( worker );

	/* worker objects */
	ARGS* 		args 		= worker->args;
	TASKS* 		tasks 		= worker->tasks;
	
	TIMES* 		times 		= worker->times;
	NAT_SCORES* scores 		= worker->scores;

	RESULTS* 	results 	= worker->results;
	RESULT* 	result 		= worker->result;

	/* set flags for pipeline tasks */
	/* TASKS */
	{
		/* linear algs */
		tasks->linear 			= true;		/* if any other linear tasks are flagged, this must be too */
		tasks->lin_fwd 			= false;	/* optional, but can't recover alignment */
		tasks->lin_bck 			= false;	/* optional, but can't recover alignment */
		tasks->lin_vit 			= false;	/* optional, but can't recover alignment */
		tasks->lin_trace 		= false;	/* optional, but can't recover alignment */
		tasks->lin_bound_fwd 	= true;
		tasks->lin_bound_bck 	= true;
		/* quadratic algs */
		tasks->quadratic 		= true;		/* if any other quadratic tasks are flagged, this must be too */
		tasks->quad_fwd 		= true;		/* optional */
		tasks->quad_bck 		= false;	/* optional */
		tasks->quad_vit 		= true;		/* viterbi required for cloud search */
		tasks->quad_trace 		= true;		/* traceback required for cloud search  */
		tasks->quad_bound_fwd 	= false;	/* required step of cloud search */
		tasks->quad_bound_bck 	= false;	/* optional */
	}

	/* load or build, then sort target and query index files */
	WORK_index( worker );

	/* set and verify index ranges */
	WORK_set_ranges( worker );

	/* report ranges */
	printf_vlo("# num_targets: %d, num_queries: %d\n", 
		worker->t_index->N, worker->q_index->N );
	printf_vlo("# target_range: (%d,%d), query_range: (%d,%d)\n", 
		args->t_range.beg, args->t_range.end,
		args->q_range.beg, args->q_range.end );
	
	/* get loop ranges */
	int i_beg, i_end, j_beg, j_end;
	int search_cnt, search_tot;
	i_beg 	= args->t_range.beg;
	i_end 	= args->t_range.end;
	j_beg 	= args->q_range.beg;
	j_end 	= args->q_range.end;
	search_cnt 	= 0;
	search_tot 	= (i_end - i_beg) * (j_end - j_beg);

	WORK_report_header( worker );

	/* loop over target range */
	for (int i = i_beg; i < i_end; i++) 
	{
		/* load in next target */
		WORK_load_target_by_id( worker, i );
		result->target_id = i;
		// HMM_PROFILE_Dump( worker->t_prof, stdout );

		/* loop over query range */
		for (int j = j_beg; j < j_end; j++) 
		{
			printf_vlo("# Running search (%d of %d) for target (%d of %d) and query (%d of %d)...\n",
				search_cnt+1, search_tot, i+1, i_end, j+1, j_end);

			/* load in next query */
			WORK_load_query_by_id( worker, j );
			result->query_id = j;

			/* resize dp matrices to new query/target */
			WORK_reuse( worker );

			const int pad = 5;
			printf_vall("# (%*d) | T_ID: (%*d/%*d) | Q_ID: (%*d/%*d) | T_LEN: %*d | Q_LEN: %*d \n", 
				pad, search_cnt, pad, i, -pad, i_end, pad, j, -pad, j_end, -pad, worker->t_prof->N, -pad, worker->q_seq->N );

			/* perform given tasks on them */
			printf_vall("# ==>viterbi...\n");
			WORK_viterbi_and_traceback( worker );
			printf_vall("# ==>forward-backward...\n");
			WORK_forward_backward( worker );
			printf_vall("# ==>cloud search...\n");
			WORK_cloud_search( worker );
			WORK_convert_scores( worker );

			/* output results to file */
			printf("# ==>report result...\n");
			WORK_report_result_current( worker );

			search_cnt++;
		}
	}

	/* report tail output */
   	WORK_report_footer( worker );

   	WORK_close( worker );
	WORK_cleanup( worker );
}

/*
 *  FUNCTION: 	
 *  SYNOPSIS: 	
 */
void main_pipeline_set_flags()
{
	
}