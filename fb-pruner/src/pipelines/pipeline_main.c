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
	/* worker objects */
	ARGS* 		args 		= worker->args;
	TASKS* 		tasks 		= worker->tasks;
	TIMES* 		times 		= worker->times;
	SCORES* 	scores 		= worker->scores;
	REPORT* 	report 		= worker->report;
	RESULTS* 	results 	= worker->results;
	RESULT* 	result 		= worker->result;

	/* set flags for pipeline tasks */
	/* linear algs */
	{
		tasks->linear 			= true;		/* if any other linear tasks are flagged, this must be too */
		tasks->lin_fwd 			= false;	/* currently not supported */
		tasks->lin_bck 			= false;	/* currently not supported */
		tasks->lin_vit 			= false;	/* currently not supported */
		tasks->lin_trace 		= false;	/* currently not supported */
		tasks->lin_bound_fwd 	= true;
		tasks->lin_bound_bck 	= true;
		/* quadratic algs */
		tasks->quadratic 		= true;		/* if any other quadratic tasks are flagged, this must be too */
		tasks->quad_fwd 		= false;	/* optional */
		tasks->quad_bck 		= false;	/* optional */
		tasks->quad_vit 		= true;		/* viterbi required for cloud search */
		tasks->quad_trace 		= true;		/* traceback required for cloud search  */
		tasks->quad_bound_fwd 	= false;	/* required step of cloud search */
		tasks->quad_bound_bck 	= false;	/* optional */
	}

	/* set flags for reporting */
	/* SCORES */
	{
		/* linear algs */
		report->lin_fwd_sc 			= true;
		report->lin_bck_sc 			= true;
		report->lin_vit_sc 			= true;
		report->lin_bound_fwd_sc 	= true;
		report->lin_bound_bck_sc 	= true;
		/* quadratic algs */
	 	report->quad_fwd_sc 		= false;
		report->quad_bck_sc 		= false;
		report->quad_vit_sc 		= false;
		report->quad_bound_fwd_sc 	= false;
		report->quad_bound_bck_sc 	= false;
	}

	/* TIMES */
	{
		/* linear algs */
		report->lin_fwd_t 			= true;
		report->lin_bck_t 			= true;
		report->lin_vit_t 			= true;
		report->lin_trace_t 		= true;
		report->lin_cloud_fwd_t 	= true;
		report->lin_cloud_bck_t 	= true;
		report->lin_merge_t 		= true;
		report->lin_reorient_t 		= true;
		report->lin_bound_fwd_t 	= true;
		report->lin_bound_bck_t 	= true;
		/* quadratic algs */
		report->quad_fwd_t 			= false;
		report->quad_bck_t 			= false;
		report->quad_vit_t 			= false;
		report->quad_trace_t 		= false;
		report->quad_cloud_fwd_t 	= false;
		report->quad_cloud_bck_t 	= false;
		report->quad_merge_t 		= false;
		report->quad_reorient_t 	= false;
		report->quad_bound_fwd_t 	= false;
		report->quad_bound_bck_t 	= false;
	}

	/* set file pointer */
	worker->out_file = fopen( args->output_filepath, "w+" );
	if ( worker->out_file == NULL ) {
		fprintf(stderr, "ERROR: Bad File Pointer.\n");
		exit(EXIT_FAILURE);
	}
	WORK_print_result_header( worker );

	/* initialize data structures needed for tasks */
	WORK_init( worker );

	/* load or build, then sort target and query index files */
	WORK_index( worker );

	/* set and verify index ranges */
	WORK_set_ranges( worker );
	/* report ranges */
	printf_vlo("target_num: %d, query_num: %d\n", 
		worker->t_index->N, worker->q_index->N );
	printf_vlo("target_range: (%d,%d), query_range: (%d,%d)\n", 
		args->t_range.beg, args->t_range.end,
		args->q_range.beg, args->q_range.end );

	int i_beg, i_end, j_beg, j_end;
	int search_cnt, search_tot;
	i_beg 	= args->t_range.beg;
	i_end 	= args->t_range.end;
	j_beg 	= args->q_range.beg;
	j_end 	= args->q_range.end;
	search_cnt 	= 0;
	search_tot 	= (i_end - i_beg) * (j_end - j_beg);

	/* loop over target range */
	for (int i = i_beg; i < i_end; i++) 
	{
		/* load in next target */
		WORK_load_target_by_id( worker, i );
		// HMM_PROFILE_Dump( worker->t_prof, stdout );
		result->target_id = i;

		/* loop over query range */
		for (int j = j_beg; j < j_end; j++) 
		{
			printf_vlo("Running search (%d of %d) for target (%d of %d) and query (%d of %d)...\n",
				search_cnt, search_tot, i, i_end, j, j_end);

			/* load in next query */
			WORK_load_query_by_id( worker, j );
			// SEQUENCE_Dump( worker->q_seq, stdout );
			result->query_id = j;

			/* resize dp matrices to new query/target */
			WORK_reuse( worker );

			const int pad = 5;
			printf_vall(" (%*d) | T_ID: (%*d/%*d) | Q_ID: (%*d/%*d) | T_LEN: %*d | Q_LEN: %*d \n", 
				pad, search_cnt, pad, i, -pad, i_end, pad, j, -pad, j_end, -pad, worker->t_prof->N, -pad, worker->q_seq->N );

			/* perform given tasks on them */
			printf_vall("viterbi...\n");
			WORK_viterbi_and_traceback( worker );
			printf_vall("forward-backward...\n");
			WORK_forward_backward( worker );
			printf_vall("cloud search...\n");
			WORK_cloud_search( worker );

			/* output results to file */
			WORK_print_result_current( worker );
			STRING_Replace( worker->t_prof->name, ' ', '_' );
			STRING_Replace( worker->q_seq->name, ' ', '_' );
			fprintf( stdout, 
				"##_SCORES_TIMES_: %d %d %s %d %d %s %d %d %f %f %d %f %f %f %f \n",
				worker->t_id, worker->t_prof->N, worker->t_prof->name, 
				worker->q_id, worker->q_seq->N, worker->q_seq->name,
				result->total_cells, result->cloud_cells, 
				args->alpha, args->alpha_max, args->beta,
				times->quad_vit, scores->quad_vit,
				times->lin_total_cloud, scores->lin_cloud_fwd );

			search_cnt++;
		}
	}
	printf("...done.\n");

	fclose( worker->out_file );
}

/*
 *  FUNCTION: 	
 *  SYNOPSIS: 	
 */
void main_pipeline_set_flags()
{
	
}