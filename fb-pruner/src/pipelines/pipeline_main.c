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
		tasks->lin_fwd 			= true;
		tasks->lin_bck 			= true;
		tasks->lin_vit 			= false;	/* currently not supported */
		tasks->lin_trace 		= false;	/* currently not supported */
		tasks->lin_bound_fwd 	= true;
		tasks->lin_bound_bck 	= true;
		/* quadratic algs */
		tasks->quadratic 		= true;		/* if any other quadratic tasks are flagged, this must be too */
		tasks->quad_fwd 		= true;		/* optional */
		tasks->quad_bck 		= true;		/* optional */
		tasks->quad_vit 		= true;		/* viterbi required for cloud search */
		tasks->quad_trace 		= true;		/* traceback required for cloud search  */
		tasks->quad_bound_fwd 	= true;		/* required step of cloud search */
		tasks->quad_bound_bck 	= false;	/* */
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
	 	report->quad_fwd_sc 		= true;
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

	/* load index files */
	WORK_load_target_index( worker );
	WORK_load_query_index( worker );
	F_INDEX_Sort_by_Id( worker->t_index );
   	F_INDEX_Sort_by_Id( worker->q_index );

	/* allocate data structs */
	worker->t_prof 	= HMM_PROFILE_Create();
	worker->t_seq	= SEQUENCE_Create();
	worker->q_seq	= SEQUENCE_Create();

	/* set and verify ranges */
	WORK_set_ranges( worker );

	printf("target_num: %d, query_num: %d\n", 
		worker->t_index->N, worker->q_index->N );
	printf("target_range: (%d,%d), query_range: (%d,%d)\n", 
		args->t_range.beg, args->t_range.end,
		args->q_range.beg, args->q_range.end );

	int search_cnt 	= 0;
	int search_tot 	= (args->t_range.end - args->t_range.beg) * (args->q_range.end - args->q_range.beg);

	int i_beg 	= args->t_range.beg;
	int i_end 	= args->t_range.end;
	int j_beg 	= args->q_range.beg;
	int j_end 	= args->q_range.end;

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
				"##_SCORES_TIMES_: %d %d %s %d %d %s %d %d %f %f %f %f %f %f ",
				worker->t_id, worker->t_prof->N, worker->t_prof->name, 
				worker->q_id, worker->q_seq->N, worker->q_seq->name,
				result->total_cells, result->cloud_cells, 
				times->quad_vit, scores->quad_vit, 
				times->quad_fwd, scores->quad_fwd, 
				times->lin_total_cloud, scores->lin_cloud_fwd );

			search_cnt++;
		}
	}
	printf("...done.\n");

	fclose( worker->out_file );
}