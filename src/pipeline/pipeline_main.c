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

/* data structures */
#include "objects/structs.h"
#include "utilities/utility.h"
#include "error_handler.h"

/* file parsers */
#include "parsers/arg_parser.h"
#include "parsers/hmm_parser.h"
#include "parsers/seq_parser.h"
#include "parsers/m8_parser.h"
#include "parsers/index_parser.h"
#include "parsers/seq_to_profile.h"

/* objects */
#include "objects/f_index.h"
#include "objects/results.h"
#include "objects/alignment.h"
#include "objects/sequence.h"
#include "objects/hmm_profile.h"
#include "objects/edgebound.h"
#include "objects/clock.h"
#include "objects/matrix/matrix_2d.h"
#include "objects/matrix/matrix_3d.h"
#include "objects/mystring.h"
#include "objects/vectors/vector_range_2d.h"

/* viterbi & fwdbck (quadratic) */
#include "algs_quad/viterbi_quad.h"
#include "algs_quad/traceback_quad.h"
#include "algs_quad/fwdback_quad.h"
/* viterbi & fwdbck (linear) */
#include "algs_linear/fwdback_linear.h"

/* cloud search (naive) */
#include "algs_naive/bounded_fwdbck_naive.h"
/* cloud search (quadratic space) */
#include "algs_quad/cloud_search_quad.h"
#include "algs_quad/merge_reorient_quad.h"
#include "algs_quad/bounded_fwdbck_quad.h"
/* cloud search (linear space) */
#include "algs_linear/cloud_search_linear.h"
#include "algs_linear/merge_reorient_linear.h"
#include "algs_linear/bounded_fwdbck_linear.h"
/* temp test */
#include "algs_linear/cloud_search_linear_rows.h"

/* debugging methods */
#include "testing.h"

/* header */
#include "pipeline.h"

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

	/* set file pointer */
	worker->out_file = fopen( args->output_filepath, "w+" );
	WORK_print_result_header( worker );

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
		tasks->quad_fwd 		= false;
		tasks->quad_bck 		= false;
		tasks->quad_vit 		= true;		/* cloud search depends on viterbi */
		tasks->quad_trace 		= true;		/* cloud search depends on traceback */
		tasks->quad_bound_fwd 	= false;
		tasks->quad_bound_bck 	= false;
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

	/* initialize data structures needed for tasks */
	WORK_init( worker );

	/* load target index file */
	WORK_load_target_index( worker );
	// F_INDEX_Dump( worker->t_index, stdout );
	/* load query index file */
	WORK_load_query_index( worker );
	// F_INDEX_Dump( worker->q_index, stdout );

	/* allocate data structs */
	worker->t_prof 	= HMM_PROFILE_Create();
	worker->t_seq	= SEQUENCE_Create();
	worker->q_seq	= SEQUENCE_Create();

	/* set and verify ranges */
	WORK_set_ranges( worker );

	/* loop over targets */
	for (int i = args->t_range.beg; i < args->t_range.end; i++) 
	{
		/* load in next target */
		WORK_load_target_by_id( worker, i );
		// HMM_PROFILE_Dump( worker->t_prof, stdout );
		result->target_id = i;

		/* loop over queries */
		for (int j = args->q_range.beg; j < args->q_range.end; j++) 
		{
			/* load in next query */
			WORK_load_query_by_id( worker, j );
			// SEQUENCE_Dump( worker->q_seq, stdout );
			result->query_id = j;

			/* resize dp matrices to new query/target */
			WORK_reuse( worker );

			/* perform given tasks on them */
			// printf("viterbi...\n");
			WORK_viterbi_and_traceback( worker );
			// printf("forward-backward...\n");
			WORK_forward_backward( worker );
			// printf("cloud search...\n");
			WORK_cloud_search( worker );

			/* output results to file */
			WORK_print_result_current( worker );
		}
	}

	fclose( worker->out_file );
}