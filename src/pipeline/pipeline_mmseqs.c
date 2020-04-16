/*******************************************************************************
 *  FILE:      pipeline_mmseqs.c
 *  PURPOSE:   Cloud Search Pipeline for MMSEQS.
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

/* mmseqs pipeline */
void mmseqs_pipeline( WORKER* worker ) 
{
	/* file pointer for writing out to file */
   FILE* fp  	= NULL;

   /* worker objects */
   ARGS*    args   	= worker->args;
   TASKS*   tasks   	= worker->tasks;
   REPORT*  report   = worker->report;
   TIMES*   times   	= worker->times;
   SCORES* 	scores 	= worker->scores;
   CLOCK*   clock 	= worker->clock;
   /* alignment for mmseqs window */
   ALIGNMENT*	tr 	= worker->traceback;

   /* input results file from MMSEQS pipeline */
   RESULTS* results_in = RESULTS_Create();
   RESULT* 	result;
   RESULTS_M8_Parse( results_in, args->mmseqs_res_filepath );
   RESULTS_M8_Dump( results_in, stdout );

   /* set tasks */
   args->is_testing 		= false;
   tasks->lin_bound_fwd = true;

   /* set values */
   bool report_all 		= true;
   float threshold_sc	= 0;

   /* initialize worker data structures */
   WORK_init( worker );

   /* === INDEX FILES === */
   WORK_load_target_index( worker );
   F_INDEX_Sort_by_Name( worker->t_index );
   WORK_load_query_index( worker );
   F_INDEX_Sort_by_Name( worker->q_index );

   /* add header to file */
   fp = fopen(args->output_filepath, "w+");
   fprintf(fp, ">");
   fprintf(fp, "{%s}\t", "t_id");
   fprintf(fp, "{%s}\t", "t_name");
   fprintf(fp, "{%s}\t", "q_id");
   fprintf(fp, "{%s}\t", "q_name");
   fprintf(fp, "{%s}\t", "cloud_sc");
   fprintf(fp, "{%s}\t", "cloud_time");
   fprintf(fp, "\n");

   /* === ITERATE OVER EACH RESULT === */
   /* Look through each input result */
   for (int i = 0; i < results_in->N; i++) 
   {
      /* get next result from list */
      result = &(results_in->data[i]);

      /* load target and query by looking them up by name */
      WORK_load_target_by_name( worker, result->target_name );
      WORK_load_query_by_name( worker, result->query_name );

      /* change sizes of data structs */
      WORK_reuse( worker );

      /* get search window from mmseqs */
	   ALIGNMENT_Reuse( tr );
	   ALIGNMENT_Pushback( tr, &((TRACE){ result->query_start, result->target_start, M_ST }) );
	   ALIGNMENT_Pushback( tr, &((TRACE){ result->query_end, result->target_end, M_ST }) );
	   tr->beg = 0;
	   tr->end = 1;
      
      /* run cloud search */
      WORK_cloud_search( worker );

      /* if it clears scoring threshold, add to results */
      if ( scores->lin_cloud_fwd > threshold_sc || report_all ) {
         RESULTS_PushBack( worker->results, result );
      }

      /* print out results */
      float cloud_tot = times->lin_cloud_fwd + times->lin_cloud_bck + 
                        times->lin_merge + times->lin_reorient + times->lin_bound_fwd;
      fprintf(fp, "%d\t", result->target_id);
      fprintf(fp, "%s\t", result->target_name);
      fprintf(fp, "%d\t", result->query_id);
      fprintf(fp, "%s\t", result->query_name);
      fprintf(fp, "%.5f\t", scores->lin_cloud_fwd);
      fprintf(fp, "%.5f\t", cloud_tot);
      fprintf(fp, "\n");
   }
   fclose(fp);

   /* final output of results */
   // fp = fopen()
   // fp = stdout;
   // RESULTS_My_Dump( worker->results, fp );
   // if (fp != stdout) fclose(fp);
}