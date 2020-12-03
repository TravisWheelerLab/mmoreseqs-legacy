/*******************************************************************************
 *  FILE:      pipeline_mmseqs.c
 *  PURPOSE:   MMORE-SEQS full pipeline.  Runs MMseqs, then passes results to MMORE-SEQS.
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
#include "../objects/structs.h"
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../parsers/parsers.h"
#include "../algs_linear/algs_linear.h"
#include "../algs_quad/algs_quad.h"
#include "../algs_naive/algs_naive.h"

/* header */
#include "pipelines.h"

/* mmseqs pipeline */
void 
mmseqs_pipeline( WORKER* worker )
{
   printf("=== MMSEQS PIPELINE ===\n");

   /* initialize worker data structures */
   WORK_init( worker );
   WORK_open( worker );

   /* file pointer for writing out to file */
   ARGS*          args        = worker->args;
   TASKS*         tasks       = worker->tasks;
   /* set tasks */
   tasks->linear              = true;
   tasks->lin_bound_fwd       = true;
   tasks->lin_bound_bck       = true;
   /* worker objects */
   TIMES*         times       = worker->times;
   NAT_SCORES*    scores      = worker->scores;
   /* timer */
   CLOCK*         clok        = worker->clok;
   /* alignment for mmseqs window */
   ALIGNMENT*     tr          = worker->traceback;
   /* input results */
   RESULTS*       results_in  = worker->results_in;
   RESULT*        result_in   = NULL;
   /* input results */
   RESULTS*       results_out = NULL;
   RESULT*        result_out  = NULL;

   /* counter for skipped alignments */
   int   skips = 0;

   /* set flags for pipeline tasks */
   /* TASKS */
   {
      /* sparse algs */
      tasks->sparse           = true;
      tasks->sparse_bound_fwd = true;
      tasks->sparse_bound_bck = true;
      tasks->sparse_bias_corr = false;
      /* linear algs */
      tasks->linear           = true;     /* if any other linear tasks are flagged, this must be too */
      tasks->lin_fwd          = true;     /* optional, but can't recover alignment */
      tasks->lin_bck          = true;     /* optional, but can't recover alignment */
      tasks->lin_vit          = false;    /* optional, but can't recover alignment */
      tasks->lin_trace        = false;    /* optional, but can't recover alignment */
      tasks->lin_bound_fwd    = true;     /* */  
      tasks->lin_bound_bck    = true;     /* */
      /* quadratic algs */
      tasks->quadratic        = false;    /* if any other quadratic tasks are flagged, this must be too */
      tasks->quad_fwd         = false;    /* optional */
      tasks->quad_bck         = false;    /* optional */
      tasks->quad_vit         = false;    /* viterbi required for cloud search */
      tasks->quad_trace       = false;    /* traceback required for cloud search  */
      tasks->quad_bound_fwd   = false;    /* required step of cloud search */
      tasks->quad_bound_bck   = false;    /* optional */
      tasks->quad_bias_corr   = false;
   }

   if ( args->is_compo_bias == BIAS_CORR_QUAD ) {
      tasks->quad_bias_corr = true;
   } 
   else if ( args->is_compo_bias == BIAS_CORR_SPARSE ) {
      tasks->sparse_bias_corr = true;
   }

   /* get result range */
   int i_rng, i_cnt, i_beg, i_end;
   i_cnt = 0;

   /* m8+ file contains target_id, query_id, and result_id fields */
   RESULTS_M8_Parse( results_in, args->mmseqs_res_filepath, args->list_range.beg, args->list_range.end );
   // RESULTS_M8_Dump( results_in, stdout );

   /* current hit */
   int      res_id            = -1;
   int      t_id              = -1;
   int      q_id              = -1;
   char*    q_name            = NULL;
   char*    t_name            = NULL;
   /* previous hit */
   int      t_id_prv          = -1;
   int      q_id_prv          = -1;
   char*    t_name_prv        = NULL;
   char*    q_name_prv        = NULL;
   /* window start and end points */
   TRACE*   beg               = NULL;
   TRACE*   end               = NULL;
   /* TODO: set values for reporting results */
   bool     report_all        = true;
   float    threshold_sc      = 0;

   /* load indexes  */
   WORK_load_index_by_name( worker );
   /* get result range */
   WORK_load_mmseqs_list( worker );
   worker->num_searches = args->list_range.end - args->list_range.beg;

   i_beg = args->list_range.beg;
   i_end = args->list_range.end;
   i_rng = i_end - i_beg;

   /* add header to all reports */
   WORK_report_header( worker );

   printf("# Beginning search through mmseqs-m8 list on range (%d,%d)...\n", i_beg, i_end);

   /* === ITERATE OVER EACH RESULT === */
   /* Look through each input result */
   for (int i = i_beg; i < i_end; i++, i_cnt++)
   {
      printf_vlo("\n# (%d/%d): Running cloud search for result (%d of %d)...\n", 
         i_cnt, i_rng, i+1, i_end );

      result_in = &worker->results_in->data[i];

      if ( args->verbose_level >= VERBOSE_LOW ) 
      {
         fprintf( stdout, "=== M8 Entry : [%d] ===\n", i);
         RESULT_M8_Dump( result_in, stdout );
      }

      /* name from results */
      int pad     = -40;
      t_name      = result_in->target_name;
      q_name      = result_in->query_name;
      printf_vhi("# T_NAME:  [cur] %*s =>  [prv] %*s\n", pad, t_name, pad, t_name_prv );
      printf_vhi("# Q_NAME:  [cur] %*s =>  [prv] %*s\n", pad, q_name, pad, q_name_prv );

      /* find target/query in index and load from database */
      if ( STRING_Equal( t_name, t_name_prv ) == false ) {
         WORK_load_target_by_name( worker, t_name );
         // HMM_PROFILE_Dump( worker->t_prof, stdout );
      }

      /* start time for current */
      worker->result->time = CLOCK_Get_RealTime();

      if ( STRING_Equal( q_name, q_name_prv ) == false ) {
         WORK_load_query_by_name( worker, q_name );
         // SEQUENCE_Dump( worker->q_seq, stdout );
      }

      t_name_prv  = t_name;
      q_name_prv  = q_name;

      /* clear old data and change sizes of data structs */
      WORK_reuse( worker );

      /* get search window from mmseqs results */
      ALIGNMENT_Reuse( tr, worker->q_seq->N, worker->t_prof->N );
      ALIGNMENT_Pushback( tr, &((TRACE) { result_in->target_start, result_in->query_start, M_ST }) );
      ALIGNMENT_Pushback( tr, &((TRACE) { result_in->target_end, result_in->query_end, M_ST }) );
      tr->beg = 0;
      tr->end = 1;
      printf_vhi("DIM:: TARGET: {%d}, QUERY: {%d}\n", worker->q_seq->N, worker->t_prof->N );
      printf_vhi("VIT_TRACEBACK:: {%6d,%6d}->{%6d,%6d}\n", result_in->target_start, result_in->query_start, result_in->target_end, result_in->query_end );

      // #if DEBUG
      // {
      //    /* get search window by running viterbi (debug only) */
      //    tr = worker->traceback;
      //    beg = &(tr->traces[tr->beg]);
      //    end = &(tr->traces[tr->end]);
      //    printf_vall("MY TRACEBACK: (%d,%d) -> (%d,%d)\n", beg->i, beg->j, end->i, end->j);
      // }
      // #else 
      // {
      //    /* get search window from mmseqs results */
      //    ALIGNMENT_Reuse( tr, worker->q_seq->N, worker->t_prof->N );
      //    ALIGNMENT_Pushback( tr, &((TRACE) { result_in->target_start, result_in->query_start, M_ST }) );
      //    ALIGNMENT_Pushback( tr, &((TRACE) { result_in->target_end, result_in->query_end, M_ST }) );
      //    tr->beg = 0;
      //    tr->end = 1;
      // }
      // #endif

      // /* if fixing mmseqs alignment, we need to trim in the alignment if it exceeds the bounds of the profile or sequence */
      // if ( args->adjust_mmseqs_alns == true )
      // {
      //    int fixed_aln  = false;
      //    int t_fix = false;
      //    int q_fix = false;

      //    if ( result_in->target_start > worker->q_seq->N - 1 ) {
      //       result_in->target_start = worker->q_seq->N - 1;
      //       fixed_aln = true;
      //       t_fix = true;
      //    } 
      //    if ( result_in->target_end > worker->q_seq->N - 1 ) {
      //       result_in->target_end = worker->q_seq->N - 1;
      //       fixed_aln = true;
      //       t_fix = true;
      //    }
      //    if ( result_in->query_start > worker->t_prof->N - 1 ) {
      //       result_in->query_start = worker->t_prof->N - 1;
      //       fixed_aln = true;
      //       q_fix = true;
      //    } 
      //    if ( result_in->query_end > worker->t_prof->N - 1 ) {
      //       result_in->query_end = worker->t_prof->N - 1;
      //       fixed_aln = true;
      //       q_fix = true;
      //    }
      //    /* log out viterbi alignment adjustments to error file */ 
      //    if (fixed_aln == true) {
      //       FILE* ferr = fopen("fixed_aln.err", "a");
      //       char* t_str = "<target>";
      //       char* q_str = "<query>";
      //       int pad = -6;
      //       /* log to error file */
      //       fprintf(ferr, "%*sERROR occurred on mmseqs results (%d). %s %s was out-of-bounds. Target_size = %d, Query_size = %d\n%*s", 
      //          pad, ">>", i, (t_fix ? t_str : ""), (q_fix ? q_str : ""), worker->t_prof->N, worker->q_seq->N, pad, "" );
      //       RESULT_M8_Dump( result_in, ferr );
      //       fprintf(ferr, "%*sMMseqs coords: (%d,%d)=>(%d,%d) || Adjusted coords: (%d,%d)=>(%d,%d)\n",
      //          pad, "", q_mmseqs.i, t_mmseqs.i, q_mmseqs.j, t_mmseqs.j, result_in->query_start, result_in->target_start, result_in->query_end, result_in->target_end );
      //       fclose(ferr);
      //       /* log to standard error */
      //       fprintf(stderr, "%*sERROR occurred on mmseqs results (%d). %s %s was out-of-bounds.  Query_size = %d, Target_size = %d\n", 
      //          pad, "#>>", i, (q_fix ? q_str : ""), (t_fix ? t_str : ""), worker->q_seq->N, worker->t_prof->N );
      //       RESULT_M8_Dump( result_in, ferr );
      //       fprintf(stderr, "%*sMMseqs coords: (%d,%d)=>(%d,%d) || Adjusted coords: (%d,%d)=>(%d,%d)\n",
      //          pad, "#",  q_mmseqs.i, t_mmseqs.i, q_mmseqs.j, t_mmseqs.j, result_in->query_start, result_in->target_start, result_in->query_end, result_in->target_end );
      //    }
      //    // if (fixed_aln == true) {
      //    //    skips++;
      //    //    printf("# out-of-bounds. skipping this alignment (%d)...\n", skips);
      //    //    continue;
      //    // }
      // }

      #if DEBUG
      {
         // WORK_forward_backward( worker );
      }
      #endif

      /* run cloud search */
      WORK_cloud_search( worker );
      /* capture alignment */
      // WORK_capture_alignment( worker );
      /* convert bitscore to eval */
      WORK_convert_scores( worker );

      /* end time */
      worker->result->time = CLOCK_Get_RealTime() - worker->clok->start;
      printf("MYTIME: %f\n", worker->result->time);

      if ( args->verbose_level >= VERBOSE_LOW || true  ) 
      {
         RESULT* result_out = worker->result;

         float percent_cells = (float) result_out->cloud_cells / (float) result_out->total_cells;
         float cloud_tot = times->lin_cloud_fwd + times->lin_cloud_bck + times->lin_merge + times->lin_reorient + times->lin_bound_fwd;
         float speedup = cloud_tot / times->quad_fwd; 

         printf("PRUNING =>  cloud_cells: %5d, total_cells: %5d, percent_cells: %2.3f\n", 
            result_out->cloud_cells, result_out->total_cells, percent_cells );
         printf("TIMES   =>  cloud_time: %2.3f, fwd_time: %2.3f, speed ratio: %2.3f\n", 
            cloud_tot, times->quad_fwd, speedup );
         printf("SCORES  =>  lin_bound_fwd_sc: %2.3f, sparse_bound_fwd_sc: %2.3f, fwd_sc: %2.3f, \n", 
            scores->lin_bound_fwd, scores->sparse_bound_fwd, scores->quad_fwd );
         printf("SCORES  =>  lin_bound_bck_sc: %2.3f, sparse_bound_bck_sc: %2.3f, bck_sc: %2.3f, \n", 
            scores->lin_bound_bck, scores->sparse_bound_bck, scores->quad_bck );
      }

      /* print results */
      WORK_report_result_current( worker );
   }

   /* results footer */
   WORK_report_footer( worker ); 
   /* close all file pointers */
   WORK_close( worker );
   /* free work data */
   WORK_cleanup( worker );
}