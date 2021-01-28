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
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../parsers/_parsers.h"
#include "../algs_linear/_algs_linear.h"
#include "../algs_quad/_algs_quad.h"
#include "../algs_naive/_algs_naive.h"
#include "../reporting/_reporting.h"
#include "../work/_work.h"

/* header */
#include "_pipelines.h"

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
   /* worker objects */
   TIMES*         times       = worker->times;
   TIMES*         times_tot   = worker->times_totals;
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

   /* set flags for pipeline tasks */
   /* TASKS */
   {
      /* sparse algs */
      tasks->sparse           = true;
      tasks->sparse_bound_fwd = false;
      tasks->sparse_bound_bck = false;
      tasks->sparse_bias_corr = true;
      /* linear algs */
      tasks->linear           = true;     /* if any other linear tasks are flagged, this must be too */
      tasks->lin_fwd          = false;    /* optional, but can't recover alignment */
      tasks->lin_bck          = false;    /* optional, but can't recover alignment */
      tasks->lin_vit          = false;    /* optional, but can't recover alignment */
      tasks->lin_trace        = false;    /* optional, but can't recover alignment */
      tasks->lin_cloud_fwd    = true;     /* required for sparse and linear bound fwdbck */
      tasks->lin_cloud_bck    = true;     /* required for sparse and linear bound fwdbck */
      tasks->lin_bound_fwd    = true;     /* can't be used to recover alignment */  
      tasks->lin_bound_bck    = false;    /* can't be used to recover alignment */
      /* quadratic algs */
      tasks->quadratic        = false;    /* if any other quadratic tasks are flagged, this must be too */
      tasks->quad_fwd         = false;    /* optional */
      tasks->quad_bck         = false;    /* optional */
      tasks->quad_vit         = false;    /* viterbi required for cloud search */
      tasks->quad_trace       = false;    /* traceback required for cloud search  */
      tasks->quad_bound_fwd   = false;    /* required step of cloud search */
      tasks->quad_bound_bck   = false;    /* optional */
      tasks->quad_bias_corr   = false;    /* optional: requires quadratic forward backward */
   }


   /* get result range */
   int i_rng, i_cnt, i_beg, i_end;
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
   /* scores (for use in cutoff thresholds) */
   float    vit_sc;
   float    cloud_sc;
   float    bound_sc;
   /* timers */
   float    iter_start;

   /* load indexes  */
   WORK_load_index_by_name( worker );

   /* get result range */
   WORK_load_mmseqs_list( worker );
   worker->num_searches = args->list_range.end - args->list_range.beg;

   i_beg = args->list_range.beg;
   i_end = args->list_range.end;
   i_rng = i_end - i_beg;
   i_cnt = 0;

   /* add header to all reports */
   WORK_report_header( worker );
   printf_vlo("# Beginning search through mmseqs-m8 list on range (%d,%d)...\n", i_beg, i_end);

   /* === ITERATE OVER EACH RESULT === */
   /* Look through each input result (i = index in full list, i_cnt = index relative to search range) */
   for (int i = i_beg; i < i_end; i++, i_cnt++)
   {
      /* start timer for iteration */
      iter_start = CLOCK_Get_Time( clok );
      /* reset timers */
      WORK_times_init( worker, worker->times );
      /* reset scores */
      WORK_scores_init( worker, worker->scores );

      printf_vlo("\n# (%d/%d): Running cloud search for result (%d of %d)...\n", 
         i_cnt, i_rng, i+1, i_end );

      /* results should only contain entries which we wish to compute */
      result_in = &VEC_X( worker->results_in, i_cnt );
      /* result_id contains the index of the corresponding entry in mmseq's output */
      worker->result->result_id = i;
       /* record the score reported by mmseqs */
      worker->result->vit_natsc = result_in->bit_score;

      /* report the input */
      if ( args->verbose_level >= VERBOSE_LOW ) 
      {
         fprintf( stdout, "=== M8 Entry : [%d] ===\n", i);
         RESULT_M8_Dump( result_in, stdout );
      }

      /* mmseqs viterbi scoring filter */
      vit_sc = result_in->vit_natsc;
      if ( args->filter_on && vit_sc > args->threshold_vit ) {
         printf("MMSEQS viterbi score does not meet threshold. ");
         continue;
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

      /* if the query is not the same as last on list, then get new query */
      if ( STRING_Equal( q_name, q_name_prv ) == false ) {
         WORK_load_query_by_name( worker, q_name );
         // SEQUENCE_Dump( worker->q_seq, stdout );
      }

      /* lookback 1 name */
      t_name_prv  = t_name;
      q_name_prv  = q_name;

      /* clear old data and change sizes of data structs */
      WORK_reuse( worker );

      // WORK_load_m8();

      /* get search window from mmseqs results */
      ALIGNMENT_Reuse( tr, worker->q_seq->N, worker->t_prof->N );
      ALIGNMENT_Pushback( tr, &((TRACE) { result_in->target_start, result_in->query_start, M_ST }) );
      ALIGNMENT_Pushback( tr, &((TRACE) { result_in->target_end, result_in->query_end, M_ST }) );
      tr->beg = 0;
      tr->end = 1;
      printf_vhi("DIM:: TARGET: {%d}, QUERY: {%d}\n", worker->q_seq->N, worker->t_prof->N );
      printf_vhi("VIT_TRACEBACK:: {%6d,%6d}->{%6d,%6d}\n", result_in->target_start, result_in->query_start, result_in->target_end, result_in->query_end );
      
      /* run cloud search */
      WORK_cloud_search( worker );

      /* cloud search scoring filter */
      cloud_sc = worker->scores->threshold_cloud_compo;
      if ( args->filter_on && cloud_sc > args->threshold_cloud ) {
         printf("MMSEQS viterbi score does not meet threshold. ");
         continue;
      }

      /* run bounded forward/backward */
      WORK_bound_forward_backward( worker );

      /* bound forward scoring filter */
      bound_sc = worker->result->bound_fwd_natsc;
      if ( args->filter_on && bound_sc > args->threshold_bound_fwd ) {
         printf("MMSEQS viterbi score does not meet threshold. ");
         continue;
      }

      /* compute posterior and bias, find domains, compute domain-specific posterior and bias */
      WORK_posterior( worker );

      /* capture total runtime for current iteration */
      worker->times->total = CLOCK_Get_Diff( clok, iter_start, CLOCK_Get_Time( clok ) );
      printf("RESULT_TIME: %f\n\n", worker->times->total);

      /* add current iteration times to totals */
      WORK_times_add( worker );

      /* print results */
      WORK_report_result_current( worker );
   }

   /* capture total-ish program runtime (from clock init to end of main loop) */
   worker->times_totals->program_runtime = CLOCK_Get_Diff( clok, clok->program_start, CLOCK_Get_Time( clok ) );

   /* report total times */
   REPORT_mytimeout_totals( worker, worker->result, stdout );

   /* results footer */
   WORK_report_footer( worker ); 
   /* close all file pointers */
   WORK_close( worker );

   /* free work data */
   WORK_cleanup( worker );
}

STATUS_FLAG
mmseqs_pipeline_Set_Default_Tasks()
{

}

STATUS_FLAG
mmseqs_pipeline_Set_Default_Args()
{

}