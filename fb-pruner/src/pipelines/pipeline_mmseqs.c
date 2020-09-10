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

/* mmseqs pipeline */
void 
mmseqs_pipeline( WORKER* worker )
{
   /* initialize worker data structures */
   WORK_init( worker );
   WORK_open( worker );

   /* file pointer for writing out to file */
   FILE*    fp       = NULL;
   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;

   /* set tasks */
   tasks->linear        = true;
   tasks->lin_bound_fwd = true;
   // #if DEBUG
   // {
   //    tasks->quadratic  = true;
   //    tasks->quad_vit   = true;
   //    tasks->quad_trace = true;
   //    tasks->quad_fwd   = true;
   //    tasks->quad_bck   = true;
   // }
   // #endif

   /* worker objects */
   TIMES*         times       = worker->times;
   NAT_SCORES*    scores      = worker->scores;
   /* timer */
   CLOCK*         clok        = worker->clok;
   /* alignment for mmseqs window */
   ALIGNMENT*     tr          = worker->traceback;
   /* input results file from MMSEQS pipeline */
   worker->results_in         = RESULTS_Create();
   RESULTS*       results_in  = worker->results_in;
   RESULT*        result_in   = NULL;
   /* output results from fb-pruner */
   RESULT*        result      = worker->result;
   RESULT*        result_out  = worker->result;
   /* 

   /* m8+ file contains target_id, query_id, and result_id fields */
   // RESULTS_M8_Plus_Parse( results_in, args->mmseqs_res_filepath );
   RESULTS_M8_Parse( results_in, args->mmseqs_res_filepath );
   // RESULTS_M8_Dump( results_in, stdout );

   /* target and query ids */
   int      res_id   = -1;
   int      t_id     = -1;
   int      q_id     = -1;
   char*    q_name   = NULL;
   char*    t_name   = NULL;
   int      t_len    = -1;
   int      q_len    = -1;
   /* previous hit ids */
   int      t_id_prv    = -1;
   int      q_id_prv    = -1;
   char*    t_name_prv  = NULL;
   char*    q_name_prv  = NULL;
   /* window start and end points */
   TRACE*   beg;
   TRACE*   end;

   /* TODO: set values for reporting results */
   bool     report_all        = true;
   float    threshold_sc      = 0;

   /* === INDEX FILES === */
   printf_vlo("# Loading index file at: '%s'\n", args->t_indexpath );
   worker->t_index = F_INDEX_Load( worker->t_index, args->t_indexpath );

   printf_vlo("# Loading index file at: '%s'\n", args->q_indexpath );
   worker->q_index = F_INDEX_Load( worker->q_index, args->q_indexpath );
   /* sort indexes by id */
   F_INDEX_Sort_by_Name( worker->t_index );
   F_INDEX_Sort_by_Name( worker->q_index );
   /* output index */
   if ( args->verbose_level >= VERBOSE_ALL ) {
      F_INDEX_Dump( worker->t_index, stdout );
      F_INDEX_Dump( worker->q_index, stdout );
   }

   /* get result range */
   int i_rng = 0;
   int i_cnt = 0;
   int i_beg = 0;
   int i_end = results_in->N;
   if ( args->mmseqs_range.beg >= 0 ) {
      i_beg = args->mmseqs_range.beg;
      i_end = MIN(args->mmseqs_range.end, i_end);
   }
   i_rng = i_end - i_beg;
   printf_vhi("# MMSEQS RESULTS RANGE: (%d,%d)\n", args->mmseqs_range.beg, args->mmseqs_range.end );

   /* add header to all reports */
   WORK_report_header( worker );

   /* === ITERATE OVER EACH RESULT === */
   /* Look through each input result */
   for (int i = i_beg; i < i_end; i++, i_cnt++)
   {
      result_out->result_id = i;
      printf_vlo("# (%d/%d): Running cloud search for result (%d of %d)...\n", 
         i_cnt, i_rng, i, i_end );

      /* get next result from list */
      result_in   = &(results_in->data[i]);
      /* result id */
      res_id      = result_in->result_id;

      printf_vhi("=== M8 Entry : %d ===\n", i);
      if ( args->verbose_level >= VERBOSE_HIGH ) {
         RESULT_M8_Dump( result_in, stdout );
      }

      /* name from results */
      t_name      = result_in->target_name;
      q_name      = result_in->query_name;
      printf_vhi("# T_NAME   (INPUT):\t%s\n", t_name );
      printf_vhi("# Q_NAME   (INPUT):\t%s\n", q_name );

      /* find target/query in index and load from database */
      if ( STRING_Equal( t_name, t_name_prv ) == false ) {
         WORK_load_target_by_name( worker, t_name );
      }
      if ( STRING_Equal( t_name, t_name_prv ) == false )  {
         WORK_load_query_by_name( worker, q_name );
      }

      printf_vhi("# T_NAME (OUTPUT):\t[%d] %s\n", worker->t_id, t_name );
      printf_vhi("# Q_NAME (OUTPUT):\t[%d] %s\n", worker->q_id, q_name );

      t_name_prv  = t_name;
      q_name_prv  = q_name;
      t_len       = worker->t_prof->N;
      q_len       = worker->q_seq->N;

      /* clear old data and change sizes of data structs */
      WORK_reuse( worker );

      // #if DEBUG
      // {
      //    WORK_viterbi_and_traceback( worker );
      //    WORK_forward_backward( worker );
      // }
      // #endif

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

      /* get search window from mmseqs results */
      ALIGNMENT_Reuse( tr, worker->q_seq->N, worker->t_prof->N );
      ALIGNMENT_Pushback( tr, &((TRACE) { result_in->target_start, result_in->query_start, M_ST }) );
      ALIGNMENT_Pushback( tr, &((TRACE) { result_in->target_end, result_in->query_end, M_ST }) );
      tr->beg = 0;
      tr->end = 1;

      /* run cloud search */
      WORK_cloud_search( worker );

      /* convert bitscore to eval */
      WORK_convert_scores( worker );

      #if DEBUG
      {
         float percent_cells = (float) result_out->cloud_cells / (float) result_out->total_cells;
         printf("PRUNING => cloud_cells: %d, total_cells: %d, percent_cells: %.3f\n", 
            result_out->cloud_cells, result_out->total_cells, percent_cells );

         float cloud_tot = times->lin_cloud_fwd + times->lin_cloud_bck +
                     times->lin_merge + times->lin_reorient + times->lin_bound_fwd;
         float speedup = cloud_tot / times->quad_fwd; 
         printf("TIMES => cloud_time: %.3f, fwd_time: %.3f, speed ratio: %.3f\n", 
            cloud_tot, times->quad_fwd, speedup );
         printf("SCORES => cloud_sc: %.3f, fwd_sc: %.3f\n", 
            scores->lin_bound_fwd, scores->quad_fwd );
      }
      #endif

      /* print results */
      WORK_report_result_current( worker );
   }

   WORK_report_footer( worker ); 

   WORK_close( worker );
   WORK_cleanup( worker );
}