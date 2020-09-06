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

/* private functions */
void 
print_result_mmseqs(    WORKER*  worker,
                        ARGS*    args,
                        RESULT*  result,
                        SCORES*  scores,
                        TIMES*   times );

void 
print_header_mmseqs(  WORKER*  worker );

/* mmseqs pipeline */
void 
mmseqs_pipeline( WORKER* worker )
{
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

   /* initialize worker data structures */
   WORK_init( worker );

   /* worker objects */
   TIMES*         times       = worker->times;
   NAT_SCORES*    scores      = worker->scores;

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

   /* m8+ file contains target_id, query_id, and result_id fields */
   RESULTS_M8_Plus_Parse( results_in, args->mmseqs_res_filepath );
   // RESULTS_M8_Dump( results_in, stdout );

   /* target and query ids */
   int      res_id   = -1;
   int      t_cid    = -1;
   int      q_cid    = -1;
   int      t_mid    = -1;
   int      q_mid    = -1;
   char*    q_name   = NULL;
   char*    t_name   = NULL;
   int      t_len    = -1;
   int      q_len    = -1;
   /* previous hit ids */
   int      t_cid_prv    = -1;
   int      q_cid_prv    = -1;
   /* window start and end points */
   TRACE*   beg;
   TRACE*   end;

   /* TODO: set values for reporting results */
   bool     report_all        = true;
   float    threshold_sc      = 0;

   /* === INDEX FILES === */
   printf("# Loading index file at: '%s'\n", args->t_indexpath );
   worker->t_index = F_INDEX_Load( worker->t_index, args->t_indexpath );

   printf("# Loading index file at: '%s'\n", args->q_indexpath );
   worker->q_index = F_INDEX_Load( worker->q_index, args->q_indexpath );
   /* sort indexes by id */
   F_INDEX_Sort_by_Id( worker->t_index );
   F_INDEX_Sort_by_Id( worker->q_index );
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
   printf_vhi("# MMSEQS PLUS RANGE: (%d,%d)\n", args->mmseqs_range.beg, args->mmseqs_range.end );

   /* open outfile and add header to file */
   worker->output_fp = fopen(args->output_filepath, "w+");
   fp = worker->output_fp;
   
   WORK_report_result_current( worker );

   /* === ITERATE OVER EACH RESULT === */
   /* Look through each input result */
   for (int i = i_beg; i < i_end; i++, i_cnt++)
   {
      result_out->result_id = i;
      printf_vlo("# (%d/%d): Running cloud search for result (%d of %d)...\n", 
         i_cnt, i_rng, i, i_end );

      /* NOTE: Do not need to keep record of mmseqs id */

      /* get next result from list */
      result_in   = &(results_in->data[i]);
      /* result id */
      res_id      = result_in->result_id;
      /* mmseqs ids from result */
      t_mid       = result_in->target_mid;
      q_mid       = result_in->query_mid;
      /* cloud ids from result */
      t_cid       = result_in->target_id;
      q_cid       = result_in->query_id;
      /* name from results */
      q_name      = result_in->query_name;
      t_name      = result_in->target_name;

      /* list mmseqs and cloud ids */
      printf_vhi("# t_mid: %d, q_mid: %d, t_cid: %d, q_cid: %d\n",
         t_mid, q_mid, t_cid, q_cid );

      /* TODO: should swap query and target program-wide? */
      /* NOTE: query and target are cross-labeled compared to mmseqs */
      int swp  = q_cid;
      q_cid    = t_cid;
      t_cid    = swp;

      /* load target and query by looking them up by id (if we aren't using the same from last search) */
      printf_vhi("# t_cid = %d -> %d, q_cid = %d -> %d\n", 
         t_cid, t_cid_prv, q_cid, q_cid_prv);

      if ( t_cid != t_cid_prv ) {
         WORK_load_target_by_id( worker, t_cid );
      }
      if (q_cid != q_cid_prv) {
         WORK_load_query_by_id( worker, q_cid );
      }
      t_cid_prv = t_cid;
      q_cid_prv = q_cid;

      printf_vhi("# T_NAME   (LOAD):\t%s \n", worker->t_prof->name );
      printf_vhi("# Q_NAME   (LOAD):\t%s \n", worker->q_seq->name );
      printf_vhi("# T_NAME (RESULT):\t%s \n", t_name );
      printf_vhi("# Q_NAME (RESULT):\t%s \n", q_name );

      t_len    = worker->t_prof->N;
      q_len    = worker->q_seq->N;

      /* change sizes of data structs */
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

      /* print search terms */
      const int pad = 10;
      printf_vhi(" (%d/%d) | T_ID: %*d | Q_ID: %*d | T_LEN: %*d | Q_LEN: %*d | BEG: (%d,%d) | END: (%d,%d) \n",
            i, i_end, -pad, t_cid, -pad, q_cid, -pad, t_len, -pad, q_len,
            result_in->query_start, result_in->target_start, result_in->query_end, result_in->target_end );

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

      /* NOTE: report bias */
      // fprintf( stdout, 
      //       "##_SCORES_TIMES_: %d %d | %d %d %25s | %d %d %25s | %d %d | %4.2f %4.2f %d | %d %9.2e | %7.4f %7.4f %9.2e\n",
      //       i, i_end,
      //       worker->t_id, worker->t_prof->N, worker->t_prof->name, 
      //       worker->q_id, worker->q_seq->N, worker->q_seq->name,
      //       result->total_cells, result->cloud_cells, 
      //       args->alpha, args->beta, args->gamma,
      //       result_in->bit_score, result_in->e_value,
      //       times->lin_bound_total, scores->lin_cloud_fwd, evals->lin_cloud_fwd );

      /* print results */
      WORK_report_result_current( worker );
   }

   WORK_report_footer( worker ); 
   WORK_cleanup( worker );
}