/*******************************************************************************
 *  FILE:      pipeline_time.c
 *  PURPOSE:   Time Trial Cloud Search Pipeline.
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

/* data structures and file parsers */
#include "objects/structs.h"
#include "utility.h"
#include "seq_parser.h"
#include "hmm_parser.h"

/* objects */
#include "objects/sequence.h"
#include "objects/hmm_profile.h"
#include "objects/alignment.h"
#include "objects/edgebound.h"
#include "objects/clock.h"
#include "objects/matrix/matrix_2d.h"
#include "objects/matrix/matrix_3d.h"
#include "objects/vectors/vector_range_2d.h"

/* viterbi & fwdbck (quadratic) */
#include "viterbi_quad.h"
#include "traceback_quad.h"
#include "fwdback_quad.h"

/* cloud search (naive) */
#include "bounded_fwdbck_naive.h"
/* cloud search (quadratic space) */
#include "cloud_search_quad.h"
#include "merge_reorient_quad.h"
#include "bounded_fwdbck_quad.h"
/* cloud search (linear space) */
#include "cloud_search_linear.h"
#include "merge_reorient_linear.h"
#include "bounded_fwdbck_linear.h"
/* temp test */
#include "cloud_search_linear_rows.h"

/* debugging methods */
#include "testing.h"

/* header */
#include "pipeline.h"

/* ****************************************************************************************** *
 *  
 *  FUNCTION:  time_pipeline()
 *  SYNOPSIS:  Runs a workflow pipeline. 
 *             Takes in a single target/query pair.  
 *             Runs generic Forward-Backward algorithm
 *
 *  ARGS:      <args>     parsed commandline arguments
 *
 *  RETURN:    No Return.
 *
/* ****************************************************************************************** */
void time_pipeline(ARGS* args) 
{
   /* Set Commandline Arguments */
   float          alpha             = args->alpha;
   int            beta              = args->beta;

   char*          hmm_file          = args->target_filepath;
   char*          fasta_file        = args->query_filepath;

   int            t_filetype        = args->target_filetype;
   int            q_filetype        = args->query_filetype;

   char*          out_file          = args->output_filepath;

   /* Set Mode */
   // int            mode           = MODE_MULTILOCAL;    /* HMMER standard mode (allows jumps) */
   int            mode              = MODE_UNILOCAL;      /* Cloud mode (prohibiits jumps) */

   /* Target & Query */
   HMM_PROFILE*   target_prof       = NULL;
   SEQUENCE*      query_seq         = NULL;

   /* Times & Scores */
   CLOCK*         cl                = CLOCK_Create();
   SCORES*        scores            = NULL;
   TIMES*         times             = NULL;

   bool           is_testing        = false;
   bool           test              = false;
   time_t         t                 = 0; 
   time_t         tot_t             = 0;

   int            T                 = 0;
   int            Q                 = 0;

   float          sc                = 0.f; 
   float          perc_cells        = 0.f;
   int            num_cells         = 0;  
   int            window_cells      = 0; 
   int            tot_cells         = 0;

   /* FILE OUTPUT */
   FILE*          fp                = NULL;
   int            pad               = 0;   

   /* ----------------------------------------------------------------------------------------- */

   /* Allocate memory for time & scores */
   scores = (SCORES*) malloc( sizeof(SCORES) );
   if (scores == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc for SCORES.\n");
      exit(EXIT_FAILURE);
   }
   times = (TIMES*) malloc( sizeof(TIMES) );
   if (times == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc for TIMES.\n");
      exit(EXIT_FAILURE);
   }

   /* Print Arguments */
   fp = stdout;
   pad = 15;

   fprintf(fp, "=== ARGS ====================\n");
   fprintf(fp, "%*s:\t%s\n",     pad, "HMM_FILENAME",       hmm_file);
   fprintf(fp, "%*s:\t%s\n",     pad, "QUERY FILETYPE",     FILE_TYPE_NAMES[q_filetype]);
   fprintf(fp, "%*s:\t%s\n",     pad, "FASTA_FILENAME",     fasta_file);
   fprintf(fp, "%*s:\t%s\n",     pad, "TARGET FILETYPE",    FILE_TYPE_NAMES[t_filetype]);
   fprintf(fp, "%*s:\t%.3f\n",   pad, "ALPHA",              alpha);
   fprintf(fp, "%*s:\t%d\n",     pad, "BETA",               beta);
   fprintf(fp, "%*s:\t%s\n",     pad, "MODE",               MODE_NAMES[mode]);
   fprintf(fp, "=============================\n\n");

   if (fp != stdout) fclose(fp);

   // printf("=== LOAD HMM_PROFILE / QUERY -> START ===\n");
   CLOCK_Start(cl);

   target_prof = HMM_PROFILE_Parse(hmm_file, 0);
   T = target_prof->N;
   HMM_PROFILE_Config(target_prof, mode);
   // HMM_PROFILE_Dump(target_prof, stdout);

   query_seq = SEQUENCE_Fasta_Parse(fasta_file, 0);
   Q = query_seq->N;
   // SEQUENCE_Dump(query_seq, stdout);

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->load_hmm = (float)t;
   // printf("=== LOAD HMM_PROFILE / QUERY -> END ===\n\n");

   tot_cells = (Q + 1) * (T + 1);

   ALIGNMENT*  tr           = ALIGNMENT_Create();
   EDGEBOUNDS* edg_fwd_lin  = EDGEBOUNDS_Create();
   EDGEBOUNDS* edg_bck_lin  = EDGEBOUNDS_Create();
   EDGEBOUNDS* edg_row_lin  = NULL;
   EDGEBOUNDS* edg_diag_lin = NULL;

   /* allocate memory for quadratic algs (for DEBUGGING) */
   MATRIX_3D*  st_MATRIX    = MATRIX_3D_Create(NUM_NORMAL_STATES,  Q+1, T+1);
   float*      st_MX        = st_MATRIX->data;
   MATRIX_2D*  sp_MATRIX    = MATRIX_2D_Create(NUM_SPECIAL_STATES, Q+1);
   float*      sp_MX        = sp_MATRIX->data;
   /* allocate memory for linear algs */
   MATRIX_3D*  st_MATRIX3   = MATRIX_3D_Create(NUM_NORMAL_STATES,  Q+1, T+1);
   float*      st_MX3       = st_MATRIX3->data;

   // /* run viterbi algorithm */
   // printf("=== VITERBI -> START ===\n");
   CLOCK_Start(cl);

   viterbi_Quad(query_seq, target_prof, Q, T, st_MX, sp_MX, &sc);
   scores->viterbi_sc = sc;

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->viterbi = (float)t;
   // printf("=== VITERBI -> END ===\n\n");

   /* run traceback of viterbi */
   // printf("=== ALIGNMENT -> START ===\n");
   CLOCK_Start(cl);

   traceback_Build(query_seq, target_prof, Q, T, st_MX, sp_MX, tr);

   TRACE beg = tr->traces[tr->beg];
   TRACE end = tr->traces[tr->end];
   window_cells = (end.i - beg.i) * (end.j - beg.j);

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->traceback = (float)t;
   // printf("=== ALIGNMENT -> END ===\n\n");

   /* run forward/backward algorithms */
   init_Logsum();

   // printf("=== FORWARD/BACKWARD -> START ===\n");
   CLOCK_Start(cl);

   forward_Quad(query_seq, target_prof, Q, T, st_MX, sp_MX, &sc);
   // printf("Forward Score: %f\n", sc);
   scores->fwd_sc = sc;

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->fwd = (float)t;
   CLOCK_Start(cl);

   backward_Quad(query_seq, target_prof, Q, T, st_MX, sp_MX, &sc);
   // printf("Backward Score: %f\n", sc);
   scores->bck_sc = sc;

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->bck = (float)t;
   // printf("=== FORWARD/BACKWARD -> END ===\n\n");

   /* cloud forward (linear) */
   // printf("=== CLOUD FORWARD/BACKWARD (Linear) -> START ===\n");
   CLOCK_Start(cl);

   cloud_Forward_Linear(query_seq, target_prof, Q, T, st_MX, st_MX3, sp_MX, tr, edg_fwd_lin, alpha, beta, is_testing);
   // EDGEBOUNDS_Save(edg_fwd_lin, "output/my.cloud_fwd.lin.diags.pipe.edg");

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->cloud_fwd = (float)t;
   CLOCK_Start(cl);

   cloud_Backward_Linear(query_seq, target_prof, Q, T, st_MX, st_MX3, sp_MX, tr, edg_bck_lin, alpha, beta, is_testing);
   // EDGEBOUNDS_Save(edg_bck_lin, "output/my.cloud_bck.lin.diags.pipe.edg");

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->cloud_bck = (float)t;
   // printf("=== CLOUD FORWARD/BACKWARD (Linear) -> END ===\n\n");

   /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
   // printf("=== MERGE & REORIENT CLOUD (Linear) -> START ===\n");
   CLOCK_Start(cl);

   edg_diag_lin = EDGEBOUNDS_Merge(Q, T, edg_fwd_lin, edg_bck_lin);

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->merge = (float)t;
   CLOCK_Start(cl);

   edg_row_lin = EDGEBOUNDS_Reorient(Q, T, edg_diag_lin);

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->reorient = (float)t;
   // printf("=== MERGE & REORIENT CLOUD (Linear) -> END ===\n\n");

   // EDGEBOUNDS_Save(edg_row_lin, "output/my.pipeline.row.edg");
   // EDGEBOUNDS_Print(edg_row_lin);

   num_cells            = EDGEBOUNDS_Count(edg_row_lin);
   scores->perc_cells   = (float)num_cells/(float)tot_cells;
   scores->perc_window  = (float)num_cells/(float)window_cells;

   /* output number of cells computed */
   fp = stdout;
   pad = 15;

   fprintf(fp, "=== STATS ===================\n");
   fprintf(fp, "%*s:\t%d\n",           pad, "QUERY LENGTH",    Q);
   fprintf(fp, "%*s:\t%d\n",           pad, "TARGET LENGTH",   T);
   fprintf(fp, "%*s:\t(%d,%d)\n",      pad, "START WINDOW",    beg.i, beg.j);
   fprintf(fp, "%*s:\t(%d,%d)\n",      pad, "END WINDOW",      end.i, end.j);
   fprintf(fp, "%*s:\t%d/%d = %f\n",   pad, "TOTAL PERC",      num_cells, tot_cells, scores->perc_cells);
   fprintf(fp, "%*s:\t%d/%d = %f\n",   pad, "WINDOW PERC",     num_cells, window_cells, scores->perc_window);
   fprintf(fp, "=============================\n\n");

   if (fp != stdout) fclose(fp);
   
   // printf("=== BOUNDED FORWARD/BACKWARD (Linear) -> START ===\n");
   CLOCK_Start(cl);

   bound_Forward_Linear(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row_lin, is_testing, &sc);
   // printf("Bound Forward Score (Linear): %f\n", sc);
   scores->cloud_fwd_sc = sc;

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->bound_fwd = (float)t;
   CLOCK_Start(cl);

   bound_Backward_Linear(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row_lin, is_testing, &sc);
   // printf("Bound Backward Score (Linear): %f\n", sc);
   scores->cloud_bck_sc = sc;

   CLOCK_Stop(cl);
   t = CLOCK_Ticks(cl);
   times->bound_bck = (float)t;
   // printf("=== BOUNDED FORWARD/BACKWARD (Linear) -> END ===\n\n");

   /* total times */
   tot_t = 0;
   tot_t += times->load_hmm;
   tot_t += times->fwd + times->bck;
   tot_t += times->viterbi + times->traceback;
   tot_t += times->merge + times->reorient;
   tot_t += times->cloud_fwd + times->cloud_bck;
   tot_t += times->bound_fwd + times->bound_bck;

   /* output times */
   fp = stdout;
   pad = 15;

   fprintf(fp, "=== TIMES ===================\n");
   fprintf(fp, "%*s:\t%f\n", pad, "FWD",         times->fwd/tot_t );
   fprintf(fp, "%*s:\t%f\n", pad, "BCK",         times->bck/tot_t );
   fprintf(fp, "%*s:\t%f\n", pad, "VIT",         times->viterbi/tot_t );
   fprintf(fp, "%*s:\t%f\n", pad, "TRB",         times->traceback/tot_t );
   fprintf(fp, "%*s:\t%f\n", pad, "CLOUD_FWD",   times->cloud_fwd/tot_t );
   fprintf(fp, "%*s:\t%f\n", pad, "CLOUD_BCK",   times->cloud_bck/tot_t );
   fprintf(fp, "%*s:\t%f\n", pad, "MERGE",       times->merge/tot_t );
   fprintf(fp, "%*s:\t%f\n", pad, "REORIENT",    times->reorient/tot_t );
   fprintf(fp, "%*s:\t%f\n", pad, "BOUND_FWD",   times->bound_fwd/tot_t );
   fprintf(fp, "%*s:\t%f\n", pad, "BOUND_BCK",   times->bound_bck/tot_t );
   fprintf(fp, "\n");
   float full_fwdbck = times->fwd + times->bck;
   float full_cloud  = times->cloud_fwd + times->cloud_bck + 
                       times->merge + times->reorient + 
                       times->bound_fwd + times->bound_bck;
   fprintf(fp, "%*s:\t%f\n", pad, "FULL_FWDBCK",   full_fwdbck/tot_t );
   fprintf(fp, "%*s:\t%f\n", pad, "FULL_CLOUD",    full_cloud/tot_t );
   fprintf(fp, "\n");
   fprintf(fp, "%*s:\t%f\n", pad, "SPEEDUP",       full_fwdbck/full_cloud);
   fprintf(fp, "%*s:\t%f\n", pad, "TIME RATIO",    full_cloud/full_fwdbck);
   fprintf(fp, "============================\n\n");

   fprintf(fp, "=== SCORES =================\n");
   fprintf(fp, "%*s:\t%f\n", pad, "VITERBI",    scores->viterbi_sc );
   fprintf(fp, "%*s:\t%f\n", pad, "FWD",        scores->fwd_sc );
   fprintf(fp, "%*s:\t%f\n", pad, "BCK",        scores->bck_sc );
   fprintf(fp, "%*s:\t%f\n", pad, "BOUND_FWD",  scores->cloud_fwd_sc );
   fprintf(fp, "%*s:\t%f\n", pad, "BOUND_BCK",  scores->cloud_bck_sc );
   fprintf(fp, "============================\n\n");

   if (fp != stdout) fclose(fp);

   /* report stats */
   printf("Writing results to: '%s'\n", out_file);
   fp = fopen(out_file, "a+");

   if (fp == NULL) {
      fprintf(stderr, "ERROR: Cannot open file => %s\n", args->output_filepath);
      exit(EXIT_FAILURE);
   }
   fprintf(fp, "%s\t", hmm_file);
   fprintf(fp, "%s\t", fasta_file);
   fprintf(fp, "%f\t", scores->viterbi_sc);
   fprintf(fp, "%f\t", scores->fwd_sc);
   fprintf(fp, "%f\t", scores->bck_sc);
   fprintf(fp, "%f\t", scores->cloud_fwd_sc);
   fprintf(fp, "%f\t", scores->cloud_bck_sc);
   fprintf(fp, "%f\t", alpha);
   fprintf(fp, "%d\t", beta);
   fprintf(fp, "%f\t", scores->perc_cells);
   fprintf(fp, "%f\t", scores->perc_window);
   fprintf(fp, "%d\t", Q);
   fprintf(fp, "%d\n", T);

   if (fp != stdout) fclose(fp);

   /* free matrices */
   MATRIX_3D_Destroy(st_MATRIX);
   MATRIX_2D_Destroy(sp_MATRIX);
   MATRIX_3D_Destroy(st_MATRIX3);

   /* free hmm_profile and sequence*/
   HMM_PROFILE_Destroy(target_prof);
   SEQUENCE_Destroy(query_seq);

   /* free alignment */
   ALIGNMENT_Destroy(tr);

   /* free edgebounds */
   EDGEBOUNDS_Destroy(edg_fwd_lin);
   EDGEBOUNDS_Destroy(edg_bck_lin);
   EDGEBOUNDS_Destroy(edg_diag_lin);
   EDGEBOUNDS_Destroy(edg_row_lin);

   printf("...finished. \n");
}