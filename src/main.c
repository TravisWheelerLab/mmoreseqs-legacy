/*******************************************************************************
 *  @file main.c
 *  @brief Main Method, Argument Parser, Unit Tests. 
 *
 *  @author Dave Rich
 *  @bug Lots.
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
#include "structs.h"
#include "misc.h"
#include "hmm_parser.h"

/* objects */
#include "objects/edgebound.h"

/* viterbi & fwdbck (quadratic) */
#include "viterbi.h"
#include "forward_backward.h"

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

/* debugging methods */
#include "testing.h"

/* macros */
#define DEBUG false

void parse_args(int argc, char *argv[], ARGS *args);
void test(ARGS *args, char *hmm_file, char *fasta_file, float alpha, int beta);

/* MAIN */
int main (int argc, char *argv[])
{
   printf("begin main...\n");
   ARGS *args = (ARGS *)malloc( sizeof(ARGS) );

   char *hmm_file, *fasta_file;
   char *arg;
   int  i, j;

   /* DEFAULT TESTS */
   args->target_hmm_file = "data/test1_2.hmm";
   args->query_fasta_file = "data/test1_1.fa";
   args->alpha = 20.0f;
   args->beta = 5;
   args->search_mode = MODE_NAIVE;
   args->outfile_name = "stats/cloud_stats.tsv";
   args->test = false;
   // args->outfile = stdout;

   printf("parsing args...\n");
   parse_args(argc, argv, args);

   if (args->test) {
      printf("running test...\n");
      // test(args, args->target_hmm_file, args->query_fasta_file, args->alpha, args->beta);
   }
   else
   {
      printf("running main pipeline...\n");
      cloud_search_pipeline(args, args->target_hmm_file, args->query_fasta_file, args->alpha, args->beta);
   }
}


/* Parses Arguments from the command line */
void parse_args (int argc, char *argv[], ARGS *args)
{
   int i, len;
   int num_main_args, max_main_args;
   char *arg_tmp;
   
   num_main_args = 0;
   max_main_args = 2;

   for (i = 1; i < argc; ++i)
   {
      if ( argv[i][0] == '-' )
      {
         switch (argv[i][1]) {
            case 'a':   /* alpha value */
               i++;
               if (i < argc) {
                  args->alpha = atof(argv[i]);
               } else {
                  fprintf(stderr, "Error: -a flag requires argument.\n");
               }
               break;
            case 'b':   /* beta value */
               i++;
               if (i < argc) {
                  args->beta = atoi(argv[i]);
               } else {
                  fprintf(stderr, "Error: -b flag requires argument.\n");
               }
               break;
            case 'o':   /* append to outfile */
               i++;
               if (i < argc) {
                  args->outfile_name = argv[i];
                  // args->outfile = fopen(argv[i], "w+");
               } else {
                  fprintf(stderr, "Error: -o flag requires argument.\n");
               } 
               break;
            // case 'q':   /* window range for query */
            //    i++;
            //    if (i+1 < argc) {
            //       args->range_query = malloc( sizeof(RANGE) );
            //       args->range_query->st = atoi(argv[i]);
            //       i++;
            //       args->range_query->end = atoi(argv[i]);
            //    } else {
            //       fprintf(stderr, "Error: -q flag requires two arguments.\n");
            //    }
            //    break;
            // case 't':   /* window range for target */
            //    i++;
            //    if (i+1 < argc) {
            //       args->range_target = malloc( sizeof(RANGE) );
            //       args->range_target->st = atoi(argv[i]);
            //       i++;
            //       args->range_target->end = atoi(argv[i]);
            //    } else {
            //       fprintf(stderr, "Error: -t flag requires two arguments.\n");
            //    } 
            //    break;
            // case 'm':   /* mode */
            //    i++;
            //    if (i < argc) {
            //       arg_tmp = argv[i];
            //    } else {
            //       fprintf(stderr, "Error: -m flag requires argument.\n");
            //    }  
            //    break; 
            case '+':   /* run test */
               i++;
               args->test = true; 
               break;
            default:
               fprintf(stderr, "Error: -%c is not a valid flag.\n", argv[i][1]);
               exit(1);
               break;
         }
      }
      else
      {
         len = get_str_len(argv[i]);
         if (num_main_args == 0)
         {
            args->target_hmm_file = argv[i];
         }
         else if (num_main_args == 1)
         {
            args->query_fasta_file = argv[i];
         }
         else
         {
            fprintf(stderr, "Error: Too many main arguments.\n");
         }
         num_main_args++;
      }
   }
}


/* unit test */
void test(ARGS *args, char *hmm_file, char *fasta_file, float alpha, int beta)
{
   /* PRINT ARGS */
   printf("HMM_FILENAME: %s\n", hmm_file);
   printf("FASTA_FILENAME: %s\n", fasta_file);
   printf("ALPHA: %f\n", alpha);
   printf("BETA: %d\n\n", beta);
   SCORES *scores = (SCORES*)malloc( sizeof(SCORES) );

   printf("building hmm profile...\n");
   /* get target profile */
   HMM_PROFILE *target_prof = (HMM_PROFILE *)malloc( sizeof(HMM_PROFILE) );
   hmmprofile_Create(target_prof, hmm_file);
   // hmmprofile_Display(target_prof);
   hmmprofile_Save(target_prof, "output/my.pre-profile.tsv");

   printf("configuring...\n");
   /* mode choices: MODE_NONE, MODE_MULTILOCAL, MODE_MULTIGLOCAL, MODE_UNILOCAL, MODE_UNIGLOCAL */
   char* modes[] = { "None", "Multi-local", "Multi-glocal", "Uni-local", "Uni-glocal" }; 
   int mode; 
   // mode = MODE_MULTILOCAL;    /* HMMER standard mode (allows jumps) */
   mode = MODE_UNILOCAL;   /* Cloud mode (prohibiits jumps) */
   printf("MODE: %s\n", modes[mode]);
   hmmprofile_Config(target_prof, mode);
   // hmmprofile_Display(target_prof);
   hmmprofile_Save(target_prof, "output/my.post-profile.tsv");
   int T = target_prof->leng;

   printf("building query sequence...\n");
   /* get query sequence */
   SEQ *query_seq = (SEQ *)malloc( sizeof(SEQ) );
   seq_Create(query_seq, fasta_file);
   seq_Display(query_seq);
   int Q = query_seq->leng;

   printf("I/O was successful!\n\n");

   /* allocate memory to store results */
   float sc, perc_cells;
   int num_cells, window_cells; 
   int tot_cells = (Q + 1) * (T + 1);
   TRACEBACK *tr = (TRACEBACK *)malloc( sizeof(TRACEBACK) );
   EDGEBOUNDS *edg_fwd = edgebounds_Create();
   EDGEBOUNDS *edg_bck = edgebounds_Create();
   EDGEBOUNDS *edg = edgebounds_Create();
   EDGEBOUNDS *edg_tmp = edgebounds_Create();
   EDGEBOUNDS *edg_fwd_tmp = edgebounds_Create();
   EDGEBOUNDS *edg_bck_tmp = edgebounds_Create();
   EDGEBOUNDS *edg_row = edgebounds_Create();
   EDGEBOUNDS *edg_quad_row = edgebounds_Create();
   EDGEBOUNDS *edg_quad_diag = edgebounds_Create();
   EDGEBOUNDS *edg_lin_row = edgebounds_Create();
   EDGEBOUNDS *edg_lin_diag = edgebounds_Create();

   /* allocate memory for quadratic algs (for DEBUGGING) */
   float *st_MX = (float *) malloc( sizeof(float) * (NUM_NORMAL_STATES * (Q+1) * (T+1)) );
   float *sp_MX = (float *) malloc( sizeof(float) * (NUM_SPECIAL_STATES * (Q+1)) );
   /* allocate memory for comparing quadratic algs (for DEBUGGING) */
   float *st_MX_tmp = (float *) malloc( sizeof(float) * (NUM_NORMAL_STATES * (Q+1) * (T+1)) );
   /* allocate memory for cloud matrices (for DEBUGGING) */
   float *st_MX_cloud = (float *) malloc( sizeof(float) * (NUM_NORMAL_STATES * (Q+1) * (T+1)) );
   float *sp_MX_cloud = (float *) malloc( sizeof(float) * (NUM_SPECIAL_STATES * (Q+1)) );
   /* allocate memory for linear algs */
   float *st_MX3 = (float *) malloc( sizeof(float) * (NUM_NORMAL_STATES * ((T+1)+(Q+1)) * 3) );

   bool test;

   /* run viterbi algorithm */
   printf("=== VITERBI -> START ===\n");
   sc = viterbi_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, tr);
   printf("Viterbi Score: %f\n", sc);
   scores->viterbi_sc = sc;
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_Save(Q, T, st_MX, sp_MX, "output/my.viterbi.mx");
   printf("=== VITERBI -> END ===\n\n");

   /* run traceback of viterbi */
   printf("=== TRACEBACK -> START ===\n");
   traceback_Build(query_seq, target_prof, Q, T, st_MX, sp_MX, tr);
   // traceback_Print(tr);
   traceback_Save(tr, "output/my.traceback_list.tsv");
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   traceback_Show(Q, T, st_MX, sp_MX, tr);
   printf("START: (%d,%d) -> END: (%d,%d)\n", tr->first_m.i, tr->first_m.j, tr->last_m.i, tr->last_m.j);
   window_cells = (tr->last_m.i - tr->first_m.i) * (tr->last_m.j - tr->first_m.j);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.traceback.mx");
   printf("=== TRACEBACK -> END ===\n\n");

   /* run forward/backward algorithms */
   init_Logsum();

   printf("=== FORWARD -> START ===\n");
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   sc = forward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX);
   printf("Forward Score: %f\n", sc);
   scores->fwd_sc = sc;
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.fwd.mx");
   printf("=== FORWARD -> END ===\n\n");

   printf("=== BACKWARD -> START ===\n");
   sc = backward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX);
   printf("Backward Score: %f\n", sc);
   scores->bck_sc = sc;
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bck.mx");
   printf("=== BACKWARD -> END ===\n\n");

   /* TEST */
   printf("=== TEST -> START ===\n");
   fwd_test_cycle(Q, T, st_MX, sp_MX, tr);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.test_fwd.quad.mx");
   bck_test_cycle(Q, T, st_MX, sp_MX, tr);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.test_bck.quad.mx");
   fwd_test_cycle3(Q, T, st_MX, st_MX3, sp_MX, tr);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.test_fwd.lin.mx");
   bck_test_cycle3(Q, T, st_MX, st_MX3, sp_MX, tr);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.test_bck.lin.mx");
   printf("=== TEST -> END ===\n\n");

   /* cloud forward (quadratic) */
   printf("=== CLOUD FORWARD (Quadratic) -> START ===\n");
   cloud_forward_Run(query_seq, target_prof, Q, T, st_MX_tmp, sp_MX, tr, edg_fwd_tmp, alpha, beta);
   dp_matrix_trace_Save(Q, T, st_MX_tmp, sp_MX, tr, "output/my.cloud_fwd_vals.quad.mx");
   dp_matrix_Clear_X(Q, T, st_MX_tmp, sp_MX, 0);
   cloud_Fill(Q, T, st_MX_tmp, sp_MX, edg_fwd_tmp, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX_tmp, sp_MX, tr, "output/my.cloud_fwd.quad.mx");
   edgebounds_Save(edg_fwd_tmp, "output/my.cloud_fwd.quad.diags.edg");
   printf("=== CLOUD FORWARD (Quadratic) -> END ===\n\n");

   /* cloud forward (linear) */
   printf("=== CLOUD FORWARD (Linear) -> START ===\n");
   cloud_forward_Run3(query_seq, target_prof, Q, T, st_MX, st_MX3, sp_MX, tr, edg_fwd, alpha, beta, true);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud_fwd_vals.lin.mx");
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_fwd, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud_fwd.lin.mx");
   edgebounds_Save(edg_fwd, "output/my.cloud_fwd.lin.diags.edg");
   test = dp_matrix_Compare(Q, T, st_MX, sp_MX, st_MX_tmp, sp_MX);
   printf("Cloud Forward - Lin vs Quad? %s\n", test ? "TRUE" : "FALSE" );
   printf("=== CLOUD FORWARD (Linear) -> END ===\n\n");

   /* cloud backward (quadratic) */
   printf("=== CLOUD BACKWARD (Quadratic) -> START ===\n");
   cloud_backward_Run(query_seq, target_prof, Q, T, st_MX_tmp, sp_MX, tr, edg_bck_tmp, alpha, beta);
   dp_matrix_trace_Save(Q, T, st_MX_tmp, sp_MX, tr, "output/my.cloud_bck_vals.quad.mx");
   dp_matrix_Clear_X(Q, T, st_MX_tmp, sp_MX, 0);
   cloud_Fill(Q, T, st_MX_tmp, sp_MX, edg_bck_tmp, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX_tmp, sp_MX, tr, "output/my.cloud_bck.quad.mx");
   edgebounds_Save(edg_bck_tmp, "output/my.cloud_bck.quad.diags.edg");
   printf("=== CLOUD BACKWARD (Quadratic) -> END ===\n\n");

   /* cloud backward (linear) */
   printf("=== CLOUD BACKWARD (Linear) -> START ===\n");
   cloud_backward_Run3(query_seq, target_prof, Q, T, st_MX, st_MX3, sp_MX, tr, edg_bck, alpha, beta, true);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud_bck_vals.lin.mx");
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_bck, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud_bck.lin.mx");
   edgebounds_Save(edg_bck, "output/my.cloud_bck.lin.diags.edg");
   test = dp_matrix_Compare(Q, T, st_MX, sp_MX, st_MX_tmp, sp_MX);
   printf("Cloud Backward - Lin vs Quad? %s\n", test ? "TRUE" : "FALSE" );
   printf("=== CLOUD BACKWARD (Linear) -> END ===\n\n");

   /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
   printf("=== MERGE & REORIENT CLOUD (Naive) -> START ===\n");
   edgebounds_Merge_Reorient_Naive(edg_fwd, edg_bck, edg_quad_diag, edg_quad_row, Q, T, st_MX, sp_MX);
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_quad_diag, 1, MODE_DIAG);
   edgebounds_Save(edg_quad_diag, "output/my.cloud.quad.diags.edg");
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud.quad.diags.mx");
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_quad_row, 1, MODE_ROW);
   edgebounds_Save(edg_quad_row, "output/my.cloud.quad.rows.edg");
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud.quad.rows.mx");
   dp_matrix_Copy(Q, T, st_MX, sp_MX, st_MX_cloud, sp_MX_cloud);
   printf("=== MERGE & REORIENT CLOUD (Naive) -> END ===\n\n");

   /* stats */
   num_cells = cloud_Cell_Count(Q, T, st_MX, sp_MX);
   scores->perc_cells = (float)num_cells/(float)tot_cells;
   printf("Perc. Total Cells Computed = %d/%d = %f\n", num_cells, tot_cells, scores->perc_cells);
   scores->perc_window = (float)num_cells/(float)window_cells;
   printf("Perc. Window Cells Computed = %d/%d = %f\n", num_cells, window_cells, scores->perc_window);

   /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
   printf("=== MERGE & REORIENT CLOUD (Linear) -> START ===\n");
   edgebounds_Merge(Q, T, edg_fwd, edg_bck, edg_lin_diag);
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_lin_diag, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud.lin.diags.mx");
   edgebounds_Save(edg_lin_diag, "output/my.cloud.lin.diags.edg");
   edgebounds_Reorient(Q, T, edg_lin_diag, edg_lin_row);
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_lin_row, 1, MODE_ROW);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud.lin.rows.mx");
   edgebounds_Save(edg_lin_row, "output/my.cloud.lin.rows.edg");
   printf("=== MERGE & REORIENT CLOUD (Linear) -> END ===\n\n");

   /* create cloud that covers entire matrix (full fwd/bck) */
   // printf("=== TEST CLOUD -> START ===\n");
   // dp_matrix_Clear_X(Q, T, st_MX_cloud, sp_MX_cloud, 1);
   // edgebounds_Build_From_Cloud(edg, Q, T, st_MX_cloud, MODE_ROW);
   // edgebounds_Save(edg, "output/my.full_fwdbck.edg");
   // dp_matrix_trace_Save(Q, T, st_MX_cloud, sp_MX_cloud, tr, "output/my.cloud.full_fwdbck.mx");
   // printf("=== TEST CLOUD -> END ===\n\n");

   /* bounded forward */
   printf("=== BOUNDED FORWARD -> START ===\n");
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   dp_matrix_Clear_X(Q, T, st_MX_cloud, sp_MX_cloud, 0);
   cloud_Fill(Q, T, st_MX_cloud, sp_MX_cloud, edg_quad_row, 1, MODE_ROW);
   forward_Bounded_Naive_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, st_MX_cloud, &sc);
   printf("Bounded Forward Score (Naive): %f\n", sc);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_fwd.naive.mx");
   
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   forward_bounded_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, edg_quad_row, &sc);
   printf("Bounded Forward Score (Quadratic): %f\n", sc);
   scores->cloud_fwd_sc = sc;
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_fwd.quad.mx");

   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   forward_bounded_Run3(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_quad_row, &sc, true);
   printf("Bound Forward Score (Linear): %f\n", sc);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_fwd.lin.mx");
   printf("=== BOUNDED FORWARD -> END ===\n\n");

   /* bounded backward */
   printf("=== BOUNDED BACKWARD -> START ===\n");
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   backward_Bounded_Naive_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, st_MX_cloud, &sc);
   printf("Bounded Backward Score (Naive): %f\n", sc);
   scores->cloud_bck_sc = sc;
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_bck.naive.mx");

   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   backward_bounded_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, edg_quad_row, &sc);
   printf("Bounded Backward Score (Quadratic): %f\n", sc);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_bck.quad.mx");

   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   backward_bounded_Run3(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_quad_row, &sc, true);
   printf("Bound Backward Score (Linear): %f\n", sc);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_bck.lin.mx");
   printf("=== BOUNDED BACKWARD -> END ===\n\n");

   /* sample stats */
   char *fileout = args->outfile_name;
   printf("Writing results to: '%s'\n", fileout);
   FILE *fp = fopen(fileout, "a+");
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
   fprintf(fp, "%f\n", scores->perc_window);
   fclose(fp);

   /* TODO: fix free memory! */
   free(tr);
   free(edg_fwd);
   free(edg_bck);
   free(edg);
   free(st_MX);
   free(sp_MX);
   free(st_MX_cloud);
   free(sp_MX_cloud);
   free(st_MX3);

   printf("...test finished. \n");
}


/* standard pipeline */
void cloud_search_pipeline(ARGS *args, char *hmm_file, char *fasta_file, float alpha, int beta) 
{
   /* PRINT ARGS */
   printf("HMM_FILENAME: %s\n", hmm_file);
   printf("FASTA_FILENAME: %s\n", fasta_file);
   printf("ALPHA: %f\n", alpha);
   printf("BETA: %d\n\n", beta);
   SCORES *scores = (SCORES*)malloc( sizeof(SCORES) );
   float sc;
   clock_t time_st, time_end, duration;

   HMM_PROFILE *target_prof = (HMM_PROFILE *)malloc( sizeof(HMM_PROFILE) );
   SEQ *query_seq = (SEQ *)malloc( sizeof(SEQ) );
   TRACEBACK *tr = (TRACEBACK *)malloc( sizeof(TRACEBACK) );

   printf("building hmm profile...\n");
   int mode = MODE_UNILOCAL;   /* Cloud mode (prohibiits jumps) */
   hmmprofile_Create(target_prof, hmm_file);
   hmmprofile_Config(target_prof, mode);
   int T = target_prof->leng;

   printf("building query sequence...\n");
   /* get query sequence */
   seq_Create(query_seq, fasta_file);
   seq_Display(query_seq);
   int Q = query_seq->leng;

   float *st_MX = (float *) malloc( sizeof(float) * (NUM_NORMAL_STATES * (Q+1) * (T+1)) );
   float *sp_MX = (float *) malloc( sizeof(float) * (NUM_SPECIAL_STATES * (Q+1)) );
   float *st_MX3 = (float *) malloc( sizeof(float) * (NUM_NORMAL_STATES * ((T+1)+(Q+1)) * 3) );
   EDGEBOUNDS *edg_fwd = edgebounds_Create();
   EDGEBOUNDS *edg_bck = edgebounds_Create();
   EDGEBOUNDS *edg_diag = edgebounds_Create();
   EDGEBOUNDS *edg_row = edgebounds_Create();

   /* run viterbi algorithm */
   printf("=== VITERBI/TRACE -> START ===\n");
   time_st = clock();
   sc = viterbi_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, tr);
   scores->viterbi_sc = sc;
   traceback_Build(query_seq, target_prof, Q, T, st_MX, sp_MX, tr);
   time_end = clock();
   duration = (time_end - time_st) * 1000 / CLOCKS_PER_SEC;
   printf("VITERBI/TRACE took %d ms\n", duration);
   printf("=== VITERBI/TRACE -> END ===\n\n");

   /* cloud forward/backward (linear) */
   printf("=== CLOUD FORWARD/BACKWARD (Linear) -> START ===\n");
   time_st = clock();
   sc = viterbi_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, tr);
   cloud_forward_Run3(query_seq, target_prof, Q, T, st_MX, st_MX3, sp_MX, tr, edg_fwd, alpha, beta, false);
   cloud_backward_Run3(query_seq, target_prof, Q, T, st_MX, st_MX3, sp_MX, tr, edg_bck, alpha, beta, false);
   time_end = clock();
   duration = (time_end - time_st) * 1000 / CLOCKS_PER_SEC;
   printf("CLOUD FWD/BCK took %d ms\n", duration);
   printf("=== CLOUD FORWARD/BACKWARD (Linear) -> END ===\n\n");

   /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
   printf("=== MERGE & REORIENT CLOUD (Linear) -> START ===\n");
   time_st = clock();
   edgebounds_Merge(Q, T, edg_fwd, edg_bck, edg_diag);
   edgebounds_Reorient(Q, T, edg_diag, edg_row);
   time_end = clock();
   duration = (time_end - time_st) * 1000 / CLOCKS_PER_SEC;
   printf("MERGE/REORIENT took %d ms\n", duration);
   printf("=== MERGE & REORIENT CLOUD (Linear) -> END ===\n\n");

   /* bounded forward/backward */
   printf("=== BOUNDED FORWARD/BACKWARD (Linear) -> START ===\n"); 
   time_st = clock();
   forward_bounded_Run3(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row, &sc, false);
   scores->cloud_fwd_sc = sc;
   backward_bounded_Run3(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row, &sc, false);
   scores->cloud_bck_sc = sc;
   time_end = clock();
   duration = (time_end - time_st) * 1000 / CLOCKS_PER_SEC;
   printf("MERGE/REORIENT took %d ms\n", duration);
   printf("=== BOUNDED FORWARD/BACKWARD (Linear) -> END ===\n\n");   

   printf("Cloud Forward: %f\n", scores->cloud_fwd_sc);
   printf("Cloud Backward: %f\n", scores->cloud_bck_sc);

   edgebounds_Destroy(edg_fwd);
   edgebounds_Destroy(edg_bck);
   edgebounds_Destroy(edg_diag);
   edgebounds_Destroy(edg_row);
   free(st_MX);
   free(st_MX3);
   free(sp_MX);
   free(tr);
   free(scores);
}