/*******************************************************************************
 *  FILE:      main.c
 *  PURPOSE:   Main Method, Cloud Fwd/Bck Pipeline, Argument Parser, Unit Tests
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
#include "hmm_parser.h"

/* objects */
#include "objects/alignment.h"
#include "objects/sequence.h"
#include "objects/hmm_profile.h"
#include "objects/edgebound.h"
#include "objects/clock.h"
#include "objects/matrix_2d.h"
#include "objects/matrix_3d.h"
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

/* debug macros */
#ifndef DEBUG
#define DEBUG false
#endif

/* === OBJECTS === */
typedef struct {
   float  alpha; 
   int    beta;
   int    search_mode;
   char*  target_hmm_file;
   char*  query_fasta_file;
   RANGE* range_query; 
   RANGE* range_target;
   char*  outfile_name;
   int    test;
   FILE*  outfile;
} ARGS;

typedef struct {
   float load_hmm;
   float viterbi, traceback;
   float fwd, bck;
   float cloud_fwd, cloud_bck;
   float merge, reorient;
   float bound_fwd, bound_bck;
} TIMES;

/* === HEADERS === */
void parse_args(int argc, char *argv[], ARGS *args);
void test_pipeline(ARGS *args, char *hmm_file, char *fasta_file, float alpha, int beta);
void cloud_search_pipeline(ARGS *args, char *hmm_file, char *fasta_file, float alpha, int beta);


/* === MAIN === */
int main (int argc, char *argv[])
{
   printf("begin main...\n");
   ARGS *args = (ARGS *)malloc( sizeof(ARGS) );

   char *hmm_file, *fasta_file;
   char *arg;
   int  i, j;

   /* DEFAULT ARGUMENTS */
   args->target_hmm_file  = "data/test1_2.hmm";
   args->query_fasta_file = "data/test1_1.fa";
   args->alpha            = 20.0f;
   args->beta             = 5;
   args->search_mode      = MODE_NAIVE;
   args->outfile_name     = "stats/cloud_stats.tsv";
   args->test             = false;
   // args->outfile          = stdout;

   printf("parsing args...\n");
   parse_args(argc, argv, args);

   if (args->test) {
      printf("running test...\n\n");
      test_pipeline(args, args->target_hmm_file, args->query_fasta_file, args->alpha, args->beta);
   }
   else
   {
      printf("running main pipeline...\n\n");
      cloud_search_pipeline(args, args->target_hmm_file, args->query_fasta_file, args->alpha, args->beta);
   }
   exit(EXIT_SUCCESS);
}


/* Parses Arguments from the command line */
void parse_args (int argc, char *argv[], ARGS *args)
{
   int i, len;
   int num_main_args, max_main_args;
   char *arg_tmp;
   
   num_main_args = 0;
   max_main_args = 2;

   if (argc == 1) {
      printf("Usage: <hmm_file> <fasta_file>\n");
      printf("Opts: -a <alpha>, -b <beta>, -o <output_file>, -T\n");
   }

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
            case 'T':   /* run test */
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
         len = strlen(argv[i]) + 1;
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


/* pipeline run optimized and unoptimized versions of search algs */
void test_pipeline(ARGS *args, char *hmm_file, char *fasta_file, float alpha, int beta)
{
   /* PROFILE & SEQUENCE */
   HMM_PROFILE*   target_prof  = NULL;
   SEQUENCE*      query_seq    = NULL;

   /* EDGEBOUNDS & ALIGNMENTS */
   ALIGNMENT  *tr              = ALIGNMENT_Create();

   EDGEBOUNDS *edg_fwd_lin     = EDGEBOUNDS_Create();
   EDGEBOUNDS *edg_bck_lin     = EDGEBOUNDS_Create();
   EDGEBOUNDS *edg_row_lin     = EDGEBOUNDS_Create();
   EDGEBOUNDS *edg_diag_lin    = EDGEBOUNDS_Create();

   EDGEBOUNDS *edg_fwd_quad    = EDGEBOUNDS_Create();
   EDGEBOUNDS *edg_bck_quad    = EDGEBOUNDS_Create();
   EDGEBOUNDS *edg_row_quad    = EDGEBOUNDS_Create();
   EDGEBOUNDS *edg_diag_quad   = EDGEBOUNDS_Create();

   /* SCORES => stores result scores */
   SCORES*        scores       = (SCORES*)malloc( sizeof(SCORES) );

   /* TEST => embed into quadratic matrix */
   bool           is_testing   = true;
   int            test         = 0;

   /* records percentage of space computed */
   int            T            = 0;
   int            Q            = 0;
   float          sc           = 0.f; 
   float          perc_cells   = 0.f;
   int            num_cells    = 0;  
   int            window_cells = 0; 
   int            tot_cells    = (Q + 1) * (T + 1);

   /* MODE => determines whether
   /* mode choices: MODE_NONE, MODE_MULTILOCAL, MODE_MULTIGLOCAL, MODE_UNILOCAL, MODE_UNIGLOCAL */
   char* modes[] = { "None", "Multi-local", "Multi-glocal", "Uni-local", "Uni-glocal" }; 
   int   mode; 
         // mode = MODE_MULTILOCAL;    /* HMMER standard mode (allows jumps) */
         mode = MODE_UNILOCAL;      /* Cloud mode (prohibiits jumps) */
   printf("MODE: %s\n", modes[mode]);

   /* PRINT ARGS */
   printf("HMM_FILENAME: %s\n", hmm_file);
   printf("FASTA_FILENAME: %s\n", fasta_file);
   printf("ALPHA: %f\n", alpha);
   printf("BETA: %d\n\n", beta);

   printf("=== BUILD HMM_PROFILE / QUERY -> START ===\n");
   
   /* build target profile */
   printf("building hmm profile...\n");
   target_prof = HMM_PROFILE_Parse(hmm_file);
   HMM_PROFILE_Config(target_prof, mode);
   HMM_PROFILE_Dump(target_prof, fopen("output/my.post-profile.tsv", "w") );
   T = target_prof->N;

   /* build query sequence */
   printf("building query sequence...\n");
   query_seq = SEQUENCE_Parse(fasta_file);
   Q = query_seq->N;
   printf("=== BUILD HMM_PROFILE / QUERY -> END ===\n");

   /* allocate memory for quadratic algs (for DEBUGGING) */
   MATRIX_3D* st_MATRIX       = MATRIX_3D_Create(NUM_NORMAL_STATES,  Q+1, T+1);
   float*     st_MX           = st_MATRIX->data;
   MATRIX_2D* sp_MATRIX       = MATRIX_2D_Create(NUM_SPECIAL_STATES, Q+1);
   float*     sp_MX           = sp_MATRIX->data;
   /* allocate memory for comparing quadratic algs (for DEBUGGING) */
   MATRIX_3D* st_MATRIX_tmp   = MATRIX_3D_Create(NUM_NORMAL_STATES,  Q+1, T+1);
   float*     st_MX_tmp       = st_MATRIX_tmp->data;
   /* allocate memory for cloud matrices (for DEBUGGING) */
   MATRIX_3D* st_MATRIX_cloud = MATRIX_3D_Create(NUM_NORMAL_STATES,  Q+1, T+1);
   float*     st_MX_cloud     = st_MATRIX_cloud->data;
   MATRIX_2D* sp_MATRIX_cloud = MATRIX_2D_Create(NUM_SPECIAL_STATES, Q+1);
   float*     sp_MX_cloud     = sp_MATRIX_cloud->data;
   /* allocate memory for linear algs */
   MATRIX_3D* st_MATRIX3      = MATRIX_3D_Create(NUM_NORMAL_STATES,  Q+1, T+1);
   float*     st_MX3          = st_MATRIX3->data;

   /* run viterbi algorithm */
   printf("=== VITERBI -> START ===\n");
   sc = viterbi_Quad(query_seq, target_prof, Q, T, st_MX, sp_MX, &sc);
   printf("Viterbi Score: %f\n", sc);
   scores->viterbi_sc = sc;
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_Save(Q, T, st_MX, sp_MX, "output/my.viterbi.mx");
   printf("=== VITERBI -> END ===\n\n");

   /* run traceback of viterbi */
   printf("=== TRACEBACK -> START ===\n");
   traceback_Build(query_seq, target_prof, Q, T, st_MX, sp_MX, tr);
   // ALIGNMENT_Dump(tr, stdout);
   ALIGNMENT_Save(tr, "output/my.traceback.tsv");
   TRACE beg = tr->traces[tr->beg];
   TRACE end = tr->traces[tr->end];
   printf("START: (%d,%d) -> END: (%d,%d)\n", beg.i, beg.j, end.i, end.j);
   window_cells = (end.i - beg.i) * (end.j - beg.j);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.traceback.mx");
   printf("=== TRACEBACK -> END ===\n\n");

   /* run forward/backward algorithms */
   init_Logsum();

   printf("=== FORWARD -> START ===\n");
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   sc = forward_Quad(query_seq, target_prof, Q, T, st_MX, sp_MX, &sc);
   printf("Forward Score: %f\n", sc);
   scores->fwd_sc = sc;
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.fwd.mx");
   printf("=== FORWARD -> END ===\n\n");

   printf("=== BACKWARD -> START ===\n");
   sc = backward_Quad(query_seq, target_prof, Q, T, st_MX, sp_MX, &sc);
   printf("Backward Score: %f\n", sc);
   scores->bck_sc = sc;
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bck.mx");
   printf("=== BACKWARD -> END ===\n\n");

   // /* TEST */
   // printf("=== TEST -> START ===\n");
   // fwd_test_cycle(Q, T, st_MX, sp_MX, tr);
   // dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.test_fwd.quad.mx");
   // bck_test_cycle(Q, T, st_MX, sp_MX, tr);
   // dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.test_bck.quad.mx");
   // fwd_test_cycle3(Q, T, st_MX, st_MX3, sp_MX, tr);
   // dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.test_fwd.lin.mx");
   // bck_test_cycle3(Q, T, st_MX, st_MX3, sp_MX, tr);
   // dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.test_bck.lin.mx");
   // printf("=== TEST -> END ===\n\n");

   /* cloud forward (quadratic) */
   printf("=== CLOUD FORWARD (Quadratic) -> START ===\n");
   cloud_Forward_Quad(query_seq, target_prof, Q, T, st_MX_tmp, sp_MX, tr, edg_fwd_quad, alpha, beta);
   dp_matrix_trace_Save(Q, T, st_MX_tmp, sp_MX, tr, "output/my.cloud_fwd_vals.quad.mx");
   dp_matrix_Clear_X(Q, T, st_MX_tmp, sp_MX, 0);
   cloud_Fill(Q, T, st_MX_tmp, sp_MX, edg_fwd_quad, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX_tmp, sp_MX, tr, "output/my.cloud_fwd.quad.mx");
   EDGEBOUNDS_Save(edg_fwd_quad, "output/my.cloud_fwd.quad.diags.edg");
   printf("=== CLOUD FORWARD (Quadratic) -> END ===\n\n");

   /* cloud forward (linear) */
   printf("=== CLOUD FORWARD (Linear) -> START ===\n");
   cloud_Forward_Linear(query_seq, target_prof, Q, T, st_MX, st_MX3, sp_MX, tr, edg_fwd_lin, alpha, beta, is_testing);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud_fwd_vals.lin.mx");
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_fwd_lin, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud_fwd.lin.mx");
   EDGEBOUNDS_Save(edg_fwd_lin, "output/my.cloud_fwd.lin.diags.edg");
   test = dp_matrix_Compare(Q, T, st_MX, sp_MX, st_MX_tmp, sp_MX);
   printf("Cloud Forward - Lin vs Quad? %s\n", test ? "TRUE" : "FALSE" );
   printf("=== CLOUD FORWARD (Linear) -> END ===\n\n");

   /* cloud backward (quadratic) */
   printf("=== CLOUD BACKWARD (Quadratic) -> START ===\n");
   cloud_Backward_Quad(query_seq, target_prof, Q, T, st_MX_tmp, sp_MX, tr, edg_bck_quad, alpha, beta);
   dp_matrix_trace_Save(Q, T, st_MX_tmp, sp_MX, tr, "output/my.cloud_bck_vals.quad.mx");
   dp_matrix_Clear_X(Q, T, st_MX_tmp, sp_MX, 0);
   cloud_Fill(Q, T, st_MX_tmp, sp_MX, edg_bck_quad, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX_tmp, sp_MX, tr, "output/my.cloud_bck.quad.mx");
   EDGEBOUNDS_Save(edg_bck_quad, "output/my.cloud_bck.quad.diags.edg");
   printf("=== CLOUD BACKWARD (Quadratic) -> END ===\n\n");

   /* cloud backward (linear) */
   printf("=== CLOUD BACKWARD (Linear) -> START ===\n");
   cloud_Backward_Linear(query_seq, target_prof, Q, T, st_MX, st_MX3, sp_MX, tr, edg_bck_lin, alpha, beta, is_testing);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud_bck_vals.lin.mx");
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_bck_lin, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud_bck.lin.mx");
   EDGEBOUNDS_Save(edg_bck_lin, "output/my.cloud_bck.lin.diags.edg");
   test = dp_matrix_Compare(Q, T, st_MX, sp_MX, st_MX_tmp, sp_MX);
   printf("Cloud Backward - Lin vs Quad? %s\n", test ? "TRUE" : "FALSE" );
   printf("=== CLOUD BACKWARD (Linear) -> END ===\n\n");

   /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
   printf("=== MERGE & REORIENT CLOUD (Naive) -> START ===\n");
   EDGEBOUNDS_Merge_Reorient_Naive(edg_fwd_quad, edg_bck_quad, edg_diag_quad, edg_row_quad, Q, T, st_MX, sp_MX);
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_diag_quad, 1, MODE_DIAG);
   EDGEBOUNDS_Save(edg_diag_quad, "output/my.cloud.quad.diags.edg");
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud.quad.diags.mx");
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_row_quad, 1, MODE_ROW);
   EDGEBOUNDS_Save(edg_row_quad, "output/my.cloud.quad.rows.edg");
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
   edg_diag_lin = EDGEBOUNDS_Merge(Q, T, edg_fwd_lin, edg_bck_lin);
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_diag_lin, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud.lin.diags.mx");
   EDGEBOUNDS_Save(edg_diag_lin, "output/my.cloud.lin.diags.edg");
   edg_row_lin = EDGEBOUNDS_Reorient(Q, T, edg_diag_lin);
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_row_lin, 1, MODE_ROW);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.cloud.lin.rows.mx");
   EDGEBOUNDS_Save(edg_row_lin, "output/my.cloud.lin.rows.edg");
   printf("=== MERGE & REORIENT CLOUD (Linear) -> END ===\n\n");

   // /* create cloud that covers entire matrix (full fwd/bck) */
   // // printf("=== TEST CLOUD -> START ===\n");
   // // dp_matrix_Clear_X(Q, T, st_MX_cloud, sp_MX_cloud, 1);
   // // EDGEBOUNDS_Build_From_Cloud(edg, Q, T, st_MX_cloud, MODE_ROW);
   // // EDGEBOUNDS_Save(edg, "output/my.full_fwdbck.edg");
   // // dp_matrix_trace_Save(Q, T, st_MX_cloud, sp_MX_cloud, tr, "output/my.cloud.full_fwdbck.mx");
   // // printf("=== TEST CLOUD -> END ===\n\n");

   /* fill cloud for naive algs */
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   dp_matrix_Clear_X(Q, T, st_MX_cloud, sp_MX_cloud, 0);
   cloud_Fill(Q, T, st_MX_cloud, sp_MX_cloud, edg_row_lin, 1, MODE_ROW);

   /* bounded forward */
   printf("=== BOUNDED FORWARD -> START ===\n");

   bound_Forward_Naive(query_seq, target_prof, Q, T, st_MX, sp_MX, st_MX_cloud, &sc);
   printf("Bounded Forward Score (Naive): %f\n", sc);
   scores->cloud_fwd_sc = sc;
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_fwd.naive.mx");
   
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   bound_Forward_Quad(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row_lin, is_testing, &sc);
   printf("Bounded Forward Score (Quadratic): %f\n", sc);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_fwd.quad.mx");

   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   dp_matrix_Clear3(Q, T, st_MX3, sp_MX);
   bound_Forward_Linear(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row_lin, is_testing, &sc);
   printf("Bound Forward Score (Linear): %f\n", sc);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_fwd.lin.mx");

   printf("=== BOUNDED FORWARD -> END ===\n\n");

   /* bounded backward */
   printf("=== BOUNDED BACKWARD -> START ===\n");

   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   bound_Backward_Naive(query_seq, target_prof, Q, T, st_MX, sp_MX, st_MX_cloud, &sc);
   printf("Bounded Backward Score (Naive): %f\n", sc);
   scores->cloud_bck_sc = sc;
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_bck.naive.mx");

   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   dp_matrix_Clear(Q, T, st_MX_cloud, sp_MX);
   bound_Backward_Quad(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row_lin, is_testing, &sc);
   printf("Bounded Backward Score (Quadratic): %f\n", sc);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/my.bounded_bck.quad.mx");
   dp_matrix_trace_Save(Q, T, st_MX_cloud, sp_MX, tr, "output/my.bounded_bck.test.mx");

   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   dp_matrix_Clear_X3(Q, T, st_MX3, sp_MX, 0);
   bound_Backward_Linear(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row_lin, is_testing, &sc);
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
   fprintf(fp, "%f\t", scores->perc_window);
   fprintf(fp, "%d\t", Q);
   fprintf(fp, "%d\n", T);
   fclose(fp);

   /* FREE MEMORY */
   /* TODO: fix free memory! */

   /* free matrices */
   MATRIX_3D_Destroy(st_MATRIX);
   MATRIX_2D_Destroy(sp_MATRIX);
   MATRIX_3D_Destroy(st_MATRIX3);

   MATRIX_3D_Destroy(st_MATRIX_tmp);

   MATRIX_3D_Destroy(st_MATRIX_cloud);
   MATRIX_2D_Destroy(sp_MATRIX_cloud);

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

   // printf("...test finished. \n");
}

/* standard pipeline */
void cloud_search_pipeline(ARGS *args, char *hmm_file, char *fasta_file, float alpha, int beta) 
{
   /* Target & Query */
   HMM_PROFILE* target_prof = NULL;
   SEQUENCE*    query_seq   = NULL;

   /* PRINT ARGS */
   printf("HMM_FILENAME: %s\n", hmm_file);
   printf("FASTA_FILENAME: %s\n", fasta_file);
   printf("ALPHA: %f\n", alpha);
   printf("BETA: %d\n\n", beta);

   /* Timing & Scoring */
   SCORES *scores      = (SCORES*)malloc( sizeof(SCORES) );
   bool   is_testing   = false;
   bool   test         = false;
   CLOCK* cl           = clock_Create();
   TIMES* times        = (TIMES *) malloc( sizeof(TIMES) );
   time_t t            = 0; 
   time_t tot_t        = 0;

   /* records percentage of space computed */
   int    T            = 0;
   int    Q            = 0;
   float  sc           = 0.f; 
   float  perc_cells   = 0.f;
   int    num_cells    = 0;  
   int    window_cells = 0; 
   int    tot_cells    = (Q + 1) * (T + 1);

   /* MODE => determines whether
   /* mode choices: MODE_NONE, MODE_MULTILOCAL, MODE_MULTIGLOCAL, MODE_UNILOCAL, MODE_UNIGLOCAL */
   char* modes[] = { "None", "Multi-local", "Multi-glocal", "Uni-local", "Uni-glocal" }; 
   int   mode; 
         // mode = MODE_MULTILOCAL;    /* HMMER standard mode (allows jumps) */
         mode = MODE_UNILOCAL;      /* Cloud mode (prohibiits jumps) */
   printf("MODE: %s\n", modes[mode]);

   // printf("=== LOAD HMM_PROFILE / QUERY -> START ===\n");
   clock_Start(cl);

   target_prof = HMM_PROFILE_Parse(hmm_file);
   HMM_PROFILE_Config(target_prof, mode);
   T = target_prof->N;
   // HMM_PROFILE_Dump(target_prof, stdout);

   query_seq = SEQUENCE_Parse(fasta_file);
   Q = query_seq->N;
   // SEQUENCE_Dump(query_seq, stdout);

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->load_hmm = (float)t;
   // printf("=== LOAD HMM_PROFILE / QUERY -> END ===\n\n");

   printf(" QUERY LENGTH: %d\n", Q);
   printf("TARGET LENGTH: %d\n\n", T);

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
   clock_Start(cl);

   viterbi_Quad(query_seq, target_prof, Q, T, st_MX, sp_MX, &sc);
   printf("Viterbi Score: %f\n", sc);
   scores->viterbi_sc = sc;

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->viterbi = (float)t;
   // printf("=== VITERBI -> END ===\n\n");

   /* run traceback of viterbi */
   // printf("=== ALIGNMENT -> START ===\n");
   clock_Start(cl);

   traceback_Build(query_seq, target_prof, Q, T, st_MX, sp_MX, tr);

   TRACE beg = tr->traces[tr->beg];
   TRACE end = tr->traces[tr->end];
   printf("START: (%d,%d) -> END: (%d,%d)\n", beg.i, beg.j, end.i, end.j);
   window_cells = (end.i - beg.i) * (end.j - beg.j);

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->traceback = (float)t;
   // printf("=== ALIGNMENT -> END ===\n\n");

   /* run forward/backward algorithms */
   init_Logsum();

   // printf("=== FORWARD/BACKWARD -> START ===\n");
   clock_Start(cl);

   forward_Quad(query_seq, target_prof, Q, T, st_MX, sp_MX, &sc);
   // printf("Forward Score: %f\n", sc);
   scores->fwd_sc = sc;

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->fwd = (float)t;
   clock_Start(cl);

   backward_Quad(query_seq, target_prof, Q, T, st_MX, sp_MX, &sc);
   // printf("Backward Score: %f\n", sc);
   scores->bck_sc = sc;

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->bck = (float)t;
   // printf("=== FORWARD/BACKWARD -> END ===\n\n");

   /* cloud forward (linear) */
   // printf("=== CLOUD FORWARD/BACKWARD (Linear) -> START ===\n");
   clock_Start(cl);

   cloud_Forward_Linear(query_seq, target_prof, Q, T, st_MX, st_MX3, sp_MX, tr, edg_fwd_lin, alpha, beta, is_testing);
   // EDGEBOUNDS_Save(edg_fwd_lin, "output/my.cloud_fwd.lin.diags.pipe.edg");

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->cloud_fwd = (float)t;
   clock_Start(cl);

   cloud_Backward_Linear(query_seq, target_prof, Q, T, st_MX, st_MX3, sp_MX, tr, edg_bck_lin, alpha, beta, is_testing);
   // EDGEBOUNDS_Save(edg_bck_lin, "output/my.cloud_bck.lin.diags.pipe.edg");

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->cloud_bck = (float)t;
   // printf("=== CLOUD FORWARD/BACKWARD (Linear) -> END ===\n\n");

   /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
   // printf("=== MERGE & REORIENT CLOUD (Linear) -> START ===\n");
   clock_Start(cl);

   edg_diag_lin = EDGEBOUNDS_Merge(Q, T, edg_fwd_lin, edg_bck_lin);

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->merge = (float)t;
   clock_Start(cl);

   edg_row_lin = EDGEBOUNDS_Reorient(Q, T, edg_diag_lin);

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->reorient = (float)t;
   // printf("=== MERGE & REORIENT CLOUD (Linear) -> END ===\n\n");

   // EDGEBOUNDS_Save(edg_row_lin, "output/my.pipeline.row.edg");
   // EDGEBOUNDS_Print(edg_row_lin);

   /* stats */
   num_cells = EDGEBOUNDS_Count(edg_row_lin);
   scores->perc_cells = (float)num_cells/(float)tot_cells;
   printf("PERCENT TOTAL CELLS COMPUTED: %d/%d = %f\n\n", num_cells, tot_cells, scores->perc_cells);
   // scores->perc_window = (float)num_cells/(float)window_cells;
   // printf("Perc. Window Cells Computed = %d/%d = %f\n\n", num_cells, window_cells, scores->perc_window);

   // printf("=== BOUNDED FORWARD/BACKWARD (Linear) -> START ===\n");
   clock_Start(cl);

   bound_Forward_Linear(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row_lin, is_testing, &sc);
   // printf("Bound Forward Score (Linear): %f\n", sc);
   scores->cloud_fwd_sc = sc;

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->bound_fwd = (float)t;
   clock_Start(cl);

   bound_Backward_Linear(query_seq, target_prof, Q, T, st_MX3, st_MX, sp_MX, edg_row_lin, is_testing, &sc);
   // printf("Bound Backward Score (Linear): %f\n", sc);
   scores->cloud_bck_sc = sc;

   clock_Stop(cl);
   t = clock_Ticks(cl);
   times->bound_bck = (float)t;
   // printf("=== BOUNDED FORWARD/BACKWARD (Linear) -> END ===\n\n");

   /* report times */
   tot_t = 0;
   tot_t += times->load_hmm;
   tot_t += times->fwd + times->bck;
   tot_t += times->viterbi + times->traceback;
   tot_t += times->merge + times->reorient;
   tot_t += times->cloud_fwd + times->cloud_bck;
   tot_t += times->bound_fwd + times->bound_bck;

   printf("=== TIMES ===\n");
   printf("FWD: %f\n", times->fwd/tot_t );
   printf("BCK: %f\n", times->bck/tot_t );
   printf("VIT: %f\n", times->viterbi/tot_t );
   printf("TRB: %f\n", times->traceback/tot_t );
   printf("CLOUD_FWD: %f\n", times->cloud_fwd/tot_t );
   printf("CLOUD_BCK: %f\n", times->cloud_bck/tot_t );
   printf("MERGE: %f\n", times->merge/tot_t );
   printf("REORIENT: %f\n", times->reorient/tot_t );
   printf("BOUND_FWD: %f\n", times->bound_fwd/tot_t );
   printf("BOUND_BCK: %f\n", times->bound_bck/tot_t );
   printf("\n");
   float full_fwdbck = times->fwd + times->bck;
   float full_cloud  = times->cloud_fwd + times->cloud_bck + 
                       times->merge + times->reorient + 
                       times->bound_fwd + times->bound_bck;
   printf("FULL FWDBCK: %f\n", full_fwdbck/tot_t );
   printf("FULL CLOUD : %f\n", full_cloud/tot_t );
   printf("\n");
   printf("SPEEDUP: %f\n", full_cloud/full_fwdbck);
   printf("=============\n\n");

   printf("=== SCORES ===\n");
   printf("VITERBI: %f\n", scores->viterbi_sc );
   printf("FWD: %f\n", scores->fwd_sc );
   printf("BCK: %f\n", scores->bck_sc );
   printf("BOUND_FWD: %f\n", scores->cloud_fwd_sc );
   printf("BOUND_BCK: %f\n", scores->cloud_bck_sc );
   printf("==============\n\n");

   /* report stats */
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
   fprintf(fp, "%f\t", scores->perc_window);
   fprintf(fp, "%d\t", Q);
   fprintf(fp, "%d\n", T);
   fclose(fp);

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