/*******************************************************************************
 *  FILE:      pipeline_test.c
 *  PURPOSE:   Test Cloud Search Pipeline.
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

/* pipeline run optimized and unoptimized versions of search algs */
void test_pipeline( WORKER* worker )
{
   /* Commandline Arguments */
   FILE*          fp             = NULL;

   ARGS*          args           = worker->args;

   float          alpha          = args->alpha;
   int            beta           = args->beta;

   char*          t_filepath     = args->t_filepath;
   char*          q_filepath     = args->q_filepath;

   int            t_filetype     = args->t_filetype;
   int            q_filetype     = args->q_filetype;

   char*          output_file    = args->output_filepath;

   /* TODO: multi mode is unneccesary (only support for UNILOCAL) */
   /* mode types: MODE_NONE, MODE_MULTILOCAL, MODE_MULTIGLOCAL, MODE_UNILOCAL, MODE_UNIGLOCAL */

   /* Cloud Search mode (prohibits jumps) */
   int            mode           = MODE_UNILOCAL;      

   /* PROFILE & SEQUENCE */
   HMM_PROFILE*   t_prof         = HMM_PROFILE_Create();
   SEQUENCE*      t_seq          = SEQUENCE_Create();  /* only used if target is a fasta file */
   SEQUENCE*      q_seq          = SEQUENCE_Create();

   /* EDGEBOUNDS & ALIGNMENTS */
   ALIGNMENT*     tr             = ALIGNMENT_Create();

   EDGEBOUNDS*    edg            = EDGEBOUNDS_Create();

   EDGEBOUNDS*    edg_fwd_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_bck_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_row_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_diag_lin   = EDGEBOUNDS_Create();

   EDGEBOUNDS*    edg_fwd_quad   = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_bck_quad   = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_row_quad   = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_diag_quad  = EDGEBOUNDS_Create();

   /* SCORES => stores result scores */
   TIMES*         times          = worker->times;
   SCORES*        scores         = worker->scores;

   /* TEST => embed into quadratic matrix? */
   bool           is_testing     = true;
   int            test           = 0;

   /* records percentage of space computed */
   int            T              = 0;
   int            Q              = 0;
   float          sc             = 0.f; 
   float          perc_cells     = 0.f;
   int            num_cells      = 0;  
   int            window_cells   = 0; 
   int            tot_cells      = 0;

   /* ====================================================================================== */

   /* PRINT ARGS */
   int pad = 20;
   printf("%*s: %s\n", pad, "MODE",             MODE_NAMES[mode]);
   printf("%*s: %s\n", pad, "HMM_FILENAME",     t_filepath);
   printf("%*s: %s\n", pad, "FASTA_FILENAME",   q_filepath);
   printf("%*s: %f\n", pad, "ALPHA",            alpha);
   printf("%*s: %d\n", pad, "BETA",             beta);
   printf("\n");

   printf("=== BUILD HMM_PROFILE / QUERY -> START ===\n");

   /* build q_seq sequence */
   printf("building q_seq sequence...\n");
   if ( q_filetype == FILE_FASTA ) 
   {
      SEQUENCE_Fasta_Parse( q_seq, q_filepath, 0 );
   }
   else 
   {
      fprintf(stderr, "ERROR: Only FASTA filetypes are supported for queries.\n");
      exit(EXIT_FAILURE);
   }
   SEQUENCE_Dump( q_seq, stdout );

   /* build t_prof profile */
   printf("building hmm profile...\n");
   if ( t_filetype == FILE_HMM || true ) 
   {
      HMM_PROFILE_Parse( t_prof, t_filepath, 0 );
      HMM_PROFILE_Convert_NegLog_To_Real( t_prof );
      HMM_PROFILE_Config( t_prof, mode );
   }
   else if ( t_filetype == FILE_FASTA )
   {
      SEQUENCE_Fasta_Parse( t_seq, t_filepath, 0 );
      SEQUENCE_to_HMM_PROFILE( t_seq, t_prof );
   }
   else
   {
      fprintf(stderr, "ERROR: Only HMM and FASTA filetypes are supported for t_profs.\n");
      exit(EXIT_FAILURE);
   }
   // HMM_PROFILE_Dump( t_prof, stdout );

   printf("=== BUILD HMM_PROFILE / QUERY -> END ===\n\n");

   Q = q_seq->N;
   T = t_prof->N;
   tot_cells = (T+1) * (Q+1);

   /* allocate memory for quadratic algs (for DEBUGGING) */
   MATRIX_3D* st_MATRIX       = MATRIX_3D_Create( NUM_NORMAL_STATES,  Q+1, T+1 );
   /* allocate memory for comparing quadratic algs (for DEBUGGING) */
   MATRIX_3D* st_MATRIX_tmp   = MATRIX_3D_Create( NUM_NORMAL_STATES, Q+1, T+1 );
   /* allocate memory for naive cloud matrices (for DEBUGGING) */
   MATRIX_3D* st_MATRIX_cloud = MATRIX_3D_Create( NUM_NORMAL_STATES, Q+1, T+1 );
   MATRIX_2D* sp_MATRIX_cloud = MATRIX_2D_Create( NUM_SPECIAL_STATES, Q+1 );
   /* allocate memory for linear algs */
   MATRIX_3D* st_MATRIX3      = MATRIX_3D_Create( NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
   MATRIX_2D* sp_MATRIX       = MATRIX_2D_Create( NUM_SPECIAL_STATES, Q+1 );

   /* run viterbi algorithm */
   printf("=== VITERBI -> START ===\n");
   viterbi_Quad(q_seq, t_prof, Q, T, st_MATRIX, sp_MATRIX, &sc);
   printf("Viterbi Score: %f\n", sc);
   scores->quad_vit = sc;
   dp_matrix_Save(Q, T, st_MATRIX, sp_MATRIX, "output/my.viterbi.mx");
   printf("=== VITERBI -> END ===\n\n");

   /* run traceback of viterbi */
   printf("=== TRACEBACK -> START ===\n");
   traceback_Build(q_seq, t_prof, Q, T, st_MATRIX, sp_MATRIX, tr);
   ALIGNMENT_Save(tr, "output/my.traceback.tsv");
   TRACE beg = tr->traces[tr->beg];
   TRACE end = tr->traces[tr->end];
   printf("START: (%d,%d) -> END: (%d,%d)\n", beg.i, beg.j, end.i, end.j);
   window_cells = (end.i - beg.i) * (end.j - beg.j);
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.traceback.mx");
   printf("=== TRACEBACK -> END ===\n\n");

   /* run forward/backward algorithms */
   init_Logsum();

   printf("=== FORWARD -> START ===\n");
   /* ==> forward (quadratic) */
   forward_Quad(q_seq, t_prof, Q, T, st_MATRIX, sp_MATRIX, &sc);
   printf("Forward Score (quadratic): %f\n", sc);
   scores->quad_fwd = sc;
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.fwd.quad.mx");
   /* ==> forward (linear) */
   forward_Linear(q_seq, t_prof, Q, T, st_MATRIX3, st_MATRIX, sp_MATRIX, &sc);
   printf("Forward Score (linear): %f\n", sc);
   scores->lin_fwd = sc;
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.fwd.lin.mx");
   printf("=== FORWARD -> END ===\n\n");

   printf("=== BACKWARD -> START ===\n");
   /* ==> backward (quadratic) */
   backward_Quad(q_seq, t_prof, Q, T, st_MATRIX, sp_MATRIX, &sc);
   printf("Backward Score (quadratic): %f\n", sc);
   scores->quad_bck = sc;
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.bck.quad.mx");
   /* ==> backward (linear) */
   backward_Linear(q_seq, t_prof, Q, T, st_MATRIX3, st_MATRIX, sp_MATRIX, &sc);
   printf("Backward Score (linear): %f\n", sc);
   scores->lin_bck = sc;
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.bck.lin.mx");
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
   cloud_Forward_Quad(q_seq, t_prof, Q, T, st_MATRIX_tmp, sp_MATRIX, tr, edg_fwd_quad, alpha, beta);
   dp_matrix_trace_Save(Q, T, st_MATRIX_tmp, sp_MATRIX, tr, "output/my.cloud_fwd_vals.quad.mx");
   dp_matrix_Clear_X(Q, T, st_MATRIX_tmp, sp_MATRIX, 0);
   cloud_Fill(Q, T, st_MATRIX_tmp, sp_MATRIX, edg_fwd_quad, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MATRIX_tmp, sp_MATRIX, tr, "output/my.cloud_fwd.quad.mx");
   EDGEBOUNDS_Save(edg_fwd_quad, "output/my.cloud_fwd.quad.diags.edg");
   printf("=== CLOUD FORWARD (Quadratic) -> END ===\n\n");

   /* cloud forward (linear) */
   printf("=== CLOUD FORWARD (Linear) -> START ===\n");
   cloud_Forward_Linear(q_seq, t_prof, Q, T, st_MATRIX3, st_MATRIX, sp_MATRIX, tr, edg_fwd_lin, alpha, beta, is_testing);
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.cloud_fwd_vals.lin.mx");
   dp_matrix_Clear_X(Q, T, st_MATRIX, sp_MATRIX, 0);
   cloud_Fill(Q, T, st_MATRIX, sp_MATRIX, edg_fwd_lin, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.cloud_fwd.lin.mx");
   EDGEBOUNDS_Save(edg_fwd_lin, "output/my.cloud_fwd.lin.diags.edg");
   test = dp_matrix_Compare(Q, T, st_MATRIX, sp_MATRIX, st_MATRIX_tmp, sp_MATRIX);
   printf("Cloud Forward - Lin vs Quad? %s\n", test ? "TRUE" : "FALSE" );
   printf("=== CLOUD FORWARD (Linear) -> END ===\n\n");

   /* cloud backward (quadratic) */
   printf("=== CLOUD BACKWARD (Quadratic) -> START ===\n");
   cloud_Backward_Quad(q_seq, t_prof, Q, T, st_MATRIX_tmp, sp_MATRIX, tr, edg_bck_quad, alpha, beta);
   dp_matrix_trace_Save(Q, T, st_MATRIX_tmp, sp_MATRIX, tr, "output/my.cloud_bck_vals.quad.mx");
   dp_matrix_Clear_X(Q, T, st_MATRIX_tmp, sp_MATRIX, 0);
   cloud_Fill(Q, T, st_MATRIX_tmp, sp_MATRIX, edg_bck_quad, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MATRIX_tmp, sp_MATRIX, tr, "output/my.cloud_bck.quad.mx");
   EDGEBOUNDS_Save(edg_bck_quad, "output/my.cloud_bck.quad.diags.edg");
   printf("=== CLOUD BACKWARD (Quadratic) -> END ===\n\n");

   /* cloud backward (linear) */
   printf("=== CLOUD BACKWARD (Linear) -> START ===\n");
   cloud_Backward_Linear(q_seq, t_prof, Q, T, st_MATRIX3, st_MATRIX, sp_MATRIX, tr, edg_bck_lin, alpha, beta, is_testing);
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.cloud_bck_vals.lin.mx");
   dp_matrix_Clear_X(Q, T, st_MATRIX, sp_MATRIX, 0);
   cloud_Fill(Q, T, st_MATRIX, sp_MATRIX, edg_bck_lin, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.cloud_bck.lin.mx");
   EDGEBOUNDS_Save(edg_bck_lin, "output/my.cloud_bck.lin.diags.edg");
   test = dp_matrix_Compare(Q, T, st_MATRIX, sp_MATRIX, st_MATRIX_tmp, sp_MATRIX);
   printf("Cloud Backward - Lin vs Quad? %s\n", test ? "TRUE" : "FALSE" );
   printf("=== CLOUD BACKWARD (Linear) -> END ===\n\n");

   /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
   printf("=== MERGE & REORIENT CLOUD (Naive) -> START ===\n");
   EDGEBOUNDS_Merge_Reorient_Naive(edg_fwd_quad, edg_bck_quad, edg_diag_quad, edg_row_quad, Q, T, st_MATRIX, sp_MATRIX);

   dp_matrix_Clear_X(Q, T, st_MATRIX, sp_MATRIX, 0);
   cloud_Fill(Q, T, st_MATRIX, sp_MATRIX, edg_diag_quad, 1, MODE_DIAG);
   EDGEBOUNDS_Save(edg_diag_quad, "output/my.cloud.quad.diags.edg");
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.cloud.quad.diags.mx");
   dp_matrix_Clear_X(Q, T, st_MATRIX, sp_MATRIX, 0);
   cloud_Fill(Q, T, st_MATRIX, sp_MATRIX, edg_row_quad, 1, MODE_ROW);
   EDGEBOUNDS_Save(edg_row_quad, "output/my.cloud.quad.rows.edg");
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.cloud.quad.rows.mx");
   dp_matrix_Copy(Q, T, st_MATRIX, sp_MATRIX, st_MATRIX_cloud, sp_MATRIX_cloud);
   printf("=== MERGE & REORIENT CLOUD (Naive) -> END ===\n\n");

   /* stats */
   num_cells = cloud_Cell_Count(Q, T, st_MATRIX, sp_MATRIX);
   scores->perc_cells = (float)num_cells/(float)tot_cells;
   printf("Perc. Total Cells Computed = %d/%d = %f\n", num_cells, tot_cells, scores->perc_cells);
   scores->perc_window = (float)num_cells/(float)window_cells;
   printf("Perc. Window Cells Computed = %d/%d = %f\n", num_cells, window_cells, scores->perc_window);
   printf("\n");

   /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
   printf("=== MERGE & REORIENT CLOUD (Linear) -> START ===\n");
   EDGEBOUNDS_Merge(Q, T, edg_fwd_lin, edg_bck_lin, edg_diag_lin);

   dp_matrix_Clear_X(Q, T, st_MATRIX, sp_MATRIX, 0);
   cloud_Fill(Q, T, st_MATRIX, sp_MATRIX, edg_diag_lin, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.cloud.lin.diags.mx");
   EDGEBOUNDS_Save(edg_diag_lin, "output/my.cloud.lin.diags.edg");
   EDGEBOUNDS_Reorient(Q, T, edg_diag_lin, edg_row_lin);
   dp_matrix_Clear_X(Q, T, st_MATRIX, sp_MATRIX, 0);
   cloud_Fill(Q, T, st_MATRIX, sp_MATRIX, edg_row_lin, 1, MODE_ROW);
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.cloud.lin.rows.mx");
   EDGEBOUNDS_Save(edg_row_lin, "output/my.cloud.lin.rows.edg");
   printf("=== MERGE & REORIENT CLOUD (Linear) -> END ===\n\n");

   /* create cloud that covers entire matrix (full fwd/bck) */
   // printf("=== TEST CLOUD -> START ===\n");
   // dp_matrix_Clear_X(Q, T, st_MATRIX_cloud, sp_MATRIX_cloud, 1);
   // EDGEBOUNDS_Build_From_Cloud(edg, Q, T, st_MATRIX_cloud, MODE_DIAG);
   // EDGEBOUNDS_Save(edg, "output/my.full_fwdbck.edg");
   // EDGEBOUNDS_Dump(edg, stdout);
   // dp_matrix_trace_Save(Q, T, st_MATRIX_cloud, sp_MATRIX_cloud, tr, "output/my.cloud.full_fwdbck.mx");
   // printf("=== TEST CLOUD -> END ===\n\n");

   /* fill cloud for naive algs */
   dp_matrix_Clear(Q, T, st_MATRIX, sp_MATRIX);
   dp_matrix_Clear_X(Q, T, st_MATRIX_cloud, sp_MATRIX_cloud, 0);
   cloud_Fill(Q, T, st_MATRIX_cloud, sp_MATRIX_cloud, edg_row_lin, 1, MODE_ROW);

   /* bounded forward */
   printf("=== BOUNDED FORWARD -> START ===\n");

   bound_Forward_Naive(q_seq, t_prof, Q, T, st_MATRIX, sp_MATRIX, st_MATRIX_cloud, &sc);
   printf("Bounded Forward Score (Naive): %f\n", sc);
   scores->naive_cloud_fwd = sc;
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.bounded_fwd.naive.mx");
   
   bound_Forward_Quad(q_seq, t_prof, Q, T, st_MATRIX, sp_MATRIX, edg_row_lin, is_testing, &sc);
   printf("Bounded Forward Score (Quadratic): %f\n", sc);
   scores->quad_cloud_fwd = sc;
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.bounded_fwd.quad.mx");

   bound_Forward_Linear(q_seq, t_prof, Q, T, st_MATRIX3, st_MATRIX, sp_MATRIX, edg_row_lin, is_testing, &sc);
   printf("Bound Forward Score (Linear): %f\n", sc);
   scores->lin_cloud_fwd = sc;
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.bounded_fwd.lin.mx");

   printf("=== BOUNDED FORWARD -> END ===\n\n");

   /* bounded backward */
   printf("=== BOUNDED BACKWARD -> START ===\n");

   bound_Backward_Naive(q_seq, t_prof, Q, T, st_MATRIX, sp_MATRIX, st_MATRIX_cloud, &sc);
   printf("Bounded Backward Score (Naive): %f\n", sc);
   scores->naive_cloud_bck = sc;
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.bounded_bck.naive.mx");

   bound_Backward_Quad(q_seq, t_prof, Q, T, st_MATRIX, sp_MATRIX, edg_row_lin, is_testing, &sc);
   printf("Bounded Backward Score (Quadratic): %f\n", sc);
   scores->quad_cloud_fwd = sc;
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.bounded_bck.quad.mx");
   dp_matrix_trace_Save(Q, T, st_MATRIX_cloud, sp_MATRIX, tr, "output/my.bounded_bck.test.mx");

   bound_Backward_Linear(q_seq, t_prof, Q, T, st_MATRIX3, st_MATRIX, sp_MATRIX, edg_row_lin, is_testing, &sc);
   printf("Bound Backward Score (Linear): %f\n", sc);
   scores->lin_cloud_fwd = sc;
   dp_matrix_trace_Save(Q, T, st_MATRIX, sp_MATRIX, tr, "output/my.bounded_bck.lin.mx");

   printf("=== BOUNDED BACKWARD -> END ===\n\n");

   /* sample stats */
   printf("Writing results to: '%s'\n", output_file);
   fp = fopen(output_file, "a+");

   fprintf(fp, "%s\t", t_filepath);
   fprintf(fp, "%s\t", q_filepath);
   fprintf(fp, "%f\t", scores->quad_vit);
   fprintf(fp, "%f\t", scores->quad_fwd);
   fprintf(fp, "%f\t", scores->quad_bck);
   fprintf(fp, "%f\t", scores->naive_cloud_fwd);
   fprintf(fp, "%f\t", scores->naive_cloud_bck);
   fprintf(fp, "%f\t", alpha);
   fprintf(fp, "%d\t", beta);
   fprintf(fp, "%f\t", scores->perc_cells);
   fprintf(fp, "%f\t", scores->perc_window);
   fprintf(fp, "%d\t", Q);
   fprintf(fp, "%d\n", T);

   if ( fp != stdout ) fclose( fp );
}