/*******************************************************************************
 *  FILE:      pipeline_int_test.c
 *  PURPOSE:   Test Cloud Search Pipeline.
 *             Requires -DDEBUG=1 for full functionality.
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

/*
 *  FUNCTION:  itest_pipeline()
 *  SYNOPSIS:  Pipeline runs integration tests. 
 *             Runs optimized and unoptimized versions of search algs and compares results.
 *             For full functionality, must be compiled in DEBUG mode.
 */
void itest_pipeline( WORKER* worker )
{
   #if DEBUG
   {
      fprintf(stdout, "Running integration test...\n");
      debugger->verbose_level = VERBOSE_LOW;
   }
   #else
   {
      fprintf(stdout, "Test pipeline only supported in DEBUG build, not PRODUCTION build.  Recompile with BUILD=DEBUG.\n" );
      exit(EXIT_FAILURE);
   }
   #endif 

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
   /* mode types: EDG_NONE, MODE_MULTILOCAL, MODE_MULTIGLOCAL, MODE_UNILOCAL, MODE_UNIGLOCAL */

   /* Cloud Search mode (prohibits jumps) */
   int            mode           = MODE_UNILOCAL;      

   /* PROFILE & SEQUENCE */
   HMM_PROFILE*   t_prof         = HMM_PROFILE_Create();
   SEQUENCE*      t_seq          = SEQUENCE_Create();  /* only used if target is a fasta file */
   SEQUENCE*      q_seq          = SEQUENCE_Create();

   /* EDGEBOUNDS & ALIGNMENTS */
   ALIGNMENT*     tr             = ALIGNMENT_Create();

   EDGEBOUNDS*    edg_fwd_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_bck_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_row_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_diag_lin   = EDGEBOUNDS_Create();

   EDGEBOUNDS*    edg_fwd_quad   = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_bck_quad   = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_row_quad   = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_diag_quad  = EDGEBOUNDS_Create();

   /* temporary edgebound object for storing row-wise edgebounds during cloud search */
   EDGEBOUND_ROWS* edg_row_tmp   = EDGEBOUND_ROWS_Create();

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
   printf("loading query...\n");
   if ( q_filetype == FILE_FASTA ) 
   {
      SEQUENCE_Fasta_Parse( q_seq, q_filepath, 0 );
   }
   else 
   {
      fprintf(stderr, "ERROR: Only FASTA filetypes are supported for queries.\n");
      exit(EXIT_FAILURE);
   }
   // SEQUENCE_Dump( q_seq, stdout );

   /* build t_prof profile */
   printf("loading target...\n");
   if ( t_filetype == FILE_HMM || true ) 
   {
      HMM_PROFILE_Parse( t_prof, t_filepath, 0 );
      HMM_PROFILE_Convert_NegLog_To_Real( t_prof );
      HMM_PROFILE_Config( t_prof, mode );
      HMM_PROFILE_ReconfigLength( t_prof, q_seq->N );
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

   printf("TARGET LEN:\t%d\n", t_prof->N);
   printf(" QUERY LEN:\t%d\n", q_seq->N);
   printf("=== BUILD HMM_PROFILE / QUERY -> END ===\n\n");

   Q = q_seq->N;
   T = t_prof->N;
   tot_cells = (T+1) * (Q+1);

   /* resize edgebounds */
   EDGEBOUNDS_Reuse( edg_fwd_lin, Q, T );
   EDGEBOUNDS_Reuse( edg_bck_lin, Q, T );
   EDGEBOUNDS_Reuse( edg_row_lin, Q, T );
   EDGEBOUNDS_Reuse( edg_diag_lin, Q, T );

   EDGEBOUNDS_Reuse( edg_fwd_quad, Q, T );
   EDGEBOUNDS_Reuse( edg_bck_quad, Q, T );
   EDGEBOUNDS_Reuse( edg_row_quad, Q, T );
   EDGEBOUNDS_Reuse( edg_diag_quad, Q, T );

   EDGEBOUND_ROWS_Reuse( edg_row_tmp, Q, T );

   /* allocate memory for quadratic algs */
   MATRIX_3D*  st_MX_naive    = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES,  Q+1, T+1 );
   MATRIX_2D*  sp_MX_naive    = MATRIX_2D_Create_Clean( NUM_SPECIAL_STATES, Q+1 );
   MATRIX_2D*  cloud_MX_naive = MATRIX_2D_Create_Clean( Q+1, T+1 );
   /* allocate memory for quadratic algs */
   MATRIX_3D*  st_MX_quad     = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES,  Q+1, T+1 );
   MATRIX_2D*  sp_MX_quad     = MATRIX_2D_Create_Clean( NUM_SPECIAL_STATES, Q+1 );
   MATRIX_2D*  cloud_MX_quad  = MATRIX_2D_Create_Clean( Q+1, T+1 );
   /* allocate memory for comparing row-wise algs */
   MATRIX_2D*  cloud_MX_rows  = MATRIX_2D_Create_Clean( Q+1, T+1 );
   /* allocate memory for linear algs */
   MATRIX_3D*  st_MX_lin      = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES,  Q+1, T+1 );
   MATRIX_3D*  st_MX3_lin     = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
   MATRIX_2D*  sp_MX_lin      = MATRIX_2D_Create_Clean( NUM_SPECIAL_STATES, Q+1 );
   MATRIX_2D*  cloud_MX_lin   = MATRIX_2D_Create_Clean( Q+1, T+1 );
   /* allocate memory for testing */
   MATRIX_3D*  st_MX_diff     = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES,  Q+1, T+1 );
   MATRIX_3D*  st_MX3_diff    = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
   MATRIX_2D*  sp_MX_diff     = MATRIX_2D_Create_Clean( NUM_SPECIAL_STATES, Q+1 );
   MATRIX_2D*  cloud_MX_diff  = MATRIX_2D_Create_Clean( Q+1, T+1 );

   /* allocate full search cloud */
   MATRIX_2D*  full_cloud_MX  = MATRIX_2D_Create_Clean( Q+1, T+1 );
   EDGEBOUNDS* full_cloud_edg = EDGEBOUNDS_Create();
   MATRIX_2D_Fill( full_cloud_MX, 1 );
   EDGEBOUNDS_Build_From_Cloud( Q, T, full_cloud_edg, full_cloud_MX, EDG_ROW );

   /* debug matrix */
   #if DEBUG
   {
      MATRIX_3D_Reuse_Clean( debugger->test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_2D_Reuse_Clean( debugger->cloud_MX, Q+1, T+1 );
   }
   #endif

   #if DEBUG
   {
      // printf("=== TEST -> START ===\n");
      // fwd_test_cycle(Q, T, st_MX, sp_MX, tr);
      // bck_test_cycle(Q, T, st_MX, sp_MX, tr);
      // fwd_test_cycle3(Q, T, st_MX, st_MX3, sp_MX, tr);
      // bck_test_cycle3(Q, T, st_MX, st_MX3, sp_MX, tr);
      // printf("=== TEST -> END ===\n\n"); 
   }
   #endif

   /* run viterbi algorithm */
   {
      printf("=== VITERBI -> START ===\n");
      /* ==> viterbi (quadratic) */
      viterbi_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, &sc);
      printf("Viterbi Score (quad):\t%f\n", sc);
      scores->quad_vit = sc;
      DP_MATRIX_Save(Q, T, st_MX_quad, sp_MX_quad, "test_output/my.viterbi.quad.mx");
      /* ==> viterbi (linear) */
      /* TODO */
      printf("=== VITERBI -> END ===\n\n");
   }

   /* run traceback of viterbi */
   {
      printf("=== TRACEBACK -> START ===\n");
      /* ==> viterbi (quadratic) */
      traceback_Build(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr);
      ALIGNMENT_Save(tr, "test_output/my.traceback.tsv");
      TRACE* beg = &(tr->traces[tr->beg]);
      TRACE* end = &(tr->traces[tr->end]);
      printf("START: (%d,%d) -> END: (%d,%d)\n", beg->i, beg->j, end->i, end->j);
      #if DEBUG
      {
         if ( debugger->verbose_level >= VERBOSE_ALL ) {
            DP_MATRIX_VIZ_Dump( cloud_MX_quad, debugout );
         }
      }
      #endif
      /* ==> viterbi (linear) */
      /* TODO */
      printf("=== TRACEBACK -> END ===\n\n");
   }

   /* run forward/backward algorithms */
   logsum_Init();

   /* run forward */
   {
      printf("=== FORWARD -> START ===\n");
      /* ==> forward (quadratic) */
      forward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, &sc);
      printf("Forward Score (quad):\t%f\n", sc);
      scores->quad_fwd = sc;
      DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.fwd.quad.mx");
      /* ==> forward (linear) */
      forward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, &sc);
      printf("Forward Score  (lin):\t%f\n", sc);
      scores->lin_fwd = sc;
      MATRIX_3D_Copy( debugger->test_MX, st_MX_lin );
      // DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.fwd.lin.mx");
      /* ==> forward (comparison) */
      if ( debugger->verbose_level >= VERBOSE_ALL ) {
         printf("forward (quad):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_quad, debugout );
         printf("forward (lin):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_lin, debugout );
      }
      int dp_cmp  = DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin );
      int cmp = ( dp_cmp == 0 ) ? true : false;
      printf("Forward:\tVALUES?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      printf("=== FORWARD -> END ===\n\n");
   }

   /* run backward */
   {
      printf("=== BACKWARD -> START ===\n");
      /* ==> backward (quadratic) */
      backward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, &sc);
      printf("Backward Score (quad):\t%f\n", sc);
      scores->quad_bck = sc;
      DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.bck.quad.mx");
      /* ==> backward (linear) */
      backward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, &sc);
      printf("Backward Score  (lin):\t%f\n", sc);
      scores->lin_bck = sc;
      MATRIX_3D_Copy( debugger->test_MX, st_MX_lin );
      // DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.bck.lin.mx");
      /* ==> backward (comparison) */
      if ( debugger->verbose_level >= VERBOSE_ALL ) {
         printf("backward (quad):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_quad, debugout );
         printf("backward (lin):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_lin, debugout );
      }
      int dp_cmp  = ( DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin ) == 0 );
      printf("Backward:\tVALUES?\t\t%s\n", dp_cmp ? "PASS" : "FAIL" );
      printf("=== BACKWARD -> END ===\n\n");
   }

   /* need to clean matrix after fwd/bck */
   DP_MATRIX_Fill( Q, T, st_MX3_lin, sp_MX_lin, -INF );

   /* run cloud forward */
   {
      /* cloud forward (quadratic) */
      printf("=== CLOUD FORWARD (quadratic) -> START ===\n");
      cloud_Forward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr, edg_fwd_quad, alpha, beta);
      if ( debugger->verbose_level >= VERBOSE_ALL ) {
         MATRIX_2D_Copy( cloud_MX_quad, debugger->cloud_MX );
         DP_MATRIX_VIZ_Dump( cloud_MX_quad, debugout );
      }
      DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.cloud_fwd.quad.mx");
      EDGEBOUNDS_Save(edg_fwd_quad, "test_output/my.cloud_fwd.quad.diags.edg");
      printf("=== CLOUD FORWARD (quadratic) -> END ===\n\n");
      /* cloud forward (linear) */
      printf("=== CLOUD FORWARD (linear) -> START ===\n");
      // cloud_Forward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, tr, edg_fwd_lin, alpha, beta);
      cloud_Forward_Linear_Rows(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, tr, edg_row_tmp, edg_fwd_lin, &(worker->cloud_params) );
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
      if ( debugger->verbose_level >= VERBOSE_ALL ) {
         MATRIX_2D_Copy( cloud_MX_lin, debugger->cloud_MX );
         DP_MATRIX_VIZ_Dump( cloud_MX_lin, debugout );
      }
      #if ( CLOUD_METHOD == CLOUD_DIAGS )
      DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.cloud_fwd.lin.diags.mx");
      EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_fwd.lin.diags.edg");
      #elif ( CLOUD_METHOD == CLOUD_ROWS )
      DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.cloud_fwd.lin.rows.mx");
      EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_fwd.lin.rows.edg");
      #endif
      printf("=== CLOUD FORWARD (linear) -> END ===\n\n");

      /* cloud forward (comparison) */
      printf("=== CLOUD FORWARD (lin v quad - comparison test) -> START ===\n");
      /* compare dp matrix values */      
      int dp_cmp = ( DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin ) == 0 );
      printf("Cloud Forward:\tMATRIX?\t\t%s\n", ( dp_cmp == true ) ? "PASS" : "FAIL" );
      if ( ( dp_cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) ) {
         DP_MATRIX_Diff( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin, st_MX_diff, sp_MX_diff );
         fprintf( debugout, "=> DP MATRIX DIFF:\n" );
         DP_MATRIX_Dump( Q, T, st_MX_diff, sp_MX_diff, debugout );
      }  
      /* compare cloud */
      int cloud_cmp = ( EDGEBOUNDS_Compare_by_Cloud( edg_fwd_quad, cloud_MX_quad, edg_fwd_lin, cloud_MX_lin ) == 0 );
      printf("Cloud Forward:\tCLOUD?\t\t%s\n", ( cloud_cmp == true ) ? "PASS" : "FAIL" );
      if ( cloud_cmp == false ) {
         DP_MATRIX_VIZ_Compare( cloud_MX_diff, edg_fwd_lin, edg_fwd_quad );
         FILE* fp = fopen( "test_output/my.quadvlin.cloud_fwd.viz", "w+" );
         DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, fp );
         fclose(fp);
      }
      /* compare edgebounds */
      int edg_cmp = ( EDGEBOUNDS_Compare( edg_fwd_quad, edg_fwd_lin ) == 0 );
      printf("Cloud Forward:\tEDGES?\t\t%s\n", ( edg_cmp == true ) ? "PASS" : "FAIL" );
      if ( ( edg_cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) ) {
         fprintf( debugout, "=> LINEAR EDGEBOUNDS:\n" );
         EDGEBOUNDS_Dump( edg_fwd_lin, debugout );
         fprintf( debugout, "=> QUADRATIC EDGEBOUNDS:\n" );
         EDGEBOUNDS_Dump( edg_fwd_quad, debugout );
      }
      printf("=== CLOUD FORWARD (comparison test) -> END ===\n\n");
   }

   /* run cloud backward */
   {
      /* cloud backward (quadratic) */
      printf("=== CLOUD BACKWARD (quadratic) -> START ===\n");
      cloud_Backward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr, edg_bck_quad, alpha, beta);
      if ( debugger->verbose_level >=  VERBOSE_ALL ) {
         MATRIX_2D_Copy( cloud_MX_quad, debugger->cloud_MX );
         DP_MATRIX_VIZ_Dump( cloud_MX_quad, debugout );
      }
      DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.cloud_bck.quad.mx");
      EDGEBOUNDS_Save(edg_bck_quad, "test_output/my.cloud_bck.quad.diags.edg");
      printf("=== CLOUD BACKWARD (quadratic) -> END ===\n\n");
      /* cloud backward (linear) */
      printf("=== CLOUD BACKWARD (linear) -> START ===\n");
      // cloud_Backward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, tr, edg_bck_lin, alpha, beta);
      cloud_Backward_Linear_Rows(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, tr, edg_row_tmp, edg_bck_lin, &(worker->cloud_params) );
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
      if ( debugger->verbose_level >=  VERBOSE_ALL ) {
         MATRIX_2D_Copy( cloud_MX_lin, debugger->cloud_MX );
         // DP_MATRIX_VIZ_Dump( cloud_MX_lin, debugout );
      }
      #if ( CLOUD_METHOD == CLOUD_DIAGS )
      DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.cloud_bck.lin.diags.mx");
      EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_bck.lin.diags.edg");
      #elif ( CLOUD_METHOD == CLOUD_ROWS )
      DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.cloud_bck.lin.rows.mx");
      EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_bck.lin.rows.edg");
      #endif
      printf("=== CLOUD BACKWARD (linear) -> END ===\n\n");

      /* cloud backward (comparison) */
      printf("=== CLOUD BACKWARD (lin v quad - comparison test) -> START ===\n");
      /* BUG: cloud_backward_linear() and cloud_backward_quad() have different scores at left edge */
      /* compare dp matrix values */
      int dp_cmp = ( DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin ) == 0 );
      printf("Cloud Backward:\tMATRIX?\t\t%s\n", ( dp_cmp == true ) ? "PASS" : "FAIL" );
      if ( ( dp_cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) ) {
         DP_MATRIX_Diff( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin, st_MX_diff, sp_MX_diff );
         fprintf( debugout, "=> DP MATRIX DIFF:\n" );
         DP_MATRIX_Dump( Q, T, st_MX_diff, sp_MX_diff, debugout );
      } 
      /* compare cloud cells */
      int cloud_cmp = ( EDGEBOUNDS_Compare_by_Cloud( edg_bck_quad, cloud_MX_quad, edg_bck_lin, cloud_MX_lin ) == 0 );
      printf("Cloud Backward:\tCLOUD?\t\t%s\n", ( cloud_cmp == true ) ? "PASS" : "FAIL" );
      if ( ( cloud_cmp == false ) ) {
         DP_MATRIX_VIZ_Compare( cloud_MX_diff, edg_bck_lin, edg_bck_quad );
         FILE* fp = fopen( "test_output/my.quadvlin.cloud_bck.viz", "w+" );
         DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, fp );
         fclose(fp);
      }
      /* compare edgebounds data */
      int edg_cmp = ( EDGEBOUNDS_Compare( edg_bck_quad, edg_bck_lin ) == 0 );
      printf("Cloud Backward:\tEDGES?\t\t%s\n", ( edg_cmp == true ) ? "PASS" : "FAIL" );
      if ( ( edg_cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) ) {
         fprintf( debugout, "=> LINEAR EDGEBOUNDS:\n");
         EDGEBOUNDS_Dump( edg_bck_lin, debugout );
         fprintf( debugout, "=> QUADRATIC EDGEBOUNDS:\n");
         EDGEBOUNDS_Dump( edg_bck_quad, debugout );
      }
      printf("=== CLOUD BACKWARD (comparison test) -> END ===\n\n");
   }

   /* set edgebounds equivalent for merge/reorient testing */
   EDGEBOUNDS_Copy( edg_fwd_lin, edg_fwd_quad );
   EDGEBOUNDS_Copy( edg_bck_lin, edg_bck_quad );
   /* visualize forward and backward clouds */
   {
      DP_MATRIX_VIZ_Compare( cloud_MX_diff, edg_fwd_lin, edg_bck_quad );
      DP_MATRIX_VIZ_Trace( cloud_MX_diff, tr );
      // DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, stdout );
      FILE* fp = fopen( "test_output/my.fwdbck_cloud.viz", "w+" );
      DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, fp );
      fclose(fp);
   }

   /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
   {
      printf("=== MERGE & REORIENT CLOUD (naive) -> START ===\n");
      EDGEBOUNDS_Merge_Reorient_Naive(Q, T, edg_fwd_quad, edg_bck_quad, edg_diag_quad, edg_row_quad, cloud_MX_quad);
      EDGEBOUNDS_Save(edg_diag_quad, "test_output/my.cloud.quad.diags.edg");
      EDGEBOUNDS_Save(edg_row_quad, "test_output/my.cloud.quad.rows.edg");
      int edg_cmp_quad = ( EDGEBOUNDS_Compare_by_Cloud( edg_diag_quad, cloud_MX_quad, edg_row_quad, cloud_MX_diff ) == 0 );
      printf("Rows vs Diags:\tEDGES?\t\t%s\n", ( edg_cmp_quad == true ) ? "PASS" : "FAIL" );
      printf("=== MERGE & REORIENT CLOUD (naive) -> END ===\n\n");

      /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
      printf("=== MERGE & REORIENT CLOUD (linear) -> START ===\n");
      EDGEBOUNDS_Merge(Q, T, edg_fwd_lin, edg_bck_lin, edg_diag_lin);
      EDGEBOUNDS_Reorient_to_Row(Q, T, edg_diag_lin, edg_row_lin);
      EDGEBOUNDS_Save(edg_diag_lin, "test_output/my.cloud.lin.diags.edg");
      EDGEBOUNDS_Save(edg_row_lin, "test_output/my.cloud.lin.rows.edg");
      int edg_cmp_lin = ( EDGEBOUNDS_Compare_by_Cloud( edg_diag_lin, cloud_MX_lin, edg_row_lin, cloud_MX_diff ) == 0 );
      printf("Rows vs Diags:\tEDGES?\t\t%s\n", ( edg_cmp_lin == true ) ? "PASS" : "FAIL" );
      printf("=== MERGE & REORIENT CLOUD (linear) -> END ===\n\n");

      printf("=== MERGE & REORIENT CLOUD ( lin v quad - comparison test) -> START ===\n");
      /* test post-merge */
      int diag_cmp = ( EDGEBOUNDS_Compare_by_Cloud( edg_diag_lin, cloud_MX_lin, edg_row_quad, cloud_MX_quad ) == 0 );
      printf("Linear vs Quad:\t    POST-MERGE?\t\t%s\n", ( diag_cmp == true ) ? "PASS" : "FAIL" );
      if ( ( diag_cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) || true ) {
         DP_MATRIX_VIZ_Compare( cloud_MX_diff, edg_diag_lin, edg_diag_quad );
         DP_MATRIX_VIZ_Trace( cloud_MX_diff, tr );
         // DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, stdout );
         FILE* fp = fopen( "test_output/my.quadvlin.diag.viz", "w+" );
         DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, fp );
         fclose(fp);
      } 
      /* test post-reorient */
      int row_cmp = ( EDGEBOUNDS_Compare_by_Cloud( edg_row_lin, cloud_MX_lin, edg_row_quad, cloud_MX_quad ) == 0 );
      printf("Linear vs Quad:\t POST-REORIENT?\t\t%s\n", ( row_cmp == true ) ? "PASS" : "FAIL" );
      if ( ( diag_cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) || true ) {
         DP_MATRIX_VIZ_Compare( cloud_MX_diff, edg_row_lin, edg_row_quad );
         DP_MATRIX_VIZ_Trace( cloud_MX_diff, tr );
         // DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, stdout );
         FILE* fp = fopen( "test_output/my.quadvlin.diag.viz", "w+" );
         DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, fp );
         fclose(fp);
      } 
      printf("=== MERGE & REORIENT CLOUD (lin v quad - comparison test) -> END ===\n\n");
   }

   // /* stats */
   // num_cells = cloud_Cell_Count(Q, T, st_MX, sp_MX);
   // scores->perc_cells = (float)num_cells/(float)tot_cells;
   // printf("Perc. Total Cells Computed = %d/%d = %f\n", num_cells, tot_cells, scores->perc_cells);
   // scores->perc_window = (float)num_cells/(float)window_cells;
   // printf("Perc. Window Cells Computed = %d/%d = %f\n", num_cells, window_cells, scores->perc_window);
   // printf("\n");

   // /* create cloud that covers entire matrix (full fwd/bck) */
   // // printf("=== TEST CLOUD -> START ===\n");
   // // DP_MATRIX_Clear_X(Q, T, st_MX_cloud, sp_MX_cloud, 1);
   // // EDGEBOUNDS_Build_From_Cloud(edg, Q, T, st_MX_cloud, EDG_DIAG);
   // // EDGEBOUNDS_Save(edg, "test_output/my.full_fwdbck.edg");
   // // EDGEBOUNDS_Dump(edg, stdout);
   // // DP_MATRIX_Trace_Save(Q, T, st_MX_cloud, sp_MX_cloud, tr, "test_output/my.cloud.full_fwdbck.mx");
   // // printf("=== TEST CLOUD -> END ===\n\n");

   /* fill cloud for naive algs */
   MATRIX_2D_Fill( cloud_MX_naive, 0 );
   MATRIX_2D_Cloud_Fill( cloud_MX_naive, edg_row_lin, 1 );
   // DP_MATRIX_VIZ_Dump( cloud_MX_naive, stdout );
   // EDGEBOUNDS_Dump( edg_row_lin, stdout );

   /* bounded forward */
   printf("=== BOUNDED FORWARD -> START ===\n");
   /* bounded forward (naive) */
   bound_Forward_Naive(q_seq, t_prof, Q, T, st_MX_naive, sp_MX_naive, cloud_MX_naive, &sc);
   printf("Bounded Forward Score (naive):\t%f\n", sc);
   scores->naive_cloud_fwd = sc;
   DP_MATRIX_Trace_Save(Q, T, st_MX_naive, sp_MX_naive, tr, "test_output/my.bound_fwd.naive.mx");
   /* bounded forward (quadratic) */
   bound_Forward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, edg_row_lin, &sc);
   printf("Bounded Forward Score  (quad):\t%f\n", sc);
   scores->quad_cloud_fwd = sc;
   DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.bound_fwd.quad.mx");
   // DP_MATRIX_VIZ_Dump( debugger->cloud_MX, stdout );
   /* bounded forward (comparison) */
   bound_Forward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, edg_row_lin, &sc);
   printf("Bounded Forward Score   (lin):\t%f\n", sc);
   scores->lin_cloud_fwd = sc;
   MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
   DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.bound_fwd.lin.mx");
   // DP_MATRIX_VIZ_Dump( debugger->cloud_MX, stdout );
   /* bounded forward (comparison) */
   int sc_cmp = ( cmp_tol_cust( scores->naive_cloud_fwd, scores->quad_cloud_fwd, 1e-3 ) && 
                  cmp_tol_cust( scores->naive_cloud_fwd, scores->lin_cloud_fwd, 1e-3 ) );
   printf("Bounded Forward:\tSCORES?\t\t%s\n", sc_cmp ? "PASS" : "FAIL" );
   int dp_cmp  = ( DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin ) == 0 );
   printf("Bounded Forward:\tVALUES?\t\t%s\n", dp_cmp ? "PASS" : "FAIL" );
   if ( dp_cmp == false ) {
      DP_MATRIX_Diff( st_MX_lin, sp_MX_lin, st_MX_quad, sp_MX_quad, st_MX_diff, sp_MX_diff );
      // DP_MATRIX_Dump( Q, T, st_MX_diff, sp_MX_diff, stdout );
   }
   printf("=== BOUNDED FORWARD -> END ===\n\n");

   /* bounded backward */
   printf("=== BOUNDED BACKWARD -> START ===\n");
   /* bounded backward (naive) */
   bound_Backward_Naive(q_seq, t_prof, Q, T, st_MX_naive, sp_MX_naive, cloud_MX_naive, &sc);
   printf("Bounded Backward Score (naive):\t%f\n", sc);
   scores->naive_cloud_bck = sc;
   DP_MATRIX_Trace_Save(Q, T, st_MX_naive, sp_MX_naive, tr, "test_output/my.bound_bck.naive.mx");
   /* bounded backward (quadratic) */
   bound_Backward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, edg_row_lin, &sc);
   printf("Bounded Backward Score  (quad):\t%f\n", sc);
   scores->quad_cloud_fwd = sc;
   DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.bound_bck.quad.mx");
   /* bounded backward (linear) */
   bound_Backward_Linear(q_seq, t_prof, Q, T, st_MX_lin, sp_MX_lin, edg_row_lin, &sc);
   printf("Bounded Backward Score   (lin):\t%f\n", sc);
   scores->lin_cloud_fwd = sc;
   MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
   DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.bound_bck.lin.mx");
   /* bounded backward (comparison) */
   #if DEBUG
   {
      int sc_cmp = ( cmp_tol_cust( scores->naive_cloud_bck, scores->quad_cloud_bck, 1e-3 ) && 
                     cmp_tol_cust( scores->naive_cloud_bck, scores->lin_cloud_bck, 1e-3 ) );
      printf("Bounded Backward:\tSCORES?\t\t%s\n", sc_cmp ? "PASS" : "FAIL" );
      int dp_cmp  = ( DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin ) == 0 );
      printf("Bounded Backward:\tVALUES?\t\t%s\n", dp_cmp ? "PASS" : "FAIL" );
      if ( ( dp_cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) ) {
         DP_MATRIX_Diff( st_MX_lin, sp_MX_lin, st_MX_quad, sp_MX_quad, st_MX_diff, sp_MX_diff );
         // DP_MATRIX_Dump( Q, T, st_MX_diff, sp_MX_diff, stdout );
      }
   }
   #endif
   printf("=== BOUNDED BACKWARD -> END ===\n\n");

   // /* sample stats */
   // printf("Writing results to: '%s'\n", output_file);
   // fp = fopen(output_file, "a+");

   // fprintf(fp, "%s\t", t_filepath);
   // fprintf(fp, "%s\t", q_filepath);
   // fprintf(fp, "%f\t", scores->quad_vit);
   // fprintf(fp, "%f\t", scores->quad_fwd);
   // fprintf(fp, "%f\t", scores->quad_bck);
   // fprintf(fp, "%f\t", scores->naive_cloud_fwd);
   // fprintf(fp, "%f\t", scores->naive_cloud_bck);
   // fprintf(fp, "%f\t", alpha);
   // fprintf(fp, "%d\t", beta);
   // fprintf(fp, "%f\t", scores->perc_cells);
   // fprintf(fp, "%f\t", scores->perc_window);
   // fprintf(fp, "%d\t", Q);
   // fprintf(fp, "%d\n", T);

   // if ( fp != stdout ) fclose( fp );
}