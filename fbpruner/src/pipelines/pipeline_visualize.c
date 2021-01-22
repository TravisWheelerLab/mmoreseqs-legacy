/*******************************************************************************
 *  FILE:      pipeline_vizualize.c
 *  PURPOSE:   Vizualization Cloud Search Pipeline.
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
#include "../objects/structs.h"
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../parsers/parsers.h"
#include "../algs_linear/algs_linear.h"
#include "../algs_quad/algs_quad.h"
#include "../algs_naive/algs_naive.h"

/* header */
#include "pipelines.h"

/*
 *  FUNCTION:  vizualization_pipeline()
 *  SYNOPSIS:  Pipeline runs viterbi, forward/backward, and cloud search.  
 *             Output visualizations for python scripts.
 */
void vizualization_pipeline( WORKER* worker )
{
   printf("VIZ MODE:\n");
   #if DEBUG
   {
      fprintf(stdout, "Running vizualization...\n");
      debugger->verbose_level = VERBOSE_LOW;
   }
   #else
   {
      fprintf(stdout, "Test pipeline only supported in DEBUG build, not PRODUCTION build.  Recompile with -BUILD=DEBUG.\n" );
      exit(EXIT_FAILURE);
   }
   #endif 

   /* Commandline Arguments */
   FILE*          fp             = NULL;

   ARGS*          args           = worker->args;

   float          alpha          = args->alpha;
   float          beta           = args->beta;
   int            gamma          = args->gamma;

   char*          t_filepath     = args->t_filepath;
   char*          q_filepath     = args->q_filepath;

   int            t_filetype     = args->t_filetype;
   int            q_filetype     = args->q_filetype;

   char*          output_file    = args->output_filepath;

   CLOUD_PARAMS*   cloud_params  = &(worker->cloud_params);

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

   /* temporary cloud search objects */
   EDGEBOUND_ROWS*   edg_row_tmp   = EDGEBOUND_ROWS_Create();
   VECTOR_INT*       lb_vec[3];
   VECTOR_INT*       rb_vec[3];
   for (int i = 0; i < 3; i++ ) {
      VECTOR_INT_Create( lb_vec[i] );
      VECTOR_INT_Create( rb_vec[i] );
   }

   /* SCORES => stores result scores */
   TIMES*         times          = worker->times;
   NAT_SCORES*    scores         = worker->scores;

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
   printf("%*s: %f\n", pad, "BETA",             beta);
   printf("%*s: %d\n", pad, "GAMMA",            gamma);
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
   MATRIX_3D_Reuse_Clean( debugger->test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
   MATRIX_2D_Reuse_Clean( debugger->cloud_MX, Q+1, T+1 );

   /* run viterbi algorithm */
   {
      printf("=== VITERBI -> START ===\n");
      /* ==> viterbi (quadratic) */
      run_Viterbi_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, &sc);
      printf("Viterbi Score (quad):\t%f\n", sc);
      scores->quad_vit = sc;
      DP_MATRIX_Save(Q, T, st_MX_quad, sp_MX_quad, "test_output/my.viterbi.quad.mx");
      printf("=== VITERBI -> END ===\n\n");
   }

   /* run traceback of viterbi */
   {
      printf("=== TRACEBACK -> START ===\n");
      /* ==> viterbi (quadratic) */
      run_Traceback_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr);
      ALIGNMENT_Save(tr, "test_output/my.traceback.tsv");
      printf("=== TRACEBACK -> END ===\n\n");
   }

   /* run forward/backward algorithms */
   logsum_Init();

   /* run forward */
   {
      printf("=== FORWARD -> START ===\n");
      /* ==> forward (quadratic) */
      run_Forward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, &sc);
      printf("Forward Score (quad):\t%f\n", sc);
      scores->quad_fwd = sc;
      DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.fwd.quad.mx");
      /* ==> forward (linear) */
      run_Forward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, &sc);
      printf("Forward Score  (lin):\t%f\n", sc);
      scores->lin_fwd = sc;
      MATRIX_3D_Copy( debugger->test_MX, st_MX_lin );
      DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.fwd.lin.mx");
      printf("=== FORWARD -> END ===\n\n");
   }

   /* run backward */
   {
      // printf("=== BACKWARD -> START ===\n");
      // /* ==> backward (quadratic) */
      // run_Backward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, &sc);
      // printf("Backward Score (quad):\t%f\n", sc);
      // scores->quad_bck = sc;
      // DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.bck.quad.mx");
      /* ==> backward (linear) */
      // run_Backward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, &sc);
      // printf("Backward Score  (lin):\t%f\n", sc);
      // scores->lin_bck = sc;
      // MATRIX_3D_Copy( debugger->test_MX, st_MX_lin );
      // DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.bck.lin.mx");
      // printf("=== BACKWARD -> END ===\n\n");
   }

   /* need to clean matrix after fwd/bck */
   DP_MATRIX_Fill( Q, T, st_MX3_lin, sp_MX_lin, -INF );

   /* run cloud forward */
   {
      /* cloud forward (quadratic) */
      printf("=== CLOUD FORWARD (quadratic) -> START ===\n");
      run_Cloud_Forward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr, edg_row_tmp, lb_vec, rb_vec, edg_fwd_quad, cloud_params);
      if ( debugger->verbose_level >= VERBOSE_ALL ) {
         MATRIX_2D_Copy( cloud_MX_quad, debugger->cloud_MX );
         DP_MATRIX_VIZ_Dump( cloud_MX_quad, debugout );
      }
      DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.cloud_fwd.quad.mx");
      EDGEBOUNDS_Save(edg_fwd_quad, "test_output/my.cloud_fwd.quad.diags.edg");
      printf("=== CLOUD FORWARD (quadratic) -> END ===\n\n");
      /* cloud forward (linear) */
      printf("=== CLOUD FORWARD (linear) -> START ===\n");
      printf("# OUTSIDE: alpha = %9.4f, beta = %9.4f\n", worker->cloud_params.alpha, worker->cloud_params.beta );
      run_Cloud_Forward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, tr, edg_row_tmp, edg_fwd_lin, cloud_params);
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
      #if ( CLOUD_METHOD == CLOUD_DIAGS )
      {
         DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.cloud_fwd.lin.diags.mx");
         EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_fwd.lin.diags.edg");
      }
      #elif ( CLOUD_METHOD == CLOUD_ROWS )
      {
         DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.cloud_fwd.lin.rows.mx");
         EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_fwd.lin.rows.edg");
      }
      #endif
      printf("=== CLOUD FORWARD (linear) -> END ===\n\n");
   }

   /* run cloud backward */
   {
      /* cloud backward (quadratic) */
      printf("=== CLOUD BACKWARD (quadratic) -> START ===\n");
      run_Cloud_Backward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr, edg_row_tmp, lb_vec, rb_vec, edg_bck_quad, cloud_params);
      if ( debugger->verbose_level >=  VERBOSE_ALL ) {
         MATRIX_2D_Copy( cloud_MX_quad, debugger->cloud_MX );
         DP_MATRIX_VIZ_Dump( cloud_MX_quad, debugout );
      }
      DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.cloud_bck.quad.mx");
      EDGEBOUNDS_Save(edg_bck_quad, "test_output/my.cloud_bck.quad.diags.edg");
      printf("=== CLOUD BACKWARD (quadratic) -> END ===\n\n");
      /* cloud backward (linear) */
      printf("=== CLOUD BACKWARD (linear) -> START ===\n");
      run_Cloud_Backward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, tr, edg_row_tmp, edg_bck_lin, cloud_params);
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
      #if ( CLOUD_METHOD == CLOUD_DIAGS )
      DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.cloud_bck.lin.diags.mx");
      EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_bck.lin.diags.edg");
      #elif ( CLOUD_METHOD == CLOUD_ROWS )
      DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.cloud_bck.lin.rows.mx");
      EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_bck.lin.rows.edg");
      #endif
      printf("=== CLOUD BACKWARD (linear) -> END ===\n\n");
   }

   /* set edgebounds equivalent for merge/reorient testing */
   edg_fwd_quad = EDGEBOUNDS_Copy( edg_fwd_quad, edg_fwd_lin );
   edg_bck_quad = EDGEBOUNDS_Copy( edg_bck_quad, edg_bck_lin );
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
      EDGEBOUNDS_Union(Q, T, edg_fwd_lin, edg_bck_lin, edg_diag_lin);
      EDGEBOUNDS_Reorient_to_Row(Q, T, edg_diag_lin, edg_row_lin);
      EDGEBOUNDS_Save(edg_diag_lin, "test_output/my.cloud.lin.diags.edg");
      EDGEBOUNDS_Save(edg_row_lin, "test_output/my.cloud.lin.rows.edg");
      int edg_cmp_lin = ( EDGEBOUNDS_Compare_by_Cloud( edg_diag_lin, cloud_MX_lin, edg_row_lin, cloud_MX_diff ) == 0 );
      printf("Rows vs Diags:\tEDGES?\t\t%s\n", ( edg_cmp_lin == true ) ? "PASS" : "FAIL" );
      printf("=== MERGE & REORIENT CLOUD (linear) -> END ===\n\n");
   }
   /* fill cloud for naive algs */
   MATRIX_2D_Fill( cloud_MX_naive, 0 );
   MATRIX_2D_Cloud_Fill( cloud_MX_naive, edg_row_lin, 1 );
   // DP_MATRIX_VIZ_Dump( cloud_MX_naive, stdout );
   // EDGEBOUNDS_Dump( edg_row_lin, stdout );

   /* bounded forward */
   printf("=== BOUNDED FORWARD -> START ===\n");
   /* bounded forward (naive) */
   // run_Bound_Forward_Naive(q_seq, t_prof, Q, T, st_MX_naive, sp_MX_naive, cloud_MX_naive, &sc);
   // printf("Bounded Forward Score (naive):\t%f\n", sc);
   // scores->naive_cloud_fwd = sc;
   // DP_MATRIX_Trace_Save(Q, T, st_MX_naive, sp_MX_naive, tr, "test_output/my.bound_fwd.naive.mx");
   // /* bounded forward (quadratic) */
   // run_Bound_Forward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, edg_row_lin, &sc);
   // printf("Bounded Forward Score  (quad):\t%f\n", sc);
   // scores->quad_cloud_fwd = sc;
   // DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.bound_fwd.quad.mx");
   /* bounded forward (comparison) */
   run_Bound_Forward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, edg_row_lin, &sc);
   printf("Bounded Forward Score   (lin):\t%f\n", sc);
   scores->lin_bound_fwd = sc;
   MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
   DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.bound_fwd.lin.mx");
   printf("=== BOUNDED FORWARD -> END ===\n\n");

   // /* bounded backward */
   // printf("=== BOUNDED BACKWARD -> START ===\n");
   // /* bounded backward (naive) */
   // run_Bound_Backward_Naive(q_seq, t_prof, Q, T, st_MX_naive, sp_MX_naive, cloud_MX_naive, &sc);
   // printf("Bounded Backward Score (naive):\t%f\n", sc);
   // scores->naive_cloud_bck = sc;
   // DP_MATRIX_Trace_Save(Q, T, st_MX_naive, sp_MX_naive, tr, "test_output/my.bound_bck.naive.mx");
   // /* bounded backward (quadratic) */
   // run_Bound_Backward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, edg_row_lin, &sc);
   // printf("Bounded Backward Score  (quad):\t%f\n", sc);
   // scores->quad_cloud_fwd = sc;
   // DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.bound_bck.quad.mx");
   // /* bounded backward (linear) */
   // run_Bound_Backward_Linear(q_seq, t_prof, Q, T, st_MX_lin, sp_MX_lin, edg_row_lin, &sc);
   // printf("Bounded Backward Score   (lin):\t%f\n", sc);
   // scores->lin_cloud_fwd = sc;
   // MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
   // DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.bound_bck.lin.mx");
   // printf("=== BOUNDED BACKWARD -> END ===\n\n");


   /* output results */
   STRING_Replace( t_prof->name, ' ', '_' );
   STRING_Replace( q_seq->name, ' ', '_' );
   printf("## %s %s %s %s %5.1f %5.1f %d %9.4f %9.4f %9.4f\n", 
      t_filepath, q_filepath, t_prof->name, q_seq->name,
      alpha, beta, gamma,
      scores->quad_vit, scores->quad_fwd, scores->lin_bound_fwd);


   /* CLEAN-UP */
   /* PROFILE & SEQUENCE */
   t_prof         = HMM_PROFILE_Destroy( t_prof );
   t_seq          = SEQUENCE_Destroy( t_seq );  /* only used if target is a fasta file */
   q_seq          = SEQUENCE_Destroy( q_seq );

   /* EDGEBOUNDS & ALIGNMENTS */
   tr             = ALIGNMENT_Destroy( tr );

   edg_fwd_lin    = EDGEBOUNDS_Destroy( edg_fwd_lin );
   edg_bck_lin    = EDGEBOUNDS_Destroy( edg_bck_lin );
   edg_row_lin    = EDGEBOUNDS_Destroy( edg_row_lin );
   edg_diag_lin   = EDGEBOUNDS_Destroy( edg_diag_lin );

   edg_fwd_quad   = EDGEBOUNDS_Destroy( edg_fwd_quad );
   edg_bck_quad   = EDGEBOUNDS_Destroy( edg_bck_quad );
   edg_row_quad   = EDGEBOUNDS_Destroy( edg_row_quad );
   edg_diag_quad  = EDGEBOUNDS_Destroy( edg_diag_quad );

   /* temporary edgebound object for storing row-wise edgebounds during cloud search */
   edg_row_tmp   = EDGEBOUND_ROWS_Destroy( edg_row_tmp );

   /* te memory for quadratic algs */
   st_MX_naive    = MATRIX_3D_Destroy( st_MX_naive );
   sp_MX_naive    = MATRIX_2D_Destroy( sp_MX_naive );
   cloud_MX_naive = MATRIX_2D_Destroy( cloud_MX_naive );
   /* allocate memory for quadratic algs */
   st_MX_quad     = MATRIX_3D_Destroy( st_MX_quad );
   sp_MX_quad     = MATRIX_2D_Destroy( sp_MX_quad );
   cloud_MX_quad  = MATRIX_2D_Destroy( cloud_MX_quad );
   /* allocate memory for comparing row-wise algs */
   cloud_MX_rows  = MATRIX_2D_Destroy( cloud_MX_rows );
   /* allocate memory for linear algs */
   st_MX_lin      = MATRIX_3D_Destroy( st_MX_lin );
   st_MX3_lin     = MATRIX_3D_Destroy( st_MX3_lin );
   sp_MX_lin      = MATRIX_2D_Destroy( sp_MX_lin );
   cloud_MX_lin   = MATRIX_2D_Destroy( cloud_MX_lin );
   /* allocate memory for testing */
   st_MX_diff     = MATRIX_3D_Destroy( st_MX_diff );
   st_MX3_diff    = MATRIX_3D_Destroy( st_MX3_diff );
   sp_MX_diff     = MATRIX_2D_Destroy( sp_MX_diff );
   cloud_MX_diff  = MATRIX_2D_Destroy( cloud_MX_diff );

   /* allocate full search cloud */
   full_cloud_MX  = MATRIX_2D_Destroy( full_cloud_MX );
   full_cloud_edg = EDGEBOUNDS_Destroy( full_cloud_edg );

   WORK_cleanup( worker );
}