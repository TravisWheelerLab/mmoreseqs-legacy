/*******************************************************************************
 *  FILE:      pipeline_int_test.c
 *  PURPOSE:   Test Cloud Search Pipeline.
 *             Requires -BUILD=DEBUG for full functionality.
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
#include "algs_sparse.h"

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

   /* vars for testing and comparing */
   int   cmp;
   int   dp_cmp;
   int   edg_cmp;
   int   cloud_cmp;
   int   row_cmp;
   int   diag_cmp;

   /* Commandline Arguments */
   FILE*    fp             = NULL;
   ARGS*    args           = worker->args;

   float    alpha          = args->alpha;
   float    beta           = args->beta;
   int      gamma          = args->gamma;

   char*    t_filepath     = args->t_filepath;
   char*    q_filepath     = args->q_filepath;

   int      t_filetype     = args->t_filetype;
   int      q_filetype     = args->q_filetype;

   char*    output_file    = args->output_filepath;

   /* TODO: multi mode is unneccesary (only support for UNILOCAL) */
   /* mode types: EDG_NONE, MODE_MULTILOCAL, MODE_MULTIGLOCAL, MODE_UNILOCAL, MODE_UNIGLOCAL */

   /* Cloud Search mode (prohibits jumps) */
   int            mode           = MODE_UNILOCAL;      

   /* PROFILE & SEQUENCE */
   HMM_PROFILE*   t_prof         = HMM_PROFILE_Create();
   SEQUENCE*      t_seq          = SEQUENCE_Create();  /* only used if target is a fasta file */
   SEQUENCE*      q_seq          = SEQUENCE_Create();

   /* EDGEBOUNDS & ALIGNMENTS */
   ALIGNMENT*     tr_quad        = ALIGNMENT_Create();
   ALIGNMENT*     tr_quad_2      = ALIGNMENT_Create();
   ALIGNMENT*     tr_quad_3      = ALIGNMENT_Create();
   ALIGNMENT*     tr_lin         = ALIGNMENT_Create();
   ALIGNMENT*     tr_sparse      = ALIGNMENT_Create();
   ALIGNMENT*     tr             = tr_quad;

   EDGEBOUNDS*    edg_fwd_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_bck_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_row_lin    = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_diag_lin   = EDGEBOUNDS_Create();

   EDGEBOUNDS*    edg_fwd_quad   = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_bck_quad   = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_row_quad   = EDGEBOUNDS_Create();
   EDGEBOUNDS*    edg_diag_quad  = EDGEBOUNDS_Create();

   /* temporary edgebound object for storing row-wise edgebounds during cloud search */
   EDGEBOUND_ROWS*   edg_row_tmp    = EDGEBOUND_ROWS_Create();
   CLOUD_PARAMS*     cloud_params   = &(worker->cloud_params);

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

   /* worker params */
   cloud_params->alpha  = alpha;
   cloud_params->beta   = beta;
   cloud_params->gamma  = gamma;

   /* PRINT ARGS */
   int pad = 20;
   printf("%*s: %s\n",     pad, "MODE",               MODE_NAMES[mode]);
   printf("%*s: %s\n",     pad, "HMM_FILENAME",       t_filepath);
   printf("%*s: %s\n",     pad, "FASTA_FILENAME",     q_filepath);
   printf("%*s: %.3f\n",   pad, "ALPHA",              cloud_params->alpha);
   printf("%*s: %.3f\n",   pad, "BETA",               cloud_params->beta);
   printf("%*s: %d\n",     pad, "GAMMA",              cloud_params->gamma);
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
   SEQUENCE_Dump( q_seq, stdout );

   /* build t_prof profile */
   printf("loading target...\n");
   if ( t_filetype == FILE_HMM || true ) 
   {
      HMM_PROFILE_Parse( t_prof, t_filepath, 0 );
      HMM_PROFILE_Convert_NegLog_To_Real( t_prof );
      HMM_PROFILE_Config( t_prof, mode );
      HMM_PROFILE_Dump( t_prof, stdout );
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
   HMM_PROFILE_Dump( t_prof, stdout );

   printf("TARGET LEN:\t%d\n", t_prof->N);
   printf(" QUERY LEN:\t%d\n", q_seq->N);
   printf("=== BUILD HMM_PROFILE / QUERY -> END ===\n\n");

   /* initialize logsum table */
   logsum_Init();

   /* get dimensions */
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
   MATRIX_3D*           st_MX_naive       = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES,  Q+1, T+1 );
   MATRIX_2D*           sp_MX_naive       = MATRIX_2D_Create_Clean( NUM_SPECIAL_STATES, Q+1 );
   MATRIX_2D*           cloud_MX_naive    = MATRIX_2D_Create_Clean( Q+1, T+1 );
   /* allocate memory for quadratic algs */
   MATRIX_3D*           st_MX_quad        = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES,  Q+1, T+1 );
   MATRIX_2D*           sp_MX_quad        = MATRIX_2D_Create_Clean( NUM_SPECIAL_STATES, Q+1 );
   MATRIX_2D*           cloud_MX_quad     = MATRIX_2D_Create_Clean( Q+1, T+1 );
   /* allocate memory for comparing row-wise algs */
   MATRIX_2D*           cloud_MX_rows     = MATRIX_2D_Create_Clean( Q+1, T+1 );
   /* allocate memory for linear algs */
   MATRIX_3D*           st_MX3_lin        = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
   MATRIX_3D*           st_MX_lin         = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES,  Q+1, T+1 );
   MATRIX_2D*           sp_MX_lin         = MATRIX_2D_Create_Clean( NUM_SPECIAL_STATES, Q+1 );
   MATRIX_2D*           cloud_MX_lin      = MATRIX_2D_Create_Clean( Q+1, T+1 );
   /* allocate memory for sparse algs */
   MATRIX_3D_SPARSE*    st_SMX            = MATRIX_3D_SPARSE_Create(  );
   MATRIX_3D*           st_MX_sparse      = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES,  Q+1, T+1 );
   MATRIX_2D*           sp_MX_sparse      = MATRIX_2D_Create_Clean( NUM_SPECIAL_STATES, Q+1 );
   MATRIX_2D*           cloud_MX_sparse   = MATRIX_2D_Create_Clean( Q+1, T+1 );
   /* allocate memory for testing */
   MATRIX_3D*           st_MX_diff        = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES,  Q+1, T+1 );
   MATRIX_3D*           st_MX3_diff       = MATRIX_3D_Create_Clean( NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
   MATRIX_2D*           sp_MX_diff        = MATRIX_2D_Create_Clean( NUM_SPECIAL_STATES, Q+1 );
   MATRIX_2D*           cloud_MX_diff     = MATRIX_2D_Create_Clean( Q+1, T+1 );

   /* allocate full search cloud for  */
   MATRIX_2D*           full_cloud_MX     = MATRIX_2D_Create_Clean( Q+1, T+1 );
   EDGEBOUNDS*          full_cloud_edg    = EDGEBOUNDS_Create();

   /* vectors for edgebounds in cloud search */
   VECTOR_INT* lb_vec[3];
   VECTOR_INT* rb_vec[3];
   for ( int i = 0; i < 3; i++ ) {
      lb_vec[i] = VECTOR_INT_Create();
      rb_vec[i] = VECTOR_INT_Create();
   }

   /* debug matrix */
   {
      debugger->cloud_MX   = MATRIX_2D_Create_Clean( 1, 1 );
      debugger->cloud_MX3  = MATRIX_2D_Create_Clean( 1, 1 );
      debugger->test_MX    = MATRIX_3D_Create_Clean( 1, 1, 1 );
      debugger->test_MX3   = MATRIX_3D_Create_Clean( 1, 1, 1 );
      debugger->test_edg   = EDGEBOUNDS_Create();

      MATRIX_3D_Reuse_Clean( debugger->test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_2D_Reuse_Clean( debugger->cloud_MX, Q+1, T+1 );
      MATRIX_3D_Reuse_Clean( debugger->test_MX3, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_2D_Reuse_Clean( debugger->cloud_MX3, Q+1, T+1 );
   }

   /* build edgebounds which cover the entire search space (to compare against full fwd/back) */
   {
      MATRIX_2D_Fill( full_cloud_MX, 1 );
      EDGEBOUNDS_Build_From_Cloud( Q, T, full_cloud_edg, full_cloud_MX, EDG_ROW );
   }

   /* build sparse matrix from full edgebounds */
   {
      MATRIX_3D_SPARSE_Shape_Like_Edgebounds( st_SMX, full_cloud_edg );
      // EDGEBOUNDS_Dump( st_SMX->edg_outer, stdout );
      EDGEBOUNDS_Save( st_SMX->edg_outer, "test_output/my.fwd.full.sparse.outer.edg" );
      EDGEBOUNDS_Save( st_SMX->edg_inner, "test_output/my.fwd.full.sparse.inner.edg" );
   }

   /* run viterbi algorithm */
   {
      printf("=== FULL VITERBI -> START ===\n");
      /* ==> viterbi (quadratic) */
      printf("==> viterbi quadratic\n");
      run_Viterbi_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, &sc);
      printf("Viterbi Score (quad):\t%f\n", sc);
      scores->quad_vit = sc;
      DP_MATRIX_Save(Q, T, st_MX_quad, sp_MX_quad, "test_output/my.viterbi.quad.mx");
      /* ==> viterbi (linear) */
      printf("==> viterbi linear\n");
      run_Viterbi_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, &sc);
      printf("Viterbi Score (lin):\t%f\n", sc);
      scores->quad_vit = sc;
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
      DP_MATRIX_Save(Q, T, st_MX_quad, sp_MX_quad, "test_output/my.viterbi.lin.mx");
      /* ==> viterbi (comparison) */
      printf("==> viterbi comparison: linear v. quadratic\n");
      dp_cmp  = DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin );
      cmp = ( dp_cmp == 0 ) ? true : false;
      printf("MATRIX VALUES?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      if ( cmp == false ) {
         printf("viterbi (quad):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_quad, stdout );
         DP_MATRIX_Dump( Q, T, st_MX_quad, sp_MX_quad, stdout );
         printf("viterbi (lin):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_lin, stdout );
         DP_MATRIX_Dump( Q, T, st_MX_lin, sp_MX_quad, stdout );
      }
      printf("=== FULL VITERBI -> END ===\n\n");
   }

   /* run traceback of viterbi */
   {
      printf("=== TRACEBACK -> START ===\n");
      /* ==> traceback (quadratic) */
      printf("==> traceback quadratic (via HMMER)\n");
      ALIGNMENT_Reuse(tr_quad, Q, T);
      run_Traceback_Quad_via_hmmer(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr_quad);
      ALIGNMENT_Save(tr_quad, "test_output/my.traceback.quad.hmmer.tsv");
      /* ==> traceback (quadratic via compare) */
      printf("==> traceback quadratic (via compare)\n");
      ALIGNMENT_Reuse(tr_quad_2, Q, T);
      run_Traceback_Quad_via_cmp(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr_quad_2);
      ALIGNMENT_Save(tr_quad_2, "test_output/my.traceback.quad.cmp.tsv");
      /* ==> traceback comparison */
      cmp = ( ALIGNMENT_Compare( tr_quad, tr_quad_2 ) == 0 ) ? true : false;
      printf("ALIGNMENTS (HMMER vs Compare)?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      if ( cmp == false ) {
         printf("traceback quadratic (via hmmer):\n");
         ALIGNMENT_Dump(tr_quad, stdout);
         printf("traceback quadratic (via compare):\n");
         ALIGNMENT_Dump(tr_quad_2, stdout);
      }
      /* ==> traceback (quadratic via maximum) */
      printf("==> traceback quadratic (via max)\n");
      ALIGNMENT_Reuse(tr_quad_3, Q, T);
      run_Traceback_Quad_via_max(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr_quad_3);
      ALIGNMENT_Save(tr_quad_3, "test_output/my.traceback.quad.max.tsv");
      /* ==> traceback comparison */
      cmp = ( ALIGNMENT_Compare( tr_quad, tr_quad_3 ) == 0 ) ? true : false;
      printf("ALIGNMENTS (HMMER vs Max)?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      if ( cmp == false ) {
         printf("traceback quadratic (via hmmer):\n");
         ALIGNMENT_Dump(tr_quad, stdout);
         printf("traceback quadratic (via max):\n");
         ALIGNMENT_Dump(tr_quad_3, stdout);
      }
      /* ==> traceback (linear) */
      printf("==> traceback linear\n");
      /* TODO: implment and test linear traceback */
      /* ==> print begin and end points in traceback */
      TRACE* beg = &(tr->traces->data[tr->beg]);
      TRACE* end = &(tr->traces->data[tr->end]);
      printf("START: (%d,%d) -> END: (%d,%d)\n", beg->i, beg->j, end->i, end->j);
      printf("=== TRACEBACK -> END ===\n\n");
   }

   DP_MATRIX_Clean(Q, T, st_MX3_lin, sp_MX_lin);

   /* run sparse forward (compare to full forward) */
   {
      printf("=== FULL VITERBI SPARSE -> START ===\n");
      /* ==> forward (quadratic) */
      printf("==> bound viterbi sparse\n");
      run_Bound_Viterbi_Sparse(q_seq, t_prof, Q, T, st_SMX, sp_MX_sparse, st_SMX->edg_inner, &sc); 
      printf("Bound Viterbi Score (sparse):\t%f\n", sc);
      scores->quad_cloud_fwd = sc;
      MATRIX_3D_Copy( st_MX_sparse, debugger->test_MX );
      DP_MATRIX_Save(Q, T, st_MX_sparse, sp_MX_sparse, "test_output/my.vit.full.sparse.mx");
      /* ==> forward (comparison) */
      printf("==> bound forward comparison: sparse v. normal matrix\n");
      dp_cmp = 0;
      dp_cmp += DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_sparse, sp_MX_sparse );
      cmp = ( dp_cmp == 0 ) ? true : false;
      if ( cmp == false ) {
         printf("bound forward (sparse):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_sparse, stdout );
         printf("bound forward (quad):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_quad, stdout );
      }
      printf("MATRIX VALUES?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      printf("=== FULL VITERBI SPARSE -> END ===\n\n");
   }

   /* run traceback of viterbi */
   {
      printf("=== TRACEBACK SPARSE -> START ===\n");
      /* ==> traceback (quadratic) */
      printf("==> traceback sparse\n");
      ALIGNMENT_Reuse(tr_quad, Q, T);
      run_Traceback_Sparse(q_seq, t_prof, Q, T, st_SMX, sp_MX_sparse, st_SMX->edg_inner, tr_sparse);
      ALIGNMENT_Save(tr_sparse, "test_output/my.traceback.sparse.tsv");
      cmp = ( ALIGNMENT_Compare( tr_quad, tr_sparse ) == 0 ) ? true : false;
      printf("ALIGNMENTS (quad vs sparse)?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      if ( cmp == false ) {
         printf("traceback quadratic:\n");
         ALIGNMENT_Dump(tr_quad, stdout);
         printf("traceback sparse:\n");
         ALIGNMENT_Dump(tr_sparse, stdout);
      }
      TRACE* beg = &(tr->traces->data[tr->beg]);
      TRACE* end = &(tr->traces->data[tr->end]);
      printf("START: (%d,%d) -> END: (%d,%d)\n", beg->i, beg->j, end->i, end->j);
      printf("=== TRACEBACK SPARSE -> END ===\n\n");
   }

   /* run forward */
   {
      printf("=== FORWARD -> START ===\n");
      /* ==> forward (quadratic) */
      printf("==> forward quadratic\n");
      run_Forward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, &sc);
      printf("Forward Score (quad):\t%f\n", sc);
      scores->quad_fwd = sc;
      DP_MATRIX_Save(Q, T, st_MX_quad, sp_MX_quad, "test_output/my.fwd.quad.mx");
      /* ==> forward (linear) */
      printf("==> forward linear\n");
      run_Forward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, &sc);
      printf("Forward Score  (lin):\t%f\n", sc);
      scores->lin_fwd = sc;  
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX ); 
      DP_MATRIX_Save(Q, T, st_MX_lin, sp_MX_lin, "test_output/my.fwd.lin.mx");
      /* ==> forward (comparison) */
      printf("==> forward comparison: quadratic v. linear\n");
      dp_cmp  = DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin );
      cmp = ( dp_cmp == 0 ) ? true : false;
      if ( cmp == false ) {
         printf("forward (quad):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_quad, stdout );
         printf("forward (lin):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_lin, stdout );
      }
      printf("MATRIX VALUES?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      printf("=== FORWARD -> END ===\n\n");
   }

   DP_MATRIX_Clean(Q, T, st_MX3_lin, sp_MX_lin);

   /* run bounded forward (compare to full forward) */
   {
      printf("=== FULL FORWARD BOUND -> START ===\n");
      /* ==> forward (quadratic) */
      printf("==> bound forward quadratic\n");
      DP_MATRIX_Clean(Q, T, st_MX_quad, sp_MX_quad);
      run_Bound_Forward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, full_cloud_edg, &sc);
      printf("Bound Forward Score (quad):\t%f\n", sc);
      scores->quad_cloud_fwd = sc;
      DP_MATRIX_Save(Q, T, st_MX_quad, sp_MX_quad, "test_output/my.fwd.full.cloud.quad.mx");
      /* ==> forward (linear) */
      printf("==> bound forward linear\n");
      DP_MATRIX_Clean(Q, T, st_MX3_lin, sp_MX_lin);
      run_Bound_Forward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, full_cloud_edg, &sc);
      printf("Bound Forward Score  (lin):\t%f\n", sc);
      scores->lin_cloud_fwd = sc;  
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX ); 
      DP_MATRIX_Save(Q, T, st_MX_lin, sp_MX_lin, "test_output/my.fwd.full.cloud.lin.mx");
      /* ==> forward (comparison) */
      printf("==> bound forward comparison: quadratic v. linear\n");
      dp_cmp = 0;
      dp_cmp += DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin );
      cmp = ( dp_cmp == 0 ) ? true : false;
      if ( cmp == false ) {
         printf("bound forward (quad):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_quad, stdout );
         printf("bound forward (lin):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_lin, stdout );
      }
      printf("MATRIX VALUES?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      printf("=== FULL FORWARD BOUND -> END ===\n\n");
   }

   // DP_MATRIX_Clean(Q, T, st_MX3_lin, sp_MX_lin);

   /* run sparse forward (compare to full forward) */
   {
      printf("=== FULL FORWARD SPARSE -> START ===\n");
      /* ==> forward (quadratic) */
      printf("==> bound forward sparse\n");

      run_Bound_Forward_Sparse(q_seq, t_prof, Q, T, st_SMX, sp_MX_lin, st_SMX->edg_inner, &sc); 
      printf("Bound Forward Score (sparse):\t%f\n", sc);
      scores->quad_cloud_fwd = sc;
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
      DP_MATRIX_Save(Q, T, st_MX_lin, sp_MX_lin, "test_output/my.fwd.full.sparse.mx");
      /* ==> forward (comparison) */
      printf("==> bound forward comparison: sparse v. normal matrix\n");
      dp_cmp = 0;
      dp_cmp += DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin );
      cmp = ( dp_cmp == 0 ) ? true : false;
      if ( cmp == false ) {
         printf("bound forward (sparse):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_lin, stdout );
         printf("bound forward (quad):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_quad, stdout );
      }
      printf("MATRIX VALUES?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      printf("=== FULL FORWARD SPARSE -> END ===\n\n");
   }

   /* run backward */
   {
      printf("=== BACKWARD -> START ===\n");
      /* ==> backward (quadratic) */
      printf("==> backward quadratic\n");
      run_Backward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, &sc);
      printf("Backward Score (quad):\t%f\n", sc);
      scores->quad_bck = sc;
      DP_MATRIX_Save(Q, T, st_MX_quad, sp_MX_quad, "test_output/my.bck.quad.mx");
      /* ==> backward (linear) */
      printf("==> backward linear\n");
      run_Backward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, &sc);
      printf("Backward Score  (lin):\t%f\n", sc);
      scores->lin_bck = sc;
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
      DP_MATRIX_Save(Q, T, st_MX_lin, sp_MX_lin, "test_output/my.bck.lin.mx");
      /* ==> abackward (comparison) */
      printf("==> backward comparison: quadratic v. linear\n");
      dp_cmp  = DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin );
      cmp = ( dp_cmp == 0 ) ? true : false;
      if ( cmp == false ) {
         printf("backward (quad):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_quad, stdout );
         printf("backward (lin):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_lin, stdout );
      }
      printf("MATRIX VALUES?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      printf("=== BACKWARD -> END ===\n\n");
   }

   /* run bounded backward (compare to full backward) */
   {
      printf("=== FULL BACKWARD BOUND -> START ===\n");
      /* ==> bound backward (quadratic) */
      printf("==> bound forward quadratic\n");
      DP_MATRIX_Clean(Q, T, st_MX_quad, sp_MX_quad);
      run_Bound_Backward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, full_cloud_edg, &sc);
      printf("Bound Backward Score (quad):\t%f\n", sc);
      scores->quad_cloud_fwd = sc;
      DP_MATRIX_Save(Q, T, st_MX_quad, sp_MX_quad, "test_output/my.bck.cloud.quad.mx");
      /* ==> backward (linear) */
      printf("==> bound forward linear\n");
      DP_MATRIX_Clean(Q, T, st_MX3_lin, sp_MX_lin);
      run_Bound_Backward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, full_cloud_edg, &sc);
      printf("Bound Backward Score  (lin):\t%f\n", sc);
      scores->lin_cloud_fwd = sc;  
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX ); 
      DP_MATRIX_Save(Q, T, st_MX_lin, sp_MX_lin, "test_output/my.bck.cloud.lin.mx");
      /* ==> backward (comparison) */
      printf("==> bound backward comparison: quadratic v. linear\n");
      dp_cmp = 0;
      dp_cmp += DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin );
      cmp = ( dp_cmp == 0 ) ? true : false;
      if ( cmp == false ) {
         printf("bound backward (quad):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_quad, stdout );
         printf("bound backward (lin):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_lin, stdout );
      }
      printf("MATRIX VALUES?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      printf("=== FULL BACKWARD BOUND -> END ===\n\n");
   }

   /* run sparse forward (compare to full forward) */
   {
      printf("=== FULL BACKWARD SPARSE -> START ===\n");
      /* ==> forward (quadratic) */
      printf("==> bound backward sparse\n");
      run_Bound_Backward_Sparse(q_seq, t_prof, Q, T, st_SMX, sp_MX_lin, st_SMX->edg_inner, &sc); 
      printf("Bound Backward Score (sparse):\t%f\n", sc);
      scores->quad_cloud_fwd = sc;
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
      DP_MATRIX_Save(Q, T, st_MX_lin, sp_MX_lin, "test_output/my.bck.full.sparse.mx");
      /* ==> forward (comparison) */
      printf("==> bound backward comparison: sparse v. normal matrix\n");
      dp_cmp = 0;
      dp_cmp += DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin );
      cmp = ( dp_cmp == 0 ) ? true : false;
      if ( cmp == false ) {
         printf("bound backward (sparse):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_lin, stdout );
         printf("bound backward (quad):\n");
         DP_MATRIX_MAT_Dump( Q, T, st_MX_quad, stdout );
      }
      printf("MATRIX VALUES?\t\t%s\n", cmp ? "PASS" : "FAIL" );
      printf("=== FULL BACKWARD SPARSE -> END ===\n\n");
   }

   /* need to clean matrix after fwd/bck */
   DP_MATRIX_Fill( Q, T, st_MX_quad, sp_MX_quad, -INF );
   // DP_MATRIX_Fill( Q, T, st_MX3_lin, sp_MX_lin, -INF );

   /* run cloud forward */
   {
      printf("=== CLOUD FORWARD -> START ===\n");
      /* cloud forward (quadratic) */
      printf("==> cloud forward quadratic\n");
      run_Cloud_Forward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr, edg_row_tmp, lb_vec, rb_vec, edg_fwd_quad, cloud_params );
      if ( debugger->verbose_level >= VERBOSE_ALL ) {
         MATRIX_2D_Copy( cloud_MX_quad, debugger->cloud_MX );
         DP_MATRIX_VIZ_Dump( cloud_MX_quad, stdout );
      }
      DP_MATRIX_Save(Q, T, st_MX_quad, sp_MX_quad, "test_output/my.cloud_fwd.quad.mx");
      EDGEBOUNDS_Save(edg_fwd_quad, "test_output/my.cloud_fwd.quad.diags.edg");
      /* cloud forward (linear) */
      printf("==> cloud forward linear\n");
      run_Cloud_Forward_Linear(
         q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, tr, edg_row_tmp, edg_fwd_lin, cloud_params );
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
      if ( debugger->verbose_level >= VERBOSE_ALL ) {
         MATRIX_2D_Copy( cloud_MX_lin, debugger->cloud_MX );
         DP_MATRIX_VIZ_Dump( cloud_MX_lin, stdout );
      }
      #if ( CLOUD_METHOD == CLOUD_DIAGS )
      {
         DP_MATRIX_Save(Q, T, st_MX_lin, sp_MX_lin, "test_output/my.cloud_fwd.lin.diags.mx");
         EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_fwd.lin.diags.edg");
      }
      #elif ( CLOUD_METHOD == CLOUD_ROWS )
      {
         DP_MATRIX_Save(Q, T, st_MX_lin, sp_MX_lin, "test_output/my.cloud_fwd.lin.rows.mx");
         EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_fwd.lin.rows.edg");
      }
      #endif
      /* cloud forward (comparison) */
      printf("==> cloud forward comparison: quadratic v. linear\n");
      /* compare dp matrix values */      
      cmp = ( DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin ) == 0 );
      printf("MATRIX VALUES?\t\t%s\n", ( cmp == true ) ? "PASS" : "FAIL" );
      if ( ( cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) ) {
         DP_MATRIX_Diff( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin, st_MX_diff, sp_MX_diff );
         fprintf( stdout, "=> DP MATRIX DIFF:\n" );
         DP_MATRIX_Dump( Q, T, st_MX_diff, sp_MX_diff, stdout );
      }  
      /* compare cloud */
      cmp = ( EDGEBOUNDS_Compare_by_Cloud( edg_fwd_quad, cloud_MX_quad, edg_fwd_lin, cloud_MX_lin ) == 0 );
      printf("CLOUD SHAPE?\t\t%s\n", ( cmp == true ) ? "PASS" : "FAIL" );
      if ( cmp == false ) {
         DP_MATRIX_VIZ_Compare( cloud_MX_diff, edg_fwd_lin, edg_fwd_quad );
         FILE* fp = fopen( "test_output/my.quadvlin.cloud_fwd.viz", "w+" );
         DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, fp );
         fclose(fp);
      }
      /* compare edgebounds */
      cmp = ( EDGEBOUNDS_Compare( edg_fwd_quad, edg_fwd_lin ) == 0 );
      printf("EDGEBOUNDS?\t\t%s\n", ( cmp == true ) ? "PASS" : "FAIL" );
      if ( ( cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) ) {
         fprintf( stdout, "=> LINEAR EDGEBOUNDS:\n" );
         EDGEBOUNDS_Dump( edg_fwd_lin, stdout );
         fprintf( stdout, "=> QUADRATIC EDGEBOUNDS:\n" );
         EDGEBOUNDS_Dump( edg_fwd_quad, stdout );
      }
      printf("=== CLOUD FORWARD -> END ===\n\n");
   }

   /* need to clean matrix after fwd/bck */
   DP_MATRIX_Fill( Q, T, st_MX_quad, sp_MX_quad, -INF );
   // DP_MATRIX_Fill( Q, T, st_MX3_lin, sp_MX_lin, -INF );

   /* run cloud backward */
   {
      printf("=== CLOUD BACKWARD -> START ===\n");
      /* cloud backward (quadratic) */
      printf("==> cloud backward quadratic\n");
      run_Cloud_Backward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, tr, edg_row_tmp, lb_vec, rb_vec, edg_bck_quad, cloud_params);
      if ( debugger->verbose_level >=  VERBOSE_ALL ) {
         MATRIX_2D_Copy( cloud_MX_quad, debugger->cloud_MX );
         DP_MATRIX_VIZ_Dump( cloud_MX_quad, stdout );
      }
      DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.cloud_bck.quad.mx");
      EDGEBOUNDS_Save(edg_bck_quad, "test_output/my.cloud_bck.quad.diags.edg");
      /* cloud backward (linear) */
      printf("==> cloud backward linear\n");
      run_Cloud_Backward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, tr, edg_row_tmp, edg_bck_lin, cloud_params );
      MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
      if ( debugger->verbose_level >=  VERBOSE_ALL ) {
         MATRIX_2D_Copy( cloud_MX_lin, debugger->cloud_MX );
         // DP_MATRIX_VIZ_Dump( cloud_MX_lin, stdout );
      }
      #if ( CLOUD_METHOD == CLOUD_DIAGS )
      {
         DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.cloud_bck.lin.diags.mx");
         EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_bck.lin.diags.edg");
      }
      #elif ( CLOUD_METHOD == CLOUD_ROWS )
      {
         DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.cloud_bck.lin.rows.mx");
         EDGEBOUNDS_Save(edg_bck_lin, "test_output/my.cloud_bck.lin.rows.edg");
      }
      #endif

      /* cloud backward (comparison) */
      printf("==> cloud backward comparison: quadratic v. linear\n");
      /* BUG: cloud_backward_linear() and cloud_backward_quad() have different scores at left edge */
      /* compare dp matrix values */
      cmp = ( DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin ) == 0 );
      printf("MATRIX VALUES?\t\t%s\n", ( cmp == true ) ? "PASS" : "FAIL" );
      if ( ( cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) ) {
         DP_MATRIX_Diff( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin, st_MX_diff, sp_MX_diff );
         fprintf( stdout, "=> DP MATRIX DIFF:\n" );
         DP_MATRIX_Dump( Q, T, st_MX_diff, sp_MX_diff, stdout );
      } 
      /* compare cloud cells */
      cmp = ( EDGEBOUNDS_Compare_by_Cloud( edg_bck_quad, cloud_MX_quad, edg_bck_lin, cloud_MX_lin ) == 0 );
      printf("CLOUD SHAPE?\t\t%s\n", ( cmp == true ) ? "PASS" : "FAIL" );
      if ( ( cmp == false ) ) {
         DP_MATRIX_VIZ_Compare( cloud_MX_diff, edg_bck_lin, edg_bck_quad );
         FILE* fp = fopen( "test_output/my.cloud_bck.quadvlin.viz", "w+" );
         DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, fp );
         fclose(fp);
      }
      /* compare edgebounds data */
      cmp = ( EDGEBOUNDS_Compare( edg_bck_quad, edg_bck_lin ) == 0 );
      printf("EDGEBOUNDS?\t\t%s\n", ( cmp == true ) ? "PASS" : "FAIL" );
      if ( ( cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) ) {
         fprintf( stdout, "=> LINEAR EDGEBOUNDS:\n");
         EDGEBOUNDS_Dump( edg_bck_lin, stdout );
         fprintf( stdout, "=> QUADRATIC EDGEBOUNDS:\n");
         EDGEBOUNDS_Dump( edg_bck_quad, stdout );
      }
      printf("=== CLOUD BACKWARD -> END ===\n\n");
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
      printf("=== MERGE & REORIENT CLOUD -> START ===\n");
     /* merge and reorient (naive) */
      printf("==> merge & reorient (naive)\n");
      EDGEBOUNDS_Merge_Reorient_Naive(Q, T, edg_fwd_quad, edg_bck_quad, edg_diag_quad, edg_row_quad, cloud_MX_quad);
      EDGEBOUNDS_Save(edg_diag_quad, "test_output/my.cloud.quad.diags.edg");
      EDGEBOUNDS_Save(edg_row_quad, "test_output/my.cloud.quad.rows.edg");
      cmp = ( EDGEBOUNDS_Compare_by_Cloud( edg_diag_quad, cloud_MX_quad, edg_row_quad, cloud_MX_diff ) == 0 );
      printf("Rows vs Diags:\tEDGES?\t\t%s\n", ( cmp == true ) ? "PASS" : "FAIL" );
     /* merge and reorient (naive) */
      printf("==> merge & reorient (naive)\n");
      EDGEBOUNDS_Merge_Together(Q, T, edg_fwd_lin, edg_bck_lin, edg_diag_lin);
      EDGEBOUNDS_Reorient_to_Row(Q, T, edg_diag_lin, edg_row_lin);
      EDGEBOUNDS_Save(edg_diag_lin, "test_output/my.cloud.lin.diags.edg");
      EDGEBOUNDS_Save(edg_row_lin, "test_output/my.cloud.lin.rows.edg");
      cmp = ( EDGEBOUNDS_Compare_by_Cloud( edg_diag_lin, cloud_MX_lin, edg_row_lin, cloud_MX_diff ) == 0 );
      printf("Rows vs Diags:\tEDGES?\t\t%s\n", ( cmp == true ) ? "PASS" : "FAIL" );
      /* test post-merge */
      printf("==> merge & reorient comparison: naive v. linear\n");
      int cmp_1 = ( EDGEBOUNDS_Compare_by_Cloud( edg_diag_lin, cloud_MX_lin, edg_row_quad, cloud_MX_quad ) == 0 );
      printf("POST-MERGE?\t\t%s\n", ( cmp_1 == true ) ? "PASS" : "FAIL" );
      if ( ( cmp_1 == false ) && ( debugger->verbose_level >= VERBOSE_ALL || true ) ) 
      {
         DP_MATRIX_VIZ_Compare( cloud_MX_diff, edg_diag_lin, edg_diag_quad );
         DP_MATRIX_VIZ_Trace( cloud_MX_diff, tr );
         DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, stdout );
         FILE* fp = fopen( "test_output/my.quadvlin.diag.viz", "w+" );
         DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, fp );
         fclose(fp);
      } 
      /* test post-reorient */
      int cmp_2 = ( EDGEBOUNDS_Compare_by_Cloud( edg_row_lin, cloud_MX_lin, edg_row_quad, cloud_MX_quad ) == 0 );
      printf("POST-REORIENT?\t\t%s\n", ( cmp_2 == true ) ? "PASS" : "FAIL" );
      if ( ( cmp_2 == false ) && ( debugger->verbose_level >= VERBOSE_ALL || true ) ) 
      {
         DP_MATRIX_VIZ_Compare( cloud_MX_diff, edg_row_lin, edg_row_quad );
         DP_MATRIX_VIZ_Trace( cloud_MX_diff, tr );
         DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, stdout );
         FILE* fp = fopen( "test_output/my.quadvlin.diag.viz", "w+" );
         DP_MATRIX_VIZ_Color_Dump( cloud_MX_diff, fp );
         fclose(fp);
      } 
      if ( ( cmp_1 == false || cmp_2 == false ) && ( debugger->verbose_level >= VERBOSE_ALL || true ) )
      {
         printf("Edgebounds (quadratic, diag-wise):\n");
         EDGEBOUNDS_Dump( edg_diag_quad, stdout );
         printf("Edgebounds (linear, diag-wise):\n");
         EDGEBOUNDS_Dump( edg_diag_lin, stdout );
         printf("Edgebounds (quadratic, row-wise):\n");
         EDGEBOUNDS_Dump( edg_row_quad, stdout );
         printf("Edgebounds (linear, row-wise):\n");
         EDGEBOUNDS_Dump( edg_row_lin, stdout );
      }
      printf("=== MERGE & REORIENT CLOUD -> END ===\n\n");
   }

   /* build sparse matrix around edgebounds */
   {
      MATRIX_3D_SPARSE_Reuse( st_SMX );
      MATRIX_3D_SPARSE_Shape_Like_Edgebounds( st_SMX, edg_row_lin );
      // EDGEBOUNDS_Dump( st_SMX->edg_outer, stdout );
      EDGEBOUNDS_Save( st_SMX->edg_outer, "test_output/my.sparse.outer.edg" );
      EDGEBOUNDS_Save( st_SMX->edg_inner, "test_output/my.sparse.inner.edg" );
   }

   // /* stats */
   // num_cells = cloud_Cell_Count(Q, T, st_MX, sp_MX);
   // scores->perc_cells = (float)num_cells/(float)tot_cells;
   // printf("Perc. Total Cells Computed = %d/%d = %f\n", num_cells, tot_cells, scores->perc_cells);
   // scores->perc_window = (float)num_cells/(float)window_cells;
   // printf("Perc. Window Cells Computed = %d/%d = %f\n", num_cells, window_cells, scores->perc_window);
   // printf("\n");

   /* fill cloud for naive algs */
   MATRIX_2D_Fill( cloud_MX_naive, 0 );
   MATRIX_2D_Cloud_Fill( cloud_MX_naive, edg_row_lin, 1 );
   // DP_MATRIX_VIZ_Dump( cloud_MX_naive, stdout );
   // EDGEBOUNDS_Dump( edg_row_lin, stdout );

   /* need to clean matrix after fwd/bck */
   DP_MATRIX_Fill( Q, T, st_MX_quad, sp_MX_quad, -INF );
   // DP_MATRIX_Fill( Q, T, st_MX3_lin, sp_MX_lin, -INF );

   /* bounded forward */
   printf("=== BOUND FORWARD -> START ===\n");
   /* bounded forward (naive) */
   printf("==> bound forward (naive)\n");
   run_Bound_Forward_Naive(q_seq, t_prof, Q, T, st_MX_naive, sp_MX_naive, cloud_MX_naive, &sc);
   printf("Bounded Forward Score (naive):\t%f\n", sc);
   scores->naive_cloud_fwd = sc;
   DP_MATRIX_Trace_Save(Q, T, st_MX_naive, sp_MX_naive, tr, "test_output/my.bound_fwd.naive.mx");
   /* bounded forward (quadratic) */
   printf("==> bound forward (quad)\n");
   run_Bound_Forward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, edg_row_lin, &sc);
   printf("Bound Forward Score  (quad):\t%f\n", sc);
   scores->quad_cloud_fwd = sc;
   DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.bound_fwd.quad.mx");
   /* bounded forward (linear) */
   printf("==> bound forward (lin)\n");
   run_Bound_Forward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, edg_row_lin, &sc);
   printf("Bound Forward Score   (lin):\t%f\n", sc);
   scores->lin_cloud_fwd = sc;
   MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
   DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.bound_fwd.lin.mx");
   /* bounded forward (comparison) */
   printf("==> bound forward comparison: naive vs. quad vs. lin\n");
   cmp = ( cmp_tol_cust( scores->naive_cloud_fwd, scores->quad_cloud_fwd, 1e-3 ) && 
            cmp_tol_cust( scores->naive_cloud_fwd, scores->lin_cloud_fwd, 1e-3 ) );
   printf("Bound Forward:\tSCORES?\t\t%s\n", cmp ? "PASS" : "FAIL" );
   dp_cmp  = ( DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin ) == 0 );
   printf("Bound Forward:\tVALUES?\t\t%s\n", dp_cmp ? "PASS" : "FAIL" );
   if ( dp_cmp == false ) {
      DP_MATRIX_Diff( st_MX_lin, sp_MX_lin, st_MX_quad, sp_MX_quad, st_MX_diff, sp_MX_diff );
      DP_MATRIX_Dump( Q, T, st_MX_naive, sp_MX_naive, stdout );
      DP_MATRIX_Dump( Q, T, st_MX_quad, sp_MX_quad, stdout );
      DP_MATRIX_Dump( Q, T, st_MX_lin, sp_MX_lin, stdout );
   }
   printf("=== BOUNDED FORWARD -> END ===\n\n");

   /* need to clean matrix after fwd/bck */
   DP_MATRIX_Fill( Q, T, st_MX_quad, sp_MX_quad, -INF );
   // DP_MATRIX_Fill( Q, T, st_MX3_lin, sp_MX_lin, -INF );

   /* bounded backward */
   printf("=== BOUND BACKWARD -> START ===\n");
   /* bounded backward (naive) */
   printf("==> bound backward (naive)\n");
   run_Bound_Backward_Naive(q_seq, t_prof, Q, T, st_MX_naive, sp_MX_naive, cloud_MX_naive, &sc);
   printf("Bound Backward Score (naive):\t%f\n", sc);
   scores->naive_cloud_bck = sc;
   // MATRIX_3D_Copy( st_MX_naive, debugger->test_MX );
   DP_MATRIX_Trace_Save(Q, T, st_MX_naive, sp_MX_naive, tr, "test_output/my.bound_bck.naive.mx");
   /* bounded backward (quadratic) */
   printf("==> bound backward (quad)\n");
   run_Bound_Backward_Quad(q_seq, t_prof, Q, T, st_MX_quad, sp_MX_quad, edg_row_lin, &sc);
   printf("Bound Backward Score  (quad):\t%f\n", sc);
   scores->quad_cloud_fwd = sc;
   // MATRIX_3D_Copy( st_MX_quad, debugger->test_MX );
   DP_MATRIX_Trace_Save(Q, T, st_MX_quad, sp_MX_quad, tr, "test_output/my.bound_bck.quad.mx");
   /* bounded backward (linear) */
   printf("==> bound backward (lin)\n");
   run_Bound_Backward_Linear(q_seq, t_prof, Q, T, st_MX3_lin, sp_MX_lin, edg_row_lin, &sc);
   printf("Bound Backward Score   (lin):\t%f\n", sc);
   scores->lin_cloud_fwd = sc;
   MATRIX_3D_Copy( st_MX_lin, debugger->test_MX );
   DP_MATRIX_Trace_Save(Q, T, st_MX_lin, sp_MX_lin, tr, "test_output/my.bound_bck.lin.mx");
   /* bounded backward (comparison) */
   printf("==> bound backward comparison: naive vs. quad vs. lin\n");
   #if DEBUG
   {
      int sc_cmp = ( cmp_tol_cust( scores->naive_cloud_bck, scores->quad_cloud_bck, 1e-3 ) && 
                     cmp_tol_cust( scores->naive_cloud_bck, scores->lin_cloud_bck, 1e-3 ) );
      printf("Bound Backward:\tSCORES?\t\t%s\n", sc_cmp ? "PASS" : "FAIL" );
      int dp_cmp  = ( DP_MATRIX_Compare( st_MX_quad, sp_MX_quad, st_MX_lin, sp_MX_lin ) == 0 );
      printf("Bound Backward:\tVALUES?\t\t%s\n", dp_cmp ? "PASS" : "FAIL" );
      if ( ( dp_cmp == false ) && ( debugger->verbose_level >= VERBOSE_ALL ) ) {
         DP_MATRIX_Diff( st_MX_lin, sp_MX_lin, st_MX_quad, sp_MX_quad, st_MX_diff, sp_MX_diff );
         DP_MATRIX_Dump( Q, T, st_MX_naive, sp_MX_naive, stdout );
         DP_MATRIX_Dump( Q, T, st_MX_quad, sp_MX_quad, stdout );
         DP_MATRIX_Dump( Q, T, st_MX_lin, sp_MX_lin, stdout );
      }
   }
   #endif
   printf("=== BOUND BACKWARD -> END ===\n\n");

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