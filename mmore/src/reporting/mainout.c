/*******************************************************************************
 *  FILE:      mainout.c
 *  PURPOSE:   Reporting for standard output for main pipeline.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       - 
 *  NOTES:
 *    - Need to buffer output!
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

/* header */
#include "_reporting.h"
#include "mainout.h"

/* === STDOUT OUTPUT === */
/* EXAMPLE:
 *
   [header/]
   # hmmsearch :: search profile(s) against a sequence database
   # HMMER 3.3 (Nov 2019); http://hmmer.org/
   # Copyright (C) 2019 Howard Hughes Medical Institute.
   # Freely distributed under the BSD open source license.
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   # query HMM file:                  test-input/3-PAP.hmm
   # target sequence database:        test-input/3-PAP.fa
   # per-seq hits tabular output:     tblout.csv
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   [entry/]
   Query:       3-PAP  [M=133]
   Accession:   PF12578.1
   Description: Myotubularin-associated protein
   Scores for complete sequences (score includes all domains):
      --- full sequence ---   --- best 1 domain ---    -#dom-
       E-value  score  bias    E-value  score  bias    exp  N  Sequence                 Description
       ------- ------ -----    ------- ------ -----   ---- --  --------                 -----------
       3.2e-12   34.5   0.0    1.5e-08   22.6   0.0    2.5  2  3-PAP/16/510-647/718-827  domains: MTMRA_DANRE/548-685 C3Z9W9
       2.7e-11   31.5   0.0    3.9e-10   27.7   0.0    2.4  2  3-PAP/13/1-136/374-513    domains: MTMRB_MOUSE/553-688 B7Q8P0
       1.1e-07   19.9   0.0    2.1e-07   18.9   0.0    1.5  1  3-PAP/14/86-218/365-501   domains: MTMRC_PONAB/559-691 A4HUS9

   Domain annotation for each sequence (and alignments):
   >> 3-PAP/16/510-647/718-827  domains: MTMRA_DANRE/548-685 C3Z9W9_BRAFL/506-615
      #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
    ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
      1 !    9.1   0.0   0.00022   0.00022      68     106 ..     582     618 ..     557     637 .. 0.78
      2 !   22.6   0.0   1.5e-08   1.5e-08      75     110 ..     769     803 ..     719     824 .. 0.73

     Alignments for each domain:
     == domain 1  score: 9.1 bits;  conditional E-value: 0.00022
                        3-PAP  68 LePqcrilDlevWdqCYfRWlPvLeikgGGqpqvDLfnR 106
                                  L Pq     l+vW+  + RW+P  +i +GG   v  f++
     3-PAP/16/510-647/718-827 582 LLPQLLPSHLSVWKLYFLRWVPEAQIPHGGP--VTAFHK 618
                                  667777789********************95..344444 PP

     == domain 2  score: 22.6 bits;  conditional E-value: 1.5e-08
                        3-PAP  75 lDlevWdqCYfRWlPvLeikgGGqpqvDLfnRllLs 110
                                    l++W qCY RW+P     gGG p  + f+  lL 
     3-PAP/16/510-647/718-827 769 AGLKLWTQCYMRWIPWAHTVGGGPPS-EYFHQCLLV 803
                                  5799************9998888665.777777664 PP

      ...

      [footer/]

      Internal pipeline statistics summary:
      -------------------------------------
      Query model(s):                            1  (133 nodes)
      Target sequences:                          3  (2475 residues searched)
      Passed MSV filter:                         3  (1); expected 0.1 (0.02)
      Passed bias filter:                        3  (1); expected 0.1 (0.02)
      Passed Vit filter:                         3  (1); expected 0.0 (0.001)
      Passed Fwd filter:                         3  (1); expected 0.0 (1e-05)
      Initial search space (Z):                  3  [actual number of targets]
      Domain search space  (domZ):               3  [number of targets reported over threshold]
      # CPU time: 0.01u 0.00s 00:00:00.01 Elapsed: 00:00:00.01
      # Mc/sec: 24.80
      //
      [ok]
 *
 */

/* FUNCTION:   REPORT_stdout_header()
 * SYNOPSIS:   Print Header to Output <fp>.
 *             (modeled after HMMER, see example)
 */
void REPORT_stdout_header(    WORKER*     worker,
                              FILE*       fp )
{
   ARGS* args     = worker->args;
   int   left_pad = 30;

   /* header info */
   REPORT_horizontal_rule( fp );
   fprintf( fp, "# %s :: %s :: %s\n", 
      BUILD_PROGRAM,
      PIPELINES[args->pipeline_mode].name,
      BUILD_DESCRIPT );
   fprintf( fp, "# %s %s (%s): http://github.com/TravisWheelerLab/fb-pruner/\n",
      BUILD_PROGRAM, 
      BUILD_NAME, 
      BUILD_DATE );
   fprintf( fp, "# Copyright (C) 2020 Travis Wheeler Lab, University of Montana.\n" );
   /* input files */
   REPORT_horizontal_rule( fp );
   fprintf( fp, "# %*s: %s\n", 
      left_pad, "HMM file (Target)", args->t_filepath );
   fprintf( fp, "# %*s: %s\n", 
      left_pad, "SEQ file (Query)", args->q_filepath );
   fprintf( fp, "# %*s: %s\n", 
      left_pad, "Target Index", args->t_indexpath );
   fprintf( fp, "# %*s: %s\n", 
      left_pad, "Query Index", args->q_indexpath );

   if ( args->pipeline_mode == PIPELINE_MMSEQS ) {
      fprintf( fp, "# %*s: %s\n", 
         left_pad, "MMSEQS RESULT file", args->mmseqs_m8_filepath );
   }
   if ( args->pipeline_mode == PIPELINE_MAIN ) {

   }

   REPORT_horizontal_rule( fp );
   fprintf( fp, "\n");
}

/*    FUNCTION:   REPORT_stdout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_stdout_entry(  WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp )
{
   ARGS*          args           = worker->args;
   BUFFER*        buffer         = worker->buffer;
   ALIGNMENT*     aln            = worker->trace_vit;
   HMM_PROFILE*   t_prof         = worker->t_prof;
   SEQUENCE*      q_seq          = worker->q_seq;
   ALL_SCORES*    scores         = &result->scores;
   SCORES*        finalsc        = &result->final_scores;

   /* if running an alignment */
   bool           is_run_aln     = false;

   /* short names for alignment */
   STR            t_name         = "=T=";      /* short target name for alignment */
   STR            q_name         = "=Q=";      /* short query name for alignment */
   STR            c_name         = "=X=";      /* short center name for alignment */

   /* alignment strings */
   STR            cigar_aln      = NULL;
   STR            target_aln     = NULL;
   STR            query_aln      = NULL;
   STR            center_aln     = NULL;
   STR            state_aln      = NULL;

   /* some should be passed as arguments (reporter object?) */
   int            pad_left       = 3;
   int            ind_left       = 6;
   int            field_width    = 15;       
   int            name_width     = 10;       /* number of characters allowed in name field */
   /* number of residues per window in alignment window */
   int            def_width      = 100;      /* default window size */
   int            aln_width      = 100;      /* current window size (can change to fit alignment) */

   /* type of alignment, if any */
   if ( args->is_run_postaln == true ) {
      is_run_aln = true;
      aln = worker->trace_post;
   }
   elif ( args->is_run_vitaln == true ) {
      is_run_aln = true;
      aln = worker->trace_vit;
   } 
   elif ( args->is_run_mmseqsaln == true ) {
      is_run_aln = true;
      fprintf(stderr, "ERROR: MMSEQS Alignment not supported.\n");
      ERRORCHECK_exit(EXIT_FAILURE);
   }
   else {
      is_run_aln = false;
   }

   /* if alignemnt strings have not been produced yet, do it now */
   if ( aln->is_cigar_aln == false ) {
      ALIGNMENT_Build_MMSEQS_Style( worker->trace_vit, worker->q_seq, worker->t_prof );
   }
   cigar_aln = VECTOR_CHAR_GetArray( aln->cigar_aln );
   cigar_aln = ( STR_GetLength(cigar_aln) > 0 ? cigar_aln : "--" );

   if ( aln->is_hmmer_aln == false ) {
      ALIGNMENT_Build_HMMER_Style( worker->trace_vit, worker->q_seq, worker->t_prof );
   }
   target_aln  = VECTOR_CHAR_GetArray( aln->target_aln );
   query_aln   = VECTOR_CHAR_GetArray( aln->query_aln );
   center_aln  = VECTOR_CHAR_GetArray( aln->center_aln );
   state_aln   = VECTOR_CHAR_GetArray( aln->state_aln );

   /* Meta Data */
   REPORT_horizontal_rule( fp );
   fprintf( fp, "%*s %s [L=%d]\n", 
      -field_width, "Query:", t_prof->name, t_prof->N );
   fprintf( fp, "%*s %s [L=%d]\n", 
      -field_width, "Target:", q_seq->name, q_seq->N );
   fprintf( fp, "%*s %s\n", 
      -field_width, "Accession:", (t_prof->acc ? t_prof->acc : "--" ) );
   fprintf( fp, "%*s %s\n", 
      -field_width, "Description:", (t_prof->desc ? t_prof->desc : "--" ) );
   /* Scores Header */
   fprintf( fp, "== %*s\n", 
      0, "Scores for complete sequences:");
   
   /* Domain Table */
   int N_regions = ALIGNMENT_GetNumRegions( aln );
   for ( int i_domain = 0; i_domain < N_regions; i_domain++ ) 
   {
      /* Alignment Header */
      fprintf( fp, "==== %*s %d :: %s %3.2f %s | %s %-3.2e %s\n",
         0, "Domain",      i_domain,                        /* domain number */
         "Score:",         finalsc->seq_sc,     "bits",     /* bit-score */
         "E-value:",       finalsc->eval,       ""          /* e-value (raw) */
      );

      /* MMSEQS-style Cigar Alignment */
      if ( is_run_aln == true )
      {
         fprintf( fp, "%*s\n", 
            0, "==== MMSEQS Cigar Alignment ===:" );
         fprintf( fp, "%*.*s %s\n\n",
            name_width, name_width, "", cigar_aln );
      }

      /* HMMER-style Pairwise Alignment */
      if ( is_run_aln == true )
      {
         int      offset   = 0; 
         TRACE    tr_beg, tr_end;
         {
            printf( "target: %d, query: %d, center: %d\n", 
               VECTOR_CHAR_GetSize( aln->target_aln ), VECTOR_CHAR_GetSize( aln->query_aln ), VECTOR_CHAR_GetSize( aln->center_aln ) );

            fprintf( fp, "==== HMMER-Style Alignment === | Domain: %d\n",
               i_domain );
            /* create alignment rows */
            for ( int i = aln->beg; i <= aln->end; i += aln_width, offset += aln_width ) 
            {
               
               /* if remaining alignment exceeds window size, constrain it */
               aln_width = MIN( aln->end - aln->beg, def_width );

               fprintf( fp, "%*.*s %5d %.*s %-5d\n",
                  name_width, name_width,                         /* padding */
                  t_name,                                         /* name */
                  VEC_X( aln->traces, i ).t_0,                    /* starting index */
                  aln_width,                                      /* number of residues per line */
                  &VEC_X( aln->target_aln, offset ),              /* alignment residues */
                  VEC_X( aln->traces, i + aln_width - 1 ).t_0     /* ending index */
               );
               fprintf( fp, "%*.*s %5d %.*s %-5d\n",
                  name_width, name_width,                         /* padding */
                  c_name,                                         /* name */
                  i,                                              /* starting index */
                  aln_width,                                      /* number of residues per line */
                  &VEC_X( aln->center_aln, offset ),              /* alignment residues */
                  i + aln_width - 1                               /* ending index */
               );
               fprintf( fp, "%*.*s %5d %.*s %-5d\n",
                  name_width, name_width,                         /* padding */
                  q_name,                                         /* name */
                  VEC_X( aln->traces, i ).q_0,                    /* starting index */
                  aln_width,                                      /* number of residues per line */
                  &VEC_X( aln->query_aln, offset ),               /* alignment residues */
                  VEC_X( aln->traces, i + aln_width - 1 ).q_0     /* ending index */
               );
               fprintf( fp, "\n" );
            }
         }
      }
      
   }

   REPORT_horizontal_rule( fp );
   fprintf( fp, "\n");
}

/*!  FUNCTION:   REPORT_stdout_align()
 *   SYNOPSIS:   Print <i>th alignment.
 */
void 
REPORT_stdout_domtbl(   WORKER*  worker,
                        RESULT*  result,
                        int      i_domain,
                        FILE*    fp )
{
   
}

/*!  FUNCTION:   REPORT_stdout_align()
 *   SYNOPSIS:   Print <i>th alignment.
 */
void 
REPORT_stdout_hmmer_align(    WORKER*  worker,
                              int      i_domain,
                              FILE*    fp )
{

}

/*!  FUNCTION:   REPORT_stdout_footer()
 *   SYNOPSIS:   Print Summary Statistics after all searches completed.
 *                (modeled after HMMER, see example)
 */
void 
REPORT_stdout_footer(   WORKER*  worker,
                        FILE*    fp )
{
   REPORT_stdout_footer_search_summary( worker, fp );
   REPORT_stdout_footer_time_summary( worker, fp );

   /* success */
   fprintf( fp, "\n# [ok.]\n" );
}

/*!  FUNCTION:   REPORT_stdout_footer_search_summary()
 *   SYNOPSIS:   Print Summary Statistics after all searches completed.
 *                (modeled after HMMER, see example)
 */
void 
REPORT_stdout_footer_search_summary(   WORKER*  worker,
                                       FILE*    fp )
{
   STATS*   stats    = worker->stats;

   /* formatting settings */
   char        str[128];
   const int   left_pad    = -45;
   const int   center_pad  = 0;

   /* statistics summary */
   fprintf( fp, "\nInternal pipeline statistics summary:\n" );
   fprintf( fp, "----------------------------------------\n" );
   fprintf( fp, "%*s %*d   %s\n", 
      left_pad, "Target models:", 
      center_pad, stats->n_target_db, "targets" );
   fprintf( fp, "%*s %*d   %s\n", 
      left_pad, "Query sequences:", 
      center_pad, stats->n_query_db, "queries" );
   /* TODO: Add breakdown of number of queries passing mmseqs prefilter */
   fprintf( fp, "%*s %*d   %s\n", 
      left_pad, "Total Number searches:", 
      center_pad, stats->n_searches, "searches" );
   fprintf( fp, "%*s %*d   %s\n", 
      left_pad, "Number Passed MMSEQS Prefilter:", 
      center_pad, stats->n_passed_prefilter, "searches" );
   fprintf( fp, "%*s %*d   %s\n", 
      left_pad, "Number Passed MMSEQS Viterbi Filter:", 
      center_pad, stats->n_passed_viterbi, "searches" );
   fprintf( fp, "%*s %*d   %s\n", 
      left_pad, "Number Passed MMORE Cloud Filter:", 
      center_pad, stats->n_passed_cloud, "searches" );
   fprintf( fp, "%*s %*d   %s\n", 
      left_pad, "Number Passed MMORE Fwdback Filter:", 
      center_pad, stats->n_passed_fwdback, "searches" );
   fprintf( fp, "%*s %*d   %s\n", 
      left_pad, "Initial search space (Z):", 
      center_pad, stats->n_query_db, "[actual number of targets]" );
   fprintf( fp, "%*s %*d   %s\n", 
      left_pad, "Domain search space (Z):", 
      center_pad, stats->n_reported_domains, "[number of targets reported over threshold]" );
   fprintf( fp, "\n" );
}

/*!  FUNCTION:   REPORT_stdout_footer_time_summary()
 *   SYNOPSIS:   Print Runtime Summary.
 */
void 
REPORT_stdout_footer_time_summary(     WORKER*  worker,
                                       FILE*    fp )
{
   TIMES*   t_times  = worker->times_totals;

   /* formatting settings */
   char        str[128];
   const int   left_pad = -45;
   const int   center_pad = 0;

   /* runtime breakdown */
   fprintf( fp, "\nRuntime breakdown:\n" );
   fprintf( fp, "----------------------------------------\n" );
   fprintf( fp, "%*s %*.3f %s\n", 
      left_pad, "Total Runtime (secs):", 
      center_pad, t_times->program, "secs" );
}
