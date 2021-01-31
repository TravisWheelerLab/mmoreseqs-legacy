/*******************************************************************************
 *  FILE:      mainout.c
 *  PURPOSE:   Reporting for standard output for main pipeline.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       - 
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
         left_pad, "MMSEQS RESULT file", (args->mmseqs_res_filepath ? args->mmseqs_res_filepath : "--") );
   }
   if ( args->pipeline_mode == PIPELINE_MAIN ) {

   }

   REPORT_horizontal_rule( fp );
   printf("\n");
}

/*    FUNCTION:   REPORT_stdout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_stdout_entry(  WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp )
{
   HMM_PROFILE_Set_Consensus( worker->t_prof );

   ALIGNMENT*     aln            = worker->traceback;
   HMM_PROFILE*   t_prof         = worker->t_prof;
   SEQUENCE*      q_seq          = worker->q_seq;

   /* TODO: some should be passed as arguments (reporter object?) */
   int            pad_left       = 3;
   int            ind_left       = 6;
   int            field_width    = 15;
   int            name_width     = 10;
   int            aln_width      = 100;   /* row width of alignments */

   /* if alignemnt strings have not been produced yet, do it now */
   if ( aln->cigar_aln == NULL ) {
      ALIGNMENT_Build_MMSEQS_Style( worker->traceback, worker->q_seq, worker->t_prof );
   }
   if ( aln->center_aln == NULL ) {
      ALIGNMENT_Build_HMMER_Style( worker->traceback, worker->q_seq, worker->t_prof );
   }

   /* Meta Data */
   REPORT_horizontal_rule( fp );
   fprintf( fp, "%*s %s [L=%d]\n", 
      field_width, "Query:", t_prof->name, t_prof->N );
   fprintf( fp, "%*s %s [L=%d]\n", 
      field_width, "Target:", q_seq->name, q_seq->N );
   // fprintf( fp, "%*s %s\n", 
   //    field_width, "Accession:", (t_prof->acc ? t_prof->acc : "--" ) );
   // fprintf( fp, "%*s %s\n", 
   //    field_width, "Description:", (t_prof->desc ? t_prof->desc : "--" ) );
   /* Scores Header */
   fprintf( fp, "== %*s\n", 
      0, "Scores for complete sequences:");
   /* Alignment Header */
   fprintf( fp, "== %*s\n", 
      0, "Alignment:");
   fprintf( fp, "==== %*s %d :: %s %3.2f | %s %-3.2e | %s %-3.2e\n",
      0, "Domain", 1,                                    /* domain number */
      "Score (bits):", result->final_scores.seq_sc,         /* bit-score */
      "Raw E-value:", result->final_scores.eval,            /* e-value (raw) */
      "Cond. E-value:", result->final_scores.eval           /* e-value (conditional) */
   );
   fprintf( fp, "%*s %s\n", 
      0, "==== Cigar:", (aln->cigar_aln != NULL ? aln->cigar_aln : "--") );
   /* create alignment rows */
   int offset = 0; 
   for ( int i = aln->beg; i < aln->end; i += aln_width, offset += aln_width ) 
   {
      if ( aln->end - i < aln_width ) {
         aln_width = aln->end - i;
      }

      fprintf( fp, "%*.*s %5d %.*s %-5d\n",
         name_width, name_width,                   /* padding */
         t_prof->name,                             /* name */
         aln->traces->data[i].t_0,                 /* starting index */
         aln_width,                                /* number of residues per line */
         &aln->target_aln[offset],                 /* alignment residues */
         aln->traces->data[i + aln_width - 1].t_0  /* ending index */
      );
      fprintf( fp, "%*.*s %5d %.*s %-5d\n",
         name_width, name_width,                   /* padding */
         "",                                       /* name */
         i,                                        /* starting index */
         aln_width,                                /* number of residues per line */
         &aln->center_aln[offset],                 /* alignment residues */
         i + aln_width - 1                         /* ending index */
      );
      fprintf( fp, "%*.*s %5d %.*s %-5d\n",
         name_width, name_width,                   /* padding */
         q_seq->name,                              /* name */
         aln->traces->data[i].q_0,                 /* starting index */
         aln_width,                                /* number of residues per line */
         &aln->query_aln[offset],                  /* alignment residues */
         aln->traces->data[i + aln_width - 1].q_0  /* ending index */
      );
   }

   REPORT_horizontal_rule( fp );
   fprintf( fp, "\n");
}

/*    FUNCTION:   REPORT_stdout_footer()
 *    SYNOPSIS:   Print Summary Statistics after all searches completed.
 *                (modeled after HMMER, see example)
 */
void REPORT_stdout_footer(    WORKER*  worker,
                              FILE*    fp )
{
   char  str[50];
   int   left_pad = -35;
   int   center_pad = 0;

   /* statistics summary */
   fprintf( fp, "\nInternal pipeline statistics summary:\n" );
   fprintf( fp, "----------------------------------------\n" );
   fprintf( fp, "%*s %*d %s\n", 
      left_pad, "Target models:", center_pad, worker->t_index->N, "" );
   fprintf( fp, "%*s %*d %s\n", 
      left_pad, "Query sequences:", center_pad, worker->q_index->N, "" );
   /* TODO: Add breakdown of number of queries passing mmseqs prefilter */
   fprintf( fp, "%*s %*d %s\n", 
      left_pad, "Number searches:", center_pad, worker->num_searches, "" );
   fprintf( fp, "%*s %*d %s\n", 
      left_pad, "Passed FB-Pruner Filter:", center_pad, 0, "" );
   fprintf( fp, "%*s %*d %s\n", 
      left_pad, "Initial search space (Z):", center_pad, worker->q_index->N, "" );
   fprintf( fp, "%*s %*d %s\n", 
      left_pad, "Domain search space (Z):", center_pad, 0, "" );
   fprintf( fp, "\n" );
   /* runtime breakdown */
   fprintf( fp, "\nRuntime breakdown:\n" );
   fprintf( fp, "----------------------------------------\n" );
   fprintf( fp, "%*s %*.3f %s\n", 
      left_pad, "Total Runtime (secs):", center_pad, 0.0f, "secs" );
   fprintf( fp, "%*s %*.3f %s\n", 
      left_pad, "Sum of Runtime (secs):", center_pad, 0.0f, "secs" );
   /* success */
   fprintf( fp, "\n# [ok.]\n" );
}
