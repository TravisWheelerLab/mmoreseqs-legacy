/*******************************************************************************
 *  FILE:      report.c
 *  PURPOSE:   Reporting Subroutines for generating output.
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
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../parsers/parsers.h"
#include "../algs_linear/algs_linear.h"
#include "../algs_quad/algs_quad.h"
#include "../algs_naive/algs_naive.h"

/* header */
#include "report.h"

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
      PIPELINE_NAMES[args->pipeline_mode],
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
   worker->runtime = CLOCK_Get_Total_Runtime( worker->clok );

   // /* test that clock works */
   // double test_beg = CLOCK_Get_RealTime();
   // my_delay(3000); 
   // double test_end = CLOCK_Get_RealTime();

   char str[50];
   int left_pad = -35;
   int center_pad = 0;

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
      left_pad, "Total Runtime (secs):", center_pad, worker->runtime * 1000, "secs" );
   fprintf( fp, "%*s %*.3f %s\n", 
      left_pad, "Sum of Runtime (secs):", center_pad, 0.0f, "secs" );
   /* success */
   fprintf( fp, "\n# [ok.]\n" );
}

/* === TBLOUT OUTPUT === */
/* EXAMPLE:
 *
      [header/]
      #                                                                   --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
      # target name            accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
      #    ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
      [entry/]
      3-PAP/16/510-647/718-827 -          3-PAP                PF12578.1    3.2e-12   34.5   0.0   1.5e-08   22.6   0.0   2.5   2   0   0   2   2   2   2 domains: MTMRA_DANRE/548-685 C3Z9W9_BRAFL/506-615
      3-PAP/13/1-136/374-513   -          3-PAP                PF12578.1    2.7e-11   31.5   0.0   3.9e-10   27.7   0.0   2.4   2   0   0   2   2   2   1 domains: MTMRB_MOUSE/553-688 B7Q8P0_IXOSC/504-643
      3-PAP/14/86-218/365-501  -          3-PAP                PF12578.1    1.1e-07   19.9   0.0   2.1e-07   18.9   0.0   1.5   1   0   0   1   1   1   1 domains: MTMRC_PONAB/559-691 A4HUS9_LEIIN/8-144
      [footer/]
      #
      # Program:         hmmsearch
      # Version:         3.3 (Nov 2019)
      # Pipeline mode:   SEARCH
      # Query file:      test-input/3-PAP.hmm
      # Target file:     test-input/3-PAP.fa
      # Option settings: hmmsearch --tblout tblout.csv test-input/3-PAP.hmm test-input/3-PAP.fa 
      # Current dir:     /home/devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/fb-pruner
      # Date:            Fri Aug 28 00:42:33 2020
      # [ok]
 *
 */

/*    FUNCTION:   REPORT_tblout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_tblout_header( WORKER*  worker,
                           FILE*    fp )
{
   int qnamew = 20;
   int tnamew = 20;
   int qaccw  = 10;
   int taccw  = 10;
   int posw   = 7;

   fprintf( fp, "#%*s %22s %22s %33s\n",
      tnamew + qnamew + taccw + qaccw + 2, 
      "", 
      "--- full sequence ----", 
      "--- best 1 domain ----", 
      "--- domain number estimation ----" 
   );

   fprintf( fp, "#%-*s %-*s %-*s %-*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
      tnamew - 1, 
      " target name",        
      taccw, 
      "accession",  
      qnamew, 
      "query name",           
      qaccw, 
      "accession",  
      "  E-value", 
      " score", 
      " bias", 
      "  E-value", 
      " score", 
      " bias", 
      "exp", 
      "reg", 
      "clu", 
      " ov", 
      "env", 
      "dom", 
      "rep", 
      "inc", 
      "description of target" 
   );

   fprintf( fp, "#%*s %*s %*s %*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
      tnamew - 1, 
      "-------------------", 
      taccw, 
      "----------", 
      qnamew, 
      "--------------------", 
      qaccw, 
      "----------", 
      "---------", 
      "------", 
      "-----", 
      "---------", 
      "------", 
      "-----", 
      "---", 
      "---", 
      "---", 
      "---", 
      "---", 
      "---", 
      "---", 
      "---", 
      "---------------------" 
   );
}

/*    FUNCTION:   REPORT_tblout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_tblout_entry(  WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp )
{
   int qnamew = 20;
   int tnamew = 20;
   int qaccw  = 10;
   int taccw  = 10;
   int posw   = 7;

   /* computed output */
   float bias_correction = result->final_scores.pre_sc - result->final_scores.seq_sc;

   fprintf( fp, "%-*s %-*s %-*s %-*s %9.2g %6.1f %5.1f %9.2g %6.1f %5.1f %5.1f %3d %3d %3d %3d %3d %3d %3d %s\n", 
      /* query / target data */
      qnamew, worker->q_seq->name,                          /* query name */
      qaccw,  (NULL ? worker->q_seq->acc : "-"),            /* query accession */
      tnamew, worker->t_prof->name,                         /* target name */
      taccw,  (NULL ? worker->t_prof->acc : "-"),           /* target accession */
      /* full sequence */
      result->final_scores.eval,                            /* evalue */
      result->final_scores.nat_sc,                          /* score */
      bias_correction,                                      /* bias correction (pre_score - score) */
      /* best domain */
      0.0,                                                  /* evalue */
      0.0,                                                  /* score */
      0.0,                                                  /* bias correction */
      /* domain number estimation */
      0.0,                                                  /* number expected */
      0,                                                    /* number regions */
      0,                                                    /* number clustered */
      0,                                                    /* number overlapped */
      0,                                                    /* number envelopes */
      0,                                                    /* number domains */
      0,                                                    /* number reported */
      0,                                                    /* number included */
      /* target description */
      (NULL ? worker->t_prof->desc : "--")                  /* query description */
   );
}

/*    FUNCTION:   REPORT_tblout_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
void REPORT_tblout_footer(    WORKER*  worker,
                              FILE*    fp )
{
   ARGS* args = worker->args;

   int left_pad = 20;

   /* get current directory */
   char cwd[256];
   getcwd(cwd, sizeof(cwd));
   /* get current date/time */
   char* time_str = NULL;
   time_str = CLOCK_Get_DateTimeString( NULL );

   fprintf( fp, "# \n");
   fprintf( fp, "# %*s %s\n",       left_pad,   "Program:",          BUILD_PROGRAM);
   fprintf( fp, "# %*s %s (%s)\n",  left_pad,   "Version:",          BUILD_VERSION, BUILD_DATE );
   fprintf( fp, "# %*s %s\n",       left_pad,   "Pipeline mode:",    PIPELINE_NAMES[args->pipeline_mode] );
   fprintf( fp, "# %*s %s\n",       left_pad,   "Query file:",       args->q_filepath );
   fprintf( fp, "# %*s %s\n",       left_pad,   "Target file:",      args->t_filepath );
   fprintf( fp, "# %*s %s\n",       left_pad,   "Option settings:",  "" );
   if (cwd)       fprintf( fp, "# %*s %s\n",       left_pad,   "Current dir:",      cwd );
   // if (time_str)  fprintf( fp, "# %*s %s\n",       left_pad,   "Date:",             time_str );
   fprintf( fp, "# %*s\n",          left_pad,   "[ok]" );
}

/* === M8 OUTPUT === */
/* BLAST-style output 
   Description: The file is formatted as a tab-separated list with 12 columns: 
   - (1,2) identifiers query and target sequences/profiles, 
   - (3) sequence identity, 
   - (4) alignment leng
   - (5) number of mismatches, 
   - (6) number of gap openings, 
   - (7,8) domain (start,end)-position in query and in target, 
   - (9,10) domain
   - (11) E-value, and 
   - (12) bit score.
 */

/* EXAMPLE:
 *
   [header/]
   #  qid         tid                                  %id    aln   mis  gap  qb   qe   tb    te    eval       bit 
   #  ----------  -------------------------------      ---    ---   ---  ---  ---  ---  ---   ---   ----       ---
   [entry/]
      Aldolase_II  Aldolase_II/544/29-203/211-382      0.252  151   105  0    21   171  51    191   2.631E-07  54
      Aldolase_II  decoy94676                          0.230  91    69   0    66   156  164   253   2.572E+01  29
      Aldolase_II  decoy95461                          0.553  24    10   0    16   39   102   124   4.583E+01  29
   [footer/]
   # [ok]
 *
 */

/*    FUNCTION:   REPORT_m8out_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_m8out_header(  WORKER*  worker,
                           FILE*    fp )
{
   fprintf( fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
      "query",
      "target",
      "pident",
      "alnlen",
      "mismatch",
      "gapopen",
      "qstart",
      "qend",
      "tstart",
      "tend",
      "eval",
      "bits",
      "cigar"
   );
   fprintf( fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
      "-----",
      "------",
      "------",
      "------",
      "--------",
      "------",
      "------",
      "-----",
      "-------",
      "----",
      "----",
      "----",
      "-----"
   );
}

/*    FUNCTION:   REPORT_m8out_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_m8out_entry(   WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp )
{
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      q_seq    = worker->q_seq;
   ALIGNMENT*     aln      = worker->traceback;

   int aln_len = aln->end - aln->beg + 1;
   TRACE* beg  = &aln->traces->data[aln->beg];
   TRACE* end  = &aln->traces->data[aln->end];

   fprintf( fp, "%s\t%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%9.2g\t%6.1f\t%s\n", 
      t_prof->name,                    /* target name */
      q_seq->name,                     /* query name */
      aln->perc_id,                    /* percent id (number matches) */
      aln_len,                         /* alignment length */
      aln->num_misses,                 /* number of mismatches */    
      aln->num_gaps,                   /* number of gap openings */
      beg->t_0,                        /* query start */
      end->t_0,                        /* query end */
      beg->q_0,                        /* target start */
      end->q_0,                        /* target end */
      result->final_scores.eval,       /* evalue */
      result->final_scores.nat_sc,     /* bits */
      ( aln->cigar_aln != NULL ? aln->cigar_aln : "--" )
   );
}

/*    FUNCTION:   REPORT_m8out_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
void REPORT_m8out_footer(    WORKER*  worker,
                              FILE*    fp )
{
   fprintf( fp, "# [ok] [m8out]\n" );
}

/* === MYOUT FUNCTIONS === */

/* === MY OUTPUT === */
/* Custom-style output 
   Description: The file is formatted as a tab-separated list with 12 columns: 
   - (1) query profile identifier
   - (2) target sequences identifier 
   - (3) e-value, 
   - (4) bitscore before composition bias correction
   - (5) composition bias score
   - (6) bitscore after composition bias correction
   - (7) mmseqs viterbi bitscore
   - (8) total DP matrix cells used by full viterbi/forward-backward
   - (9) DP matrix cloud of cells used by MMORE
   - (10) start-end range of query in cloud
   - (11) start-end range of target in cloud
   - (12) time to run given search
 */

/* EXAMPLE:
 *
   [header/]
   #  qid          tid                                 %id    aln   mis  gap  qb   qe   tb    te    eval       bit 
   #  ---          ---                                 ---    ---   ---  ---  ---  ---  ---   ---   ----       ---
   [entry/]
      Aldolase_II  Aldolase_II/544/29-203/211-382      0.252  151   105  0    21   171  51    191   2.631E-07  54
      Aldolase_II  decoy94676                          0.230  91    69   0    66   156  164   253   2.572E+01  29
      Aldolase_II  decoy95461                          0.553  24    10   0    16   39   102   124   4.583E+01  29
   [footer/]
   # [ok]
 *
 */

/*    FUNCTION:   REPORT_myout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_myout_header(  WORKER*  worker,
                           FILE*    fp )
{
   fprintf( fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
      "query",
      "target",
      "evalue",
      "pre-sc",
      "comp-bias",
      "seq-sc",
      "mmseqs-vit-sc",
      "total-cells",
      "MMORE-cells",
      "perc-cells",
      "q-range",
      "t-range",
      "time"
   );
   fprintf( fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
      "------",
      "------",
      "------",
      "------",
      "--------",
      "------",
      "-------------",
      "----------",
      "----------",
      "----------",
      "-------",
      "-------",
      "----"
   );
}

/*    FUNCTION:   REPORT_myout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_myout_entry(   WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp )
{
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      q_seq    = worker->q_seq;
   ALIGNMENT*     aln      = worker->traceback;

   int aln_len = aln->end - aln->beg + 1;
   TRACE* beg  = &aln->traces->data[aln->beg];
   TRACE* end  = &aln->traces->data[aln->end];

   float percent_cells = (float)result->cloud_cells / (float)result->total_cells;

   fprintf( fp, "%s\t%s\t%.3g\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.5f\t%d-%d\t%d-%d\t%.5f\n", 
      t_prof->name,                    /* target name */
      q_seq->name,                     /* query name */
      result->final_scores.eval,       /* evalue */
      result->final_scores.pre_sc,     /* seq scores before correction */
      result->final_scores.seq_bias,   /* seq bias */
      result->final_scores.seq_sc,     /* seq score after correction */
      result->vit_natsc,               /* viterbi score (in mmore, this comes from mmseqs) */
      result->total_cells,             /* total number of cells computed by full viterbi */
      result->cloud_cells,             /* total number of cells computed by mmore */
      percent_cells,                   /* percent of total cells computed by mmore */
      result->target_start,            /* start of target range */
      result->target_end,              /* end of target range */
      result->query_start,             /* start of query range */
      result->query_end,               /* end of query range */
      result->time                     /* time for entire iteration */
   );
}

/*    FUNCTION:   REPORT_myout_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
void REPORT_myout_footer(    WORKER*  worker,
                              FILE*    fp )
{
   fprintf( fp, "# [ok] [myout]\n" );
}

/* === MYDOMOUT FUNCTIONS === */

/* === MYDOM OUTPUT === */
/* Custom-style output 
   Description: The file is formatted as a tab-separated list with 12 columns: 
   - (1) query profile identifier
   - (2) target sequences identifier 
   - (3) e-value, 
   - (4) bitscore before composition bias correction
   - (5) composition bias score
   - (6) bitscore after composition bias correction
   - (7) mmseqs viterbi bitscore
   - (8) total DP matrix cells used by full viterbi/forward-backward
   - (9) DP matrix cloud of cells used by MMORE
   - (10) start-end range of query in cloud
   - (11) start-end range of target in cloud
   - (12) time to run given search
 */

/* EXAMPLE:
 *
   [header/]
   #  qid          tid                                 %id    aln   mis  gap  qb   qe   tb    te    eval       bit 
   #  ---          ---                                 ---    ---   ---  ---  ---  ---  ---   ---   ----       ---
   [entry/]
      Aldolase_II  Aldolase_II/544/29-203/211-382      0.252  151   105  0    21   171  51    191   2.631E-07  54
      Aldolase_II  decoy94676                          0.230  91    69   0    66   156  164   253   2.572E+01  29
      Aldolase_II  decoy95461                          0.553  24    10   0    16   39   102   124   4.583E+01  29
   [footer/]
   # [ok]
 *
 */

/*    FUNCTION:   REPORT_mydomout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_mydomout_header(  WORKER*  worker,
                           FILE*    fp )
{
   fprintf( fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
      "query",
      "target",
      "evalue",
      "pre-sc",
      "comp-bias",
      "seq-sc",
      "mmseqs-vit-sc",
      "total-cells",
      "MMORE-cells",
      "perc-cells",
      "q-range",
      "t-range",
      "time"
   );
   fprintf( fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
      "------",
      "------",
      "------",
      "------",
      "--------",
      "------",
      "-------------",
      "----------",
      "----------",
      "----------",
      "-------",
      "-------",
      "----"
   );
}

/*    FUNCTION:   REPORT_mydomout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_mydomout_entry(   WORKER*  worker,
                              RESULT*  result,
                              FILE*    fp )
{
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      q_seq    = worker->q_seq;
   ALIGNMENT*     aln      = worker->traceback;
   DOMAIN_DEF*    dom_def  = worker->dom_def;

   int      best_idx    = dom_def->best;
   RANGE    dom_rng     = VEC_X( dom_def->dom_ranges, best_idx );
   float    dom_fwdsc   = VEC_X( dom_def->dom_fwdsc, best_idx );
   float    dom_bias    = VEC_X( dom_def->dom_bias, best_idx);
   float    dom_sc      = dom_def->best_sc;

   float percent_cells = (float)result->cloud_cells / (float)result->total_cells;

   fprintf( fp, "%s\t%s\t%.3g\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.5f\t%d-%d\t%d-%d\t%.5f\n", 
      t_prof->name,                    /* target name */
      q_seq->name,                     /* query name */
      result->final_scores.eval,       /* evalue */
      dom_fwdsc,     /* seq scores before correction */
      dom_bias,   /* seq bias */
      dom_sc,     /* seq score after correction */
      result->vit_natsc,               /* viterbi score (in mmore, this comes from mmseqs) */
      result->total_cells,             /* total number of cells computed by full viterbi */
      result->cloud_cells,             /* total number of cells computed by mmore */
      percent_cells,                   /* percent of total cells computed by mmore */
      result->target_start,            /* start of target range */
      result->target_end,              /* end of target range */
      dom_rng.beg,                     /* start of query range */
      dom_rng.end,                     /* end of query range */
      result->time                     /* time for entire iteration */
   );
}

/*    FUNCTION:   REPORT_mydomout_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
void REPORT_mydomout_footer(    WORKER*  worker,
                              FILE*    fp )
{
   fprintf( fp, "# [ok] [myout]\n" );
}


/* === UTILTITY FUNCTIONS === */

/*    FUNCTION:   REPORT_horizontal_rule()
 *    SYNOPSIS:   Print a horizontal rule. 
 */
inline
void REPORT_horizontal_rule( FILE* fp )
{
   fprintf( fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n" );
}

/*    FUNCTION:   REPORT_horizontal_rule()
 *    SYNOPSIS:   Print a horizontal rule of specified length.
 */
void REPORT_hr_size( FILE* fp, 
                     int   length )
{
   char* base_hr = "-";

   fprintf( fp, "#");
   for (int i = 0; i < length; i++) {
      fprintf( fp, "%s", base_hr);
   }
   fprintf( fp, "\n");
}