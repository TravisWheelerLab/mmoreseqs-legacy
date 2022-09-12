/*******************************************************************************
 *  - FILE:  hmmerout.c
 *  - DESC:  Reporting HMMER-style output.
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

/* === SD2DOUT OUTPUT === */
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
       3.2e-12   34.5   0.0    1.5e-08   22.6   0.0    2.5  2  3-PAP/16/510-647/718-827  domains: MD2MRA_DANRE/548-685 C3Z9W9
       2.7e-11   31.5   0.0    3.9e-10   27.7   0.0    2.4  2  3-PAP/13/1-136/374-513    domains: MD2MRB_MOUSE/553-688 B7Q8P0
       1.1e-07   19.9   0.0    2.1e-07   18.9   0.0    1.5  1  3-PAP/14/86-218/365-501   domains: MD2MRC_PONAB/559-691 A4HUS9

   Domain annotation for each sequence (and alignments):
   >> 3-PAP/16/510-647/718-827  domains: MD2MRA_DANRE/548-685 C3Z9W9_BRAFL/506-615
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

/*! FUNCTION:   REPORT_hmmerout_header()
 *  SYNOPSIS:   Print Header to Output <fp>.
 *              (modeled after HMMER output, see example)
 */
void REPORT_hmmerout_header(WORKER* worker, FILE* fp) {
  ARGS* args = worker->args;
  int left_pad = 30;

  /* header info */
  REPORT_horizontal_rule(fp);
  fprintf(fp, "# %s :: %s :: %s\n",
          BUILD_PROGRAM,
          PIPELINES[args->pipeline_mode].name,
          BUILD_DESC);
  fprintf(fp, "# %s (%s): %s\n",
          BUILD_PROGRAM,
          BUILD_DATE,
          BUILD_REPO);
  fprintf(fp, "# %s.\n",
          BUILD_COPYRIGHT);
  /* input files */
  REPORT_horizontal_rule(fp);
  fprintf(fp, "# %*s: %s\n",
          left_pad, "HMM file (Target)", args->t_filein);
  fprintf(fp, "# %*s: %s\n",
          left_pad, "SEQ file (Query)", args->q_filein);
  fprintf(fp, "# %*s: %s\n",
          left_pad, "Target Index", args->t_index_filein);
  fprintf(fp, "# %*s: %s\n",
          left_pad, "Query Index", args->q_index_filein);

  if (args->adjust_mmseqs_alns) {
    fprintf(fp, "# %*s: %s\n",
            left_pad, "MMSEQS RESULT file", args->mmseqs_m8_filein);
  } else {
  }

  REPORT_horizontal_rule(fp);
  fprintf(fp, "\n");
}

/*! FUNCTION:   REPORT_hmmerout_entry()
 *  SYNOPSIS:   Print all alignment data for current search
 *              (modeled after HMMER, see example)
 */
void REPORT_hmmerout_entry(WORKER* worker, RESULT* result, FILE* fp) {
  TASKS* tasks = worker->tasks;
  ARGS* args = worker->args;
  BUFFER* buffer = worker->buffer;
  ALIGNMENT* aln;
  ALIGNMENT* aln_vit = worker->trace_vit;
  ALIGNMENT* aln_post = worker->trace_post;
  HMM_PROFILE* t_prof = worker->t_prof;
  SEQUENCE* q_seq = worker->q_seq;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;

  /* if running an alignment */
  bool is_run_aln = false;

  /* short names for alignment */
  STR t_name = "=T=";   /* short target name for alignment */
  STR q_name = "=Q=";   /* short query name for alignment */
  STR c_name = "=X=";   /* residue name for alignment */
  STR st_name = "===";  /* state name for alignment */

  /* alignment strings */
  STR cigar_aln = NULL;
  STR target_aln = NULL;
  STR query_aln = NULL;
  STR center_aln = NULL;
  STR state_aln = NULL;

  /* some should be passed as arguments (reporter object?) */
  int pad_left = 3;
  int ind_left = 6;
  int field_width = 15;
  int name_width = 10; /* number of characters allowed in name field */
  /* number of residues per window in alignment window */
  int def_width = 100; /* default window size */
  int aln_width = 100; /* current window size (can change to fit alignment) */

  /* type of alignment, if any */
  if (args->is_run_postaln == true) {
    is_run_aln = true;
    aln = worker->trace_post;
  }
  elif (args->is_run_vitaln == true) {
    is_run_aln = true;
    aln = worker->trace_vit;
  }
  elif (args->is_run_mmseqsaln == true) {
    is_run_aln = true;
    ERRORCHECK_exit(EXIT_FAILURE);
  }
  else {
    is_run_aln = false;
  }

  if (is_run_aln == true) {
    /* if alignemnt strings have not been produced yet, do it now */
    if (aln->is_cigar_aln == false) {
      ALIGNMENT_Build_MMSEQS_Style(aln, worker->q_seq, worker->t_prof);
    }
    cigar_aln = VECTOR_CHAR_GetArray(aln->cigar_aln);
    cigar_aln = (STR_GetLength(cigar_aln) > 0 ? cigar_aln : "--");

    if (aln->is_hmmer_aln == false) {
      ALIGNMENT_Build_HMMER_Style(aln, worker->q_seq, worker->t_prof);
    }
    target_aln = VECTOR_CHAR_GetArray(aln->target_aln);
    query_aln = VECTOR_CHAR_GetArray(aln->query_aln);
    center_aln = VECTOR_CHAR_GetArray(aln->center_aln);
    state_aln = VECTOR_CHAR_GetArray(aln->state_aln);
  }

  /* Meta Data */
  REPORT_horizontal_rule(fp);
  fprintf(fp, "%*s %s [L=%d]\n",
          -field_width, "Query:", t_prof->name, t_prof->N);
  fprintf(fp, "%*s %s [L=%d]\n",
          -field_width, "Target:", q_seq->name, q_seq->N);
  fprintf(fp, "%*s %s\n",
          -field_width, "Accession:", (t_prof->acc ? t_prof->acc : "--"));
  fprintf(fp, "%*s %s\n",
          -field_width, "Description:", (t_prof->desc ? t_prof->desc : "--"));
  /* Scores Header */
  fprintf(fp, "== %*s\n",
          0, "Scores for complete sequences:");

  /* Domain Table */
  if (args->is_run_domains == true && is_run_aln == true) {
    int N_regions = ALIGNMENT_GetNumRegions(aln);
    fprintf(fp, "NUMBER DOMAINS: %d\n", N_regions);
    for (int i_domain = 0; i_domain < N_regions; i_domain++) {
      /* Alignment Header */
      fprintf(fp, "==== %*s %d :: %s %3.2f %s | %s %-3.2e %s | %s %d\n",
              0, "Domain", i_domain,             /* domain number */
              "Score:", finalsc->seq_sc, "bits", /* bit-score */
              "E-value:", finalsc->eval, "",     /* e-value (raw) */
              "Length:", aln->aln_len - 2
      );

      /* MMSEQS-style Cigar Alignment */
      if (is_run_aln == true) {
        fprintf(fp, "%*s\n",
                0, "==== MMSEQS Cigar Alignment ===");
        fprintf(fp, "%*.*s %s\n",
                name_width, name_width, "", cigar_aln);
      }

      /* HMMER-style Pairwise Alignment */
      if (is_run_aln == true) {
        TRACE tr_beg, tr_end;
        {
          fprintf(fp, "==== HMMER-Style Alignment ===\n");
          /* create alignment rows */
          int offset_length = VECTOR_CHAR_GetSize(aln->state_aln) - 1;
          int aln_length = aln->end - aln->beg - 1;
          int beg_idx, end_idx;
          int beg_offset, end_offset;
          end_idx = 1;
          for (int offset = 1; offset <= offset_length; offset += def_width) {
            /* if remaining alignment exceeds window size, constrain it */
            beg_offset = offset;
            end_offset = MIN(beg_offset + def_width - 1, offset_length - 1);
            aln_width = end_offset - beg_offset - 1;
            beg_idx = aln->beg + beg_offset;
            end_idx = aln->beg + end_offset - 1;
            if (aln_width < 1) break;

            /* query */
            fprintf(fp, "%*.*s %5d %.*s %-5d\n",
                    name_width, name_width,                   /* padding */
                    t_name,                                   /* name */
                    VEC_X(aln->traces, beg_idx).t_0,          /* starting index */
                    aln_width,                                /* number of residues per line */
                    &VEC_X(aln->query_aln, beg_offset),       /* alignment residues */
                    VEC_X(aln->traces, end_idx - 1).t_0       /* ending index */
            );
            /* alignment */
            fprintf(fp, "%*.*s %5d %.*s %-5d\n",
                    name_width, name_width,                   /* padding */
                    c_name,                                   /* name */
                    beg_offset,                               /* starting index */
                    aln_width,                                /* number of residues per line */
                    &VEC_X(aln->center_aln, beg_offset),      /* alignment residues */
                    end_offset                                /* ending index */
            );
            /* target */
            fprintf(fp, "%*.*s %5d %.*s %-5d\n",
                    name_width, name_width,                   /* padding */
                    q_name,                                   /* name */
                    VEC_X(aln->traces, beg_idx).q_0,          /* starting index */
                    aln_width,                                /* number of residues per line */
                    &VEC_X(aln->target_aln, beg_offset),      /* alignment residues */
                    VEC_X(aln->traces, end_idx).q_0           /* ending index */
            );
            // /* state */
            // fprintf(fp, "%*.*s %5d %.*s %-5d\n",
            //         name_width, name_width,                   /* padding */
            //         st_name,                                  /* name */
            //         beg_offset,                               /* starting index */
            //         aln_width,                                /* number of residues per line */
            //         &VEC_X(aln->state_aln, beg_offset),       /* alignment residues */
            //         end_offset                                /* ending index */
            // );
            fprintf(fp, "\n");

            if (aln_width < def_width - 2) break;
          }
        }
      }
    }
  }

  REPORT_horizontal_rule(fp);
  fprintf(fp, "\n");
}

/*!   FUNCTION:   REPORT_hmmerout_footer()
 *    SYNOPSIS:   Print Summary Statistics after all searches completed.
 *                (modeled after HMMER, see example)
 */
void REPORT_hmmerout_footer(WORKER* worker, FILE* fp) {
  REPORT_hmmerout_footer_search_summary(worker, fp);
  REPORT_hmmerout_footer_time_summary(worker, fp);

  /* success */
  fprintf(fp, "\n# [ok.]\n");
}

/*!   FUNCTION:   REPORT_hmmerout_footer_search_summary()
 *    SYNOPSIS:   Print Summary Statistics after all searches completed.
 *                (modeled after HMMER, see example)
 */
void REPORT_hmmerout_footer_search_summary(WORKER* worker, FILE* fp) {
  STATS* stats = worker->stats;

  /* formatting settings */
  char str[128];
  const int left_pad = -45;
  const int center_pad = 0;
  int n_searches;
  n_searches = stats->n_target_db * stats->n_query_db;

  /* statistics summary */
  fprintf(fp, "\nInternal pipeline statistics summary:\n");
  fprintf(fp, "----------------------------------------\n");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Target models:",
          center_pad, stats->n_target_db, "targets");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Query sequences:",
          center_pad, stats->n_query_db, "queries");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Total Number searches:",
          center_pad, n_searches, "searches");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Number Passed MMSEQS Prefilter:",
          center_pad, stats->n_passed_prefilter, "searches");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Number Passed MMSEQS Viterbi Filter:",
          center_pad, stats->n_passed_viterbi, "searches");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Number Passed MMORE Viterbi Filter:",
          center_pad, stats->n_passed_viterbi, "searches");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Number Passed MMORE Cloud Filter:",
          center_pad, stats->n_passed_cloud, "searches");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Number Passed MMORE Fwdback Filter:",
          center_pad, stats->n_passed_fwdback, "searches");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Number Passed MMORE Reporting Filter:",
          center_pad, stats->n_passed_report, "searches");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Initial search space (Z):",
          center_pad, stats->n_query_db, "[actual number of targets]");
  fprintf(fp, "%*s %*d   %s\n",
          left_pad, "Domain search space (Z):",
          center_pad, stats->n_reported_domains, "[number of targets reported over threshold]");
  fprintf(fp, "\n");
}

/*!  FUNCTION:   REPORT_stdout_footer_time_summary()
 *   SYNOPSIS:   Print Runtime Summary.
 */
void REPORT_hmmerout_footer_time_summary(WORKER* worker,
                                         FILE* fp) {
  TIMES* t_times = worker->times_totals;

  /* formatting settings */
  char str[128];
  const int left_pad = -45;
  const int center_pad = 0;

  /* runtime breakdown */
  fprintf(fp, "\nRuntime breakdown:\n");
  fprintf(fp, "----------------------------------------\n");
  fprintf(fp, "%*s %*.3f %s\n",
          left_pad, "Total Runtime (secs):",
          center_pad, t_times->program, "secs");
}
