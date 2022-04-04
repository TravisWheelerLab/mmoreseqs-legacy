/*******************************************************************************
 *  - FILE:      myout.c
 *  - DESC:    Reporting Subroutines for generating output.
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
#include "myout.h"

/* === MYOUT FUNCTIONS === */

/* === MY OUTPUT === */
/* Custom-style output, similar to BLAST-style .m8
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

/*!   FUNCTION:   REPORT_myout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after m8, see example)
 */
void REPORT_myout_header(WORKER* worker,
                         FILE* fp) {
  const int num_fields = 21;
  const char* headers[] = {
      "result-id",
      "target-hmm",
      "query-seq",
      "t-len",
      "q-len",
      "evalue",
      "pre-sc",
      "comp-bias",
      "seq-sc",
      "dom-sum-sc",
      "mmseqs-vit-sc",
      "mmore-vit-sc",
      "total-cells",
      "MMORE-cells",
      "perc-cells",
      "t-aln",
      "q-aln",
      "t-bounds",
      "q-bounds",
      "time",
      "time-noload"};

  REPORT_header(fp, headers, num_fields);
}

/*!   FUNCTION:   REPORT_myout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_myout_entry(WORKER* worker,
                        RESULT* result,
                        FILE* fp) {
  TIMES* times = worker->times;
  HMM_PROFILE* t_prof = worker->t_prof;
  SEQUENCE* q_seq = worker->q_seq;
  ALIGNMENT* aln;
  ALIGNMENT* aln_mm;
  ALIGNMENT* aln_vit = worker->trace_vit;
  ALIGNMENT* aln_fwd = worker->trace_post;
  DOMAIN_DEF* dom_def = worker->dom_def;
  ALL_SCORES* scores = &result->scores;
  SCORES* final = &result->final_scores;

  if (aln_fwd->is_hmmer_aln == true) {
    aln = aln_fwd;
  }
  elif (aln_vit->is_hmmer_aln == true) {
    aln = aln_vit;
  }
  else {
    aln = aln_vit;
  }

  int aln_len = aln->end - aln->beg + 1;
  TRACE aln_beg = VEC_X(aln->traces, aln->beg);
  TRACE aln_end = VEC_X(aln->traces, aln->end);

  RANGE q_bounds, t_bounds;

  float time_noload = times->loop - times->load_query - times->load_target;

  if (result->target_range.end - result->target_range.beg <= 0) {
    EDGEBOUNDS_Find_BoundingBox(worker->edg_row, &q_bounds, &t_bounds);
  }

  fprintf(fp, "%d\t%s\t%s\t%d\t%d\t%.2e\t%.3f\t%.3f\t%.3f\t%.3f\t%.2e\t%.2e\t%d\t%d\t%.5f\t%d-%d\t%d-%d\t%d-%d\t%d-%d\t%.5f\t%.5f\n",
          worker->mmseqs_id,           /* id index in mmseqs list */
          t_prof->name,                /* target name */
          q_seq->name,                 /* query name */
          worker->t_prof->N,           /* target length */
          worker->q_seq->N,            /* query length */
          final->eval,                 /* evalue */
          final->pre_sc,               /* seq scores before bias correction (in bits) */
          final->null2_seq_bias_bitsc, /* seq bias (in bits) */
          final->seq_sc,               /* seq score after correction (in bits) */
          dom_def->dom_sumsc,          /* sum of all domain scores; if domains were not computed, zero */
          final->viterbi_eval,         /* viterbi eval (in mmore, this comes from mmseqs) */
          final->viterbi_mmore_eval,   /* viterbi eval (in mmore, this comes from mmseqs) */
          result->total_cells,         /* total number of cells computed by full viterbi */
          result->cloud_cells,         /* total number of cells computed by mmore */
          result->perc_cells,          /* percent of total cells computed by mmore */
          aln_beg.t_0, aln_end.t_0,    /* target alignment range */
          aln_beg.q_0, aln_end.q_0,    /* query alignment range */
          t_bounds.beg, t_bounds.end,  /* target bounds of cloud search area */
          q_bounds.beg, q_bounds.end,  /* query bounds of cloud search area */
          times->loop,                 /* time for entire iteration */
          time_noload                  /* time with load times */
  );

  /* TODO: WIP */
  // const int num_fields = 12;
  // const GEN fields[] = {
  //    GEN_Wrap( &t_prof->name,                    DATATYPE_STRING,     sizeof(char*) ),
  //    GEN_Wrap( &q_seq->name,                     DATATYPE_STRING,     sizeof(char*) ),
  //    GEN_Wrap( &result->final_scores.eval,       DATATYPE_FLOAT_EXP,  sizeof(float) ),
  //    GEN_Wrap( &result->final_scores.pre_sc,     DATATYPE_FLOAT,      sizeof(float) ),
  //    GEN_Wrap( &result->final_scores.seq_bias,   DATATYPE_FLOAT,      sizeof(float) ),
  //    GEN_Wrap( &dom_def->dom_sumsc,              DATATYPE_FLOAT,      sizeof(float) ),
  //    GEN_Wrap( &result->vit_natsc,               DATATYPE_FLOAT,      sizeof(float) ),
  //    GEN_Wrap( &result->total_cells,             DATATYPE_INT,        sizeof(int) ),
  //    GEN_Wrap( &result->cloud_cells,             DATATYPE_INT,        sizeof(int) ),
  //    GEN_Wrap( &result->target_range,            DATATYPE_RANGE,      sizeof(int) ),
  //    GEN_Wrap( &result->query_range,             DATATYPE_RANGE,      sizeof(int) ),
  //    GEN_Wrap( &result->time,                    DATATYPE_FLOAT,      sizeof(float) ),
  // };

  // REPORT_entry( fp, fields, num_fields, sig_digits );
}

/*!   FUNCTION:   REPORT_myout_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)
 */
void REPORT_myout_footer(WORKER* worker,
                         FILE* fp) {
  fprintf(fp, "# [ok] [myout]\n");
}
