/*******************************************************************************
 *  - FILE:      work_loader.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Loads query sequences and target hmm profiles.
 *  NOTES:
 *    - Currently, the load_by_name() functions are the only ones in use.
 *    - The names could be stored in a map for quicker lookups.
 *  TODOS:
 *    - Eventually, plan for these functions is that state is handled entirely by the
 *      WORKER object, so other arguments are unneccessary (encapsulated by WORKER).
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
#include "../algs_sparse/_algs_sparse.h"
#include "../reporting/_reporting.h"

/* header */
#include "_work.h"
#include "work_loader.h"

/*! FUNCTION:  	WORK_load_args()
 *  SYNOPSIS:  	Updates <worker> object to reflect <args>.
 */
void WORK_load_args(WORKER* worker) {
  ARGS* args = worker->args;
  RESULT* result = worker->result;

  result->threshold_viterbi = args->threshold_vit;
  result->threshold_cloud = args->threshold_cloud;
  result->threshold_fwdback = args->threshold_fwd;
}

/*! FUNCTION:  	WORK_load_mmseqs_file()
 *  SYNOPSIS:  	Loads mmseqs input m8 file into <results_in>, located at <mmseqs_m8_filein>.
 *                Verifies that search index range does not exceed bounds of list.
 */
void WORK_load_mmseqs_file(WORKER* worker) {
  ARGS* args = worker->args;
  STATS* stats = worker->stats;

  /* If range values are negative, then range is set to full list set */
  if (args->list_range.beg < 0 && args->list_range.end < 0) {
    args->list_range.beg = 0;
    args->list_range.end = INT_MAX;
  }

  /* load .m8 file from index range (list_range.beg, list_range.end) */
  RESULTS_M8_Parse(
      worker->mmseqs_data, args->mmseqs_m8_filein, args->list_range.beg, args->list_range.end);
  /* this is a fix because query and target are cross-labeled between MMSEQS and MMORE */


  M8_RESULTS_Swap_Target_and_Query(worker->mmseqs_data);

  /* Truncate or extract valid result range */
  args->list_range.beg = MAX(args->list_range.beg, 0);
  args->list_range.end = MIN(args->list_range.end, worker->mmseqs_data->N + args->list_range.beg);
  /* Compute number of searches */
  int num_searches = args->list_range.end - args->list_range.beg;
  worker->n_searches = num_searches;
  stats->n_searches = num_searches;
}

/*! FUNCTION:  	WORK_load_mmseqs_by_id()
 *  SYNOPSIS:  	Load <i>th mmseqs input from .m8 <mmseqs_data> list into <worker>.
 */
void WORK_load_mmseqs_by_id(WORKER* worker,
                            int id) {
  ARGS* args = worker->args;
  RESULT* result = worker->result;
  ALL_SCORES* scores = &result->scores;
  SCORES* finalsc = &result->final_scores;

  /* load mmseqs data */
  worker->mmseqs_id = id;
  worker->mmseqs_prv = worker->mmseqs_cur;
  worker->mmseqs_cur = &VEC_X(worker->mmseqs_data, id);

  /* load viterbi scores */
  finalsc->viterbi_natsc = worker->mmseqs_cur->bitsc;
  finalsc->viterbi_eval = worker->mmseqs_cur->eval;
  result->score_viterbi = worker->mmseqs_cur->eval;

  /* report input */
  if (args->verbose_level >= VERBOSE_HIGH) {
    fprintf(stdout, "=== M8 Entry : [%d] ===\n", worker->mmseqs_id);
    M8_RESULT_Dump(worker->mmseqs_cur, stdout);
  }
}

/*! FUNCTION:  	WORK_load_mmseqs_alignment()
 *  SYNOPSIS:  	Load mmseq's viterbi alignment <trace_vit> into <worker>.
 *                Query and Target should already be loaded.
 */
void WORK_load_mmseqs_alignment(WORKER* worker) {
  ARGS* args = worker->args;
  SEQUENCE* q_seq = worker->q_seq;
  int Q = q_seq->N;
  HMM_PROFILE* t_prof = worker->t_prof;
  int T = t_prof->N;
  M8_RESULT* mm_m8 = worker->mmseqs_cur;

  if (args->verbose_level >= VERBOSE_ALL) {
    M8_RESULT_Dump(mm_m8, stdout);
  }

  /* load viterbi alignment */
  TRACE aln_beg, aln_end;
  ALIGNMENT_Reuse(worker->trace_vit, Q, T);
  /* if no cigar alignment, just use the beginning and end of the alignment */
  aln_beg.t_0 = mm_m8->t_beg;
  aln_beg.q_0 = mm_m8->q_beg;
  aln_beg.st = M_ST;
  ALIGNMENT_AddTrace(worker->trace_vit, aln_beg);
  aln_end.t_0 = mm_m8->t_end;
  aln_end.q_0 = mm_m8->q_end;
  aln_end.st = M_ST;
  ALIGNMENT_AddTrace(worker->trace_vit, aln_end);
  /* since only start and end point of alignment is known, set endpoints by default */
  ALIGNMENT_AddRegion(worker->trace_vit, 0, 1, mm_m8->eval);
  ALIGNMENT_SetRegion(worker->trace_vit, 0);
}

/*! FUNCTION:  	WORK_load_target()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <mmseqs_data>'s target name field.
 */
void WORK_load_target(WORKER* worker) {
  ARGS* args = worker->args;
  TIMES* times = worker->times;
  CLOCK* timer = worker->timer;
  HMM_PROFILE* t_prof = worker->t_prof;


  CLOCK_Start(timer);

  /* update target entry name */
  worker->t_name_prv = worker->t_name;
  worker->t_name = worker->mmseqs_cur->target_name;
  /* if current and previous mmseqs entry targets are not the same, then load new target */
  if (STRING_Equal(worker->t_name, worker->t_name_prv) == false) {
    /* load new target */
    WORK_load_target_by_name(worker, worker->t_name);
  }

  CLOCK_Stop(timer);
  times->load_target = CLOCK_Duration(timer);

  /* report input */
  if (args->verbose_level >= VERBOSE_HIGH) {
    fprintf(stdout, "=== TARGET : [%d] ===\n", worker->t_id);
    F_INDEX_Node_Dump(worker->t_index, worker->t_id, stdout);
    fprintf(stdout, "[%d] { NAME: %s, LENGTH: %d }\n", worker->t_id, t_prof->name, t_prof->N);
  }
}

/*! FUNCTION:  	WORK_load_query()
 *  SYNOPSIS:  	Loads <query> HMM_PROFILE by <mmseqs_data>'s query name field.
 */
void WORK_load_query(WORKER* worker) {
  ARGS* args = worker->args;
  TIMES* times = worker->times;
  CLOCK* timer = worker->timer;
  SEQUENCE* q_seq = worker->q_seq;

  CLOCK_Start(timer);

  /* update query entry name */
  worker->q_name_prv = worker->q_name;
  worker->q_name = worker->mmseqs_cur->query_name;
  /* if current and previous mmseqs entry targets are not the same */
  if (STRING_Equal(worker->q_name, worker->q_name_prv) == false) {
    /* load new target */
    WORK_load_query_by_name(worker, worker->q_name);
  }

  CLOCK_Stop(timer);
  times->load_query = CLOCK_Duration(timer);

  /* report input */
  if (args->verbose_level >= VERBOSE_HIGH) {
    fprintf(stdout, "=== QUERY  : [%d] ===\n", worker->q_id);
    F_INDEX_Node_Dump(worker->q_index, worker->q_id, stdout);
    fprintf(stdout, "[%d] { NAME: %s, LENGTH: %d }\n", worker->q_id, q_seq->name, q_seq->N);
  }
}

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <name> field.
 *                Depends on <task> settings in <worker>.
 */
void WORK_load_target_by_name(WORKER* worker,
                              STR name) {

  /* find target id by searching target name in file index */
  int findex_id = F_INDEX_Search_Name(worker->t_index, name);
  if (findex_id == -1) {
    fprintf(stderr, "ERROR: Target name '%s' not found in F_INDEX.\n", name);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
  /* use target id to load target */
  WORK_load_target_by_findex_id(worker, findex_id);
}

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <name> field.
 *                Depends on <task> settings in <worker>.
 */
void WORK_load_query_by_name(WORKER* worker,
                             STR name) {

  /* find query id by searching query name in file index */
  int findex_id = F_INDEX_Search_Name(worker->q_index, name);
  if (findex_id == -1) {
    fprintf(stderr, "ERROR: Query name '%s' not found in F_INDEX.\n", name);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
  WORK_load_query_by_findex_id(worker, findex_id);
}

/*! FUNCTION:  	WORK_load_target_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void WORK_load_target_by_findex_id(WORKER* worker,
                                   int index_id) {
  ARGS* args = worker->args;
  F_INDEX_NODE* my_idx = &worker->t_index->nodes[index_id];
  worker->t_id_prv = worker->t_id;
  worker->t_id = index_id;

  /* load target profile by file type */
  switch (args->t_filetype) {
    case FILE_HMM: {
      HMM_PROFILE_Parse(worker->t_prof, args->t_mmore_filein, my_idx->offset);
      HMM_PROFILE_Convert_NegLog_To_Real(worker->t_prof);
      HMM_PROFILE_Config(worker->t_prof, args->search_mode);
    } break;
    case FILE_FASTA: {
      SEQUENCE_Fasta_Parse(worker->t_seq, args->t_filein, my_idx->offset);
      SEQUENCE_to_HMM_PROFILE(worker->t_seq, worker->t_prof);
      HMM_PROFILE_Dump(worker->t_prof, stdout);
      ERRORCHECK_exit(EXIT_SUCCESS);
    } break;
    default: {
      fprintf(stderr, "ERROR: Only HMM and FASTA filetypes are supported for targets.\n");
      ERRORCHECK_exit(EXIT_FAILURE);
    }
  }
}

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void WORK_load_query_by_findex_id(WORKER* worker,
                                  int index_id) {
  ARGS* args = worker->args;
  HMM_PROFILE* t_prof = worker->t_prof;
  SEQUENCE* t_seq = worker->t_seq;
  SEQUENCE* q_seq = worker->q_seq;
  F_INDEX_NODE* my_idx = &worker->q_index->nodes[index_id];

  worker->q_id_prv = worker->q_id;
  worker->q_id = index_id;

  /* load query by file type */
  switch (args->q_filetype) {
    /* fasta only supported file type */
    case FILE_FASTA: {
      SEQUENCE_Fasta_Parse(worker->q_seq, args->q_mmore_filein, my_idx->offset);
      // SEQUENCE_Dump( worker->q_seq, stdout );
    } break;
    case FILE_HMM: {
    }
    default: {
      fprintf(stderr, "ERROR: Only FASTA filetypes are supported for queries.\n");
      ERRORCHECK_exit(EXIT_FAILURE);
    }
  }

  /* set special state transitions based on query sequence length */
  if (t_prof != NULL) {
    HMM_PROFILE_ReconfigLength(worker->t_prof, worker->q_seq->N);
  } else {
    fprintf(stderr, "ERROR: Target profile must be loaded before Query Sequence. Currently NULL.\n");
    ERRORCHECK_exit(EXIT_FAILURE);
  }
}
