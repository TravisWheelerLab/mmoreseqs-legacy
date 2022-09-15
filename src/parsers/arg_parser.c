/*******************************************************************************
 *  - FILE:  arg_parser.h
 *  - DESC:  Parses command line arguments.
 *    NOTES:
 *       - To be phased out in preference of cli using CLI11.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
/* local imports */
#include "../objects/structs.h"
#include "../objects/_objects.h"
/* header */
#include "arg_parser.h"

void ARGS_SetDefaults(ARGS* args) {
  /* --- PROGRAM --- */
  args->mmoreseqs_program = STR_Create("mmoreseqs");
  args->mmseqs_program = STR_Create("mmseqs");
  args->hmmer_program = STR_Create("hmmbuild");
  args->mmoreseqs_scripts = STR_Create(SCRIPT_DIR);
  /* --- PIPELINE MODE --- */
  args->pipeline_mode = PIPELINE_EASY_SEARCH;
  args->pipeline_name = NULL;
  args->verbose_level = VERBOSE_LOW;
  args->num_threads = 1;
  args->search_mode = MODE_UNILOCAL;
  args->qt_search_space = SELECT_ALL_V_ALL;
  args->tmp_folderpath = NULL;
  args->tmp_remove = false;

  /* --- PREPARATION OPTIONS --- */
  /* files/folders */
  args->is_run_prep = true;
  args->is_prep_copy = true;
  args->prep_folderpath = NULL;
  args->target_prep = NULL;
  args->query_prep = NULL;
  /* types */
  args->target_prep_type = FILE_MSA;
  args->query_prep_type = FILE_FASTA;

  /* --- PIPELINE OPTIONS --- */
  args->search_name = STR_Create("P2S");
  args->is_run_bias = true;
  args->is_run_pruned = true;
  args->is_run_full = false; /* DEBUG */
  args->is_run_domains = true;
  args->is_run_mmseqsaln = false;
  args->is_run_vit_mmore = false; 
  args->is_run_vit = false;  
  args->is_run_vitaln = true;    
  args->is_run_optacc = false;   
  args->is_run_post = false;
  args->is_run_postaln = false;
  args->is_run_mmseqs = true;
  args->is_run_mmseqs_pref = true;
  args->is_run_mmseqs_align = true;
  args->is_run_mmseqs_ungapped = true;
  args->is_run_mmseqs_sens = false;
  args->is_run_mmore = true;
  args->is_run_convert = true;
  args->use_pvals = true;

  /* --- DEBUG OPTIONS --- */
  args->is_use_local_tools = false;
  args->is_recycle_mx = false;
  args->is_debug = true;
  args->dbg_folderpath = STR_Create(DEBUG_FOLDER "/");
  args->dbg_arg_dump = false;

  /* --- INPUT --- */
  /* filepath */
  args->t_filein = NULL;
  args->q_filein = NULL;
  args->t_mmore_filein = NULL;
  args->q_mmore_filein = NULL;
  args->t_mmseqs_p_filein = NULL;
  args->t_mmseqs_s_filein = NULL;
  args->q_mmseqs_filein = NULL;
  /* filetype */
  args->is_guess_filetype = true;
  args->t_filetype = FILE_HMM;
  args->q_filetype = FILE_FASTA;
  args->t_mmore_filetype = FILE_HMM;
  args->q_mmore_filetype = FILE_FASTA;
  args->t_mmseqs_p_filetype = FILE_MMDB_P;
  args->t_mmseqs_s_filetype = FILE_MMDB_S;
  args->q_mmseqs_filetype = FILE_MMDB_S;
  /* indexes */
  args->t_index_filein = NULL;
  args->q_index_filein = NULL;
  /* db sizes */
  args->q_dbsize = -1;
  args->t_dbsize = -1;

  /* --- INTERRIM OUTPUT --- */
  args->mmseqs_m8_filein = STR_Create("mmseqs.results.m8out");
  args->is_mmseqs_m8out = false;
  args->mmseqs_m8out = STR_Create("mmseqs.results.m8out");

  /* --- OUTPUT --- */
  /* default outputs */
  args->is_redirect_stdout = false;
  args->stdout_fileout = STR_Create("mmore.results.stdout");
  args->is_redirect_stderr = false;
  args->stderr_fileout = STR_Create("mmore.results.stderr");
  /* file outputs */
  args->is_allout = false;
  args->allout_fileout = STR_Create("mmore.results");
  args->is_hmmerout = false;
  args->hmmerout_fileout = STR_Create("mmore.results.hmmerout");
  args->is_m8out = false;
  args->m8out_fileout = STR_Create("mmore.results.m8");
  args->is_myout = false;
  args->myout_fileout = STR_Create("mmore.results.myout");
  args->is_mydom = false;
  args->mydom_fileout = STR_Create("mmore.results.mydomout");
  args->is_mytimeout = false;
  args->mytime_fileout = STR_Create("mmore.results.mytimeout");
  args->is_mythreshout = false;
  args->mythresh_fileout = STR_Create("mmore.results.mythreshout");
  // args->is_customout = false;
  // args->customout_fileout = STR_Create("results.customout");

  /* --- RANGE OPTIONS --- */
  args->t_range = (RANGE){-1, -1};
  args->q_range = (RANGE){-1, -1};
  // args->use_range = false;
  args->list_range = (RANGE){-1, -1};

  /* --- MMORE / FB-PRUNER --- */
  args->alpha = 12.0f;
  args->beta = 16.0f;
  args->gamma = 5;
  args->hard_limit = -12.0f;
  args->mmore_evalue = 2e2f;
  args->mmore_pvalue = 1e-3f;

  /* --- MMSEQS --- */
  args->mmseqs_maxhits = 1000; /* mmseqs default: 300 */
  args->mmseqs_altalis = 0;
  args->mmseqs_kmer = 7;
  args->mmseqs_kscore = 80;       /* mmseqs default: 95 */
  args->mmseqs_ungapped_vit = 0;  /* mmseqs default: 15 */
  args->mmseqs_sensitivity = 7.5; /* mmseqs default: 4 (overrides kscore) */
  args->mmseqs_p2s_pvalue = 1e-2f;
  args->mmseqs_s2s_pvalue = 1e7;
  args->mmseqs_pvalue = 1e-2f; /* mmseqs default: 1E-3 */
  args->mmseqs_evalue = 1e-2f;

  /* --- SCORE FILTERS (P_VALUES) --- */
  args->is_run_filter = false;
  args->is_run_viterbi_filter = true;
  args->is_run_cloud_filter = true;
  args->is_run_boundfwd_filter = true;
  args->is_run_fwdback_filter = true;
  args->is_run_report_filter = true;
  args->threshold_pre = 80;      /* mmseqs default: 95 */
  args->threshold_ungapped = 15; /* mmseqs default: 15 */
  args->threshold_p2s = 1e-2f;   /* mmseqs default: 1E-3 */
  args->threshold_s2s = INF;
  args->threshold_vit = 1e-2f;
  args->threshold_cloud = 1e-3f;
  args->threshold_boundfwd = 1e-4f;
  args->threshold_fwd = 1e-4f;
  args->threshold_report_eval = 10;
}

void ARGS_Dump(ARGS* args,
               FILE* fp) {
  int pad = 20;
  bool align = -1; /* -1 for right alignment, 1 for left alignment */

  fprintf(fp, "# === MMORE OPTIONS =======================\n");
  /* --- PIPELINE --- */
  fprintf(fp, "# === PIPELINE ===\n");
//  fprintf(fp, "# %*s:\t%s %s\n", align * pad, "PIPELINE", "pipeline :)");
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "VERBOSITY_MODE", VERBOSITY_NAMES[args->verbose_level], args->verbose_level);
  fprintf(fp, "# === PROGRAMS ===\n");
  fprintf(fp, "# %*s:\t%s\n", align * pad, "MMORESEQS_PROGRAM", args->mmoreseqs_program);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "MMSEQS_PROGRAM", args->mmseqs_program);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "HMMER_PROGRAM", args->hmmer_program);
  fprintf(fp, "# === SCRIPTS ===\n");
  fprintf(fp, "# %*s:\t%s\n", align * pad, "PREP_SCRIPT", PREP_SCRIPT);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "PREPSEARCH_SCRIPT", PREPSEARCH_SCRIPT);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "SEARCH_SCRIPT", SEARCH_SCRIPT);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "EASYSEARCH_SCRIPT", EASYSEARCH_SCRIPT);
  fprintf(fp, "# === PIPELINE OPTIONS ===\n");
  fprintf(fp, "# %*s:\t[%d] : [%d] [%d]\n",
          align * pad, "RUN_MMSEQS", args->is_run_mmseqs, args->is_run_mmseqs_pref, args->is_run_mmseqs_align);
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "RUN_MMORE", args->is_run_mmore);
  /* --- INPUT --- */
  fprintf(fp, "# === INPUT ===\n");
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "TARGET", args->t_filein, FILETYPE_NAME_Get(args->t_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "QUERY", args->q_filein, FILETYPE_NAME_Get(args->q_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "TARGET_PREP", args->target_prep, FILETYPE_NAME_Get(args->target_prep_type));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "QUERY_PREP", args->query_prep, FILETYPE_NAME_Get(args->query_prep_type));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "TARGET_MMORE", args->t_mmore_filein, FILETYPE_NAME_Get(args->t_mmore_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "QUERY_MMORE", args->q_mmore_filein, FILETYPE_NAME_Get(args->q_mmore_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "TARGET_MMSEQS_P", args->t_mmseqs_p_filein, FILETYPE_NAME_Get(args->t_mmseqs_p_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "TARGET_MMSEQS_S", args->t_mmseqs_s_filein, FILETYPE_NAME_Get(args->t_mmseqs_s_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "QUERY_MMSEQS", args->q_mmseqs_filein, FILETYPE_NAME_Get(args->t_mmseqs_s_filetype));
  fprintf(fp, "# %*s:\t%s\n", align * pad, "T_INDEX_PATH", args->t_index_filein);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "Q_INDEX_PATH", args->q_index_filein);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "MMSEQS_M8", args->mmseqs_m8_filein);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "TMP_FOLDER", args->tmp_folderpath);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "PREP_FOLDER", args->prep_folderpath);
  fprintf(fp, "# \n");
  /* --- MMSEQS --- */
  fprintf(fp, "# === MMSEQS ===\n");
  fprintf(fp, "# %*s:\t%d\n", align * pad, "MMSEQS_KMER", args->mmseqs_kmer);
  fprintf(fp, "# %*s:\t%d [%d]\n", align * pad, "MMSEQS_KSCORE", args->mmseqs_kscore, !args->is_run_mmseqs_sens);
  fprintf(fp, "# %*s:\t%.3f [%d]\n", align * pad, "MMSEQS_SENS", args->mmseqs_sensitivity, args->is_run_mmseqs_sens);
  fprintf(fp, "# %*s:\t%d [%d]\n", align * pad, "MMSEQS_UNGAPPED", args->mmseqs_ungapped_vit, args->is_run_mmseqs_ungapped);
  fprintf(fp, "# %*s:\t%d\n", align * pad, "MMSEQS_ALTALIS", args->mmseqs_altalis);
  fprintf(fp, "# %*s:\t%.2e\n", align * pad, "MMSEQS_P2S_PVAL", args->threshold_p2s);
  fprintf(fp, "# %*s:\t%.2e\n", align * pad, "MMSEQS_S2S_PVAL", args->threshold_s2s);
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "MMSEQS_VITALN", args->is_run_mmseqsaln);
  /* --- MMORE --- */
  fprintf(fp, "# === MMORE ===\n");
  fprintf(fp, "# %*s:\t%.2f\n", align * pad, "MMORE_ALPHA", args->alpha);
  fprintf(fp, "# %*s:\t%.2f\n", align * pad, "MMORE_BETA", args->beta);
  fprintf(fp, "# %*s:\t%d\n", align * pad, "MMORE_GAMMA", args->gamma);
  fprintf(fp, "# %*s:\t%.2f\n", align * pad, "MMORE_HARD_LIMIT", args->hard_limit);
  fprintf(fp, "# %*s:\t%.2e [%d]\n", align * pad, "MMORE_VITERBI_PVAL", args->threshold_vit, args->is_run_viterbi_filter);
  fprintf(fp, "# %*s:\t%.2e [%d]\n", align * pad, "MMORE_CLOUD_PVAL", args->threshold_cloud, args->is_run_cloud_filter);
  fprintf(fp, "# %*s:\t%.2e [%d]\n", align * pad, "MMORE_BOUNDFWD_PVAL", args->threshold_boundfwd, args->is_run_boundfwd_filter);
  fprintf(fp, "# %*s:\t%.2e [%d]\n", align * pad, "MMORE_REPORT_EVAL", args->threshold_report_eval, args->is_run_report_filter);

  fprintf(fp, "# %*s:\t(%d,%d)\n", align * pad, "MMORE_RANGE", args->list_range.beg, args->list_range.end);
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "MMORE_FULL", args->is_run_full);
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "MMORE_VIT_MMORE", args->is_run_vit_mmore);
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "MMORE_DOMAINS", args->is_run_domains);
  fprintf(fp, "# %*s:\t[%d] [%d]\n", align * pad, "MMORE_VITALN", args->is_run_vitaln, args->is_run_vit);
  fprintf(fp, "# %*s:\t[%d] [%d]\n", align * pad, "MMORE_POSTALN", args->is_run_postaln, args->is_run_optacc);
  fprintf(fp, "# \n");
  /* --- OUTPUT --- */
  fprintf(fp, "# === OUTPUT ===\n");
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "OUTPUT_FILEPATH", args->stdout_fileout, args->is_redirect_stdout);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "HMMEROUT_FILEPATH", args->hmmerout_fileout, args->is_hmmerout);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "M8OUT_FILEPATH", args->m8out_fileout, args->is_m8out);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "MYOUT_FILEPATH", args->myout_fileout, args->is_myout);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "MYDOMOUT_FILEPATH", args->mydom_fileout, args->is_mydom);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "MYTHRESHOUT_FILEPATH", args->mythresh_fileout, args->is_mythreshout);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "MYTIMEOUT_FILEPATH", args->mytime_fileout, args->is_mytimeout);
  fprintf(fp, "# ==============================================\n\n");
}

