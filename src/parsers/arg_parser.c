/*******************************************************************************
 *  - FILE:  arg_parser.h
 *  - DESC:  Parses command line arguments.
 *    NOTES:
 *       - To be phased out in preference of cli using CLI11.
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
#include <argp.h>
/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../macros/_macros.h"
/* header */
#include "_parsers.h"
#include "arg_parser.h"

void ARGS_Parse(ARGS* args,
                int argc,
                char* argv[],
                COMMANDLINE* cmd,
                ARG_OPTS* arg_opts) {
  DBG_PRINTF(stdout, "# === NUM_ARGS: %d\n", argc);
  ARGS_SetDefaults(args);

  /* if no arguments given */
  if (argc <= 1) {
    ARGS_HelpInfo();
    exit(EXIT_SUCCESS);
  }

  /* check for help flag */
  if (STR_Equals(argv[1], "-h") || STR_Equals(argv[1], "--help")) {
    ARGS_HelpInfo();
    exit(EXIT_SUCCESS);
  }
  /* check for version flag */
  if (STR_Equals(argv[1], "--version")) {
    ARGS_VersionInfo();
    exit(EXIT_SUCCESS);
  }

  /* current argument index */
  int arg_cur = 1;

  /* parse command */
  ARGS_ParseCommand(args, argc, argv, &arg_cur);
  /* parse flags and options */
  ARGS_ParseOptions(args, argc, argv, &arg_cur);

  if (args->dbg_arg_dump) {
    ARGS_Dump(args, stdout);
    exit(EXIT_SUCCESS);
  }
}

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
  args->mmseqs_kscore = 75;       /* mmseqs default: 95 */
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
  args->threshold_report_eval = 2e2f;
}

void ARGS_Dump(ARGS* args,
               FILE* fp) {
  int pad = 20;
  bool align = -1; /* -1 for right alignment, 1 for left alignment */

  fprintf(fp, "# === MMORE OPTIONS =======================\n");
  /* --- PIPELINE --- */
  fprintf(fp, "# === PIPELINE ===\n");
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "PIPELINE", PIPELINES[args->pipeline_mode].name, args->pipeline_mode);
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
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "RUN_CONVERT", args->is_run_convert);
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

FILETYPE ARGS_FindFiletype(STR filename) {
  for (int i = 0; i < NUM_FILETYPE_EXTS; i++) {
    char* ext_name = FILETYPE_EXTS[i].s;
    int ext_type = FILETYPE_EXTS[i].i;
    if (STRING_EndsWith(filename, ext_name, strlen(ext_name)) == 0) {
      return ext_type;
    }
  }

  fprintf(stderr, "WARNING: '%s' filetype could not be found. Set to TYPE_NULL.\n", filename);

  return FILE_NULL;
}

void ARGS_HelpInfo() {
  /* basic usage */
  printf("Usage: <command> <args...> [options]\n\n");

  printf("COMMANDS:\n");
  for (int i = 0; i < NUM_PIPELINES; i++) {
    printf("%s\n", PIPELINES_ARG_HELP[i]);
  }
}

void ARGS_CommandHelpInfo(ARGS* args) {
  fprintf(stdout, "USAGE: %s\n", PIPELINES_ARG_HELP[args->pipeline_mode]);
}

void ARGS_VersionInfo() {
  fprintf(stdout, "MMORESeqs: %s\n", BUILD_DESC);
  fprintf(stdout, "%s\n", BUILD_COPYRIGHT);
  fprintf(stdout, "VERSION: v%s\n", BUILD_VERSION);
  fprintf(stdout, "BUILD_TYPE: %s\n", BUILD_TYPE);
  fprintf(stdout, "GIT_HASH: %s\n", BUILD_HASH);
  ERRORCHECK_exit(EXIT_SUCCESS);
}

/* CLI using argp.h */

static char doc[] = "MMOREseqs -- " BUILD_DESC;

enum OPT_KEYS {
  _KEY,

  VERBOSE_KEY,
  VERSION_KEY,
  NUM_THREADS_KEY,
  EVAL_KEY,
  USE_PVALS_KEY,

  RUN_PREP_KEY,
  RUN_MMSEQS_KEY,
  RUN_MMSEQS_PREF_KEY,
  RUN_MMSEQS_ALIGN_KEY,
  RUN_CONVERT_KEY,
  RUN_MMORE_KEY,
  RUN_VIT_KEY,
  RUN_POST_KEY,

  PREP_KEY,
  TMP_KEY,
  MMSEQS_M8_KEY,
  INDEX_KEY,

  MMSEQS_KMER_KEY,
  MMSEQS_KSCORE_KEY,
  MMSEQS_SENS_KEY,
  MMSEQS_UNGAPPED_VIT_KEY,
  MMSEQS_EVAL_KEY,
  MMSEQS_PVAL_KEY,
  MMSEQS_HITS_PER_SEARCH_KEY,
  MMSEQS_ALTALIS_KEY,
  MMSEQS_SPLIT_KEY,

  VIT_FILTER_KEY,
  CLD_FILTER_KEY,
  FWD_FILTER_KEY,
  ALPHA_KEY,
  BETA_KEY,
  GAMMA_KEY,
  HARD_LIMIT_KEY,
  RANGE_KEY,

  STDOUT_KEY,
  STDERR_KEY,
  M8OUT_KEY,
  DOMTBLOUT_KEY,
  MMSEQS_M8OUT_KEY,

  PREP_COPY_KEY,
  PREP_LINK_TARGET_MMORE_KEY,
  PREP_LINK_QUERY_MMORE_KEY,
  PREP_LINK_TARGET_MMSEQS_KEY,
  PREP_LINK_QUERY_MMSEQS_KEY,

  RUN_DOMAINS_KEY,
  RUN_BIAS_KEY,
  RUN_FILTER_KEY,
  RUN_VIT_FILTER_KEY,
  RUN_CLD_FILTER_KEY,
  RUN_FWD_FILTER_KEY,
  RUN_VIT_MMORE_KEY,
  RUN_FULL_KEY,
  RUN_MMSEQSALN_KEY,
  RUN_VITALN_KEY,
  RUN_POSTALN_KEY,

  SEARCH_MODE_KEY,
  SEARCH_TYPE_KEY,
  GUESS_FTYPE_KEY,
  MMORESEQS_FTYPE_KEY,
  MMORESEQS_MAIN_FTYPE_KEY,
  PROGRAM_MMSEQS_KEY,
  PROGRAM_HMMER_KEY,
  PROGRAM_MMORESEQS_KEY,
  LOCAL_TOOLS_KEY,
  SCRIPT_DIR_KEY,
  DBSIZES_KEY,
  MMSEQS_TIMES_KEY,
  MMSEQS_DBSIZES_KEY,

  DEBUG_KEY,
  DEBUG_VIZ_KEY,
  ALLOUT_KEY,
  MYOUT_KEY,
  MYTIMEOUT_KEY,
  MYTHRESHOUT_KEY,
  CUSTOMOUT_KEY,
  DEBUGOUT_KEY,

  NULL_KEY
};

STATUS_FLAG ARGS_ParseCommand(ARGS* args,
                              int argc, 
                              char* argv[],
                              int* arg_cur_p) {

  /* current arg */
  int arg_cur = arg_cur_p[0];
  /* args remaining */
  int args_rem = argc - arg_cur;

  /* first argument is command pipeline */
  bool found_pipeline = false;
  for (int i = 0; i < NUM_PIPELINE_MODES; i++) {
    PIPELINE* pipeline = &PIPELINES[i];
    if (STR_Compare(argv[1], pipeline->name) == 0) {
      args->pipeline_mode = i;
      args->pipeline_name = STR_Create(pipeline->name);
      found_pipeline = true;
      break;
    }
  }

  /* check that valid pipeline mode was entered */
  if (found_pipeline == false || (args->pipeline_mode < 0) || (args->pipeline_mode >= NUM_PIPELINE_MODES)) {
    fprintf(stderr, "ERROR: Invalid pipeline/command was given: %s.\n", argv[1]);
    ARGS_CommandHelpInfo(args);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
  args_rem -= 1;
  arg_cur += 1;

  /* set number of main arguments based on given pipeline */
  int num_main_args = PIPELINES[args->pipeline_mode].num_main_args;

  /* check for help flag */
  if (argc >= 3 && (STR_Equals(argv[2], "-h") || STR_Equals(argv[2], "--help"))) {
    ARGS_CommandHelpInfo(args);
    exit(EXIT_SUCCESS);
  }

  /* check proper number of main args remain */
  if (args_rem < num_main_args) {
    fprintf(stderr, "ERROR: Improper number of main arguments. [required: %d/%d]\n", args_rem, num_main_args);
    ARGS_CommandHelpInfo(args);
    ERRORCHECK_exit(EXIT_FAILURE);
  }

  /* check if proper number of args */
  for (int i = 0; i < num_main_args; i++) {
    if (STR_StartsWith(argv[i + 1], "--") == 0) {
      fprintf(stderr, "ERROR: Improper number of main arguments (must preceed options). [required: %d/%d]\n", i, num_main_args);
      ARGS_CommandHelpInfo(args);
      ERRORCHECK_exit(EXIT_FAILURE);
    }
  }

  /* parse main commands */
  if (STR_Equals(args->pipeline_name, "search")) {
    args->t_filein = STR_Set(args->t_filein, argv[2]);
    args->q_filein = STR_Set(args->q_filein, argv[3]);
    args->t_mmore_filein = STR_Set(args->t_mmore_filein, argv[2]);
    args->q_mmore_filein = STR_Set(args->q_mmore_filein, argv[3]);
    args->t_mmseqs_p_filein = STR_Set(args->t_mmseqs_p_filein, argv[4]);
    args->t_mmseqs_s_filein = STR_Set(args->t_mmseqs_s_filein, argv[5]);
    args->q_mmseqs_filein = STR_Set(args->q_mmseqs_filein, argv[6]);

    args->t_mmore_p_filetype = FILE_HMM;
    args->q_filetype = FILE_FASTA;
    args->t_mmseqs_p_filetype = FILE_MMDB_P;
    args->t_mmseqs_s_filetype = FILE_MMDB_S;
    args->q_mmseqs_filetype = FILE_MMDB_S;
  }
  elif (STR_Equals(args->pipeline_name, "mmore-search")) {
    args->t_filein = STR_Set(args->t_filein, argv[2]);
    args->q_filein = STR_Set(args->q_filein, argv[3]);
    args->t_mmore_filein = STR_Set(args->t_mmore_filein, argv[2]);
    args->q_mmore_filein = STR_Set(args->q_mmore_filein, argv[3]);
    args->mmseqs_m8_filein = STR_Set(args->mmseqs_m8_filein, argv[4]);

    args->t_filetype = FILE_HMM;
    args->q_filetype = FILE_FASTA;
  }
  elif (STR_Equals(args->pipeline_name, "mmseqs-search")) {
    args->t_mmseqs_p_filein = STR_Set(args->t_mmseqs_p_filein, argv[2]);
    args->t_mmseqs_s_filein = STR_Set(args->t_mmseqs_s_filein, argv[3]);
    args->q_mmseqs_filein = STR_Set(args->q_mmseqs_filein, argv[4]);

    args->t_mmseqs_p_filetype = FILE_MMDB_P;
    args->t_mmseqs_s_filetype = FILE_MMDB_S;
    args->q_mmseqs_filetype = FILE_MMDB_S;
  }
  elif (STR_Equals(args->pipeline_name, "prep")) {
    args->target_prep = STR_Set(args->target_prep, argv[2]);
    args->query_prep = STR_Set(args->query_prep, argv[3]);
    args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[4]);
    args->prep_folderpath = STR_Set(args->prep_folderpath, argv[4]);

    args->target_prep_type = FILE_MSA;
    args->query_prep_type = FILE_FASTA;
  }
  elif (STR_Equals(args->pipeline_name, "prep-search")) {
    args->prep_folderpath = STR_Set(args->prep_folderpath, argv[2]);
    args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[2]);
  }
  elif (STR_Equals(args->pipeline_name, "easy-search")) {
    args->target_prep = STR_Set(args->target_prep, argv[2]);
    args->query_prep = STR_Set(args->query_prep, argv[3]);
    args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[4]);
    args->prep_folderpath = STR_Set(args->prep_folderpath, argv[4]);

    args->target_prep_type = FILE_MSA;
    args->query_prep_type = FILE_FASTA;
  }
  elif (STR_Equals(args->pipeline_name, "index")) {
    args->t_filein = STR_Set(args->t_filein, argv[2]);
    args->q_filein = STR_Set(args->q_filein, argv[3]);
    args->t_index_filein = STR_Set(args->t_index_filein, argv[4]);
    args->q_index_filein = STR_Set(args->q_index_filein, argv[5]);
  }
  else {
    fprintf(stderr, "ERROR: Command '%s' is currently not supported.\n", args->pipeline_name);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
  
  /* return remaining number of args */
  args_rem -= num_main_args;
  arg_cur += num_main_args;
  arg_cur_p[0] = arg_cur;
}

STATUS_FLAG ARGS_ParseOptions(ARGS* args,
                              int argc, 
                              char* argv[],
                              int* arg_cur_p) {
  /* current arg */
  int arg_cur = arg_cur_p[0];
  /* args remaining */
  int args_rem = argc - arg_cur;
  /* required arguments */
  int req_args = 0;

  char* flag = NULL;
  for (int i = arg_cur; i < argc; ++i) {
    /* if long flag */
    if (STR_ComparePrefix(argv[i], "--", 2) == 0) {
      /* === HELP OPTIONS === */
      if (STR_Equals(argv[i], (flag = "--help"))) {
        ARGS_CommandHelpInfo(args);
        ERRORCHECK_exit(EXIT_SUCCESS);
      }
      elif (STR_Equals(argv[i], (flag = "--version"))) {
        ARGS_VersionInfo();
      }
      /* === DEBUG OPTIONS (should only affect debug builds) === */
      elif (STR_Equals(argv[i], (flag = "--debug"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          debugger->is_debugging = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--debug-viz"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          debugger->is_viz = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--debug-arg-dump"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->dbg_arg_dump = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === PIPELINE OPTIONS === */
      elif (STR_Equals(argv[i], (flag = "--verbose"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->verbose_level = atoi(argv[i]);
          if (args->verbose_level < 0 || args->verbose_level > NUM_VERBOSITY_MODES - 1) {
            fprintf(stderr, "ERROR: Verbose Level (%d) is outside acceptable range (%d,%d).\n",
                    args->verbose_level, 0, NUM_VERBOSITY_MODES - 1);
            args->verbose_level = MAX(args->verbose_level, 0);
            args->verbose_level = MIN(args->verbose_level, NUM_VERBOSITY_MODES - 1);
            fprintf(stderr, "WARNING: Verbose Level set to: %d.\n", args->verbose_level);
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--num-threads"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->num_threads = atoi(argv[i]);
          args->num_threads = atoi(argv[i]);
          if (args->num_threads < 1) {
            fprintf(stderr, "ERROR: Number of Threads (%d) is outside acceptable range (%d,%d).\n",
                    args->num_threads, 1, 16);
            args->num_threads = MAX(args->num_threads, 1);
            args->num_threads = MIN(args->num_threads, 16);
            fprintf(stderr, "WARNING: Number of Threads set to: %d.\n", args->num_threads);
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--enforce-errors"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          if (atoi(argv[i])) {
            args->enforce_warnings = false;
          }
          elif (atoi(argv[i]) == 1) {
            args->enforce_warnings = true;
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--enforce-warnings"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          if (atoi(argv[i])) {
            args->enforce_warnings = false;
          }
          elif (atoi(argv[i]) == 1) {
            args->enforce_warnings = true;
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--search-type"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          if (STR_Equals(argv[i], "P2S") == true) {
            args->search_name = STR_Set(args->search_name, argv[i]);
            args->target_prep_type = FILE_MSA;
            args->query_prep_type = FILE_FASTA;
          }
          elif (STR_Equals(argv[i], "S2S") == true) {
            args->search_name = STR_Set(args->search_name, argv[i]);
            args->target_prep_type = FILE_FASTA;
            args->query_prep_type = FILE_FASTA;
          }
          else {
            fprintf(stderr, "ERROR: '%s' is not a valid argument for --search-type.\n", argv[i]);
            ERRORCHECK_exit(EXIT_FAILURE);
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === TASKS === */
      elif (STR_Equals(argv[i], (flag = "--run-mmseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_mmseqs = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-mmseqs-pref"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_mmseqs_pref = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-mmseqs-align"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_mmseqs_align = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-convert"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_convert = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-mmore"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_mmore = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--use-pvals"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->use_pvals = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === INPUT PROGRAMS === */
      elif (STR_Equals(argv[i], (flag = "--program-mmseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_program = STR_Set(args->mmseqs_program, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--program-hmmer"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->hmmer_program = STR_Set(args->hmmer_program, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--program-mmoreseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmoreseqs_program = STR_Set(args->mmoreseqs_program, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--script-dir"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmoreseqs_scripts = STR_Set(args->mmoreseqs_scripts, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === INPUT FILES === */
      elif (STR_Equals(argv[i], (flag = "--index"))) {
        req_args = 2;
        if (i + req_args < argc) {
          i++;
          args->t_index_filein = STR_Set(args->t_index_filein, argv[i]);
          i++;
          args->q_index_filein = STR_Set(args->q_index_filein, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--local-tools"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_use_local_tools = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--guess-ftype"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_guess_filetype = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmoreseqs-ftype"))) {
        req_args = 5;
        if (i + req_args < argc) {
          i++;
          args->t_filetype = atoi(argv[i]);
          i++;
          args->q_filetype = atoi(argv[i]);
          i++;
          args->t_mmseqs_p_filetype = atoi(argv[i]);
          i++;
          args->t_mmseqs_s_filetype = atoi(argv[i]);
          i++;
          args->q_mmseqs_filetype = atoi(argv[i]);
          /* since filetype supplied, turn off guesser */
          args->is_guess_filetype = false;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmoreseqs-main-ftype"))) {
        req_args = 2;
        if (i + req_args < argc) {
          i++;
          args->t_filetype = atoi(argv[i]);
          i++;
          args->q_filetype = atoi(argv[i]);
          /* since filetype supplied, turn off guesser */
          args->is_guess_filetype = false;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--tmp"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[i]);
          args->prep_folderpath = STR_Set(args->prep_folderpath, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--prep"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[i]);
          args->prep_folderpath = STR_Set(args->prep_folderpath, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-m8"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_m8_filein = STR_Set(args->mmseqs_m8_filein, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === PREP FILES === */
      elif (STR_Equals(argv[i], (flag = "--prep-copy"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_prep_copy = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--prep-link-target-mmore"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->link_target_mmore_prep = STR_Set(args->link_target_mmore_prep, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--prep-link-query-mmore"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->link_query_mmore_prep = STR_Set(args->link_query_mmore_prep, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--prep-link-target-mmseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->link_target_mmseqs_prep = STR_Set(args->link_target_mmseqs_prep, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--prep-link-query-mmseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->link_query_mmseqs_prep = STR_Set(args->link_query_mmseqs_prep, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === MMORE METADATA === */
      elif (STR_Equals(argv[i], (flag = "--dbsizes"))) {
        req_args = 2;
        if (i + req_args < argc) {
          i++;
          args->q_dbsize = atoi(argv[i]);
          i++;
          args->t_dbsize = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: '%s' flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === MMORE PARAMETERS === */
      elif (STR_Equals(argv[i], (flag = "--alpha"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->alpha = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--beta"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->beta = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--gamma"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->gamma = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--hard-limit"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->hard_limit = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* ==== MMORE OPTIONS === */
      elif (STR_Equals(argv[i], (flag = "--run-prep"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-bias"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_bias = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-full"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_full = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-domains"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_domains = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-vit-mmore"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_vit_mmore = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-mmseqsaln"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_mmseqsaln = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-vitaln"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_vit = atoi(argv[i]);
          args->is_run_vitaln = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-vit"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_vit = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-postaln"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_optacc = atoi(argv[i]);
          args->is_run_postaln = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-post"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_optacc = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === MMORE FILTERS === */
      elif (STR_Equals(argv[i], (flag = "--run-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_filter = atoi(argv[i]);
          args->is_run_viterbi_filter = atoi(argv[i]);
          args->is_run_cloud_filter = atoi(argv[i]);
          args->is_run_boundfwd_filter = atoi(argv[i]);
          args->is_run_fwdback_filter = atoi(argv[i]);
          args->is_run_report_filter = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-vit-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_viterbi_filter = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-cld-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_cloud_filter = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-fwd-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_boundfwd_filter = atoi(argv[i]);
          args->is_run_fwdback_filter = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--vit-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->threshold_vit = atof(argv[i]);
          // args->is_run_viterbi_filter   = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--cld-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->threshold_cloud = atof(argv[i]);
          // args->is_run_cloud_filter  = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--fwd-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->threshold_boundfwd = atof(argv[i]);
          args->threshold_fwd = atof(argv[i]);
          // args->is_run_boundfwd_filter  = true;
          // args->is_run_fwdback_filter   = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--eval"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->threshold_report_eval = atof(argv[i]);
          args->mmore_evalue = atof(argv[i]);
          args->is_run_report_filter = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === MMSEQS PARAMETERS === */
      elif (STR_Equals(argv[i], (flag = "--mmseqs-split"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_split = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-kmer"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_kmer = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-kscore"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_kscore = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-sens"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_sensitivity = atof(argv[i]);
          args->is_run_mmseqs_sens = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-ungapped-score"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_ungapped_vit = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-eval"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_evalue = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-pval"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_pvalue = atof(argv[i]);
          args->mmseqs_p2s_pvalue = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-maxseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_maxhits = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-altalis"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_altalis = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === MMSEQS METADATA === */
      elif (STR_Equals(argv[i], (flag = "--mmseqs-times"))) {
        req_args = 5;
        if (i + req_args < argc) {
          i++;
          args->prep_time = atof(argv[i]);
          i++;
          args->mmseqs_time = atof(argv[i]);
          i++;
          args->mmseqs_prefilter_time = atof(argv[i]);
          i++;
          args->mmseqs_p2s_time = atof(argv[i]);
          i++;
          args->mmseqs_s2s_time = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-dbsizes"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->prefilter_dbsize = atoi(argv[i]);
          i++;
          args->mmseqs_dbsize = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === SEARCH/RANGE OPTIONS === */
      elif (STR_Equals(argv[i], (flag = "--range"))) {
        req_args = 2;
        if (i + req_args < argc) {
          i++;
          args->list_range.beg = atoi(argv[i]);
          i++;
          args->list_range.end = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: '%s' flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--search-mode"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->search_mode = atoi(argv[i]);
          if (args->search_mode < 0 || args->search_mode > NUM_SELECT_SEARCHES) {
            fprintf(stderr, "ERROR: Verbose level (%d) is outside acceptable range (%d,%d).\n",
                    args->search_mode, 0, NUM_SELECT_SEARCHES);
            args->search_mode = MAX(args->search_mode, 0);
            args->search_mode = MIN(args->search_mode, NUM_SELECT_SEARCHES - 1);
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === INTERRIM OUTPUT === */
      elif (STR_Equals(argv[i], (flag = "--mmseqs-m8out"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          args->mmseqs_m8out = STR_Set(args->mmseqs_m8out, argv[i]);
          // args->mmseqs_m8_filein = STR_Set(args->mmseqs_m8_filein, argv[i]);
          args->is_mmseqs_m8out = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === OUTPUT === */
      elif (STR_Equals(argv[i], (flag = "--stderr"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          fclose(stdout);
          free(args->stderr_fileout);
          args->stderr_fileout = STR_Create(argv[i]);
          args->is_redirect_stderr = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--stdout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          free(args->stdout_fileout);
          args->stdout_fileout = STR_Create(argv[i]);
          args->is_redirect_stdout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--allout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          STR_Destroy(args->hmmerout_fileout);
          args->hmmerout_fileout = STR_Concat(argv[i], ".hmmerout");
          STR_Destroy(args->m8out_fileout);
          args->m8out_fileout = STR_Concat(argv[i], ".m8out");
          args->is_m8out = true;
          STR_Destroy(args->myout_fileout);
          args->myout_fileout = STR_Concat(argv[i], ".myout");
          args->is_myout = true;
          STR_Destroy(args->mydom_fileout);
          args->mydom_fileout = STR_Concat(argv[i], ".mydomout");
          args->is_mydom = true;
          STR_Destroy(args->mytime_fileout);
          args->mytime_fileout = STR_Concat(argv[i], ".mytimeout");
          args->is_mytimeout = true;
          STR_Destroy(args->mythresh_fileout);
          args->mythresh_fileout = STR_Concat(argv[i], ".mythreshout");
          args->is_mythreshout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--domtblout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->mydom_fileout);
          args->mydom_fileout = STR_Create(argv[i]);
          args->is_mydom = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--hmmerout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->hmmerout_fileout);
          args->hmmerout_fileout = STR_Create(argv[i]);
          args->is_hmmerout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--m8out"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          STR_Destroy(args->m8out_fileout);
          args->m8out_fileout = STR_Create(argv[i]);
          args->is_m8out = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--myout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->myout_fileout);
          args->myout_fileout = STR_Create(argv[i]);
          args->is_myout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mydomtblout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->mydom_fileout);
          args->mydom_fileout = STR_Create(argv[i]);
          args->is_mydom = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mytimeout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->mytime_fileout);
          args->mytime_fileout = STR_Create(argv[i]);
          args->is_mytimeout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mythreshout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->mythresh_fileout);
          args->mythresh_fileout = STR_Create(argv[i]);
          args->is_mythreshout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--customout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->hmmerout_fileout);
          args->hmmerout_fileout = STR_Create(argv[i]);
          args->is_hmmerout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--debugout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->dbg_folderpath);
          args->dbg_folderpath = STR_Create(argv[i]);
          args->is_debug = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === NOT PROPER FLAG === */
      else {
        fprintf(stderr, "ERROR: '%s' is not a recognized flag.\n", argv[i]);
        ERRORCHECK_exit(EXIT_FAILURE);
      }
    } else {
      fprintf(stderr, "ERROR: '%s' is not associated with an argument or flag\n", argv[i]);
      ERRORCHECK_exit(EXIT_FAILURE);
    }
  }
}

STATUS_FLAG ARGS_ParseOptionsArgp(ARGS* args,
                                  int argc,
                                  char* argv[]) {
  /* add all commandline options */
  static struct argp_option options[] = {
    /* General options */
    {"verbose", VERBOSE_KEY, NULL, OPTION_ARG_OPTIONAL, 
        "The amount of output.." },
    {"version", VERSION_KEY, NULL, OPTION_ARG_OPTIONAL, 
        "Produce version information."},
    {"num-threads", NUM_THREADS_KEY, NULL, OPTION_ARG_OPTIONAL, 
        "The number of parallel threads to run."},
    {"eval", EVAL_KEY, NULL, OPTION_ARG_OPTIONAL, 
        "Set E-value filter threshold cutoff score for reporting."},
    {"use-pvals", USE_PVALS_KEY, NULL, OPTION_ARG_OPTIONAL, 
        "Whether to use P-value (as opposed to using the default E-value) for assessing reporting and filtering thresholds."},
    /* Pipeline options */
    {"run-prep", RUN_PREP_KEY, NULL, OPTION_ARG_OPTIONAL, 
        "Run file preparation stage of pipeline."},
    {"run-mmseqs", RUN_MMSEQS_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-mmseqs-pref", RUN_MMSEQS_PREF_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-mmseqs-align", RUN_MMSEQS_ALIGN_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-convert", RUN_CONVERT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-mmore", RUN_MMORE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-vit", RUN_VIT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-post", RUN_POST_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    /* Input file options */
    {"prep", PREP_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"tmp", TMP_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-m8", MMSEQS_M8_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"index", INDEX_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    /* MMSeqs options */
    {"mmseqs-kmer", MMSEQS_KMER_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-kscore", MMSEQS_KSCORE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-sens", MMSEQS_SENS_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-ungapped-vit", MMSEQS_UNGAPPED_VIT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-eval", MMSEQS_EVAL_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-pval", MMSEQS_PVAL_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-hits-per-search", MMSEQS_HITS_PER_SEARCH_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-altalis", MMSEQS_ALTALIS_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-split", MMSEQS_SPLIT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    /* MMORE options */
    {"vit-filter", VIT_FILTER_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"cld-filter", CLD_FILTER_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"fwd-filter", FWD_FILTER_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"alpha", ALPHA_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"beta", BETA_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"gamma", GAMMA_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"hard-limit", HARD_LIMIT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"range", RANGE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    /* Output file options */
    {"stdout", STDOUT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"stderr", STDERR_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"m8out", M8OUT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"domtblout", DOMTBLOUT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-m8out", MMSEQS_M8OUT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    /* Uncommon options */
    {"run-domains", RUN_DOMAINS_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-bias", RUN_BIAS_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-filter", RUN_FILTER_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-vit-filter", RUN_VIT_FILTER_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-cld-filter", RUN_CLD_FILTER_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-fwd-filter", RUN_FWD_FILTER_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-vit-mmore", RUN_VIT_MMORE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-full", RUN_FULL_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-mmseqsaln", RUN_MMSEQSALN_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-vitaln", RUN_VITALN_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"run-postaln", RUN_POSTALN_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    /* Prep options */
    {"prep-copy", PREP_COPY_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"prep-link-target-mmore", PREP_LINK_TARGET_MMORE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"prep-link-query-mmore", PREP_LINK_QUERY_MMORE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"prep-link-target-mmseqs", PREP_LINK_TARGET_MMSEQS_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"prep-link-query-mmore", PREP_LINK_QUERY_MMSEQS_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    /* Forwarding options */
    {"search-mode", SEARCH_MODE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"search-type", SEARCH_TYPE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"guess-ftype", GUESS_FTYPE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmoreseqs-ftype", MMORESEQS_FTYPE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmoreseqs-main-ftype", MMORESEQS_MAIN_FTYPE_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"program-mmseqs", PROGRAM_MMSEQS_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"program-hmmer", PROGRAM_HMMER_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"program-mmoreseqs", PROGRAM_MMORESEQS_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"local-tools", LOCAL_TOOLS_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"script-dir", SCRIPT_DIR_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"dbsizes", DBSIZES_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-times", MMSEQS_TIMES_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mmseqs-dbsizes", MMSEQS_DBSIZES_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    /* Developer options */
    {"debug", DEBUG_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"debug-viz", DEBUG_VIZ_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"allout", ALLOUT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"myout", MYOUT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mytimeout", MYTIMEOUT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"mythreshout", MYTHRESHOUT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"customout", CUSTOMOUT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    {"debugout", DEBUGOUT_KEY, NULL, OPTION_ARG_OPTIONAL, ""},
    /* terminator option */
    { 0 }
  };
  
  // struct argp argp = {options, parse_opt, 0, 0};
  // return argp_parse(&argp, argc, argv, 0, 0, &args);
}

static error_t
parse_opt(int key, char *arg, struct argp_state *state) {
  ARGS* args;
  int req_args;

  switch(key) {
    /* General options */
    case VERBOSE_KEY:
      args->verbose_level = atoi(arg);
      break;
    case VERSION_KEY:
      ARGS_VersionInfo();
      break;
    case NUM_THREADS_KEY:
      args->num_threads = atoi(arg);
      break;
    case EVAL_KEY:
      break;
    case USE_PVALS_KEY:
      break;
    /* Pipeline options */
    case RUN_PREP_KEY:
      break;
    case RUN_MMSEQS_KEY:
      break;
    case RUN_MMSEQS_PREF_KEY:
      break;
    case RUN_MMSEQS_ALIGN_KEY:
      break;
    case RUN_CONVERT_KEY:
      break;
    case RUN_MMORE_KEY:
      break;
    case RUN_VIT_KEY:
      break;
    case RUN_POST_KEY:
      break;
    /* Input file options */
    case PREP_KEY:
      break;
    case TMP_KEY:
      break;
    case MMSEQS_M8_KEY:
      break;
    case INDEX_KEY:
      break;
    /* MMSeqs options */
    case MMSEQS_KMER_KEY:
      break;
    case MMSEQS_KSCORE_KEY:
      break;
    case MMSEQS_SENS_KEY:
      break;
    case MMSEQS_UNGAPPED_VIT_KEY:
      break;
    case MMSEQS_EVAL_KEY:
      break;
    case MMSEQS_PVAL_KEY:
      break;
    case MMSEQS_HITS_PER_SEARCH_KEY:
      break;
    case MMSEQS_ALTALIS_KEY:
      break;
    case MMSEQS_SPLIT_KEY:
      break;
    /* MMORE options */
    case VIT_FILTER_KEY:
      break;
    case CLD_FILTER_KEY:
      break;
    case FWD_FILTER_KEY:
      break;
    case ALPHA_KEY:
      break;
    case BETA_KEY:
      break;
    case GAMMA_KEY:
      break;
    case HARD_LIMIT_KEY:
      break;
    case RANGE_KEY:
      break;
    /* Output file options */
    case STDOUT_KEY:
      break;
    case STDERR_KEY:
      break;
    case M8OUT_KEY:
      break;
    case DOMTBLOUT_KEY:
      break;
    case MMSEQS_M8OUT_KEY:
      break;
    /* */
    case RUN_DOMAINS_KEY:
      break;
    case RUN_BIAS_KEY:
      break;
    case RUN_FILTER_KEY:
      break;
    case RUN_VIT_FILTER_KEY:
      break;
    case RUN_CLD_FILTER_KEY:
      break;
    case RUN_FWD_FILTER_KEY:
      break;
    case RUN_VIT_MMORE_KEY:
      break;
    case RUN_FULL_KEY:
      break;
    case RUN_MMSEQSALN_KEY:
      break;
    case RUN_VITALN_KEY:
      break;
    case RUN_POSTALN_KEY:
      break;
    /* */
    case PREP_COPY_KEY:
      break;
    case PREP_LINK_TARGET_MMORE_KEY:
      break;
    case PREP_LINK_QUERY_MMORE_KEY:
      break;
    case PREP_LINK_TARGET_MMSEQS_KEY:
      break;
    case PREP_LINK_QUERY_MMSEQS_KEY:
      break;
    /* */
    case SEARCH_MODE_KEY:
      break;
    case SEARCH_TYPE_KEY:
      break;
    case GUESS_FTYPE_KEY:
      break;
    case MMORESEQS_FTYPE_KEY:
      break;
    case MMORESEQS_MAIN_FTYPE_KEY:
      break;
    case PROGRAM_MMSEQS_KEY:
      break;
    case PROGRAM_HMMER_KEY:
      break;
    case PROGRAM_MMORESEQS_KEY:
      break;
    case LOCAL_TOOLS_KEY:
      break;
    case SCRIPT_DIR_KEY:
      break;
    case DBSIZES_KEY:
      break;
    case MMSEQS_TIMES_KEY:
      break;
    case MMSEQS_DBSIZES_KEY:
      break;
    /* Developer options */
    case DEBUG_KEY:
      break;
    case DEBUG_VIZ_KEY:
      break;
    case ALLOUT_KEY:
      break;
    case MYOUT_KEY:
      break;
    case MYTIMEOUT_KEY:
      break;
    case MYTHRESHOUT_KEY:
      break;
    case CUSTOMOUT_KEY:
      break;
    case DEBUGOUT_KEY:
      break;
    /* Errors: too many or too few args,  */
    case ARGP_KEY_ARG:
      if (state->arg_num >= 2) {
        argp_usage(state);
      }
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 2) {
        argp_usage(state);
      }
      break;
    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}
