/*******************************************************************************
 *  - FILE:  arg_parser.h
 *  - DESC:  Parses command line arguments.
 *    NOTES:
 *       - To be phased out in preference of cli using CLI11.
 *******************************************************************************/

/* imports */
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

  /* --- RANGE OPTIONS --- */
  args->t_range = (RANGE){-1, -1};
  args->q_range = (RANGE){-1, -1};
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
