/*******************************************************************************
 *  FILE:      work_loop.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for inner search loop.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>

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
#include "work_scripting.h"

/*! FUNCTION:  	WORK_load_script_env_args()
 *  SYNOPSIS:  	Sets all environmental variables for values from <args>.
 */
void 
WORK_load_script_env_args( WORKER* worker )
{
   ARGS*          args     = worker->args;
   SCRIPTRUNNER*  runner   = worker->runner;
   char           buffer[256];

   /* --- ENVIRONMENTAL ARGS --- */
   /* MMORE ARGS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMORE",         args->t_filein );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_MMORE",          args->q_filein );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMORE_TYPE",    FILETYPE_NAME_Get( args->t_filetype ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_MMORE_TYPE",     FILETYPE_NAME_Get( args->q_filetype ) );
   /* MMSEQS ARGS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_P",      args->t_mmseqs_p_filein );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_S",      args->t_mmseqs_s_filein );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_MMSEQS",         args->q_mmseqs_filein );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_P_TYPE", FILETYPE_NAME_Get( args->t_mmseqs_p_filetype ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_S_TYPE", FILETYPE_NAME_Get( args->t_mmseqs_s_filetype ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_MMSEQS_TYPE",    FILETYPE_NAME_Get( args->q_mmseqs_filetype ) );
   /* PREP ARGS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "PREP_DIR",             args->prep_folderpath );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET",               args->target_prep );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY",                args->query_prep );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_TYPE",          FILETYPE_NAME_Get( args->target_prep_type ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_TYPE",           FILETYPE_NAME_Get( args->query_prep_type ) );
   /* OTHER ARGS */
   /* PIPELINE OPTIONS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "VERBOSE",              INT_ToString( args->verbose_level,          buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "NUM_THREADS",          INT_ToString( args->num_threads,            buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_RM_TEMP",           INT_ToString( args->tmp_remove,             buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_PREP",              INT_ToString( args->is_run_prep,            buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_COPY",              INT_ToString( args->is_prep_copy,           buffer ) );
   // SCRIPTRUNNER_Add_Env_Variable( runner, "DO_OVERWRITE",         INT_ToString( args->is_overwrite,           buffer ) );
   // SCRIPTRUNNER_Add_Env_Variable( runner, "DO_IGNORE_WARNINGS",   INT_ToString( args->is_ignore_warnings,     buffer ) );
   /* SEARCH OPTIONS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "SEARCH_TYPE",          args->search_name );
   // SCRIPTRUNNER_Add_Env_Variable( runner, "SEARCH_TYPE",          INT_ToString( args->search_type,            buffer ) );
   
   /* MMSEQS OPTIONS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_DO_UNGAPPED",   INT_ToString( args->is_run_mmseqs_ungapped, buffer ) );
   /* MMSEQS PARAMETERS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_MAXSEQS",       INT_ToString( args->mmseqs_maxhits,         buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_ALTALIS",       INT_ToString( args->mmseqs_altalis,         buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_KMER",          INT_ToString( args->mmseqs_kmer,            buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_KSCORE",        INT_ToString( args->mmseqs_prefilter,       buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_UNGAPPED",      INT_ToString( args->mmseqs_ungapped_vit,    buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_EVAL",          INT_ToString( args->mmseqs_evalue,          buffer ) );
   // SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_SENS",          INT_ToString( args->mmseqs_sensitivity,     buffer ) );
   /* MMORE PARAMETERS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_ALPHA",          FLT_ToString( args->alpha,                  buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_BETA",           FLT_ToString( args->beta,                   buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_GAMMA",          INT_ToString( args->gamma,                  buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_VITERBI",        FLT_ToString( args->threshold_vit,          buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_CLOUD",          FLT_ToString( args->threshold_cloud,        buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_BOUNDFWD",       FLT_ToString( args->threshold_boundfwd,    buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_FWDBACK",        FLT_ToString( args->threshold_fwd,          buffer ) );
   /* MMORE OPTIONS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_DO_MMORE",       INT_ToString( args->is_run_mmore,           buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_DO_FILTER",      INT_ToString( args->is_run_filter,          buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_DO_BIAS",        INT_ToString( args->is_run_bias,            buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_DO_DOMAIN",      INT_ToString( args->is_run_domains,         buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_DO_FULL",        INT_ToString( args->is_run_full,            buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_DO_ALLOUT",      INT_ToString( args->is_run_full,            buffer ) );
   /* INTERRIM FILES (FILE CREATED DURING PROCESS) */
   // SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_S_OUT",             args->t_mmseqs_s_fileout );
   // SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_P_OUT",             args->t_mmseqs_s_fileout );
   /* OUTPUT FILES */
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MMORE_SD2DOUT_OUT",        args->stdout_fileout,            args->is_redirect_stdout );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MMORE_SD2DERR_OUT",        args->stderr_fileout,            args->is_redirect_stderr );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MMORE_M8OUT_OUT",         args->m8out_fileout,             args->is_m8out );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MMORE_MYOUT_OUT",         args->myout_fileout,             args->is_myout );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MMORE_MYDOM_OUT",         args->mydom_fileout,             args->is_mydom );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MMORE_MYTIME_OUT",        args->mytime_fileout,            args->is_mytimeout );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MMORE_MYTHRESH_OUT",      args->mythresh_fileout,          args->is_mythreshout );
}

/*! FUNCTION:  	WORK_load_script_tools()
 *  SYNOPSIS:  	Sets all environmental variables for tool locations (mmseqs, hmmer, esl, python, ...).
 */
void 
WORK_load_script_tools( WORKER* worker )
{
   // ARGS*          args     = worker->args;
   // SCRIPTRUNNER*  runner   = worker->runner;

   // /* TOOLS */
   // SCRIPTRUNNER_Add_Env_Variable( runner, "ROOT_DIR",          project_path );
   // SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_DIR",         mmore_path );
   // SCRIPTRUNNER_Add_Env_Variable( runner, "HMMER_DIR",         hmmer_path );
   // SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_DIR",        mmseqs_path ); 
   // SCRIPTRUNNER_Add_Env_Variable( runner, "USE_LOCAL_TOOLS",   INT_ToString( args->is_use_local_tools, buffer ) );
}