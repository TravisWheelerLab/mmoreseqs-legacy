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

   /* ENVIRONMENTAL ARGS */
   /* pass main args without type appended */
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMORE",         args->t_filepath );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_MMORE",          args->q_filepath );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_P",      args->t_mmseqs_p_filepath );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_S",      args->t_mmseqs_s_filepath );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_MMSEQS",         args->q_mmseqs_filepath );
   /* MAIN ARG TYPES */
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMORE_TYPE",    FILE_TYPE_NAMES[args->t_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_MMORE_TYPE",     FILE_TYPE_NAMES[args->q_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_P_TYPE", FILE_TYPE_NAMES[args->t_mmseqs_p_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_S_TYPE", FILE_TYPE_NAMES[args->t_mmseqs_s_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_MMSEQS_TYPE",    FILE_TYPE_NAMES[args->q_mmseqs_filetype] );
   /* PREPARATION FILES */
   SCRIPTRUNNER_Add_Env_Variable( runner, "PREP_DIR",             args->prep_folder );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_PREP",          args->target_prep );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_PREP",           args->target_prep );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_PREP_TYPE",     args->target_prep );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_PREP",           args->target_prep );
   /* INTERRIM FILES (FILE CREATED DURING PROCESS) */
   
   SCRIPTRUNNER_Add_Env_Variable( runner, "TEMP_DIR",             args->tmp_folderpath );
   
   /* OUTPUT FILES */
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "OUTPUT_PATH",       args->output_filepath,                       args->is_redirect_stdout );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "ERROR_PATH",        args->output_filepath,                       args->is_redirect_stderr );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "M8OUT_PATH",        args->m8out_filepath,                        args->is_m8out );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MYOUT_PATH",        args->myout_filepath,                        args->is_myout );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MYDOMOUT_PATH",     args->mydomout_filepath,                     args->is_mydomout );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MYTIMEOUT_PATH",    args->mytimeout_filepath,                    args->is_mytimeout );
   SCRIPTRUNNER_If_Add_Env_Variable( runner, "MYTHRESHOUT_PATH",  args->mythreshout_filepath,                  args->is_mythreshout );
   /* MMSEQS PARAMETERS */
   /* PREFILTER */
   /* P2S SEARCH */
   /* S2S SEARCH */
   SCRIPTRUNNER_Add_Env_Variable( runner, "HITS_PER_SEARCH",      INT_To_String( args->mmseqs_hits_per_search, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "KMER",                 INT_To_String( args->mmseqs_kmer,            buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "PREFILTER_THRESH",     INT_To_String( args->mmseqs_prefilter,       buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "UNGAPPEDVIT_THRESH",   INT_To_String( args->mmseqs_ungapped_vit,    buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "GAPPEDVIT_THRESH",     INT_To_String( args->mmseqs_evalue,          buffer ) );
   /* MMORE PARAMETERS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "ALPHA",                FLT_To_String( args->alpha,                  buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "BETA",                 FLT_To_String( args->beta,                   buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "GAMMA",                INT_To_String( args->gamma,                  buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "VITERBI_THRESH",       FLT_To_String( args->threshold_vit,          buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "CLOUD_THRESH",         FLT_To_String( args->threshold_cloud,        buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "BOUNDFWD_THRESH",      FLT_To_String( args->threshold_bound_fwd,    buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "FWDBACK_THRESH",       FLT_To_String( args->threshold_fwd,          buffer ) );
   /* TASK OPTIONS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_PREP",              INT_To_String( args->is_run_prep,            buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_FILTER",            INT_To_String( args->is_run_filter,          buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_BIAS",              INT_To_String( args->is_run_bias,            buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_DOMAIN",            INT_To_String( args->is_run_domains,         buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_FULL",              INT_To_String( args->is_run_full,            buffer ) );
   /* OPTIONS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "VERBOSE",              INT_To_String( args->verbose_level,          buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "RM_TEMP",              INT_To_String( args->tmp_remove,             buffer ) );
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
   // SCRIPTRUNNER_Add_Env_Variable( runner, "USE_LOCAL_TOOLS",   INT_To_String( args->is_use_local_tools, buffer ) );
}