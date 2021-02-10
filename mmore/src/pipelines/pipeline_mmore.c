/*******************************************************************************
 *  FILE:      pipeline_mmore.c
 *  PURPOSE:   Cloud Search Pipeline for MMORE pipeline.
 *             Set arguments and issues command to bash script.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>

/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../work/_work.h"

/* header */
#include "_pipelines.h"

/*! FUNCTION:  mmore_pipeline()
 *  SYNOPSIS:  Overarching MMORE pipeline. 
 *             This pipeline runs the full 
 */
STATUS_FLAG 
mmore_pipeline( WORKER* worker )
{
   printf("=== MMORE PIPELINE ===\n");
   
   ARGS*          args     = worker->args;
   TASKS*         tasks    = worker->tasks;
   SCRIPTRUNNER*  runner   = worker->runner;

   /* get location of tools */
   STR project_path     = ROOT_DIR;
   // printf("# PROJECT_LOCATION: %s\n", project_path );
   STR mmseqs_path      = MMSEQS_BIN;
   // printf("# MMSEQS_LOCATION: %s\n", mmseqs_path );
   STR hmmer_path       = HMMER_BIN;
   // printf("# HMMER_LOCATION: %s\n", mmseqs_path );
   STR mmore_path       = MMORE_BIN;
   // printf("# MMORE_LOCATION: %s\n", mmore_path );
   STR script_relpath   = "/scripts/workflows/mmore.sh";
   STR script_path      = STR_Concat( project_path, script_relpath );
   // printf("# SCRIPT_LOCATION: %s\n\n", script_path );

   /* buffer for converting args to string */
   STR   str = NULL;
   char  buffer[256];

   /* SCRIPT */
   SCRIPTRUNNER_SetScript( runner, script_path );
   SCRIPTRUNNER_Add_Script_Argument( runner, NULL, args->t_filepath );
   SCRIPTRUNNER_Add_Script_Argument( runner, NULL, args->q_filepath );
   SCRIPTRUNNER_Add_Script_Argument( runner, NULL, args->t_mmseqs_filepath );
   /* TOOLS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "ROOT_DIR",          project_path );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_DIR",         mmore_path );
   SCRIPTRUNNER_Add_Env_Variable( runner, "HMMER_DIR",         hmmer_path );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_DIR",        mmseqs_path ); 
   SCRIPTRUNNER_Add_Env_Variable( runner, "USE_LOCAL_TOOLS",   INT_To_String( args->is_use_local_tools, buffer ) ); 
   /* MAIN ARGS */
   /* pass main args with type appended */
   str = STR_Set( str, "TARGET_");
   str = STR_Append( str, FILE_TYPE_NAMES[args->t_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner,    str,     args->t_filepath );
   str = STR_Set( str, "QUERY_");
   str = STR_Append( str, FILE_TYPE_NAMES[args->q_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner,    str,     args->q_filepath );
   str = STR_Set( str, "TARGET_MMSEQS_" );
   str = STR_Append( str, FILE_TYPE_NAMES[args->t_mmseqs_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner,    str,     args->t_mmseqs_filepath );
   /* pass main args without type appended */
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET",               FILE_TYPE_NAMES[args->t_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY",                FILE_TYPE_NAMES[args->q_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS",        FILE_TYPE_NAMES[args->t_mmseqs_filetype] );
   /* MAIN ARG TYPES */
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_TYPE",          FILE_TYPE_NAMES[args->t_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner, "QUERY_TYPE",           FILE_TYPE_NAMES[args->q_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner, "TARGET_MMSEQS_TYPE",   FILE_TYPE_NAMES[args->t_mmseqs_filetype] );
   /* INTERRIM FILES */
   if ( args->is_redirect_stdout ) {
      SCRIPTRUNNER_Add_Env_Variable( runner, "OUTPUT_PATH",  args->output_filepath );
   }
   /* OUTPUT FILES */
   if ( args->is_redirect_stdout ) {
      SCRIPTRUNNER_Add_Env_Variable( runner, "OUTPUT_PATH",  args->output_filepath );
   }
   if ( args->is_redirect_stderr ) {
      SCRIPTRUNNER_Add_Env_Variable( runner, "ERROR_PATH",  args->output_filepath );
   }
   if ( args->is_m8out ) {
      SCRIPTRUNNER_Add_Env_Variable( runner, "M8OUT_PATH",  args->m8out_filepath );
   }
   if ( args->is_myout ) {
      SCRIPTRUNNER_Add_Env_Variable( runner, "MYOUT_PATH",  args->myout_filepath );
   }
   if ( args->is_mydomout ) {
      SCRIPTRUNNER_Add_Env_Variable( runner, "MYDOMOUT_PATH",  args->mydomout_filepath );
   }
   if ( args->is_mytimeout ) {
      SCRIPTRUNNER_Add_Env_Variable( runner, "MYTIMEOUT_PATH",  args->mytimeout_filepath );
   }
   if ( args->is_mythreshout ) {
      SCRIPTRUNNER_Add_Env_Variable( runner, "MYTHRESHOUT_PATH",  args->mythreshout_filepath );
   }
   /* MMSEQS PARAMETERS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "HITS_PER_SEARCH",      INT_To_String( args->mmseqs_hits_per_search, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "KMER",                 INT_To_String( args->mmseqs_kmer, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "PREFILTER_THRESH",     INT_To_String( args->mmseqs_prefilter, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "UNGAPPEDVIT_THRESH",   INT_To_String( args->mmseqs_ungapped_vit, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "GAPPEDVIT_THRESH",     INT_To_String( args->mmseqs_evalue, buffer ) );
   /* MMORE PARAMETERS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "ALPHA",                FLT_To_String( args->alpha, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "BETA",                 FLT_To_String( args->beta, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "GAMMA",                INT_To_String( args->gamma, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "VITERBI_THRESH",       FLT_To_String( args->threshold_vit, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "CLOUD_THRESH",         FLT_To_String( args->threshold_cloud, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "BOUNDFWD_THRESH",      FLT_To_String( args->threshold_bound_fwd, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "FWDBACK_THRESH",       FLT_To_String( args->threshold_fwd, buffer ) );
   /* TASK OPTIONS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_FILTER",   INT_To_String( args->is_run_filter, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_BIAS",     INT_To_String( args->is_compo_bias, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_DOMAIN",   INT_To_String( args->is_run_domains, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "DO_FULL",     INT_To_String( args->is_run_full, buffer ) );
   /* OPTIONS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "VERBOSE",     INT_To_String( args->verbose_level, buffer ) );
   SCRIPTRUNNER_Add_Env_Variable( runner, "RM_TEMP",     INT_To_String( args->tmp_remove, buffer ) );

   /* EXECUTE SCRIPT */
   SCRIPTRUNNER_Execute( runner );

   STR_Destroy(script_path);
   STR_Destroy(str);

   /* code should not reach here */
   fprintf( stderr, "# MMORE PIPELINE FAILED.\n" );
   exit(EXIT_FAILURE);
}

