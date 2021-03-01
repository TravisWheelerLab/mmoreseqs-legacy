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
#include "../parsers/_parsers.h"
#include "../algs_linear/_algs_linear.h"
#include "../algs_quad/_algs_quad.h"
#include "../algs_naive/_algs_naive.h"
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
   SCRIPTRUNNER_Add_Script_Argument( runner, NULL, args->t_mmseqs_p_filepath );
   SCRIPTRUNNER_Add_Script_Argument( runner, NULL, args->q_mmseqs_filepath );
   /* TOOLS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "ROOT_DIR",          project_path );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_DIR",         mmore_path );
   SCRIPTRUNNER_Add_Env_Variable( runner, "HMMER_DIR",         hmmer_path );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_DIR",        mmseqs_path ); 
   SCRIPTRUNNER_Add_Env_Variable( runner, "USE_LOCAL_TOOLS",   INT_To_String( args->is_use_local_tools, buffer ) ); 
   
   /* MAIN COMMANDLINE ARGS */
   /* pass main args with type appended */
   str = STR_Set( str, "TARGET_MMORE_");
   str = STR_Append( str, FILE_TYPE_NAMES[args->t_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner,    str,     args->t_filepath );
   str = STR_Set( str, "QUERY_MMORE_");
   str = STR_Append( str, FILE_TYPE_NAMES[args->q_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner,    str,     args->q_filepath );
   str = STR_Set( str, "TARGET_MMSEQS_P_" );
   str = STR_Append( str, FILE_TYPE_NAMES[args->t_mmseqs_p_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner,    str,     args->t_mmseqs_p_filepath );
   str = STR_Set( str, "TARGET_MMSEQS_S_" );
   str = STR_Append( str, FILE_TYPE_NAMES[args->t_mmseqs_s_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner,    str,     args->t_mmseqs_s_filepath );
   str = STR_Set( str, "QUERY_MMSEQS_" );
   str = STR_Append( str, FILE_TYPE_NAMES[args->q_mmseqs_filetype] );
   SCRIPTRUNNER_Add_Env_Variable( runner,    str,     args->q_mmseqs_filepath );

   /* ENVIRONMENTAL ARGS */
   WORK_load_script_env_args( worker );

   /* EXECUTE SCRIPT */
   SCRIPTRUNNER_Execute( runner );

   STR_Destroy(script_path);
   STR_Destroy(str);

   /* code should not reach here */
   fprintf( stderr, "# MMORE PIPELINE FAILED.\n" );
   ERRORCHECK_exit(EXIT_FAILURE);
}

