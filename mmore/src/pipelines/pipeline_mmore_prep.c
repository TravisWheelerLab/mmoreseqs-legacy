/*******************************************************************************
 *  FILE:      pipeline_prep.c
 *  PURPOSE:   Cloud Search Pipeline for MMORE pipeline.
 *             Preparation pipeline, converts files to formats to be used by 
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

/*! FUNCTION:  prep_pipeline()
 *  SYNOPSIS:  Helper for overarching MMORE-SEQS pipeline. 
 *             Prepares files for use in the MMORE-SEQS Pipeline
 */
STATUS_FLAG 
mmore_prep_pipeline( WORKER* worker )
{
   printf("=== MMORE PREP PIPELINE ===\n");
   
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
   STR script_relpath   = "/scripts/workflows/mmore-prep.sh";
   STR script_path      = STR_Concat( project_path, script_relpath );
   // printf("# SCRIPT_LOCATION: %s\n\n", script_path );

   /* buffer for converting args to string */
   STR   str = NULL;
   char  buffer[256];

   /* COMMANDLINE ARGS */
   SCRIPTRUNNER_SetScript( runner, script_path );
   SCRIPTRUNNER_Add_Script_Argument( runner, NULL, args->query_prep );
   SCRIPTRUNNER_Add_Script_Argument( runner, NULL, args->target_prep );
   SCRIPTRUNNER_Add_Script_Argument( runner, NULL, FILE_TYPE_NAMES[args->target_prep_type] );
   SCRIPTRUNNER_Add_Script_Argument( runner, NULL, FILE_TYPE_NAMES[args->query_prep_type] );

   
   /* COMMANDLINE ENVIRONMENTAL VARIABLES */
   /* pass main args with type appended */
   str = STR_Set( str, "TARGET_PREP_");
   str = STR_Append( str, FILE_TYPE_NAMES[args->target_prep_type] );
   SCRIPTRUNNER_Add_Env_Variable( runner,    str,     args->target_prep );
   str = STR_Set( str, "QUERY_PREP_");
   str = STR_Append( str, FILE_TYPE_NAMES[args->query_prep_type] );
   SCRIPTRUNNER_Add_Env_Variable( runner,    str,     args->query_prep );

   /* TOOLS */
   SCRIPTRUNNER_Add_Env_Variable( runner, "ROOT_DIR",          project_path );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMORE_DIR",         mmore_path );
   SCRIPTRUNNER_Add_Env_Variable( runner, "HMMER_DIR",         hmmer_path );
   SCRIPTRUNNER_Add_Env_Variable( runner, "MMSEQS_DIR",        mmseqs_path ); 
   SCRIPTRUNNER_Add_Env_Variable( runner, "USE_LOCAL_TOOLS",   INT_To_String( args->is_use_local_tools, buffer ) ); 
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

