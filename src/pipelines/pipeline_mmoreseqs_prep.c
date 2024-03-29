/*******************************************************************************
 *  - FILE:      pipeline_prep.c
 *  - DESC:    Cloud Search Pipeline for MMORE pipeline.
 *             Preparation pipeline, converts files to formats to be used by
 *             Set arguments and issues command to bash script.
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
mmoreseqs_prep_pipeline(WORKER* worker) {
  printf("=== MMORESEQS: PREP PIPELINE ===\n");

  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  SCRIPTRUNNER* runner = worker->runner;

  /* get location of tools */
  STR project_path = ROOT_DIR;
  // printf("# PROJECT_LOCATION: %s\n", project_path );
  STR mmseqs_path = MMSEQS_BIN;
  // printf("# MMSEQS_LOCATION: %s\n", mmseqs_path );
  STR hmmer_path = HMMER_BIN;
  // printf("# HMMER_LOCATION: %s\n", mmseqs_path );
  STR mmore_path = MMORE_BIN;
  // printf("# MMORE_LOCATION: %s\n", mmore_path );
  STR script_relpath = "/scripts/workflows/mmoreseqs-prep-prepare.sh";
  STR script_path = STR_Concat(project_path, script_relpath);
  printf_vhi("# SCRIPT_LOCATION: %s\n\n", script_path);

  /* buffer for converting args to string */
  STR str = NULL;
  char buffer[256];

  /* COMMANDLINE ARGS */
  SCRIPTRUNNER_SetScript(runner, script_path);
  SCRIPTRUNNER_Add_Script_Argument(runner, NULL, args->target_prep);
  SCRIPTRUNNER_Add_Script_Argument(runner, NULL, args->query_prep);
  SCRIPTRUNNER_Add_Script_Argument(runner, NULL, args->prep_folderpath);
  SCRIPTRUNNER_Add_Script_Argument(runner, NULL, FILETYPE_NAME_Get(args->target_prep_type));
  SCRIPTRUNNER_Add_Script_Argument(runner, NULL, FILETYPE_NAME_Get(args->query_prep_type));

  /* COMMANDLINE ENVIRONMENTAL VARIABLES */
  /* pass main args with type appended */
  str = STR_Concat("TARGET_PREP_", FILETYPE_NAME_Get(args->target_prep_type));
  SCRIPTRUNNER_Add_Env_Variable(runner, str, args->target_prep);
  str = STR_Concat("QUERY_PREP_", FILETYPE_NAME_Get(args->query_prep_type));
  SCRIPTRUNNER_Add_Env_Variable(runner, str, args->query_prep);

  /* TOOLS */
  SCRIPTRUNNER_Add_Env_Variable(runner, "ROOT_DIR", project_path);
  SCRIPTRUNNER_Add_Env_Variable(runner, "MMORE_DIR", mmore_path);
  SCRIPTRUNNER_Add_Env_Variable(runner, "HMMER_DIR", hmmer_path);
  SCRIPTRUNNER_Add_Env_Variable(runner, "MMSEQS_DIR", mmseqs_path);
  SCRIPTRUNNER_Add_Env_Variable(runner, "USE_LOCAL_TOOLS", INT_ToString(args->is_use_local_tools, buffer));
  SCRIPTRUNNER_Add_Env_Variable(runner, "LINK_QUERY_MMORE", args->link_query_mmore_prep);
  SCRIPTRUNNER_Add_Env_Variable(runner, "LINK_TARGET_MMORE", args->link_target_mmore_prep);
  SCRIPTRUNNER_Add_Env_Variable(runner, "LINK_QUERY_MMSEQS", args->link_query_mmseqs_prep);
  SCRIPTRUNNER_Add_Env_Variable(runner, "LINK_TARGET_MMSEQS", args->link_target_mmseqs_prep);

  /* ENVIRONMENTAL ARGS */
  WORK_load_script_env_args(worker);
  SCRIPTRUNNER_Add_Env_Variable(runner, "VERBOSE", INT_ToString(args->verbose_level, buffer));

  /* EXECUTE SCRIPT */
  SCRIPTRUNNER_Execute(runner);

  /* code only reaches here in the result of an error */
  fprintf(stderr, "# MMORE PIPELINE FAILED.\n");
  SCRIPTRUNNER_Dump(runner, stderr);
  STR_Destroy(script_path);
  STR_Destroy(str);
  ERRORCHECK_exit(EXIT_FAILURE);
}
