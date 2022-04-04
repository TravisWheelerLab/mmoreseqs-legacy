/*******************************************************************************
 *  - FILE:      pipeline_null.c
 *  - DESC:    Empty Pipeline (for testing only).
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
#include "../work/_work.h"

/* header */
#include "_pipelines.h"

/*! FUNCTION:  null_pipeline()
 *  SYNOPSIS:  Pipeline does nothing.
 *             For testing.
 */
STATUS_FLAG
null_pipeline(WORKER* worker) {
  printf("=== NULL PIPELINE ===\n");

  SCRIPTRUNNER* runner = worker->runner;

  STR project_path = ROOT_DIR;
  printf("# PROJECT_LOCATION: %s\n", project_path);
  STR mmseqs_path = MMSEQS_BIN;
  printf("# MMSEQS_LOCATION: %s\n", mmseqs_path);
  STR hmmer_path = HMMER_BIN;
  printf("# HMMER_LOCATION: %s\n", hmmer_path);
  STR script_relpath = "/scripts/workflows/test_workflow.sh";
  STR script_path = STR_Concat(project_path, script_relpath);
  printf("# SCRIPT_LOCATION: %s\n\n", script_path);

  SCRIPTRUNNER_SetScript(runner, script_path);
  SCRIPTRUNNER_Add_Env_Variable(runner, "PROJECT_DIR", project_path);
  SCRIPTRUNNER_Add_Env_Variable(runner, "HMMER_DIR", hmmer_path);
  SCRIPTRUNNER_Add_Env_Variable(runner, "MMSEQS_DIR", mmseqs_path);
  SCRIPTRUNNER_Add_Script_Argument(runner, NULL, "target.hmm");
  SCRIPTRUNNER_Add_Script_Argument(runner, NULL, "query.fasta");
  SCRIPTRUNNER_Execute(runner);

  STR_Destroy(script_path);

  return STATUS_SUCCESS;
}
