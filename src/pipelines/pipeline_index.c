/*******************************************************************************
 *  - FILE:      pipeline_index.c
 *  - DESC:    Index pipeline.
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

/*! FUNCTION:  index_pipeline()
 *  SYNOPSIS:  Index Pipeline: Indexes FASTA or HMM files.
 */
STATUS_FLAG
index_pipeline(WORKER* worker) {
  printf("=== INDEX PIPELINE ===\n");

  FILE* fp = NULL;
  ARGS* args = worker->args;
  TASKS* tasks = worker->tasks;
  CLOCK* timer = worker->timer;

  CLOCK_Start(timer);

  /* initialize data structures needed for tasks */
  WORK_init(worker);
  WORK_open(worker);

  /* building, sorting, and outputting target index */
  printf_vhi("# building target index from:\t%s\n",
             args->t_filein);
  WORK_build_target_index(worker);
  F_INDEX_Sort_by_Id(worker->t_index);
  printf_vhi("# outputting target index to:\t%s\n",
             args->t_index_filein);
  WORK_output_target_index(worker);

  /* building, sorting, and outputting query index */
  printf_vhi("# building query index from:\t%s\n",
             args->q_filein);
  WORK_build_query_index(worker);
  F_INDEX_Sort_by_Id(worker->q_index);
  printf_vhi("# outputting query index to:\t%s\n",
             args->q_index_filein);
  WORK_output_query_index(worker);

  /* clean up worker data structs */
  WORK_close(worker);
  WORK_cleanup(worker);
}
