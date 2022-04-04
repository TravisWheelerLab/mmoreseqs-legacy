/*******************************************************************************
 *  - FILE:      pipeline_int_test.c
 *  - DESC:    Integrated Testing Pipeline.
 *             Requires -BUILD=DEBUG for full functionality.
 *  NOTES:
 *    - WIP.
 *    - Many unit tests still need to be built.
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

/*! FUNCTION:  itest_pipeline()
 *  SYNOPSIS:  Pipeline runs integration tests.
 *             Runs optimized and unoptimized versions of search algs and compares results.
 *             For full functionality, must be compiled in DEBUG mode.
 */
STATUS_FLAG
itest_pipeline(WORKER* worker) {
  printf("=== INTEGRATION TEST PIPELINE ===\n");
}
