/*******************************************************************************
 *  - FILE:     pipeline_main.c
 *  - DESC:   Main Cloud Search Pipeline.
 * 				    Runs Viterbi and Pruned Forward-Backward internally.
 *  NOTES:    Phasing out
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

/* standard pipeline */
STATUS_FLAG
generic_pipeline(WORKER* worker) {
}

/*! FUNCTION:
 *  SYNOPSIS:
 */
void main_pipeline_set_flags() {
}
