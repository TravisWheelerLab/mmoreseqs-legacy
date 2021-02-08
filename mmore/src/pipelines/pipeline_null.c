/*******************************************************************************
 *  FILE:      pipeline_null.c
 *  PURPOSE:   Empty Pipeline
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
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
null_pipeline( WORKER* worker )
{
   printf("=== NULL PIPELINE ===\n");

   char* command[] = { MACRO_XSTR(PROJECT_LOCATION) "./scripts/workflows/mmseqs_plus_easy.sh", "target.hmm", "query.fasta", "", NULL };

   printf("# EXECUTING: %s\n", command[0] );
   printf("# MMSEQS_PLUS_SCRIPT: %s\n", MMSEQS_PLUS_SCRIPT );
   printf("# FASTA_TO_HMM_SCRIPT: %s\n", FASTA_TO_HMM_SCRIPT );
   int exit_code = execvp( command[0], command );

   return STATUS_SUCCESS;
}