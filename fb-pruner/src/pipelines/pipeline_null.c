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
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../parsers/parsers.h"
#include "../algs_linear/algs_linear.h"
#include "../algs_quad/algs_quad.h"
#include "../algs_naive/algs_naive.h"

/* header */
#include "pipelines.h"

/*  FUNCTION:  null_pipeline()
 *  SYNOPSIS:  Pipeline that does nothing.
 */
void 
null_pipeline( WORKER* worker )
{
   char* command[] = { MACRO_XSTR(PROJECT_LOCATION) "./scripts/workflows/mmseqs_plus_easy.sh", "target.hmm", "query.fasta", "", NULL };

   printf("# EXECUTING: %s\n", command[0] );
   printf("# MMSEQS_PLUS_SCRIPT: %s\n", MMSEQS_PLUS_SCRIPT );
   printf("# FASTA_TO_HMM_SCRIPT: %s\n", FASTA_TO_HMM_SCRIPT );
   int exit_code = execvp( command[0], command );

   return;
}