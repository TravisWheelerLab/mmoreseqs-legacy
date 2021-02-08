/*******************************************************************************
 *  FILE:      pipeline_interactive.c
 *  PURPOSE:   Pipeline for building hmm's from single fasta sequence.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:      
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

/*! FUNCTION:  	interactive_pipeline()
 *  SYNOPSIS:  	Interactive, tutorialized pipeline.
 *                An interactive pipeline which shows users how to work with tool and purpose of options.
 */
STATUS_FLAG 
interactive_pipeline( WORKER* worker );
STATUS_FLAG 
interactive_pipeline( WORKER* worker )
{
	printf("=== INTERACTIVE PIPELINE ===\n");
}