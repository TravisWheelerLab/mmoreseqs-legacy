/*******************************************************************************
 *  FILE:      pipeline_hmmbuild.c
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

/* header */
#include "_pipelines.h"

/* mmseqs pipeline */
void hmmbuild_pipeline( WORKER* worker )
{
	printf("=== HMMBUILD PIPELINE ===\n");

	/* initialize worker */
	WORK_init( worker );

   /* worker objects */
   ARGS*       args        = worker->args;
   TASKS*      tasks       = worker->tasks;
   TIMES*      times       = worker->times;

   /* load or build index */
   WORK_index( worker );

   // if ( args->t_filetype == FILE_FASTA ) {
   // 	WORK_convert_target_to_hmm( WORKER* worker );
   // }
   // if ( args->q_filetype == FILE_FASTA ) {
   // 	WORK_convert_query_to_hmm( WORKER* worker );
   // }

   WORK_load_target_index( worker );

   SEQUENCE_to_HMM_PROFILE( worker->t_seq, worker->t_prof );

   WORK_cleanup( worker );
}