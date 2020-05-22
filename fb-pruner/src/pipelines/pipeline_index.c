/*******************************************************************************
 *  FILE:      pipeline_index.c
 *  PURPOSE:   Index pipeline...
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
#include "structs.h"
#include "utilities.h"
#include "objects.h"
#include "parsers.h"
#include "algs_linear.h"
#include "algs_quad.h"
#include "algs_naive.h"

/* header */
#include "pipelines.h"

/* ****************************************************************************************** *
 *  FUNCTION:  index_pipeline()
 *  SYNOPSIS:  Indexing workflow pipeline.
 *             Indexes a target and query file, then saves those indexes to file.
/* ****************************************************************************************** */
void index_pipeline( WORKER* worker )
{
   FILE*       fp       = NULL;

   ARGS*       args     = worker->args;
   TASKS*      tasks    = worker->tasks;
   CLOCK*      clok     = worker->clok;

   // printf("load target index...\n");
   // WORK_load_target_index( worker );
   // printf("load query index...\n");
   // WORK_load_query_index( worker );

   // printf("output target index...\n");
   // WORK_output_target_index( worker );
   // printf("output query index...\n");
   // WORK_output_query_index( worker );

   /* begin time */
   CLOCK_Start(clok);

   /* create and save target index */
   {
      /* if no index location given, use default (same as main file but with .idx extension) */
      if ( args->t_indexpath == NULL )
      {
         char* ext = ".idx";
         args->t_indexpath = (char*) malloc( sizeof(char) * (strlen(args->t_filepath) + strlen(ext) + 1) );
         if ( args->t_indexpath == NULL ) {
            fprintf(stderr, "ERROR: malloc failed.\n");
         }
         strcpy( args->t_indexpath, args->t_filepath );
         strcat( args->t_indexpath, ext );
      }

      /* build index on the fly */
      printf("building index of file...\n");
      if ( args->t_filetype == FILE_HMM ) {
         worker->t_index = F_INDEX_Hmm_Build(args->t_filepath);
      }
      else if (args->t_filetype == FILE_FASTA) {
         worker->t_index = F_INDEX_Fasta_Build(args->t_filepath);
      }
      else {
         fprintf(stderr, "ERROR: target filetype is not supported.\n" );
         exit(EXIT_FAILURE);
      }
      /* identify the target file being indexed */
      worker->t_index->source_path = strdup(args->t_filepath);

      /* save index file */
      F_INDEX_Save( worker->t_index, args->t_indexpath );
   }

   /* create and save target index */
   {
      /* if no index location given, use default (same as main file but with .idx extension) */
      if ( args->q_indexpath == NULL )
      {
         char* ext = ".idx";
         args->q_indexpath = (char*) malloc( sizeof(char) * (strlen(args->q_filepath) + strlen(ext) + 1) );
         if ( args->q_indexpath == NULL ) {
            fprintf(stderr, "ERROR: malloc failed.\n");
         }
         strcpy( args->q_indexpath, args->q_filepath );
         strcat( args->q_indexpath, ext );
      }

      /* build index on the fly */
      printf("building index of file...\n");
      if ( args->q_filetype == FILE_HMM ) {
         worker->q_index = F_INDEX_Hmm_Build(args->q_filepath);
      }
      else if (args->q_filetype == FILE_FASTA) {
         worker->q_index = F_INDEX_Fasta_Build(args->q_filepath);
      }
      else {
         fprintf(stderr, "ERROR: target filetype is not supported.\n" );
         exit(EXIT_FAILURE);
      }
      /* identify the target file being indexed */
      worker->q_index->source_path = strdup(args->q_filepath);

      /* save index file */
      F_INDEX_Save( worker->q_index, args->q_indexpath );
   }

   // printf("\nTARGET:\n");
   // F_INDEX_Dump( worker->t_index, stdout );
   // printf("QUERY:\n");
   // F_INDEX_Dump( worker->q_index, stdout );
}
