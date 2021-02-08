/*******************************************************************************
 *  FILE:      work_index.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Build, load and retrieve file indexes.
 *
 *  AUTHOR:    Dave Rich
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
#include "../algs_sparse/_algs_sparse.h"
#include "../reporting/_reporting.h"

/* header */
#include "_work.h"
#include "work_index.h"

/*! FUNCTION:  	WORK_load_indexes()
 *  SYNOPSIS:  	Load or build target and query index files <t_index> for <q_index>.
 *                Stored in <worker>.
 */
void 
WORK_load_indexes( WORKER* worker )
{
   STATS* stats = worker->stats;

   /* load by id */
   WORK_load_indexes_by_name( worker );

   /* pull database size from index */
   stats->n_query_db = worker->q_index->N;
   printf("INDEX SIZES = %d, %d\n", worker->q_index->N, worker->t_index->N );
}

/*! FUNCTION:  	WORK_load_indexes_by_id()
 *  SYNOPSIS:  	Load or build target and query index files <t_index> for <q_index>.
 */
void 
WORK_load_indexes_by_id( WORKER* worker )
{
   /* load target index */
   CLOCK_Start( worker->timer );
   WORK_load_target_index( worker );
   F_INDEX_Sort_by_Id( worker->t_index );
   CLOCK_Stop( worker->timer ); 
   worker->times->load_target_index = CLOCK_Duration( worker->timer );

   /* load query index */
   CLOCK_Start( worker->timer );
   WORK_load_query_index( worker );
   F_INDEX_Sort_by_Id( worker->q_index );
   CLOCK_Stop( worker->timer ); 
   worker->times->load_query_index = CLOCK_Duration( worker->timer );
}

/*! FUNCTION:  	WORK_load_indexes_by_name()
 *  SYNOPSIS:  	Load or build target index <t_index> and query index <q_index> for <worker>.
 *                Then sorts index by name.
 */
void 
WORK_load_indexes_by_name( WORKER* worker )
{
   ARGS* args  = worker->args;

   /* load target index */
   CLOCK_Start( worker->timer );
   WORK_load_target_index( worker );
   F_INDEX_Sort_by_Name( worker->t_index );
   CLOCK_Stop( worker->timer ); 
   worker->times->load_target_index = CLOCK_Duration( worker->timer );
   /* report index */
   // if ( args->verbose_level >= VERBOSE_HIGH ) 
   // {
   //    printf("=== TARGET F_INDEX ===\n");
   //    F_INDEX_Dump( worker->t_index, stdout );
   // }

   /* load query index */
   CLOCK_Start( worker->timer );
   WORK_load_query_index( worker );
   F_INDEX_Sort_by_Name( worker->q_index );
   CLOCK_Stop( worker->timer ); 
   worker->times->load_query_index = CLOCK_Duration( worker->timer );
   /* report index */
   // if ( args->verbose_level >= VERBOSE_HIGH ) 
   // {
   //    printf("=== QUERY F_INDEX ===\n");
   //    F_INDEX_Dump( worker->q_index, stdout );
   // }
}

/*! FUNCTION:  	WORK_load_target_index()
 *  SYNOPSIS:  	Load target index <t_index> for <worker>.
 */
void 
WORK_load_target_index( WORKER* worker )
{
   FILE*    fp       = NULL;
   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;

   CLOCK_Start( worker->timer );

   /* default index location  (same as main file but with .idx extension) */
   char* t_indexpath_tmp = NULL;
   if (args->t_indexpath == NULL) {
      /* index file extension */
      char* ext = ".idx";
      t_indexpath_tmp = STR_Concat( args->t_filepath, ext );
   }

   /* if target index path was given at the command line, load that */
   if (args->t_indexpath != NULL) {
      printf_vhi("# loading indexpath from commandline: '%s'...\n", args->t_indexpath );
      worker->t_index = F_INDEX_Load( worker->t_index, args->t_indexpath );
   }
   /* else, check if index exists at default file location */
   else if ( access( t_indexpath_tmp, F_OK ) == 0 ) {
      printf_vhi("# found index at database location: '%s'...\n", t_indexpath_tmp );
      args->t_indexpath = STR_Create( t_indexpath_tmp );
      worker->t_index   = F_INDEX_Load( worker->t_index, t_indexpath_tmp );
   }
   /* else, build index on the fly */
   else {
      printf_vhi("# building index of file...\n");
      if (args->t_filetype == FILE_HMM) {
         worker->t_index = F_INDEX_Hmm_Build( worker->t_index, args->t_filepath );
      }
      elif (args->t_filetype == FILE_FASTA) {
         worker->t_index = F_INDEX_Fasta_Build( worker->t_index, args->t_filepath );
      }
      else {
         fprintf(stderr, "ERROR: target filetype is not supported.\n" );
         exit(EXIT_FAILURE);
      }
      /* identify the query file being indexed */
      worker->t_index->source_path  = STR_Create( args->t_filepath );
      args->t_indexpath             = STR_Create( t_indexpath_tmp );

      /* save index file */
      fp = fopen( t_indexpath_tmp, "w+" );
      F_INDEX_Dump( worker->t_index, fp );
      fclose(fp);
   }
   STR_Destroy(t_indexpath_tmp);

   CLOCK_Stop( worker->timer );
   worker->times->load_target_index = CLOCK_Duration( worker->timer );
}

/*! FUNCTION:  	WORK_load_query_index()
 *  SYNOPSIS:  	Load query index <q_index> for <worker>.
 */
void 
WORK_load_query_index( WORKER* worker )
{
   FILE*    fp       = NULL;
   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;

   CLOCK_Start( worker->timer );

   /* default index location (same as main file but with .idx extenstion) */
   char* q_indexpath_tmp = NULL;
   if ( args->q_indexpath == NULL ) {
      /* index file extension */
      char* ext = ".idx";
      q_indexpath_tmp = STR_Concat( args->q_filepath, ext );
   }

   /* if query index path was given at the command line, load that */
   if ( args->q_indexpath != NULL ) {
      /* load file passed by commandline */
      printf_vhi("# loading indexpath from commandline: '%s'...\n", args->q_indexpath );
      worker->q_index = F_INDEX_Load( worker->q_index, args->q_indexpath );
   } 
   /* else, check if index exists at default file location */
   elif ( access( q_indexpath_tmp, F_OK ) == 0 ) {
      printf_vhi("# found index at database location: '%s'...\n", q_indexpath_tmp );
      args->q_indexpath    = STR_Create( q_indexpath_tmp );
      worker->q_index      = F_INDEX_Load( worker->q_index, q_indexpath_tmp );
   }  
   /* else, build index on the fly */
   else {
      printf_vhi("# building index of file...\n");
      if (args->q_filetype == FILE_HMM) {
         worker->q_index = F_INDEX_Hmm_Build( worker->q_index, args->q_filepath );
      }
      elif ( args->q_filetype == FILE_FASTA ) {
         worker->q_index = F_INDEX_Fasta_Build( worker->q_index, args->q_filepath );
      }
      else {
         fprintf(stderr, "ERROR: query filetype is not supported.\n" );
         exit(EXIT_FAILURE);
      }
      /* identify the query file being indexed */
      worker->q_index->source_path  = STR_Create(args->q_filepath);
      args->q_indexpath             = STR_Create(q_indexpath_tmp);

      /* save index file */
      fp = fopen( q_indexpath_tmp, "w+" );
      F_INDEX_Dump( worker->q_index, fp );
      fclose(fp);
   }
   STR_Destroy(q_indexpath_tmp);

   
   CLOCK_Stop( worker->timer );
   worker->times->load_query_index = CLOCK_Duration( worker->timer );
}

/*! FUNCTION:  	WORK_build_target_index()
 *  SYNOPSIS:  	Build target index <t_index> for <worker>.
 */
void 
WORK_build_target_index( WORKER* worker )
{
   FILE*    fp       = NULL;
   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;

   CLOCK_Start( worker->timer );

   /* build index on the fly */
   if ( args->t_filetype == FILE_HMM ) {
      worker->t_index = F_INDEX_Hmm_Build( worker->t_index, args->t_filepath );
   }
   else if (args->t_filetype == FILE_FASTA) {
      worker->t_index = F_INDEX_Fasta_Build( worker->t_index, args->t_filepath );
   }
   else {
      fprintf(stderr, "ERROR: target filetype is not supported.\n" );
      exit(EXIT_FAILURE);
   }
   /* identify the query file being indexed */
   worker->t_index->source_path = STR_Create(args->t_filepath);

   /* if output location not specified, then use default naming scheme */
   /* default index location  (same as main file but with .idx extension) */
   if ( args->t_indexpath == NULL ) {
      char* ext = ".idx";
      args->t_indexpath = STR_Concat( args->t_filepath, ext );
   }

   CLOCK_Stop( worker->timer );
   worker->times->load_target_index = CLOCK_Duration( worker->timer );
}

/*! FUNCTION:  	WORK_build_query_index()
 *  SYNOPSIS:  	Build query index <q_index> for <worker>.
 */
void 
WORK_build_query_index( WORKER* worker )
{
   FILE*    fp       = NULL;
   ARGS*    args     = worker->args;
   TASKS*   tasks    = worker->tasks;

   CLOCK_Start( worker->timer );

   /* build index on the fly */
   if ( args->q_filetype == FILE_HMM ) {
      worker->q_index = F_INDEX_Hmm_Build( worker->q_index, args->q_filepath );
   }
   else if ( args->q_filetype == FILE_FASTA ) {
      worker->q_index = F_INDEX_Fasta_Build( worker->q_index, args->q_filepath );
   }
   else {
      fprintf(stderr, "ERROR: query filetype is not supported.\n" );
      exit(EXIT_FAILURE);
   }
   /* identify the query file being indexed */
   worker->q_index->source_path = STR_Create(args->q_filepath);

   /* if output location not specified, then use default naming scheme */
   /* default index location  (same as main file but with .idx extension) */
   if ( args->q_indexpath == NULL ) {
      char* ext = ".idx";
      args->q_indexpath = STR_Concat( args->q_filepath, ext );
   }

   CLOCK_Stop( worker->timer );
   worker->times->load_query_index = CLOCK_Duration( worker->timer );
}

/*! FUNCTION:  	WORK_output_target_index()
 *  SYNOPSIS:  	Write target index out to <t_indexpath>.
 */
void 
WORK_output_target_index( WORKER* worker )
{
   FILE*    fp    = NULL;
   ARGS*    args  = worker->args;

   /* determine the output file to save to */
   if ( args->t_indexpath == NULL ) {
      /* if no output name is given, append ".idx" and save in same directory */
      const char* ext = ".idx";
      args->t_indexpath = STR_Concat( args->t_filepath, ext );
   }
   /* output target index */
   fp = fopen( args->t_indexpath, "w" );
   F_INDEX_Dump( worker->t_index, fp );
   fclose(fp);
}

/*! FUNCTION:  	WORK_output_target_index()
 *  SYNOPSIS:  	Write query index out to <q_indexpath>.
 */
void 
WORK_output_query_index( WORKER* worker )
{
   FILE*    fp    = NULL;
   ARGS*    args  = worker->args;

   /* determine the output file to save to */
   if ( args->q_indexpath == NULL ) {
      /* if no output name is given, append ".idx" and save in same directory */
      const char* ext = ".idx";
      args->q_indexpath = STR_Concat( args->q_filepath, ext );
   }
   /* output query index */
   fp = fopen( args->q_indexpath, "w" );
   F_INDEX_Dump( worker->q_index, fp );
   fclose(fp);
}

