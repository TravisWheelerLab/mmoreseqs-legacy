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

/*! FUNCTION:  	WORK_index()
 *  SYNOPSIS:  	Load or build target and query index files <t_index> for <q_index>.
 *                Stored in <worker>.
 */
void 
WORK_index( WORKER* worker )
{
   /* load by id */
   WORK_load_indexes_by_id( worker );
}

/*! FUNCTION:  	WORK_load_indexes_by_id()
 *  SYNOPSIS:  	Load or build target and query index files <t_index> for <q_index>.
 */
void 
WORK_load_indexes_by_id( WORKER* worker )
{
   /* load target index */
   CLOCK_Start( worker->clok );
   WORK_load_target_index( worker );
   F_INDEX_Sort_by_Id( worker->t_index );
   CLOCK_Stop( worker->clok ); 
   worker->times->load_target_index = CLOCK_Duration( worker->clok );

   /* load query index */
   CLOCK_Start( worker->clok );
   WORK_load_query_index( worker );
   F_INDEX_Sort_by_Id( worker->q_index );
   CLOCK_Stop( worker->clok ); 
   worker->times->load_query_index = CLOCK_Duration( worker->clok );
}

/*! FUNCTION:  	WORK_load_indexes_by_name()
 *  SYNOPSIS:  	Load or build target index <t_index> and query index <q_index> for <worker>.
 *                Then sorts index by name.
 */
void 
WORK_load_indexes_by_name( WORKER* worker )
{
   /* load target index */
   CLOCK_Start( worker->clok );
   WORK_load_target_index( worker );
   F_INDEX_Sort_by_Name( worker->t_index );
   CLOCK_Stop( worker->clok ); 
   worker->times->load_target_index = CLOCK_Duration( worker->clok );

   /* load query index */
   CLOCK_Start( worker->clok );
   WORK_load_query_index( worker );
   F_INDEX_Sort_by_Name( worker->q_index );
   CLOCK_Stop( worker->clok ); 
   worker->times->load_query_index = CLOCK_Duration( worker->clok );
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

   /* begin timer */
   CLOCK_Start( worker->clok );

   /* default index location  (same as main file but with .idx extension) */
   char* t_indexpath_tmp = NULL;
   if (args->t_indexpath == NULL) {
      char* ext = ".idx";
      t_indexpath_tmp = ERROR_malloc( sizeof(char) * (strlen(args->t_filepath) + strlen(ext) + 1) );
      if (t_indexpath_tmp == NULL) {
           fprintf(stderr, "ERROR: malloc failed.\n");
       }
      strcpy( t_indexpath_tmp, args->t_filepath );
      strcat( t_indexpath_tmp, ext );
   }

   /* load or build target file index */
   if (args->t_indexpath != NULL) {
      /* load file passed by commandline */
      printf_vhi("# loading indexpath from commandline: '%s'...\n", args->t_indexpath );
      worker->t_index = F_INDEX_Load( worker->t_index, args->t_indexpath );
   }
   else if ( access( t_indexpath_tmp, F_OK ) == 0 ) {
      /* check if standard extension index file exists */
      printf_vhi("# found index at database location: '%s'...\n", t_indexpath_tmp );
      args->t_indexpath = strdup( t_indexpath_tmp );
      worker->t_index = F_INDEX_Load( worker->t_index, t_indexpath_tmp );
   }
   else {
      /* build index on the fly */
      printf("# building index of file...\n");
      if (args->t_filetype == FILE_HMM) {
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
      worker->t_index->source_path = strdup(args->t_filepath);
      args->t_indexpath = strdup(t_indexpath_tmp);

      /* save index file */
      fp = fopen( t_indexpath_tmp, "w+" );
      F_INDEX_Dump( worker->t_index, fp );
      fclose(fp);
   }
   ERROR_free(t_indexpath_tmp);

   /* end and save time */
   CLOCK_Stop( worker->clok );
   worker->times->load_target_index = CLOCK_Duration( worker->clok );
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

   /* begin time */
   CLOCK_Start( worker->clok );

   /* default index location (same as main file but with .idx extenstion) */
   char* q_indexpath_tmp = NULL;
   if (args->q_indexpath == NULL) {
      char* ext = ".idx";
      q_indexpath_tmp = ERROR_malloc( sizeof(char) * (strlen(args->q_filepath) + strlen(ext) + 1) );
      if (q_indexpath_tmp == NULL) {
           fprintf(stderr, "ERROR: malloc failed.\n");
       }
      strcpy( q_indexpath_tmp, args->q_filepath );
      strcat( q_indexpath_tmp, ext );
   }

   /* load or build target file index */
   if (args->q_indexpath != NULL) {
      /* load file passed by commandline */
      printf_vhi("# loading indexpath from commandline: '%s'...\n", args->q_indexpath );
      worker->q_index = F_INDEX_Load( worker->q_index, args->q_indexpath );
   } 
   else if ( access( q_indexpath_tmp, F_OK ) == 0 ) {
      /* check if standard extension index file exists */
      printf_vhi("# found index at database location: '%s'...\n", q_indexpath_tmp );
      args->q_indexpath = strdup( q_indexpath_tmp );
      worker->q_index = F_INDEX_Load( worker->q_index, q_indexpath_tmp );
   }  
   else {
      /* build index on the fly */
      printf_vhi("# building index of file...\n");
      if (args->q_filetype == FILE_HMM) {
         worker->q_index = F_INDEX_Hmm_Build( worker->q_index, args->q_filepath );
      }
      else if (args->q_filetype == FILE_FASTA) {
         worker->q_index = F_INDEX_Fasta_Build( worker->q_index, args->q_filepath );
      }
      else {
         fprintf(stderr, "ERROR: query filetype is not supported.\n" );
         exit(EXIT_FAILURE);
      }
      /* identify the query file being indexed */
      worker->q_index->source_path = strdup(args->q_filepath);
      args->q_indexpath = strdup(q_indexpath_tmp);

      /* save index file */
      fp = fopen( q_indexpath_tmp, "w+" );
      F_INDEX_Dump( worker->q_index, fp );
      fclose(fp);
   }
   ERROR_free(q_indexpath_tmp);

   /* end and save time */
   CLOCK_Stop( worker->clok );
   worker->times->load_query_index = CLOCK_Duration( worker->clok );
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

   /* begin time */
   CLOCK_Start( worker->clok );

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
   worker->t_index->source_path = strdup(args->t_filepath);

   /* if output location not specified, then use default naming scheme */
   /* default index location  (same as main file but with .idx extension) */
   if ( args->t_indexpath == NULL ) {
      char* t_indexpath_tmp = NULL;
      char* ext = ".idx";
      t_indexpath_tmp = ERROR_malloc( sizeof(char) * (strlen(args->t_filepath) + strlen(ext) + 1) );
      if (t_indexpath_tmp == NULL) {
           fprintf(stderr, "ERROR: malloc failed.\n");
      }
      strcpy( t_indexpath_tmp, args->t_filepath );
      strcat( t_indexpath_tmp, ext );
      args->t_indexpath = strdup(t_indexpath_tmp);
      ERROR_free(t_indexpath_tmp);
   }

   /* end and save time */
   CLOCK_Stop( worker->clok );
   worker->times->load_target_index = CLOCK_Duration( worker->clok );
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

   /* begin time */
   CLOCK_Start( worker->clok );

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
   worker->q_index->source_path = strdup(args->q_filepath);

   /* if output location not specified, then use default naming scheme */
   /* default index location  (same as main file but with .idx extension) */
   if ( args->q_indexpath == NULL ) {
      char* q_indexpath_tmp = NULL;
      char* ext = ".idx";
      q_indexpath_tmp = ERROR_malloc( sizeof(char) * (strlen(args->q_filepath) + strlen(ext) + 1) );
      if (q_indexpath_tmp == NULL) {
           fprintf(stderr, "ERROR: malloc failed.\n");
      }
      strcpy( q_indexpath_tmp, args->q_filepath );
      strcat( q_indexpath_tmp, ext );
      args->q_indexpath = strdup(q_indexpath_tmp);
      ERROR_free(q_indexpath_tmp);
   }

   /* end and save time */
   CLOCK_Stop( worker->clok );
   worker->times->load_query_index = CLOCK_Duration( worker->clok );
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
      args->t_indexpath = ERROR_malloc( sizeof(char) * (strlen(args->t_filepath) + strlen(ext) + 1) );
      if ( args->t_indexpath == NULL ) {
         fprintf(stderr, "ERROR: malloc failed.\n");
      }
      strcpy( args->t_indexpath, args->t_filepath );
      strcat( args->t_indexpath, ext );
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
      args->q_indexpath = ERROR_malloc( sizeof(char) * (strlen(args->q_filepath) + strlen(ext) + 1) );
      if ( args->q_indexpath == NULL ) {
         fprintf(stderr, "ERROR: malloc failed.\n");
      }
      strcpy( args->q_indexpath, args->q_filepath );
      strcat( args->q_indexpath, ext );
   }
   /* output query index */
   fp = fopen( args->q_indexpath, "w" );
   F_INDEX_Dump( worker->q_index, fp );
   fclose(fp);
}

