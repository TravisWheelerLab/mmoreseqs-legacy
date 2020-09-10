/*******************************************************************************
 *  FILE:      results.c
 *  PURPOSE:   RESULTS object
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

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "results.h"

/*
 *  FUNCTION:  RESULTS_Create()
 *  SYNOPSIS:
 */
RESULTS* RESULTS_Create()
{
   RESULTS*    results  = NULL;;
   const int   min_size = 16;

   results = (RESULTS*) ERROR_malloc( sizeof(RESULTS) );

   results->filepath  = NULL;
   results->N         = 0;
   results->Nalloc    = 0;
   results->data      = NULL;

   RESULTS_Resize( results, min_size );

   return results;
}

/*
 *  FUNCTION:  RESULTS_Destroy()
 *  SYNOPSIS:
 */
void* RESULTS_Destroy( RESULTS* results )
{
   if ( results == NULL ) return results;

   for (int i = 0; i < results->N; i++) {
      ERROR_free(results->data[i].target_name);
      ERROR_free(results->data[i].query_name);
   }

   ERROR_free(results->data);
   results->data = NULL;
   ERROR_free(results->filepath);
   results->filepath = NULL;
   
   ERROR_free(results);
   results = NULL;
   return results;
}

/* add result to results */
void RESULTS_Pushback( RESULTS*  results, 
                       RESULT*   res )
{
   results->data[results->N] = *res;

   results->N++;
   if ( results->N >= results->Nalloc ) {
      RESULTS_Resize( results, results->Nalloc * 2 );
   }
}

/* resize results */
void RESULTS_Resize( RESULTS* results,
                     int      size )
{
   results->Nalloc = size;
   results->data   = (RESULT*) ERROR_realloc( results->data, sizeof(RESULT) * size );
}

/* output results to file pointer */
void RESULTS_M8_Dump( RESULTS*   results,
                           FILE*      fp )
{
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to open file.\n");
      exit(EXIT_FAILURE);
   }

   int max_len = 0;
   for (int i = 0; i < results->N; i++) {
      RESULT* res = &(results->data[i]);

      if ( max_len < strlen(res->query_name) ) {
         max_len = strlen(res->query_name);
      }

      if ( max_len < strlen(res->target_name) ) {
         max_len = strlen(res->target_name);
      }
   }

   /* print headers */
   fprintf( fp, "# M8+ RESULTS:\n");

   fprintf( fp, "# " );
   fprintf( fp, "%s\n",    "res_id" );
   fprintf( fp, "%s\t",    "q_id" );
   fprintf( fp, "%s\t",    "t_id" );

   fprintf( fp, "%s\t",    "query" );
   fprintf( fp, "%s\t",    "target" );

   fprintf( fp, "%s\t",    "perc_id" );
   fprintf( fp, "%s\t",    "aln_len" );
   fprintf( fp, "%s\t",    "mismatch" );
   fprintf( fp, "%s\t",    "gaps" );

   fprintf( fp, "%s\t",    "q_st" );
   fprintf( fp, "%s\t",    "q_end");
   fprintf( fp, "%s\t",    "t_st" );
   fprintf( fp, "%s\t",    "t_end" );

   fprintf( fp, "%s\t",    "eval" );
   fprintf( fp, "%s\t",    "bit-sc" );
   fprintf( fp, "%s\n",    "cloud-fwd-sc");

   /* print each result */
   for (int i = 0; i < results->N; i++) {
      RESULT* res = &(results->data[i]);

      fprintf( fp, "%d\t",    res->result_id );
      fprintf( fp, "%d\t",    res->query_id );
      fprintf( fp, "%d\t",    res->target_id );

      fprintf( fp, "%s\t",    res->query_name );
      fprintf( fp, "%s\t",    res->target_name );

      fprintf( fp, "%.3f\t",  res->perc_id );
      fprintf( fp, "%d\t",    res->aln_len );
      fprintf( fp, "%d\t",    res->mismatch );
      fprintf( fp, "%d\t",    res->gap_openings );

      fprintf( fp, "%d\t",    res->query_start );
      fprintf( fp, "%d\t",    res->query_end );
      fprintf( fp, "%d\t",    res->target_start );
      fprintf( fp, "%d\t",    res->target_end );

      fprintf( fp, "%.3e\t",  res->e_value );
      fprintf( fp, "%d\t",    res->bit_score );
      
      fprintf( fp, "%f\n",    res->cloud_fwd_sc);
   }
}

/* output results to file pointer */
void RESULTS_My_Dump( RESULTS*   results,
                      FILE*      fp )
{
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to open file.\n");
      exit(EXIT_FAILURE);
   }

   int max_len = 0;
   for (int i = 0; i < results->N; i++) {
      RESULT* res = &(results->data[i]);

      if ( max_len < strlen(res->query_name) ) {
         max_len = strlen(res->query_name);
      }

      if ( max_len < strlen(res->target_name) ) {
         max_len = strlen(res->target_name);
      }
   }

   /* headers */
   fprintf( fp, "# M8 RESULTS:\n");

   fprintf( fp, "# " );
   fprintf( fp, "%s\t",    "query" );
   fprintf( fp, "%s\t",    "target" );

   fprintf( fp, "%s\t",    "perc_id" );
   fprintf( fp, "%s\t",    "aln_len" );
   fprintf( fp, "%s\t",    "mismatch" );
   fprintf( fp, "%s\t",    "gaps" );

   fprintf( fp, "%s\t",    "q_st" );
   fprintf( fp, "%s\t",    "q_end");
   fprintf( fp, "%s\t",    "t_st" );
   fprintf( fp, "%s\t",    "t_end" );

   fprintf( fp, "%s\t",    "eval" );
   fprintf( fp, "%s\t",    "bit_sc" );

   fprintf( fp, "%s\n",    "cloud_fwd_sc" );

   for (int i = 0; i < results->N; i++) {
      RESULT* res = &(results->data[i]);

      fprintf( fp, "%s\t",    res->query_name );
      fprintf( fp, "%s\t",    res->target_name );

      fprintf( fp, "%.3f\t",  res->perc_id );
      fprintf( fp, "%d\t",    res->aln_len );
      fprintf( fp, "%d\t",    res->mismatch );
      fprintf( fp, "%d\t",    res->gap_openings );

      fprintf( fp, "%d\t",    res->query_start );
      fprintf( fp, "%d\t",    res->query_end );
      fprintf( fp, "%d\t",    res->target_start );
      fprintf( fp, "%d\t",    res->target_end );

      fprintf( fp, "%.3e\t",  res->e_value );
      fprintf( fp, "%d\t",    res->bit_score );

      fprintf( fp, "%f\n",    res->cloud_fwd_sc );
   }
}

/* TODO: run print commands on single line */
/* output single results to file pointer */
void RESULT_M8_Dump( RESULT*  res,
                     FILE*    fp )
{   
   printf("# M8 SINGLE:\t");

   fprintf( fp, "%s\t",    res->query_name );
   fprintf( fp, "%s\t",    res->target_name );

   fprintf( fp, "%.3f\t",  res->perc_id );
   fprintf( fp, "%d\t",    res->aln_len );
   fprintf( fp, "%d\t",    res->mismatch );
   fprintf( fp, "%d\t",    res->gap_openings );

   fprintf( fp, "%d\t",    res->query_start );
   fprintf( fp, "%d\t",    res->query_end );
   fprintf( fp, "%d\t",    res->target_start );
   fprintf( fp, "%d\t",    res->target_end );

   fprintf( fp, "%.3e\t",  res->e_value );
   fprintf( fp, "%d\n",    res->bit_score );
}
