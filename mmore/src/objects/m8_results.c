/*******************************************************************************
 *  FILE:      m8_results.c
 *  PURPOSE:   M8_RESULTS object.
 *             Stores results in .m8 format.
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
#include "../utilities/_utilities.h"

/* header */
#include "_objects.h"
#include "m8_results.h"

/*! FUNCTION:  M8_RESULTS_Create()
 *  SYNOPSIS:  Create <results>, allocate memory and return pointer.
 */
M8_RESULTS* 
M8_RESULTS_Create()
{
   M8_RESULTS*    results  = NULL;;
   const int      min_size = 16;

   results = ERROR_malloc( sizeof(M8_RESULTS) );

   results->N         = 0;
   results->Nalloc    = 0;
   results->data      = NULL;

   M8_RESULTS_Resize( results, min_size );

   return results;
}

/*! FUNCTION:  M8_RESULTS_Destroy()
 *  SYNOPSIS:  Destroy <results>, free memory and return NULL pointer.
 */
M8_RESULTS* 
M8_RESULTS_Destroy( M8_RESULTS* results )
{
   if ( results == NULL ) return results;

   for (int i = 0; i < results->N; i++) {
      ERROR_free(results->data[i].target_name);
      ERROR_free(results->data[i].query_name);
   }

   results->data  = ERROR_free(results->data);
   results        = ERROR_free(results);
   return results;
}

/*! FUNCTION:  M8_RESULTS_Pushback()
 *  SYNOPSIS:  Add <res> to <results> list, resize array if full.
 */
void 
M8_RESULTS_Pushback(    M8_RESULTS*    results, 
                        M8_RESULT*     res )
{
   results->data[results->N] = *res;

   results->N++;
   if ( results->N >= results->Nalloc ) {
      M8_RESULTS_Resize( results, results->Nalloc * 2 );
   }
}

/*! FUNCTION:  M8_RESULTS_Pushback()
 *  SYNOPSIS:  Resize <results> list to be able to store <size> entries.
 */
void 
M8_RESULTS_Resize(   M8_RESULTS* 	results,
                     size_t   	   size )
{
   results->Nalloc = size;
   results->data   = ERROR_realloc( results->data, sizeof(RESULT) * size );
}

/*! FUNCTION:  M8_RESULTS_GetX()
 *  SYNOPSIS:  Get <i>th entry in <results> list.
 */
M8_RESULT* 
M8_RESULTS_GetX(  M8_RESULTS* 	results,
                  int   	      i )
{
   return &(results->data[i]);
}

/*! FUNCTION:  M8_RESULTS_Swap_Target_and_Query()
 *  SYNOPSIS:  The target and query are cross-labeled between MMSEQS and MMORE.
 *             Ideally, this should be remedied and MMORE should be swapped (todo list).
 *             For now, this will swap the fields so that m8 results are corrected for MMORE.
 */
M8_RESULT* 
M8_RESULTS_Swap_Target_and_Query(   M8_RESULTS* 	results )
{
   M8_RESULT*  result;
   int         tmp_int;
   STR         tmp_str;

   for (int i = 0; i < results->N; i++) {
      result = M8_RESULTS_GetX( results, i );
      /* swap all fields */
      INT_Swap( &result->query_id, &result->target_id );
      STR_Swap( &result->query_name, &result->target_name );
      INT_Swap( &result->q_beg, &result->t_beg );
      INT_Swap( &result->q_end, &result->t_end );
   }
}

/*! FUNCTION:  M8_RESULTS_Dump()
 *  SYNOPSIS:  Output all entries in <results> in .m8 format to file <fp>.
 */
void 
M8_RESULTS_Dump(  M8_RESULTS*    results,
                  FILE*          fp )
{
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to open file.\n");
      ERRORCHECK_exit(EXIT_FAILURE);
   }

   int max_len = 0;
   for (int i = 0; i < results->N; i++) 
   {
      M8_RESULT* res = M8_RESULTS_GetX( results, i );

      if ( max_len < strlen(res->query_name) ) {
         max_len = strlen(res->query_name);
      }

      if ( max_len < strlen(res->target_name) ) {
         max_len = strlen(res->target_name);
      }
   }

   /* print headers */
   fprintf( fp, "# M8+ M8_RESULTS:\n");

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

   /* print each result */
   for (int i = 0; i < results->N; i++) 
   {
      M8_RESULT* res = M8_RESULTS_GetX( results, i );

      fprintf( fp, "%d\t",    res->result_id );
      fprintf( fp, "%d\t",    res->query_id );
      fprintf( fp, "%d\t",    res->target_id );

      fprintf( fp, "%s\t",    res->query_name );
      fprintf( fp, "%s\t",    res->target_name );

      fprintf( fp, "%.3f\t",  res->perc_id );
      fprintf( fp, "%d\t",    res->aln_len );
      fprintf( fp, "%d\t",    res->mismatch );
      fprintf( fp, "%d\t",    res->gap_openings );

      fprintf( fp, "%d\t",    res->q_beg );
      fprintf( fp, "%d\t",    res->q_end );
      fprintf( fp, "%d\t",    res->t_beg );
      fprintf( fp, "%d\t",    res->t_end );

      fprintf( fp, "%.3e\t",  res->eval );
      fprintf( fp, "%.1f\t",  res->bitsc );
   }
}

/*! FUNCTION:  M8_RESULTS_Dump()
 *  SYNOPSIS:  Output single <result> to file <fp>.
 */
void 
M8_RESULT_Dump(   M8_RESULT*  res,
                  FILE*       fp )
{  
   fprintf( fp, "{ ");

   fprintf( fp, "%s: %s | ",    "QUERY",    res->query_name );
   fprintf( fp, "%s: %s | ",    "TARGET",   res->target_name );
   fprintf( fp, "%s: %.3f | ",  "PERCID",   res->perc_id );
   fprintf( fp, "%s: %d | ",    "ALN_LEN",  res->aln_len );
   fprintf( fp, "%s: %d | ",    "MISMAT",   res->mismatch );
   fprintf( fp, "%s: %d | ",    "GAPS",     res->gap_openings );

   fprintf( fp, "%s: %d-%d | ",  "Q_RANGE",  res->q_beg, res->q_end );
   fprintf( fp, "%s: %d-%d | ",  "T_RANGE",  res->t_beg, res->t_end );

   fprintf( fp, "%s: %.3e | ",  "EVAL",     res->eval );
   fprintf( fp, "%s: %.1f   ",  "BITSC",    res->bitsc );

   fprintf( fp, "}\n");
}

/*! FUNCTION:  M8_RESULTS_Dump()
 *  SYNOPSIS:  Output single <result> to file <fp>.
 */
void 
M8_RESULT_Dump_old(   M8_RESULT*  res,
                      FILE*       fp )
{   
   fprintf( fp, "{ ");

   fprintf( fp, "%s ",    res->query_name );
   fprintf( fp, "%s ",    res->target_name );

   fprintf( fp, "%.3f ",  res->perc_id );
   fprintf( fp, "%d ",    res->aln_len );
   fprintf( fp, "%d ",    res->mismatch );
   fprintf( fp, "%d ",    res->gap_openings );

   fprintf( fp, "%d ",    res->q_beg );
   fprintf( fp, "%d ",    res->q_end );
   fprintf( fp, "%d ",    res->t_beg );
   fprintf( fp, "%d ",    res->t_end );

   fprintf( fp, "%.3e ",  res->eval );
   fprintf( fp, "%.1f ",  res->bitsc );

   fprintf( fp, "}\n");
}