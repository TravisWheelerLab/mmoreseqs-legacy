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
#include "../utilities/_utilities.h"
#include "_objects.h"

/* header */
#include "results.h"

/*! FUNCTION:  RESULTS_Create()
 *  SYNOPSIS:  
 */
RESULTS* 
RESULTS_Create()
{
   RESULTS*    results  = NULL;;
   const int   min_size = 16;

   results = (RESULTS*) ERROR_malloc( sizeof(RESULTS) );

   results->N         = 0;
   results->Nalloc    = 0;
   results->data      = NULL;
   
   results->max_in_queue   = 0;
   results->num_hits       = 0;
   results->num_searches   = 0;

   RESULTS_Resize( results, min_size );

   return results;
}

/*! FUNCTION:  RESULTS_Destroy()
 *  SYNOPSIS:
 */
RESULTS* 
RESULTS_Destroy( RESULTS* results )
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

/* add result to results */
void 
RESULTS_Pushback(    RESULTS*  results, 
                     RESULT*   res )
{
   results->data[results->N] = *res;

   results->N++;
   if ( results->N >= results->Nalloc ) {
      RESULTS_Resize( results, results->Nalloc * 2 );
   }
}

/* resize results */
void 
RESULTS_Resize(   RESULTS*    results,
                  size_t      size )
{
   results->Nalloc = size;
   results->data   = (RESULT*) ERROR_realloc( results->data, sizeof(RESULT) * size );
}

/* output results to file pointer */
void 
RESULTS_Dump(  RESULTS*   results,
               FILE*      fp )
{
   
}

/* TODO: run print commands on single line */
/* output single results to file pointer */
void 
RESULT_Dump(   RESULT*  res,
               FILE*    fp )
{   
   
}
