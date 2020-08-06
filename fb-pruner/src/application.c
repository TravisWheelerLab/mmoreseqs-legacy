/*******************************************************************************
 *  FILE:      application.c
 *  PURPOSE:   Entry Point to Application, Argument Parsing, Choose Pipeline
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
#include <sys/types.h>
#include <sys/stat.h>

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "parsers.h"
#include "pipelines.h"
#include "objects.h"

/* === HEADER === */

/* === MAIN ENTRY-POINT TO PROGRAM === */
int main ( int argc, char *argv[] )
{
   /* initialize random number generator */
   RNG_Init();

   /* initialize debugging toolkit */
   #if DEBUG
   {
      printf("=== IN DEBUG MODE ===\n");

      debugger             = (DEBUG_KIT*) ERRORCHECK_malloc( sizeof(DEBUG_KIT), __FILE__, __LINE__, __FUNCTION__ );
      /* debugging options */
      debugger->is_embed   = true;
      debugger->is_viz     = true;

      /* create output folder for testing */
      char* test_folder = "test_output/";
      if ( mkdir(test_folder, 0777) == -1 ) {
         fprintf(stderr, "Could not create folder: %s\n", test_folder);
         // exit(EXIT_FAILURE);
      }   
   }
   #endif

   /* parse command line arguments */
   args = ARGS_Create();
   ARGS_Parse( args, argc, argv );

   /* output arguments */
   ARGS_Dump( args, stdout );

   /* build worker */
   WORKER* worker = WORKER_Create_with_Args( args );

   /* jumps to pipeline based on -p flag */
   printf_vlo("> Running %s...\n\n", PIPELINE_NAMES[args->pipeline_mode] );
   PIPELINES[ args->pipeline_mode ]( worker );

   /* free debugging toolkit */
   #if DEBUG 
   {
      // fclose( debugout );
      free( debugger );
   }
   #endif

   /* clean up allocated data */
   WORKER_Destroy( worker );

   printf("...exited successfully.\n");
   exit(EXIT_SUCCESS);
}
