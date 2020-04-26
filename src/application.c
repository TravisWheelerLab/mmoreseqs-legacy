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
   DBG_PRINTF("In DEBUG mode...\n\n");

   /* initialize debugging toolkit */
   #if DEBUG
   {
      debugger             = (DEBUG_KIT*) ERRORCHECK_malloc( sizeof(DEBUG_KIT), __FILE__, __LINE__, __FUNCTION__ );
      debugger->dbfp_path  = DEBUGOUT;
      debugger->dbfp       = fopen( debugger->dbfp_path, "w+" );
      debugger->test_MX    = MATRIX_3D_Create( NUM_NORMAL_STATES, 1, 1 );
      debugger->cloud_MX   = MATRIX_2D_Create( 1, 1 ); 
      /* debugging options */
      debugger->is_embed   = false;
      debugger->is_viz     = false;
   }
   #endif

   /* parse command line arguments */
   args = ARGS_Create();
   ARGS_Parse( args, argc, argv );

   /* output arguments */
   ARGS_Dump( args, stdout );

   /* build worker */
   WORKER* worker = WORKER_Create();
   worker->args = args;
   worker->tasks = (TASKS*) malloc( sizeof(TASKS) );

   /*

   /* jumps to pipeline based on -p flag */
   printf_vall("> Running %s...\n\n", PIPELINE_NAMES[args->pipeline_mode] );
   PIPELINES[ args->pipeline_mode ]( worker );

   /* free debugging toolkit */
   #if DEBUG 
      MATRIX_3D_Destroy( debugger->test_MX );
      MATRIX_2D_Destroy( debugger->cloud_MX );
      fclose( debugout );
      free( debugger );
   #endif

   /* clean up allocated data */
   WORKER_Destroy( worker );

   exit(EXIT_SUCCESS);
}
