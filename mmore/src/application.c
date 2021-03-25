/*******************************************************************************
 *  FILE:      application.c
 *  PURPOSE:   Entry Point to Application: Worker/Debugger Initialization, Argument Parsing, Chooses Pipeline
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

/* import stdlib */
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
// #include <dos.h> /* date function */

/* import local libraries */
#include "easel.h"

/* include local files */
#include "objects/structs.h"
#include "utilities/_utilities.h"
#include "parsers/_parsers.h"
#include "pipelines/_pipelines.h"
#include "objects/_objects.h"

/* === HEADER === */

/* === MAIN ENTRY-POINT TO PROGRAM === */
STATUS_FLAG 
main ( int argc, char *argv[] )
{
   /* full program runtimes */
   float program_start, program_end, program_runtime;

   /* initialize debugging toolkit */
   #if DEBUG
   {
      printf("===> IN DEBUG MODE <===\n");

      debugger = DEBUGGER_Create(DEBUG_FOLDER "");
      /* debugging options */
      debugger->is_embed   = true;
      debugger->is_viz     = true;
   }
   #endif

   /* initialize worker and args object */
   WORKER* worker = WORKER_Create();
   WORKER_Init( worker );

   /* ideally, the clock would start before all, but this will work */
   program_start = CLOCK_GetTime( worker->timer );

   /* initialize random number generator */
   RNG_Init();

   /* create commandline object */
   COMMANDLINE_Load( worker->cmd, argc, argv );
   COMMANDLINE_SimpleDump( worker->cmd, stdout );

   /* parse command line arguments */
   ARGS_Parse( worker->args, argc, argv, worker->cmd, worker->arg_opts );

   /* output arguments */
   if ( worker->args->verbose_level >= VERBOSE_LOW ) {
      ARGS_Dump( worker->args, stdout );
   }

   /* Run pipeline determined by args */
   PIPELINES[ worker->args->pipeline_mode ].pipeline_main( worker );

   /* free debugging toolkit */
   #if DEBUG 
   {
      DEBUGGER_Destroy( debugger );
   }
   #endif

   program_end       = CLOCK_GetTime( worker->timer );
   program_runtime   = CLOCK_GetDiff( worker->timer, program_start, program_end );

   /* clean up allocated data */
   WORKER_Destroy( worker );

   printf("# Program_Runtime: %f\n", program_runtime );
   printf("# Completed Successfully.\n");
   ERRORCHECK_exit(EXIT_SUCCESS);
}
