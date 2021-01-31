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
int 
main ( int argc, char *argv[] )
{
   /* times */
   float program_start, program_end, program_runtime;

   /* initialize debugging toolkit */
   #if DEBUG
   {
      printf("===> IN DEBUG MODE <===\n");

      debugger = DEBUGGER_Create("test-output");
      /* debugging options */
      debugger->is_embed   = true;
      debugger->is_viz     = true;
   }
   #endif

   /* initialize worker and args object */
   ARGS*    args;
   WORKER*  worker;

   args     = ARGS_Create();
   worker   = WORKER_Create_with_Args( args );

   /* ideally, the clock would start before all, but this will work */
   program_start = CLOCK_Get_Time( worker->clok );

   /* initialize random number generator */
   RNG_Init();

   /* parse command line arguments */
   ARGS_Parse( args, argc, argv );

   /* output arguments */
   ARGS_Dump( args, stdout );
   
   /* Run pipeline determined by args */
   PIPELINES[ args->pipeline_mode ].func( worker );

   /* free debugging toolkit */
   #if DEBUG 
   {
      DEBUGGER_Destroy( debugger );
   }
   #endif

   program_end = CLOCK_Get_Time( worker->clok );
   program_runtime = CLOCK_Get_Diff( worker->clok, program_start, program_end );

   /* clean up allocated data */
   WORKER_Destroy( worker );

   printf("# Program_Runtime: %f\n", program_runtime );
   printf("# Completed Successfully.\n");
   exit(EXIT_SUCCESS);
}
