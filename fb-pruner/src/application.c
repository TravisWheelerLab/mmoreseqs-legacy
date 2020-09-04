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
#include "structs.h"
#include "utilities.h"
#include "parsers.h"
#include "pipelines.h"
#include "objects.h"

/* === HEADER === */

/* === MAIN ENTRY-POINT TO PROGRAM === */
int main ( int argc, char *argv[] )
{
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

   /* initialize random number generator */
   RNG_Init();

   /* parse command line arguments */
   ARGS_Parse( args, argc, argv );

   /* report initial header */
   REPORT_stdout_header( worker, stdout );

   /* output arguments */
   ARGS_Dump( args, stdout );
   
   /* jumps to pipeline based on -p flag */
   PIPELINES[ args->pipeline_mode ]( worker );

   /* free debugging toolkit */
   #if DEBUG 
   {
      DEBUGGER_Destroy( debugger );
   }
   #endif

   /* clean up allocated data */
   WORKER_Destroy( worker );

   exit(EXIT_SUCCESS);
}
