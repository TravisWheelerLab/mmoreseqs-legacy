/*******************************************************************************
 *  FILE:      workerer.c
 *  PURPOSE:   WORKER object.c
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
#include "worker.h"

/* constructor */
WORKER* WORKER_Create()
{
   WORKER* worker = NULL;
   
   worker = (WORKER*) malloc( sizeof(WORKER) );
   if (worker == NULL) {
      fprintf(stderr, "ERROR: Failed to malloc WORKER.\n");
      exit(EXIT_FAILURE);
   }

   /* set all pointers null */
   worker->args      = NULL;
   worker->tasks     = NULL;
   worker->report    = NULL;

   worker->t_file    = NULL;
   worker->q_file    = NULL;

   worker->t_seq     = NULL;
   worker->t_prof    = NULL;
   worker->q_seq     = NULL;

   worker->t_index   = NULL;
   worker->q_index   = NULL;

   worker->edg_fwd      = NULL;
   worker->edg_bck      = NULL;
   worker->edg_diag     = NULL;
   worker->edg_row      = NULL;
   worker->edg_rows_tmp = NULL;

   worker->traceback = NULL;

   worker->st_MX     = NULL;
   worker->sp_MX     = NULL;
   worker->st_MX3    = NULL;

   worker->times        = NULL;
   worker->scores       = NULL;
   worker->results      = NULL;
   worker->results_in   = NULL;
   worker->result       = NULL;
   worker->clok         = NULL;

   /* create complex data structs */
   worker->clok      = CLOCK_Create();

   /* malloc all basic data structures */
   worker->tasks        = (TASKS*) calloc( 1, sizeof(TASKS) );    /* sets all tasks to false */
   worker->report       = (REPORT*) calloc( 1, sizeof(REPORT) );  /* sets all report fields to false */
   worker->times        = (TIMES*) malloc( sizeof(TIMES) );       
   worker->scores       = (SCORES*) malloc( sizeof(SCORES) );     

   if ( worker->tasks == NULL || worker->report == NULL || worker->times == NULL || worker->scores == NULL ) 
   {
      fprintf(stderr, "ERROR: Failed to malloc WORKER.\n");
      exit(EXIT_FAILURE);
   }

   return worker;
}

/* constructor with args supplied */
WORKER* WORKER_Create_with_Args( ARGS* args )
{
   WORKER* worker = NULL;
   worker = WORKER_Create();
   worker->args = args;

   return worker;
}

/* destructor */
void* WORKER_Destroy( WORKER* worker )
{
   if (worker == NULL) return worker;

   ARGS_Destroy( worker->args );
   worker->args = NULL;

   worker->clok      = CLOCK_Destroy( worker->clok );

   free( worker->tasks   );
   free( worker->report  );
   free( worker->times   );
   free( worker->scores  );

   free( worker );
   worker = NULL;
   return worker;
}