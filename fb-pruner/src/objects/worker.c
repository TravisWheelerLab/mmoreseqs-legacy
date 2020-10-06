/*******************************************************************************
 *  FILE:      worker.c
 *  PURPOSE:   WORKER object.
 *             Master worker object. Maintains memory for all data structures used in pipeline.
 *             Contains data shared by all worker threads.
 *             TODO: Worker will have thread workers.
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

/** FUNCTION:  WORKER_Create()
 *  SYNOPSIS:  Create new WORKER object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
WORKER* 
WORKER_Create()
{
   WORKER* worker       = NULL;
   
   worker = (WORKER*) ERROR_malloc( sizeof(WORKER) );

   /* set all pointers null */
   worker->args         = NULL;
   worker->tasks        = NULL;

   /* file pointers are all initially NULL */
   worker->output_fp    = NULL;
   worker->tblout_fp    = NULL;
   worker->m8out_fp     = NULL;
   worker->myout_fp     = NULL;

   worker->t_file       = NULL;
   worker->q_file       = NULL;

   worker->t_seq        = NULL;
   worker->t_prof       = NULL;
   worker->q_seq        = NULL;

   worker->t_index      = NULL;
   worker->q_index      = NULL;

   worker->q_seq        = NULL;

   worker->edg_fwd      = NULL;
   worker->edg_bck      = NULL;
   worker->edg_diag     = NULL;
   worker->edg_row      = NULL;

   worker->dom_def      = NULL;

   worker->edg_rows_tmp = NULL;

   for ( int i = 0; i < 3; i++ ) 
   {
      worker->lb_vec[i] = NULL;
      worker->rb_vec[i] = NULL;
   }

   worker->traceback    = NULL;

   worker->st_MX        = NULL;
   worker->st_MX3       = NULL;

   worker->sp_MX        = NULL;
   worker->sp_MX_fwd    = NULL;
   worker->sp_MX_bck    = NULL;

   worker->times        = NULL;
   worker->times_raw    = NULL;

   worker->results      = NULL;
   worker->results_in   = NULL;
   worker->result       = NULL;

   worker->dom_def      = NULL;

   worker->num_searches = 0;

   /* create complex data structs */
   worker->clok         = CLOCK_Create();

   /* malloc all basic data structures */
   worker->tasks        = (TASKS*) ERROR_malloc( sizeof(TASKS) );
   worker->times        = (TIMES*) ERROR_malloc( sizeof(TIMES) );
   worker->times_totals = (TIMES*) ERROR_malloc( sizeof(TIMES) );
   worker->scores       = (NAT_SCORES*) ERROR_malloc( sizeof(NAT_SCORES) );

   /* initialize all values to zero */
   memset( worker->tasks, 0, sizeof(TASKS) ); 
   memset( worker->times, 0, sizeof(TIMES) );
   memset( worker->times_totals, 0, sizeof(TIMES) );
   memset( worker->scores, 0, sizeof(NAT_SCORES ) );

   /* get start time */
   worker->times->program_start = CLOCK_Get_RealTime();

   return worker;
}

/** FUNCTION:  WORKER_Create_with_Args()
 *  SYNOPSIS:  Create new WORKER object and adds <args> to worker.
 */
WORKER* 
WORKER_Create_with_Args( ARGS* args )
{
   WORKER* worker = NULL;
   worker = WORKER_Create();
   worker->args = args;

   return worker;
}

/** FUNCTION:  WORKER_Create_Threads()
 *  SYNOPSIS:  Creates {N_threads} WORKER_THREAD objects for {worker}.
 *             Stored in {worker->theads}.
 *             Exact number of threads are allocated (should not change during program lifetime).
 */
void
WORKER_Create_Threads(  WORKER*  worker,
                        int      N_threads )
{
   WORKER_THREAD* threads = worker->threads;
   threads = (WORKER_THREAD*) ERROR_realloc( threads, sizeof(WORKER_THREAD) );
}

/** FUNCTION:  WORKER_Destroy()
 *  SYNOPSIS:  Frees WORKER object and returns pointer.
 */
void* 
WORKER_Destroy( WORKER* worker )
{
   if (worker == NULL) return worker;

   ARGS_Destroy( worker->args );
   worker->args = NULL;

   worker->clok = CLOCK_Destroy( worker->clok );

   ERROR_free( worker->tasks );
   ERROR_free( worker->times );
   ERROR_free( worker->times_totals );
   ERROR_free( worker->scores );

   worker = ERROR_free( worker );
   return worker;
}