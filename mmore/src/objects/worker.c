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
#include "../utilities/_utilities.h"
#include "_objects.h"

/* header */
#include "worker.h"

/*! FUNCTION:  WORKER_Create()
 *  SYNOPSIS:  Create new WORKER object and returns pointer.
 *             Most data is left NULL to be supplied by WORKER_init() and WORK_init().
 */
WORKER* 
WORKER_Create()
{
   /* allocate main object data */
   WORKER* worker;
   worker                  = ERROR_malloc( sizeof(WORKER) );

   /* --- pipeline --- */
   // worker->pipeline        = NULL;
   /* pipeline settings */
   worker->args            = NULL;
   worker->tasks           = NULL;
   /* pipeline tools */
   worker->timer           = CLOCK_Create();
   /* script running tool */
   worker->runner          = SCRIPTRUNNER_Create();

   /* --- reader/writer & buffer --- */
   worker->reader          = NULL;
   // worker->writer          = NULL;
   worker->buffer          = NULL;

   /* --- input files --- */
   worker->q_seq_file      = NULL;
   worker->t_prof_file     = NULL;
   worker->q_index_file    = NULL;
   worker->t_index_file    = NULL;
   worker->mmseqs_file     = NULL;
   worker->hitlist_file    = NULL;

   /* --- output files --- */
   worker->output_file     = NULL;
   worker->tblout_file     = NULL;
   worker->m8out_file      = NULL;
   worker->mydomout_file   = NULL;
   worker->mytimeout_file  = NULL;
   worker->mythreshout_file = NULL;

   /* --- input data --- */
   /* m8 results from mmseqs */
   worker->mmseqs_data     = NULL;
   worker->hitlist_data    = NULL;
   /* indexes of query and target data files */
   worker->q_index         = NULL;
   worker->t_index         = NULL;

   /* --- output data --- */
   /* stats */
   worker->stats           = NULL;
   /* results */
   worker->results         = NULL;
   worker->result          = NULL;
   /* times */
   worker->times           = NULL;
   worker->times_totals    = NULL;
   /* scores */
   // worker->scores          = NULL;
   // worker->evals           = NULL;

   /* --- working data --- */
   /* query and target structs */
   worker->q_seq           = NULL;
   worker->t_seq           = NULL;
   worker->t_prof          = NULL;
   worker->hmm_bg          = NULL;
   /* edgebounds for cloud search */
   worker->edg_fwd         = NULL;
   worker->edg_bck         = NULL;
   worker->edg_diag        = NULL;
   worker->edg_row         = NULL;
   worker->edg_rows_tmp    = NULL;
   /* left and right bound vectors for cloud search */
   for ( int i = 0; i < 3; i++ ) 
   {
      worker->lb_vec[i]    = NULL;
      worker->rb_vec[i]    = NULL;
   }
   /* cloud pruning parameters */
   worker->cloud_params    = (CLOUD_PARAMS) { -1, -1 }; 
   /* tracebacks */
   worker->trace_vit       = NULL;
   worker->trace_post      = NULL;
   /* quadratic space dp matrices */
   worker->st_MX           = NULL;
   worker->st_MX_fwd       = NULL;
   worker->st_MX_bck       = NULL;
   worker->st_MX_post      = NULL;
   worker->st_MX_optacc    = NULL;
   /* linear space dp matrices */
   worker->st_MX3          = NULL;
   worker->st_MX3_fwd      = NULL;
   worker->st_MX3_bck      = NULL;
   /* sparse dp matrices */
   worker->st_SMX          = NULL;
   worker->st_SMX_fwd      = NULL;
   worker->st_SMX_bck      = NULL;
   worker->st_SMX_post     = NULL;
   worker->st_SMX_optacc   = NULL;
   /* special state matrices */
   worker->sp_MX           = NULL;
   worker->sp_MX_fwd       = NULL;
   worker->sp_MX_bck       = NULL;
   worker->sp_MX_post      = NULL;
   worker->sp_MX_optacc    = NULL;
   /* domain definitions */
   worker->dom_def         = NULL;

   /* --- loop data & variables --- */
   /* search id */
   worker->search_rng      = (RANGE) { -1, -1 };  
   worker->search_id       = -1;  
   worker->n_searches      = -1;  
   worker->n_searches_run  = -1;  
   /* mmseqs variables */
   worker->mmseqs_id       = -1; 
   worker->mmseqs_cur      = NULL; 
   worker->mmseqs_prv      = NULL; 
   /* hitlist variables */
   worker->hitlist_id      = -1; 
   worker->hitlist_cur     = NULL;  
   worker->hitlist_prv     = NULL;  
   /* current iteration query/target id */
   worker->q_id            = -1;
   worker->t_id            = -1; 
   worker->q_id_prv        = -1; 
   worker->t_id_prv        = -1;  
   /* current iteration query/target name */
   worker->q_name          = NULL;
   worker->t_name          = NULL; 
   worker->q_name_prv      = NULL; 
   worker->t_name_prv      = NULL;  

   /* --- multi-threading --- */
   worker->N_threads       = 0;
   worker->Nalloc_threads  = 0;
   worker->threads         = NULL;

   return worker;
}

/* TODO: Set so it won't overwrite non-NULL data */
/*! FUNCTION:  WORKER_Init()
 *  SYNOPSIS:  Initialize data pipeline-agnostic data structures.
 */
WORKER* 
WORKER_Init( WORKER* worker )
{
   /* malloc all basic data structures */
   if ( worker->args == NULL) {
      worker->args         = ERROR_malloc( sizeof(ARGS) );
   }
   worker->stats           = ERROR_malloc( sizeof(STATS) );
   worker->tasks           = ERROR_malloc( sizeof(TASKS) );
   worker->times           = ERROR_malloc( sizeof(TIMES) );
   worker->times_totals    = ERROR_malloc( sizeof(TIMES) );
   // worker->scores          = ERROR_malloc( sizeof(ALL_SCORES) );
   /* initialize all values to zero */
   memset( worker->stats,        0, sizeof(STATS) );
   memset( worker->tasks,        0, sizeof(TASKS) ); 
   memset( worker->times,        0, sizeof(TIMES) );
   memset( worker->times_totals, 0, sizeof(TIMES) );
   // memset( worker->scores,       0, sizeof(ALL_SCORES ) );
}

/*! FUNCTION:  WORKER_Create_with_Args()
 *  SYNOPSIS:  Create new WORKER object and adds <args> to worker.
 */
WORKER* 
WORKER_Create_with_Args( ARGS* args )
{
   WORKER* worker;
   worker         = WORKER_Create();
   worker->args   = args;

   return worker;
}

/*! FUNCTION:  WORKER_Create_Threads()
 *  SYNOPSIS:  Creates {N_threads} WORKER_THREAD objects for {worker}.
 *             Stored in {worker->theads}.
 *             Exact number of threads are allocated (should not change during program lifetime).
 */
void
WORKER_Create_Threads(  WORKER*  worker,
                        int      N_threads )
{
   WORKER_THREAD* threads;
   threads = worker->threads;
   threads = ERROR_realloc( threads, sizeof(WORKER_THREAD) * N_threads );
}

/*! FUNCTION:  WORKER_Destroy()
 *  SYNOPSIS:  Frees WORKER object and returns NULL pointer.
 */
WORKER* 
WORKER_Destroy( WORKER* worker )
{
   if (worker == NULL) {
      return worker;
   }

   worker->args            = ARGS_Destroy( worker->args );
   worker->tasks           = ERROR_free( worker->tasks );
   worker->timer           = CLOCK_Destroy( worker->timer );
   worker->runner          = SCRIPTRUNNER_Destroy( worker->runner );

   worker->stats           = ERROR_free( worker->stats );
   worker->times           = ERROR_free( worker->times );
   worker->times_totals    = ERROR_free( worker->times_totals );
   // worker->scores          = ERROR_free( worker->scores );

   worker                  = ERROR_free( worker );
   
   return worker;
}