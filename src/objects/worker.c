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

/* objects */
#include "objects/f_index.h"
#include "objects/results.h"
#include "objects/alignment.h"
#include "objects/sequence.h"
#include "objects/hmm_profile.h"
#include "objects/edgebound.h"
#include "objects/clock.h"
#include "objects/matrix/matrix_2d.h"
#include "objects/matrix/matrix_3d.h"
#include "objects/mystring.h"
#include "objects/vectors/vector_range_2d.h"

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

   worker->args      = NULL;
   worker->tasks     = NULL;

   worker->q_index   = NULL;
   worker->t_index   = NULL;

   worker->q_file    = NULL;
   worker->t_file    = NULL;

   worker->q_seq     = NULL;
   worker->t_seq     = NULL;
   worker->t_prof    = NULL;

   worker->edg_fwd   = NULL;
   worker->edg_bck   = NULL;
   worker->edg_diag  = NULL;
   worker->edg_row   = NULL;

   worker->traceback = NULL;

   worker->st_MX     = NULL;
   worker->sp_MX     = NULL;
   worker->st_MX3    = NULL;

   worker->times     = NULL;
   worker->scores    = NULL;
   worker->results   = NULL;
   worker->clock     = NULL;

   worker->tasks     = (TASKS*) malloc( sizeof(TASKS) );
   worker->times     = (TIMES*) malloc( sizeof(TIMES) );
   worker->scores    = (SCORES*) malloc( sizeof(SCORES) );
   worker->results   = (RESULTS*) malloc( sizeof(RESULTS) );

   if ( worker->tasks == NULL || worker->times == NULL || worker->scores == NULL || worker->results == NULL ) {
      fprintf(stderr, "ERROR: Failed to malloc WORKER.\n");
      exit(EXIT_FAILURE);
   }

   worker->clock     = CLOCK_Create();
}

/* destructor */
void WORKER_Destroy( WORKER* worker )
{
   free( worker->args    );
   free( worker->tasks   );

   free( worker->q_index );
   free( worker->t_index );

   free( worker->q_file  );
   free( worker->t_file  );

   free( worker->q_seq   );
   free( worker->t_seq   );
   free( worker->t_prof  );

   free( worker->edg_fwd  );
   free( worker->edg_bck  );
   free( worker->edg_diag );
   free( worker->edg_row  );
   free( worker->traceback  );

   free( worker->st_MX  );
   free( worker->sp_MX  );
   free( worker->st_MX3 );

   free( worker->times   );
   free( worker->scores  );
   free( worker->results );
   free( worker->clock   );

   worker->args      = NULL;
   worker->tasks     = NULL;

   worker->q_index   = NULL;
   worker->t_index   = NULL;

   worker->q_file    = NULL;
   worker->t_file    = NULL;

   worker->q_seq     = NULL;
   worker->t_seq     = NULL;
   worker->t_prof    = NULL;

   worker->edg_fwd   = NULL;
   worker->edg_bck   = NULL;
   worker->edg_diag  = NULL;
   worker->edg_row   = NULL;

   worker->traceback = NULL;

   worker->st_MX     = NULL;
   worker->sp_MX     = NULL;
   worker->st_MX3    = NULL;

   worker->times     = NULL;
   worker->scores    = NULL;
   worker->results   = NULL;
   worker->clock     = NULL;

   free( worker );

   worker = NULL;
}