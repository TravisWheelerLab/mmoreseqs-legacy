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
#include "objects/args.h"

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

   worker->t_index   = NULL;
   worker->q_index   = NULL;

   worker->t_file    = NULL;
   worker->q_file    = NULL;

   worker->t_seq     = NULL;
   worker->t_prof    = NULL;
   worker->q_seq     = NULL;

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

   /* malloc all basic data structures */
   worker->tasks        = (TASKS*) calloc( 1, sizeof(TASKS) ); /* sets all tasks to false */
   worker->report       = (REPORT*) calloc( 1, sizeof(REPORT) ); /* sets all tasks to false */
   worker->times        = (TIMES*) malloc( sizeof(TIMES) );
   worker->scores       = (SCORES*) malloc( sizeof(SCORES) );
   worker->results      = (RESULTS*) malloc( sizeof(RESULTS) );
   worker->result       = (RESULT*) malloc( sizeof(RESULT) );
   if ( worker->tasks == NULL || worker->times == NULL || worker->scores == NULL || worker->results == NULL ) {
      fprintf(stderr, "ERROR: Failed to malloc WORKER.\n");
      exit(EXIT_FAILURE);
   }

   /* create edgebounds objects */
   worker->edg_fwd   = EDGEBOUNDS_Create();
   worker->edg_bck   = EDGEBOUNDS_Create();
   worker->edg_diag  = EDGEBOUNDS_Create();
   worker->edg_row   = EDGEBOUNDS_Create();

   worker->clock     = CLOCK_Create();

   return worker;
}

/* destructor */
void WORKER_Destroy( WORKER* worker )
{
   if (worker == NULL) return;

   ARGS_Destroy( worker->args );
   free( worker->tasks );

   F_INDEX_Destroy( worker->q_index );
   F_INDEX_Destroy( worker->t_index );

   free( worker->q_file  );
   free( worker->t_file  );

   SEQUENCE_Destroy( worker->q_seq   );
   SEQUENCE_Destroy( worker->t_seq   );
   HMM_PROFILE_Destroy( worker->t_prof  );

   EDGEBOUNDS_Destroy( worker->edg_fwd  );
   EDGEBOUNDS_Destroy( worker->edg_bck  );
   EDGEBOUNDS_Destroy( worker->edg_diag );
   EDGEBOUNDS_Destroy( worker->edg_row  );
   ALIGNMENT_Destroy( worker->traceback );

   MATRIX_3D_Destroy( worker->st_MX3 );
   MATRIX_3D_Destroy( worker->st_MX  );
   MATRIX_2D_Destroy( worker->sp_MX  );

   free( worker->times   );
   free( worker->scores  );
   free( worker->results );
   free( worker->result  );

   CLOCK_Destroy( worker->clock   );

   free( worker );
   worker = NULL;
}