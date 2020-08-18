/*******************************************************************************
 *  FILE:      pipeline_index.c
 *  PURPOSE:   Index pipeline...
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
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
#include "objects.h"
#include "parsers.h"
#include "algs_linear.h"
#include "algs_quad.h"
#include "algs_naive.h"

/* header */
#include "pipelines.h"

/* set default args for index pipeline */


/* ****************************************************************************************** *
 *  FUNCTION:  index_pipeline()
 *  SYNOPSIS:  Indexing workflow pipeline.
 *             Indexes a target and query file, then saves those indexes to file.
/* ****************************************************************************************** */
void index_pipeline( WORKER* worker )
{
   FILE*       fp       = NULL;

   ARGS*       args     = worker->args;
   TASKS*      tasks    = worker->tasks;
   CLOCK*      clok     = worker->clok;

   /* begin time */
   CLOCK_Start(clok);

   /* initialize data structures needed for tasks */
   WORK_init( worker );

   /* load or build file indexes */
   WORK_load_target_index( worker );
   WORK_load_query_index( worker );
   /* sort file indexes */
   F_INDEX_Sort_by_Id( worker->t_index );
   F_INDEX_Sort_by_Id( worker->q_index );

   /* print index out to their outfiles */
   WORK_output_target_index( worker );
   WORK_output_query_index( worker );
   /* clean up worker data structs */
   WORK_cleanup( worker );
}

/* set default arguments for index_pipeline */
void index_set_args( WORKER* worker )
{

}
