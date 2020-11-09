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
#include "../objects/structs.h"
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../parsers/parsers.h"
#include "../algs_linear/algs_linear.h"
#include "../algs_quad/algs_quad.h"
#include "../algs_naive/algs_naive.h"

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
   WORK_open( worker );

   /* building, sorting, and outputting target index */
   printf_vlo("# building target index from:\t%s\n", 
      args->t_filepath );
   WORK_build_target_index( worker );
   F_INDEX_Sort_by_Id( worker->t_index );
   WORK_output_target_index( worker );
   printf_vlo("# outputting target index to:\t%s\n",
      args->t_indexpath );

   /* building, sorting, and outputting target index */
   printf_vlo("# building query index from:\t%s\n", 
      args->q_filepath );
   WORK_build_query_index( worker );
   F_INDEX_Sort_by_Id( worker->q_index );
   WORK_output_query_index( worker );
   printf_vlo("# outputting query index to:\t%s\n",
      args->q_indexpath );
   
   /* clean up worker data structs */
   WORK_close( worker );
   WORK_cleanup( worker );
}

/* set default arguments for index_pipeline */
void index_set_args( WORKER* worker )
{

}
