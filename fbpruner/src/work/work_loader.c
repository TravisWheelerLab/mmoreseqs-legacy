/*******************************************************************************
 *  FILE:      work_loader.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Loads query sequences and target hmm profiles.
 *
 *  AUTHOR:    Dave Rich
 *  BUGS:   
 *    - Currently, the load_by_name() functions are the only ones in use.
 *    - The names could be stored in a map for quicker lookups.
 *  NOTES:
 *    - Eventually, plan for these functions is that state is handled entirely by the 
 *      WORKER object, so other arguments are unneccessary (encapsulated by WORKER).   
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
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../parsers/_parsers.h"
#include "../algs_linear/_algs_linear.h"
#include "../algs_quad/_algs_quad.h"
#include "../algs_naive/_algs_naive.h"
#include "../algs_sparse/_algs_sparse.h"
#include "../reporting/_reporting.h"

/* header */
#include "_work.h"
#include "work_loader.h"

/*! FUNCTION:  	WORK_load_target_by_id()
 *  SYNOPSIS:  	Loads mmseqs input m8 file into <results_in>, located at <mmseqs_res_filepath>.
 *                Verifies that search index range does not exceed bounds of list.
 */
void 
WORK_load_mmseqs_list( WORKER* worker )
{
   ARGS* args = worker->args;

   /* If range values are negative, then range is set to full list set */
   if ( args->list_range.beg < 0 && args->list_range.end < 0 ) {
      args->list_range.beg = 0;
      args->list_range.end = INT_MAX;
   }

   /* m8+ file contains target_id, query_id, and result_id fields */
   RESULTS_M8_Parse( 
      worker->results_in, args->mmseqs_res_filepath, args->list_range.beg, args->list_range.end );

   /* Truncate or extract valid result range */
   args->list_range.beg = MAX(args->list_range.beg, 0);
   args->list_range.end = MIN(args->list_range.end, worker->results_in->N + args->list_range.beg);
}

/*! FUNCTION:  	WORK_load_target_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_target_by_id(    WORKER*     worker,
                           int         id )
{
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   CLOCK*         clok           = worker->clok;
   int            index_id       = 0;
   int            index_offset   = 0;

   /* begin time */
   CLOCK_Start(clok);

   /* search and load target by offset */
   index_id = F_INDEX_Search_Id( worker->t_index, id );
   if ( index_id == -1 ) {
      fprintf(stderr, "ERROR: Target id '%d' not found in F_INDEX.\n", id );
      exit(EXIT_FAILURE);
   }
   WORK_load_target_by_index_id( worker, index_id );

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_target = CLOCK_Duration(clok);
}

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_query_by_id(  WORKER*     worker,
                        int         id )
{
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   CLOCK*         clok           = worker->clok;
   int            index_id       = 0;
   int            index_offset   = 0;

   /* begin time */
   CLOCK_Start(clok);

   /* search and load target by offset */
   index_id       = F_INDEX_Search_Id( worker->q_index, id );
   if ( index_id == -1 ) {
      fprintf(stderr, "ERROR: Query id '%d' not found in F_INDEX.\n", id );
      exit(EXIT_FAILURE);
   }
   WORK_load_query_by_index_id( worker, index_id );

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_query = CLOCK_Duration(clok);
}

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <name> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_target_by_name(  WORKER*    worker,
                           char*      name )
{
   // printf_vall("==> WORK_load_target_by_name()\n");
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   CLOCK*         clok           = worker->clok;
   int            index_id       = 0;
   int            index_offset   = 0;

   /* begin time */
   CLOCK_Start(clok);

   /* search and load target by offset */
   index_id = F_INDEX_Search_Name( worker->t_index, name );
   if ( index_id == -1 ) {
      fprintf(stderr, "ERROR: Target name '%s' not found in F_INDEX.\n", name );
      exit(EXIT_FAILURE);
   }
   WORK_load_target_by_index_id( worker, index_id );

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_target = CLOCK_Duration(clok);
}

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <name> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_query_by_name(   WORKER*     worker,
                           char*       name )
{
   // printf_vall("==> WORK_load_query_by_name()\n");
   ARGS*          args           = worker->args;
   TASKS*         tasks          = worker->tasks;
   TIMES*         times          = worker->times;
   CLOCK*         clok           = worker->clok;
   int            index_id       = 0;
   int            index_offset   = 0;

   /* begin time */
   CLOCK_Start(clok);

   /* search and load target by offset */
   index_id = F_INDEX_Search_Name( worker->q_index, name );
   if ( index_id == -1 ) {
      fprintf(stderr, "ERROR: Query name '%s' not found in F_INDEX.\n", name );
      exit(EXIT_FAILURE);
   }
   WORK_load_query_by_index_id( worker, index_id);

   /* end and save time */
   CLOCK_Stop(clok);
   times->load_query = CLOCK_Duration(clok);
}

/*! FUNCTION:  	WORK_load_target_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_target_by_index_id(    WORKER*     worker,
                                 int         index_id )
{
   ARGS*          args     = worker->args;
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      t_seq    = worker->t_seq;
   SEQUENCE*      q_seq    = worker->q_seq;
   F_INDEX_NODE*  my_idx   = &worker->t_index->nodes[index_id];

   worker->t_id = index_id;

   /* load target profile by file type */
   switch ( args->t_filetype )
   {
      case FILE_HMM: {
         HMM_PROFILE_Parse( worker->t_prof, args->t_filepath, my_idx->offset ); 
         HMM_PROFILE_Convert_NegLog_To_Real( worker->t_prof );
         HMM_PROFILE_Config( worker->t_prof, args->search_mode );
         // HMM_PROFILE_Dump( worker->t_prof, stdout );
      } break;
      case FILE_FASTA: {
         SEQUENCE_Fasta_Parse( worker->t_seq, args->t_filepath, my_idx->offset );
         SEQUENCE_to_HMM_PROFILE( worker->t_seq, worker->t_prof );
         HMM_PROFILE_Dump( worker->t_prof, stdout );
         exit(EXIT_SUCCESS);
      } break;
      default: {
         fprintf(stderr, "ERROR: Only HMM and FASTA filetypes are supported for targets.\n");
         exit(EXIT_FAILURE);
      }
   }
}

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void 
WORK_load_query_by_index_id(  WORKER*     worker,
                              int         index_id )
{
   ARGS*          args     = worker->args;
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      t_seq    = worker->t_seq;
   SEQUENCE*      q_seq    = worker->q_seq;
   F_INDEX_NODE*  my_idx   = &worker->q_index->nodes[index_id];

   worker->q_id = index_id;

   /* load query by file type */
   switch ( args->q_filetype )
   {
      /* fasta only supported file type */
      case FILE_FASTA: {
         SEQUENCE_Fasta_Parse( worker->q_seq, args->q_filepath, my_idx->offset );
         // SEQUENCE_Dump( worker->q_seq, stdout );
      } break;
      case FILE_HMM: {

      } 
      default: {
         fprintf(stderr, "ERROR: Only FASTA filetypes are supported for queries.\n");
         exit(EXIT_FAILURE);
      }
   }

   /* set special state transitions based on query sequence length */
   if ( t_prof != NULL ) 
   {
      HMM_PROFILE_ReconfigLength( worker->t_prof, worker->q_seq->N );
   } else {
      fprintf(stderr, "ERROR: Target profile must be loaded before Query Sequence. Currently NULL.\n");
      exit(EXIT_FAILURE);
   }
}

