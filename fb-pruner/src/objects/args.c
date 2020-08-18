/*******************************************************************************
 *  FILE:      alignment.c
 *  PURPOSE:   ALIGNMENT Object.
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
#include "alignment.h"

/* constructor */
ARGS* ARGS_Create()
{
   ARGS* args = NULL;

   args = (ARGS*) calloc( 1, sizeof(ARGS) );
   if (args == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc ARGS.\n");
      exit(EXIT_FAILURE);
   }

   args->t_filepath = NULL;
   args->q_filepath = NULL; 

   args->t_indexpath = NULL;  
   args->q_indexpath = NULL;

   args->output_filepath = NULL; 

   args->mmseqs_res_filepath = NULL;
   args->mmseqs_tmp_filepath = NULL;

   args->t_lookup_filepath = NULL;
   args->q_lookup_filepath = NULL;

   return args;
}

/* destructor */
void ARGS_Destroy(ARGS* args)
{
   if (args == NULL) return;
   
   /* free all strings */
   ERRORCHECK_free( args->t_filepath );
   ERRORCHECK_free( args->q_filepath );

   ERRORCHECK_free( args->t_indexpath );
   ERRORCHECK_free( args->q_indexpath );

   ERRORCHECK_free( args->output_filepath );

   ERRORCHECK_free( args->mmseqs_res_filepath );
   ERRORCHECK_free( args->mmseqs_tmp_filepath );

   ERRORCHECK_free( args->t_lookup_filepath );
   ERRORCHECK_free( args->q_lookup_filepath );

   ERRORCHECK_free( args );
   args = NULL;
}