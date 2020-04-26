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
   
   free( args->t_filepath );
   free( args->q_filepath );

   free( args->t_indexpath );
   free( args->q_indexpath );

   free( args->output_filepath );

   free( args->mmseqs_res_filepath );
   free( args->mmseqs_tmp_filepath );

   free( args->t_lookup_filepath );
   free( args->q_lookup_filepath );

   free( args );
   args = NULL;
}