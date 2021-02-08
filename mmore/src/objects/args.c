/*******************************************************************************
 *  FILE:      args.c
 *  PURPOSE:   ARGS Object. Used for Parsing Commandline Args.
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
#include "alignment.h"

/*! FUNCTION:  ARGS_Create()
 *  SYNOPSIS:  
 */
ARGS* 
ARGS_Create()
{
   ARGS* args = NULL;

   args = (ARGS*) calloc( 1, sizeof(ARGS) );
   if (args == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc ARGS.\n");
      exit(EXIT_FAILURE);
   }

   args->cmdline                 = NULL;
   args->opts                    = NULL;

   args->tmp_folderpath          = NULL;
   args->dbg_folderpath          = NULL;

   args->t_filepath              = NULL;
   args->q_filepath              = NULL; 
   args->t_indexpath             = NULL;  
   args->q_indexpath             = NULL;

   args->mmseqs_res_filepath     = NULL;
   args->hitlist_filepath        = NULL;

   args->output_filepath         = NULL; 
   args->tblout_filepath         = NULL;
   args->m8out_filepath          = NULL;
   args->myout_filepath          = NULL;
   args->mydomout_filepath       = NULL;
   args->customout_filepath      = NULL;

   return args;
}

/*! FUNCTION:  ARGS_Destroy()
 *  SYNOPSIS:  
 */
ARGS* 
ARGS_Destroy( ARGS* args )
{
   if (args == NULL) return NULL;
   
   /* free all strings */
   ERROR_free( args->cmdline );
   ERROR_free( args->opts );

   ERROR_free( args->tmp_folderpath );
   ERROR_free( args->dbg_folderpath );

   ERROR_free( args->t_filepath );
   ERROR_free( args->q_filepath );
   ERROR_free( args->t_indexpath );
   ERROR_free( args->q_indexpath );

   ERROR_free( args->mmseqs_res_filepath );
   ERROR_free( args->hitlist_filepath );

   ERROR_free( args->output_filepath );
   ERROR_free( args->tblout_filepath );
   ERROR_free( args->m8out_filepath );
   ERROR_free( args->myout_filepath );
   ERROR_free( args->mydomout_filepath );
   ERROR_free( args->customout_filepath );

   ERROR_free( args );
   args = NULL;
   return NULL;
}