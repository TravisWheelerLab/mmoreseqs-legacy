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
      ERRORCHECK_exit(EXIT_FAILURE);
   }

   args->cmdline                 = NULL;
   args->opts                    = NULL;

   args->tmp_folderpath          = NULL;
   args->dbg_folderpath          = NULL;

   args->t_filepath              = NULL;
   args->q_filepath              = NULL; 
   args->t_mmseqs_p_filepath       = NULL;

   args->t_indexpath             = NULL;  
   args->q_indexpath             = NULL;

   args->mmseqs_m8_filepath      = NULL;
   args->hitlist_filepath        = NULL;

   args->output_filepath         = NULL; 
   args->tblout_filepath         = NULL;
   args->m8out_filepath          = NULL;
   args->myout_filepath          = NULL;
   args->mydomout_filepath       = NULL;
   args->mytimeout_filepath      = NULL;
   args->mythreshout_filepath    = NULL;
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
   STR_Destroy( args->cmdline );
   STR_Destroy( args->opts );

   STR_Destroy( args->tmp_folderpath );
   STR_Destroy( args->dbg_folderpath );

   STR_Destroy( args->t_filepath );
   STR_Destroy( args->q_filepath );
   STR_Destroy( args->t_mmseqs_p_filepath );

   STR_Destroy( args->t_indexpath );
   STR_Destroy( args->q_indexpath );

   STR_Destroy( args->mmseqs_m8_filepath );
   STR_Destroy( args->hitlist_filepath );

   STR_Destroy( args->output_filepath );
   STR_Destroy( args->tblout_filepath );
   STR_Destroy( args->m8out_filepath );
   STR_Destroy( args->myout_filepath );
   STR_Destroy( args->mydomout_filepath );
   STR_Destroy( args->mytimeout_filepath );
   STR_Destroy( args->mythreshout_filepath );
   STR_Destroy( args->customout_filepath );

   ERROR_free( args );
   args = NULL;
   return NULL;
}