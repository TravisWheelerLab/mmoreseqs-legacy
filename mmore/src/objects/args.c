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
   args = ERROR_malloc( sizeof(ARGS) );

   args->cmdline                 = NULL;
   args->opts                    = NULL;

   // args->pipeline                = NULL;
   args->pipeline_mode           = -1;
   args->pipeline_name           = NULL;

   args->search_name             = NULL;

   args->tmp_folderpath          = NULL;
   args->dbg_folderpath          = NULL;
   args->prep_folderpath         = NULL;

   args->t_filepath              = NULL;
   args->q_filepath              = NULL;

   args->t_mmore_filepath        = NULL;
   args->t_mmore_p_filepath      = NULL;
   args->t_mmore_s_filepath      = NULL;
   args->q_mmore_filepath        = NULL;

   args->t_mmseqs_filepath       = NULL;
   args->t_mmseqs_p_filepath     = NULL;
   args->t_mmseqs_s_filepath     = NULL;
   args->q_mmseqs_filepath       = NULL;

   args->t_indexpath             = NULL;  
   args->q_indexpath             = NULL;

   args->mmseqs_m8_filepath      = NULL;
   args->hitlist_filepath        = NULL;

   args->allout_fileout          = NULL;
   args->stdout_fileout          = NULL; 
   args->tblout_fileout          = NULL;
   args->m8out_fileout           = NULL;
   args->myout_fileout           = NULL;
   args->mydom_fileout           = NULL;
   args->mytime_fileout          = NULL;
   args->mythresh_fileout        = NULL;
   args->customout_filepath      = NULL;

   args->alpha                   = 12.0f;

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

   STR_Destroy( args->pipeline_name );

   STR_Destroy( args->search_name );

   STR_Destroy( args->tmp_folderpath );
   STR_Destroy( args->dbg_folderpath );
   STR_Destroy( args->prep_folderpath );

   STR_Destroy( args->t_filepath );
   STR_Destroy( args->q_filepath );

   STR_Destroy( args->t_mmore_filepath );
   STR_Destroy( args->t_mmore_p_filepath );
   STR_Destroy( args->t_mmore_s_filepath );
   STR_Destroy( args->q_mmore_filepath );

   STR_Destroy( args->t_mmseqs_filepath );
   STR_Destroy( args->t_mmseqs_p_filepath );
   STR_Destroy( args->t_mmseqs_s_filepath );
   STR_Destroy( args->q_mmseqs_filepath );

   STR_Destroy( args->t_indexpath );
   STR_Destroy( args->q_indexpath );

   STR_Destroy( args->mmseqs_m8_filepath );
   STR_Destroy( args->hitlist_filepath );

   STR_Destroy( args->stdout_fileout );
   STR_Destroy( args->tblout_fileout );
   STR_Destroy( args->m8out_fileout );
   STR_Destroy( args->myout_fileout );
   STR_Destroy( args->mydom_fileout );
   STR_Destroy( args->mytime_fileout );
   STR_Destroy( args->mythresh_fileout );
   STR_Destroy( args->customout_filepath );

   ERROR_free( args );
   args = NULL;
   return NULL;
}