/*******************************************************************************
 *  FILE:      cmd_opts.c
 *  PURPOSE:   CMD_OPTS Object. Used for Parsing Commandline Args.
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
#include "_utilities.h"
#include "_objects.h"

/* header */
#include "opts.h"

/*
 *  FUNCTION:  OPTS_Create()
 *  SYNOPSIS:  
 */
OPTS* OPTS_Create()
{
   OPTS* opts = (OPTS*) ERROR_malloc( sizeof(OPTS) );

   return opts;
}

/*
 *  FUNCTION:  OPTS_Destroy()
 *  SYNOPSIS:  Returns Null pointer.
 */
OPTS* OPTS_Destroy( OPTS* opts )
{
   if (opts == NULL) return NULL;
   
   ERROR_free( opts );

   return NULL;
}