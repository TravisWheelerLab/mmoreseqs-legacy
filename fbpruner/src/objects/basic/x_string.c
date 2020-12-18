/*******************************************************************************
 *  FILE:      x_string.c
 *  PURPOSE:   X_STRING object
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
#include "../structs.h"
#include "../../utilities/utilities.h"
#include "../objects.h"

/* header */
#include "x_string.h"

/*
 *  FUNCTION:  X_STRING_Create()
 *  SYNOPSIS:
 */
 X_STRING* 
 X_STRING_Create( char* chars )
{
   if ( chars == NULL ) return NULL;

   X_STRING* str    	= ERROR_malloc( sizeof(X_STRING) );
   str->data      	= strdup( chars );
   str->N         	= strlen( chars );

   return str;
}

/*
 *  FUNCTION:  X_STRING_Destroy()
 *  SYNOPSIS:
 */
X_STRING* 
X_STRING_Destroy( X_STRING* str )
{
   ERROR_free( str->data );
   ERROR_free( str );
   return NULL;
}