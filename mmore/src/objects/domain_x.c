/*******************************************************************************
 *  FILE:      dom_def.c
 *  PURPOSE:   DOMAIN_X object.
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
#include "domain_x.h"

/*! FUNCTION:  DOMAIN_X_Create()
 *  SYNOPSIS:  Create new DOMAIN_X object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
DOMAIN_X* 
DOMAIN_X_Create()
{
   DOMAIN_X* dom = NULL;
   dom = (DOMAIN_X*) ERROR_malloc( sizeof(DOMAIN_X) );

   return dom;
}

/*! FUNCTION:  DOMAIN_X_Destroy()
 *  SYNOPSIS:  Frees DOMAIN_X object and returns NULL pointer.
 */
DOMAIN_X* 
DOMAIN_X_Destroy( DOMAIN_X* dom )
{
   if (dom == NULL) return dom;

   dom = ERROR_free( dom );

   return dom;
}

/*! FUNCTION:  DOMAIN_X_Reuse()
 *  SYNOPSIS:  
 */
int 
DOMAIN_X_Reuse(   DOMAIN_X*    dom )
{
   
   return STATUS_SUCCESS;
}