/*******************************************************************************
 *  FILE:      debugger.c
 *  PURPOSE:   DEBUGGER Object.
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
#include "debugger.h"

/* constructor */
DEBUG_KIT* DEBUGGER_Create()
{
   DEBUG_KIT* debugger = (DEBUG_KIT*) malloc( sizeof(DEBUG_KIT) );

   return debugger;
}

/* destructor */
void DEBUGGER_Destroy( DEBUG_KIT* debugger )
{
   free(debugger);
}