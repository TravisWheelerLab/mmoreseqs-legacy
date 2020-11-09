/*******************************************************************************
 *  FILE:      workerer.c
 *  PURPOSE:   WORKER_THREAD object.
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
#include "../utilities/utilities.h"
#include "objects.h"

/* header */
#include "worker.h"

/* constructor */
WORKER_THREAD* 
WORKER_THREAD_Create()
{
   WORKER_THREAD* worker = NULL;
   
   worker = (WORKER_THREAD*) ERROR_malloc( sizeof(WORKER_THREAD) );

   return worker;
}

/* destructor */
void* 
WORKER_THREAD_Destroy( WORKER_THREAD* worker )
{
   ERROR_free( worker );
   return NULL;
}