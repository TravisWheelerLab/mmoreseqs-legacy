/*******************************************************************************
 *  @file traceback.c
 *  @brief Traceback Datatype
 *
 *  @author Dave Rich
 *  @bug Lots.
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
#include "traceback.h"


TRACEBACK* traceback_Create()
{
   const int min_size = 16;
   TRACEBACK*tr = (TRACEBACK*) malloc( sizeof(TRACEBACK) );
   tr->traces = (TRACE*) malloc( sizeof(TRACE) * min_size );
}

