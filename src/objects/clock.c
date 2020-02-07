/*******************************************************************************
 *  @file clock.c
 *  @brief Clock object
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
#include <ctype.h>
#include <time.h>

/* local imports */
#include "clock.h"

CLOCK* clock_Create()
{
   CLOCK*cl = (CLOCK*) malloc( sizeof(CLOCK) );
   cl->start = 0;
   cl->stop = 0;
   cl->duration = 0;
}

void clock_Destroy(CLOCK*cl)
{
   free(cl);
}

time_t clock_Start(CLOCK*cl)
{
   cl->start = clock();
   return cl->start;
}

time_t clock_Stop(CLOCK*cl)
{
   cl->stop = clock();
   return cl->stop;
}

time_t clock_pTicks(CLOCK*cl, char*str)
{
   cl->stop = clock();
   cl->duration = cl->stop - cl->start;
   printf("%s took %d ticks\n", str, cl->duration);
   return cl->duration;
}

float ticks_to_msec(time_t t)
{
   float new_t;
   new_t = t * 1000.0 / CLOCKS_PER_SEC;
   return new_t;
}