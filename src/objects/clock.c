/*******************************************************************************
 *  @file clock.c
 *  @brief CLOCK object
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
#include "structs.h"

/* header */
#include "clock.h"

CLOCK* clock_Create()
{
   const int min_size = 16;

   CLOCK*cl = (CLOCK*) malloc( sizeof(CLOCK) );
   cl->start = 0;
   cl->stop = 0;

   cl->N = 0;
   cl->N = min_size;
   cl->stamps = (float*) malloc( sizeof(float) * min_size );
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

time_t clock_Ticks(CLOCK*cl)
{
   cl->duration = cl->stop - cl->start;
   return cl->duration;
}

time_t clock_pTicks(CLOCK*cl, char*str)
{
   printf("%s took %d ticks\n", str, clock_Ticks(cl) );
   return cl->duration;
}

float clock_Secs(CLOCK*cl)
{
   cl->duration = cl->stop - cl->start;
   return ticks_to_msec(cl->duration);
}

time_t clock_pSecs(CLOCK*cl, char*str)
{
   printf("%s took %d msecs\n", str, clock_Secs(cl) );
   return cl->duration;
}

float ticks_to_msec(time_t t)
{
   return ((float)t * 1000.f) / ((float)CLOCKS_PER_SEC);
}