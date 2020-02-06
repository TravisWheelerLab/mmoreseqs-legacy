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

}

time_t clock_Start(CLOCK* cl)
{

}

time_t clock_Stop(CLOCK* cl)
{

}

float clocktime_to_msec(time_t t)
{
   float new_t;
   new_t = t * 1000.0 / CLOCKS_PER_SEC;
   return new_t;
}