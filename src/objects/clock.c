/*******************************************************************************
 *  FILE:      clock.c
 *  PURPOSE:   CLOCK Object for testing runtimes. 
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Unknown.
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

/* constructor */
CLOCK* CLOCK_Create()
{
   CLOCK*      cl       = NULL;   
   const int   min_size = 16; 

   cl = (CLOCK*) malloc( sizeof(CLOCK) );
   if (cl == NULL) {
      const char* obj_name = "CLOCK";
      fprintf(stderr, "ERROR: Unable to malloc object %s.\n", obj_name);
      exit(EXIT_FAILURE);
   }

   cl->start   = 0;
   cl->stop    = 0;
   cl->N       = 0;
   cl->Nalloc  = min_size;

   cl->stamps  = (float*) malloc( sizeof(float) * min_size );
   if (cl->stamps == NULL) {
      const char* obj_name = "CLOCK STAMPS";
      fprintf(stderr, "ERROR: Unable to malloc object %s: <%p>.\n", obj_name, cl);
      exit(EXIT_FAILURE);
   }

   return cl;
}

/* destructor */
void CLOCK_Destroy(CLOCK*cl)
{
   free(cl->stamps);
   free(cl);
}

/* */
time_t CLOCK_Start(CLOCK*cl)
{
   cl->start = clock();
   return cl->start;
}

/* */
time_t CLOCK_Stop(CLOCK*cl)
{
   cl->stop = clock();
   return cl->stop;
}

/* */
time_t CLOCK_Ticks(CLOCK*cl)
{
   cl->duration = cl->stop - cl->start;
   return cl->duration;
}

/* */
time_t CLOCK_pTicks(CLOCK*cl, char*str)
{
   printf("%s took %d ticks\n", str, CLOCK_Ticks(cl) );
   return cl->duration;
}

/* */
float CLOCK_Secs(CLOCK*cl)
{
   cl->duration = cl->stop - cl->start;
   return ticks_to_msec(cl->duration);
}

/* */
time_t CLOCK_Print_Secs(CLOCK*cl, char*str)
{
   printf("%s took %d msecs\n", str, CLOCK_Secs(cl) );
   return cl->duration;
}

/* */
float ticks_to_msec(time_t t)
{
   return ((float)t * 1000.f) / ((float)CLOCKS_PER_SEC);
}
