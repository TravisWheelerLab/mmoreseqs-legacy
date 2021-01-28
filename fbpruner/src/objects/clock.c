/*******************************************************************************
 *  FILE:      clok.c
 *  PURPOSE:   CLOCK Object for testing runtimes. 
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *    - None.
 *  NOTES:
 *    - After implementing a simple clock, it did not run correctly on the cluster.
 *      So, I have wrapped the Easel stopwatch in this implementation.
 *      The implementation version is selected by the CLOCK_TYPE macro at compile time.
 *  ISSUES:
 *    - Easel's stopwatch doesn't allow for timing multiple events at the same time,
 *      To capture times, my code just quickly stops easel stopwatch, then starts it again, 
 *      storing the accumulating offset value. Functional, if a little hacky.
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
#include <sys/times.h>

/* easel imports */
#include "../../lib/easel/esl_stopwatch.h"

/* local imports */
#include "structs.h"
#include "../utilities/_utilities.h"

/* header */
#include "_objects.h"
#include "clock.h"

/* private functions */
double 
CLOCK_Get_MyTime();
double 
CLOCK_Get_EaselTime();
/* hack to get esl's realtime function? */
static double 
stopwatch_getRealTime(void);

/*! FUNCTION:  CLOCK_Create()
 *  SYNOPSIS:  Create new CLOCK object and returns pointer.
 *  RETURN:    Pointer to CLOCK object.
 */
CLOCK* 
CLOCK_Create()
{
   CLOCK*      cl       = NULL;   
   const int   min_size = 16; 

   cl = (CLOCK*) ERROR_malloc( sizeof(CLOCK) );

   cl->start   = 0;
   cl->stop    = 0;
   cl->stamps  = VECTOR_DBL_Create_by_Size( min_size );

   #if ( CLOCK_TYPE == CLOCK_MYCLOCK )
   {
      cl->program_start = CLOCK_Get_Time( cl );
   }
   #elif ( CLOCK_TYPE == CLOCK_EASEL )
   {
      cl->esl_stopwatch = esl_stopwatch_Create();
      esl_stopwatch_Start( cl->esl_stopwatch );
      cl->program_start = 0.0f;
      cl->esl_total     = 0.0f; 
   }
   #endif

   return cl;
}

/*! FUNCTION:  CLOCK_Destroy()
 *  SYNOPSIS:  Destroy CLOCK object.
 *  RETURN:    NULL pointer.
 */
CLOCK* 
CLOCK_Destroy( CLOCK* cl )
{
   if (cl == NULL) return NULL;

   #if ( CLOCK_TYPE == CLOCK_EASEL )
   {
      esl_stopwatch_Destroy( cl->esl_stopwatch );
   }
   #endif
   VECTOR_DBL_Destroy( cl->stamps );
   ERROR_free(cl);

   return NULL;
}

/*! FUNCTION:  CLOCK_Start()
 *  SYNOPSIS:  Set start point for <cl> clock stopwatch timer.
 *  RETURN:    Start time.
 */
double 
CLOCK_Start( CLOCK* cl )
{
   cl->start = CLOCK_Get_Time( cl );

   return cl->start;
}

/*! FUNCTION:  CLOCK_Stop()
 *  SYNOPSIS:  Set stop point for <cl> clock timer.
 *  RETURN:    Stop time.
 */
double 
CLOCK_Stop( CLOCK* cl )
{
   cl->stop = CLOCK_Get_Time( cl );

   return cl->stop;
}

/*! FUNCTION:  CLOCK_Duration()
 *  SYNOPSIS:  Get time <duration> between CLOCK_Start() and CLOCK_Stop().
 *  RETURN:    Duration time.
 */
double 
CLOCK_Duration( CLOCK* cl )
{
   cl->duration = cl->stop - cl->start;

   return (double)(cl->duration);
}

/*! FUNCTION:  CLOCK_Get_Diff()
 *  SYNOPSIS:  Get time <duration> between <start> and <stop>, acquired from CLOCK_Get_Time().
 *  RETURN:    Duration time.
 */
double 
CLOCK_Get_Diff(   CLOCK*   cl,
                  double   start,
                  double   stop )
{
   float duration;

   duration = stop - start;

   return duration;
}

/*! FUNCTION:  CLOCK_Get_Time()
 *  SYNOPSIS:  Get current time.
 *  RETURN:    Current time.
 */
double 
CLOCK_Get_Time( CLOCK* cl )
{
   double time;
   #if   ( CLOCK_TYPE == CLOCK_MYCLOCK )
   {
      time = CLOCK_Get_MyTime();
   }
   #elif ( CLOCK_TYPE == CLOCK_EASEL )
   {
      time = CLOCK_Get_EaselTime( cl );
   }
   #endif

   return time;
}

/*! FUNCTION:  CLOCK_Get_MyTime()
 *  SYNOPSIS:  Get current time (using simple method).
 *  RETURN:    Current time.
 */
double 
CLOCK_Get_MyTime()
{
   long  ltime = clock();
   float ftime = (double)ltime / (double)CLOCKS_PER_SEC;
   return ftime;
}

/*! FUNCTION:  CLOCK_Get_EaselTime()
 *  SYNOPSIS:  Capture time (from Easel stopwatch implementation)
 *  RETURN:    Current time.
 */
double 
CLOCK_Get_EaselTime( CLOCK* cl )
{
   double time;

   /* we briefly stop the stopwatch, add that time elapsed to a running total, then restart watch */
   esl_stopwatch_Stop( cl->esl_stopwatch );
   cl->esl_total += esl_stopwatch_GetElapsed( cl->esl_stopwatch );
   esl_stopwatch_Start( cl->esl_stopwatch );
   time = cl->esl_total;

   return time;
}

/*! FUNCTION:  CLOCK_Get_DateTimeString()
 *  SYNOPSIS:  Get string representation of <time>, stores in <str_buf>.
 *             <str_buf> must be allocated by user.
 *  RETURN:    Pointer to <str_buf>.
 */
char* 
CLOCK_Get_DateTimeString(  double   time,
                           char*    str_buf )
{
      // time_t   t     = time(NULL);
      // struct   tm tm = *localtime(&t);
      // char*    tm_str[256];
      // sprintf( tm_str "%d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
      // return tm_str;
   return NULL;
}

/*! FUNCTION:  CLOCK_Get_Total_Runtime()
 *  SYNOPSIS:  Get lifetime, the total time passed since clock created.
 *  RETURN:    Pointer to <str_buf>.
 */
double 
CLOCK_Get_Total_Runtime( CLOCK* cl )
{
   double current_time  = CLOCK_Get_Time( cl );
   double duration      = (cl->stop - cl->start);
   return (double)(cl->duration);
}

/*! FUNCTION:  CLOCK_Delay()
 *  SYNOPSIS:  Delay (sleep) program for <sec> seconds.
 */
STATUS_FLAG
CLOCK_Delay( int sec )
{
   STATUS_FLAG status = sleep(sec);
   return status;
}

/*! FUNCTION:  CLOCK_UnitTest()
 *  SYNOPSIS:  UnitTest for CLOCK object.
 */
STATUS_FLAG
CLOCK_UnitTest( )
{
   
}
