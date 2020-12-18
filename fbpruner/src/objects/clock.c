/*******************************************************************************
 *  FILE:      clok.c
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
#include <sys/times.h>

/* local imports */
#include "structs.h"
#include "../utilities/utilities.h"
#include "objects.h"

/* header */
#include "clock.h"


double CLOCK_Get_RealTime();

/* constructor */
CLOCK* CLOCK_Create()
{
   CLOCK*      cl       = NULL;   
   const int   min_size = 16; 

   cl = (CLOCK*) ERROR_malloc( sizeof(CLOCK) );

   cl->start   = 0;
   cl->stop    = 0;
   cl->stamps  = VECTOR_DBL_Create_by_Size( min_size );

   /* set program start time */
   cl->program_start = CLOCK_Get_RealTime();

   return cl;
}

/* destructor */
void* CLOCK_Destroy( CLOCK* cl )
{
   if (cl == NULL) return NULL;

   VECTOR_DBL_Destroy( cl->stamps );
   ERROR_free(cl);

   return NULL;
}

/* Set the Stopwatch */
double CLOCK_Start(CLOCK*cl)
{
   cl->start = CLOCK_Get_RealTime();
   return cl->start;
}

/* Stop the Stopwatch */
double CLOCK_Stop(CLOCK* cl)
{
   cl->stop = CLOCK_Get_RealTime();
   return cl->stop;
}

/* Convert duration from ticks to msec */
double CLOCK_Secs(CLOCK*cl)
{
   cl->duration = cl->stop - cl->start;
   return (double)(cl->duration);
}

/* dependencies for Get_RealTime() */
#if defined(_WIN32)
#include <Windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>  /* POSIX flags */
#include <time.h> /* clock_gettime(), time() */
#include <sys/time.h>   /* gethrtime(), gettimeofday() */

#if defined(__MACH__) && defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

#else
#error "Unable to define getRealTime( ) for an unknown OS."
#endif

/* capture time based on system (pulled directly from easel) */
double CLOCK_Get_RealTime()
{
   #if defined(_WIN32)
   FILETIME tm;
   ULONGLONG t;
#if defined(NTDDI_WIN8) && NTDDI_VERSION >= NTDDI_WIN8
   /* Windows 8, Windows Server 2012 and later. ---------------- */
   GetSystemTimePreciseAsFileTime( &tm );
#else
   /* Windows 2000 and later. ---------------------------------- */
   GetSystemTimeAsFileTime( &tm );
#endif
   t = ((ULONGLONG)tm.dwHighDateTime << 32) | (ULONGLONG)tm.dwLowDateTime;
   return (double)t / 10000000.0;

#elif (defined(__hpux) || defined(hpux)) || ((defined(__sun__) || defined(__sun) || defined(sun)) && (defined(__SVR4) || defined(__svr4__)))
   /* HP-UX, Solaris. ------------------------------------------ */
   return (double)gethrtime( ) / 1000000000.0;

#elif defined(__MACH__) && defined(__APPLE__)
   /* OSX. ----------------------------------------------------- */
   static double timeConvert = 0.0;
   if ( timeConvert == 0.0 )
   {
      mach_timebase_info_data_t timeBase;
      (void)mach_timebase_info( &timeBase );
      timeConvert = (double)timeBase.numer /
         (double)timeBase.denom /
         1000000000.0;
   }
   return (double)mach_absolute_time( ) * timeConvert;

#elif defined(_POSIX_VERSION)
   /* POSIX. --------------------------------------------------- */
#if defined(_POSIX_TIMERS) && (_POSIX_TIMERS > 0)
   {
      struct timespec ts;
#if defined(CLOCK_MONOTONIC_PRECISE)
      /* BSD. --------------------------------------------- */
      const clockid_t id = CLOCK_MONOTONIC_PRECISE;
#elif defined(CLOCK_MONOTONIC_RAW)
      /* Linux. ------------------------------------------- */
      const clockid_t id = CLOCK_MONOTONIC_RAW;
#elif defined(CLOCK_HIGHRES)
      /* Solaris. ----------------------------------------- */
      const clockid_t id = CLOCK_HIGHRES;
#elif defined(CLOCK_MONOTONIC)
      /* AIX, BSD, Linux, POSIX, Solaris. ----------------- */
      const clockid_t id = CLOCK_MONOTONIC;
#elif defined(CLOCK_REALTIME)
      /* AIX, BSD, HP-UX, Linux, POSIX. ------------------- */
      const clockid_t id = CLOCK_REALTIME;
#else
      const clockid_t id = (clockid_t)-1; /* Unknown. */
#endif /* CLOCK_* */
      if ( id != (clockid_t)-1 && clock_gettime( id, &ts ) != -1 )
         return (double)ts.tv_sec +
            (double)ts.tv_nsec / 1000000000.0;
      /* Fall thru. */
   }
#endif /* _POSIX_TIMERS */

   /* AIX, BSD, Cygwin, HP-UX, Linux, OSX, POSIX, Solaris. ----- */
   struct timeval tm;
   gettimeofday( &tm, NULL );
   return (double)tm.tv_sec + (double)tm.tv_usec / 1000000000.0;
#else
   return -1.0;      /* Failed. */
#endif
}

/* returns string containing current date and time */
char* CLOCK_Get_DateTimeString( CLOCK* cl )
{
      // time_t   t     = time(NULL);
      // struct   tm tm = *localtime(&t);
      // char*    tm_str[256];
      // sprintf( tm_str "%d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
      // return tm_str;
   return NULL;
}

/* return total runtime since clock created */
double CLOCK_Get_Total_Runtime( CLOCK* cl )
{
   double current_time  = CLOCK_Get_RealTime();
   double duration      = (cl->stop - cl->start);
   return (double)(cl->duration);
}

/* Delay for <sec> milliseconds */
int CLOCK_Delay( int milliseconds )
{

}