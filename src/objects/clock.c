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
#include "utilities.h"
#include "objects.h"

/* header */
#include "clock.h"


double CLOCK_Get_RealTime();

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

   cl->stamps  = (double*) malloc( sizeof(double) * min_size );
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
   if (cl == NULL) return;

   free(cl->stamps);
   free(cl);
   cl = NULL;
}

/* Set the Stopwatch */
double CLOCK_Start(CLOCK*cl)
{
   cl->start = CLOCK_Get_RealTime();
   return cl->start;
}

/* Stop the Stopwatch */
double CLOCK_Stop(CLOCK*cl)
{
   cl->stop = CLOCK_Get_RealTime();
   return cl->stop;
}

/* Convert duration from ticks to msec */
float CLOCK_Secs(CLOCK*cl)
{
   cl->duration = cl->stop - cl->start;
   return (float)(cl->duration);
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
   return (double)tm.tv_sec + (double)tm.tv_usec / 1000000.0;
#else
   return -1.0;      /* Failed. */
#endif
}
