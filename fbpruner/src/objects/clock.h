/*******************************************************************************
 *  FILE:      clok.c
 *  PURPOSE:   CLOCK Object    
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _CLOCK_H
#define _CLOCK_H

/*! FUNCTION:  CLOCK_Create()
 *  SYNOPSIS:  Create new CLOCK object and returns pointer.
 *  RETURN:    Pointer to CLOCK object.
 */
CLOCK* 
CLOCK_Create();

/*! FUNCTION:  CLOCK_Destroy()
 *  SYNOPSIS:  Destroy CLOCK object.
 *  RETURN:    NULL pointer.
 */
CLOCK* 
CLOCK_Destroy( CLOCK* cl );

/*! FUNCTION:  CLOCK_Start()
 *  SYNOPSIS:  Set start point for <cl> clock stopwatch timer.
 *  RETURN:    Start time.
 */
double 
CLOCK_Start( CLOCK* cl );

/*! FUNCTION:  CLOCK_Stop()
 *  SYNOPSIS:  Set stop point for <cl> clock timer.
 *  RETURN:    Stop time.
 */
double 
CLOCK_Stop( CLOCK* cl );

/*! FUNCTION:  CLOCK_Duration()
 *  SYNOPSIS:  Get time <duration> between CLOCK_Start() and CLOCK_Stop().
 *  RETURN:    Duration time.
 */
double 
CLOCK_Duration( CLOCK* cl );

/*! FUNCTION:  CLOCK_Get_Diff()
 *  SYNOPSIS:  Get time <duration> between <start> and <stop>, acquired from CLOCK_Get_Time().
 *  RETURN:    Duration time.
 */
double 
CLOCK_Get_Diff(   CLOCK*   cl,
                  double   start,
                  double   stop );

/*! FUNCTION:  CLOCK_Get_Time()
 *  SYNOPSIS:  Get current time.
 *  RETURN:    Current time.
 */
double 
CLOCK_Get_Time( CLOCK* cl );

/*! FUNCTION:  CLOCK_Get_MyTime()
 *  SYNOPSIS:  Get current time (using my simple method)
 *  RETURN:    Current time.
 */
double 
CLOCK_Get_MyTime();

/*! FUNCTION:  CLOCK_Get_EaselTime()
 *  SYNOPSIS:  Capture time based on system (using easel method)
 *  RETURN:    Current time.
 */
double 
CLOCK_Get_EaselTime();

/*! FUNCTION:  CLOCK_Get_EaselTime()
 *  SYNOPSIS:  Get string representation of <time>, stores in <str_buf>.
 *             <str_buf> must be allocated by user.
 *  RETURN:    Pointer to <str_buf>.
 */
char* 
CLOCK_Get_DateTimeString(  double   time,
                           char*    str_buf );

/*! FUNCTION:  CLOCK_Get_Total_Runtime()
 *  SYNOPSIS:  Get lifetime, the total time passed since clock created.
 *  RETURN:    Pointer to <str_buf>.
 */
double 
CLOCK_Get_Total_Runtime( CLOCK* cl );

/*! FUNCTION:  CLOCK_Delay()
 *  SYNOPSIS:  Delay (sleep) program for <sec> seconds.
 */
STATUS_FLAG
CLOCK_Delay( int sec );

#endif /* _CLOCK_H */
