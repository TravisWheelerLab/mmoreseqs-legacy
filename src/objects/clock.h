/*******************************************************************************
 *  FILE:      clok.c
 *  PURPOSE:   CLOCK Object    
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Unknown.
 *******************************************************************************/

#ifndef _CLOCK_H
#define _CLOCK_H

/*
 *  FUNCTION:  CLOCK_Create()
 *  SYNOPSIS:  Create new CLOCK object and returns pointer.
 */
CLOCK* CLOCK_Create();

/*
 *  FUNCTION:  CLOCK_Destroy()
 *  SYNOPSIS:  Destroy CLOCK object.
 */
void CLOCK_Destroy( CLOCK* cl );

/* start timer */
double CLOCK_Start( CLOCK* cl );

/* stop timer  */
double CLOCK_Stop( CLOCK* cl );

/* Convert duration from ticks to msec */
float CLOCK_Secs( CLOCK* cl );


/* capture time based on system (pulled directly from easel) */
double CLOCK_Get_RealTime(void);

#endif /* _CLOCK_H */
