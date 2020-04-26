/*******************************************************************************
 *  FILE:      clok.c
 *  PURPOSE:   CLOCK Object    
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Unknown.
 *******************************************************************************/

#ifndef _CLOCK_H
#define _CLOCK_H

/* constructor */
CLOCK* CLOCK_Create();

/* destructor */
void CLOCK_Destroy(CLOCK* cl);

/* start timer */
time_t CLOCK_Start(CLOCK* cl);

/* stop timer  */
time_t CLOCK_Stop(CLOCK* cl);

/* get duration in ticks */
time_t CLOCK_Ticks(CLOCK* cl);

/* get duration in msecs */
float CLOCK_Secs(CLOCK* cl);

/* (test) print duration in ticks */
void CLOCK_pTicks(CLOCK* 	cl, 
				  char*		str);

/* convert time in ticks to milliseconds */
float ticks_to_msec(time_t t);

#endif /* _CLOCK_H */
