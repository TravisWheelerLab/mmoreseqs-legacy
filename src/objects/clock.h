/*******************************************************************************
 *  FILE:      clock.c
 *  PURPOSE:   CLOCK Object    
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Unknown.
 *******************************************************************************/

#ifndef _CLOCK_H
#define _CLOCK_H

// typedef struct {
//    time_t   start;
//    time_t   stop;
//    time_t   duration;

//    float    N;
//    float    Nalloc;
//    float*   stamps;
// } CLOCK;


CLOCK* CLOCK_Create();

void CLOCK_Destroy(CLOCK*cl);

time_t CLOCK_Start(CLOCK*cl);

time_t CLOCK_Stop(CLOCK*cl);

time_t CLOCK_pTicks(CLOCK*cl, char*str);

time_t CLOCK_Ticks(CLOCK*cl);

float ticks_to_msec(time_t t);

#endif /* _CLOCK_H */
