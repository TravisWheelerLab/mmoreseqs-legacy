/*******************************************************************************
 *  @file clock.h
 *  @brief CLOCK object.
 *
 *  @author Dave Rich
 *  @bug Lots.
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


CLOCK* clock_Create();

void clock_Destroy(CLOCK*cl);

time_t clock_Start(CLOCK*cl);

time_t clock_Stop(CLOCK*cl);

time_t clock_pTicks(CLOCK*cl, char*str);

time_t clock_Ticks(CLOCK*cl);

float ticks_to_msec(time_t t);

#endif /* _CLOCK_H */
