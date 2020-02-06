/*******************************************************************************
 *  @file clock.h
 *  @brief Clock object.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _CLOCK_H
#define _CLOCK_H

typedef struct {
   time_t start;
   time_t stop;
   time_t duration;
} CLOCK;

CLOCK* clock_Create();

time_t clock_Start(CLOCK* cl);

time_t clock_Stop(CLOCK* cl);

#endif /* _CLOCK_H */
