/*
 *  timer.c
 *  
 *  Common timer functionality
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "timer.h"

float elapsed_time()
{
	static int sec = -1;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	if(sec < 0) sec = tv.tv_sec;
	return (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}

inline void timer_start(volatile stopwatch_t* t)
{
	t->start = elapsed_time();
}

inline void timer_stop(volatile stopwatch_t* t)
{
    t->elapsed += elapsed_time() - t->start;
}

void metrics_init(volatile metrics_t* m, char* name)
{
    int i;
    
    stopwatch_t* tptr = (stopwatch_t*)m;
    
    for(i=0; i<NUM_TIMERS; i++, tptr++)
    {
        tptr->start = 0.0;
        tptr->elapsed = 0.0;
    }
    
    sprintf(m->name, "%s\0", name);
}

inline int64_t year2sec(int32_t years)
{
    return years * 31556926;
}

inline int32_t day2sec(int32_t days)
{
    return days * 86400;
}

inline int32_t hour2sec(int32_t hours)
{
    return hours * 3600;
}

inline int32_t minute2sec(int32_t minutes)
{
    return minutes * 60;
}

inline int32_t sec2year(int64_t seconds)
{
    return seconds / 31556926;
}

inline int32_t sec2day(int32_t seconds)
{
    return seconds / 86400;
}

inline int32_t sec2hour(int32_t seconds)
{
    return seconds / 3600;
}

inline int32_t sec2minute(int32_t seconds)
{
    return seconds / 60;
}
