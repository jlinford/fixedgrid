/*
 *  timer.c
 *  
 *  Common timer functionality
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <string.h>
#include <sys/time.h>

#include "timer.h"

float elapsed_time()
{
	static int sec = -1;
	struct timeval tv;
	gettimeofday(&tv, 0);
	if(sec < 0) sec = tv.tv_sec;
	return (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
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
    
    strncpy(m->name, name, 56);
}
