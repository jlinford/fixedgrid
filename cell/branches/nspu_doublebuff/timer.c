#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "timer.h"

/* Global timers */
static fixedgrid_timer_t timers[MAX_TIMERS];

static char *timer_names[MAX_TIMERS] = 
{
    "Wallclock   ",
    "Array Init  ",
    "Array Copy  ",
    "File I/O    ",
    "Discretize  ",
    "Advec / Diff",
    "Row discret.",
    "Col discret."
};

double elapsed_time()
{
	static int sec = -1;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	if(sec < 0) sec = tv.tv_sec;
	return (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}

inline void timer_clear(int t)
{
	timers[t].elapsed = 0.0;
}

inline void timer_start(int t)
{
	timers[t].start = elapsed_time();
}

inline void timer_stop(int t)
{
	timers[t].elapsed += elapsed_time() - timers[t].start;
}

inline double timer_read(int t)
{
	return timers[t].elapsed;
}

void print_timer_summary(char *header)
{
    int i;
    
    printf("\n%s\n", header);
    for(i=0; i<MAX_TIMERS; i++)
    {
        printf("  %d: %s: %f\n", i, timer_names[i], timer_read(i));
    }

}

inline long int year2sec(int years)
{
    return years * 31556926;
}

inline int day2sec(int days)
{
    return days * 86400;
}

inline int hour2sec(int hours)
{
    return hours * 3600;
}

inline int minute2sec(int minutes)
{
    return minutes * 60;
}

inline int sec2year(long int seconds)
{
    return seconds / 31556926;
}

inline int sec2day(int seconds)
{
    return seconds / 86400;
}

inline int sec2hour(int seconds)
{
    return seconds / 3600;
}

inline int sec2minute(int seconds)
{
    return seconds / 60;
}
