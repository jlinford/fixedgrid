#include <stdio.h>
#include <sys/time.h>

#include "timer.h"

static struct timers
{
	double start[MAX_TIMERS];
	double elapsed[MAX_TIMERS];
	char names[MAX_TIMERS][255];
} tt;

static char *timer_names[] = 
{
    "Wallclock   ",
    "Array Init  ",
    "Array Copy  ",
    "File I/O    ",
    "Discretize  ",
    "Advec / Diff"
};


double elapsed_time()
{
	static int sec = -1;
	struct timeval tv;
	gettimeofday(&tv, (void*)0);
	if(sec < 0) sec = tv.tv_sec;
	return (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}

inline void timer_clear(int n)
{
	tt.elapsed[n] = 0.0;
}

inline void timer_start(int n)
{
	tt.start[n] = elapsed_time();
}

inline void timer_stop(int n)
{
	double t, now;
	now = elapsed_time();
	t = now - tt.start[n];
	tt.elapsed[n] += t;
}

inline double timer_read(int n)
{
	return tt.elapsed[n];
}

void initialize_timers()
{
    int i;
    
    for(i=0; i<MAX_TIMERS; i++)
        timer_clear(i);
}

void print_timer_summary()
{
    int i;
    
    printf("\nTimers:\n");
    
    for(i=0; i<MAX_TIMERS; i++)
    {
        printf("  %d: %s: %f\n", i, timer_names[i], timer_read(i));
        //printf("  %d: %f\n", i, timer_read(i));
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
