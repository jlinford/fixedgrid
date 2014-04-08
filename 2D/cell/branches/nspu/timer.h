#ifndef __TIMER_H__
#define __TIMER_H__

#define TIMER_WALLCLOCK     0
#define TIMER_ARRAY_INIT    1
#define TIMER_ARRAY_COPY    2
#define TIMER_FIO           3
#define TIMER_DISCRETIZE    4
#define TIMER_ADVEC_DIFF    5
#define MAX_TIMERS          6

void timer_clear(int n);

void timer_start(int n);

void timer_stop(int n);

double timer_read(int n);

void initialize_timers();

void print_timer_summary();
            
long int year2sec(int years);

int day2sec(int days);

int hour2sec(int hours);

int minute2sec(int minutes);

int sec2year(long int seconds);

int sec2day(int seconds);

int sec2hour(int seconds);

int sec2minute(int seconds);

#endif
