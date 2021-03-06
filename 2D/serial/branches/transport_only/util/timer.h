/*
 *  timer.h
 *  
 *  Common timer functionality
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __TIMER_H__
#define __TIMER_H__

#include <stdint.h>

#define NUM_TIMERS 9

/* Stopwtch for gathering metrics */
typedef struct stopwatch
{
    float start;
    float elapsed;
} stopwatch_t;

/* Thread metrics */
typedef struct metrics
{
    stopwatch_t wallclock;
    stopwatch_t array_init;
    stopwatch_t array_copy;
    stopwatch_t file_io;
    stopwatch_t row_discret;
    stopwatch_t col_discret;
    stopwatch_t tot_discret;
    stopwatch_t chem;
    stopwatch_t comm;
    char name[56];
} metrics_t;

/**************************************************
 * Function Prototypes                            *
 **************************************************/

void timer_start(volatile stopwatch_t* n);

void timer_stop(volatile stopwatch_t* n);

void metrics_init(volatile metrics_t* m, char* name);

int64_t year2sec(int32_t years);

int32_t day2sec(int32_t days);

int32_t hour2sec(int32_t hours);

int32_t minute2sec(int32_t minutes);

int32_t sec2year(int64_t seconds);

int32_t sec2day(int32_t seconds);

int32_t sec2hour(int32_t seconds);

int32_t sec2minute(int32_t seconds);

#endif
