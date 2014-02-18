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

/**************************************************
 * Includes                                       *
 **************************************************/

#include "common.h"

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
