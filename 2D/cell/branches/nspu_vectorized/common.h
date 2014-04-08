/*
 *  common.h
 *  
 *  Definitions needed by both PPU and SPU code.
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __COMMON_H__
#define __COMMON_H__

/**************************************************
 * Includes                                       *
 **************************************************/

#include <stdint.h>

/**************************************************
 * Macros                                         *
 **************************************************/

#define TRUE  1
#define FALSE 0

#define SPE_STATUS_WAITING  0
#define SPE_STATUS_ROW      1
#define SPE_STATUS_COL      2
#define SPE_STATUS_STOPPED  100

#define NUM_TIMERS 9

/**************************************************
 * Data types                                     *
 **************************************************/

/* Boolean type */
typedef short bool;

/* SPE program arguments */
typedef struct spe_argv
{
    union
    {
        double   dbl;
        uint64_t u64;
        float    flt[2];
        uint32_t u32[2];
    } 
    arg[16];
} spe_argv_t;

/* SPE program environment variables */
typedef struct spe_env
{
    uint32_t speid;
    uint32_t nprocs;
    uint64_t metptr;
} spe_env_t;

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

#endif