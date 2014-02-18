/*
 *  fixedgrid.h
 *  
 *  PPU-only definitions.
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __FIXEDGRID_H__
#define __FIXEDGRID_H__

/**************************************************
 * Includes                                       *
 **************************************************/

#include <stdint.h>
#include <pthread.h>
#include <libspe2.h>

#include "params.h"
#include "common.h"
#include "timer.h"

/**************************************************
 * Macros                                         *
 **************************************************/

#define SPE_MAX_THREADS 16

//#define conc(spec, row, col) conc[NROWS*NCOLS*spec + NCOLS*row + col]
#define conc(spec, row, col) conc[NCOLS*row + col]
#define wind_u(row, col) wind_u[row*NCOLS+col]
#define wind_v(row, col) wind_u[row*NCOLS+col]
#define diff(row, col) wind_u[row*NCOLS+col]

/**************************************************
 * Data types                                     *
 **************************************************/

/* Thread data datatype */
typedef struct spe_pthread_data 
{
    uint32_t status;
    spe_context_ptr_t speid;
    pthread_t pthread;
    
    volatile metrics_t metrics __attribute__((aligned(128)));
    volatile spe_argv_t argv __attribute__((aligned(128)));
    volatile spe_env_t env __attribute__((aligned(128)));
} spe_pthread_data_t;

/* Program state (global variables) */
typedef struct fixedgrid
{        
    /* 2D concentration data: NROWS x NCOLS x NSPEC */
    //volatile double* conc;
    volatile double conc[NROWS * NCOLS * NSPEC] __attribute__((aligned(128)));
    
    /* 2D wind field data: NROWS x NCOLS */
    //volatile double* wind_u;
    //volatile double* wind_v;
    volatile double wind_u[NROWS * NCOLS] __attribute__((aligned(128)));
    volatile double wind_v[NROWS * NCOLS] __attribute__((aligned(128)));
    
    /* 2D diffusion tensor data: NROWS x NCOLS */
    //volatile double* diff;
    volatile double diff[NROWS * NCOLS] __attribute__((aligned(128)));
    
    /* SPE threads */
    spe_pthread_data_t threads[SPE_MAX_THREADS];
    
    /* Time (seconds) */
    double end_time;
    double dt;
    uint32_t steps;

    /* Environment data */
    double temp;

    /* Parallelization */
    int nprocs;
    
    /* Metrics */
    metrics_t ppe_metrics;

} fixedgrid_t;

#endif
