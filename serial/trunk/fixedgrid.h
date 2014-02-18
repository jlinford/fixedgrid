/*
 *  fixedgrid.h
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

#include "params.h"
#include "timer.h"

/**************************************************
 * Macros                                         *
 **************************************************/

#define TRUE  1
#define FALSE 0

/**************************************************
 * Data types                                     *
 **************************************************/

typedef short bool;

/* Program state (global variables) */
typedef struct fixedgrid
{
    /* Concentration field */
    real_t conc[NSPEC][NROWS][NCOLS];
    
    /* Wind vector field */
    real_t wind_u[NROWS][NCOLS];
    real_t wind_v[NROWS][NCOLS];
    
    /* Diffusion tensor field */
    real_t diff[NROWS][NCOLS];    
    
    /* Time (seconds) */
    real_t time;
    real_t tstart;
    real_t tend;
    real_t dt;
    
    /* Parallelization */
    /* This is always == 1 for serial code */
    uint32_t nprocs;
    
    /* Metrics */
    metrics_t metrics;
    
} fixedgrid_t;


#endif
