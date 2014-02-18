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
    double conc[NSPEC][NROWS][NCOLS];
    
    /* Wind vector field */
    double wind_u[NROWS][NCOLS];
    double wind_v[NROWS][NCOLS];
    
    /* Diffusion tensor field */
    double diff[NROWS][NCOLS];    
    
    /* Time (seconds) */
    double end_time;
    double dt;
    uint32_t steps;
    
    /* Environment data */
    double temp;

    uint32_t nprocs;
    
    /* Metrics */
    metrics_t metrics;
    
} fixedgrid_t;


#endif
