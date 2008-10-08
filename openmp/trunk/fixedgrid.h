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

#define conc(x, y, z, s) __conc[s][z][y][x]
#define wind_u(x, y, z)  __wind_u[z][y][x]
#define wind_v(x, y, z)  __wind_v[z][y][x]
#define wind_w(x, y, z)  __wind_w[z][y][x]
#define diff_h(x, y, z)  __diff_h[z][y][x]
#define diff_v(x, y, z)  __diff_v[z][y][x]
#define   temp(x, y, z)    __temp[z][y][x]

/**************************************************
 * Data types                                     *
 **************************************************/

typedef short bool;

/* Program state (global variables) */
typedef struct fixedgrid
{
    /* Concentration field */
    real_t __conc[NSPEC][NZ][NY][NX];
    
    /* Wind vector field */
    real_t __wind_u[NZ][NY][NX];
    real_t __wind_v[NZ][NY][NX];
    real_t __wind_w[NZ][NY][NX];
    
    /* Diffusion tensor field */
    real_t __diff_h[NZ][NY][NX];
    real_t __diff_v[NZ][NY][NX];
    
    /* Temperature field */
    real_t __temp[NZ][NY][NX];
    
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
