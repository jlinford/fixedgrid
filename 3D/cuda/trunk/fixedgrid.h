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
#include <cuda.h>

#include "params.h"
#include "common.h"
#include "timer.h"

/**************************************************
 * Macros                                         *
 **************************************************/

#define conc(x, y, z, s) conc[s][z][y][x]
#define wind_u(x, y, z) wind_u[z][y][x]
#define wind_v(x, y, z) wind_v[z][y][x]
#define wind_w(x, y, z) wind_w[z][y][x]
#define diff_h(x, y, z) diff_h[z][y][x]
#define diff_v(x, y, z) diff_v[z][y][x]
#define   temp(x, y, z)   temp[z][y][x]

/**************************************************
 * Data types                                     *
 **************************************************/

/* Program state (global variables) */
typedef struct fixedgrid
{        
    /* 3D concentration data */
    real_t conc[NSPEC][NZ][NY][NX];
    
    /* 3D wind field data */
    real_t wind_u[NZ][NY][NX];
    real_t wind_v[NZ][NY][NX];
    real_t wind_w[NZ][NY][NX];
    
    /* 3D diffusion data */
    real_t diff_h[NZ][NY][NX];
    real_t diff_v[NZ][NY][NX];
    
    /* 3D temperature data */
    real_t temp[NZ][NY][NX];
    
    /* Device data pointers */
    CUdeviceptr dev_conc;
    CUdeviceptr dev_wind;
    CUdeviceptr dev_diff;
    CUdeviceptr dev_buff;
    CUdeviceptr dev_conc_out;
    
    /* Device function pointers */
    CUfunction discretize_init;
    CUfunction discretize_final;
    CUfunction advec_diff_x;
    CUfunction advec_diff_y;
    CUfunction advec_diff_z;
        
    /* Time (seconds) */
    real_t time;
    real_t tstart;
    real_t tend;
    real_t dt;
    
    /* Metrics */
    metrics_t metrics;

} fixedgrid_t;

/**************************************************
 * Inline fuctions                                *
 **************************************************/

#endif
