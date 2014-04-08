/*
 *  discretize.c
 *  
 *  Transport module.
 *  Main kernel.
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include "discretize.h"
#include "timer.h"
#include "common.h"
#include "params.h"
#include "cuda_grid.h"

typedef struct neighbors
{
    real_t left[2];
    real_t right[2];
} neighbors_t;

/**
 * Upwinded advection/diffusion equation
 */
extern "C" __device__
real_t advec_diff(real_t cell_size,
                  real_t c2l, real_t w2l, real_t d2l, 
                  real_t c1l, real_t w1l, real_t d1l, 
                  real_t   c, real_t   w, real_t   d, 
                  real_t c1r, real_t w1r, real_t d1r, 
                  real_t c2r, real_t w2r, real_t d2r)
{
    real_t wind, diff_term, advec_term;
    real_t advec_termR, advec_termL;

    wind = (w1l + w) / 2.0;
    if(wind >= 0.0) advec_termL = (1.0/6.0) * ( -c2l + 5.0*c1l + 2.0*c );
    else advec_termL = (1.0/6.0) * ( 2.0*c1l + 5.0*c - c1r );
    advec_termL *= wind;
    wind = (w1r + w) / 2.0;
    if(wind >= 0.0) advec_termR = (1.0/6.0) * ( -c1l + 5.0*c + 2.0*c1r );
    else advec_termR = (1.0/6.0) * ( 2.0*c + 5.0*c1r - c2r );
    advec_termR *= wind;
    advec_term = (advec_termL - advec_termR) / cell_size;
    diff_term = ( ((d1l+d)/2)*(c1l-c) - ((d+d1r)/2)*(c-c1r) ) / (cell_size * cell_size);
    return advec_term + diff_term;
}

/**
 * Collects neighborhood data into simple struct
 */
extern "C" __device__
void get_x_neighbors(neighbors_t* n, real_t field[NZ][NY][NX])
{
    int x = blockIdx.x*BLOCK_X + threadIdx.x;
    int y = blockIdx.y*BLOCK_Y + threadIdx.y;
    int z = blockIdx.z*BLOCK_Z + threadIdx.z;
    
    if(blockIdx.x == 0 && threadIdx.x == 0)
    {
        n->left[0] = field[z][y][NX-2];
        n->left[1] = field[z][y][NX-1];
    }
    if(blockIdx.x == 0 && threadIdx.x == 1)
    {
        n->left[0] = field[z][y][NX-1];
        n->left[1] = field[z][y][0];
    }
    if(blockIdx.x > 0 || threadIdx.x > 1)
    {
        n->left[0] = field[z][y][x-2];
        n->left[1] = field[z][y][x-1];
    }
    
    if(blockIdx.x == gridDim.x-1 && threadIdx.x == BLOCK_X-1)
    {
        n->right[0] = field[z][y][1];
        n->right[1] = field[z][y][0];
    }
    if(blockIdx.x == gridDim.x-1 && threadIdx.x == BLOCK_X-2)
    {
        n->right[0] = field[z][y][0];
        n->right[1] = field[z][y][NX-1];
    }
    if(blockIdx.x < gridDim.x-1 || threadIdx.x < BLOCK_X-2)
    {
        n->right[0] = field[z][y][x+2];
        n->right[1] = field[z][y][x+1];
    }
}

/**
 * Collects neighborhood data into simple struct
 */
extern "C" __device__
void get_y_neighbors(neighbors_t* n, real_t field[NZ][NY][NX])
{
    int x = blockIdx.x*BLOCK_X + threadIdx.x;
    int y = blockIdx.y*BLOCK_Y + threadIdx.y;
    int z = blockIdx.z*BLOCK_Z + threadIdx.z;
    
    if(blockIdx.y == 0 && threadIdx.y == 0)
    {
        n->left[0] = field[z][NY-2][x];
        n->left[1] = field[z][NY-1][x];
    }
    if(blockIdx.y == 0 && threadIdx.y == 1)
    {
        n->left[0] = field[z][NY-1][x];
        n->left[1] = field[z][0][x];
    }
    if(blockIdx.y > 0 || threadIdx.y > 1)
    {
        n->left[0] = field[z][y-2][x];
        n->left[1] = field[z][y-1][x];
    }
    
    if(blockIdx.y == gridDim.y-1 && threadIdx.y == BLOCK_Y-1)
    {
        n->right[0] = field[z][1][x];
        n->right[1] = field[z][0][x];
    }
    if(blockIdx.y == gridDim.y-1 && threadIdx.y == BLOCK_Y-2)
    {
        n->right[0] = field[z][0][x];
        n->right[1] = field[z][NX-1][x];
    }
    if(blockIdx.y < gridDim.y-1 || threadIdx.y < BLOCK_Y-2)
    {
        n->right[0] = field[z][y+2][x];
        n->right[1] = field[z][y+1][x];
    }
}

/**
 * Collects neighborhood data into simple struct
 */
extern "C" __device__
void get_z_neighbors(neighbors_t* n, real_t field[NZ][NY][NX])
{
    int x = blockIdx.x*BLOCK_X + threadIdx.x;
    int y = blockIdx.y*BLOCK_Y + threadIdx.y;
    int z = blockIdx.z*BLOCK_Z + threadIdx.z;
    
    if(blockIdx.z == 0 && threadIdx.z == 0)
    {
        n->left[0] = field[NZ-2][y][x];
        n->left[1] = field[NZ-1][y][x];
    }
    if(blockIdx.z == 0 && threadIdx.z == 1)
    {
        n->left[0] = field[NZ-1][y][x];
        n->left[1] = field[0][y][x];
    }
    if(blockIdx.z > 0 || threadIdx.z > 1)
    {
        n->left[0] = field[z-2][y][x];
        n->left[1] = field[z-1][y][x];
    }
    
    if(blockIdx.z == gridDim.z-1 && threadIdx.z == BLOCK_Z-1)
    {
        n->right[0] = field[1][y][x];
        n->right[1] = field[0][y][x];
    }
    if(blockIdx.z == gridDim.z-1 && threadIdx.z == BLOCK_Z-2)
    {
        n->right[0] = field[0][y][x];
        n->right[1] = field[NZ-1][y][x];
    }
    if(blockIdx.z < gridDim.z-1 || threadIdx.z < BLOCK_Z-2)
    {
        n->right[0] = field[z+2][y][x];
        n->right[1] = field[z+1][y][x];
    }
}

/**
 * Initialize discretization kernel data
 */
extern "C" __global__
void discretize_init(real_t c_in[NZ][NY][NX], 
                     real_t buff[NZ][NY][NX],
                     real_t c_out[NZ][NY][NX])
{
    // Data index
    int x = blockIdx.x*BLOCK_X + threadIdx.x;
    int y = blockIdx.y*BLOCK_Y + threadIdx.y;
    int z = blockIdx.z*BLOCK_Z + threadIdx.z;
    
    buff[z][y][x] = c_out[z][y][x] = c_in[z][y][x];
}

/**
 * Finalize discretization
 */
extern "C" __global__
void discretize_final(real_t c_in[NZ][NY][NX],
                      real_t buff[NZ][NY][NX], 
                      real_t c_out[NZ][NY][NX])
{
    // Data index
    int x = blockIdx.x*BLOCK_X + threadIdx.x;
    int y = blockIdx.y*BLOCK_Y + threadIdx.y;
    int z = blockIdx.z*BLOCK_Z + threadIdx.z;

    // Average results into c_out
    c_out[z][y][x] = 0.5 * (c_out[z][y][x] + buff[z][y][x]);
    if(c_out[z][y][x] < 0.0) c_out[z][y][x] = 0.0;
    
    // Update original conc data for next discretization
    c_in[z][y][x] = c_out[z][y][x];
}

/**
 * X-discretization
 */
extern "C" __global__
void advec_diff_x(real_t cell_size, real_t dt,
                  real_t c_in[NZ][NY][NX], 
                  real_t wind[NZ][NY][NX], 
                  real_t diff[NZ][NY][NX], 
                  real_t buff[NZ][NY][NX])
{
    // Data index
    int x = blockIdx.x*BLOCK_X + threadIdx.x;
    int y = blockIdx.y*BLOCK_Y + threadIdx.y;
    int z = blockIdx.z*BLOCK_Z + threadIdx.z;
    
    // Change in conc
    real_t dcdx;
    
    // Data
    neighbors_t conc_n;
    neighbors_t wind_n;
    neighbors_t diff_n;
    
    // Prepare for discretization
    get_x_neighbors(&conc_n, c_in);
    get_x_neighbors(&wind_n, wind);
    get_x_neighbors(&diff_n, diff);
        
    // Discretize
    dcdx = advec_diff(cell_size, 
                      conc_n.left[0], wind_n.left[0], diff_n.left[0],
                      conc_n.left[1], wind_n.left[1], diff_n.left[1],
                       c_in[z][y][x],  wind[z][y][x],  diff[z][y][x],
                      conc_n.right[1], wind_n.right[1], diff_n.right[1],
                      conc_n.right[0], wind_n.right[0], diff_n.right[0]);
    buff[z][y][x] += dt*dcdx;
}

/**
 * X-discretization
 */
extern "C" __global__
void advec_diff_y(real_t cell_size, real_t dt,
                  real_t c_in[NZ][NY][NX], 
                  real_t wind[NZ][NY][NX], 
                  real_t diff[NZ][NY][NX], 
                  real_t buff[NZ][NY][NX])
{
    // Data index
    int x = blockIdx.x*BLOCK_X + threadIdx.x;
    int y = blockIdx.y*BLOCK_Y + threadIdx.y;
    int z = blockIdx.z*BLOCK_Z + threadIdx.z;
    
    // Change in conc
    real_t dcdx;
    
    // Data
    neighbors_t conc_n;
    neighbors_t wind_n;
    neighbors_t diff_n;
    
    // Prepare for discretization
    get_y_neighbors(&conc_n, c_in);
    get_y_neighbors(&wind_n, wind);
    get_y_neighbors(&diff_n, diff);
        
    // Discretize
    dcdx = advec_diff(cell_size, 
                      conc_n.left[0], wind_n.left[0], diff_n.left[0],
                      conc_n.left[1], wind_n.left[1], diff_n.left[1],
                       c_in[z][y][x],  wind[z][y][x],  diff[z][y][x],
                      conc_n.right[1], wind_n.right[1], diff_n.right[1],
                      conc_n.right[0], wind_n.right[0], diff_n.right[0]);
    buff[z][y][x] += dt*dcdx;
}

/**
 * X-discretization
 */
extern "C" __global__
void advec_diff_z(real_t cell_size, real_t dt,
                  real_t c_in[NZ][NY][NX], 
                  real_t wind[NZ][NY][NX], 
                  real_t diff[NZ][NY][NX], 
                  real_t buff[NZ][NY][NX])
{
    // Data index
    int x = blockIdx.x*BLOCK_X + threadIdx.x;
    int y = blockIdx.y*BLOCK_Y + threadIdx.y;
    int z = blockIdx.z*BLOCK_Z + threadIdx.z;
    
    // Change in conc
    real_t dcdx;
    
    // Data
    neighbors_t conc_n;
    neighbors_t wind_n;
    neighbors_t diff_n;
    
    // Prepare for discretization
    get_z_neighbors(&conc_n, c_in);
    get_z_neighbors(&wind_n, wind);
    get_z_neighbors(&diff_n, diff);
        
    // Discretize
    dcdx = advec_diff(cell_size, 
                      conc_n.left[0], wind_n.left[0], diff_n.left[0],
                      conc_n.left[1], wind_n.left[1], diff_n.left[1],
                       c_in[z][y][x],  wind[z][y][x],  diff[z][y][x],
                      conc_n.right[1], wind_n.right[1], diff_n.right[1],
                      conc_n.right[0], wind_n.right[0], diff_n.right[0]);
    buff[z][y][x] += dt*dcdx;
}