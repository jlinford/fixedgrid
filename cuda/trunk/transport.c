/*
 *  transport.c
 *  fixedgrid_serial
 *
 *  Created by John Linford on 6/23/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include "cuda_host.h"
#include "transport.h"
#include "discretize.h"

/**
 * Discretize conc data on device along x dimension 
 */
void discretize_all_x(fixedgrid_t* G, real_t dt)
{
#if DO_X_DISCRET == 1
    timer_start(&G->metrics.x_discret);
    
    // Configure parameters
    CU_SAFE_CALL(cuParamSetf(G->advec_diff_x, 1*sizeof(real_t), dt));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_x, 2*sizeof(real_t), G->dev_conc));
    
    // Copy atmospheric data to device
    CU_SAFE_CALL(cuMemcpyHtoD(G->dev_wind, G->wind_u, NX*NY*NZ*sizeof(real_t)));
    CU_SAFE_CALL(cuMemcpyHtoD(G->dev_diff, G->diff_h, NX*NY*NZ*sizeof(real_t)));
    
    // Run kernel functions
    CU_SAFE_CALL(cuLaunchGrid(G->discretize_init, NX/BLOCK_X, NY/BLOCK_Y));
    CU_SAFE_CALL(cuLaunchGrid(G->advec_diff_x, NX/BLOCK_X, NY/BLOCK_Y));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_x, 2*sizeof(real_t), G->dev_buff));
    CU_SAFE_CALL(cuLaunchGrid(G->advec_diff_x, NX/BLOCK_X, NY/BLOCK_Y));
    CU_SAFE_CALL(cuLaunchGrid(G->discretize_final, NX/BLOCK_X, NY/BLOCK_Y));
    
    timer_stop(&G->metrics.x_discret);
#endif
}

/**
 * Discretize conc data on device along y dimension 
 */
void discretize_all_y(fixedgrid_t* G, real_t dt)
{
#if DO_Y_DISCRET == 1
    timer_start(&G->metrics.y_discret);
    
    // Configure parameters
    CU_SAFE_CALL(cuParamSetf(G->advec_diff_y, 1*sizeof(real_t), dt));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_y, 2*sizeof(real_t), G->dev_conc));
    
    // Copy atmospheric data to device
    CU_SAFE_CALL(cuMemcpyHtoD(G->dev_wind, G->wind_v, NX*NY*NZ*sizeof(real_t)));
    CU_SAFE_CALL(cuMemcpyHtoD(G->dev_diff, G->diff_h, NX*NY*NZ*sizeof(real_t)));
    
    // Run kernel functions
    CU_SAFE_CALL(cuLaunchGrid(G->discretize_init, NX/BLOCK_X, NY/BLOCK_Y));
    CU_SAFE_CALL(cuLaunchGrid(G->advec_diff_y, NX/BLOCK_X, NY/BLOCK_Y));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_y, 2*sizeof(real_t), G->dev_buff));
    CU_SAFE_CALL(cuLaunchGrid(G->advec_diff_y, NX/BLOCK_X, NY/BLOCK_Y));
    CU_SAFE_CALL(cuLaunchGrid(G->discretize_final, NX/BLOCK_X, NY/BLOCK_Y));
    
    timer_stop(&G->metrics.y_discret);
#endif    
}

/**
 * Discretize conc data on device along z dimension 
 */
void discretize_all_z(fixedgrid_t* G, real_t dt)
{
#if DO_Z_DISCRET == 1
    timer_start(&G->metrics.z_discret);
    
    // Configure parameters
    CU_SAFE_CALL(cuParamSetf(G->advec_diff_z, 1*sizeof(real_t), dt));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_z, 2*sizeof(real_t), G->dev_conc));
    
    // Copy atmospheric data to device
    CU_SAFE_CALL(cuMemcpyHtoD(G->dev_wind, G->wind_w, NX*NY*NZ*sizeof(real_t)));
    CU_SAFE_CALL(cuMemcpyHtoD(G->dev_diff, G->diff_v, NX*NY*NZ*sizeof(real_t)));
    
    // Run kernel functions
    CU_SAFE_CALL(cuLaunchGrid(G->discretize_init, NX/BLOCK_X, NY/BLOCK_Y));
    CU_SAFE_CALL(cuLaunchGrid(G->advec_diff_z, NX/BLOCK_X, NY/BLOCK_Y));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_z, 2*sizeof(real_t), G->dev_buff));
    CU_SAFE_CALL(cuLaunchGrid(G->advec_diff_z, NX/BLOCK_X, NY/BLOCK_Y));
    CU_SAFE_CALL(cuLaunchGrid(G->discretize_final, NX/BLOCK_X, NY/BLOCK_Y));
    
    timer_stop(&G->metrics.z_discret);
#endif    
}

