/*
 *  cuda_grid.h
 *  fixedgrid
 *
 *  Created by John Linford on 8/1/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __CUDA_GRID_H__
#define __CUDA_GRID_H__

/* Set these according to the number of processing
 * elements available on the GPU:
 * BLOCK_X * BLOCK_Y * BLOCK_Z <= # PEs
 * Otherwise CUDA error code 1.
 */
#define BLOCK_X 8
#define BLOCK_Y 8
#define BLOCK_Z 4

#endif