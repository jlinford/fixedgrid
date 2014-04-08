/*
 *  gpu_host.h
 *  add2matrix
 *
 *  Created by John Linford on 4/5/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __GPU_HOST_H__
#define __GPU_HOST_H__

#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

#include "fixedgrid.h"
#include "cuda_grid.h"

#define CU_SAFE_CALL_NO_SYNC( call ) do {                                    \
    CUresult err = call;                                                     \
    if( CUDA_SUCCESS != err) {                                               \
        fprintf(stderr, "Cuda driver error %x in file '%s' in line %i.\n",   \
                err, __FILE__, __LINE__ );                                   \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#define CU_SAFE_CALL( call ) do {                                            \
    CU_SAFE_CALL_NO_SYNC(call);                                              \
    CUresult err = cuCtxSynchronize();                                       \
    if( CUDA_SUCCESS != err) {                                               \
        fprintf(stderr, "Cuda driver error %x in file '%s' in line %i.\n",   \
                err, __FILE__, __LINE__ );                                   \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

/* Globals */
extern CUdevice cu_device_global;
extern CUcontext cu_context_global;
extern CUmodule cu_module_global;

void init_device(void);
CUresult init_cuda_driver(const char *modPath);
CUresult get_device_function(const char* func_name, CUfunction* fptr);
CUresult init_discretization_kernel(fixedgrid_t* G);
#endif