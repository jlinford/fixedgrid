/*
 *  gpu_host.c
 *  add2matrix
 *
 *  Created by John Linford on 4/5/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include "cuda_host.h"
#include "fixedgrid.h"

/* Globals */
CUdevice cu_device_global;
CUcontext cu_context_global;
CUmodule cu_module_global;

void init_device(void)
{
    int dev;
    int devCount;
    int major, minor;
    CUresult err;
    
    cu_device_global = 0;
    
    err = cuInit(0);
    if (err == CUDA_SUCCESS)
    {
        CU_SAFE_CALL_NO_SYNC(cuDeviceGetCount(&devCount));
    }
    
    if (devCount == 0) 
    {
        fprintf(stderr, "No GPU device found.\n");
        exit(EXIT_FAILURE);
    }
    
    for(dev=0; dev<devCount; ++dev)
    {
        CU_SAFE_CALL_NO_SYNC(cuDeviceComputeCapability(&major, &minor, dev));
        if (major >= 1) break;
    }
    
    if (dev == devCount)
    {
        fprintf(stderr, "A GPU was found, but it does not support CUDA.\n"); 
        exit(EXIT_FAILURE);
    }
    else
    {
        CU_SAFE_CALL_NO_SYNC(cuDeviceGet(&cu_device_global, dev));
    }        
}

CUresult init_cuda_driver(const char *modPath)
{
    CUresult status = 0;
    
    init_device();
    
    status = cuCtxCreate(&cu_context_global, 0, cu_device_global);
    if(status != CUDA_SUCCESS)
    {
        cuCtxDetach(cu_context_global);
        return status;
    }

    status = cuModuleLoad(&cu_module_global, modPath);
    if(status != CUDA_SUCCESS)
    {
        cuCtxDetach(cu_context_global);
        return status;
    }
        
    return CUDA_SUCCESS;
}

CUresult get_device_function(const char* func_name, CUfunction* fptr)
{
    CUresult status = 0;
    CUfunction cu_func = NULL;

    status = cuModuleGetFunction(&cu_func, cu_module_global, func_name);
    if(status != CUDA_SUCCESS)
    {
        cuCtxDetach(cu_context_global);
        return status;
    }
    
    *fptr = cu_func;
    
    return CUDA_SUCCESS;
}

CUresult init_discretization_kernel(fixedgrid_t* G)
{
    // Get device function pointers
    CU_SAFE_CALL(get_device_function("discretize_init", &G->discretize_init));
    CU_SAFE_CALL(get_device_function("discretize_final", &G->discretize_final));
    CU_SAFE_CALL(get_device_function("advec_diff_x", &G->advec_diff_x));
    CU_SAFE_CALL(get_device_function("advec_diff_y", &G->advec_diff_y));
    CU_SAFE_CALL(get_device_function("advec_diff_z", &G->advec_diff_z));
    
    // Set block size for all functions
    CU_SAFE_CALL(cuFuncSetBlockShape(G->discretize_init, BLOCK_X, BLOCK_Y, BLOCK_Z));
    CU_SAFE_CALL(cuFuncSetBlockShape(G->discretize_final, BLOCK_X, BLOCK_Y, BLOCK_Z));
    CU_SAFE_CALL(cuFuncSetBlockShape(G->advec_diff_x, BLOCK_X, BLOCK_Y, BLOCK_Z));
    CU_SAFE_CALL(cuFuncSetBlockShape(G->advec_diff_y, BLOCK_X, BLOCK_Y, BLOCK_Z));
    CU_SAFE_CALL(cuFuncSetBlockShape(G->advec_diff_z, BLOCK_X, BLOCK_Y, BLOCK_Z));
    
    // Configure discretize_init parameters
    CU_SAFE_CALL(cuParamSeti(G->discretize_init, 0, G->dev_conc));
    CU_SAFE_CALL(cuParamSeti(G->discretize_init, 1*sizeof(real_t*), G->dev_buff));
    CU_SAFE_CALL(cuParamSeti(G->discretize_init, 2*sizeof(real_t*), G->dev_conc_out));
    CU_SAFE_CALL(cuParamSetSize(G->discretize_init, 3*sizeof(real_t*)));
    
    // Configure discretize_final parameters
    CU_SAFE_CALL(cuParamSeti(G->discretize_final, 0, G->dev_conc));
    CU_SAFE_CALL(cuParamSeti(G->discretize_final, 1*sizeof(real_t*), G->dev_buff));
    CU_SAFE_CALL(cuParamSeti(G->discretize_final, 2*sizeof(real_t*), G->dev_conc_out));
    CU_SAFE_CALL(cuParamSetSize(G->discretize_final, 3*sizeof(real_t*)));    
    
    // Configure advec_diff_x parameters
    CU_SAFE_CALL(cuParamSetf(G->advec_diff_x, 0, (real_t)DX));
    CU_SAFE_CALL(cuParamSetf(G->advec_diff_x, 1*sizeof(real_t), G->dt*0.5));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_x, 2*sizeof(real_t), G->dev_conc));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_x, 2*sizeof(real_t)+1*sizeof(real_t*), G->dev_wind));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_x, 2*sizeof(real_t)+2*sizeof(real_t*), G->dev_diff));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_x, 2*sizeof(real_t)+3*sizeof(real_t*), G->dev_buff));
    CU_SAFE_CALL(cuParamSetSize(G->advec_diff_x, 2*sizeof(real_t)+4*sizeof(real_t*)));

    // Configure advec_diff_y parameters
    CU_SAFE_CALL(cuParamSetf(G->advec_diff_y, 0, (real_t)DY));
    CU_SAFE_CALL(cuParamSetf(G->advec_diff_y, 1*sizeof(real_t), G->dt*0.5));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_y, 2*sizeof(real_t), G->dev_conc));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_y, 2*sizeof(real_t)+1*sizeof(real_t*), G->dev_wind));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_y, 2*sizeof(real_t)+2*sizeof(real_t*), G->dev_diff));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_y, 2*sizeof(real_t)+3*sizeof(real_t*), G->dev_buff));
    CU_SAFE_CALL(cuParamSetSize(G->advec_diff_y, 2*sizeof(real_t)+4*sizeof(real_t*)));

    // Configure advec_diff_z parameters
    CU_SAFE_CALL(cuParamSetf(G->advec_diff_z, 0, (real_t)DZ));
    CU_SAFE_CALL(cuParamSetf(G->advec_diff_z, 1*sizeof(real_t), G->dt));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_z, 2*sizeof(real_t), G->dev_conc));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_z, 2*sizeof(real_t)+1*sizeof(real_t*), G->dev_wind));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_z, 2*sizeof(real_t)+2*sizeof(real_t*), G->dev_diff));
    CU_SAFE_CALL(cuParamSeti(G->advec_diff_z, 2*sizeof(real_t)+3*sizeof(real_t*), G->dev_buff));
    CU_SAFE_CALL(cuParamSetSize(G->advec_diff_z, 2*sizeof(real_t)+4*sizeof(real_t*)));
    
    return CUDA_SUCCESS;
}