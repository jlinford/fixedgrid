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

/* 
 * Macros to convert matrix coordinates (row/col) to
 * cartesian coordinates (x/y), while keeping rows
 * of data contiguous 
 */
#define conc(x, y, z, s) __conc[s][z][y][x]
#define wind_u(x, y, z) __wind_u[z][y][x]
#define wind_v(x, y, z) __wind_v[z][y][x]
#define wind_w(x, y, z) __wind_w[z][y][x]
#define diff_h(x, y, z) __diff_h[z][y][x]
#define diff_v(x, y, z) __diff_v[z][y][x]
#define   temp(x, y, z)   __temp[z][y][x]

/**************************************************
 * Data types                                     *
 **************************************************/

/* Thread data datatype */
typedef struct spe_pthread_data 
{
    spe_context_ptr_t speid;
    pthread_t pthread;
    volatile spe_status_t status  __attribute__((aligned(128)));
    volatile spe_argv_t   argv    __attribute__((aligned(128)));
    volatile spe_envv_t   envv    __attribute__((aligned(128)));
    volatile metrics_t    metrics __attribute__((aligned(128)));
} spe_pthread_data_t;

/* Program state (global variables) */
typedef struct fixedgrid
{        
    /* 3D concentration data */
    volatile real_t __conc[NSPEC][NZ][NY][NX_ALIGNED] __attribute__((aligned(128)));
    
    /* 3D wind field data */
    volatile real_t __wind_u[NZ][NY][NX_ALIGNED] __attribute__((aligned(128)));
    volatile real_t __wind_v[NZ][NY][NX_ALIGNED] __attribute__((aligned(128)));
    volatile real_t __wind_w[NZ][NY][NX_ALIGNED] __attribute__((aligned(128)));
    
    /* 3D diffusion data */
    volatile real_t __diff_h[NZ][NY][NX_ALIGNED] __attribute__((aligned(128)));
    volatile real_t __diff_v[NZ][NY][NX_ALIGNED] __attribute__((aligned(128)));
    
    /* 3D temperature data */
    volatile real_t __temp[NZ][NY][NX_ALIGNED] __attribute__((aligned(128)));
    
    /* SPE threads */
    spe_pthread_data_t threads[SPE_MAX_THREADS];
    
    /* Time (seconds) */
    real_t time;
    real_t tstart;
    real_t tend;
    real_t dt;

    /* Parallelization */
    uint32_t nprocs;
    
    /* Metrics */
    metrics_t ppe_metrics;

} fixedgrid_t;

/**************************************************
 * Inline fuctions                                *
 **************************************************/

//static inline uint32_t calc_block(fixedgrid_t* G, uint32_t dim, uint64_t id)
//{
//    uint32_t block = dim / G->nprocs;
//    if(id < dim % G->nprocs) ++block;
//    return block;
//}

/**
 * Gets the status of an spe
 */
static inline uint64_t spe_get_status(fixedgrid_t* G, uint32_t id)
{
    return G->threads[id].status.value;
}

/**
 * Sets the status of an spe
 */
static inline void spe_set_status(fixedgrid_t* G, uint32_t id, uint32_t status)
{
    G->threads[id].status.value = status;
    spe_mfcio_get(G->threads[id].speid, G->threads[id].envv.ls_status,
                  (void*)G->threads[id].envv.ea_status, 16, 5, 0, 0);
}

/**
 * Waits for one spe to have SPE_WAITING_STATUS status
 */
static inline void wait_for_spe(fixedgrid_t* G, uint32_t id)
{
    while(spe_get_status(G, id) != SPE_STATUS_WAITING) ; // Intentional wait
}

/**
 * Waits for all spes to have SPE_WAITING_STATUS status
 */
static inline void wait_all_spes(fixedgrid_t* G)
{
    uint32_t i, t;
    
    do
    {
        for(i=0, t=0; i<G->nprocs; ++i)
        {
            t += spe_get_status(G, i);
        }
    }
    while(t > 0);
}

#endif
