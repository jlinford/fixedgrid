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

#define conc(spec, row, col) conc[ROW_LENGTH*NROWS*spec + ROW_LENGTH*row + col]
#define wind_u(row, col) wind_u[row*ROW_LENGTH+col]
#define wind_v(row, col) wind_v[row*ROW_LENGTH+col]
#define diff(row, col) diff[row*ROW_LENGTH+col]

/**************************************************
 * Data types                                     *
 **************************************************/

/* Thread data datatype */
typedef struct spe_pthread_data 
{
    spe_context_ptr_t speid;
    pthread_t pthread;
    volatile spe_status_t status __attribute__((aligned(128)));
    volatile spe_argv_t argv     __attribute__((aligned(128)));
    volatile spe_envv_t envv     __attribute__((aligned(128)));
    volatile metrics_t metrics   __attribute__((aligned(128)));
} spe_pthread_data_t;

/* Program state (global variables) */
typedef struct fixedgrid
{        
    /* 2D concentration data: NROWS x NCOLS x NSPEC */
    volatile real_t conc[ROW_LENGTH * NCOLS * NSPEC] __attribute__((aligned(128)));
    
    /* 2D wind field data: NROWS x NCOLS */
    volatile real_t wind_u[ROW_LENGTH * NCOLS] __attribute__((aligned(128)));
    volatile real_t wind_v[ROW_LENGTH * NCOLS] __attribute__((aligned(128)));
    
    /* 2D diffusion tensor data: NROWS x NCOLS */
    volatile real_t diff[ROW_LENGTH * NCOLS] __attribute__((aligned(128)));
    
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
                  G->threads[id].envv.ea_status, 16, 5, 0, 0);
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
