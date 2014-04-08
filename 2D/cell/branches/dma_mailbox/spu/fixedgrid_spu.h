/*
 *  fixedgrid_spu.h
 *  
 *  Prototype chemical transport model.
 *  SPE header file.
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __FIXEDGRID_SPU_H__
#define __FIXEDGRID_SPU_H__

#include <stdio.h>

/**************************************************
 * Includes                                       *
 **************************************************/

#include <malloc_align.h>
#include <free_align.h>
#include <spu_mfcio.h>
#include <stdint.h>
#include "params.h"
#include "common.h"
#include "timer.h"

/**************************************************
 * Data types                                     *
 **************************************************/

/* DMA list element */
typedef struct dma_list_element
{
    union
    {
        unsigned int all32;
        struct
        {
            unsigned nbytes: 31;
            unsigned stall: 1;
        } bits;
    } size;
    unsigned int ea_low;
} dma_list_element_t;

/* DMA list */
typedef struct spe_dma_list
{
    dma_list_element_t* data;
    uint32_t length;
} spe_dma_list_t;

/* Single buffer of data from main memory */
typedef struct spe_buffer
{
    real_t* data;
    uint64_t ea_base;
    uint32_t length;
} spe_buffer_t;


/**************************************************
 * Globals                                        *
 **************************************************/

extern volatile metrics_t metrics;      /* Metrics */
extern volatile spe_argv_t argv;        /* Program arguments */
extern volatile spe_envv_t envv;        /* Program environment variables */
extern volatile spe_status_t status;    /* This SPU's status code */

/**************************************************
 * Inline fuctions                                *
 **************************************************/

/**
 * Waits until DMA requests with given tag-mask complete.
 */
static inline void wait_for_dma(uint32_t mask)
{
    mfc_write_tag_mask(1<<mask);
    spu_mfcstat(MFC_TAG_UPDATE_ALL);
}

/**
 * Set status code
 */
static inline void set_status(uint32_t s)
{
    status.value = s;
    mfc_put(&status, envv.ea_status, 16, 5, 0, 0);
}

/**
 * Allocates a DMA list on the heap
 */
static inline void allocate_dma_list(spe_dma_list_t* list, uint32_t length)
{
    list->data = (dma_list_element_t*)_malloc_align(length*sizeof(dma_list_element_t), 7);
}

/**
 * Frees a DMA list from the heap
 */
static inline void free_dma_list(spe_dma_list_t* list)
{
    _free_align(list->data);
}

/**
 * Allocates a buffer on the heap
 */
static inline void allocate_buffer(spe_buffer_t* buffer, uint32_t length)
{
    buffer->data = (real_t*)_malloc_align(length*sizeof(real_t), 7);
}

/**
 * Frees a buffer from the heap
 */
static inline void free_buffer(spe_buffer_t* buffer)
{
    _free_align(buffer->data);
}


#endif
