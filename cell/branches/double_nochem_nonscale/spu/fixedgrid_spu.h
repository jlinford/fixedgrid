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

/**************************************************
 * Includes                                       *
 **************************************************/

#include <stdint.h>

#include "params.h"

/**************************************************
 * Macros                                         *
 **************************************************/

#define MAX_DMA_DBL 2048

/**************************************************
 * Data types                                     *
 **************************************************/

/* Single buffer of data from main memory */
typedef struct spe_buffer
{
    uint64_t ea_base;
    uint32_t length;
	double data[MAX_DMA_DBL] __attribute__((aligned(128)));
} spe_buffer_t;

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
    uint32_t length;
    dma_list_element_t data[NROWS] __attribute__((aligned(128)));
} spe_dma_list_t;

#endif