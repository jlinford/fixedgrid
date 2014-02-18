/*
 *  chemistry.h
 *  fixedgrid
 *
 *  Created by John Linford on 6/18/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __CHEMISTRY_H__
#define __CHEMISTRY_H__

/**************************************************
 * Includes                                       *
 **************************************************/

#include <stdint.h>
#include "fixedgrid_spu.h"

/**************************************************
 * Data types                                     *
 **************************************************/

/* DMA list */
typedef struct spe_chem_dma_list
{
    dma_list_element_t data[NSPEC*VECTOR_LENGTH] __attribute__((aligned(128)));
    uint32_t length;
} spe_chem_dma_list_t;

/* Single buffer of data from main memory */
typedef struct spe_chem_buffer
{
    real_t data[NSPEC*VECTOR_LENGTH] __attribute__((aligned(128)));
    uint64_t ea_base;
    uint32_t length;
} spe_chem_buffer_t;

/**************************************************
 * Function prototypes                            *
 **************************************************/

void saprc99_chem_block(uint64_t argvp);

#endif

