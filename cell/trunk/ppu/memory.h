/*
 *  memory.h
 *  
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __MEMORY_H__
#define __MEMORY_H__

/**************************************************
 * Includes                                       *
 **************************************************/

#include "fixedgrid.h"

/**************************************************
 * Function Prototypes                            *
 **************************************************/

/**
 * Allocates memory for...
 * Concentration data
 * Wind vector field
 * Diffusion tensor
 */
void allocate_global_memory(fixedgrid_t* G);

/**
 * Frees memory allocated for...
 * Concentration data
 * Wind vector field
 * Diffusion tensor
 */
void free_global_memory(fixedgrid_t* G);

#endif
