/*
 *  memory.c
 *  
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <libmisc.h>

#include "memory.h"

/**
 * Allocates memory for...
 * Concentration data
 * Wind vector field
 * Diffusion tensor
 */
void allocate_global_memory(fixedgrid_t* G)
{
    /* Allocate concentration memory */
    //G->conc = (double*)malloc_align(NROWS*NCOLS*NSPEC*sizeof(double), 7);
    
    /* Allocate wind vector filed memory */
    //G->wind_u = (double*)malloc_align(NROWS*NCOLS*sizeof(double), 7);
    //G->wind_v = (double*)malloc_align(NROWS*NCOLS*sizeof(double), 7);
    
    /* Allocate diffusion tensor memory */
    //G->diff = (double*)malloc_align(NROWS*NCOLS*sizeof(double), 7);    
}

/**
 * Frees memory allocated for...
 * Concentration data
 * Wind vector field
 * Diffusion tensor
 */
void free_global_memory(fixedgrid_t* G)
{
    //free_align((void*)G->conc);
    //free_align((void*)G->wind_u);
    //free_align((void*)G->wind_v);
    //free_align((void*)G->diff);
}
