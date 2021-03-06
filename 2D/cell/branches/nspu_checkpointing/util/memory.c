/*
 *  memory.c
 *  
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

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
    //G->conc = (value_t*)malloc_align(NROWS*NCOLS*NSPEC*sizeof(value_t), 7);
    
    /* Allocate wind vector filed memory */
    //G->wind_u = (value_t*)malloc_align(NROWS*NCOLS*sizeof(value_t), 7);
    //G->wind_v = (value_t*)malloc_align(NROWS*NCOLS*sizeof(value_t), 7);
    
    /* Allocate diffusion tensor memory */
    //G->diff = (value_t*)malloc_align(NROWS*NCOLS*sizeof(value_t), 7);    
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
