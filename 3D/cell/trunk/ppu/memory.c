/*
 *  memory.c
 *  
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <malloc_align.h>
#include <free_align.h>

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
    //G->conc = (real_t*)_malloc_align(NX_ALIGNED*NX*NSPEC*sizeof(real_t), 7);
    
    /* Allocate wind vector filed memory */
    //G->wind_u = (real_t*)malloc_align(NY*NX*sizeof(real_t), 7);
    //G->wind_v = (real_t*)malloc_align(NY*NX*sizeof(real_t), 7);
    
    /* Allocate diffusion tensor memory */
    //G->diff = (real_t*)malloc_align(NY*NX*sizeof(real_t), 7);    
}

/**
 * Frees memory allocated for...
 * Concentration data
 * Wind vector field
 * Diffusion tensor
 */
void free_global_memory(fixedgrid_t* G)
{
    //_free_align((void*)G->conc);
    //free_align((void*)G->wind_u);
    //free_align((void*)G->wind_v);
    //free_align((void*)G->diff);
}
