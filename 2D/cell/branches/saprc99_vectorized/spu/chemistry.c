/*
 *  chemistry.c
 *  cell_double
 *
 *  Created by John Linford on 6/13/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <stdio.h>

#include "chemistry.h"
#include "fixedgrid_spu.h"
#include "saprc99_Global.h"
#include "saprc99_Monitor.h"
#include "common.h"
#include "params.h"

/* DMA lists */
spe_dma_list_t clist[2];  /* DMA list for conc data */

/* Buffered local storage */
spe_buffer_t conc[2];     /* in_x */

int Rosenbrock(volatile real_t Y[], real_t Tstart, real_t Tend,
               volatile real_t AbsTol[], volatile real_t RelTol[],
               real_t RPAR[], int IPAR[]);

/**
 * Get arguments from main memory synchronously 
 */
void get_chemistry_argv(uint64_t argvp, uint32_t* rows)
{
    timer_start(&metrics.comm);
    
    mfc_get(&argv, argvp, sizeof(spe_argv_t), 31, 0, 0);
    wait_for_dma(31);
    
    conc[0].ea_base = argv.arg[0].u64;
    conc[0].length  = NSPEC;
    
    conc[1].ea_base = conc[0].ea_base;
    conc[1].length  = conc[0].length;
    
    TIME  = argv.arg[1].dbl;
    DT    = argv.arg[2].dbl;
    *rows = argv.arg[3].u32[0];
    
    timer_stop(&metrics.comm);
}

/**
 * Fetches chemistry data from main memory
 */
void fetch_chem_buffer(uint32_t i, uint64_t ea_off)
{
    uint32_t j;
    
    /* Build DMA lists */
    timer_start(&metrics.array_copy);
    
    clist[i].length = conc[i].length;
    
    for(j=0; j<clist[i].length; ++j)
    {
        clist[i].data[j].size.all32 = 16;
        clist[i].data[j].ea_low = conc[i].ea_base + ea_off + j*ROW_SIZE*NROWS;
    }
    
    timer_stop(&metrics.array_copy);
    
    /* Fetch data */
    spu_mfcdma32(conc[i].data, (unsigned int)clist[i].data, clist[i].length*sizeof(dma_list_element_t), i, MFC_GETL_CMD);
}

/**
 * Runs KPP-generated SAPRC'99 chemical process in vector data buffer
 */
void saprc99_buffer(uint32_t i)
{
    /* Iterators */
    uint32_t j, k;
    
    /* Integration method statistics.
     * (Used to detect limit violation.) */
    static int stats[8];    
    
    /* Integration method parameters */
    real_t RPAR[20];
    int    IPAR[20];
    int    IERR;
    
    /* Initalize parameters */
    for(j=0; j<20; j++)
    {
        IPAR[j] = 0;
        RPAR[j] = 0.0;
    }
    IPAR[0] = 0;        /* non-autonomous */
    IPAR[1] = 1;        /* scalar tolerances */
    RPAR[2] = STEPMIN;  /* starting step */
    IPAR[3] = 5;        /* method selection: Rodas4 */
    
    wait_for_dma(i);
    
    for(j=0; j<VECTOR_LENGTH; j++)
    {
        /* Deinterlace data */
        for(k=0; k<NSPEC; k++)
        {
            conc[1].data[k] = conc[0].data[k*VECTOR_LENGTH+j];
        }
        
        /* Point method at current data */
        C   = &(conc[1].data[0]);
        VAR = &(conc[1].data[0]);
        FIX = &(conc[1].data[NFIXST]);
        
        /* Reset statistics for each integration */
        for(k=0; k<8; k++)
        {
            IPAR[10+k] = stats[k];
        }
        IPAR[12] = 0;
        
        /* Integrate */
        #if DO_CHEMISTRY == 1   /* Need this here for linker */
        IERR = Rosenbrock(VAR, TIME, TIME+DT, ATOL, RTOL, RPAR, IPAR);
        #endif
        
        if(IERR < 0)
        {
            printf("\n Rosenbrock: Unsucessful step at T=%g: IERR=%d\n", TIME, IERR);
            set_status(SPE_STATUS_STOPPED);
            exit(1);
        }
         
        /* Reinterlace data */
        for(k=0; k<NSPEC; k++)
        {
            conc[0].data[k*VECTOR_LENGTH+j] = conc[1].data[k];
        }
    }
    
    /* Record final statistics */
    for(k=0; k<8; k++)
    {
        stats[k] = IPAR[10+k];
    }
    
    /* Record last step for next method invocation */
    STEPMIN = RPAR[11];
}

/**
 * Writes chemistry data to main memory
 */
void write_chem_buffer(uint32_t i)
{
    spu_mfcdma32(conc[i].data, (unsigned int)clist[i].data, clist[i].length*sizeof(dma_list_element_t), i, MFC_PUTL_CMD);
    wait_for_dma(i);
}

/**
 * Performs chemistry along the depth of a row-block of the concentration matrix
 */
void saprc99_chem_block(uint64_t argvp)
{
    /* Iterators */
    uint32_t i, j;
    
    uint32_t rows;
    
    timer_start(&metrics.chem);
    
    allocate_dma_list(&clist[0], NSPEC);
    allocate_dma_list(&clist[1], NSPEC);
    allocate_buffer(&conc[0], NSPEC*VECTOR_LENGTH);
    allocate_buffer(&conc[1], NSPEC*VECTOR_LENGTH);
    
    /* Get arguments from main memory synchronously */
    get_chemistry_argv(argvp, &rows);
    
    for(i=0; i<rows; i++)
    {
        printf("        SPE %d: Chemistry on row %d of %d...\n", envv.speid, i, rows);
        for(j=0; j<(ROW_LENGTH/VECTOR_LENGTH); j++)
        {
            fetch_chem_buffer(0, j*16);
            
            saprc99_buffer(0);
            
            write_chem_buffer(0);
        }
    }
    
    free_dma_list(&clist[0]);
    free_dma_list(&clist[1]);
    free_buffer(&conc[0]);
    free_buffer(&conc[1]);

    set_status(SPE_STATUS_WAITING);
}
