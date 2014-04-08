/*
 *  transport.c
 *  fixedgrid
 *
 *  Created by John Linford on 6/17/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <spu_mfcio.h>
#include <spu_intrinsics.h>

#include "fixedgrid_spu.h"
#include "transport.h"
#include "discretize.h"
#include "timer.h"

/* DMA lists */
spe_dma_list_t clist[2];  /* DMA list for conc data */
spe_dma_list_t wlist[2];  /* DMA list for wind data */
spe_dma_list_t dlist[2];  /* DMA list for diffusion data */

/* Buffered local storage */
/* (Note that transport calculations are actually tripple-buffered via buff) */
spe_buffer_t conc[2];     /* in_x */
spe_buffer_t wind[2];     /* in_x */
spe_buffer_t diff[2];     /* in_x */
spe_buffer_t buff[2];     /* out_x */

/* Boundary values */
real_t cbound[4][4];    /* Concentration boundary */
real_t wbound[4][4];    /* Wind boundary */
real_t dbound[4][4];    /* Diffusion boundary */

/**
 * Get arguments from main memory synchronously 
 */
void get_transport_argv(uint64_t argvp, real_t *dt, real_t *size,
                        uint32_t *block, real_t *temp, real_t *time)
{
    timer_start(&metrics.comm);
    
    mfc_get(&argv, argvp, sizeof(spe_argv_t), 5, 0, 0);
    
    wait_for_dma(5);
    
    conc[0].ea_base = argv.arg[0].u64;
    wind[0].ea_base = argv.arg[1].u64;
    diff[0].ea_base = argv.arg[2].u64;
    buff[0].ea_base = argv.arg[0].u64;
    conc[0].length  = argv.arg[5].u32[0];
    wind[0].length  = argv.arg[5].u32[0];
    diff[0].length  = argv.arg[5].u32[0];
    buff[0].length  = argv.arg[5].u32[0];
    
    conc[1].ea_base = conc[0].ea_base;
    wind[1].ea_base = wind[0].ea_base;
    diff[1].ea_base = diff[0].ea_base;
    buff[1].ea_base = buff[0].ea_base;
    conc[1].length  = conc[0].length;
    wind[1].length  = wind[0].length;
    diff[1].length  = diff[0].length;
    buff[1].length  = buff[0].length;
    
    *dt     = argv.arg[3].dbl;
    *size   = argv.arg[4].dbl;
    *block  = argv.arg[5].u32[1];
    *temp   = argv.arg[6].dbl;
    *time   = argv.arg[7].dbl;
    
    timer_stop(&metrics.comm);
}

/**
 * Initializes boundary value arrays.
 */
void create_shift_boundary(int i, real_t *c, real_t *w, real_t *d, int length)
{
    timer_start(&metrics.array_copy);
    cbound[i][0] = c[length-2];
    cbound[i][1] = c[length-1];
    cbound[i][2] = c[0];
    cbound[i][3] = c[1];
    wbound[i][0] = w[length-2];
    wbound[i][1] = w[length-1];
    wbound[i][2] = w[0];
    wbound[i][3] = w[1];
    dbound[i][0] = d[length-2];
    dbound[i][1] = d[length-1];
    dbound[i][2] = d[0];
    dbound[i][3] = d[1];
    timer_stop(&metrics.array_copy);
}

/**
 * Fetches row data from main memory
 */
void fetch_row_buffer(uint32_t i, uint64_t ea_off) 
{
    mfc_get(conc[i].data, conc[i].ea_base+ea_off, ROW_SIZE, i, 0, 0);
    mfc_get(wind[i].data, wind[i].ea_base+ea_off, ROW_SIZE, i, 0, 0);
    mfc_getb(diff[i].data, diff[i].ea_base+ea_off, ROW_SIZE, i, 0, 0);
}

/**
 * Fetches column data from main memory
 */
void fetch_col_buffer(uint32_t i, uint64_t ea_off)
{
    uint32_t j;
    
    /* Build DMA lists */
    timer_start(&metrics.array_copy);
    
    clist[i].length = conc[i].length;
    wlist[i].length = wind[i].length;
    dlist[i].length = diff[i].length;
    
    //init_dmalist_element_size(16);
    
    for(j=0; j<clist[i].length; ++j)
    {
        clist[i].data[j].ea_low = conc[i].ea_base + ea_off + j*ROW_SIZE;
        wlist[i].data[j].ea_low = wind[i].ea_base + ea_off + j*ROW_SIZE;
        dlist[i].data[j].ea_low = diff[i].ea_base + ea_off + j*ROW_SIZE;
        clist[i].data[j].size.all32 = 16;
        wlist[i].data[j].size.all32 = 16;
        dlist[i].data[j].size.all32 = 16;
    }
    
    timer_stop(&metrics.array_copy);    
    
    /* Fetch data */
    spu_mfcdma32(conc[i].data, (unsigned int)clist[i].data, clist[i].length*sizeof(dma_list_element_t), i, MFC_GETL_CMD);
    spu_mfcdma32(wind[i].data, (unsigned int)wlist[i].data, wlist[i].length*sizeof(dma_list_element_t), i, MFC_GETL_CMD);
    spu_mfcdma32(diff[i].data, (unsigned int)dlist[i].data, dlist[i].length*sizeof(dma_list_element_t), i, MFC_GETL_CMD);    
}

/**
 * Applies transport processes to row data
 */
void transport_row_buffer(uint32_t i, uint32_t length, real_t size, real_t dt)
{
    /* Wait for buffer */
    wait_for_dma(i);
    
    /* Build boundary conditions */
    create_shift_boundary(i, conc[i].data, wind[i].data, diff[i].data, length);
    
    /* Transport in buffer */
    discretize_row(length, 
                   conc[i].data, wind[i].data, diff[i].data, 
                   cbound[i], wbound[i], dbound[i], 
                   size, dt, 
                   buff[i].data);
}

/**
 * Applies transport processes to column data
 */
void transport_col_buffer(uint32_t i, uint32_t length, real_t size, real_t dt)
{
    /* Reinterpret data as vectors */
    vector real_t* vconc = (vector real_t*)conc[i].data;
    vector real_t* vwind = (vector real_t*)wind[i].data;
    vector real_t* vdiff = (vector real_t*)diff[i].data;
    vector real_t* vbuff = (vector real_t*)buff[i].data;
    
    /* Splat scalars to vectors */
    vector real_t vsize = spu_splats(size);
    vector real_t vdt = spu_splats(dt);
    
    /* Boundary values */
    vector real_t vcbound[4];
    vector real_t vwbound[4];
    vector real_t vdbound[4];
    
    /* Wait for buffer */
    wait_for_dma(i);
    
    /* Prepare boundary values */
    timer_start(&metrics.array_copy);
    vcbound[0] = vconc[length-2];
    vcbound[1] = vconc[length-1];
    vcbound[2] = vconc[0];
    vcbound[3] = vconc[1];
    vwbound[0] = vwind[length-2];
    vwbound[1] = vwind[length-1];
    vwbound[2] = vwind[0];
    vwbound[3] = vwind[1];
    vdbound[0] = vdiff[length-2];
    vdbound[1] = vdiff[length-1];
    vdbound[2] = vdiff[0];
    vdbound[3] = vdiff[1];
    timer_stop(&metrics.array_copy);
    
    /* Transport in buffer */
    discretize_col(length, vconc, vwind, vdiff, vcbound, vwbound, vdbound, vsize, vdt, vbuff);
}

/**
 * Writes row data to main memory
 */
void write_row_buffer(uint32_t i, uint64_t ea_off)
{
    mfc_put(buff[i].data, buff[i].ea_base+ea_off, ROW_SIZE, i, 0, 0);
    wait_for_dma(i);
}

/**
 * Writes column data to main memory
 */
void write_col_buffer(uint32_t i)
{
    spu_mfcdma32(buff[i].data, (unsigned int)clist[i].data, clist[i].length*sizeof(dma_list_element_t), i, MFC_PUTL_CMD);
    wait_for_dma(i);
}

/**
 * Performs transport in a column-block of the concentration matrix
 */
void discretize_col_block(uint64_t argvp)
{
    /* Iterators */
    uint32_t i, b;
    
    /* Arguments */
    real_t dt;
    real_t size;
    uint32_t block;
    real_t temp;
    real_t time;
    
    timer_start(&metrics.col_discret);
    
    allocate_dma_list(&clist[0], NROWS);
    allocate_dma_list(&clist[1], NROWS);
    allocate_dma_list(&wlist[0], NROWS);
    allocate_dma_list(&wlist[1], NROWS);
    allocate_dma_list(&dlist[0], NROWS);
    allocate_dma_list(&dlist[1], NROWS);

    allocate_buffer(&conc[0], NROWS*VECTOR_LENGTH);
    allocate_buffer(&conc[1], NROWS*VECTOR_LENGTH);
    allocate_buffer(&wind[0], NROWS*VECTOR_LENGTH);
    allocate_buffer(&wind[1], NROWS*VECTOR_LENGTH);
    allocate_buffer(&diff[0], NROWS*VECTOR_LENGTH);
    allocate_buffer(&diff[1], NROWS*VECTOR_LENGTH);
    allocate_buffer(&buff[0], NROWS*VECTOR_LENGTH);
    allocate_buffer(&buff[1], NROWS*VECTOR_LENGTH);    
    
    
    /* Get arguments from main memory synchronously */
    get_transport_argv(argvp, &dt, &size, &block, &temp, &time);
    
    if(block == 1)
    {
        /* Start in_0 transfer */
        fetch_col_buffer(0, 0);
        
        /* Process in_0 */
        transport_col_buffer(0, NROWS, size, dt);
        
        /* Write out_0 back to main memory */
        write_col_buffer(0);
    }
    else
    {
        /* Start in_0 transfer */
        fetch_col_buffer(0, 0);
        
        /* Start in_1 transfer */
        fetch_col_buffer(1, 16);
        
        /* Process in_0 */
        transport_col_buffer(0, NROWS, size, dt);
        
        for(i=0; i<block-2; i++)
        {
            b = i % 2;
            
            /* Write buffer back to main memory */
            write_col_buffer(b);
            
            /* Start buffer transfer */
            fetch_col_buffer(b, (i+2)*16);
            
            /* Process buffer */
            transport_col_buffer((b+1)%2, NROWS, size, dt);
        }
        
        /* Discretize final column vector */
        b = (block - 2) % 2;
        
        /* Write out_b back to main memory */
        write_col_buffer(b);
        
        /* Process in_(b+1) */
        transport_col_buffer((b+1)%2, NROWS, size, dt);
        
        /* Write out_(b+1) back to main memory */
        write_col_buffer((b+1)%2);
    }
    
    free_dma_list(&clist[0]);
    free_dma_list(&clist[1]);
    free_dma_list(&wlist[0]);
    free_dma_list(&wlist[1]);
    free_dma_list(&dlist[0]);
    free_dma_list(&dlist[1]);
    
    free_buffer(&conc[0]);
    free_buffer(&conc[1]);
    free_buffer(&wind[0]);
    free_buffer(&wind[1]);
    free_buffer(&diff[0]);
    free_buffer(&diff[1]);
    free_buffer(&buff[0]);
    free_buffer(&buff[1]);
    
    timer_stop(&metrics.col_discret);
    
    set_status(SPE_STATUS_WAITING);
}

/**
 * Performs transport in a row-block of the concentration matrix
 */
void discretize_row_block(uint64_t argvp)
{
    uint32_t i, b;
    
    /* Arguments */
    real_t dt;
    real_t size;
    uint32_t block;
    real_t temp;
    real_t time;
    
    timer_start(&metrics.row_discret);
    
    allocate_buffer(&conc[0], ROW_LENGTH);
    allocate_buffer(&conc[1], ROW_LENGTH);
    allocate_buffer(&wind[0], ROW_LENGTH);
    allocate_buffer(&wind[1], ROW_LENGTH);
    allocate_buffer(&diff[0], ROW_LENGTH);
    allocate_buffer(&diff[1], ROW_LENGTH);
    allocate_buffer(&buff[0], ROW_LENGTH);
    allocate_buffer(&buff[1], ROW_LENGTH);
    
    /* Get arguments from main memory synchronously */
    get_transport_argv(argvp, &dt, &size, &block, &temp, &time);
    
    if(block == 1)
    {
        /* Start in_0 transfer */
        fetch_row_buffer(0, 0);
        
        /* Process in_0 */
        transport_row_buffer(0, NCOLS, size, dt);
        
        /* Write out_0 back to main memory */
        write_row_buffer(0, 0);
    }
    else
    {
        /* Start in_0 transfer */
        fetch_row_buffer(0, 0);
        
        /* Start in_1 transfer */
        fetch_row_buffer(1, ROW_SIZE);
        
        /* Process in_0 */
        transport_row_buffer(0, NCOLS, size, dt);
        
        /* Loop over rows in this block */
        for(i=0; i<block-2; i++)
        {
            b = i % 2;
            
            /* Write out_b back to main memory */
            write_row_buffer(b, i*ROW_SIZE);
            
            /* Start in_b transfer */
            fetch_row_buffer(b, (i+2)*ROW_SIZE);
            
            /* Process in_(b+1) */
            transport_row_buffer((b+1)%2, NCOLS, size, dt);
        }
        
        /* Discretize final row */
        b = (block - 2) % 2;
        
        /* Write out_b back to main memory */
        write_row_buffer(b, i*ROW_SIZE);
        
        /* Process in_(b+1) */
        transport_row_buffer((b+1)%2, NCOLS, size, dt);
        
        /* Write out_(b+1) back to main memory */
        write_row_buffer((b+1)%2, (i+1)*ROW_SIZE);
    }
    
    free_buffer(&conc[0]);
    free_buffer(&conc[1]);
    free_buffer(&wind[0]);
    free_buffer(&wind[1]);
    free_buffer(&diff[0]);
    free_buffer(&diff[1]);
    free_buffer(&buff[0]);
    free_buffer(&buff[1]);
    
    timer_stop(&metrics.row_discret);
    
    /* Signal PPE */
    set_status(SPE_STATUS_WAITING);
}

