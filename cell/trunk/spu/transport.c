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

/* Shuffle buffers for row discretization */
spe_buffer_t shuffle[4];

/* DMA lists */
spe_dma_list_t clist[3];  /* DMA list for conc data */
spe_dma_list_t wlist[3];  /* DMA list for wind data */
spe_dma_list_t dlist[3];  /* DMA list for diffusion data */

/* Buffered local storage */
spe_buffer_t conc[3];     /* in_x */
spe_buffer_t wind[3];     /* in_x */
spe_buffer_t diff[3];     /* in_x */
spe_buffer_t buff[3];     /* out_x */

/* Boundary values */
vector real_t cbound[3][4];    /* Concentration boundary */
vector real_t wbound[3][4];    /* Wind boundary */
vector real_t dbound[3][4];    /* Diffusion boundary */

inline void copy_as_vector(real_t* s, vector real_t* v, uint32_t n)
{
    uint32_t i = 0;
    
    //timer_start(&metrics.array_copy);
    
    while(n > 8)
    {
        v[i] = spu_splats(s[i]); ++i;
        v[i] = spu_splats(s[i]); ++i;
        v[i] = spu_splats(s[i]); ++i;
        v[i] = spu_splats(s[i]); ++i;
        v[i] = spu_splats(s[i]); ++i;
        v[i] = spu_splats(s[i]); ++i;
        v[i] = spu_splats(s[i]); ++i;
        v[i] = spu_splats(s[i]); ++i;
        n -= 8;
    }
    while(n > 4)
    {
        v[i] = spu_splats(s[i]); ++i;
        v[i] = spu_splats(s[i]); ++i;
        v[i] = spu_splats(s[i]); ++i;
        v[i] = spu_splats(s[i]); ++i;
        n -= 4;
    }
    while(n > 0)
    {
        v[i] = spu_splats(s[i]); ++i;
        --n;
    }
    
    //timer_stop(&metrics.array_copy);
}

inline void copy_as_scalar(vector real_t* v, real_t* s, uint32_t n)
{
    uint32_t i = 0;
    
    //timer_start(&metrics.array_copy);
    
    while(n > 8)
    {
        s[i] = spu_extract(v[i], 0); ++i;
        s[i] = spu_extract(v[i], 0); ++i;
        s[i] = spu_extract(v[i], 0); ++i;
        s[i] = spu_extract(v[i], 0); ++i;
        s[i] = spu_extract(v[i], 0); ++i;
        s[i] = spu_extract(v[i], 0); ++i;
        s[i] = spu_extract(v[i], 0); ++i;
        s[i] = spu_extract(v[i], 0); ++i;        
        n -= 8;
    }
    while(n > 4)
    {
        s[i] = spu_extract(v[i], 0); ++i;
        s[i] = spu_extract(v[i], 0); ++i;
        s[i] = spu_extract(v[i], 0); ++i;
        s[i] = spu_extract(v[i], 0); ++i;
        n -= 4;
    }
    while(n > 0)
    {
        s[i] = spu_extract(v[i], 0); ++i;
        --n;
    }
    
    //timer_stop(&metrics.array_copy);
}

/**
 * Get arguments from main memory synchronously 
 */
void get_transport_argv(uint64_t argvp, real_t *dt, real_t *size, uint32_t *block)
{
    mfc_get(&argv, argvp, sizeof(spe_argv_t), GET_ARG_TAG_MASK, 0, 0);
    wait_for_dma(GET_ARG_TAG_MASK);
    
    conc[0].ea_base = argv.arg[0].u64;
    wind[0].ea_base = argv.arg[1].u64;
    diff[0].ea_base = argv.arg[2].u64;
    buff[0].ea_base = argv.arg[0].u64;
    conc[0].length  = argv.arg[5].u32[0];
    wind[0].length  = conc[0].length;
    diff[0].length  = conc[0].length;
    buff[0].length  = conc[0].length;
    
    conc[1].ea_base = conc[0].ea_base;
    wind[1].ea_base = wind[0].ea_base;
    diff[1].ea_base = diff[0].ea_base;
    buff[1].ea_base = buff[0].ea_base;
    conc[1].length  = conc[0].length;
    wind[1].length  = wind[0].length;
    diff[1].length  = diff[0].length;
    buff[1].length  = buff[0].length;
    
    conc[2].ea_base = conc[0].ea_base;
    wind[2].ea_base = wind[0].ea_base;
    diff[2].ea_base = diff[0].ea_base;
    buff[2].ea_base = buff[0].ea_base;
    conc[2].length  = conc[0].length;
    wind[2].length  = wind[0].length;
    diff[2].length  = diff[0].length;
    buff[2].length  = buff[0].length;
    
    clist[0].length = conc[0].length;
    wlist[0].length = wind[0].length;
    dlist[0].length = diff[0].length;    
    
    clist[1].length = conc[1].length;
    wlist[1].length = wind[1].length;
    dlist[1].length = diff[1].length;
    
    clist[2].length = conc[1].length;
    wlist[2].length = wind[1].length;
    dlist[2].length = diff[1].length;    
    
    shuffle[0].length = conc[0].length;
    shuffle[1].length = wind[0].length;
    shuffle[2].length = diff[0].length;
    shuffle[3].length = buff[0].length;
    
    *dt     = argv.arg[3].dbl;
    *size   = argv.arg[4].dbl;
    *block  = argv.arg[5].u32[1];
}

/**
 * Initializes boundary value arrays.
 */
void create_shift_boundary(uint32_t i, vector real_t *c, vector real_t *w, vector real_t *d, uint32_t length)
{
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
}

/**
 * Fetches row data from main memory
 */
void fetch_row_buffer(uint32_t i, uint64_t ea_off) 
{
    const uint32_t rtag1 = i+9;
    const uint32_t rtag2 = i+10;
    const uint32_t rtag3 = i+11;
        
    mfc_get(shuffle[0].data, conc[i].ea_base+ea_off, 
            conc[i].length*sizeof(real_t), rtag1, 0, 0);
    mfc_get(shuffle[1].data, wind[i].ea_base+ea_off, 
            wind[i].length*sizeof(real_t), rtag2, 0, 0);
    mfc_get(shuffle[2].data, diff[i].ea_base+ea_off, 
            diff[i].length*sizeof(real_t), rtag3, 0, 0);
    
    wait_for_dma(rtag1);
    copy_as_vector(shuffle[0].data, (vector real_t*)conc[i].data, conc[i].length);
    wait_for_dma(rtag2);
    copy_as_vector(shuffle[1].data, (vector real_t*)wind[i].data, wind[i].length);
    wait_for_dma(rtag3);
    copy_as_vector(shuffle[2].data, (vector real_t*)diff[i].data, diff[i].length);
}

/**
 * Fetches column data from main memory
 */
void fetch_col_buffer(uint32_t i, uint64_t ea_off)
{
    uint32_t j, x;
    
    /* Build DMA lists */
    //timer_start(&metrics.array_copy);
    
    #define UNROLL_ELEMENT \
    clist[i].data[j].eal = conc[i].ea_base + ea_off + j*ROW_SIZE; \
    wlist[i].data[j].eal = wind[i].ea_base + ea_off + j*ROW_SIZE; \
    dlist[i].data[j].eal = diff[i].ea_base + ea_off + j*ROW_SIZE; \
    clist[i].data[j].size = 16; \
    wlist[i].data[j].size = 16; \
    dlist[i].data[j].size = 16
    
    j=0; x = clist[i].length;
    while(x > 8)
    {
        UNROLL_ELEMENT; ++j;
        UNROLL_ELEMENT; ++j;
        UNROLL_ELEMENT; ++j;
        UNROLL_ELEMENT; ++j;
        UNROLL_ELEMENT; ++j;
        UNROLL_ELEMENT; ++j;
        UNROLL_ELEMENT; ++j;
        UNROLL_ELEMENT; ++j;       
        x -= 8;
    }
    while(x > 4)
    {
        UNROLL_ELEMENT; ++j;
        UNROLL_ELEMENT; ++j;
        UNROLL_ELEMENT; ++j;
        UNROLL_ELEMENT; ++j;
        x -= 4;
    }
    while(x > 0)
    {
        UNROLL_ELEMENT; ++j;
        --x;
    }
    
    #undef UNROLL_ELEMENT
    
    //timer_stop(&metrics.array_copy);    
    
    /* Fetch data */
    mfc_getlb(conc[i].data, conc[i].ea_base, clist[i].data, 
              clist[i].length*sizeof(mfc_list_element_t), i, 0, 0);
    mfc_getl(wind[i].data, wind[i].ea_base, wlist[i].data, 
             wlist[i].length*sizeof(mfc_list_element_t), i, 0, 0);
    mfc_getl(diff[i].data, diff[i].ea_base, dlist[i].data, 
             dlist[i].length*sizeof(mfc_list_element_t), i, 0, 0);
}

/**
 * Applies transport processes to row data
 */
void transport_buffer(uint32_t i, uint32_t length, real_t size, real_t dt)
{
    /* Reinterpret data as vectors */
    vector real_t* vconc = (vector real_t*)conc[i].data;
    vector real_t* vwind = (vector real_t*)wind[i].data;
    vector real_t* vdiff = (vector real_t*)diff[i].data;
    vector real_t* vbuff = (vector real_t*)buff[i].data;
    
    /* Splat scalars to vectors */
    vector real_t vsize = spu_splats(size);
    vector real_t vdt = spu_splats(dt);
        
    /* Wait for input buffer */
    wait_for_dma(i);
    
    /* Prepare boundary values */
    create_shift_boundary(i, vconc, vwind, vdiff, conc[i].length);
    
    /* Transport in buffer */
    discretize(length, vconc, vwind, vdiff, 
               cbound[i], wbound[i], dbound[i], vsize, vdt, vbuff);
}

/**
 * Writes row data to main memory
 */
void write_row_buffer(uint32_t i, uint64_t ea_off)
{
    copy_as_scalar((vector real_t*)buff[i].data, shuffle[3].data, buff[i].length);
    
    mfc_put(shuffle[3].data, buff[i].ea_base+ea_off, 
            buff[i].length*sizeof(real_t), i, 0, 0);
}

/**
 * Writes column data to main memory
 */
void write_col_buffer(uint32_t i)
{
    mfc_putl(buff[i].data, conc[i].ea_base, clist[i].data, 
             clist[i].length*sizeof(mfc_list_element_t), i, 0, 0);
}

/**
 * Performs transport in a column-block of the concentration matrix
 */
void discretize_col_block(uint64_t argvp)
{
    /* Iterators */
    uint32_t i, f, p, w;
    
    /* Arguments */
    real_t dt;
    real_t size;
    uint32_t block;
    
    timer_start(&metrics.col_discret);
    
    /* Get arguments from main memory synchronously */
    get_transport_argv(argvp, &dt, &size, &block);
    
    allocate_dma_list(&clist[0]);
    allocate_dma_list(&clist[1]);
    allocate_dma_list(&clist[2]);
    allocate_dma_list(&wlist[0]);
    allocate_dma_list(&wlist[1]);
    allocate_dma_list(&wlist[2]);
    allocate_dma_list(&dlist[0]);
    allocate_dma_list(&dlist[1]);
    allocate_dma_list(&dlist[2]);

    allocate_buffer(&conc[0]);
    allocate_buffer(&conc[1]);
    allocate_buffer(&conc[2]);
    allocate_buffer(&wind[0]);
    allocate_buffer(&wind[1]);
    allocate_buffer(&wind[2]);
    allocate_buffer(&diff[0]);
    allocate_buffer(&diff[1]);
    allocate_buffer(&diff[2]);
    allocate_buffer(&buff[0]);
    allocate_buffer(&buff[1]);
    allocate_buffer(&buff[2]);
    
    if(block == 1)
    {
        /* Start in_0 transfer */
        fetch_col_buffer(0, 0);
        
        /* Process in_0 */
        transport_buffer(0, conc[0].length, size, dt);
        
        /* Write out_0 back to main memory */
        write_col_buffer(0);
        
        wait_for_dma(0);
    }
    else
    {
        /* Start in_0 transfer */
        fetch_col_buffer(0, 0);
        
        /* Start in_1 transfer */
        fetch_col_buffer(1, 16);
        
        /* Process in_0 */
        transport_buffer(0, conc[0].length, size, dt);
        
        for(i=0; i<block-2; i++)
        {
            w = i % 3;
            p = (i+1) % 3;
            f = (i+2) % 3;
            
            /* Write buffer back to main memory */
            write_col_buffer(w);
            
            /* Start buffer transfer */
            fetch_col_buffer(f, (i+2)*16);
            
            /* Process buffer */
            transport_buffer(p, conc[p].length, size, dt);
        }
        
        /* Discretize final column vector */
        w = i % 3;
        p = (i+1) % 3;
        
        /* Write out_b back to main memory */
        write_col_buffer(w);
        
        /* Process in_(b+1) */
        transport_buffer(p, conc[p].length, size, dt);
        
        /* Write out_(b+1) back to main memory */
        write_col_buffer(p);
        
        /* Make sure DMA is complete before we exit */
        mfc_write_tag_mask( (1<<w) | (1<<p) );
        spu_mfcstat(MFC_TAG_UPDATE_ALL);
    }
    
    free_dma_list(&clist[0]);
    free_dma_list(&clist[1]);
    free_dma_list(&clist[2]);
    free_dma_list(&wlist[0]);
    free_dma_list(&wlist[1]);
    free_dma_list(&wlist[2]);
    free_dma_list(&dlist[0]);
    free_dma_list(&dlist[1]);
    free_dma_list(&dlist[2]);
    
    free_buffer(&conc[0]);
    free_buffer(&conc[1]);
    free_buffer(&conc[2]);
    free_buffer(&wind[0]);
    free_buffer(&wind[1]);
    free_buffer(&wind[2]);
    free_buffer(&diff[0]);
    free_buffer(&diff[1]);
    free_buffer(&diff[2]);
    free_buffer(&buff[0]);
    free_buffer(&buff[1]);
    free_buffer(&buff[2]);
    
    timer_stop(&metrics.col_discret);
    
    set_status(SPE_STATUS_WAITING);
}

/**
 * Performs transport in a row-block of the concentration matrix
 */
void discretize_row_block(uint64_t argvp)
{
    uint32_t i, f, p, w;
    
    /* Arguments */
    real_t dt;
    real_t size;
    uint32_t block;
    
    timer_start(&metrics.row_discret);
        
    /* Get arguments from main memory synchronously */
    get_transport_argv(argvp, &dt, &size, &block);
    
    allocate_buffer(&conc[0]);
    allocate_buffer(&conc[1]);
    allocate_buffer(&conc[2]);
    allocate_buffer(&wind[0]);
    allocate_buffer(&wind[1]);
    allocate_buffer(&wind[2]);
    allocate_buffer(&diff[0]);
    allocate_buffer(&diff[1]);
    allocate_buffer(&diff[2]);
    allocate_buffer(&buff[0]);
    allocate_buffer(&buff[1]);
    allocate_buffer(&buff[2]);
    allocate_buffer(&shuffle[0]);
    allocate_buffer(&shuffle[1]);
    allocate_buffer(&shuffle[2]);
    allocate_buffer(&shuffle[3]);
    
    if(block == 1)
    {
        /* Start in_0 transfer */
        fetch_row_buffer(0, 0);
        
        /* Process in_0 */
        transport_buffer(0, conc[0].length, size, dt);
        
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
        transport_buffer(0, conc[0].length, size, dt);
        
        /* Loop over rows in this block */
        for(i=0; i<block-2; i++)
        {
            w = i % 3;
            p = (i+1) % 3;
            f = (i+2) % 3;
            
            /* Write out_b back to main memory */
            write_row_buffer(w, i*ROW_SIZE);
            
            /* Start in_b transfer */
            fetch_row_buffer(f, (i+2)*ROW_SIZE);
            
            /* Process in_(b+1) */
            transport_buffer(p, conc[p].length, size, dt);
        }
        
        /* Discretize final row */
        w = i % 3;
        p = (i+1) % 3;
        
        /* Write out_b back to main memory */
        write_row_buffer(w, i*ROW_SIZE);
        
        /* Process in_(b+1) */
        transport_buffer(p, conc[p].length, size, dt);
        
        /* Write out_(b+1) back to main memory */
        write_row_buffer(p, (i+1)*ROW_SIZE);

        /* Make sure DMA is complete before we exit */
        mfc_write_tag_mask( (1<<w) | (1<<p) );
        spu_mfcstat(MFC_TAG_UPDATE_ALL);
    }
    
    free_buffer(&conc[0]);
    free_buffer(&conc[1]);
    free_buffer(&conc[2]);
    free_buffer(&wind[0]);
    free_buffer(&wind[1]);
    free_buffer(&wind[2]);
    free_buffer(&diff[0]);
    free_buffer(&diff[1]);
    free_buffer(&diff[2]);
    free_buffer(&buff[0]);
    free_buffer(&buff[1]);
    free_buffer(&buff[2]);
    free_buffer(&shuffle[0]);
    free_buffer(&shuffle[1]);
    free_buffer(&shuffle[2]);
    free_buffer(&shuffle[3]);

    timer_stop(&metrics.row_discret);
    
    /* Signal PPE */
    set_status(SPE_STATUS_WAITING);
}

