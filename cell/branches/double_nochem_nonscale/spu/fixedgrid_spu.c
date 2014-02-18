/*
 *  fixedgrid_spu.c
 *  
 *  Prototype chemical transport model.
 *  SPE program file.
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <spu_mfcio.h>

#include "fixedgrid_spu.h"
#include "common.h"
#include "params.h"
#include "discretize.h"
#include "timer.h"

/* Metrics */
volatile metrics_t metrics __attribute__((aligned(128)));

/* Program arguments */
volatile spe_argv_t argv __attribute__((aligned(128)));

/* Program environment variables */
volatile spe_env_t env __attribute__((aligned(128)));

/* DMA lists */
spe_dma_list_t clist[2];
spe_dma_list_t wlist[2];
spe_dma_list_t dlist[2];

/* Double-buffered local storage */
spe_buffer_t conc[2];
spe_buffer_t wind[2];
spe_buffer_t diff[2];
spe_buffer_t buff[2];

/* Boundary values */
double cbound[4][4];
double wbound[4][4];
double dbound[4][4];

/* This SPU's status code */
int status;

/**
 * Set status code
 */
inline void set_status(uint32_t s)
{
    status = s;
    spu_write_out_mbox(status);
}

/**
 * Update and get status code
 */
inline uint32_t get_status()
{
    if(spu_stat_in_mbox() > 0)
    {
        status = spu_read_in_mbox();
    }
    return status;
}

/**
 * Waits until DMA requests with given tag-mask complete.
 */
inline void wait_for_dma(uint32_t mask)
{
    mfc_write_tag_mask(1<<mask);
    spu_mfcstat(MFC_TAG_UPDATE_ALL);
}

/**
 * Get arguments from main memory synchronously 
 */
inline void get_spe_argv(uint64_t argvp, double *dt, double *size, 
                         uint32_t *block, double *temp, double *time)
{
    timer_start(&metrics.comm);
    
    mfc_get(&argv, argvp, sizeof(spe_argv_t), 31, 0, 0);
    
    wait_for_dma(31);
    
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
inline void create_shift_boundary(int i, double *c, double *w, double *d, int length)
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
inline void fetch_row_buffer(uint32_t i, uint64_t ea_off) 
{
    mfc_getb(conc[i].data, conc[i].ea_base+ea_off, conc[i].length*sizeof(double), i, 0, 0);
    mfc_getb(wind[i].data, wind[i].ea_base+ea_off, wind[i].length*sizeof(double), i, 0, 0);
    mfc_getb(diff[i].data, diff[i].ea_base+ea_off, diff[i].length*sizeof(double), i, 0, 0);
}

/**
 * Fetches column data from main memory
 */
inline void fetch_col_buffer(uint32_t i, uint64_t ea_off)
{
    const uint32_t colsize = NCOLS*sizeof(double);
    uint32_t j;
    
    /* Build DMA lists */
    timer_start(&metrics.array_copy);
    
    clist[i].length = conc[i].length;
    wlist[i].length = wind[i].length;
    dlist[i].length = diff[i].length;    
    
    for(j=0; j<clist[i].length; ++j)
    {
        clist[i].data[j].ea_low = conc[i].ea_base + ea_off + j*colsize;
        wlist[i].data[j].ea_low = wind[i].ea_base + ea_off + j*colsize;
        dlist[i].data[j].ea_low = diff[i].ea_base + ea_off + j*colsize;
    }
    
    timer_stop(&metrics.array_copy);    
    
    /* Fetch data */
    spu_mfcdma32(conc[i].data, (unsigned int)clist[i].data, clist[i].length*sizeof(dma_list_element_t), i, MFC_GETL_CMD);
    spu_mfcdma32(wind[i].data, (unsigned int)wlist[i].data, wlist[i].length*sizeof(dma_list_element_t), i, MFC_GETL_CMD);
    spu_mfcdma32(diff[i].data, (unsigned int)dlist[i].data, dlist[i].length*sizeof(dma_list_element_t), i, MFC_GETL_CMD);    
}

/**
 * Applies transport and chemistry processes to row data
 */
inline void process_row_buffer(uint32_t i, uint32_t length, double size, double dt)
{
    /* Wait for buffer */
    wait_for_dma(i);
    
    /* Build boundary conditions */
    create_shift_boundary(i, conc[i].data, wind[i].data, diff[i].data, length);
    
    /* Chemistry in buffer 0 */
    /* .... */
    
    /* Transport in buffer 0 */
    discretize(length, conc[i].data, wind[i].data, diff[i].data, 
               cbound[i], wbound[i], dbound[i], size, dt, buff[i].data);
}

/**
 * Applies transport and chemistry processes to column data
 */
void process_col_buffer(uint32_t i, uint32_t length, double size, double dt)
{
    uint32_t j;
    
    /* Temp variables */
    uint32_t t0, t1;
    
    /* Extra column buffers */
    double c2[NROWS];
    double w2[NROWS];
    double d2[NROWS];
    double b2[NROWS];
    
    /* Wait for buffer */
    wait_for_dma(i);
    
    /* Deinterlace buffers (FIXME: unroll) */
    timer_start(&metrics.array_copy);
    for(j=0; j<length; j++)
    {
        t1 = 2*j;
        t0 = t1+1;
        c2[j] = conc[i].data[t0];
        w2[j] = wind[i].data[t0];
        d2[j] = diff[i].data[t0];
        conc[i].data[j] = conc[i].data[t1];
        wind[i].data[j] = wind[i].data[t1];
        diff[i].data[j] = diff[i].data[t1];
    }
    timer_stop(&metrics.array_copy);
    
    /* Prepare boundary values */
    create_shift_boundary(i, conc[i].data, wind[i].data, diff[i].data, length);
    create_shift_boundary(i+1, c2, w2, d2, length);
    
    /* Transport in buffer */
    discretize(length, conc[i].data, wind[i].data, diff[i].data, cbound[i], wbound[i], dbound[i], size, dt, buff[i].data);
    discretize(length, c2, w2, d2, cbound[i+1], wbound[i+1], dbound[i+1], size, dt, b2);
    
    /* Reinterlace buffers (FIXME: unroll) */
    timer_start(&metrics.array_copy);
    for(j=length-1; j>0; j--)
    {
        t1 = 2*j;
        t0 = t1+1;
        buff[i].data[t0] = b2[j];
        buff[i].data[t1] = buff[i].data[j];
    }
    buff[i].data[1] = b2[0];
    timer_stop(&metrics.array_copy);
}

/**
 * Writes row data to main memory
 */
inline void write_row_buffer(uint32_t i, uint64_t ea_off)
{
    mfc_put(buff[i].data, buff[i].ea_base+ea_off, buff[i].length*sizeof(double), i, 0, 0);
}

/**
 * Writes column data to main memory
 */
inline void write_col_buffer(uint32_t i)
{
    spu_mfcdma32(buff[i].data, (unsigned int)clist[i].data, clist[i].length*sizeof(dma_list_element_t), i, MFC_PUTL_CMD);
    wait_for_dma(i);
}

/**
 * FIXME: Function assumes an even number of columns in domain
 */
void discretize_col_block(uint64_t argvp)
{
    /* Iterators */
    uint32_t i, b;
        
    /* Arguments */
    double dt;
    double size;
    uint32_t block;
    double temp;
    double time;
    
    timer_start(&metrics.col_discret);
    
    /* Get arguments from main memory synchronously */
    get_spe_argv(argvp, &dt, &size, &block, &temp, &time);
    
    //printf("concp: 0x%llx, windp: 0x%llx, diffp: 0x%llx, dt: %f, size: %f, length: %d, block: %d\n", concp, windp, diffp, dt, size, length, block);
    
    /* Start buffer 0 transfer */
    fetch_col_buffer(0, 0);
    
    /* Start buffer 1 transfer */
    fetch_col_buffer(1, 2*sizeof(double));
    
    for(i=0; i<block - 2; i+=2)
    {
        b = (i>>1) % 2;
        
        /* Process buffer */
        process_col_buffer(b, NROWS, size, dt);
        
        /* Write buffer back to main memory */
        write_col_buffer(b);
        
        /* Start buffer transfer */
        fetch_col_buffer(b, (i+4)*sizeof(double));
    }
    
    /* Discretize final two columns */
    b = (block - 2) % 2;
    
    /* Process buffer */
    process_col_buffer(b, NROWS, size, dt);
    
    /* Write buffer back to main memory */
    write_col_buffer(b);
    
    timer_stop(&metrics.col_discret);
    
    set_status(SPE_STATUS_WAITING);
}

void discretize_row_block(uint64_t argvp)
{
    const uint32_t rowsize = NCOLS*sizeof(double);

    uint32_t i, b;
    
    /* Arguments */
    double dt;
    double size;
    uint32_t block;
    double temp;
    double time;
    
    timer_start(&metrics.row_discret);
    
    /* Get arguments from main memory synchronously */
    get_spe_argv(argvp, &dt, &size, &block, &temp, &time);
    
    //printf("dt: %f, size: %f, block: %d, temp: %f, time: %f\n", dt, size, block, temp, time);
    //printf("concp: 0x%llx, windp: 0x%llx, diffp: 0x%llx, buffp: 0x%llx, length: %d\n", conc[0].ea_base, wind[0].ea_base, diff[0].ea_base, buff[0].ea_base, conc[0].length);
    
    /* Start buffer 0 transfer */
    fetch_row_buffer(0, 0);
    
    /* Start buffer 1 transfer */
    fetch_row_buffer(1, rowsize);
    
    /* Loop over rows in this block */
    for(i=0; i<block-1; i++)
    {
        b = i % 2;
        
        /* Process buffer */
        process_row_buffer(b, NCOLS, size, dt);
        
        /* Write buffer back to main memory */
        write_row_buffer(b, i*rowsize);
        
        /* Start buffer transfer */
        fetch_row_buffer(b, (i+2)*rowsize);
    }
    
    /* Discretize final row */
    b = (block - 1) % 2;
    
    /* Process buffer */
    process_row_buffer(b, NCOLS, size, dt);
    
    /* Write buffer back to main memory */
    write_row_buffer(b, i*rowsize);
    
    timer_stop(&metrics.row_discret);
    
    /* Signal PPE */
    set_status(SPE_STATUS_WAITING);
}

/**
 * SPE program entry point
 */
int main(uint64_t id, uint64_t argvp, uint64_t envp)
{
    int k, t0;
    char sbuff[20];
    
    timer_start(&metrics.wallclock);
    
    status = SPE_STATUS_WAITING;
    
    /* Get environment information */
    timer_start(&metrics.comm);
    mfc_get(&env, envp, sizeof(spe_env_t), 31, 0, 0);
    wait_for_dma(31);
    timer_stop(&metrics.comm);
    
    /* Initialize metrics */
    sprintf(sbuff, "SPE %d of %d", env.speid, env.nprocs);
    metrics_init(&metrics, sbuff);

    /* Unrolled DMA list initialization */
    timer_start(&metrics.array_copy);
    
    t0 = NROWS;
    k = 0;
    
#define SIZE_ELEMENT \
clist[0].data[k].size.all32 = 16; \
wlist[0].data[k].size.all32 = 16; \
dlist[0].data[k].size.all32 = 16; \
clist[1].data[k].size.all32 = 16; \
wlist[1].data[k].size.all32 = 16; \
dlist[1].data[k].size.all32 = 16; \
++k;

    while( (t0 / 8) > 0 )
    {
        SIZE_ELEMENT
        SIZE_ELEMENT
        SIZE_ELEMENT
        SIZE_ELEMENT
        SIZE_ELEMENT
        SIZE_ELEMENT
        SIZE_ELEMENT
        SIZE_ELEMENT
        t0 -= 8;
    }
    while( (t0 / 4) > 0 )
    {
        SIZE_ELEMENT
        SIZE_ELEMENT
        SIZE_ELEMENT
        SIZE_ELEMENT
        t0 -= 4;
    }
    while( (t0 / 2) > 0 )
    {
        SIZE_ELEMENT
        SIZE_ELEMENT
        t0 -= 2;
    }
    if(t0 == 1)
    {
        SIZE_ELEMENT
    }
    
#undef SIZE_ELEMENT
    
    timer_stop(&metrics.array_copy);
    /* End DMA list initialization */
    
    /* Wait for events */
    while((status = get_status()) != SPE_STATUS_STOPPED)
    {
        switch(status)
        {
            case SPE_STATUS_ROW:
                discretize_row_block(argvp);
                break;
                
            case SPE_STATUS_COL:
                discretize_col_block(argvp);
                break;
                
            case SPE_STATUS_WAITING:
            case SPE_STATUS_STOPPED:
                break;
            
        default:
            fprintf(stderr, "Unknown SPE status: %d.\n", status);
            exit(1);
        }
    }
    
    timer_stop(&metrics.wallclock);
    
    /* Send metrics back to PPU */
    mfc_put(&metrics, env.metptr, sizeof(metrics_t), 31, 0, 0);
    wait_for_dma(31);
    
    return 0;
}
