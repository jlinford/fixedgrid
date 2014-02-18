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
#include "timer.h"
#include "chemistry.h"
#include "transport.h"
#include "saprc99_Global.h"

/* Metrics */
volatile metrics_t metrics __attribute__((aligned(128)));

/* Program arguments */
volatile spe_argv_t argv __attribute__((aligned(128)));

/* Program environment variables */
volatile spe_envv_t envv __attribute__((aligned(128)));

/* This SPU's status code */
int status;

/* KPP-generated SAPRC'99 mechanism data */
volatile double * C;        /* Concentration of all species */
volatile double * VAR;      /* First variable species */
volatile double * FIX;      /* First fixed species */
double RCONST[NREACT];      /* Rate constants (global) */
double ATOL[NVAR];          /* Absolute tolerance */
double RTOL[NVAR];          /* Relative tolerance */
double TIME;                /* Current integration time */
double DT;                  /* Integration step */
double SUN;                 /* Sunlight intensity between [0,1] */
double TEMP;                /* Temperature */
double STEPMIN;             /* Lower bound for integration step */

///**
// * Sets size.all32 = size for all dma lists
// */
//void init_dmalist_element_size(uint32_t size)
//{
//    /* Unrolled DMA list initialization */
//    timer_start(&metrics.array_copy);
//    
//    int t0 = NROWS;
//    int k = 0;
//    
//#define SIZE_ELEMENT \
//clist[0].data[k].size.all32 = size; \
//wlist[0].data[k].size.all32 = size; \
//dlist[0].data[k].size.all32 = size; \
//clist[1].data[k].size.all32 = size; \
//wlist[1].data[k].size.all32 = size; \
//dlist[1].data[k].size.all32 = size; \
//++k;
//    
//    while( (t0 / 8) > 0 )
//    {
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        t0 -= 8;
//    }
//    while( (t0 / 4) > 0 )
//    {
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        t0 -= 4;
//    }
//    while( (t0 / 2) > 0 )
//    {
//        SIZE_ELEMENT
//        SIZE_ELEMENT
//        t0 -= 2;
//    }
//    if(t0 == 1)
//    {
//        SIZE_ELEMENT
//    }
//    
//#undef SIZE_ELEMENT
//    
//    timer_stop(&metrics.array_copy);
//    /* End DMA list initialization */
//}

//void init_dmalist_element_size(uint32_t size)
//{
//    uint32_t k;
//    for(k=0; k<clist[0].length; ++k)
//    {
//        clist[0].data[k].size.all32 = size;
//        wlist[0].data[k].size.all32 = size;
//        dlist[0].data[k].size.all32 = size;
//        clist[1].data[k].size.all32 = size;
//        wlist[1].data[k].size.all32 = size;
//        dlist[1].data[k].size.all32 = size;
//    }
//}

/**
 * Initalizes program environment variables and metrics.
 */
void init_environment(uint64_t envvp)
{
    char sbuff[20];
    
    timer_start(&metrics.comm);
    mfc_get(&envv, envvp, sizeof(spe_envv_t), 31, 0, 0);
    wait_for_dma(31);
    timer_stop(&metrics.comm);
    
    /* Initialize metrics */
    sprintf(sbuff, "SPE %d of %d", envv.speid, envv.nprocs);
    metrics_init(&metrics, sbuff);    
}

/**
 * Initializes model parameters.
 */
void init_model(uint64_t argvp)
{
    uint32_t i;
    
    /* Get model parameters from main memory */
    mfc_get(&argv, argvp, sizeof(spe_argv_t), 31, 0, 0);
    wait_for_dma(31);
    
    timer_start(&metrics.array_init);
    
    /* constant rate coefficients */
    RCONST[149] = 1.5e-11;
    /* END constant rate coefficients */
    
    /* Initialize tolerances */
    for( i = 0; i < NVAR; i++ ) 
    {
        ATOL[i] = 1.0;
        RTOL[i] = 1.0e-3;
    }
    
    /* Set saprc'99 parameters */
    TIME = (real_t)argv.arg[4].dbl;
    DT   = (real_t)argv.arg[6].dbl;
    SUN  = (real_t)0.0;
    TEMP = (real_t)argv.arg[7].dbl;
    STEPMIN = 0.01;
    
    timer_stop(&metrics.array_init);
}


/**
 * SPE program entry point
 */
int main(uint64_t id __attribute__((unused)), uint64_t argvp, uint64_t envvp)
{
    timer_start(&metrics.wallclock);
    
    status = SPE_STATUS_INIT;
    
    init_environment(envvp);
    
    init_model(argvp);
    
    set_status(SPE_STATUS_WAITING);
    
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
                
            case SPE_STATUS_CHEM:
                saprc99_chem_block(argvp);
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
    mfc_put(&metrics, envv.metptr, sizeof(metrics_t), 31, 0, 0);
    wait_for_dma(31);
    
    return 0;
}
