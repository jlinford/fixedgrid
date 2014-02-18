/*
 *  common.h
 *  
 *  Definitions which must always agree on both PPU and SPU.
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __COMMON_H__
#define __COMMON_H__

/**************************************************
 * Includes                                       *
 **************************************************/

#include <stdint.h>
#include "params.h"

/**************************************************
 * Macros                                         *
 **************************************************/

#define TRUE  1
#define FALSE 0

#define SPE_STATUS_WAITING  0
#define SPE_STATUS_INIT     1
#define SPE_STATUS_ROW      2
#define SPE_STATUS_COL      3
#define SPE_STATUS_CHEM     4
#define SPE_STATUS_STOPPED  100

/**************************************************
 * Data types                                     *
 **************************************************/

/* Boolean type */
typedef short bool;

/* SPE program arguments */
typedef struct spe_argv
{
    union
    {
        double   dbl;
        uint64_t u64;
        float    flt[2];
        uint32_t u32[2];
    } 
    arg[16];
} spe_argv_t;

/* SPE program environment variables */
typedef struct spe_envv
{
    uint32_t speid;
    uint32_t nprocs;
    uint32_t ls_status;
    uint64_t ea_status;
    uint64_t ea_metrics;
} spe_envv_t;

typedef union spe_status
{
    uint32_t value;
    uint32_t all[4];
} spe_status_t;

#endif
