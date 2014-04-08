#ifndef __FIXEDGRID_H__
#define __FIXEDGRID_H__

#include <stdint.h>

#define MAX_THREADS 8
#define MAX_DMA_DBL 2048

#define SPE_STATUS_WAITING  0
#define SPE_STATUS_WORKING  1
#define SPE_STATUS_STOPPED  100

typedef short bool;
#define TRUE  1
#define FALSE 0

typedef union
{
    double   dbl;
    uint64_t u64;
    uint32_t u32[2];
} arg_t;

typedef struct
{
    arg_t arg[16];
} fixedgrid_spe_argv_t __attribute__((aligned(128)));


#endif
