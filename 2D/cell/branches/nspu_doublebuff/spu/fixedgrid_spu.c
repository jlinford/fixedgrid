#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <spu_mfcio.h>

#include "params.h"
#include "fixedgrid.h"
#include "discretize.h"
#include "timer.h"

volatile fixedgrid_spe_argv_t argv __attribute__((aligned(128)));

/* Local storage */
volatile double conc0[MAX_DMA_DBL] __attribute__((aligned(128)));
volatile double wind0[MAX_DMA_DBL] __attribute__((aligned(128)));
volatile double diff0[MAX_DMA_DBL] __attribute__((aligned(128)));
volatile double buff0[MAX_DMA_DBL] __attribute__((aligned(128)));

/* Double buffering */
volatile double conc1[MAX_DMA_DBL] __attribute__((aligned(128)));
volatile double wind1[MAX_DMA_DBL] __attribute__((aligned(128)));
volatile double diff1[MAX_DMA_DBL] __attribute__((aligned(128)));
volatile double buff1[MAX_DMA_DBL] __attribute__((aligned(128)));

int status;

void set_status(int s)
{
    status = s;
    spu_write_out_mbox(status);
}

int get_status()
{
    if(spu_stat_in_mbox() > 0)
    {
        status = spu_read_in_mbox();
    }
    return status;
}

void begin_discretize(uint64_t argvp)
{
    uint32_t i;
    
    uint64_t concp;
    uint64_t windp;
    uint64_t diffp;
    double dt;
    double size;
    uint32_t length;
    uint32_t block;
    
    /* Get arguments from main memory synchronously */
    mfc_get(&argv, argvp, sizeof(fixedgrid_spe_argv_t), 31, 0, 0);
    mfc_write_tag_mask(1<<31);
    mfc_read_tag_status_all();
    
    concp  = argv.arg[0].u64;
    windp  = argv.arg[1].u64;
    diffp  = argv.arg[2].u64;
    dt     = argv.arg[3].dbl;
    size   = argv.arg[4].dbl;
    length = argv.arg[5].u32[0];
    block  = argv.arg[5].u32[1];
    
    //printf("concp: 0x%llx, windp: 0x%llx, diffp: 0x%llx, dt: %f, size: %f, length: %d, block: %d\n", concp, windp, diffp, dt, size, length, block);
    
    /* Fill buffer 0 */
    mfc_get(conc0, concp, length*sizeof(double), 30, 0, 0);
    mfc_get(wind0, windp, length*sizeof(double), 30, 0, 0);
    mfc_get(diff0, diffp, length*sizeof(double), 30, 0, 0);
    
    for(i=1; i<block; i++)
    {
        /* Fill buffer 1 */
        mfc_getb(conc1, concp+(i*length*sizeof(double)), length*sizeof(double), 31, 0, 0);
        mfc_getb(wind1, windp+(i*length*sizeof(double)), length*sizeof(double), 31, 0, 0);
        mfc_getb(diff1, diffp+(i*length*sizeof(double)), length*sizeof(double), 31, 0, 0);

        /* Wait for buffer 0 */
        mfc_write_tag_mask(1<<30);
        mfc_read_tag_status_all();
        
        /* Process buffer 0 */
        discretize(length, conc0, wind0, diff0, size, dt, buff0);
        
        /* Write buffer 0 results */
        mfc_put(buff0, concp+((i-1)*length*sizeof(double)), length*sizeof(double), 30, 1, 1);
        
        if(++i >= block) break;
        
        /* Fill buffer 0 */
        mfc_getb(conc0, concp+(i*length*sizeof(double)), length*sizeof(double), 30, 0, 0);
        mfc_getb(wind0, windp+(i*length*sizeof(double)), length*sizeof(double), 30, 0, 0);
        mfc_getb(diff0, diffp+(i*length*sizeof(double)), length*sizeof(double), 30, 0, 0);
        
        /* Wait for buffer 1 */
        mfc_write_tag_mask(1<<31);
        mfc_read_tag_status_all();
        
        /* Process buffer 1 */
        discretize(length, conc1, wind1, diff1, size, dt, buff1);
        
        /* Write buffer 1 results */
        mfc_put(buff1, concp+((i-1)*length*sizeof(double)), length*sizeof(double), 31, 1, 1);
    }
    
    /* Process remaining buffer */
    if(i % 2 == 0)
    {
        /* Wait for buffer 1 */
        mfc_write_tag_mask(1<<31);
        mfc_read_tag_status_all();
        
        /* Process buffer 1 */
        discretize(length, conc1, wind1, diff1, size, dt, buff1);
        
        /* Write buffer 1 results */
        mfc_put(buff1, concp+((block-1)*length*sizeof(double)), length*sizeof(double), 31, 1, 1);
    }
    else
    {
        /* Wait for buffer 0 */
        mfc_write_tag_mask(1<<30);
        mfc_read_tag_status_all();
        
        /* Process buffer 0 */
        discretize(length, conc0, wind0, diff0, size, dt, buff0);
        
        /* Write buffer 0 results */
        mfc_put(buff0, concp+((block-1)*length*sizeof(double)), length*sizeof(double), 30, 1, 1);
    }
    
    set_status(SPE_STATUS_WAITING);
}

int main(uint64_t id, uint64_t argvp, uint64_t envp)
{
    status = SPE_STATUS_WAITING;
    
    while((status = get_status()) != SPE_STATUS_STOPPED)
    {
        switch(status)
        {
        case SPE_STATUS_WORKING:
            begin_discretize(argvp);
            break;
            
        case SPE_STATUS_WAITING:
        case SPE_STATUS_STOPPED:
            break;
            
        default:
            fprintf(stderr, "Unknown SPE status: %d.\n", status);
            exit(1);
        }
    }
    
    print_timer_summary("---SPU Timers---");
    
    return 0;
}
