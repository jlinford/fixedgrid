#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <spu_mfcio.h>

#include "params.h"
#include "fixedgrid.h"
#include "discretize.h"

volatile fixedgrid_spe_argv_t argv __attribute__((aligned(128)));

volatile double conc[MAX_DMA_DBL] __attribute__((aligned(128)));
volatile double wind[MAX_DMA_DBL] __attribute__((aligned(128)));
volatile double diff[MAX_DMA_DBL] __attribute__((aligned(128)));
volatile double buff[MAX_DMA_DBL] __attribute__((aligned(128)));

int status;

void begin_discretize(uint64_t argvp)
{
    int i;
    uint64_t concp;
    uint64_t windp;
    uint64_t diffp;
    double dt;
    uint32_t length;
    double size;

    /* Get arguments from main memory synchronously */
    mfc_get(&argv, argvp, sizeof(fixedgrid_spe_argv_t), 31, 0, 0);
    mfc_write_tag_mask(1<<31);
    mfc_read_tag_status_all();
    
    concp  = argv.arg[0].u64;
    windp  = argv.arg[1].u64;
    diffp  = argv.arg[2].u64;
    dt     = argv.arg[3].dbl;
    length = argv.arg[4].u32[0];
    size   = argv.arg[5].dbl;
    
    //printf("concp: 0x%llx, windp: 0x%llx, diffp: 0x%llx, dt: %f, length: %d, size: %f\n", concp, windp, diffp, dt, length, size); 
    
    /* Read in data from main memory */
    mfc_get(conc, concp, length*sizeof(double), 31, 0, 0);
    mfc_get(wind, windp, length*sizeof(double), 31, 0, 0);
    mfc_get(diff, diffp, length*sizeof(double), 31, 0, 0);
    mfc_read_tag_status_all();

    discretize(length, conc, wind, diff, size, dt, buff);
    
    /* Write data to memory */
    mfc_put(buff, concp, length*sizeof(double), 30, 0, 0);
    mfc_write_tag_mask(1<<30);
    mfc_read_tag_status_all();
}

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

int main(uint64_t id, uint64_t argvp, uint64_t envp)
{
    status = SPE_STATUS_WAITING;
    
    while((status = get_status()) != SPE_STATUS_STOPPED)
    {
        switch(status)
        {
        case SPE_STATUS_WORKING:
            begin_discretize(argvp);
            set_status(SPE_STATUS_WAITING);
            break;
            
        case SPE_STATUS_WAITING:
        case SPE_STATUS_STOPPED:
            break;
            
        default:
            fprintf(stderr, "Unknown SPE status: %d.\n", status);
            exit(1);
        }
    }
    
    return 0;
}
