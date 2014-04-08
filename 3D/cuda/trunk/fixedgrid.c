/*
 *  fixedgrid.c
 *  
 *  Prototype chemical transport model.
 *  Main program file.
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

#include "fileio.h"
#include "fixedgrid.h"
#include "transport.h"
#include "saprc99_Monitor.h"
#include "cuda_host.h"

/* Program data */
fixedgrid_t G_GLOBAL;

/**
 * Applies KPP-generated SAPRC99 chemical mechanism
 * to all concentration data.
 */
void start_saprc99(fixedgrid_t* G)
{
#if DO_CHEMISTRY == 1
    timer_start(&G->ppe_metrics.chem);
    timer_stop(&G->ppe_metrics.chem);
#endif
}

/**
 * Initializes the model
 */
void init_model(fixedgrid_t* G)
{
    uint32_t x, y, z, s;
    
    /* Chemistry buffer */
    real_t chemBuff[NSPEC];
    
    /* Initialize metrics */
    metrics_init(&G->metrics, "PPE");
        
    /* Initialize time frame */
    /* FIXME: year is ignored */
    G->tstart = day2sec(START_DOY) + hour2sec(START_HOUR) + minute2sec(START_MIN);
    G->tend   = day2sec(END_DOY)   + hour2sec(END_HOUR)   + minute2sec(END_MIN);
    G->dt = STEP_SIZE;
    G->time = G->tstart;
        
    /* Initialize chemistry and concentration data */
    printf("Loading chemistry and concentration data... ");
    timer_start(&G->metrics.array_init);

    #if DO_CHEMISTRY == 1
    
    saprc99_Initialize(chemBuff);
    
    for(s=0; s<NSPEC; s++)
    {
        for(z=0; z<NZ; z++)
        {
            for(y=0; y<NY; y++)
            {
                for(x=0; x<NX; x++)
                {
                    G->conc(x, y, z, s) = chemBuff[s];
                }
            }
        }
    }
    
    #else
    
    for(z=0; z<NZ; z++)
    {
        for(y=0; y<NY; y++)
        {
            for(x=0; x<NX; x++)
            {
                G->conc(x, y, z, ind_O3) = O3_INIT;
            }
        }
    }
    
    #endif
    
    timer_stop(&G->metrics.array_init);
    printf("done.\n");
    
    /* Initialize wind field */
    printf("Loading wind field data... ");
    timer_start(&G->metrics.array_init);
    
    for(z=0; z<NZ; z++)
    {
        for(y=0; y<NY; y++)
        {
            for(x=0; x<NX; x++)
            {
                G->wind_u(x, y, z) = WIND_U_INIT;
                G->wind_v(x, y, z) = WIND_V_INIT;
                G->wind_w(x, y, z) = WIND_W_INIT;
            }
        }
    }
    
    timer_stop(&G->metrics.array_init);
    printf("done.\n");
    
    /* Initialize diffusion field */
    printf("Loading diffusion field data... ");
    timer_start(&G->metrics.array_init);
    
    for(z=0; z<NZ; z++)
    {
        for(y=0; y<NY; y++)
        {
            for(x=0; x<NX; x++)
            {
                G->diff_h(x, y, z) = DIFF_H_INIT;
                G->diff_v(x, y, z) = DIFF_V_INIT;
            }
        }
    }
    
    timer_stop(&G->metrics.array_init);
    printf("done.\n");
    
    /* Initialize temperature field */
    printf("Loading temperature field data... ");
    timer_start(&G->metrics.array_init);
    
    for(z=0; z<NZ; z++)
    {
        for(y=0; y<NY; y++)
        {
            for(x=0; x<NX; x++)
            {
                G->temp(x, y, z) = TEMP_INIT;
            }
        }
    }
    
    timer_stop(&G->metrics.array_init);
    printf("done.\n");
    
    /* Initialize CUDA kernel and device memory */
    printf("Initializing CUDA driver interface... ");
    //FIXME: Start a timer here
    
    CU_SAFE_CALL(init_cuda_driver("data/discretize.cubin"));
    CU_SAFE_CALL(cuMemAlloc(&G->dev_conc, NZ*NY*NX*sizeof(real_t)));
    CU_SAFE_CALL(cuMemAlloc(&G->dev_wind, NZ*NY*NX*sizeof(real_t)));
    CU_SAFE_CALL(cuMemAlloc(&G->dev_diff, NZ*NY*NX*sizeof(real_t)));
    CU_SAFE_CALL(cuMemAlloc(&G->dev_buff, NZ*NY*NX*sizeof(real_t)));
    CU_SAFE_CALL(cuMemAlloc(&G->dev_conc_out, NZ*NY*NX*sizeof(real_t)));
    
    init_discretization_kernel(G);
        
    //FIXME: Stop a timer here
    printf("done.\n");
}

void update_model(fixedgrid_t* G)
{
    /*
     * Could update wind field here...
     */
    
    /*
     * Could update diffusion tensor here...
     */
    
    /*
     * Could update emissions here...
     */
}

/**
 * Processes emission sources
 */
void process_emissions(fixedgrid_t* G)
{
    /* Add O3 plume */
    G->conc(SOURCE_X, SOURCE_Y, SOURCE_Z, ind_O3) += SOURCE_RATE / (DX * DY * DZ);
}

/**
 * Shows all emission sources.
 */
void print_emission_sources()
{
    printf("\n");
    printf("Emission sources (SPEC, X, Y, Z, RATE):\n");
    printf("    (%s, %f, %f, %f, %E)\n", SPC_NAMES[ind_O3], (DX*SOURCE_X + DX*0.5), (DY*SOURCE_Y + DY*0.5), 0.0, SOURCE_RATE);
    printf("\n");
}

/**
 * Prints program banner.
 */
void print_start_banner(fixedgrid_t* G)
{
    int i;
    uint32_t steps;
    
    steps = (G->tend - G->tstart) / G->dt;
    
    printf("\n");
    printf("CONFIGURATION:\n");
    printf("    X DISCRETIZATION:   %s\n", DO_X_DISCRET == TRUE ? "TRUE" : "FALSE");
    printf("    Y DISCRETIZATION:   %s\n", DO_Y_DISCRET == TRUE ? "TRUE" : "FALSE");
    printf("    Z DISCRETIZATION:   %s\n", DO_Z_DISCRET == TRUE ? "TRUE" : "FALSE");
    printf("    SAPRC99 CHEMISTRY:  %s\n", DO_CHEMISTRY == TRUE ? "TRUE" : "FALSE");
    printf("    DOUBLE PRECISION:   %s\n", DOUBLE_PRECISION == TRUE ? "TRUE" : "FALSE");
    printf("\n");
    printf("SPACE DOMAIN:\n");
    printf("    LENGTH (X): %f meters\n", NX*DX);
    printf("    WIDTH  (Y): %f meters\n", NY*DY);
    printf("    DEPTH  (Z): %f meters\n", NZ*DZ);
    printf("\n");
    printf("TIME DOMAIN:\n");
    printf("    FROM  %d:%d.00 on day %d of year %d\n", START_HOUR, START_MIN, START_DOY, START_YEAR);
    printf("    TO    %d:%d.00 on day %d of year %d\n", END_HOUR, END_MIN, END_DOY, END_YEAR);
    printf("    TOTAL %d seconds (%d timesteps of %d seconds)\n", (int)(G->tend-G->tstart), steps, (int)G->dt);
    printf("\n");
    printf("CHEMICAL SPECIES:\n");
    printf("    TOTAL:    %d\n", NSPEC);
    printf("    EXAMINED: ");
    for(i=0; i<NLOOKAT-1; ++i)
    {
        printf("%s, ", SPC_NAMES[ LOOKAT[i] ]);
    }
    printf("%s\n", SPC_NAMES[ LOOKAT[NLOOKAT-1] ]);
    printf("    MONITORED: ");
    for(i=0; i<NMONITOR-1; ++i)
    {
        printf("%s, ", SPC_NAMES[ MONITOR[i] ]);
    }
    printf("%s\n", SPC_NAMES[ MONITOR[NMONITOR-1] ]);
    
    print_emission_sources();
    
    printf("\n");
}

/**
 * PPU program entry point.
 */
int main(int argc, char** argv)
{
    /* Get global memory pointer */
    fixedgrid_t* const G = &G_GLOBAL;
    
    /* Iterators */
    uint32_t k, iter;
    
    /* Start wall clock timer */
    timer_start(&G->metrics.wallclock);

    /* Check dimensions */
    if(NX % BLOCK_X != 0)
    {
        fprintf(stderr, "NX must be a multiple of %d\n", BLOCK_X);
        exit(1);
    }    
    if(NY % BLOCK_Y != 0)
    {
        fprintf(stderr, "NY must be a multiple of %d\n", BLOCK_Y);
        exit(1);
    }
    if(NZ % BLOCK_Z != 0)
    {
        fprintf(stderr, "NZ must be a multiple of %d\n", BLOCK_Z);
        exit(1);
    }
    
    /* Initialize the model parameters */
    init_model(G);
    
    /* Add emissions */
    process_emissions(G);
    
    /* Print startup banner */
    print_start_banner(G);
    
    /* Store initial concentration */
    printf("Writing initial concentration data... ");
    write_conc(G, 0, 0);
    printf("done.\n");    
        
    printf("\n!!!!FIXME: Report # FPEs\n");
        
    /* BEGIN CALCULATIONS */
    for(iter=1, G->time = G->tstart; G->time < G->tend; G->time += G->dt, ++iter)
    {
        start_saprc99(G);
        
        for(k=0; k<NLOOKAT; k++)
        {
            // Copy concentration data to device
            CU_SAFE_CALL(cuMemcpyHtoD(G->dev_conc, &G->conc(0, 0, 0, MONITOR[k]), NX*NY*NZ*sizeof(real_t)));
            
            discretize_all_x(G, G->dt*0.5);
            
            discretize_all_y(G, G->dt*0.5);
            
            discretize_all_z(G, G->dt);
            
            discretize_all_y(G, G->dt*0.5);
            
            discretize_all_x(G, G->dt*0.5);
            
            // Copy updated concentrations back to host
            CU_SAFE_CALL(cuMemcpyDtoH((void*)&G->conc(0, 0, 0, MONITOR[k]), G->dev_conc_out, NX*NY*NZ*sizeof(real_t)));            
        }

        update_model(G);
        
        #if WRITE_EACH_ITER == 1
        write_conc(G, iter, 0);
        #endif

        printf("  After iteration %02d: Model time = %07.2f sec.\n", iter, iter*G->dt);
    }
    /* END CALCULATIONS */
    
    /* Store concentration */
    #if WRITE_EACH_ITER != 1
    write_conc(G, iter-1, 0);
    #endif
    
    /* Show final time */
    printf("\nFinal time: %f seconds.\n", (iter-1)*G->dt);
    
    timer_stop(&G->metrics.wallclock);
    
    /* Write metrics to CSV file */
    write_metrics_as_csv(G, "NVidia CUDA");
    
    /* Cleanup and exit */

    CU_SAFE_CALL(cuMemFree(G->dev_conc));
    CU_SAFE_CALL(cuMemFree(G->dev_wind));
    CU_SAFE_CALL(cuMemFree(G->dev_diff));
    CU_SAFE_CALL(cuMemFree(G->dev_buff));
    CU_SAFE_CALL(cuMemFree(G->dev_conc_out));
    CU_SAFE_CALL_NO_SYNC(cuCtxDetach(cu_context_global));
    
    return 0;
}

