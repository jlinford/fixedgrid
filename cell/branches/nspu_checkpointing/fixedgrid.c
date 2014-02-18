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
#include <libspe2.h>

#include "fixedgrid.h"
#include "saprc99.h"
#include "spe_pthread.h"

/* Program data */
static fixedgrid_t G_GLOBAL;

/**
 * Process rows 1/2 timestep 
 */
void discretize_row(fixedgrid_t* G, uint32_t spec, value_t dt)
{
#ifdef DO_ROW_DISCRET
    uint32_t i;
    
    /* Loop blocking */
    uint32_t block;

    timer_start(&G->ppe_metrics.comm);
    
    block = NROWS / G->nprocs;
    for(i=0; i<G->nprocs; i++)
    {
        /* Configure SPE arguments */
        G->threads[i].argv.arg[0].u64 = (uint64_t)(&G->conc(spec, i*block, 0));
        G->threads[i].argv.arg[1].u64 = (uint64_t)(&G->wind_u[i*block*NCOLS]);
        G->threads[i].argv.arg[2].u64 = (uint64_t)(&G->diff[i*block*NCOLS]);
        G->threads[i].argv.arg[3].dbl = dt;
        G->threads[i].argv.arg[4].dbl = DX;
        G->threads[i].argv.arg[5].u32[0] = NCOLS;
        G->threads[i].argv.arg[5].u32[1] = (i == G->nprocs - 1 ? block + NROWS - G->nprocs*block : block);  //FIXME: Potential load imbalance
        G->threads[i].argv.arg[6].dbl = G->temp;
        G->threads[i].argv.arg[7].dbl = 0.0;
        
        /* Signal SPE */
        spe_set_status(G, i, SPE_STATUS_ROW);
    }
    
    timer_stop(&G->ppe_metrics.comm);
    
    /* Wait for SPEs to finish */
    wait_all_spes(G);
    
#endif
}

/**
 * Process columns 1 timestep 
 */
void discretize_col(fixedgrid_t* G, uint32_t spec, value_t dt)
{
#ifdef DO_COL_DISCRET
    uint32_t i;
    
    /* Loop blocking */
    uint32_t block;

    timer_start(&G->ppe_metrics.comm);
    
    /* Multiples of block must give addresses on a 16-byte boundary */
    block = NCOLS / G->nprocs;
    block -= block % (16 / sizeof(value_t)); //FIXME: broken for value_t != double
    
    for(i=0; i<G->nprocs; i++)
    {
        /* Configure SPE arguments */
        G->threads[i].argv.arg[0].u64 = (uint64_t)(&G->conc(spec, 0, i*block));
        G->threads[i].argv.arg[1].u64 = (uint64_t)(&G->wind_v[i*block]);
        G->threads[i].argv.arg[2].u64 = (uint64_t)(&G->diff[i*block]);
        G->threads[i].argv.arg[3].dbl = dt;
        G->threads[i].argv.arg[4].dbl = DY;
        G->threads[i].argv.arg[5].u32[0] = NROWS;
        G->threads[i].argv.arg[5].u32[1] = (i == G->nprocs - 1 ? block + NCOLS - G->nprocs*block : block);  //FIXME: Potential load imbalance
        G->threads[i].argv.arg[6].dbl = G->temp;
        G->threads[i].argv.arg[7].dbl = 0.0;
        
        /* Signal SPE */
        spe_set_status(G, i, SPE_STATUS_COL);
    }
    
    timer_stop(&G->ppe_metrics.comm);
    
    /* Wait for SPEs to finish */
    wait_all_spes(G);

#endif
}

/**
 * Initializes the model
 */
void init_model(fixedgrid_t* G)
{
    /* Iterators */
    uint32_t i, j, k;
    
    /* Chemistry buffer */
    value_t chemBuff[NSPEC];
    
    /* Initialize metrics */
    metrics_init(&G->ppe_metrics, "PPE");
    
    /* Allocate memory */
    printf("Allocating memory... ");
    allocate_global_memory();
    printf("done.\n");
    
    /* Initialize chemistry and concentration data */
    printf("Loading chemistry and concentration data... ");
    timer_start(&G->ppe_metrics.array_init);
    init_saprc99();
    for(i=0; i<NROWS; i++)
    {
        for(j=0; j<NCOLS; j++)
        {
            #ifdef DO_CHEMISTRY
            
            for(k=0; k<NSPEC; k++)
            {
                chemBuff[k] = G->conc(k, i, j);
            }
            
            init_conc(chemBuff);
            
            #else
            
            for(k=0; k<NSPEC; k++)
            {
                chemBuff[k] = O3_INIT;
            }
            
            #endif
            
            for(k=0; k<NSPEC; k++)
            {
                G->conc(k, i, j) = chemBuff[k];
            }
        }
    }
    timer_stop(&G->ppe_metrics.array_init);
    printf("done.\n");
    
    /* Initialize wind field */
    printf("Loading wind field data... ");
    timer_start(&G->ppe_metrics.array_init);
    for(i=0; i<NROWS; i++)
    {
        for(j=0; j<NCOLS; j++)
        {
            G->wind_u[i*NCOLS + j] = WIND_U_INIT;
            G->wind_v[i*NCOLS + j] = WIND_V_INIT;
        }
    }
    timer_stop(&G->ppe_metrics.array_init);
    printf("done.\n");
    
    /* Initialize diffusion field */
    printf("Loading diffusion field data... ");
    timer_start(&G->ppe_metrics.array_init);
    for(i=0; i<NROWS; i++)
    {
        for(j=0; j<NCOLS; j++)
        {
            G->diff[i*NCOLS + j] = DIFF_INIT;
        }
    }
    timer_stop(&G->ppe_metrics.array_init);
    printf("done.\n");
    
    /* Initialize time */
    G->end_time =  year2sec(END_YEAR - START_YEAR) + day2sec(END_DOY - START_DOY) + 
                   hour2sec(END_HOUR - START_HOUR) + minute2sec(END_MIN - START_MIN);
    G->dt = STEP_SIZE;
    G->steps = (uint32_t)(G->end_time/G->dt);
    
    /* Initialize environment */
    G->temp = TEMP_INIT;
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
    G->conc(ind_O3, SOURCE_Y, SOURCE_X) += SOURCE_RATE / (DX * DY * 1000.0);
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
void print_start_banner(const value_t x, const value_t y, const value_t z, 
                        const uint32_t t_end, const uint32_t steps)
{    
    printf("\n");
    printf("SPACE DOMAIN:\n");
    printf("    LENGTH (X): %f meters\n", x);
    printf("    WIDTH  (Y): %f meters\n", y);
    printf("    DEPTH  (Z): %f meters\n", z);
    printf("     # Species: %d\n", NSPEC);
    printf("\n");
    printf("\n");
    printf("TIME DOMAIN:\n");
    printf("    FROM  %d:%d.00 on day %d of year %d\n", START_HOUR, START_MIN, START_DOY, START_YEAR);
    printf("    TO    %d:%d.00 on day %d of year %d\n", END_HOUR, END_MIN, END_DOY, END_YEAR);
    printf("    TOTAL %ld seconds (%ld timesteps of %d seconds)\n", t_end, steps, STEP_SIZE);
    
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
    uint32_t i, k, iter;
    
    /* Start wall clock timer */
    timer_start(&G->ppe_metrics.wallclock);
    
    /* Calculate available SPEs */
    G->nprocs = spe_cpu_info_get(SPE_COUNT_USABLE_SPES, -1);
    G->nprocs = G->nprocs > SPE_MAX_THREADS ? SPE_MAX_THREADS : G->nprocs;
    
    /* Parse command line arguments */
    if(argc > 1)
    {
        i = atoi(argv[1]);
        if(i < 1)
        {
            fprintf(stderr, "Invalid number of SPUs: %d < 1.\n", i);
            exit(1);
        }
        
        if(i < G->nprocs)
        {
            G->nprocs = i;
        }
        else 
        {
            printf("%d SPUs unavailable.  Using %d instead.\n", i, G->nprocs);
        }
    }
    
    /* Create SPE threads */
    create_spe_pthreads(G);
    
    printf("\nRunning %d threads (%d SPU + 1 PPU).\n", (G->nprocs+1), G->nprocs);
    
    /* Initialize the model parameters */
    init_model(G);

    /* Add emissions */
    process_emissions(G);
    
    /* Print startup banner */
    print_start_banner(NX*DX, NY*DY, 0.0, G->end_time, G->steps);
    
    /* Store initial concentration */
    printf("Writing initial concentration data... ");
    write_conc(G, ind_O3, "O3", 0, 0);
    printf("done.\n");
    
    /* BEGIN CALCULATIONS */
    for(iter = 1; iter <= G->steps; iter++)
    {
        k = ind_O3;
        //for(k=0; k<NSPEC; k++)
        //{
            discretize_row(G, k, G->dt/2.0);
          
            discretize_col(G, k, G->dt);
            
            discretize_row(G, k, G->dt/2.0);
        //}

        update_model(G);
        
        #ifdef WRITE_EACH_ITER
        write_conc(G, ind_O3, "O3", iter, 0);
        #endif

        printf("Iteration %d of %d.  Time = %f seconds.\n", iter, G->steps, iter*G->dt);
    }
    /* END CALCULATIONS */
    
    /* Wait for SPU-thread to complete execution. */
    join_all_spes(G);
    
    /* Store concentration */
    #ifndef WRITE_EACH_ITER
    write_conc(G, ind_O3, "O3", iter-1, 0);
    #endif
    
    /* Show final time */
    printf("\nFinal time: %f seconds.\n", (iter-1)*G->dt);
    
    timer_stop(&G->ppe_metrics.wallclock);
    
    /* Print metrics */
    print_metrics(&G->ppe_metrics);
    
    metrics_t avg;
    get_avg_spe_metrics(G, &avg);
    print_metrics(&avg);
    
    for(i=0; i<G->nprocs; i++)
    {
        print_metrics(&G->threads[i].metrics);
    }
    
    /* Write metrics to CSV file */
    write_metrics_as_csv(G, "Cell B.E.");
    
    /* Cleanup and exit */
    free_global_memory(G);
    return 0;
}

