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

#include "fileio.h"
#include "fixedgrid.h"
#include "memory.h"
#include "spe_pthread.h"
#include "saprc99_Monitor.h"

/* Program data */
fixedgrid_t G_GLOBAL;

/**
 * Applies KPP-generated SAPRC99 chemical mechanism
 * to all concentration data.
 */
void start_saprc99(fixedgrid_t* G)
{
#if DO_CHEMISTRY == 1
    uint32_t i;
    
    /* Loop blocking */
    uint32_t block;
    
    timer_start(&G->ppe_metrics.chem);
    
    block = NROWS / G->nprocs;
    for(i=0; i<G->nprocs; i++)
    {
        /* Configure SPE arguments */
        G->threads[i].argv.arg[0].u64 = (uint32_t)(&G->conc(0, i*block, 0));
        G->threads[i].argv.arg[1].dbl = G->time;
        G->threads[i].argv.arg[2].dbl = G->dt;
        G->threads[i].argv.arg[3].u32[0] = (i == G->nprocs - 1 ? block + NROWS - G->nprocs*block : block);  //FIXME: Potential load imbalance
        
        /* Signal SPE */
        spe_set_status(G, i, SPE_STATUS_CHEM);
    }
    
    /* Wait for SPEs to finish */
    wait_all_spes(G);
    
    timer_stop(&G->ppe_metrics.chem);
#endif
}

/**
 * Process rows 1/2 timestep 
 */
void start_discretize_row(fixedgrid_t* G, uint32_t spec, real_t dt)
{
#if DO_ROW_DISCRET == 1
    uint32_t i;
    
    /* Loop blocking */
    uint32_t block;
    
    timer_start(&G->ppe_metrics.row_discret);
    
    block = NROWS / G->nprocs;
    for(i=0; i<G->nprocs; i++)
    {
        /* Configure SPE arguments */
        G->threads[i].argv.arg[0].u64 = (uint32_t)(&G->conc(spec, i*block, 0));
        G->threads[i].argv.arg[1].u64 = (uint32_t)(&G->wind_u(i*block, 0));
        G->threads[i].argv.arg[2].u64 = (uint32_t)(&G->diff(i*block, 0));
        G->threads[i].argv.arg[3].dbl = dt;
        G->threads[i].argv.arg[4].dbl = DX;
        G->threads[i].argv.arg[5].u32[0] = NCOLS;
        G->threads[i].argv.arg[5].u32[1] = (i == G->nprocs - 1 ? block + NROWS - G->nprocs*block : block);  //FIXME: Potential load imbalance
        
        /* Signal SPE */
        spe_set_status(G, i, SPE_STATUS_ROW);
    }
    
    /* Wait for SPEs to finish */
    wait_all_spes(G);
    
    timer_stop(&G->ppe_metrics.row_discret);
#endif
}

/**
 * Process columns 1 timestep 
 */
void start_discretize_col(fixedgrid_t* G, uint32_t spec, real_t dt)
{
#if DO_COL_DISCRET == 1
    uint32_t i;
    uint32_t col, spe;
    
    timer_start(&G->ppe_metrics.col_discret);
    
    for(col=0; col<NCOLS; col++)
    {
        spe = col % G->nprocs;
        
        /* Buffer has been filled with results by now, so write back */
        if(col >= G->nprocs)
        {
            wait_for_spe(G, spe);
            
            for(i=0; i<NROWS; i++)
            {
                G->conc(spec, i, col-G->nprocs) = G->conc_T(spe, i);
            }
        }
        
        /* Transpose data so SPE doesn't use dma lists */
        for(i=0; i<NROWS; i++)
        {
            G->conc_T(spe, i) = G->conc(spec, i, col);
            G->wind_T(spe, i) = G->wind_v(i, col);
            G->diff_T(spe, i) = G->diff(i, col);
        }
                
        /* Configure SPE arguments */
        G->threads[spe].argv.arg[0].u64 = (uint32_t)(&G->conc_T(spe, 0));
        G->threads[spe].argv.arg[1].u64 = (uint32_t)(&G->wind_T(spe, 0));
        G->threads[spe].argv.arg[2].u64 = (uint32_t)(&G->diff_T(spe, 0));
        G->threads[spe].argv.arg[3].dbl = dt;
        G->threads[spe].argv.arg[4].dbl = DY;
        G->threads[spe].argv.arg[5].u32[0] = NROWS;
        G->threads[spe].argv.arg[5].u32[1] = 1;
        
        /* Do transport as rows, now columns are transposed */
        spe_set_status(G, spe, SPE_STATUS_ROW);
    }

    /* Flush buffers */
    for(col=NCOLS-G->nprocs; col<NCOLS; col++)
    {
        spe = col % G->nprocs;
        
        wait_for_spe(G, spe);
    
        for(i=0; i<NROWS; i++)
        {
            G->conc(spec, i, col) = G->conc_T(spe, i);
        }
    }
        
    timer_stop(&G->ppe_metrics.col_discret);
#endif
}

/**
 * Initializes the model
 */
void init_model(fixedgrid_t* G)
{
    /* Iterators */
    uint32_t i, j, k;
    
    /* Loop blocking */
    uint32_t block;
    
    /* Chemistry buffer */
    real_t chemBuff[NSPEC];
    
    /* Initialize metrics */
    metrics_init(&G->ppe_metrics, "PPE");
    
    /* Allocate memory */
    printf("\nAllocating memory... ");
    allocate_global_memory(G);
    printf("done.\n");
    
    /* Initialize time frame */
    /* FIXME: year is ignored */
    G->tstart = day2sec(START_DOY) + hour2sec(START_HOUR) + minute2sec(START_MIN);
    G->tend   = day2sec(END_DOY)   + hour2sec(END_HOUR)   + minute2sec(END_MIN);
    G->dt = STEP_SIZE;
    G->time = G->tstart;
    
    /* Pack model parameters to be collected by SPEs */
    block = NROWS / G->nprocs;
    for(i=0; i<G->nprocs; i++)
    {
        G->threads[i].argv.arg[0].u64 = (uint32_t)(&G->conc(0, 0, 0));
        G->threads[i].argv.arg[1].u64 = (uint32_t)(&G->wind_u(0, 0));
        G->threads[i].argv.arg[2].u64 = (uint32_t)(&G->wind_v(0, 0));
        G->threads[i].argv.arg[3].u64 = (uint32_t)(&G->diff(0, 0));
        G->threads[i].argv.arg[4].dbl = G->tstart;
        G->threads[i].argv.arg[5].dbl = G->tend;
        G->threads[i].argv.arg[6].dbl = G->dt;
        G->threads[i].argv.arg[7].dbl = TEMP_INIT;
        G->threads[i].argv.arg[8].u32[0] = (i == G->nprocs - 1 ? block + NROWS - G->nprocs*block : block);  //FIXME: Potential load imbalance
    }
    
    /* Initialize chemistry and concentration data */
    printf("Loading chemistry and concentration data... ");
    timer_start(&G->ppe_metrics.array_init);

    #if DO_CHEMISTRY == 1
    
    saprc99_Initialize(chemBuff);
    
    for(k=0; k<NSPEC; k++)
    {
        for(i=0; i<NROWS; i++)
        {
            for(j=0; j<NCOLS; j++)
            {
                G->conc(k, i, j) = chemBuff[k];
            }
        }
    }
    
    #else
    
    for(i=0; i<NROWS; i++)
    {
        for(j=0; j<NCOLS; j++)
        {
            G->conc(ind_O3, i, j) = O3_INIT;
        }
    }
    
    #endif
    
    timer_stop(&G->ppe_metrics.array_init);
    printf("done.\n");
    
    /* Initialize wind field */
    printf("Loading wind field data... ");
    timer_start(&G->ppe_metrics.array_init);
    
    for(i=0; i<NROWS; i++)
    {
        for(j=0; j<NCOLS; j++)
        {
            G->wind_u(i, j) = WIND_U_INIT;
            G->wind_v(i, j) = WIND_V_INIT;
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
            G->diff(i, j) = DIFF_INIT;
        }
    }
    
    timer_stop(&G->ppe_metrics.array_init);
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
void print_start_banner(fixedgrid_t* G)
{
    int i;
    uint32_t steps;
    
    steps = (G->tend - G->tstart) / G->dt;
    
    printf("\n");
    printf("CONFIGURATION:\n");
    printf("    ROW DISCRETIZATION:    %s\n", DO_ROW_DISCRET == TRUE ? "TRUE" : "FALSE");
    printf("    COLUMN DISCRETIZATION: %s\n", DO_COL_DISCRET == TRUE ? "TRUE" : "FALSE");
    printf("    SAPRC99 CHEMISTRY:     %s\n", DO_CHEMISTRY == TRUE ? "TRUE" : "FALSE");
    printf("    DOUBLE PRECISION:      %s\n", DOUBLE_PRECISION == TRUE ? "TRUE" : "FALSE");
    printf("\n");
    printf("SPACE DOMAIN:\n");
    printf("    LENGTH (X): %f meters\n", NCOLS*DX);
    printf("    WIDTH  (Y): %f meters\n", NROWS*DY);
    printf("    DEPTH  (Z): %f meters\n", 1.0);
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
        
        if(i > G->nprocs)
        {
            printf("%d SPUs unavailable.  Using %d instead.\n", i, G->nprocs);
        }
        else 
        {
            G->nprocs = i;
        }
    }
    
    /* Check dimensions */
    if(NROWS < 5)
    {
        fprintf(stderr, "%d rows < 5 rows is too small for discretization.\n", NROWS);
    }
    if(NCOLS < 5)
    {
        fprintf(stderr, "%d columns < 5 columns is too small for discretization.\n", NCOLS);
    }
    
    /* Don't use more SPEs than there are rows or columns */
    if(NROWS < G->nprocs)
    {
        printf("%d SPUs available, but only %d rows, so using %d SPUs\n", G->nprocs, NROWS, NROWS);
        G->nprocs = NROWS;
    }
    if(NCOLS / VECTOR_LENGTH < G->nprocs)
    {
        printf("%d SPUs available, but only %d column vectors of size %d, so using %d SPUs\n", G->nprocs, (NCOLS/VECTOR_LENGTH), VECTOR_LENGTH, (NCOLS/VECTOR_LENGTH));
        G->nprocs = (NCOLS/VECTOR_LENGTH);
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
    
    /* Create SPE threads */
    create_spe_pthreads(G);
    
    /* Wait for SPEs to finish initialization */
    wait_all_spes(G);
    
    printf("\nRunning %d threads (%d SPU + 1 PPU).\n\n", (G->nprocs+1), G->nprocs);    
        
    /* BEGIN CALCULATIONS */
    for(iter=1, G->time = G->tstart; G->time < G->tend; G->time += G->dt, ++iter)
    {
        start_saprc99(G);
        
        for(k=0; k<NLOOKAT; k++)
        {
            start_discretize_row(G, LOOKAT[k], G->dt/2.0);
          
            start_discretize_col(G, LOOKAT[k], G->dt);
            
            start_discretize_row(G, LOOKAT[k], G->dt/2.0);
        }

        update_model(G);
        
        #if WRITE_EACH_ITER == 1
        write_conc(G, iter, 0);
        #endif

        printf("  After iteration %02d: Model time = %07.2f sec.\n", iter, iter*G->dt);
    }
    /* END CALCULATIONS */
    
    /* Wait for SPU-thread to complete execution. */
    join_all_spes(G);
    
    /* Store concentration */
    #if WRITE_EACH_ITER != 1
    write_conc(G, iter-1, 0);
    #endif
    
    /* Show final time */
    printf("\nFinal time: %f seconds.\n", (iter-1)*G->dt);
    
    timer_stop(&G->ppe_metrics.wallclock);
    
    /* Write metrics to CSV file */
    write_metrics_as_csv(G, "Cell B.E.");
    
    /* Cleanup and exit */
    free_global_memory(G);
    return 0;
}

