/*
 *  fileio.c
 *  
 *  File I/O functions.
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "fileio.h"
#include "timer.h"
#include "params.h"
#include "common.h"
#include "saprc99_Monitor.h"

static char* timer_names[NUM_TIMERS] = 
{
    "Wallclock    ",
    "Array Init   ",
    "Array Copy   ",
    "File I/O     ",
    "Row discret  ",
    "Col discret  ",
    "Total discret",
    "Chemistry    ",
    "PPU/SPU Comm "
};

void get_avg_spe_metrics(fixedgrid_t* G, metrics_t* avg)
{
    uint32_t i, j;
    
    stopwatch_t *tptr1, *tptr2;
    
    metrics_init(avg, "SPE Average");

    for(i=0; i<G->nprocs; i++)
    {
        tptr1 = (stopwatch_t*)(avg);
        tptr2 = (stopwatch_t*)(&G->threads[i].metrics);
        for(j=0; j<NUM_TIMERS; j++, tptr1++, tptr2++)
        {
            tptr1->elapsed += tptr2->elapsed;
        }
    }
    
    tptr1 = (stopwatch_t*)(avg);
    for(i=0; i<NUM_TIMERS; i++, tptr1++)
    {
        tptr1->elapsed /= G->nprocs;
    }
}

void write_conc(fixedgrid_t* G, uint32_t iter, uint32_t proc)
{
    uint32_t i, j, k;
    uint32_t spc_ind;
    float x, y;
    FILE *fptr;
    char fname[255];
    
    timer_start(&G->ppe_metrics.file_io);
    
    for(k=0; k<NMONITOR; k++)
    {
        spc_ind = MONITOR[k];
        
        /* Build file name */
        sprintf(fname, "%s/OUT_solution_%s_%02d_%05d.%03d", OUTPUT_DIR, SPC_NAMES[spc_ind], G->nprocs, iter, proc);

        /* Write to new file */
        if((fptr = (FILE*)fopen(fname, "w")) != NULL)
        {
            //printf("Saving %s concentrations to file \"%s\"...", SPC_NAMES[spc_ind], fname);
            for(i=0; i<NROWS; i++)
            {
                for(j=0; j<NCOLS; j++)
                {
                    x = DX*j + DX*0.5;
                    y = DY*i + DY*0.5;
                    fprintf(fptr, "%22.16E  %22.16E %22.16E\n", x, y, G->conc(MONITOR[k], i, j));
                }
            }
            fclose(fptr);
            //printf(" done.\n");
        }
        else
        {
            fprintf(stderr, "Couldn't open file \"%s\" for writing.", fname);
            exit(1);
        }
    }    
    timer_stop(&G->ppe_metrics.file_io);
}

void write_metrics_to_csv_file(volatile metrics_t* m, FILE* fptr)
{
    uint32_t i;
    
    stopwatch_t* tptr = (stopwatch_t*)m;
    
    fprintf(fptr, "Timer,%s,\n", m->name);
    
    for(i=0; i<NUM_TIMERS; i++, tptr++)
    {
        fprintf(fptr, "%s,%f,\n", timer_names[i], tptr->elapsed);
    }
    fprintf(fptr, ",\n,\n");
}

void write_metrics_as_csv(fixedgrid_t* G, char* platform)
{
    uint32_t i;
    uint32_t steps;
    
    FILE* fptr;
    char fname[255];
    metrics_t avg;
    
    steps = (G->tend - G->tstart) / G->dt;
    
    // Build file name
    sprintf(fname, "%s/METRICS_%03d_%02d.csv", OUTPUT_DIR, RUN_ID, G->nprocs);
    
    // Write to new file
    if((fptr = (FILE*)fopen(fname, "w")) != NULL)
    {
        // Write header
        fprintf(fptr, ",\n");
        fprintf(fptr, "Platform:,%s,\n", platform);
        fprintf(fptr, "NPROCS:,%d (%d SPE + %d PPE),\n", G->nprocs+1, G->nprocs, 1);
        fprintf(fptr, "NROWS:,%d,\n", NROWS);
        fprintf(fptr, "NCOLS:,%d,\n", NCOLS);
        fprintf(fptr, "NSPEC:,%d,\n", NSPEC);
        fprintf(fptr, "Steps:,%d,\n", steps);
        fprintf(fptr, ",\n");
        
        // Write PPE metrics
        write_metrics_to_csv_file(&G->ppe_metrics, fptr);
        
        // Write averages
        get_avg_spe_metrics(G, &avg);
        write_metrics_to_csv_file(&avg, fptr);
                
        // Write SPE metrics
        for(i=0; i<G->nprocs; i++)
        {
            write_metrics_to_csv_file(&G->threads[i].metrics, fptr);
        }
        
        fclose(fptr);
    }
    else
    {
        fprintf(stderr, "Couldn't open file \"%s\" for writing.", fname);
        exit(1);
    }
    
    printf("Metrics stored to file: %s\n", fname);
}

void print_metrics(volatile metrics_t* m)
{
    int i;
    stopwatch_t* tptr = (stopwatch_t*)m;
    
    printf("\n===== %s =====\n", m->name);
    
    for(i=0; i<NUM_TIMERS; i++, tptr++)
    {
        printf("%s: %f\n", timer_names[i], tptr->elapsed);
    }
}