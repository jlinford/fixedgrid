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
    int i, j;
    
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

void write_conc(fixedgrid_t* G, uint32_t spec, char* sname, uint32_t iter, uint32_t proc)
{
    int i, j;
    float x, y;
    FILE *fptr;
    char fname[255];
    
    timer_start(&G->ppe_metrics.file_io);
    
    /* Build file name */
    sprintf(fname, "%s/OUT_solution_%s_%03d_%05ld.%03d", OUTPUT_DIR, sname, RUN_ID, iter, proc);

    /* Write to new file */
    if((fptr = (FILE*)fopen(fname, "w")) != NULL)
    {
        for(i=0; i<NROWS; i++)
        {
            for(j=0; j<NCOLS; j++)
            {
                x = DX*j + DX*0.5;
                y = DY*i + DY*0.5;
                fprintf(fptr, "%22.16E  %22.16E %22.16E\n", x, y, G->conc(spec, i, j));
            }
        }
        fclose(fptr);
    }
    else
    {
        fprintf(stderr, "Couldn't open file \"%s\" for writing.", fname);
        exit(1);
    }
    
    timer_stop(&G->ppe_metrics.file_io);
}

void write_metrics_to_csv_file(volatile metrics_t* m, FILE* fptr)
{
    int i;
    
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
    int i;
    
    FILE* fptr;
    char fname[255];
    metrics_t avg;
    
    // Build file name
    sprintf(fname, "%s/METRICS_%03d.csv", OUTPUT_DIR, RUN_ID);
    
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
        fprintf(fptr, "Steps:,%d,\n", G->steps);
        fprintf(fptr, ",\n");
        
        // Write PPE metrics
        write_metrics_to_csv_file(&G->ppe_metrics, fptr);
        
        // Write averages
        get_avg_spe_metrics(G, &avg);
        write_metrics_to_csv_file(&avg, fptr);
                
        // Write SPE metrics
        for(i=0; i<G->nprocs; i++)
        {
            printf("Writing metrics %d...", i);
            write_metrics_to_csv_file(&G->threads[i].metrics, fptr);
            printf(" done.\n");
        }
        
        printf("Flushing...");
        fclose(fptr);
        printf(" done.\n");
    }
    else
    {
        fprintf(stderr, "Couldn't open file \"%s\" for writing.", fname);
        exit(1);
    }
    
    printf("Metrics have been stored to file: %s\n", fname);
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
