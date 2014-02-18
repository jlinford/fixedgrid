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

void write_conc(fixedgrid_t* G, uint32_t spec, char* sname, uint32_t iter, uint32_t proc)
{
    int i, j;
    double x, y;
    FILE *fptr;
    char fname[255];
    
    timer_start(&G->metrics.file_io);
    
    /* Build file name */
    sprintf(fname, "%s/OUT_solution_%s_%03d_%05d.%03d", OUTPUT_DIR, sname, RUN_ID, iter, proc);

    /* Write to new file */
    if((fptr = (FILE*)fopen(fname, "w")) != NULL)
    {
        for(i=0; i<NROWS; i++)
        {
            for(j=0; j<NCOLS; j++)
            {
                x = DX*j + DX*0.5;
                y = DY*i + DY*0.5;
                fprintf(fptr, "%22.16E  %22.16E %22.16E\n", x, y, G->conc[spec][i][j]);
            }
        }
        fclose(fptr);
    }
    else
    {
        fprintf(stderr, "Couldn't open file \"%s\" for writing.", fname);
        exit(1);
    }
    
    timer_stop(&G->metrics.file_io);
}

void write_metrics_to_csv_file(metrics_t* m, FILE* fptr)
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
    FILE* fptr;
    char fname[255];
    
    // Build file name
    sprintf(fname, "%s/METRICS_%03d.csv", OUTPUT_DIR, RUN_ID);
    
    // Write to new file
    if((fptr = (FILE*)fopen(fname, "w")) != NULL)
    {
        // Write header
        fprintf(fptr, ",\n");
        fprintf(fptr, "Platform:,%s,\n", platform);
        fprintf(fptr, "NPROCS:,%d,\n", G->nprocs);
        fprintf(fptr, "NROWS:,%d,\n", NROWS);
        fprintf(fptr, "NCOLS:,%d,\n", NCOLS);
        fprintf(fptr, "NSPEC:,%d,\n", NSPEC);
        fprintf(fptr, "Steps:,%d,\n", G->steps);
        fprintf(fptr, ",\n");
        
        // Write PPE metrics
        write_metrics_to_csv_file(&G->metrics, fptr);
        
        fclose(fptr);
    }
    else
    {
        fprintf(stderr, "Couldn't open file \"%s\" for writing.", fname);
        exit(1);
    }
    
    printf("Metrics have been stored to file: %s\n", fname);
}

void print_metrics(metrics_t* m)
{
    int i;
    stopwatch_t* tptr = (stopwatch_t*)m;
    
    printf("\n===== %s =====\n", m->name);
    
    for(i=0; i<NUM_TIMERS; i++, tptr++)
    {
        printf("%s: %f\n", timer_names[i], tptr->elapsed);
    }
}
