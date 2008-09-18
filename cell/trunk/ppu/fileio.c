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
#include "saprc99_Monitor.h"

void write_conc(fixedgrid_t* G, uint32_t iter, uint32_t proc)
{
    uint32_t x, y, z, s;
    uint32_t spc_ind;
    float coord_x, coord_y, coord_z;
    FILE *fptr;
    char fname[255];
    
    timer_start(&G->ppe_metrics.file_io);
    
    for(s=0; s<NMONITOR; s++)
    {
        spc_ind = MONITOR[s];
        
        /* Build file name */
        sprintf(fname, "%s/OUT_solution_%s_%02d_%05d.%03d", OUTPUT_DIR, SPC_NAMES[spc_ind], G->nprocs, iter, proc);
        
        /* Write to new file */
        if((fptr = (FILE*)fopen(fname, "w")) != NULL)
        {
            for(z=0; z<NZ; z++)
            {
                for(y=0; y<NY; y++)
                {
                    for(x=0; x<NX; x++)
                    {
                        coord_x = DX*x + DX*0.5;
                        coord_y = DY*y + DY*0.5;
                        coord_z = DZ*z + DZ*0.5;
                        fprintf(fptr, "%22.16E %22.16E %22.16E %22.16E\n", 
                                coord_x, coord_y, coord_z, 
                                G->conc(x, y, z, spc_ind));
                    }
                }
            }
            fclose(fptr);
        }
        else
        {
            fprintf(stderr, "Couldn't open file \"%s\" for writing.", fname);
            exit(1);
        }
    }    
    timer_stop(&G->ppe_metrics.file_io);
}

void write_metrics_to_csv_file( metrics_t* m, FILE* fptr)
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

metrics_t average_metrics(fixedgrid_t* G)
{
    uint32_t i, j;
    metrics_t m;
    stopwatch_t *tptr1, *tptr2;
    
    metrics_init(&m, "SPE Average");
    
    for(i=0; i<G->nprocs; i++)
    {
        tptr1 = (stopwatch_t*)(&m);
        tptr2 = (stopwatch_t*)(&G->threads[i].metrics);
        for(j=0; j<NUM_TIMERS; j++, tptr1++, tptr2++)
        {
            tptr1->elapsed += tptr2->elapsed;
        }
    }
    
    tptr1 = (stopwatch_t*)(&m);
    for(i=0; i<NUM_TIMERS; i++, tptr1++)
    {
        tptr1->elapsed /= G->nprocs;
    }
    
    return m;
}

void write_metrics_as_csv(fixedgrid_t* G, char* platform)
{
    uint32_t i, steps;
    metrics_t avg;
    
    FILE* fptr;
    char fname[255];
    
    steps = (G->tend - G->tstart) / G->dt;
    
    // Build file name
    sprintf(fname, "%s/METRICS_%03d_%02d.csv", OUTPUT_DIR, RUN_ID, G->nprocs);
    
    // Write to new file
    if((fptr = (FILE*)fopen(fname, "w")) != NULL)
    {
        // Write header
        fprintf(fptr, "Platform:,%s,\n", platform);
        fprintf(fptr, "NPROCS:,%d,\n", G->nprocs);
        fprintf(fptr, "NX:,%d,\n", NX);
        fprintf(fptr, "NY:,%d,\n", NY);
        fprintf(fptr, "NZ:,%d,\n", NZ);
        fprintf(fptr, "NSPEC:,%d,\n", NSPEC);
        fprintf(fptr, "Steps:,%d,\n", steps);
        fprintf(fptr, ",\n");
        
        // Write metrics
        write_metrics_to_csv_file(&G->ppe_metrics, fptr);
        
        avg = average_metrics(G);
        write_metrics_to_csv_file(&avg, fptr);

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

