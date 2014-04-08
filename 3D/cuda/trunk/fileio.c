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

void write_conc(fixedgrid_t* G, uint32_t iter, uint32_t proc)
{
    uint32_t x, y, z, k;
    uint32_t spc_ind;
    float coord_x, coord_y, coord_z;
    FILE *fptr;
    char fname[255];
    
    timer_start(&G->metrics.file_io);
    
    for(k=0; k<NMONITOR; k++)
    {
        spc_ind = MONITOR[k];
        
        /* Build file name */
        sprintf(fname, "%s/OUT_solution_%s_cuda_%05d.%03d", OUTPUT_DIR, SPC_NAMES[spc_ind], iter, proc);

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
                                G->conc(x, y, z, MONITOR[k]));
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
    timer_stop(&G->metrics.file_io);
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
    uint32_t steps;
    
    FILE* fptr;
    char fname[255];
    
    timer_start(&G->metrics.file_io);
    
    steps = (G->tend - G->tstart) / G->dt;
    
    // Build file name
    sprintf(fname, "%s/METRICS_%03d_cuda.csv", OUTPUT_DIR, RUN_ID);
    
    // Write to new file
    if((fptr = (FILE*)fopen(fname, "w")) != NULL)
    {
        // Write header
        fprintf(fptr, ",\n");
        fprintf(fptr, "Platform:,%s,\n", platform);
        fprintf(fptr, "NPROCS:,%d,\n", 1);
        fprintf(fptr, "NY:,%d,\n", NY);
        fprintf(fptr, "NX:,%d,\n", NX);
        fprintf(fptr, "NSPEC:,%d,\n", NSPEC);
        fprintf(fptr, "Steps:,%d,\n", steps);
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
    
    printf("Metrics stored to file: %s\n", fname);
    
    timer_stop(&G->metrics.file_io);
}
