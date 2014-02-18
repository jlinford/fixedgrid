#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "params.h"
#include "timer.h"
#include "fileio.h"

void write_conc(volatile double *conc, long iter, int proc)
{
    int i, j;
    double x, y;
    FILE *fptr;
    char fname[255];
    
    timer_start(TIMER_FIO);
    
    /* Build file name */
    sprintf(fname, "%s/OUT_solution_%03d_%05ld.%03d", OUTPUT_DIR, RUN_ID, iter, proc);

    /* Write to new file */    
    if((fptr = (FILE*)fopen(fname, "w+")) != NULL)
    {
        for(i=0; i<NROWS; i++)
        {
            for(j=0; j<NCOLS; j++)
            {
                x = DX*i + DX*0.5;
                y = DY*j + DY*0.5;
                fprintf(fptr, "%22.16E  %22.16E %22.16E\n", x, y, conc[i*NX+j]);
            }
        }
        fclose(fptr);
    }
    else
    {
        fprintf(stderr, "Couldn't open file \"%s\" for writing.", fname);
        exit(1);
    }
    
    timer_stop(TIMER_FIO);
}

