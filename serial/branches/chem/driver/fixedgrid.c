#include <stdio.h>
#include <stdlib.h>

#include "fixedgrid.h"
#include "params.h"
#include "timer.h"
#include "discretize.h"
#include "fileio.h"
#include "saprc99.h"

void double_array_init(unsigned int n, double *array, double val)
{
    int i;
    
    timer_start(TIMER_ARRAY_INIT);
    
    for(i=0; i<n; i++)
    {
        array[i] = val;
    }
    
    timer_stop(TIMER_ARRAY_INIT);
}

void double_array_copy(unsigned int n, double *src, double *dst)
{
    int i;
    
    timer_start(TIMER_ARRAY_COPY);
    
    for(i=0; i<n; i++)
    {
        dst[i] = src[i];
    }
    
    timer_stop(TIMER_ARRAY_COPY);
}


void print_emission_sources()
{
    printf("\n");
    printf("Emission sources (X, Y, Z, RATE):\n");
    printf("    (%f, %f, %f, %E)\n", (DX*SOURCE_X + DX*0.5), (DY*SOURCE_Y + DY*0.5), 0.0, SOURCE_RATE);
    printf("\n");
}

void print_start_banner(const double x, const double y, const double z, 
                        const long t_end, const long steps)
{    
    printf("\n");
    printf("SPACE DOMAIN:\n");
    printf("    LENGTH (X): %f meters\n", x);
    printf("    WIDTH  (Y): %f meters\n", y);
    printf("    DEPTH  (Z): %f meters\n", z);
    printf("\n");
    printf("\n");
    printf("TIME DOMAIN:\n");
    printf("    FROM  %d:%d.00 on day %d of year %d\n", START_HOUR, START_MIN, START_DOY, START_YEAR);
    printf("    TO    %d:%d.00 on day %d of year %d\n", END_HOUR, END_MIN, END_DOY, END_YEAR);
    printf("    TOTAL %ld seconds (%ld timesteps of %d seconds)\n", t_end, steps, STEP_SIZE);
    
    print_emission_sources();
    
    printf("\n");
}

int main(int argc, char** argv)
{
    /* Iterators */
    int i, j, k, s;
    
    /* Time (seconds) */
    long t_0;
    long t_end;
    long dt;
    long steps;
    long iter;
    double time;
    
    /* 2D concentration data */
    double conc[NROWS][NCOLS][NSPEC];
    double crow1[NX], crow2[NX];
    double ccol1[NY], ccol2[NY];
    
    /* 2D wind field data */
    double wind_u[NROWS][NCOLS];
    double wind_v[NROWS][NCOLS];
    double wrow[NX];
    double wcol[NY];

    /* 2D diffusion tensor data */
    double diff[NROWS][NCOLS];
    double drow[NX];
    double dcol[NY];
    
    /* Environment data */
    double temp;
        
    /* Emission control */
    bool emflag = TRUE;
    
    /* Start wall clock timer */
    timer_start(TIMER_WALLCLOCK);
    
    printf("\nRunning on 1 CPU\n");
    
    /* Initialize chemistry mechanism */
    init_saprc99();
    
    /* Initialize concentration data */
    for(i=0; i<NROWS; i++)
    {
        for(j=0; j<NCOLS; j++)
        {
            init_conc(conc[i][j]);
        }
    }

    /* Initialize wind field */
    double_array_init(NROWS*NCOLS, &(wind_u[0][0]), WIND_U_INIT);
    double_array_init(NROWS*NCOLS, &(wind_v[0][0]), WIND_V_INIT);
    
    /* Initialize diffusion field */
    double_array_init(NROWS*NCOLS, &(diff[0][0]), DIFF_INIT);
    
    /* Initialize time */
    t_0   = year2sec(START_YEAR) + day2sec(START_DOY) + 
            hour2sec(START_HOUR) + minute2sec(START_MIN);
            
    t_end = year2sec(END_YEAR) + day2sec(END_DOY) + 
            hour2sec(END_HOUR) + minute2sec(END_MIN);
    dt = STEP_SIZE;
    steps = (long)( (t_end - t_0)/dt );
    time = t_0;
    
    /* Initialize environment */
    temp = TEMP_INIT;
    
    /* Print startup banner */
    print_start_banner(NX*DX, NY*DY, 0.0, t_end, steps);
    
    /* Store initial concentration */
    write_conc(conc, ind_O3, "O3", 0, 0);
    
    /* BEGIN CALCULATIONS */
    for(iter = 1; iter <= steps; iter++)
    {
        emflag = iter*dt < 6*3600.0 ? TRUE : FALSE;
        
        timer_start(TIMER_CHEM);
        
        /* Chemistry */
        for(i=0; i<NROWS; i++)
        {
            for(j=0; j<NCOLS; j++)
            {
                saprc99(conc[i][j], temp, time, time+dt, dt);
            }
        }
        
        timer_stop(TIMER_CHEM);
                
        for(s=0; s<NSPEC; s++)
        {
            timer_start(TIMER_ROW_DISCRET);
        
            /* Discretize rows 1/2 timestep */
            for(i=0; i<NROWS; i++)
            {
                //double_array_copy(NX, conc[i], crow1);
                timer_start(TIMER_ARRAY_COPY);
                for(k=0; k<NX; k++)
                {
                    crow1[k] = conc[i][k][s];
                }
                timer_stop(TIMER_ARRAY_COPY);
                double_array_copy(NX, wind_u[i], wrow);
                double_array_copy(NX, diff[i], drow);
                
                discretize(NX, crow1, wrow, drow, DX, dt/2, crow2);
    
                //double_array_copy(NX, crow2, conc[i]);
                timer_start(TIMER_ARRAY_COPY);
                for(k=0; k<NX; k++)
                    conc[i][k][s] = crow2[k];
                timer_stop(TIMER_ARRAY_COPY);
            }
            
            timer_stop(TIMER_ROW_DISCRET);
            
            timer_start(TIMER_COL_DISCRET);
            
            /* Discretize colums 1 timestep */
            for(j=0; j<NCOLS; j++)
            {
                timer_start(TIMER_ARRAY_COPY);
                for(k=0; k<NY; k++)
                {
                    ccol1[k] = conc[k][j][s];
                    wcol[k]  = wind_v[k][j];
                    dcol[k]  = diff[k][j];
                }
                timer_stop(TIMER_ARRAY_COPY);
                
                discretize(NY, ccol1, wcol, dcol, DY, dt, ccol2);
                
                timer_start(TIMER_ARRAY_COPY);
                for(k=0; k<NY; k++)
                    conc[k][j][s] = ccol2[k];
                timer_stop(TIMER_ARRAY_COPY);
            }
            
            timer_stop(TIMER_COL_DISCRET);
            
            timer_start(TIMER_ROW_DISCRET);
            
            /* Discretize rows 1/2 timestep */
            for(i=0; i<NROWS; i++)
            {
                //double_array_copy(NX, conc[i], crow1);
                timer_start(TIMER_ARRAY_COPY);
                for(k=0; k<NX; k++)
                {
                    crow1[k] = conc[i][k][s];
                }
                timer_stop(TIMER_ARRAY_COPY);
                double_array_copy(NX, wind_u[i], wrow);
                double_array_copy(NX, diff[i], drow);
                
                discretize(NX, crow1, wrow, drow, DX, dt/2, crow2);
    
                //double_array_copy(NX, crow2, conc[i]);
                timer_start(TIMER_ARRAY_COPY);
                for(k=0; k<NX; k++)
                    conc[i][k][s] = crow2[k];
                timer_stop(TIMER_ARRAY_COPY);
            }
            
            timer_stop(TIMER_ROW_DISCRET);
        }
        
        /*
         * Could update wind field here...
         */
         
        /*
         * Could update diffusion tensor here...
         */
         
        /*
         * Could update environment here...
         */
        
        /* Add emissions */
        if(emflag)
        {
            conc[SOURCE_Y][SOURCE_X][ind_O3] += dt * (SOURCE_RATE) / (DX * DY * 1000.0);
        }
        
        /* Store concentration */
        #ifdef WRITE_EACH_ITER
        write_conc(conc, ind_03, "O3", iter, 0);
        #endif
        
        /* Indicate progress */
        if(iter % 10 == 0)
        {
            printf("Iteration %ld of %ld.  Time = %ld seconds.\n", iter, steps, iter*dt);
        }
        
        time += dt;
        
    }
    /* END CALCULATIONS */
    
    /* Store concentration */
    #ifndef WRITE_EACH_ITER
    write_conc(conc, ind_O3, "O3", iter-1, 0);
    #endif
    
    /* Show final time */
    printf("Final time: %ld seconds.\n", (iter-1)*dt);
    
    timer_stop(TIMER_WALLCLOCK);
    
    print_timer_summary("===Timers===");
    
    /* Cleanup and exit */
    return 0;
}
