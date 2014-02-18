#include <stdio.h>
#include <stdlib.h>

#include "fixedgrid.h"
#include "params.h"
#include "timer.h"
#include "discretize.h"
#include "fileio.h"

#ifdef DO_CHEMISTRY
#include "saprc99.h"
#endif

fixedgrid_t G_GLOBAL;

/**
 * Fills an array with a value.
 * @param n     Length of the array
 * @param array Array to initialize
 * @param val   Value to initialize with
 */
void double_array_init(fixedgrid_t* G, unsigned int n, double *array, double val)
{
    int i;
    
    timer_start(&G->metrics.array_copy);
    
    for(i=0; i<n; i++)
    {
        array[i] = val;
    }
    
    timer_stop(&G->metrics.array_copy);
}

/**
 * Copies one array to another
 * @param n     Length of the array
 * @param src   Source array
 * @param dst   Destination array
 */
void double_array_copy(fixedgrid_t* G, unsigned int n, double *src, double *dst)
{
    int i;
    
    timer_start(&G->metrics.array_copy);
    
    for(i=0; i<n; i++)
    {
        dst[i] = src[i];
    }
    
    timer_stop(&G->metrics.array_copy);
}

/**
 * Applies saprc99 chemical mechanism to all chemical species
 */
void saprc99_chem(fixedgrid_t* G)
{
#ifdef DO_CHEMISTRY
    int i, j, s;
    
    double buff[NSPEC];
    
    timer_start(&G->metrics.chem);
    for(i=0; i<NROWS; i++)
    {
        for(j=0; j<NCOLS; j++)
        {
            for(s=0; s<NSPEC; s++)
            {
                buff[s] = G->conc[s][i][j];
            }
            //printf("i: %d, j: %d, temp: %f, tin: %f, tout: %f, dt: %d\n", i, j, temp, time, time+dt, dt);
            saprc99(buff, temp, time, time+dt, dt);
            for(s=0; s<NSPEC; s++)
            {
                G->conc[s][i][j] = buff[s];
            }
        }
    }
    timer_stop(&G->metrics.chem);
#endif
}

/**
 * Discretize rows 1/2 timestep 
 * @param s     Species index
 */
void discretize_all_rows(fixedgrid_t* G, double dt)
{
#ifdef DO_ROW_DISCRET
    
    int i, j, s;
    
    /* Buffers */
    /*
    double crow1[NCOLS];
    double crow2[NCOLS];
    double wrow[NCOLS];
    double drow[NCOLS];
     */
    double buff[NCOLS];
    
    /* Boundary values */
    double cbound[4];
    double wbound[4];
    double dbound[4];    
    
    timer_start(&G->metrics.row_discret);
    
    for(i=0; i<NROWS; i++)
    {
        #ifdef DO_CHEMISTRY
        for(s=0; s<NSPEC; s++)
        {
        #else
        s = ind_O3;
        #endif

            /*
            timer_start(&G->metrics.array_copy);
            for(k=0; k<NCOLS; k++)
            {
                crow1[k] = G->conc[i][k][s];
            }
            timer_stop(&G->metrics.array_copy);
            
            double_array_copy(G, NCOLS, G->wind_u[i], wrow);
            double_array_copy(G, NCOLS, G->diff[i], drow);
             */
            
            timer_start(&G->metrics.array_copy);
            cbound[0] = G->conc[s][i][NCOLS-2];
            cbound[1] = G->conc[s][i][NCOLS-1];
            cbound[2] = G->conc[s][i][0];
            cbound[3] = G->conc[s][i][1];
            wbound[0] = G->wind_u[i][NCOLS-2];
            wbound[1] = G->wind_u[i][NCOLS-1];
            wbound[2] = G->wind_u[i][0];
            wbound[3] = G->wind_u[i][1];
            dbound[0] = G->diff[i][NCOLS-2];
            dbound[1] = G->diff[i][NCOLS-1];
            dbound[2] = G->diff[i][0];
            dbound[3] = G->diff[i][1];
            timer_stop(&G->metrics.array_copy);
            
            discretize(NCOLS, G->conc[s][i], G->wind_u[i], G->diff[i], cbound, wbound, dbound, DX, dt, buff);
            
            timer_start(&G->metrics.array_copy);
            for(j=0; j<NCOLS; j++)
                G->conc[s][i][j] = buff[j];
            timer_stop(&G->metrics.array_copy);
            
        #ifdef DO_CHEMISTRY
        }
        #endif
        
    }
    
    timer_stop(&G->metrics.row_discret);
    
#endif
}

/**
 * Discretize colums 1 timestep 
 * @param s     Species index
 */
void discretize_all_cols(fixedgrid_t* G, double dt)
{
#ifdef DO_COL_DISCRET
    
    int j, k, s;
    
    /* Buffers */
    double ccol1[NROWS];
    double ccol2[NROWS];
    double wcol[NROWS];
    double dcol[NROWS];
    
    /* Boundary values */
    double cbound[4];
    double wbound[4];
    double dbound[4];    
    
    timer_start(&G->metrics.col_discret);
    
    for(j=0; j<NCOLS; j++)
    {
        #ifdef DO_CHEMISTRY
        for(s=0; s<NSPEC; s++)
        {
        #else
        s = ind_O3;
        #endif
            
            timer_start(&G->metrics.array_copy);
            for(k=0; k<NROWS; k++)
            {
                ccol1[k] = G->conc[s][k][j];
                wcol[k]  = G->wind_v[k][j];
                dcol[k]  = G->diff[k][j];
            }
            
            cbound[0] = ccol1[NROWS-2];
            cbound[1] = ccol1[NROWS-1];
            cbound[2] = ccol1[0];
            cbound[3] = ccol1[1];
            wbound[0] = wcol[NROWS-2];
            wbound[1] = wcol[NROWS-1];
            wbound[2] = wcol[0];
            wbound[3] = wcol[1];
            dbound[0] = dcol[NROWS-2];
            dbound[1] = dcol[NROWS-1];
            dbound[2] = dcol[0];
            dbound[3] = dcol[1];
            timer_stop(&G->metrics.array_copy);
            
            discretize(NROWS, ccol1, wcol, dcol, cbound, wbound, dbound, DY, dt, ccol2);
            
            timer_start(&G->metrics.array_copy);
            for(k=0; k<NROWS; k++)
                G->conc[s][k][j] = ccol2[k];
            timer_stop(&G->metrics.array_copy);
            
        #ifdef DO_CHEMISTRY
        }
        #endif
        
    }
    
    timer_stop(&G->metrics.col_discret);
    
#endif
}


/**
 * Displays emission source locations and rates
 */
void print_emission_sources()
{
    printf("\n");
    printf("Emission sources (X, Y, Z, RATE):\n");
    printf("    (%f, %f, %f, %E)\n", (DX*SOURCE_X + DX*0.5), (DY*SOURCE_Y + DY*0.5), 0.0, SOURCE_RATE);
    printf("\n");
}

/**
 * Prints program startup banner
 */
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
    printf("    FROM  %02d:%02d.00 on day %03d of year %d\n", START_HOUR, START_MIN, START_DOY, START_YEAR);
    printf("    TO    %02d:%02d.00 on day %03d of year %d\n", END_HOUR, END_MIN, END_DOY, END_YEAR);
    printf("    TOTAL %ld seconds (%ld timesteps of %d seconds)\n", t_end, steps, STEP_SIZE);
    
    print_emission_sources();
    
    printf("\n");
}

/**
 * Program entry point.
 */
int main(int argc, char** argv)
{
    /* Global data pointer */
    fixedgrid_t* G = &G_GLOBAL;
    
    /* Iterators */
    int i, j, s, iter;
    
    double buff[NSPEC];

    metrics_init(&G->metrics, "Serial");
        
    /* Start wall clock timer */
    timer_start(&G->metrics.wallclock);

    G->nprocs = 1;
        
    printf("\nRunning %d thread.\n", G->nprocs);
    
    /* Initialize chemistry mechanism */
    #ifdef DO_CHEMISTRY
    init_saprc99(G);
    #endif
    
    /* Initialize concentration data */
    printf("Loading chemistry and concentration data... ");
    for(i=0; i<NROWS; i++)
    {
        for(j=0; j<NCOLS; j++)
        {
            #ifdef DO_CHEMISTRY
            init_conc(buff);
            #else
            double_array_init(G, NSPEC, buff, O3_INIT);
            #endif
            for(s=0; s<NSPEC; s++)
            {
                G->conc[s][i][j] = buff[s];
            }            
        }
    }
    printf("done.\n");

    /* Initialize wind field */
    printf("Loading wind field data...");
    double_array_init(G, NROWS*NCOLS, &(G->wind_u[0][0]), WIND_U_INIT);
    double_array_init(G, NROWS*NCOLS, &(G->wind_v[0][0]), WIND_V_INIT);
    printf(" done.\n");
    
    /* Initialize diffusion field */
    printf("Loading diffusion field data...");
    double_array_init(G, NROWS*NCOLS, &(G->diff[0][0]), DIFF_INIT);
    printf(" done.\n");
    
    /* Initialize time */
    G->end_time = year2sec(END_YEAR - START_YEAR) + day2sec(END_DOY - START_DOY) + hour2sec(END_HOUR - START_HOUR) + minute2sec(END_MIN - START_MIN);
    G->dt = STEP_SIZE;
    G->steps = G->end_time/G->dt;
    
    /* Initialize environment */
    G->temp = TEMP_INIT;

    /* Add O3 plume emissions */
    G->conc[ind_O3][SOURCE_Y][SOURCE_X] += (SOURCE_RATE) / (DX * DY * 1000.0);
    
    /* Print startup banner */
    print_start_banner(NCOLS*DX, NROWS*DY, 0.0, G->end_time, G->steps);
    
    /* Store initial concentration */
    printf("Writing initial concentration...");
    write_conc(G, ind_O3, "O3", 0, 0);
    printf(" done.\n");
    
    /* BEGIN CALCULATIONS */
    for(iter = 1; iter <= G->steps; iter++)
    {
        /* Chemistry */
        saprc99_chem(G);
        
        discretize_all_rows(G, G->dt/2.0);

        discretize_all_cols(G, G->dt);
        
        discretize_all_rows(G, G->dt/2.0);
        
        /*
         * Could update wind field here...
         */
         
        /*
         * Could update diffusion tensor here...
         */
         
        /*
         * Could update environment here...
         */
        
        /* Store concentration */
        #ifdef WRITE_EACH_ITER
        write_conc(G, ind_O3, "O3", iter, 0);
        #endif
        
        /* Indicate progress */
        printf("Iteration %d of %d.  Time = %f seconds.\n", iter, G->steps, iter*G->dt);
    }
    /* END CALCULATIONS */
    
    /* Store concentration */
    #ifndef WRITE_EACH_ITER
    write_conc(G, ind_O3, "O3", iter-1, 0);
    #endif
    
    /* Show final time */
    printf("Final time: %f seconds.\n", (iter-1)*G->dt);
    
    timer_stop(&G->metrics.wallclock);
    
    /* Print metrics */
    print_metrics(&G->metrics);
    
    /* Write metrics to CSV file */
    write_metrics_as_csv(G, "Serial");
    
    /* Cleanup and exit */
    return 0;
}
