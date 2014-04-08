#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "fixedgrid.h"
#include "params.h"
#include "timer.h"
#include "fileio.h"
#include "saprc99_Monitor.h"
#include "chemistry.h"
#include "transport.h"

void saprc99_Initialize(real_t C[NSPEC]);

fixedgrid_t G_GLOBAL;

/* KPP-generated SAPRC'99 mechanism data */
double * C;                 /* Concentration of all species */
double * VAR;               /* First variable species */
double * FIX;               /* First fixed species */
double RCONST[NREACT];      /* Rate constants (global) */
double ATOL[NVAR];          /* Absolute tolerance */
double RTOL[NVAR];          /* Relative tolerance */
double TIME;                /* Current integration time */
double DT;                  /* Integration step */
double SUN;                 /* Sunlight intensity between [0,1] */
double TEMP;                /* Temperature */
double STEPMIN;             /* Lower bound for integration step */

/**
 * Fills an array with a value.
 * @param n     Length of the array
 * @param array Array to initialize
 * @param val   Value to initialize with
 */
void array_init(fixedgrid_t* G, unsigned int n, real_t *array, real_t val)
{
    int i;
    
    timer_start(&G->metrics.array_init);
    
    for(i=0; i<n; i++)
    {
        array[i] = val;
    }
    
    timer_stop(&G->metrics.array_init);
}

/**
 * Processes emission sources
 */
void process_emissions(fixedgrid_t* G)
{
    /* Add O3 plume */
    G->conc(SOURCE_X, SOURCE_Y, SOURCE_Z, ind_O3) += SOURCE_RATE / (DX * DY * DZ);
}

/**
 * Initializes the model
 */
void init_model(fixedgrid_t* G)
{
    uint32_t i;
    uint32_t x, y, z, s;
    
    /* Chemistry buffer */
    real_t chemBuff[NSPEC];
    
    metrics_init(&G->metrics, "Serial");
        
    /* Initialize time frame */
    /* FIXME: year is ignored */
    G->tstart = day2sec(START_DOY) + hour2sec(START_HOUR) + minute2sec(START_MIN);
    G->tend   = day2sec(END_DOY)   + hour2sec(END_HOUR)   + minute2sec(END_MIN);
    G->dt = STEP_SIZE;
    G->time = G->tstart;
    
    /* Initialize chemistry and concentration data */
    printf("Loading chemistry and concentration data... ");
    timer_start(&G->metrics.array_init);
    
#if DO_CHEMISTRY == 1
    
    /* Set saprc'99 parameters */
    TIME = G->time;
    DT   = G->dt;
    SUN  = (real_t)0.0;
    TEMP = TEMP_INIT;   // FIXME: TEMP is now a field
    STEPMIN = 0.01;
    
    /* Initialize tolerances */
    for( i = 0; i < NVAR; i++ ) {
        RTOL[i] = 1.0e-3;
        ATOL[i] = 1.0;
    }
    
    /* Initialize concentrations */
    saprc99_Initialize(chemBuff);
    
    for(s=0; s<NSPEC; s++)
    {
        for(z=0; z<NZ; z++)
        {
            for(y=0; y<NY; y++)
            {
                for(x=0; x<NX; x++)
                {
                    G->conc(x, y, z, s) = chemBuff[s];
                }
            }
        }
    }
    
#else
    
    for(z=0; z<NZ; z++)
    {
        for(y=0; y<NY; y++)
        {
            for(x=0; x<NX; x++)
            {
                G->conc(x, y, z, ind_O3) = O3_INIT;
            }
        }
    }
    
#endif
    
    timer_stop(&G->metrics.array_init);
    printf("done.\n");
    
    /* Initialize wind field */
    printf("Loading wind field data...");
    array_init(G, NX*NY*NZ, &G->wind_u(0,0,0), WIND_U_INIT);
    array_init(G, NX*NY*NZ, &G->wind_v(0,0,0), WIND_V_INIT);
    array_init(G, NX*NY*NZ, &G->wind_w(0,0,0), WIND_W_INIT);
    printf(" done.\n");
    
    /* Initialize diffusion field */
    printf("Loading diffusion field data...");
    array_init(G, NX*NY*NZ, &G->diff_h(0,0,0), DIFF_H_INIT);
    array_init(G, NX*NY*NZ, &G->diff_v(0,0,0), DIFF_V_INIT);
    printf(" done.\n");
    
    /* Initialize diffusion field */
    printf("Loading temperature field data...");
    array_init(G, NX*NY*NZ, &G->temp(0,0,0), TEMP_INIT);
    printf(" done.\n");
}


/**
 * Displays emission source locations and rates
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
    uint32_t i;
    uint32_t steps;
    
    steps = (G->tend - G->tstart) / G->dt;
    
    printf("\n");
    printf("CONFIGURATION:\n");
    printf("    X DISCRETIZATION:   %s\n", DO_X_DISCRET == TRUE ? "TRUE" : "FALSE");
    printf("    Y DISCRETIZATION:   %s\n", DO_Y_DISCRET == TRUE ? "TRUE" : "FALSE");
    printf("    Z DISCRETIZATION:   %s\n", DO_Z_DISCRET == TRUE ? "TRUE" : "FALSE");
    printf("    SAPRC99 CHEMISTRY:  %s\n", DO_CHEMISTRY == TRUE ? "TRUE" : "FALSE");
    printf("    DOUBLE PRECISION:   %s\n", DOUBLE_PRECISION == TRUE ? "TRUE" : "FALSE");
    printf("\n");
    printf("SPACE DOMAIN:\n");
    printf("    LENGTH (X): %f meters\n", NX*DX);
    printf("    WIDTH  (Y): %f meters\n", NY*DY);
    printf("    DEPTH  (Z): %f meters\n", NZ*DZ);
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
 * Program entry point.
 */
int main(int argc, char** argv)
{
    /* Global data pointer */
    fixedgrid_t* G = &G_GLOBAL;
    
    /* Iterators */
    int i, iter;
    
    /* Start wall clock timer */
    timer_start(&G->metrics.wallclock);

    G->nprocs = omp_get_max_threads();

    /* Parse command line arguments */
    if(argc > 1)
    {
        i = atoi(argv[1]);
        if(i < 1)
        {
            fprintf(stderr, "Invalid number of threads: %d < 1.\n", i);
            exit(1);
        }

        if(i <= G->nprocs)
        {
            G->nprocs = i;
        }
        else
        {
            printf("%d threads unavailable.  Using %d instead.\n", i, G->nprocs);
        }
    }

    printf("\nRunning on %d threads\n", G->nprocs);

    omp_set_num_threads(G->nprocs);
    
    /* Initialize the model parameters */
    init_model(G);
    
    /* Add emissions */
    process_emissions(G);
    
    /* Print startup banner */
    print_start_banner(G);
    
    /* Store initial concentration */
    printf("Writing initial concentration...");
    write_conc(G, 0, 0);
    printf(" done.\n");
    
    /* BEGIN CALCULATIONS */
    for(iter=1, G->time = G->tstart; G->time < G->tend; G->time += G->dt, ++iter)
    {
        /* Chemistry */
        saprc99_chem(G);
        
	#pragma omp parallel shared(G, iter)
	{
		discretize_all_x(G, G->dt*0.5);

		discretize_all_y(G, G->dt*0.5);
		
		discretize_all_z(G, G->dt);
		
		discretize_all_y(G, G->dt*05);
		
		discretize_all_x(G, G->dt*05);
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
        
        /* Store concentration */
        #if WRITE_EACH_ITER == 1
        write_conc(G, iter, 0);
        #endif
        
        /* Indicate progress */
        printf("  After iteration %02d: Model time = %07.2f sec.\n", iter, iter*G->dt);
    }
    /* END CALCULATIONS */
    
    /* Store concentration */
    #if WRITE_EACH_ITER != 1
    write_conc(G, iter-1, 0);
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
