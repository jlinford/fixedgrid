#include <stdio.h>
#include <stdlib.h>
#include <libspe2.h>
#include <pthread.h>
#include <errno.h>
#include <string.h>
#include <malloc_align.h>

#include "fixedgrid.h"
#include "params.h"
#include "spu/timer.h"
#include "fileio.h"

/* Thread data datatype */
typedef struct ppu_pthread_data 
{
    uint32_t status;
    spe_context_ptr_t speid;
    pthread_t pthread;
    void *argp;
} ppu_pthread_data_t;

/* SPE program handle */
extern spe_program_handle_t fixedgrid_spu;

/* 2D concentration data */
volatile double conc[NROWS*NCOLS] __attribute__((aligned(128)));

/* 2D wind field data */
volatile double wind_u[NROWS*NCOLS] __attribute__((aligned(128)));
volatile double wind_v[NROWS*NCOLS] __attribute__((aligned(128)));

/* 2D diffusion tensor data */
volatile double diff[NROWS*NCOLS] __attribute__((aligned(128)));

/* Column buffers */
volatile double ccol[MAX_THREADS*NY] __attribute__((aligned(128)));
volatile double wcol[MAX_THREADS*NY] __attribute__((aligned(128)));
volatile double dcol[MAX_THREADS*NY] __attribute__((aligned(128)));

/* SPE argument vectors */
volatile fixedgrid_spe_argv_t spe_argvs[MAX_THREADS] __attribute__((aligned(128)));

/* SPE threads */
int nprocs;
ppu_pthread_data_t threads[MAX_THREADS];


void double_array_init(unsigned int n, volatile double *array, double val)
{
    unsigned int i;
    
    timer_start(TIMER_ARRAY_INIT);
    
    for(i=0; i<n; i++)
    {
        *(array + i) = val;
    }
    
    timer_stop(TIMER_ARRAY_INIT);
}

void double_array_copy(unsigned int n, volatile double *src, double *dst)
{
    unsigned int i;
    
    timer_start(TIMER_ARRAY_COPY);
    
    for(i=0; i<n; i++)
    {
        *(dst + i) = *(src + i);
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

void *ppu_pthread_function(void *arg) 
{
    ppu_pthread_data_t *data = (ppu_pthread_data_t*)arg;
    unsigned int entry = SPE_DEFAULT_ENTRY;
    int rc;
    
    if ((rc = spe_context_run(data->speid, &entry, 0, data->argp, NULL, NULL)) < 0) 
    {
        fprintf(stderr, "Failed spe_context_run(rc=%d, errno=%d)\n", rc, errno);
        exit(1);
    }
    pthread_exit(NULL);
}

int spe_get_status(int id)
{
    if(spe_out_mbox_status(threads[id].speid) > 0)
    {
        spe_out_mbox_read(threads[id].speid, &(threads[id].status), 1);
    }
    return threads[id].status;
}

void spe_set_status(int id, uint32_t status)
{
    threads[id].status = status;
    spe_in_mbox_write(threads[id].speid, (uint32_t*)(&threads[id].status), 1, SPE_MBOX_ANY_NONBLOCKING);
}

void wait_all_spes()
{
    int i, t;
    
    do
    {
        for(i=0, t=0; i<nprocs; i++)
        {
            t += spe_get_status(i);
        }
    }
    while(t > 0);
}

int main(int argc, char** argv)
{
    /* Iterators */
    int i, j, k;
    
    uint32_t block;
    
    /* Time (seconds) */
    long t_0;
    long t_end;
    long dt;
    long steps;
    long iter;
    
    /* Emission control */
    bool emflag = TRUE;
    
    /* Start wall clock timer */
    timer_start(TIMER_WALLCLOCK);
    
    /* Initialize parallelization */
    nprocs = spe_cpu_info_get(SPE_COUNT_USABLE_SPES, -1);
    nprocs = nprocs > MAX_THREADS ? MAX_THREADS : nprocs;
    
    if(argc > 1)
    {
        i = atoi(argv[1]);
        if(i < 1)
        {
            fprintf(stderr, "Invalid number of SPUs: %d < 1.\n", i);
            exit(1);
        }
        
        if(i < nprocs)
        {
            nprocs = i;
        }
        else 
        {
            printf("%d SPUs unavailable.  Using %d instead.\n", i, nprocs);
        }
    }
    
    /* Create SPE threads */
    for(i=0; i<nprocs; i++) 
    {
        threads[i].argp = (void*)(&spe_argvs[i]);
        
        /* Create context */
        if((threads[i].speid = spe_context_create(0, NULL)) == NULL) 
        {
            fprintf(stderr, "Failed spe_context_create(errno=%d strerror=%s)\n", errno, strerror(errno));
            exit(1);
        }
        
        /* Load program into context */
        if(spe_program_load(threads[i].speid, &fixedgrid_spu)) 
        {
            fprintf(stderr, "Failed spe_program_load(errno=%d strerror=%s)\n", errno, strerror(errno));
            exit(1);
        }
            
        /* Create thread for each SPE context */
        if(pthread_create(&threads[i].pthread, NULL, &ppu_pthread_function, &threads[i])) 
        {
            fprintf(stderr, "Failed pthread_create(errno=%d strerror=%s)\n", errno, strerror(errno));
            exit(1);
        }
        
        spe_set_status(i, SPE_STATUS_WAITING);
    }
    
    printf("\nRunning %d threads (%d SPU + 1 PPU).\n", (nprocs+1), nprocs);
    
    /* Allocate concentration memory */
    //conc = _malloc_align(NROWS*NCOLS*sizeof(double), 7);
    //conc_buff = (double*)_malloc_align(MAX_THREADS*NY*sizeof(double), 7);

    /* Allocation wind vector filed memory */
    //wind_u = _malloc_align(NROWS*NCOLS*sizeof(double), 7);
    //wind_v = _malloc_align(NROWS*NCOLS*sizeof(double), 7);
    //wind_u_buff = (double*)_malloc_align(MAX_THREADS*NY*sizeof(double), 7);
    //wind_v_buff = (double*)_malloc_align(MAX_THREADS*NY*sizeof(double), 7);
 
    /* Allocation diffusion tensor memory */
    //diff = _malloc_align(NROWS*NCOLS*sizeof(double), 7);
    //diff_buff = (double*)_malloc_align(MAX_THREADS*NY*sizeof(double), 7);

    /* Initialize concentration data */
    double_array_init(NROWS*NCOLS, conc, O3_INIT);
        
    /* Initialize wind field */
    double_array_init(NROWS*NCOLS, wind_u, WIND_U_INIT);
    double_array_init(NROWS*NCOLS, wind_v, WIND_V_INIT);
    
    /* Initialize diffusion field */
    double_array_init(NROWS*NCOLS, diff, DIFF_INIT);
    
    /* Initialize time */
    t_0 = 0.0;
    t_end = year2sec(END_YEAR - START_YEAR) + day2sec(END_DOY - START_DOY) + 
            hour2sec(END_HOUR - START_HOUR) + minute2sec(END_MIN - START_MIN);
    dt = STEP_SIZE;
    steps = (long)( (t_end - t_0)/dt );
    
    /* Print startup banner */
    print_start_banner(NX*DX, NY*DY, 0.0, t_end, steps);
    
    /* Store initial concentration */
    write_conc(&(conc[0]), 0, 0);
    
    /* BEGIN CALCULATIONS */
    for(iter = 1; iter <= steps; iter++)
    {
        emflag = iter*dt < 6*3600.0 ? TRUE : FALSE;
        
        timer_start(TIMER_ROW_DISCRET);
        
        /* Discretize rows 1/2 timestep */
        block = NROWS / nprocs;
        for(i=0; i<nprocs; i++)
        {
            /* Configure SPE arguments */
            spe_argvs[i].arg[0].u64 = (uint64_t)(&conc[i*block*NX]);
            spe_argvs[i].arg[1].u64 = (uint64_t)(&wind_u[i*block*NX]);
            spe_argvs[i].arg[2].u64 = (uint64_t)(&diff[i*block*NX]);
            spe_argvs[i].arg[3].dbl = dt/2;
            spe_argvs[i].arg[4].dbl = DX;
            spe_argvs[i].arg[5].u32[0] = NX;
            spe_argvs[i].arg[5].u32[1] = (i == nprocs - 1 ? block + NROWS % nprocs : block);  //FIXME
            
            /* Signal SPE */
            spe_set_status(i, SPE_STATUS_WORKING);
        }
        
        /* Wait for SPEs to finish */
        wait_all_spes();
        
        timer_stop(TIMER_ROW_DISCRET);
        
        timer_start(TIMER_COL_DISCRET);
        
        /* Discretize colums 1 timestep */
        for(i=0; i<NCOLS; i++)
        {
            k = i % nprocs;

            while(spe_get_status(k) > 0) ; //intentional wait
            
            if(i >= nprocs)
            {
                timer_start(TIMER_ARRAY_COPY);
                for(j=0; j<NY; j++)
                {
                    conc[i-nprocs + j*NX] = ccol[k*NY+j];
                }
                timer_stop(TIMER_ARRAY_COPY);
            }
            
            timer_start(TIMER_ARRAY_COPY);
            for(j=0; j<NY; j++)
            {
                ccol[k*NY + j] = conc[i + j*NX];
                wcol[k*NY + j] = wind_v[i + j*NX];
                dcol[k*NY + j] = diff[i + j*NX];
            }
            timer_stop(TIMER_ARRAY_COPY);

            // Configure SPE arguments 
            spe_argvs[k].arg[0].u64 = (uint64_t)(&ccol[k*NY]);
            spe_argvs[k].arg[1].u64 = (uint64_t)(&wcol[k*NY]);
            spe_argvs[k].arg[2].u64 = (uint64_t)(&dcol[k*NY]);
            spe_argvs[k].arg[3].dbl = dt;
            spe_argvs[k].arg[4].dbl = DY;
            spe_argvs[k].arg[5].u32[0] = NY;
            spe_argvs[k].arg[5].u32[1] = 1;

            // Signal SPE 
            spe_set_status(k, SPE_STATUS_WORKING);
        }

        /* Wait for SPEs to finish */
        wait_all_spes();
        
        timer_stop(TIMER_COL_DISCRET);
        
        timer_start(TIMER_ROW_DISCRET);
        
        /* Discretize rows 1/2 timestep */
        block = NROWS / nprocs;
        for(i=0; i<nprocs; i++)
        {
            /* Configure SPE arguments */
            spe_argvs[i].arg[0].u64 = (uint64_t)(&conc[i*block*NX]);
            spe_argvs[i].arg[1].u64 = (uint64_t)(&wind_u[i*block*NX]);
            spe_argvs[i].arg[2].u64 = (uint64_t)(&diff[i*block*NX]);
            spe_argvs[i].arg[3].dbl = dt/2;
            spe_argvs[i].arg[4].dbl = DX;
            spe_argvs[i].arg[5].u32[0] = NX;
            spe_argvs[i].arg[5].u32[1] = (i == nprocs - 1 ? block + NROWS % nprocs : block);  //FIXME
            
            /* Signal SPE */
            spe_set_status(i, SPE_STATUS_WORKING);
        }
        
        /* Wait for SPEs to finish */
        wait_all_spes();
        
        timer_stop(TIMER_ROW_DISCRET);
        
        /*
         * Could update wind field here...
         */
         
        /*
         * Could update diffusion tensor here...
         */
        
        /* Add emissions */
        if(emflag)
        {
            conc[SOURCE_Y*NX + SOURCE_X] += dt * (SOURCE_RATE) / (DX * DY * 1000.0);
        }
        
        /* Store concentration */
        #ifdef WRITE_EACH_ITER
        write_conc(conc, iter, 0);
        #endif
        
        /* Indicate progress */
        if(iter % 10 == 0)
        {
            printf("Iteration %ld of %ld.  Time = %ld seconds.\n", iter, steps, iter*dt);
        }
        
    }
    /* END CALCULATIONS */
    
    /* Wait for SPU-thread to complete execution. */
    for(i=0; i<nprocs; i++) 
    {
        spe_set_status(i, SPE_STATUS_STOPPED);
        if(pthread_join(threads[i].pthread, NULL)) 
        {
            perror("Failed pthread_join");
            exit(1);
        }
    }
    
    /* Store concentration */
    write_conc(conc, iter-1, 0);
    
    /* Show final time */
    printf("Final time: %ld seconds.\n", (iter-1)*dt);
    
    timer_stop(TIMER_WALLCLOCK);

    print_timer_summary("===PPU Timers===");    
    
    /* Cleanup and exit */
    return 0;
}

