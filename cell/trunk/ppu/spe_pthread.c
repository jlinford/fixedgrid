/*
 *  spe_pthread.c
 *  
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <stdlib.h>

#include "spe_pthread.h"

/**
 * Thread entry function for spe threads
 */
void *spe_pthread_function(void *arg) 
{
    spe_pthread_data_t *data = (spe_pthread_data_t*)arg;
    unsigned int entry = SPE_DEFAULT_ENTRY;
    int rc;
    
    if ((rc = spe_context_run(data->speid, &entry, 0, (void*)(&data->argv), (void*)(&data->envv), NULL)) < 0) 
    {
        fprintf(stderr, "Failed spe_context_run(rc=%d, errno=%d)\n", rc, errno);
        exit(1);
    }
    pthread_exit(NULL);
}

/**
 * Waits for all spes to complete execution
 */
void join_all_spes(fixedgrid_t* G)
{
    uint32_t i;
    
    for(i=0; i<G->nprocs; i++) 
    {
        spe_set_status(G, i, SPE_STATUS_STOPPED);
        if(pthread_join(G->threads[i].pthread, NULL)) 
        {
            perror("Failed pthread_join");
            exit(1);
        }
    }
}

/**
 * Create and start several threads on the SPEs
 * @param nprocs Number of threads to start
 */
void create_spe_pthreads(fixedgrid_t* G)
{
    uint32_t i;
    
    for(i=0; i<G->nprocs; i++) 
    {
        /* Configure environment */
        G->threads[i].envv.speid = i;
        G->threads[i].envv.nprocs = G->nprocs;
        G->threads[i].envv.ea_status = (uint32_t)(&G->threads[i].status);
        G->threads[i].envv.ea_metrics = (uint32_t)(&G->threads[i].metrics);
        G->threads[i].status.value = SPE_STATUS_INIT;
        
        /* Create context */
        if((G->threads[i].speid = spe_context_create(0, NULL)) == NULL) 
        {
            fprintf(stderr, "Failed spe_context_create (errno=%d)\n", errno);
            exit(1);
        }
        
        /* Load program into context */
        if(spe_program_load(G->threads[i].speid, &fixedgrid_spu)) 
        {
            fprintf(stderr, "Failed spe_program_load (errno=%d)\n", errno);
            exit(1);
        }
        
        /* Create thread for each SPE context */
        if(pthread_create(&G->threads[i].pthread, NULL, &spe_pthread_function, &G->threads[i])) 
        {
            fprintf(stderr, "Failed pthread_create (errno=%d)\n", errno);
            exit(1);
        }
    }
}
