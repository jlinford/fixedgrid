/*
 *  spe_pthread.c
 *  
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include "spe_pthread.h"

/**
 * Thread entry function for spe threads
 */
void *spe_pthread_function(void *arg) 
{
    spe_pthread_data_t *data = (spe_pthread_data_t*)arg;
    unsigned int entry = SPE_DEFAULT_ENTRY;
    int rc;
    
    if ((rc = spe_context_run(data->speid, &entry, 0, (void*)(&data->argv), (void*)(&data->env), NULL)) < 0) 
    {
        fprintf(stderr, "Failed spe_context_run(rc=%d, errno=%d)\n", rc, errno);
        exit(1);
    }
    pthread_exit(NULL);
}

/**
 * Gets the status of an spe
 */
inline int spe_get_status(fixedgrid_t* G, int id)
{
    if(spe_out_mbox_status(G->threads[id].speid) > 0)
    {
        spe_out_mbox_read(G->threads[id].speid, &(G->threads[id].status), 1);
    }
    return G->threads[id].status;
}

/**
 * Sets the status of an spe
 */
inline void spe_set_status(fixedgrid_t* G, int id, uint32_t status)
{
    G->threads[id].status = status;
    spe_in_mbox_write(G->threads[id].speid, (uint32_t*)(&G->threads[id].status), 1, SPE_MBOX_ANY_NONBLOCKING);
}

/**
 * Waits for one spe to have SPE_WAITING_STATUS status
 */
inline void wait_for_spe(fixedgrid_t* G, int id)
{
    while(spe_get_status(G, id) != SPE_STATUS_WAITING) ; // Intentional wait
}


/**
 * Waits for all spes to have SPE_WAITING_STATUS status
 */
inline void wait_all_spes(fixedgrid_t* G)
{
    int i, t;
    
    do
    {
        for(i=0, t=0; i<G->nprocs; ++i)
        {
            t += spe_get_status(G, i);
        }
    }
    while(t > 0);
}

/**
 * Waits for all spes to complete execution
 */
inline void join_all_spes(fixedgrid_t* G)
{
    int i;
    
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
inline void create_spe_pthreads(fixedgrid_t* G)
{
    uint32_t i;
    
    for(i=0; i<G->nprocs; i++) 
    {
        /* Configure environment */
        G->threads[i].env.speid = i;
        G->threads[i].env.nprocs = G->nprocs;
        G->threads[i].env.metptr = (uint64_t)(&G->threads[i].metrics);
        
        /* Create context */
        if((G->threads[i].speid = spe_context_create(0, NULL)) == NULL) 
        {
            fprintf(stderr, "Failed spe_context_create(errno=%d strerror=%s)\n", errno, strerror(errno));
            exit(1);
        }
        
        /* Load program into context */
        if(spe_program_load(G->threads[i].speid, &fixedgrid_spu)) 
        {
            fprintf(stderr, "Failed spe_program_load(errno=%d strerror=%s)\n", errno, strerror(errno));
            exit(1);
        }
        
        /* Create thread for each SPE context */
        if(pthread_create(&G->threads[i].pthread, NULL, &spe_pthread_function, &G->threads[i])) 
        {
            fprintf(stderr, "Failed pthread_create(errno=%d strerror=%s)\n", errno, strerror(errno));
            exit(1);
        }
        
        spe_set_status(G, i, SPE_STATUS_WAITING);
    }
}
