/*
 *  spe_pthread.h
 *  
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __SPE_PTHREAD_H__
#define __SPE_PTHREAD_H__

/**************************************************
 * Includes                                       *
 **************************************************/

#include <stdint.h>

#include "fixedgrid.h"

/**************************************************
 * External Data                                  *
 **************************************************/

/* SPE program handle */
extern spe_program_handle_t fixedgrid_spu;

/**************************************************
 * Function Prototypes                            *
 **************************************************/

/**
 * Thread entry function for spe threads
 */
void *spe_pthread_function(void *arg);

/**
 * Gets the status of an spe
 */
int spe_get_status(fixedgrid_t* G, int id);

/**
 * Sets the status of an spe
 */
void spe_set_status(fixedgrid_t* G, int id, uint32_t status);

/**
 * Waits for one spe to have SPE_WAITING_STATUS status
 */
void wait_for_spe(fixedgrid_t* G, int id);

/**
 * Waits for all spes to have SPE_WAITING_STATUS status
 */
void wait_all_spes(fixedgrid_t* G);

/**
 * Waits for all spes to complete execution
 */
void join_all_spes(fixedgrid_t* G);

/**
 * Create and start several threads on the SPEs
 */
void create_spe_pthreads(fixedgrid_t* G);


#endif
