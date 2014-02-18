/*
 *  discretize.h
 *  
 *  Transport module.
 *  Main kernel.
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __DISCRETIZE_H__
#define __DISCRETIZE_H__

/**************************************************
 * Includes                                       *
 **************************************************/

#include <stdint.h>
#include "params.h"

/**************************************************
 * Function prototypes                            *
 **************************************************/

void space_advec_diff_row(int n,  real_t *c,  real_t *w,  real_t *d, 
                          real_t *cb, real_t *wb, real_t *db, real_t cell_size, 
                           real_t *dcdx);

void discretize(const int n,  real_t *conc_in,  real_t *wind, 
                 real_t *diff, real_t *concbound, real_t *windbound, 
                real_t *diffbound, real_t cell_size, real_t dt, 
                 real_t *conc_out);


#endif

