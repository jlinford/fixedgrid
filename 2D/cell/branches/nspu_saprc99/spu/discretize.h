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

#include <stdint.h>
#include "params.h"

void discretize_row(const uint32_t n, 
                    volatile real_t *conc_in, 
                    volatile real_t *wind, 
                    volatile real_t *diff, 
                    real_t *concbound, 
                    real_t *windbound, 
                    real_t *diffbound, 
                    real_t cell_size, 
                    real_t dt, 
                    volatile real_t *conc_out);

void discretize_col(const uint32_t n, 
                    volatile vector real_t *conc_in, 
                    volatile vector real_t *wind, 
                    volatile vector real_t *diff, 
                    vector real_t *concbound, 
                    vector real_t *windbound, 
                    vector real_t *diffbound, 
                    vector real_t cell_size, 
                    vector real_t dt, 
                    volatile vector real_t *conc_out);

#endif

