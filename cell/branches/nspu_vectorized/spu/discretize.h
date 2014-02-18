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

#include "params.h"

void discretize_row(const int n, 
                    volatile value_t *conc_in, 
                    volatile value_t *wind, 
                    volatile value_t *diff, 
                    value_t *concbound, 
                    value_t *windbound, 
                    value_t *diffbound, 
                    value_t cell_size, 
                    value_t dt, 
                    volatile value_t *conc_out);

void discretize_col(const int n, 
                    volatile vector value_t *conc_in, 
                    volatile vector value_t *wind, 
                    volatile vector value_t *diff, 
                    vector value_t *concbound, 
                    vector value_t *windbound, 
                    vector value_t *diffbound, 
                    vector value_t cell_size, 
                    vector value_t dt, 
                    volatile vector value_t *conc_out);

#endif

