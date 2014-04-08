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

void space_advec_diff_row(int n, volatile double *c, volatile double *w, volatile double *d, 
                          double *cb, double *wb, double *db, double cell_size, 
                          volatile double *dcdx);

void discretize(const int n, volatile double *conc_in, volatile double *wind, 
                volatile double *diff, double *concbound, double *windbound, 
                double *diffbound, double cell_size, double dt, 
                volatile double *conc_out);

#endif

