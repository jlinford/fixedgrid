#ifndef __DISCRETIZE_H__
#define __DISCRETIZE_H__

void space_advec_diff(int n, volatile double *c, volatile double *w, volatile double *d, 
                      double cell_size, volatile double *dcdx);

void discretize(int n, volatile double *conc_in, volatile double *wind, 
                volatile double *diff, double cell_size, double dt, 
                volatile double *conc_out);

#endif

