#ifndef __DISCRETIZE_H__
#define __DISCRETIZE_H__

void space_advec_diff(const int n, const double *c, const double *w, const double *d, 
                      const double cell_size, double *dcdx);

void discretize(const int n, const double *conc_in, const double *wind, 
                const double *diff, const double cell_size, const double dt, 
                double *conc_out);

#endif

