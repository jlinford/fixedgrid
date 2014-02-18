/*
 *  discretize.c
 *  
 *  Transport module.
 *  Main kernel.
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include "discretize.h"
#include "timer.h"
#include "common.h"

extern metrics_t metrics;

/* 
 * Macro-define the upwinded advection/diffusion equation.
 * c = conc, w = wind, d = diff
 * []2l is the 2-left neighbor of [], etc.
 * []2r is the 2-right neighbor of [], etc.
 * out = output variable
 *
 * Origional code:
 *
 * wind = (w[i-1] + w[i]) / 2.0;
 * if(wind > 0.0)
 * {
 *     advec_termL = (1.0/6.0) * ( (-c[i-2]) + 5.0*c[i-1] + 2.0*c[i] );
 * }
 * else
 * {
 *     advec_termL = (1.0/6.0) * ( 2.0*c[i-1] + 5.0*c[i] + (-c[i+1]) );
 * }
 * advec_termL *= wind;
 * 
 * wind = (w[i+1] + w[i]) / 2.0;
 * if(wind > 0.0)
 * {
 *     advec_termR = (1.0/6.0) * ( (-c[i-1]) + 5.0*c[i] + 2.0*c[i+1] );
 * }
 * else
 * {
 *     advec_termR = (1.0/6.0) * ( 2.0*c[i] + 5.0*c[i+1] + (-c[i+2]) );
 * }
 * advec_termR *= wind;
 * 
 * advec_term = (advec_termL - advec_termR) / cell_size;
 * 
 * diff_term = ( ((d[i-1]+d[i])/2)*(c[i-1]-c[i]) - ((d[i]+d[i+1])/2)*(c[i]-c[i+1]) ) / (cell_size * cell_size);
 * 
 * dcdx[i-shift] = advec_term + diff_term;
 *
 */
#define ADVEC_DIFF(c2l, w2l, d2l, c1l, w1l, d1l, c, w, d, c1r, w1r, d1r, c2r, w2r, d2r, out) {  \
    double wind, diff_term, advec_term, advec_termL, advec_termR;                               \
    wind = (w1l + w) / 2.0;                                                                     \
    if(wind >= 0.0) advec_termL = (1.0/6.0) * ( -c2l + 5.0*c1l + 2.0*c );                       \
    else advec_termL = (1.0/6.0) * ( 2.0*c1l + 5.0*c - c1r );                                   \
    advec_termL *= wind;                                                                        \
    wind = (w1r + w) / 2.0;                                                                     \
    if(wind >= 0.0) advec_termR = (1.0/6.0) * ( -c1l + 5.0*c + 2.0*c1r );                       \
    else advec_termR = (1.0/6.0) * ( 2.0*c + 5.0*c1r - c2r );                                   \
    advec_termR *= wind;                                                                        \
    advec_term = (advec_termL - advec_termR) / cell_size;                                       \
    diff_term = ( ((d1l+d)/2)*(c1l-c) - ((d+d1r)/2)*(c-c1r) ) / (cell_size * cell_size);        \
    out = advec_term + diff_term;                                                               \
}

/*
 * Applies the advection / diffusion equation to a row of data
 */
void space_advec_diff_row(int n, volatile double *c, volatile double *w, volatile double *d, 
                          double *cb, double *wb, double *db, double cell_size, 
                          volatile double *dcdx)
{
    int i;
    
    /* Do boundary cell c[0] explicitly */
    ADVEC_DIFF(cb[0], wb[0], db[0],  /* 2-left neighbors */
               cb[1], wb[1], db[1],  /* 1-left neighbors */
               c[0], w[0], d[0],     /* Values */
               c[1], w[1], d[1],     /* 1-right neighbors */
               c[2], w[2], d[2],     /* 2-right neighbors */
               dcdx[0]);             /* Output */
    
    /* Do boundary cell c[1] explicitly */
    ADVEC_DIFF(cb[1], wb[1], db[1],  /* 2-left neighbors */
               cb[2], wb[2], db[2],  /* 1-left neighbors */
               c[1], w[1], d[1],     /* Values */
               c[2], w[2], d[2],     /* 1-right neighbors */
               c[3], w[3], d[3],     /* 2-right neighbors */
               dcdx[1]);             /* Output */
        
    for(i=2; i<n-2; i++)
    {
        ADVEC_DIFF(c[i-2], w[i-2], d[i-2],  /* 2-left neighbors */
                   c[i-1], w[i-1], d[i-1],  /* 1-left neighbors */
                   c[i],   w[i],   d[i],    /* Values */
                   c[i+1], w[i+1], d[i+1],  /* 1-right neighbors */
                   c[i+2], w[i+2], d[i+2],  /* 2-right neighbors */
                   dcdx[i]);                /* Output */
    }
    
    /* Do boundary cell c[n-2] explicitly */
    ADVEC_DIFF(c[n-4], w[n-4], d[n-4],  /* 2-left neighbors */
               c[n-3], w[n-3], d[n-3],  /* 1-left neighbors */
               c[n-2], w[n-2], d[n-2],  /* Values */
               cb[1],  wb[1],  db[1],   /* 1-right neighbors */
               cb[2],  wb[2],  db[2],   /* 2-right neighbors */
               dcdx[n-2]);              /* Output */
    
    /* Do boundary cell c[n-1] explicitly */
    ADVEC_DIFF(c[n-3], w[n-3], d[n-3],  /* 2-left neighbors */
               c[n-2], w[n-2], d[n-2],  /* 1-left neighbors */
               c[n-1], w[n-1], d[n-1],  /* Values */
               cb[2],  wb[2],  db[2],   /* 1-right neighbors */
               cb[3],  wb[3],  db[3],   /* 2-right neighbors */
               dcdx[n-1]);              /* Output */
}

void discretize(const int n, volatile double *conc_in, volatile double *wind, 
                volatile double *diff, double *concbound, double *windbound, 
                double *diffbound, double cell_size, double dt, 
                volatile double *conc_out)
{
    int i;
    double c[n];
    double dcdx[n];
    
    timer_start(&metrics.tot_discret);
    
    /* Copy original values (FIXME: SIMD, unroll) */
    for(i=0; i<n; i++)
    {
        c[i] = conc_out[i] = conc_in[i];
    }
    
    space_advec_diff_row(n, conc_in, wind, diff, concbound, windbound, diffbound, cell_size, dcdx);
    
    for(i=0; i<n; i++)
        c[i] += dt*dcdx[i];
    
    space_advec_diff_row(n, c, wind, diff, concbound, windbound, diffbound, cell_size, dcdx);
    
    for(i=0; i<n; i++)
        c[i] += dt*dcdx[i];
    
    for(i=0; i<n; i++)
    {
        conc_out[i] = 0.5 * (conc_out[i] + c[i]);
        if(conc_out[i] < 0.0)
            conc_out[i] = 0.0;
    }
    
    timer_stop(&metrics.tot_discret);
}
