#include "fixedgrid.h"
#include "discretize.h"
#include "timer.h"

extern fixedgrid_t G_GLOBAL;

/* 
 * The core upwinded advection/diffusion equation.
 * c = conc, w = wind, d = diff
 * x2l is the 2-left neighbor of x, etc.
 * x2r is the 2-right neighbor of x, etc.
 * out = output variable
 */
inline real_t 
advec_diff(real_t cell_size,
           real_t c2l, real_t w2l, real_t d2l, 
           real_t c1l, real_t w1l, real_t d1l, 
           real_t   c, real_t   w, real_t   d, 
           real_t c1r, real_t w1r, real_t d1r, 
           real_t c2r, real_t w2r, real_t d2r)
{
    real_t wind, diff_term, advec_term, advec_termL, advec_termR;
    
    wind = (w1l + w) / 2.0;
    if(wind >= 0.0) advec_termL = (1.0/6.0) * ( -c2l + 5.0*c1l + 2.0*c );
    else advec_termL = (1.0/6.0) * ( 2.0*c1l + 5.0*c - c1r );
    advec_termL *= wind;
    wind = (w1r + w) / 2.0;
    if(wind >= 0.0) advec_termR = (1.0/6.0) * ( -c1l + 5.0*c + 2.0*c1r );
    else advec_termR = (1.0/6.0) * ( 2.0*c + 5.0*c1r - c2r );
    advec_termR *= wind;
    advec_term = (advec_termL - advec_termR) / cell_size;
    diff_term = ( ((d1l+d)/2)*(c1l-c) - ((d+d1r)/2)*(c-c1r) ) / (cell_size * cell_size);
    return advec_term + diff_term;
}

/*
 * Applies the advection / diffusion equation to scalar data
 */
void space_advec_diff(const uint32_t n, 
                      real_t *c, 
                      real_t *w, 
                      real_t *d, 
                      real_t *cb, 
                      real_t *wb, 
                      real_t *db, 
                      real_t cell_size, 
                      real_t *dcdx)
{
    uint32_t i;
    
    /* Do boundary cell c[0] explicitly */
    dcdx[0] = advec_diff(cell_size,
                         cb[0], wb[0], db[0],  /* 2-left neighbors */
                         cb[1], wb[1], db[1],  /* 1-left neighbors */
                         c[0], w[0], d[0],     /* Values */
                         c[1], w[1], d[1],     /* 1-right neighbors */
                         c[2], w[2], d[2]);    /* 2-right neighbors */
    
    /* Do boundary cell c[1] explicitly */    
    dcdx[1] = advec_diff(cell_size,
                         cb[1], wb[1], db[1],  /* 2-left neighbors */
                         cb[2], wb[2], db[2],  /* 1-left neighbors */
                         c[1], w[1], d[1],     /* Values */
                         c[2], w[2], d[2],     /* 1-right neighbors */
                         c[3], w[3], d[3]);    /* 2-right neighbors */
    
    for(i=2; i<n-2; i++)
    {
        dcdx[i] = advec_diff(cell_size,
                             c[i-2], w[i-2], d[i-2],  /* 2-left neighbors */
                             c[i-1], w[i-1], d[i-1],  /* 1-left neighbors */
                             c[i],   w[i],   d[i],    /* Values */
                             c[i+1], w[i+1], d[i+1],  /* 1-right neighbors */
                             c[i+2], w[i+2], d[i+2]); /* 2-right neighbors */
    }
    
    /* Do boundary cell c[n-2] explicitly */
    dcdx[n-2] = advec_diff(cell_size,
                           c[n-4], w[n-4], d[n-4],  /* 2-left neighbors */
                           c[n-3], w[n-3], d[n-3],  /* 1-left neighbors */
                           c[n-2], w[n-2], d[n-2],  /* Values */
                           cb[1],  wb[1],  db[1],   /* 1-right neighbors */
                           cb[2],  wb[2],  db[2]);  /* 2-right neighbors */
    
    /* Do boundary cell c[n-1] explicitly */
    dcdx[n-1] = advec_diff(cell_size,
                           c[n-3], w[n-3], d[n-3],  /* 2-left neighbors */
                           c[n-2], w[n-2], d[n-2],  /* 1-left neighbors */
                           c[n-1], w[n-1], d[n-1],  /* Values */
                           cb[2],  wb[2],  db[2],   /* 1-right neighbors */
                           cb[3],  wb[3],  db[3]);  /* 2-right neighbors */
}


void discretize(const int n, real_t *conc_in, real_t *wind, 
                real_t *diff, real_t *concbound, real_t *windbound, 
                real_t *diffbound, real_t cell_size, real_t dt, 
                real_t *conc_out)
{
    int i;
    real_t c[n];
    real_t dcdx[n];
    
    for(i=0; i<n; i++)
    {
        c[i] = conc_out[i] = conc_in[i];
    }
    
    space_advec_diff(n, conc_in, wind, diff, concbound, windbound, diffbound, cell_size, dcdx);
    
    for(i=0; i<n; i++)
        c[i] += dt*dcdx[i];
    
    space_advec_diff(n, c, wind, diff, concbound, windbound, diffbound, cell_size, dcdx);
    
    for(i=0; i<n; i++)
        c[i] += dt*dcdx[i];
    
    for(i=0; i<n; i++)
    {
        conc_out[i] = 0.5 * (conc_out[i] + c[i]);
        if(conc_out[i] < 0.0)
            conc_out[i] = 0.0;
    }
}
