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
#include <simdmath.h>
#include <spu_intrinsics.h>

#include "fixedgrid_spu.h"
#include "discretize.h"
#include "timer.h"
#include "common.h"
#include "params.h"

/* 
 * The core upwinded advection/diffusion equation.
 * c = conc, w = wind, d = diff
 * x2l is the 2-left neighbor of x, etc.
 * x2r is the 2-right neighbor of x, etc.
 * out = output variable
 */
/*
#define _ADVEC_DIFF(c2l, w2l, d2l, c1l, w1l, d1l, c, w, d, c1r, w1r, d1r, c2r, w2r, d2r, out) {  \
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

/* Same equation but vectorized */
/*
#define _ADVEC_DIFF_V(c2l, w2l, d2l, c1l, w1l, d1l, c, w, d, c1r, w1r, d1r, c2r, w2r, d2r, out) { \
    acc1 = spu_add(w1l, w);                         \
    wind = spu_mul(acc1, HALF);                     \
    acc1 = spu_mul(c1l, FIVE);                      \
    acc2 = spu_mul(c, TWO);                         \
    advec_term_pos = spu_add(acc1, acc2);           \
    advec_term_pos = spu_sub(advec_term_pos, c2l);  \
    acc1 = spu_mul(c1l, TWO);                       \
    acc2 = spu_mul(c, FIVE);                        \
    advec_term_neg = spu_add(acc1, acc2);           \
    advec_term_neg = spu_sub(advec_term_neg, c1r);  \
    acc1 = (vector real_t)spu_cmpgt(wind, ZERO);   \
    acc1 = spu_and(acc1, advec_term_pos);           \
    acc2 = (vector real_t)spu_cmpgt(ZERO, wind);   \
    acc2 = spu_and(acc2, advec_term_neg);           \
    advec_termL = spu_add(acc1, acc2);              \
    advec_termL = spu_mul(advec_termL, SIXTH);      \
    advec_termL = spu_mul(advec_termL, wind);       \
    acc1 = spu_add(w1r, w);                         \
    wind = spu_mul(acc1, HALF);                     \
    acc1 = spu_mul(c, FIVE);                        \
    acc2 = spu_mul(c1r, TWO);                       \
    advec_term_pos = spu_add(acc1, acc2);           \
    advec_term_pos = spu_sub(advec_term_pos, c1l);  \
    acc1 = spu_mul(c, TWO);                         \
    acc2 = spu_mul(c1r, FIVE);                      \
    advec_term_neg = spu_add(acc1, acc2);           \
    advec_term_neg = spu_sub(advec_term_neg, c2r);  \
    acc1 = (vector real_t)spu_cmpgt(wind, ZERO);   \
    acc1 = spu_and(acc1, advec_term_pos);           \
    acc2 = (vector real_t)spu_cmpgt(ZERO, wind);   \
    acc2 = spu_and(acc2, advec_term_neg);           \
    advec_termR = spu_add(acc1, acc2);              \
    advec_termR = spu_mul(advec_termR, SIXTH);      \
    advec_termR = spu_mul(advec_termR, wind);       \
    acc1 = spu_sub(advec_termL, advec_termR);       \
    advec_term = VEC_DIVIDE(acc1, cell_size);       \
    acc1 = spu_add(d1l, d);                         \
    acc1 = spu_mul(acc1, HALF);                     \
    acc3 = spu_sub(c1l, c);                         \
    acc1 = spu_mul(acc1, acc3);                     \
    acc2 = spu_add(d, d1r);                         \
    acc2 = spu_mul(acc2, HALF);                     \
    acc3 = spu_sub(c, c1r);                         \
    acc2 = spu_mul(acc2, acc3);                     \
    acc1 = spu_sub(acc1, acc2);                     \
    acc2 = spu_mul(cell_size, cell_size);           \
    diff_term = VEC_DIVIDE(acc1, acc2);             \
    out = spu_add(advec_term, diff_term);           \
}
 */

inline vector real_t
advec_diff_v(vector real_t cell_size,
             vector real_t c2l, vector real_t w2l, vector real_t d2l, 
             vector real_t c1l, vector real_t w1l, vector real_t d1l, 
             vector real_t   c, vector real_t   w, vector real_t   d, 
             vector real_t c1r, vector real_t w1r, vector real_t d1r, 
             vector real_t c2r, vector real_t w2r, vector real_t d2r)
{
    const vector real_t FIVE  = spu_splats((real_t)5.0);
    const vector real_t TWO   = spu_splats((real_t)2.0);
    const vector real_t ZERO  = spu_splats((real_t)0.0);
    const vector real_t HALF  = spu_splats((real_t)0.5);
    const vector real_t SIXTH = spu_splats((real_t)1.0/(real_t)6.0);
    vector real_t acc1, acc2, acc3;
    vector real_t wind, diff_term, advec_term;
    vector real_t advec_term_pos, advec_term_neg;
    vector real_t advec_termR, advec_termL;
    
    acc1 = spu_add(w1l, w);
    wind = spu_mul(acc1, HALF);
    acc1 = spu_mul(c1l, FIVE);
    acc2 = spu_mul(c, TWO);
    advec_term_pos = spu_add(acc1, acc2);
    advec_term_pos = spu_sub(advec_term_pos, c2l);
    acc1 = spu_mul(c1l, TWO);
    acc2 = spu_mul(c, FIVE);
    advec_term_neg = spu_add(acc1, acc2);
    advec_term_neg = spu_sub(advec_term_neg, c1r);
    acc1 = (vector real_t)spu_cmpgt(wind, ZERO);
    acc1 = spu_and(acc1, advec_term_pos);
    acc2 = (vector real_t)spu_cmpgt(ZERO, wind);
    acc2 = spu_and(acc2, advec_term_neg);
    advec_termL = spu_add(acc1, acc2);
    advec_termL = spu_mul(advec_termL, SIXTH);
    advec_termL = spu_mul(advec_termL, wind);
    acc1 = spu_add(w1r, w);
    wind = spu_mul(acc1, HALF);
    acc1 = spu_mul(c, FIVE);
    acc2 = spu_mul(c1r, TWO);
    advec_term_pos = spu_add(acc1, acc2);
    advec_term_pos = spu_sub(advec_term_pos, c1l);
    acc1 = spu_mul(c, TWO);
    acc2 = spu_mul(c1r, FIVE);
    advec_term_neg = spu_add(acc1, acc2);
    advec_term_neg = spu_sub(advec_term_neg, c2r);
    acc1 = (vector real_t)spu_cmpgt(wind, ZERO);
    acc1 = spu_and(acc1, advec_term_pos);
    acc2 = (vector real_t)spu_cmpgt(ZERO, wind);
    acc2 = spu_and(acc2, advec_term_neg);
    advec_termR = spu_add(acc1, acc2);
    advec_termR = spu_mul(advec_termR, SIXTH);
    advec_termR = spu_mul(advec_termR, wind);
    acc1 = spu_sub(advec_termL, advec_termR);
    advec_term = VEC_DIVIDE(acc1, cell_size);
    acc1 = spu_add(d1l, d);
    acc1 = spu_mul(acc1, HALF);
    acc3 = spu_sub(c1l, c);
    acc1 = spu_mul(acc1, acc3);
    acc2 = spu_add(d, d1r);
    acc2 = spu_mul(acc2, HALF);
    acc3 = spu_sub(c, c1r);
    acc2 = spu_mul(acc2, acc3);
    acc1 = spu_sub(acc1, acc2);
    acc2 = spu_mul(cell_size, cell_size);
    diff_term = VEC_DIVIDE(acc1, acc2);
    return spu_add(advec_term, diff_term);
}


/*
 * Applies the advection / diffusion equation to scalar data
 */
void space_advec_diff_s(const uint32_t n, 
                        volatile real_t *c, 
                        volatile real_t *w, 
                        volatile real_t *d, 
                        real_t *cb, 
                        real_t *wb, 
                        real_t *db, 
                        real_t cell_size, 
                        volatile real_t *dcdx)
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

/*
 * Applies the advection / diffusion equation to vector data
 */
void space_advec_diff_v(const uint32_t n, 
                        volatile vector real_t *c, 
                        volatile vector real_t *w, 
                        volatile vector real_t *d, 
                        vector real_t *cb, 
                        vector real_t *wb, 
                        vector real_t *db, 
                        vector real_t cell_size, 
                        volatile vector real_t *dcdx)
{    
    uint32_t i;
    
    /* Do boundary cell c[0] explicitly */
    dcdx[0] = advec_diff_v(cell_size,
                         cb[0], wb[0], db[0],  /* 2-left neighbors */
                         cb[1], wb[1], db[1],  /* 1-left neighbors */
                         c[0], w[0], d[0],     /* Values */
                         c[1], w[1], d[1],     /* 1-right neighbors */
                         c[2], w[2], d[2]);    /* 2-right neighbors */
    
    /* Do boundary cell c[1] explicitly */    
    dcdx[1] = advec_diff_v(cell_size,
                         cb[1], wb[1], db[1],  /* 2-left neighbors */
                         cb[2], wb[2], db[2],  /* 1-left neighbors */
                         c[1], w[1], d[1],     /* Values */
                         c[2], w[2], d[2],     /* 1-right neighbors */
                         c[3], w[3], d[3]);    /* 2-right neighbors */
    
    for(i=2; i<n-2; i++)
    {
        dcdx[i] = advec_diff_v(cell_size,
                             c[i-2], w[i-2], d[i-2],  /* 2-left neighbors */
                             c[i-1], w[i-1], d[i-1],  /* 1-left neighbors */
                             c[i],   w[i],   d[i],    /* Values */
                             c[i+1], w[i+1], d[i+1],  /* 1-right neighbors */
                             c[i+2], w[i+2], d[i+2]); /* 2-right neighbors */
    }
    
    /* Do boundary cell c[n-2] explicitly */
    dcdx[n-2] = advec_diff_v(cell_size,
                           c[n-4], w[n-4], d[n-4],  /* 2-left neighbors */
                           c[n-3], w[n-3], d[n-3],  /* 1-left neighbors */
                           c[n-2], w[n-2], d[n-2],  /* Values */
                           cb[1],  wb[1],  db[1],   /* 1-right neighbors */
                           cb[2],  wb[2],  db[2]);  /* 2-right neighbors */
    
    /* Do boundary cell c[n-1] explicitly */
    dcdx[n-1] = advec_diff_v(cell_size,
                           c[n-3], w[n-3], d[n-3],  /* 2-left neighbors */
                           c[n-2], w[n-2], d[n-2],  /* 1-left neighbors */
                           c[n-1], w[n-1], d[n-1],  /* Values */
                           cb[2],  wb[2],  db[2],   /* 1-right neighbors */
                           cb[3],  wb[3],  db[3]);  /* 2-right neighbors */
}

void discretize_row(const uint32_t n, 
                    volatile real_t *conc_in, 
                    volatile real_t *wind, 
                    volatile real_t *diff, 
                    real_t *concbound, 
                    real_t *windbound, 
                    real_t *diffbound, 
                    real_t cell_size, 
                    real_t dt, 
                    volatile real_t *conc_out)
{
    uint32_t i;
    real_t c[n];
    real_t dcdx[n];
    
    timer_start(&metrics.tot_discret);
    
    /* Copy original values (FIXME: SIMD, unroll) */
    for(i=0; i<n; i++)
    {
        c[i] = conc_out[i] = conc_in[i];
    }
    
    space_advec_diff_s(n, conc_in, wind, diff, concbound, windbound, diffbound, cell_size, dcdx);
    
    for(i=0; i<n; i++)
    {
        c[i] += dt*dcdx[i];
    }
    
    space_advec_diff_s(n, c, wind, diff, concbound, windbound, diffbound, cell_size, dcdx);
    
    for(i=0; i<n; i++)
    {
        c[i] += dt*dcdx[i];
    }
    
    for(i=0; i<n; i++)
    {
        conc_out[i] = 0.5 * (conc_out[i] + c[i]);
        if(conc_out[i] < 0.0) 
            conc_out[i] = 0.0;
    }
    
    timer_stop(&metrics.tot_discret);
}

void discretize_col(const uint32_t n, 
                    volatile vector real_t *conc_in, 
                    volatile vector real_t *wind, 
                    volatile vector real_t *diff, 
                    vector real_t *concbound, 
                    vector real_t *windbound, 
                    vector real_t *diffbound, 
                    vector real_t cell_size, 
                    vector real_t dt, 
                    volatile vector real_t *conc_out)
{
    const vector real_t HALF = spu_splats((real_t)0.5);
    
    uint32_t i;
    vector real_t acc;
    vector real_t c[n];
    vector real_t dcdx[n];
        
    timer_start(&metrics.tot_discret);
    
    /* Copy original values (FIXME: unroll) */
    for(i=0; i<n; i++)
    {
        c[i] = conc_out[i] = conc_in[i];
    }
    
    space_advec_diff_v(n, conc_in, wind, diff, concbound, windbound, diffbound, cell_size, dcdx);
    
    for(i=0; i<n; i++)
    {
        //c[i] += dt*dcdx[i];
        c[i] = spu_madd(dt, dcdx[i], c[i]);
    }
    
    space_advec_diff_v(n, c, wind, diff, concbound, windbound, diffbound, cell_size, dcdx);
    
    for(i=0; i<n; i++)
    {
        //c[i] += dt*dcdx[i];
        c[i] = spu_madd(dt, dcdx[i], c[i]);
    }
    
    for(i=0; i<n; i++)
    {
        //conc_out[i] = 0.5 * (conc_out[i] + c[i]);
        acc = spu_add(conc_out[i], c[i]);
        conc_out[i] = spu_mul(HALF, acc);
        
        //if(conc_out[i] < 0.0) conc_out[i] = 0.0;
        acc = spu_splats((real_t)0.0);
        acc = (vector real_t)spu_cmpgt(conc_out[i], acc);
        conc_out[i] = spu_and(conc_out[i], acc);
    }
    
    timer_stop(&metrics.tot_discret);
}

#undef _ADVEC_DIFF
#undef _ADVEC_DIFF_V
