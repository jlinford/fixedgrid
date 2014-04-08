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

#include "discretize.h"
#include "timer.h"
#include "common.h"
#include "params.h"

extern metrics_t metrics;

/* 
 * The core upwinded advection/diffusion equation.
 * c = conc, w = wind, d = diff
 * x2l is the 2-left neighbor of x, etc.
 * x2r is the 2-right neighbor of x, etc.
 * out = output variable
 */
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
/* Same equation but vectorized */
#define _ADVEC_DIFF_V(c2l, w2l, d2l, c1l, w1l, d1l, c, w, d, c1r, w1r, d1r, c2r, w2r, d2r, out) { \
    acc1 = spu_add(w1l, w);                                 \
    wind = spu_mul(acc1, HALF);                             \
    acc1 = spu_splats(5.0);                                 \
    acc1 = spu_mul(acc1, c1l);                              \
    acc2 = spu_splats(2.0);                                 \
    acc2 = spu_mul(acc2, c);                                \
    advec_term_pos = spu_add(acc1, acc2);                   \
    advec_term_pos = spu_sub(advec_term_pos, c2l);          \
    advec_term_pos = spu_mul(advec_term_pos, SIXTH);        \
    acc1 = spu_splats(2.0);                                 \
    acc1 = spu_mul(acc1, c1l);                              \
    acc2 = spu_splats(5.0);                                 \
    acc2 = spu_mul(acc2, c);                                \
    advec_term_neg = spu_add(acc1, acc2);                   \
    advec_term_neg = spu_sub(advec_term_neg, c1r);          \
    advec_term_neg = spu_mul(advec_term_neg, SIXTH);        \
    acc1 = spu_splats(0.0);                                 \
    acc1 = (vector value_t)spu_cmpgt(wind, acc1);           \
    acc1 = spu_and(acc1, advec_term_pos);                   \
    acc2 = spu_splats(0.0);                                 \
    acc2 = (vector value_t)spu_cmpgt(acc2, wind);           \
    acc2 = spu_and(acc2, advec_term_neg);                   \
    advec_termL = spu_add(acc1, acc2);                      \
    advec_termL = spu_mul(advec_termL, wind);               \
    acc1 = spu_add(w1r, w);                                 \
    wind = spu_mul(acc1, HALF);                             \
    acc1 = spu_splats(5.0);                                 \
    acc1 = spu_mul(acc1, c);                                \
    acc2 = spu_splats(2.0);                                 \
    acc2 = spu_mul(acc2, c1r);                              \
    advec_term_pos = spu_add(acc1, acc2);                   \
    advec_term_pos = spu_sub(advec_term_pos, c1l);          \
    advec_term_pos = spu_mul(advec_term_pos, SIXTH);        \
    acc1 = spu_splats(2.0);                                 \
    acc1 = spu_mul(acc1, c);                                \
    acc2 = spu_splats(5.0);                                 \
    acc2 = spu_mul(acc2, c1r);                              \
    advec_term_neg = spu_add(acc1, acc2);                   \
    advec_term_neg = spu_sub(advec_term_neg, c2r);          \
    advec_term_neg = spu_mul(advec_term_neg, SIXTH);        \
    acc1 = spu_splats(0.0);                                 \
    acc1 = (vector value_t)spu_cmpgt(wind, acc1);           \
    acc1 = spu_and(acc1, advec_term_pos);                   \
    acc2 = spu_splats(0.0);                                 \
    acc2 = (vector value_t)spu_cmpgt(acc2, wind);           \
    acc2 = spu_and(acc2, advec_term_neg);                   \
    advec_termR = spu_add(acc1, acc2);                      \
    advec_termR = spu_mul(advec_termR, wind);               \
    acc1 = spu_sub(advec_termL, advec_termR);               \
    advec_term = VEC_DIVIDE(acc1, cell_size);               \
    acc1 = spu_add(d1l, d);                                 \
    acc1 = spu_mul(acc1, HALF);                             \
    acc3 = spu_sub(c1l, c);                                 \
    acc1 = spu_mul(acc1, acc3);                             \
    acc2 = spu_add(d, d1r);                                 \
    acc2 = spu_mul(acc2, HALF);                             \
    acc3 = spu_sub(c, c1r);                                 \
    acc2 = spu_mul(acc2, acc3);                             \
    acc1 = spu_sub(acc1, acc2);                             \
    acc2 = spu_mul(cell_size, cell_size);                   \
    diff_term = VEC_DIVIDE(acc1, acc2);                     \
    out = spu_add(advec_term, diff_term);                   \
}

/*
inline void advec_diff_v(const vector value_t cell_size,
                         vector value_t c2l, vector value_t w2l, vector value_t d2l,
                         vector value_t c1l, vector value_t w1l, vector value_t d1l,
                         vector value_t   c, vector value_t   w, vector value_t   d,
                         vector value_t c1r, vector value_t w1r, vector value_t d1r,
                         vector value_t c2r, vector value_t w2r, vector value_t d2r,
                         volatile vector value_t* out)
{    
    const vector value_t HALF  = spu_splats(0.5);
    const vector value_t SIXTH = spu_splats(1.0/6.0);
    vector value_t acc1, acc2, acc3;
    vector value_t wind, diff_term, advec_term;
    vector value_t advec_term_pos, advec_term_neg;
    vector value_t advec_termR, advec_termL;    
    
    //wind = (w1l + w) / 2.0;
    acc1 = spu_add(w1l, w);
    wind = spu_mul(acc1, HALF);
    
    //if(wind >= 0.0) advec_termL = (1.0/6.0) * ( 5.0*c1l + 2.0*c - c2l );
    acc1 = spu_splats(5.0);
    acc1 = spu_mul(acc1, c1l);
    acc2 = spu_splats(2.0);
    acc2 = spu_mul(acc2, c);
    advec_term_pos = spu_add(acc1, acc2);
    advec_term_pos = spu_sub(advec_term_pos, c2l);
    advec_term_pos = spu_mul(advec_term_pos, SIXTH);
    
    //else advec_termL = (1.0/6.0) * ( 2.0*c1l + 5.0*c - c1r );
    acc1 = spu_splats(2.0);
    acc1 = spu_mul(acc1, c1l);
    acc2 = spu_splats(5.0);
    acc2 = spu_mul(acc2, c);
    advec_term_neg = spu_add(acc1, acc2);
    advec_term_neg = spu_sub(advec_term_neg, c1r);
    advec_term_neg = spu_mul(advec_term_neg, SIXTH);
    
    // Mask vectors together to select correct upwinded values
    acc1 = spu_splats(0.0);
    acc1 = (vector value_t)spu_cmpgt(wind, acc1);
    acc1 = spu_and(acc1, advec_term_pos);
    acc2 = spu_splats(0.0);
    acc2 = (vector value_t)spu_cmpgt(acc2, wind);
    acc2 = spu_and(acc2, advec_term_neg);
    advec_termL = spu_add(acc1, acc2);
    
    //advec_termL *= wind;
    advec_termL = spu_mul(advec_termL, wind);
    
    //wind = (w1r + w) / 2.0;
    acc1 = spu_add(w1r, w);
    wind = spu_mul(acc1, HALF);
    
    //if(wind >= 0.0) advec_termR = (1.0/6.0) * ( 5.0*c + 2.0*c1r - c1l );
    acc1 = spu_splats(5.0);
    acc1 = spu_mul(acc1, c);
    acc2 = spu_splats(2.0);
    acc2 = spu_mul(acc2, c1r);
    advec_term_pos = spu_add(acc1, acc2);
    advec_term_pos = spu_sub(advec_term_pos, c1l);
    advec_term_pos = spu_mul(advec_term_pos, SIXTH);
    
    //else advec_termR = (1.0/6.0) * ( 2.0*c + 5.0*c1r - c2r );
    acc1 = spu_splats(2.0);
    acc1 = spu_mul(acc1, c);
    acc2 = spu_splats(5.0);
    acc2 = spu_mul(acc2, c1r);
    advec_term_neg = spu_add(acc1, acc2);
    advec_term_neg = spu_sub(advec_term_neg, c2r);
    advec_term_neg = spu_mul(advec_term_neg, SIXTH);
    
    // Mask vectors together to select correct upwinded values
    acc1 = spu_splats(0.0);
    acc1 = (vector value_t)spu_cmpgt(wind, acc1);
    acc1 = spu_and(acc1, advec_term_pos);
    acc2 = spu_splats(0.0);
    acc2 = (vector value_t)spu_cmpgt(acc2, wind);
    acc2 = spu_and(acc2, advec_term_neg);
    advec_termR = spu_add(acc1, acc2);
    
    //advec_termR *= wind;
    advec_termR = spu_mul(advec_termR, wind);
    
    //advec_term = (advec_termL - advec_termR) / cell_size;
    acc1 = spu_sub(advec_termL, advec_termR);
    advec_term = VEC_DIVIDE(acc1, cell_size);

    //diff_term = ( ((d1l+d)/2)*(c1l-c) - ((d+d1r)/2)*(c-c1r) ) / (cell_size * cell_size);
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

    //out = advec_term + diff_term;
    *out = spu_add(advec_term, diff_term);
}
*/

/*
 * Applies the advection / diffusion equation to scalar data
 */
void space_advec_diff_s(int n, 
                        volatile value_t *c, 
                        volatile value_t *w, 
                        volatile value_t *d, 
                        value_t *cb, 
                        value_t *wb, 
                        value_t *db, 
                        value_t cell_size, 
                        volatile value_t *dcdx)
{
    int i;
    value_t wind, diff_term, advec_term, advec_termL, advec_termR;
    
    /* Do boundary cell c[0] explicitly */
    _ADVEC_DIFF(cb[0], wb[0], db[0],  /* 2-left neighbors */
                cb[1], wb[1], db[1],  /* 1-left neighbors */
                c[0], w[0], d[0],     /* Values */
                c[1], w[1], d[1],     /* 1-right neighbors */
                c[2], w[2], d[2],     /* 2-right neighbors */
                dcdx[0]);             /* Output */
    
    /* Do boundary cell c[1] explicitly */
    _ADVEC_DIFF(cb[1], wb[1], db[1],  /* 2-left neighbors */
                cb[2], wb[2], db[2],  /* 1-left neighbors */
                c[1], w[1], d[1],     /* Values */
                c[2], w[2], d[2],     /* 1-right neighbors */
                c[3], w[3], d[3],     /* 2-right neighbors */
                dcdx[1]);             /* Output */
        
    for(i=2; i<n-2; i++)
    {
        _ADVEC_DIFF(c[i-2], w[i-2], d[i-2],  /* 2-left neighbors */
                    c[i-1], w[i-1], d[i-1],  /* 1-left neighbors */
                    c[i],   w[i],   d[i],    /* Values */
                    c[i+1], w[i+1], d[i+1],  /* 1-right neighbors */
                    c[i+2], w[i+2], d[i+2],  /* 2-right neighbors */
                    dcdx[i]);                /* Output */
    }
    
    /* Do boundary cell c[n-2] explicitly */
    _ADVEC_DIFF(c[n-4], w[n-4], d[n-4],  /* 2-left neighbors */
                c[n-3], w[n-3], d[n-3],  /* 1-left neighbors */
                c[n-2], w[n-2], d[n-2],  /* Values */
                cb[1],  wb[1],  db[1],   /* 1-right neighbors */
                cb[2],  wb[2],  db[2],   /* 2-right neighbors */
                dcdx[n-2]);              /* Output */
    
    /* Do boundary cell c[n-1] explicitly */
    _ADVEC_DIFF(c[n-3], w[n-3], d[n-3],  /* 2-left neighbors */
                c[n-2], w[n-2], d[n-2],  /* 1-left neighbors */
                c[n-1], w[n-1], d[n-1],  /* Values */
                cb[2],  wb[2],  db[2],   /* 1-right neighbors */
                cb[3],  wb[3],  db[3],   /* 2-right neighbors */
                dcdx[n-1]);              /* Output */
}

/*
 * Applies the advection / diffusion equation to vector data
 */
void space_advec_diff_v(int n, 
                        volatile vector value_t *c, 
                        volatile vector value_t *w, 
                        volatile vector value_t *d, 
                        vector value_t *cb, 
                        vector value_t *wb, 
                        vector value_t *db, 
                        vector value_t cell_size, 
                        volatile vector value_t *dcdx)
{
    int i;
    const vector value_t HALF  = spu_splats(0.5);
    const vector value_t SIXTH = spu_splats(1.0/6.0);
    vector value_t acc1, acc2, acc3;
    vector value_t wind, diff_term, advec_term;
    vector value_t advec_term_pos, advec_term_neg;
    vector value_t advec_termR, advec_termL;    
    
    /* Do boundary cell c[0] explicitly */
    _ADVEC_DIFF_V(cb[0], wb[0], db[0],  /* 2-left neighbors */
                  cb[1], wb[1], db[1],  /* 1-left neighbors */
                  c[0], w[0], d[0],     /* Values */
                  c[1], w[1], d[1],     /* 1-right neighbors */
                  c[2], w[2], d[2],     /* 2-right neighbors */
                  dcdx[0]);             /* Output */
    
    /* Do boundary cell c[1] explicitly */
    _ADVEC_DIFF_V(cb[1], wb[1], db[1],  /* 2-left neighbors */
                  cb[2], wb[2], db[2],  /* 1-left neighbors */
                  c[1], w[1], d[1],     /* Values */
                  c[2], w[2], d[2],     /* 1-right neighbors */
                  c[3], w[3], d[3],     /* 2-right neighbors */
                  dcdx[1]);             /* Output */
    
    for(i=2; i<n-2; i++)
    {
        _ADVEC_DIFF_V(c[i-2], w[i-2], d[i-2],  /* 2-left neighbors */
                      c[i-1], w[i-1], d[i-1],  /* 1-left neighbors */
                      c[i],   w[i],   d[i],    /* Values */
                      c[i+1], w[i+1], d[i+1],  /* 1-right neighbors */
                      c[i+2], w[i+2], d[i+2],  /* 2-right neighbors */
                      dcdx[i]);                /* Output */
    }
    
    /* Do boundary cell c[n-2] explicitly */
    _ADVEC_DIFF_V(c[n-4], w[n-4], d[n-4],  /* 2-left neighbors */
                  c[n-3], w[n-3], d[n-3],  /* 1-left neighbors */
                  c[n-2], w[n-2], d[n-2],  /* Values */
                  cb[1],  wb[1],  db[1],   /* 1-right neighbors */
                  cb[2],  wb[2],  db[2],   /* 2-right neighbors */
                  dcdx[n-2]);              /* Output */
    
    /* Do boundary cell c[n-1] explicitly */
    _ADVEC_DIFF_V(c[n-3], w[n-3], d[n-3],  /* 2-left neighbors */
                  c[n-2], w[n-2], d[n-2],  /* 1-left neighbors */
                  c[n-1], w[n-1], d[n-1],  /* Values */
                  cb[2],  wb[2],  db[2],   /* 1-right neighbors */
                  cb[3],  wb[3],  db[3],   /* 2-right neighbors */
                  dcdx[n-1]);              /* Output */
}

void discretize_row(const int n, 
                    volatile value_t *conc_in, 
                    volatile value_t *wind, 
                    volatile value_t *diff, 
                    value_t *concbound, 
                    value_t *windbound, 
                    value_t *diffbound, 
                    value_t cell_size, 
                    value_t dt, 
                    volatile value_t *conc_out)
{
    int i;
    value_t c[n];
    value_t dcdx[n];
    
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

void discretize_col(const int n, 
                    volatile vector value_t *conc_in, 
                    volatile vector value_t *wind, 
                    volatile vector value_t *diff, 
                    vector value_t *concbound, 
                    vector value_t *windbound, 
                    vector value_t *diffbound, 
                    vector value_t cell_size, 
                    vector value_t dt, 
                    volatile vector value_t *conc_out)
{
    const vector value_t HALF = spu_splats(0.5);
    
    int i;
    vector value_t acc;
    vector value_t c[n];
    vector value_t dcdx[n];
        
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
        acc = spu_splats(0.0);
        acc = (vector value_t)spu_cmpgt(conc_out[i], acc);
        conc_out[i] = spu_and(conc_out[i], acc);
    }
    
    timer_stop(&metrics.tot_discret);
}

#undef _ADVEC_DIFF
#undef _ADVEC_DIFF_V
