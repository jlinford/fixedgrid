#include "params.h"
#include "timer.h"

void space_advec_diff(const int n, const double *c, const double *w, const double *d, 
                      const double cell_size, double *dcdx)
{
    const int shift = 2;
    
    int i;
    double wind, diff_term;
    double advec_term, advec_termL, advec_termR;
    
    timer_start(TIMER_ADVEC_DIFF);
    
    for(i=shift; i<n+shift; i++)
    {
        wind = (w[i-1] + w[i]) / 2.0;
        if(wind > 0.0)
        {
            advec_termL = (1.0/6.0) * ( (-c[i-2]) + 5.0*c[i-1] + 2.0*c[i] );
        }
        else
        {
            advec_termL = (1.0/6.0) * ( 2.0*c[i-1] + 5.0*c[i] + (-c[i+1]) );
        }
        advec_termL *= wind;
            
        wind = (w[i+1] + w[i]) / 2.0;
        if(wind > 0.0)
        {
            advec_termR = (1.0/6.0) * ( (-c[i-1]) + 5.0*c[i] + 2.0*c[i+1] );
        }
        else
        {
            advec_termR = (1.0/6.0) * ( 2.0*c[i] + 5.0*c[i+1] + (-c[i+2]) );
        }
        advec_termR *= wind;
        
        advec_term = (advec_termL - advec_termR) / cell_size;
        
        diff_term = ( ((d[i-1]+d[i])/2)*(c[i-1]-c[i]) - ((d[i]+d[i+1])/2)*(c[i]-c[i+1]) ) / (cell_size * cell_size);
     
        dcdx[i-shift] = advec_term + diff_term;
    }
    
    timer_stop(TIMER_ADVEC_DIFF);
}

void discretize(const int n, const double *conc_in, const double *wind, 
                const double *diff, const double cell_size, const double dt, 
                double *conc_out)
{
    int k;
    int shift;
    double c[NX+4], c1[NX+4], w[NX+4], d[NX+4];
    double dcdx[NX];
    
    timer_start(TIMER_DISCRETIZE);
    
    shift = 2;
    
    /* Shift vectors to fill margin cells */
    for(k=0; k<shift; k++)
    {
        c[k] = conc_in[n-k-1];
        w[k] = wind[n-k-1];
        d[k] = diff[n-k-1];
    }
    for(k=0; k<n; k++)
    {
        c[k+shift] = conc_in[k];
        w[k+shift] = wind[k];
        d[k+shift] = diff[k];
    }
    for(k=0; k<shift; k++)
    {
        c[n+shift+k] = conc_in[k];
        w[n+shift+k] = wind[k];
        d[n+shift+k] = diff[k];
    }
    
    for(k=0; k<n+4; k++)
        c1[k] = c[k];
    
    space_advec_diff(n, c, w, d, cell_size, dcdx);
    
    for(k=0; k<n; k++)
        c1[k+shift] += dt*dcdx[k];
    
    space_advec_diff(n, c1, w, d, cell_size, dcdx);
    
    for(k=0; k<n; k++)
        c1[k+shift] += dt*dcdx[k];
    
    for(k=0; k<n; k++)
        c[k+shift] = 0.5 * (c[k+shift] + c1[k+shift]);
    
    for(k=0; k<n; k++)
        c[k+shift] = c[k+shift] > 0.0 ? c[k+shift] : 0.0;
    
    for(k=0; k<n; k++)
        conc_out[k] = c[k+shift];
    
    timer_stop(TIMER_DISCRETIZE);
}
