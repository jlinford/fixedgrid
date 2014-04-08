/*
 *  transport.c
 *  fixedgrid_serial
 *
 *  Created by John Linford on 6/23/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include "transport.h"
#include "discretize.h"

/**
 * Discretize rows
 */
void discretize_all_x(fixedgrid_t* G, real_t dt)
{
#if DO_X_DISCRET == 1
    
    int32_t x, y, z, s;
    
    real_t buff[NX];
    
    /* Boundary values */
    real_t cbound[4];
    real_t wbound[4];
    real_t dbound[4];
    
    timer_start(&G->metrics.x_discret);
    
    #pragma omp for private(z, y, x, s, buff, cbound, wbound, dbound)
    for(z=0; z<NZ; z++)
    {
        for(y=0; y<NY; y++)
        {
            for(s=0; s<NLOOKAT; s++)
            {
                cbound[0] = G->conc(NX-2, y, z, s);
                cbound[1] = G->conc(NX-1, y, z, s);
                cbound[2] = G->conc(0, y, z, s);
                cbound[3] = G->conc(1, y, z, s);
                wbound[0] = G->wind_u(NX-2, y, z);
                wbound[1] = G->wind_u(NX-1, y, z);
                wbound[2] = G->wind_u(0, y, z);
                wbound[3] = G->wind_u(1, y, z);
                dbound[0] = G->diff_h(NX-2, y, z);
                dbound[1] = G->diff_h(NX-1, y, z);
                dbound[2] = G->diff_h(0, y, z);
                dbound[3] = G->diff_h(1, y, z);
                
                discretize(NX, 
                           &G->conc(0, y, z, s), 
                           &G->wind_u(0, y, z),
                           &G->diff_h(0, y, z),
                           cbound, wbound, dbound, 
                           DX, dt, buff);
                
                timer_start(&G->metrics.array_copy);
                for(x=0; x<NX; x++)
                    G->conc(x, y, z, s) = buff[x];
                timer_stop(&G->metrics.array_copy);
            }
        }
    }
    
    timer_stop(&G->metrics.x_discret);
    
#endif
}

/**
 * Discretize y
 */
void discretize_all_y(fixedgrid_t* G, real_t dt)
{
#if DO_Y_DISCRET == 1
    
    int32_t x, y, z, s;
    
    /* Buffers */
    real_t ccol1[NY];
    real_t ccol2[NY];
    real_t wcol[NY];
    real_t dcol[NY];
    
    /* Boundary values */
    real_t cbound[4];
    real_t wbound[4];
    real_t dbound[4];
    
    timer_start(&G->metrics.y_discret);
    
    #pragma omp for private(z, y, x, s, ccol1, ccol2, wcol, dcol, cbound, wbound, dbound)
    for(z=0; z<NZ; z++)
    {
        for(x=0; x<NX; x++)
        {
            for(s=0; s<NLOOKAT; s++)
            {
                timer_start(&G->metrics.array_copy);
                for(y=0; y<NY; y++)
                {
                    ccol1[y] = G->conc(x, y, z, s);
                    wcol[y]  = G->wind_v(x, y, z);
                    dcol[y]  = G->diff_h(x, y, z);
                }
                timer_stop(&G->metrics.array_copy);
                
                cbound[0] = ccol1[NY-2];
                cbound[1] = ccol1[NY-1];
                cbound[2] = ccol1[0];
                cbound[3] = ccol1[1];
                wbound[0] = wcol[NY-2];
                wbound[1] = wcol[NY-1];
                wbound[2] = wcol[0];
                wbound[3] = wcol[1];
                dbound[0] = dcol[NY-2];
                dbound[1] = dcol[NY-1];
                dbound[2] = dcol[0];
                dbound[3] = dcol[1];
                
                discretize(NY, 
                           ccol1, wcol, dcol, 
                           cbound, wbound, dbound, 
                           DY, dt, ccol2);
                
                timer_start(&G->metrics.array_copy);
                for(y=0; y<NY; y++)
                    G->conc(x, y, z, s) = ccol2[y];
                timer_stop(&G->metrics.array_copy);
            }
        }
    }
    
    timer_stop(&G->metrics.y_discret);
    
#endif
}

/**
 * Discretize z
 */
void discretize_all_z(fixedgrid_t* G, real_t dt)
{
#if DO_Z_DISCRET == 1
    
    int32_t x, y, z, s;
    
    /* Buffers */
    real_t ccol1[NZ];
    real_t ccol2[NZ];
    real_t wcol[NZ];
    real_t dcol[NZ];
    
    /* Boundary values */
    real_t cbound[4];
    real_t wbound[4];
    real_t dbound[4];
    
    timer_start(&G->metrics.z_discret);
    
    #pragma omp for private(z, y, x, s, ccol1, ccol2, wcol, dcol, cbound, wbound, dbound)
    for(y=0; y<NY; y++)
    {
        for(x=0; x<NX; x++)
        {
            for(s=0; s<NLOOKAT; s++)
            {
                timer_start(&G->metrics.array_copy);
                for(z=0; z<NZ; z++)
                {
                    ccol1[z] = G->conc(x, y, z, s);
                    wcol[z]  = G->wind_v(x, y, z);
                    dcol[z]  = G->diff_h(x, y, z);
                }
                timer_stop(&G->metrics.array_copy);
                
                cbound[0] = ccol1[NZ-2];
                cbound[1] = ccol1[NZ-1];
                cbound[2] = ccol1[0];
                cbound[3] = ccol1[1];
                wbound[0] = wcol[NZ-2];
                wbound[1] = wcol[NZ-1];
                wbound[2] = wcol[0];
                wbound[3] = wcol[1];
                dbound[0] = dcol[NZ-2];
                dbound[1] = dcol[NZ-1];
                dbound[2] = dcol[0];
                dbound[3] = dcol[1];
                
                discretize(NZ, 
                           ccol1, wcol, dcol, 
                           cbound, wbound, dbound, 
                           DY, dt, ccol2);
                
                timer_start(&G->metrics.array_copy);
                for(z=0; z<NZ; z++)
                    G->conc(x, y, z, s) = ccol2[z];
                timer_stop(&G->metrics.array_copy);
            }
        }
    }
    
    timer_stop(&G->metrics.z_discret);
    
#endif
}

