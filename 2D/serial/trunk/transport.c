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
 * Discretize rows 1/2 timestep 
 * @param s     Species index
 */
void discretize_all_rows(fixedgrid_t* G, real_t dt)
{
#if DO_ROW_DISCRET == 1
    
    int i, j, k;
    
    real_t buff[NCOLS];
    
    /* Boundary values */
    real_t cbound[4];
    real_t wbound[4];
    real_t dbound[4];
    
    timer_start(&G->metrics.row_discret);
    
    for(i=0; i<NROWS; i++)
    {
        for(k=0; k<NLOOKAT; k++)
        {
            cbound[0] = G->conc[k][i][NCOLS-2];
            cbound[1] = G->conc[k][i][NCOLS-1];
            cbound[2] = G->conc[k][i][0];
            cbound[3] = G->conc[k][i][1];
            wbound[0] = G->wind_u[i][NCOLS-2];
            wbound[1] = G->wind_u[i][NCOLS-1];
            wbound[2] = G->wind_u[i][0];
            wbound[3] = G->wind_u[i][1];
            dbound[0] = G->diff[i][NCOLS-2];
            dbound[1] = G->diff[i][NCOLS-1];
            dbound[2] = G->diff[i][0];
            dbound[3] = G->diff[i][1];
            
            discretize(NCOLS, G->conc[k][i], G->wind_u[i], G->diff[i], cbound, wbound, dbound, DX, dt, buff);
            
            timer_start(&G->metrics.array_copy);
            for(j=0; j<NCOLS; j++)
                G->conc[k][i][j] = buff[j];
            timer_stop(&G->metrics.array_copy);
        }
    }
    
    timer_stop(&G->metrics.row_discret);
    
#endif
}

/**
 * Discretize colums 1 timestep 
 * @param s     Species index
 */
void discretize_all_cols(fixedgrid_t* G, real_t dt)
{
#if DO_COL_DISCRET == 1
    
    int i, j, k;
    
    /* Buffers */
    real_t ccol1[NROWS];
    real_t ccol2[NROWS];
    real_t wcol[NROWS];
    real_t dcol[NROWS];
    
    /* Boundary values */
    real_t cbound[4];
    real_t wbound[4];
    real_t dbound[4];    
    
    timer_start(&G->metrics.col_discret);
    
    for(j=0; j<NCOLS; j++)
    {
        for(k=0; k<NLOOKAT; k++)
        {
            timer_start(&G->metrics.array_copy);
            for(i=0; i<NROWS; i++)
            {
                ccol1[i] = G->conc[k][i][j];
                wcol[i]  = G->wind_v[i][j];
                dcol[i]  = G->diff[i][j];
            }
            timer_stop(&G->metrics.array_copy);
            
            cbound[0] = ccol1[NROWS-2];
            cbound[1] = ccol1[NROWS-1];
            cbound[2] = ccol1[0];
            cbound[3] = ccol1[1];
            wbound[0] = wcol[NROWS-2];
            wbound[1] = wcol[NROWS-1];
            wbound[2] = wcol[0];
            wbound[3] = wcol[1];
            dbound[0] = dcol[NROWS-2];
            dbound[1] = dcol[NROWS-1];
            dbound[2] = dcol[0];
            dbound[3] = dcol[1];
            
            discretize(NROWS, ccol1, wcol, dcol, cbound, wbound, dbound, DY, dt, ccol2);
            
            timer_start(&G->metrics.array_copy);
            for(i=0; i<NROWS; i++)
                G->conc[k][i][j] = ccol2[i];
            timer_stop(&G->metrics.array_copy);
        }
    }
    
    timer_stop(&G->metrics.col_discret);
    
#endif
}

