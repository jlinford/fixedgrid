/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Initialization File                                              */
/*                                                                  */
/* Generated by KPP-2.2 symbolic chemistry Kinetics PreProcessor    */
/*       (http://www.cs.vt.edu/~asandu/Software/KPP)                */
/* KPP is distributed under GPL, the general public licence         */
/*       (http://www.gnu.org/copyleft/gpl.html)                     */
/* (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           */
/* (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            */
/*     With important contributions from:                           */
/*        M. Damian, Villanova University, USA                      */
/*        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany */
/*                                                                  */
/* File                 : saprc99_Initialize.c                      */
/* Time                 : Wed Jun 11 10:09:48 2008                  */
/* Working directory    : /Users/jlinford/workspace/saprc99         */
/* Equation file        : saprc99.kpp                               */
/* Output root filename : saprc99                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "saprc99_Parameters.h"
#include "params.h"
#include "timer.h"


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Initialize - function to initialize concentrations               */
/*   Arguments :                                                    */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void saprc99_Initialize(real_t C[NSPEC])
{
    int i;
    real_t x;
    
    real_t * VAR = &C[0];
    real_t * FIX = &C[74];
        
    /* Initialize concentrations */
    x = (0.0e0)*CFACTOR;
    for( i = 0; i < NVAR; i++ )
        VAR[i] = x;

    x = (0.0e0)*CFACTOR;
    for( i = 0; i < NFIX; i++ )
        FIX[i] = x;
    
    VAR[1] = (6.77e-4)*CFACTOR;
    VAR[2] = (1.16e-3)*CFACTOR;
    VAR[3] = (3.92e-4)*CFACTOR;
    VAR[7] = (0.2E0)*CFACTOR;
    VAR[8] = (5.e-2)*CFACTOR;
    VAR[10] = (1.167e-2)*CFACTOR;
    VAR[18] = (1.e-3)*CFACTOR;
    VAR[19] = (1.88e-2)*CFACTOR;
    VAR[20] = (4.69e-2)*CFACTOR;
    VAR[22] = (3.06e-2)*CFACTOR;
    VAR[23] = (8.74e-3)*CFACTOR;
    VAR[28] = (5.89e-3)*CFACTOR;
    VAR[29] = (4.17e-2)*CFACTOR;
    VAR[30] = (1.18e-2)*CFACTOR;
    VAR[33] = (5.60e-4)*CFACTOR;
    VAR[37] = (7.51e-5)*CFACTOR;
    VAR[38] = (6.06e-4)*CFACTOR;
    VAR[40] = (8.37e-5)*CFACTOR;
    VAR[41] = (5.07e-3)*CFACTOR;
    VAR[43] = (1.89e-2)*CFACTOR;
    VAR[44] = (1.21e-4)*CFACTOR;
    VAR[46] = (4.33e-4)*CFACTOR;
    VAR[48] = (8.20e-4)*CFACTOR;
    VAR[49] = (1.30e-3)*CFACTOR;
    VAR[50] = (1.04e-2)*CFACTOR;
    VAR[51] = (8.93e-5)*CFACTOR;
    VAR[52] = (7.97e-3)*CFACTOR;
    VAR[54] = (2.316e-3)*CFACTOR;
    VAR[55] = (1.121e-2)*CFACTOR;
    VAR[57] = (7.843e-9)*CFACTOR;
    VAR[58] = (1.72e-3)*CFACTOR;
    VAR[59] = (3.26e-3)*CFACTOR;
    VAR[60] = (1.93e-3)*CFACTOR;
    VAR[61] = O3_INIT;
    VAR[68] = (1.0e-1)*CFACTOR;
    VAR[69] = (5.0e-2)*CFACTOR;
    FIX[0] = (1.0e+6)*CFACTOR;
    FIX[1] = (2.09e+5)*CFACTOR;
    FIX[2] = (2.0e+04)*CFACTOR;
    FIX[4] = (1.0e0)*CFACTOR;
}

/* End of Initialize function                                       */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


