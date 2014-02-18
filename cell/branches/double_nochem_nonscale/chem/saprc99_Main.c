/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Main Program File                                                */
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
/* File                 : saprc99_Main.c                            */
/* Time                 : Wed Jan  2 14:34:04 2008                  */
/* Working directory    : /home/jlinford/workspace/fixedgrid/serial/chem */
/* Equation file        : saprc99.kpp                               */
/* Output root filename : saprc99                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "saprc99_Parameters.h"
#include "saprc99_Global.h"
#include "saprc99_Sparse.h"

void Initialize();
void INTEGRATE( double TIN, double TOUT );

//double C[NSPEC];                         /* Concentration of all species */
double * C;                              /* Concentration of all species */
double * VAR;
double * FIX;
double RCONST[NREACT];                   /* Rate constants (global) */
double TIME;                             /* Current integration time */
double SUN;                              /* Sunlight intensity between [0,1] */
double TEMP;                             /* Temperature */
double RTOLS;                            /* (scalar) Relative tolerance */
double TSTART;                           /* Integration start time */
double TEND;                             /* Integration end time */
double DT;                               /* Integration step */
double ATOL[NVAR];                       /* Absolute tolerance */
double RTOL[NVAR];                       /* Relative tolerance */
double STEPMIN;                          /* Lower bound for integration step */
double STEPMAX;                          /* Upper bound for integration step */
double CFACTOR;                          /* Conversion factor for concentration units */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* MAIN - Main program - driver routine                             */
/*   Arguments :                                                    */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void init_saprc99()
{
  int i;
  
  RTOLS = 1e-3;
  for( i = 0; i < NVAR; i++ ) {
    RTOL[i] = RTOLS;
    ATOL[i] = 1.0;
  }
  STEPMIN = 0.01;
  STEPMAX = 900;
}

void init_conc(volatile double *conc)
{
  C = conc;
  VAR = &C[0];
  FIX = &C[74];
  
  Initialize();
}

void saprc99(volatile double *conc, double temp, double tstart, double tend, double dt)
{
  C = conc;
  VAR = &C[0];
  FIX = &C[74];
 
  TEMP = temp;      //236.21;
  
  TSTART = tstart;  //3600*12;
  TEND = tend;      //TSTART + 3600*24*5;
  DT = dt;          //3600.;
  
  TIME = TSTART;
  while (TIME < TEND) 
  {
    //printf("Integrating from %f to %f\n", TIME, TIME+DT);
    INTEGRATE( TIME , TIME+DT );
    TIME += DT;
    //printf("Done!\n");
  }

}
/* End of MAIN function                                             */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


