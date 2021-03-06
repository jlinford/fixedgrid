/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Global Data Header File                                          */
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
/* File                 : saprc99_Global.h                          */
/* Time                 : Wed Jun 11 10:09:48 2008                  */
/* Working directory    : /Users/jlinford/workspace/saprc99         */
/* Equation file        : saprc99.kpp                               */
/* Output root filename : saprc99                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef __SAPRC99_GLOBAL_H__
#define __SAPRC99_GLOBAL_H__

#include "saprc99_Parameters.h"

/* Declaration of global variables */

extern volatile double * C;                      /* Concentration of all species */
extern volatile double * VAR;
extern volatile double * FIX;
extern double RCONST[NREACT];                    /* Rate constants (global) */
extern double ATOL[NVAR];                        /* Absolute tolerance */
extern double RTOL[NVAR];                        /* Relative tolerance */
extern double TIME;                              /* Current integration time */
extern double DT;                                /* Integration step */
extern double SUN;                               /* Sunlight intensity between [0,1] */
extern double TEMP;                              /* Temperature */
extern double STEPMIN;                           /* Lower bound for integration step */
/* Now a preproc define */
//extern double CFACTOR;                           /* Conversion factor for concentration units */

#endif
