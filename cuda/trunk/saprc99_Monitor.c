/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Utility Data Initialization                                      */
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
/* File                 : saprc99_Monitor.c                         */
/* Time                 : Wed Jun 11 10:09:48 2008                  */
/* Working directory    : /Users/jlinford/workspace/saprc99         */
/* Equation file        : saprc99.kpp                               */
/* Output root filename : saprc99                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "params.h"

#if DO_CHEMISTRY == 1

/* Names of chemical species */
char *  SPC_NAMES[] = {
    "H2SO4","HCOOH","CCO_OH","RCO_OH","CCO_OOH","RCO_OOH","XN","XC",
    "SO2","O1D","ALK1","BACL","PAN","PAN2","PBZN","MA_PAN",
    "H2O2","N2O5","HONO","ALK2","ALK3","TBU_O","ALK5","ARO2",
    "HNO4","COOH","HOCOO","BZNO2_O","MEOH","ALK4","ARO1","DCB2",
    "DCB3","CRES","DCB1","NPHE","ROOH","BALD","PHEN","CO",
    "MGLY","ACET","HNO3","ETHENE","GLY","BZ_O","ISOPRENE","R2O2",
    "TERP","METHACRO","OLE1","ISOPROD","OLE2","MVK","CCHO","HCHO",
    "RNO3","O3P","RCHO","MEK","PROD2","O3","HO2","RO2_N",
    "MA_RCO3","C_O2","BZCO_O2","RO2_R","NO","NO2","NO3","CCO_O2",
    "RCO_O2","OH","AIR","O2","H2O","H2","CH4" 
}; 

/* Indexes of species to look at */
int  LOOKAT[NLOOKAT] = {
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
    24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
    36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
    60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
    72, 73, 74, 75, 76, 77, 78 
};

/* Indexes of species to monitor */
int  MONITOR[NMONITOR] = {
    43, 61, 68, 69 
};

#else

char * SPC_NAMES[] = {
    "O3"
};

/* Indexes of species to look at */
int  LOOKAT[NLOOKAT] = {
    ind_O3 
};

/* Indexes of species to monitor */
int  MONITOR[NMONITOR] = {
    ind_O3 
};

#endif


