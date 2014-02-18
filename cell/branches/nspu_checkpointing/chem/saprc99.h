/*
 *  saprc99.h
 *  
 *
 *  Created by John Linford on 4/11/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __SAPRC99_H__
#define __SAPRC99_H__

/**************************************************
 * Includes                                       *
 **************************************************/

#include "saprc99_Monitor.h"

/**************************************************
 * Function Prototypes                            *
 **************************************************/

void init_saprc99();
void init_conc(volatile double *conc);
void saprc99(volatile double *conc, double temp, double tstart, double tend, double dt);

#endif
