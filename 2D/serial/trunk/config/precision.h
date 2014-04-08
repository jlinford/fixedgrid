/*
 *  precision.h
 *  cell_double
 *
 *  Created by John Linford on 4/29/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __PRECISION_H__
#define __PRECISION_H__

#include "params.h"

#if DOUBLE_PRECISION == 1

    #define real_t double

#else

    #define real_t float

#endif

#endif
