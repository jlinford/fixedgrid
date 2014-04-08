/*
 *  transport.h
 *  fixedgrid
 *
 *  Created by John Linford on 6/17/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __TRANSPORT_H__
#define __TRANSPORT_H__

/**************************************************
 * Includes                                       *
 **************************************************/

#include <stdint.h>

/**************************************************
 * Function prototypes                            *
 **************************************************/

void discretize_row_block(uint64_t argvp);

void discretize_col_block(uint64_t argvp);

#endif

