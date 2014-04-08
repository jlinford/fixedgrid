#ifndef __FILEIO_H__
#define __FILEIO_H__

#include "params.h"

void write_conc(volatile double conc[NROWS][NCOLS], long iter, int proc);

#endif
