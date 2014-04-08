#ifndef __FILEIO_H__
#define __FILEIO_H__

#include "params.h"

void write_conc(double conc[NROWS][NCOLS][NSPEC], int spec, char* sname, long iter, int proc);

#endif
