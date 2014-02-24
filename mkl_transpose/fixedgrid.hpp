/*
 *  fixedgrid.hpp
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __FIXEDGRID_HPP__
#define __FIXEDGRID_HPP__

#include <cstdlib>
#include <cerrno>
#include <mkl.h>
#include "timer.hpp"


namespace fixedgrid {

/* Floating point value type */
typedef float real_t;

/* Model state variables */
template < size_t NROWS, size_t NCOLS >
class Model
{
public:

  typedef real_t (*matrix_t)[NCOLS];
  typedef real_t (*transposed_t)[NROWS];

  Model(int _run_id, size_t _dx, size_t _dy, size_t _dz,
      real_t conc_init, real_t wind_u_init, real_t wind_v_init, real_t diff_init) :
      write_each_iter(false), row_discret(true), col_discret(true),
      run_id(_run_id), dx(_dx), dy(_dy), dz(_dz), step(0), time(0),
      conc((matrix_t)alloc()),
      wind_u((matrix_t)alloc()),
      wind_v((transposed_t)alloc()),
      diff((matrix_t)alloc())
  { 
    conc[0:NROWS][:] = conc_init;
    wind_u[0:NROWS][:] = wind_u_init;
    wind_v[0:NCOLS][:] = wind_v_init;
    diff[0:NROWS][:] = diff_init;
  }

  ~Model()
  { }

  /* Add emission plume (mol/m^3) */
  void AddPlume(real_t plume, size_t row, size_t col) {
    conc[row][col] += plume / (dx * dy * dz);
  }

  bool AreRowsDiscretized() const {
    return row_discret;
  }

  bool AreColsDiscretized() const {
    return col_discret;
  }

  real_t GetTime() const {
    return time;
  }

  void Step(real_t tstart, real_t tend, real_t dt);

  void WriteConcFile(void);

private:

  void * alloc() {
    void * memptr;
    switch (posix_memalign(&memptr, 64, NROWS*NCOLS*sizeof(real_t))) {
      case EINVAL:
        std::cerr << "ERROR: Invalid argument to posix_memalign" << std::endl;
        memptr = NULL;
        break;
      case ENOMEM:
        std::cerr << "ERROR: Insufficient memory for allocation" << std::endl;
        memptr = NULL;
        break;
    }
    return memptr;
  }

  int WriteGnuplotBinaryMatrixFile(matrix_t mat, char const * fname);

  /* Write output each iteration? */
  bool write_each_iter;

  /* Discretize rows each iteration? */
  bool row_discret;

  /* Discretize columns each iteration? */
  bool col_discret;

  /* Run identifier */
  int run_id;

  /* Cell dimensions */
  real_t dx, dy, dz;

  /* Model step */
  size_t step;

  /* Current model time */
  real_t time;

  /* Concentration field [NROWS][NCOLS] */
  matrix_t conc;

  /* Wind field [NROWS][NCOLS] */
  matrix_t wind_u;
  transposed_t wind_v;

  /* Diffusion tensor field [NROWS][NCOLS] */
  matrix_t diff;
};


/*
 * Upwinded advection/diffusion stencil.
 * c = conc, w = wind, d = diff
 * x2l is the 2-left neighbor of x, etc.
 * x2r is the 2-right neighbor of x, etc.
 */
static inline real_t advec_diff(
		real_t const cell_size,
		real_t const c2l, real_t const w2l, real_t const d2l,
		real_t const c1l, real_t const w1l, real_t const d1l,
		real_t const   c, real_t const   w, real_t const   d,
		real_t const c1r, real_t const w1r, real_t const d1r,
		real_t const c2r, real_t const w2r, real_t const d2r)
{
  real_t windL = (w1l + w) / 2.0;
  real_t advec_termL;
  if (windL >= 0.0) {
    advec_termL = (1.0 / 6.0) * (-c2l + 5.0 * c1l + 2.0 * c);
  } else {
    advec_termL = (1.0 / 6.0) * (2.0 * c1l + 5.0 * c - c1r);
  }
  advec_termL *= windL;

  real_t windR = (w1r + w) / 2.0;
  real_t advec_termR;
  if (windR >= 0.0) {
    advec_termR = (1.0 / 6.0) * (-c1l + 5.0 * c + 2.0 * c1r);
  } else {
    advec_termR = (1.0 / 6.0) * (2.0 * c + 5.0 * c1r - c2r);
  }
  advec_termR *= windR;

  real_t advec_term = (advec_termL - advec_termR) / cell_size;
  real_t diff_term = (((d1l + d) / 2) * (c1l - c) - ((d + d1r) / 2) * (c - c1r)) / (cell_size * cell_size);

  return advec_term + diff_term;
}


/*
 * Applies advection / diffusion stencil
 */
template < size_t N >
static inline void space_advec_diff(
    real_t const cell_size,
	real_t const cb0, real_t const wb0, real_t const db0,
	real_t const cb1, real_t const wb1, real_t const db1,
	real_t const cb2, real_t const wb2, real_t const db2,
	real_t const cb3, real_t const wb3, real_t const db3,
	real_t const c[N], real_t const w[N], real_t const d[N],
	real_t dcdx[N])
{
  /* Do boundary cell c[0] explicitly */
  dcdx[0] = advec_diff(cell_size,
      cb0, wb0, db0,  /* 2-left neighbors */
      cb1, wb1, db1,  /* 1-left neighbors */
      c[0], w[0], d[0],     /* Values */
      c[1], w[1], d[1],     /* 1-right neighbors */
      c[2], w[2], d[2]);    /* 2-right neighbors */

  /* Do boundary cell c[1] explicitly */
  dcdx[1] = advec_diff(cell_size,
      cb1, wb1, db1,  /* 2-left neighbors */
      cb2, wb2, db2,  /* 1-left neighbors */
      c[1], w[1], d[1],     /* Values */
      c[2], w[2], d[2],     /* 1-right neighbors */
      c[3], w[3], d[3]);    /* 2-right neighbors */

  for(int i=2; i<N-2; ++i) {
      dcdx[i] = advec_diff(cell_size,
          c[i-2], w[i-2], d[i-2],  /* 2-left neighbors */
          c[i-1], w[i-1], d[i-1],  /* 1-left neighbors */
          c[i],   w[i],   d[i],    /* Values */
          c[i+1], w[i+1], d[i+1],  /* 1-right neighbors */
          c[i+2], w[i+2], d[i+2]); /* 2-right neighbors */
  }

  /* Do boundary cell c[n-2] explicitly */
  dcdx[N-2] = advec_diff(cell_size,
      c[N-4], w[N-4], d[N-4],  /* 2-left neighbors */
      c[N-3], w[N-3], d[N-3],  /* 1-left neighbors */
      c[N-2], w[N-2], d[N-2],  /* Values */
      cb1,  wb1,  db1,   /* 1-right neighbors */
      cb2,  wb2,  db2);  /* 2-right neighbors */

  /* Do boundary cell c[n-1] explicitly */
  dcdx[N-1] = advec_diff(cell_size,
      c[N-3], w[N-3], d[N-3],  /* 2-left neighbors */
      c[N-2], w[N-2], d[N-2],  /* 1-left neighbors */
      c[N-1], w[N-1], d[N-1],  /* Values */
      cb2,  wb2,  db2,   /* 1-right neighbors */
      cb3,  wb3,  db3);  /* 2-right neighbors */
}


template < size_t N >
static inline void discretize(
    real_t const cell_size, real_t const dt,
	real_t const cb0, real_t const wb0, real_t const db0,
	real_t const cb1, real_t const wb1, real_t const db1,
	real_t const cb2, real_t const wb2, real_t const db2,
	real_t const cb3, real_t const wb3, real_t const db3,
    real_t const conc_in[N], real_t const wind[N], real_t const diff[N],
    real_t conc_out[N])
{
  real_t c[N];
  real_t dcdx[N];

  conc_out[:] = conc_in[:];
  c[:] = conc_in[:];

  space_advec_diff<N>(
		  cell_size,
		  cb0, wb0, db0,
		  cb1, wb1, db1,
		  cb2, wb2, db2,
		  cb3, wb3, db3,
		  conc_in, wind, diff, dcdx);
  c[:] += dt * dcdx[:];

  space_advec_diff<N>(
		  cell_size,
		  cb0, wb0, db0,
		  cb1, wb1, db1,
		  cb2, wb2, db2,
		  cb3, wb3, db3,
		  c, wind, diff, dcdx);
  c[:] += dt * dcdx[:];

  conc_out[:] = 0.5 * (conc_out[:] + c[:]);
}


template < size_t NROWS, size_t NCOLS >
void discretize(real_t const dx, real_t const dt,
		real_t (*conc)[NCOLS], real_t (*wind)[NCOLS], real_t (*diff)[NCOLS])
{
	TIMER_START("Discret");

	#pragma omp parallel
	{
	  /* Buffers */
	  real_t buff[NCOLS];

	  #pragma omp for
	  for (int i = 0; i < NROWS; i++) {
		real_t cb0 = conc[i][NCOLS - 2];
		real_t cb1 = conc[i][NCOLS - 1];
		real_t cb2 = conc[i][0];
		real_t cb3 = conc[i][1];
		real_t wb0 = wind[i][NCOLS - 2];
		real_t wb1 = wind[i][NCOLS - 1];
		real_t wb2 = wind[i][0];
		real_t wb3 = wind[i][1];
		real_t db0 = diff[i][NCOLS - 2];
		real_t db1 = diff[i][NCOLS - 1];
		real_t db2 = diff[i][0];
		real_t db3 = diff[i][1];

		discretize<NCOLS>(dx, dt,
				cb0, wb0, db0,
				cb1, wb1, db1,
				cb2, wb2, db2,
				cb3, wb3, db3,
				conc[i], wind[i], diff[i],
				buff);

		conc[i][:] = buff[:];
	  }
	} // pragma omp parallel

	TIMER_STOP("Discret");
}

/*
  * From gnuplot/doc/html/gnuplot.html:
  * For binary input data, the first element of the first row must contain
  * the number of data columns. (This number is ignored for ascii input).
  * Both the coordinates and the data values in a binary input are treated
  * as single precision floats. Example commands for plotting non-uniform matrix data:
  *
  *     splot 'file' nonuniform matrix using 1:2:3  # ascii input
  *     splot 'file' binary matrix using 1:2:3      # binary input
  *
  * Thus the data organization for non-uniform matrix input is
  *
  *      <N+1>  <y0>   <y1>   <y2> ...  <yN>
  *      <x0> <z0,0> <z0,1> <z0,2> ... <z0,N>
  *      <x1> <z1,0> <z1,1> <z1,2> ... <z1,N>
  *       :      :      :      :   ...    :
  *
  * which is then converted into triplets:
  *
  *      <x0> <y0> <z0,0>
  *      <x0> <y1> <z0,1>
  *      <x0> <y2> <z0,2>
  *       :    :     :
  *      <x0> <yN> <z0,N>
  *      <x1> <y0> <z1,0>
  *      <x1> <y1> <z1,1>
  *       :    :     :
  *
  * These triplets are then converted into gnuplot iso-curves and then gnuplot
  * proceeds in the usual manner to do the rest of the plotting.
  */
template < size_t NROWS, size_t NCOLS >
int Model<NROWS,NCOLS>::WriteGnuplotBinaryMatrixFile(matrix_t mat, char const * fname)
{
  using namespace std;

  // Open file for writing
  FILE * fout = fopen(fname, "wb");
  if (!fout) {
    int err = errno;
    fprintf(stderr, "ERROR: Can't open '%s' for writing\n", fname);
    return err;
  }

  size_t nmemb = NCOLS + 1;
  float buffer[nmemb];

  // Write ncols and column indices as float values
  buffer[0] = (float)NCOLS;
  for (int j=0; j<NCOLS; ++j) {
    buffer[j+1] = (float)j;
  }
  if (fwrite(buffer, sizeof(float), nmemb, fout) != nmemb) {
    fprintf(stderr, "ERROR: Failed to write %ld elements to '%s'\n", nmemb, fname);
    return -1;
  }

  // Write matrix rows
  for (int i=0; i<NROWS; ++i) {
    buffer[0] = (float)i;
    buffer[1:NCOLS] = mat[i][:];
    if (fwrite(buffer, sizeof(float), nmemb, fout) != nmemb) {
      fprintf(stderr, "ERROR: Failed to write %ld elements to '%s'\n", nmemb, fname);
      return -1;
    }
  }

  // Cleanup and return
  fclose(fout);
  return 0;
}

template < size_t NROWS, size_t NCOLS >
void Model<NROWS,NCOLS>::WriteConcFile(void)
{
  TIMER_START("File I/O");

  // Build filename
  char fname[512];
  sprintf(fname, "fixedgrid_%03d_%05ld.dat", run_id, step);

  if (WriteGnuplotBinaryMatrixFile(conc, fname)) {
	  std::cerr << "ERROR: Failed to write '" << fname << "'" << std::endl;
  } else {
	  std::cout << "Wrote '" << fname << "'" << std::endl;
  }

  TIMER_STOP("File I/O");
}

template < size_t NROWS, size_t NCOLS >
void Model<NROWS,NCOLS>::Step(real_t tstart, real_t tend, real_t dt)
{
  TIMER_START("Step");

  for(time=tstart; time < tend; time += dt) {
	  //cout << "  Step " << step << ": Time = " << time << endl;

	  if (row_discret) {
		  discretize<NROWS, NCOLS>(dx, 0.5*dt, conc, wind_u, diff);
	  }

	  if (col_discret) {
		  // Transpose
		  TIMER_START("Transpose");
		  mkl_simatcopy('R', 'T', NROWS, NCOLS, 1, conc[0], NCOLS, NROWS);
		  mkl_simatcopy('R', 'T', NROWS, NCOLS, 1, diff[0], NCOLS, NROWS);
		  TIMER_STOP("Transpose");

		  discretize<NCOLS, NROWS>(dy, dt, (transposed_t)conc, wind_v, (transposed_t)diff);

		  // Transpose
		  TIMER_START("Transpose");
		  mkl_simatcopy('R', 'T', NCOLS, NROWS, 1, conc[0], NROWS, NCOLS);
		  mkl_simatcopy('R', 'T', NCOLS, NROWS, 1, diff[0], NROWS, NCOLS);
		  TIMER_STOP("Transpose");
	  }

	  if (row_discret) {
		  discretize<NROWS, NCOLS>(dx, 0.5*dt, conc, wind_u, diff);
	  }

    ++step;

    /*
     * Could update wind field here...
     */

    /*
     * Could update diffusion tensor here...
     */

    /*
     * Could update environment here...
     */

    /* Store concentration */
    if (write_each_iter) {
      WriteConcFile();
    }
  }

  TIMER_STOP("Step");
}


} // namespace fixedgrid

#endif
