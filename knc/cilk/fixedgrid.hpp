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
#include "timer.hpp"


namespace fixedgrid {

/* Floating point value type */
typedef float real_t;

/* Model state variables */
template < size_t M, size_t N >
class Model
{

public:

  enum { NROWS=M, NCOLS=N, NCOLS_ALIGNED=((NCOLS + 63) & ~63UL) };

  typedef real_t (*matrix_t)[NCOLS_ALIGNED];

  Model(int _run_id, size_t _dx, size_t _dy, size_t _dz,
      real_t conc_init, real_t wind_u_init, real_t wind_v_init, real_t diff_init) :
      write_each_iter(false), row_discret(true), col_discret(true),
      run_id(_run_id), dx(_dx), dy(_dy), dz(_dz), step(0), time(0),
      conc(init_matrix(conc_init)),
      wind_u(init_matrix(wind_u_init)),
      wind_v(init_matrix(wind_v_init)),
      diff(init_matrix(diff_init))
  { }

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

  matrix_t init_matrix(real_t const init) {
    void * memptr;
    switch (posix_memalign(&memptr, 64, NROWS*NCOLS_ALIGNED*sizeof(real_t))) {
      case EINVAL:
        std::cerr << "ERROR: Invalid argument to posix_memalign" << std::endl;
        return NULL;
      case ENOMEM:
        std::cerr << "ERROR: Insufficient memory for allocation" << std::endl;
        return NULL;
    }

    real_t * __restrict__ p = (real_t*)memptr;
    __assume_aligned(p, 64);
    p[0:NROWS*NCOLS_ALIGNED] = init;

    return (matrix_t)memptr;
  }

  int WriteGnuplotBinaryMatrixFile(matrix_t mat, char const * fname);

  /* Discretize rows a half timestep */
  void discretize_rows();

  /* Discretize cols a whole timestep */
  void discretize_cols();

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

  /* Time step size */
  real_t dt;

  /* Concentration field [NROWS][NCOLS] */
  matrix_t conc;

  /* Wind field [NROWS][NCOLS] */
  matrix_t wind_u;
  matrix_t wind_v;

  /* Diffusion tensor field [NROWS][NCOLS] */
  matrix_t diff;
};


/*
 * Upwinded advection/diffusion stencil.
 * c = conc, w = wind, d = diff
 * x2l is the 2-left neighbor of x, etc.
 * x2r is the 2-right neighbor of x, etc.
 */
static inline __declspec(vector)
real_t advec_diff(
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
static inline void space_advec_diff(
    size_t const N, real_t const cell_size,
    real_t const cb0, real_t const wb0, real_t const db0,
    real_t const cb1, real_t const wb1, real_t const db1,
    real_t const cb2, real_t const wb2, real_t const db2,
    real_t const cb3, real_t const wb3, real_t const db3,
    real_t const * __restrict__ c,
    real_t const * __restrict__ w,
    real_t const * __restrict__ d,
    real_t * __restrict__ dcdx)
{
  __assume_aligned(c, 64);
  __assume_aligned(w, 64);
  __assume_aligned(d, 64);
  __assume_aligned(dcdx, 64);
  __assume(N%16==0);

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


static inline void discretize(
    size_t const N, real_t const cell_size, real_t const dt,
    real_t const * __restrict__ c,
    real_t const * __restrict__ w,
    real_t const * __restrict__ d,
    real_t * __restrict co)
{
  __assume_aligned(c, 64);
  __assume_aligned(w, 64);
  __assume_aligned(d, 64);
  __assume_aligned(co, 64);
  __assume(N%16==0);

  real_t buff[N] __attribute__((aligned(64)));
  real_t dcdx[N] __attribute__((aligned(64)));

  real_t const cb0 = c[N-2];
  real_t const cb1 = c[N-1];
  real_t const cb2 = c[0];
  real_t const cb3 = c[1];
  real_t const wb0 = w[N-2];
  real_t const wb1 = w[N-1];
  real_t const wb2 = w[0];
  real_t const wb3 = w[1];
  real_t const db0 = d[N-2];
  real_t const db1 = d[N-1];
  real_t const db2 = d[0];
  real_t const db3 = d[1];

  co[0:N] = c[0:N];
  buff[0:N] = c[0:N];

  space_advec_diff(
      N, cell_size,
      cb0, wb0, db0,
      cb1, wb1, db1,
      cb2, wb2, db2,
      cb3, wb3, db3,
      c, w, d, dcdx);
  buff[0:N] += dt * dcdx[0:N];

  space_advec_diff(
      N, cell_size,
      cb0, wb0, db0,
      cb1, wb1, db1,
      cb2, wb2, db2,
      cb3, wb3, db3,
      buff, w, d, dcdx);
  buff[0:N] += dt * dcdx[0:N];

  co[0:N] = 0.5 * (co[0:N] + buff[0:N]);
}


/**
 * Discretize rows 1/2 timestep
 */
template < size_t M, size_t N >
void Model<M,N>::discretize_rows()
{
  if (row_discret) {
    TIMER_START("Row Discret");

    /* Buffers */
    real_t buff[NCOLS_ALIGNED] __attribute__((aligned(64)));

    #pragma omp for private(buff)
    for (int i = 0; i < NROWS; i++) {
      real_t * __restrict__ c = conc[i];
      real_t const * __restrict__ w = wind_u[i];
      real_t const * __restrict__ d = diff[i];

      __assume_aligned(c, 64);
      __assume_aligned(w, 64);
      __assume_aligned(d, 64);

      discretize(NCOLS_ALIGNED, dx, 0.5*dt, c, w, d, buff);
      c[0:NCOLS_ALIGNED] = buff[:];
    }

    TIMER_STOP("Row Discret");
  }
}

/**
 * Discretize colums 1 timestep
 */
template < size_t M, size_t N >
void Model<M,N>::discretize_cols()
{
  if (col_discret) {
    TIMER_START("Col Discret");

    /* Buffers */
    real_t ccol[NROWS] __attribute__((aligned(64)));
    real_t wcol[NROWS] __attribute__((aligned(64)));
    real_t dcol[NROWS] __attribute__((aligned(64)));
    real_t buff[NROWS] __attribute__((aligned(64)));

    #pragma omp for private(ccol, wcol, dcol, buff)
    for (int j = 0; j < NCOLS; j++) {
      ccol[:] = conc[0:NROWS][j];
      wcol[:] = wind_v[0:NROWS][j];
      dcol[:] = diff[0:NROWS][j];

      discretize(NCOLS, dy, dt, ccol, wcol, dcol, buff);
      conc[0:NROWS][j] = buff[:];
    }

    TIMER_STOP("Col Discret");
  }
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
template < size_t M, size_t N >
int Model<M,N>::WriteGnuplotBinaryMatrixFile(matrix_t mat, char const * fname)
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
    buffer[1:NCOLS] = mat[i][0:NCOLS];
    if (fwrite(buffer, sizeof(float), nmemb, fout) != nmemb) {
      fprintf(stderr, "ERROR: Failed to write %ld elements to '%s'\n", nmemb, fname);
      return -1;
    }
  }

  // Cleanup and return
  fclose(fout);
  return 0;
}

template < size_t M, size_t N >
void Model<M,N>::WriteConcFile(void)
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

template < size_t M, size_t N >
void Model<M,N>::Step(real_t tstart, real_t tend, real_t dt)
{
  TIMER_START("Step");

  this->dt = dt;
  for(time=tstart; time < tend; time += dt) {

    #pragma omp parallel
    {
      discretize_rows();
      discretize_cols();
      discretize_rows();
    }

    // Show progress
    std::cout << "  Step " << step << ": Time = " << time << std::endl;
    if ((step % 10) == 0) {
      PRINT_METRICS();
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
