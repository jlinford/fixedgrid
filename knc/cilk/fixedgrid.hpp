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
template < size_t NROWS, size_t NCOLS >
class Model
{
public:

  typedef real_t (*matrix_t)[NCOLS];

  Model(int _run_id, size_t _dx, size_t _dy, size_t _dz,
      real_t conc_init, real_t wind_u_init, real_t wind_v_init, real_t diff_init) :
      write_each_iter(false), row_discret(true), col_discret(true),
      run_id(_run_id), dx(_dx), dy(_dy), dz(_dz), step(0), time(0),
      conc(alloc()), wind_u(alloc()), wind_v(alloc()), diff(alloc())
  { 
    conc[0:NROWS][:] = conc_init;
    wind_u[0:NROWS][:] = wind_u_init;
    wind_v[0:NROWS][:] = wind_v_init;
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

  int WriteConcToFile();

  void Step(real_t tstart, real_t tend, real_t dt);

private:

  matrix_t alloc(void) {
    void * memptr;
    std::cout << NROWS << ", " << NCOLS << std::endl;
    switch (posix_memalign(&memptr, 64, NROWS*NCOLS*sizeof(real_t))) {
      case EINVAL:
        std::cerr << "ERROR: Invalid argument to posix_memalign" << endl;
        memptr = NULL;
        break;
      case ENOMEM:
        std::cerr << "ERROR: Insufficient memory for allocation" << endl;
        memptr = NULL;
        break;
    }
    return (matrix_t)memptr;
  }

  /* Discretize rows a half timestep */
  void discretize_rows(real_t dt);

  /* Discretize cols a whole timestep */
  void discretize_cols(real_t dt);

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
    real_t const c[N], real_t const w[N], real_t const d[N],
    real_t const cb[4], real_t const wb[4], real_t const db[4], 
    real_t dcdx[N])
{
  /* Do boundary cell c[0] explicitly */
  dcdx[0] = advec_diff(cell_size,
      cb[0], wb[0], db[0],  /* 2-left neighbors */
      cb[1], wb[1], db[1],  /* 1-left neighbors */
      c[0], w[0], d[0],     /* Values */
      c[1], w[1], d[1],     /* 1-right neighbors */
      c[2], w[2], d[2]);    /* 2-right neighbors */

  /* Do boundary cell c[1] explicitly */
  dcdx[1] = advec_diff(cell_size,
      cb[1], wb[1], db[1],  /* 2-left neighbors */
      cb[2], wb[2], db[2],  /* 1-left neighbors */
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
      cb[1],  wb[1],  db[1],   /* 1-right neighbors */
      cb[2],  wb[2],  db[2]);  /* 2-right neighbors */

  /* Do boundary cell c[n-1] explicitly */
  dcdx[N-1] = advec_diff(cell_size,
      c[N-3], w[N-3], d[N-3],  /* 2-left neighbors */
      c[N-2], w[N-2], d[N-2],  /* 1-left neighbors */
      c[N-1], w[N-1], d[N-1],  /* Values */
      cb[2],  wb[2],  db[2],   /* 1-right neighbors */
      cb[3],  wb[3],  db[3]);  /* 2-right neighbors */
}

template < size_t N >
static inline void discretize(
    real_t const cell_size, real_t const dt,
    real_t const conc_in[N],
    real_t const wind[N],
    real_t const diff[N],
    real_t const concbound[4],
    real_t const windbound[4],
    real_t const diffbound[4],
    real_t conc_out[N])
{
  real_t c[N];
  real_t dcdx[N];

  conc_out[:] = conc_in[:];
  c[:] = conc_in[:];

  space_advec_diff<N>(cell_size, conc_in, wind, diff, concbound, windbound, diffbound, dcdx);
  c[:] += dt * dcdx[:];

  space_advec_diff<N>(cell_size, c, wind, diff, concbound, windbound, diffbound, dcdx);
  c[:] += dt * dcdx[:];

  conc_out[:] = 0.5 * (conc_out[:] + c[:]);
}

template < size_t NROWS, size_t NCOLS >
int Model<NROWS,NCOLS>::WriteConcToFile() 
{
  TIMER_START("File I/O");

  char buff[512];
  sprintf(buff, "fixedgrid_%03d_%05ld.dat", run_id, step);
  //int retval = conc.WriteGnuplotBinaryMatrixFile(buff);
  int retval = 0;

  TIMER_STOP("File I/O");
  return retval;
}

template < size_t NROWS, size_t NCOLS >
void Model<NROWS,NCOLS>::Step(real_t tstart, real_t tend, real_t dt)
{
  TIMER_START("Step");

  for(time=tstart; time < tend; time += dt) {
    //cout << "  Step " << step << ": Time = " << time << endl;

    discretize_rows(dt);
    discretize_cols(dt);
    discretize_rows(dt);

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
      WriteConcToFile();
    }
  }

  TIMER_STOP("Step");
}

/**
 * Discretize rows 1/2 timestep
 */
template < size_t NROWS, size_t NCOLS >
void Model<NROWS,NCOLS>::discretize_rows(real_t dt)
{
  if (row_discret) {
    TIMER_START("Row Discret");

    #pragma omp parallel
    {
      /* Buffers */
      real_t buff[NCOLS];

      /* Boundary values */
      real_t cbound[4];
      real_t wbound[4];
      real_t dbound[4];

      #pragma omp for
      for (int i = 0; i < NROWS; i++) {
        TIMER_START("Buffering");
        cbound[0] = conc[i][NCOLS - 2];
        cbound[1] = conc[i][NCOLS - 1];
        cbound[2] = conc[i][0];
        cbound[3] = conc[i][1];
        wbound[0] = wind_u[i][NCOLS - 2];
        wbound[1] = wind_u[i][NCOLS - 1];
        wbound[2] = wind_u[i][0];
        wbound[3] = wind_u[i][1];
        dbound[0] = diff[i][NCOLS - 2];
        dbound[1] = diff[i][NCOLS - 1];
        dbound[2] = diff[i][0];
        dbound[3] = diff[i][1];
        TIMER_STOP("Buffering");

        discretize<NCOLS>(dx, 0.5*dt, conc[i], wind_u[i], diff[i], cbound, wbound, dbound, buff);

        conc[i][:] = buff[:];
      }
    } // pragma omp parallel

    TIMER_STOP("Row Discret");
  }
}

/**
 * Discretize colums 1 timestep
 */
template < size_t NROWS, size_t NCOLS >
void Model<NROWS,NCOLS>::discretize_cols(real_t dt)
{
  if (col_discret) {
    TIMER_START("Col Discret");

    #pragma omp parallel
    {
      /* Buffers */
      real_t ccol[NROWS];
      real_t wcol[NROWS];
      real_t dcol[NROWS];
      real_t buff[NROWS];

      /* Boundary values */
      real_t cbound[4];
      real_t wbound[4];
      real_t dbound[4];

      #pragma omp for
      for (int j = 0; j < NCOLS; j++) {
        TIMER_START("Buffering");

        ccol[:] = conc[0:NROWS][j];
        wcol[:] = wind_v[0:NROWS][j];
        dcol[:] = diff[0:NROWS][j];

        cbound[0] = ccol[NROWS - 2];
        cbound[1] = ccol[NROWS - 1];
        cbound[2] = ccol[0];
        cbound[3] = ccol[1];
        wbound[0] = wcol[NROWS - 2];
        wbound[1] = wcol[NROWS - 1];
        wbound[2] = wcol[0];
        wbound[3] = wcol[1];
        dbound[0] = dcol[NROWS - 2];
        dbound[1] = dcol[NROWS - 1];
        dbound[2] = dcol[0];
        dbound[3] = dcol[1];
        TIMER_STOP("Buffering");

        discretize<NROWS>(dy, dt, ccol, wcol, dcol, cbound, wbound, dbound, buff);

        conc[0:NROWS][j] = buff[:];

      }
    } // pragma omp parallel

    TIMER_STOP("Col Discret");
  }
}

} // namespace fixedgrid

#endif