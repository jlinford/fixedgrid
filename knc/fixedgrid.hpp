/*
 *  fixedgrid.hpp
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __FIXEDGRID_HPP__
#define __FIXEDGRID_HPP__

#include "array2.hpp"
#include "timer.hpp"

namespace fixedgrid {

/* Floating point value type */
typedef float real_t;

/* Model state variables */
class Model
{
public:

  typedef Array2<real_t> conc_t;
  typedef Array2<real_t> wind_t;
  typedef Array2<real_t> diff_t;

  Model(int _run_id, size_t _nrows, size_t _ncols, size_t _dx, size_t _dy, size_t _dz,
      real_t conc_init, real_t wind_u_init, real_t wind_v_init, real_t diff_init) :
      write_each_iter(false), row_discret(true), col_discret(true),
      run_id(_run_id), nrows(_nrows), ncols(_ncols), dx(_dx), dy(_dy), dz(_dz),
      step(0), time(0),
      conc(nrows, ncols, conc_init),
      wind_u(nrows, ncols, wind_u_init),
      wind_v(nrows, ncols, wind_v_init),
      diff(nrows, ncols, diff_init),
      metrics("fixedgrid::Model")
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

  Metrics const & GetMetrics() const {
    return metrics;
  }

  real_t GetTime() const {
    return time;
  }

  int WriteConcToFile();

  void Step(real_t tstart, real_t tend, real_t dt);

private:

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

  /* Matrix dimensions */
  size_t nrows, ncols;

  /* Cell dimensions */
  real_t dx, dy, dz;

  /* Model step */
  size_t step;

  /* Current model time */
  real_t time;

  /* Concentration field [NROWS][NCOLS] */
  conc_t conc;

  /* Wind field [NROWS][NCOLS] */
  wind_t wind_u;
  wind_t wind_v;

  /* Diffusion tensor field [NROWS][NCOLS] */
  diff_t diff;

  Metrics metrics;
};

} // namespace fixedgrid

#endif
