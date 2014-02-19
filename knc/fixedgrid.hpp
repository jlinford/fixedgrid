/*
 *  fixedgrid.hpp
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __FIXEDGRID_HPP__
#define __FIXEDGRID_HPP__

#include "timer.hpp"


namespace fixedgrid {

/* Floating point value type */
typedef float real_t;

/* Simple Wall clock date/time */
class SimpleDate
{
public:

	SimpleDate(int _year, int _day, int _hour, int _minute, int _second) :
		year(_year), day(_day), hour(_hour), minute(_minute), second(_second)
	{ }

	double seconds() const {
		return year*3.15569e7 + day*86400 + hour*3600 + minute*60 + second;
	}

private:

	int year;
	int day;
	int hour;
	int minute;
	int second;
};

template < typename T >
class Array2
{
public:

	typedef T element_type;

	Array2(size_t _nrows, size_t _ncols, element_type const & init) :
		nrows(_nrows), ncols(_ncols), count(_nrows*_ncols),
		data(new element_type[_nrows*_ncols]), rows(new element_type*[_nrows])
	{
		for (int i=0; i<count; ++i) {
			data[i] = init;
		}
		for (int i=0; i<nrows; ++i) {
			rows[i] = data + i*ncols;
		}
	}

	~Array2() {
		delete[] rows;
		delete[] data;
	}

	element_type * const operator[](size_t i) {
		return rows[i];
	}
	element_type const * const operator[](size_t i) const {
		return rows[i];
	}

private:

	size_t const nrows;
	size_t const ncols;
	size_t const count;
	element_type * const data;
	element_type ** const rows;
};

/* Model state variables */
class Model
{
public:

	typedef Array2<real_t> conc_t;
	typedef Array2<real_t> wind_t;
    typedef Array2<real_t> diff_t;

	Model(int _run_id, size_t _nrows, size_t _ncols,
		  size_t _dx, size_t _dy, size_t _dz,
		  real_t conc_init, real_t wind_u_init, real_t wind_v_init, real_t diff_init) :
			  write_each_iter(false), row_discret(true), col_discret(true),
		      run_id(_run_id), nrows(_nrows), ncols(_ncols),
		      dx(_dx), dy(_dy), dz(_dz), // dt(0), time(0), steps(0),
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

	void Step(SimpleDate const & tstart, SimpleDate const & tend, real_t dt);

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
