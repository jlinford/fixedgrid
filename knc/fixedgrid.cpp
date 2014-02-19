/*
 *  fixedgrid.cpp
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <cstdio>
#include "fixedgrid.hpp"

namespace fixedgrid {

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
static inline void space_advec_diff(
		int const n,
		real_t const cell_size,
		real_t const * const __restrict__ c,
		real_t const * const __restrict__ w,
		real_t const * const __restrict__ d,
		real_t const * const __restrict__ cb,
		real_t const * const __restrict__ wb,
		real_t const * const __restrict__ db,
		real_t * const __restrict__ dcdx)
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

    for(int i=2; i<n-2; i++)
    {
        dcdx[i] = advec_diff(cell_size,
                             c[i-2], w[i-2], d[i-2],  /* 2-left neighbors */
                             c[i-1], w[i-1], d[i-1],  /* 1-left neighbors */
                             c[i],   w[i],   d[i],    /* Values */
                             c[i+1], w[i+1], d[i+1],  /* 1-right neighbors */
                             c[i+2], w[i+2], d[i+2]); /* 2-right neighbors */
    }

    /* Do boundary cell c[n-2] explicitly */
    dcdx[n-2] = advec_diff(cell_size,
                           c[n-4], w[n-4], d[n-4],  /* 2-left neighbors */
                           c[n-3], w[n-3], d[n-3],  /* 1-left neighbors */
                           c[n-2], w[n-2], d[n-2],  /* Values */
                           cb[1],  wb[1],  db[1],   /* 1-right neighbors */
                           cb[2],  wb[2],  db[2]);  /* 2-right neighbors */

    /* Do boundary cell c[n-1] explicitly */
    dcdx[n-1] = advec_diff(cell_size,
                           c[n-3], w[n-3], d[n-3],  /* 2-left neighbors */
                           c[n-2], w[n-2], d[n-2],  /* 1-left neighbors */
                           c[n-1], w[n-1], d[n-1],  /* Values */
                           cb[2],  wb[2],  db[2],   /* 1-right neighbors */
                           cb[3],  wb[3],  db[3]);  /* 2-right neighbors */
}

static inline void discretize(
		int const n,
		real_t const cell_size, real_t const dt,
		real_t const * const __restrict__ conc_in,
		real_t const * const __restrict__ wind,
		real_t const * const __restrict__ diff,
		real_t const * const __restrict__ concbound,
		real_t const * const __restrict__ windbound,
		real_t const * const __restrict__ diffbound,
		real_t * const __restrict__ conc_out)
{
    real_t c[n];
    real_t dcdx[n];

    for(int i=0; i<n; i++) {
        c[i] = conc_out[i] = conc_in[i];
    }

    space_advec_diff(n, cell_size, conc_in, wind, diff, concbound, windbound, diffbound, dcdx);

    for(int i=0; i<n; i++) {
        c[i] += dt*dcdx[i];
    }

    space_advec_diff(n, cell_size, c, wind, diff, concbound, windbound, diffbound, dcdx);

    for(int i=0; i<n; i++) {
        c[i] += dt*dcdx[i];
    }

    for(int i=0; i<n; i++) {
        conc_out[i] = 0.5 * (conc_out[i] + c[i]);
    }
}


void Model::Step(SimpleDate const & start_time, SimpleDate const & end_time, real_t dt)
{
	Timer & step_timer = metrics["Step"];
	step_timer.start();

	/* Time span */
	real_t span = end_time.seconds() - start_time.seconds();

	size_t step = 0;
	for(real_t time=0; time < span; time+=dt) {
		printf("  Step %02d: Time = %07.2f sec.\n", step, time);

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
			//conc.SaveToFile();
		}
	}

	step_timer.stop();
}

/**
 * Discretize rows 1/2 timestep
 */
void Model::discretize_rows(real_t dt)
{
	if (row_discret) {
		Timer & discret_timer = metrics["Row Discret"];
		Timer & buffer_timer = metrics["Buffering"];
		discret_timer.start();

		/* Buffers */
		real_t buff[ncols];

		/* Boundary values */
		real_t cbound[4];
		real_t wbound[4];
		real_t dbound[4];

		for (int i = 0; i < nrows; i++) {
			buffer_timer.start();
			cbound[0] = conc[i][ncols - 2];
			cbound[1] = conc[i][ncols - 1];
			cbound[2] = conc[i][0];
			cbound[3] = conc[i][1];
			wbound[0] = wind_u[i][ncols - 2];
			wbound[1] = wind_u[i][ncols - 1];
			wbound[2] = wind_u[i][0];
			wbound[3] = wind_u[i][1];
			dbound[0] = diff[i][ncols - 2];
			dbound[1] = diff[i][ncols - 1];
			dbound[2] = diff[i][0];
			dbound[3] = diff[i][1];
			buffer_timer.stop();

			discretize(ncols, dx, 0.5*dt, conc[i], wind_u[i], diff[i],
					cbound, wbound, dbound, buff);

			for (int j = 0; j < ncols; j++) {
				conc[i][j] = buff[j];
			}
		}
		discret_timer.stop();
	}
}

/**
 * Discretize colums 1 timestep
 */
void Model::discretize_cols(real_t dt)
{
	if (col_discret) {
		Timer & discret_timer = metrics["Col Discret"];
		Timer & buffer_timer = metrics["Buffering"];
		discret_timer.start();

		/* Buffers */
		real_t ccol[nrows];
		real_t wcol[nrows];
		real_t dcol[nrows];
		real_t buff[nrows];

		/* Boundary values */
		real_t cbound[4];
		real_t wbound[4];
		real_t dbound[4];

		for(int j=0; j<ncols; j++) {
			buffer_timer.start();
			for(int i=0; i<nrows; i++) {
				ccol[i] = conc[i][j];
				wcol[i]  = wind_v[i][j];
				dcol[i]  = diff[i][j];
			}

			cbound[0] = ccol[nrows-2];
			cbound[1] = ccol[nrows-1];
			cbound[2] = ccol[0];
			cbound[3] = ccol[1];
			wbound[0] = wcol[nrows-2];
			wbound[1] = wcol[nrows-1];
			wbound[2] = wcol[0];
			wbound[3] = wcol[1];
			dbound[0] = dcol[nrows-2];
			dbound[1] = dcol[nrows-1];
			dbound[2] = dcol[0];
			dbound[3] = dcol[1];
			buffer_timer.stop();

			discretize(nrows, dy, dt, ccol, wcol, dcol, cbound, wbound, dbound, buff);

			for(int i=0; i<nrows; i++) {
				conc[i][j] = buff[i];
			}
		}
		discret_timer.stop();
	}
}

} // namespace fixedgrid
