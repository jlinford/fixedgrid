#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "fixedgrid.hpp"

using namespace std;
using namespace fixedgrid;

/* Run ID */
static int const runID = 1;

/* Matrix dimensions */
static size_t const nrows = 600;
static size_t const ncols = 600;

/* Cell dimensions */
static real_t const dx = 1000;
static real_t const dy = 1000;
static real_t const dz = 1000;

/* Timespan */
static SimpleDate const start_time(2000, 100, 0, 0, 0);
static SimpleDate const end_time(2000, 100, 1, 0, 0);

/* Timestep size (sec) */
static real_t const step_size = 50;

/* Initial O3 concentration */
static real_t const conc_init = 8.61E+09;

/* Initial wind vector (m/s) */
static real_t const wind_u_init = 5.0;
static real_t const wind_v_init = -5.0;

/* Initial diffusion (m^2/s) */
static real_t const diff_init = 100.0;



/**
 * Program entry point.
 */
int main(int argc, char** argv)
{
    /* Initialize model */
    Model m(runID, nrows, ncols, dx, dy, dz,
    		conc_init, wind_u_init, wind_v_init, diff_init);
    
    /* Add O3 plume */
    m.AddPlume(4.67E+23, 300, 300);

    /* Print startup banner */
    double span_seconds = end_time.seconds() - start_time.seconds();
    cout << "\n"
		 << "CONFIGURATION:\n"
		 << "    ROW DISCRETIZATION:    " << m.AreRowsDiscretized() << "\n"
		 << "    COLUMN DISCRETIZATION: " << m.AreColsDiscretized() << "\n"
		 << "    sizeof(real_t):        " << sizeof(real_t) << "\n"
		 << "\n"
		 << "SPACE DOMAIN:\n"
		 << "    LENGTH (X): " << ncols*dx << " meters\n"
		 << "    WIDTH  (Y): " << nrows*dy << " meters\n"
		 << "    DEPTH  (Z): " << dz << "meters\n"
		 << "\n"
		 << "TIME SPAN:\n"
		 << "    " << span_seconds << " seconds \n"
		 << "    " << (int)ceil(span_seconds / step_size) << " timesteps of " << step_size << " seconds\n"
		 << "\n";
    
    /* Store initial concentration */
    cout << "Writing initial concentration...";
    //m.conc.SaveToFile();
    cout << " done." << endl;
    
    /* Iterate */
    m.Step(start_time, end_time, step_size);

    /* Store final concentration */
    cout << "Writing final concentration...";
    //m.conc.SaveToFile();
    cout << " done." << endl;

    /* Show final time */
    //cout << "Final time: " << m.GetModelTime() << " seconds.\n" << endl;
    
    /* Show metrics */
    cout << m.GetMetrics();

    /* Cleanup and exit */
    return 0;
}
