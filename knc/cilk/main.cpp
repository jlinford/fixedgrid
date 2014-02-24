#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fixedgrid.hpp"

#define HOURS 3600
#define MINUTES 60
#define SECONDS 1

using namespace std;
using namespace fixedgrid;

/* Run ID */
static int const runID = 1;

/* Matrix dimensions */
//#define NROWS 600
//#define NCOLS 600

static real_t const nrows = 600;
static real_t const ncols = 600;

/* Cell dimensions */
static real_t const dx = 1000;
static real_t const dy = 1000;
static real_t const dz = 1000;

/* Timespan */
static real_t const tstart = 0*HOURS + 0*MINUTES + 0*SECONDS;
static real_t const tend = 12*HOURS + 0*MINUTES + 0*SECONDS;

/* Timestep size (sec) */
static real_t const dt = 50;

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
  /* Wall clock timer is always used */
  Timer wall_clock;
  wall_clock.start();

//  // Get grid dimensions
//  if (argc < 3) {
//	  cout << "Usage:" << argv[0] << " nrows ncols" << endl;
//	  exit(1);
//  }
//  size_t const nrows = atoi(argv[1]);
//  size_t const ncols = atoi(argv[2]);

  // Instantiate the template class
  typedef class Model<nrows, ncols> model_t;


  /* Initialize model */
  model_t m(runID, dx, dy, dz,
      conc_init, wind_u_init, wind_v_init, diff_init);

  /* Add O3 plume */
  m.AddPlume(4.67E+23, 300, 300);

  /* Print startup banner */
  double tspan = tend - tstart;
  cout << "\n"
       << "CONFIGURATION:\n"
       << "    ROW DISCRETIZATION:    " << m.AreRowsDiscretized() << "\n"
       << "    COLUMN DISCRETIZATION: " << m.AreColsDiscretized() << "\n"
       << "    sizeof(real_t):        " << sizeof(real_t) << "\n"
#ifdef _OPENMP
       << "    OMP_NUM_THREADS:       " << omp_get_max_threads() << "\n"
#endif
       << "\n"
       << "SPACE DOMAIN:\n"
       << "    LENGTH (X): " << nrows*dx << " meters\n"
       << "    WIDTH  (Y): " << ncols*dy << " meters\n"
       << "    DEPTH  (Z): " << dz << " meters\n"
       << "\n"
       << "TIME SPAN:\n"
       << "    " << tspan << " seconds \n"
       << "    " << (int)ceil(tspan / dt) << " timesteps of " << dt << " seconds\n"
       << endl;

  /* Store initial concentration */
  cout << "Writing initial concentration...";
  m.WriteConcToFile();
  cout << " done." << endl;

  /* Iterate */
  m.Step(tstart, tend, dt);

  /* Show final time */
  cout << "Final time: " << m.GetTime() << " seconds.\n" << endl;

  /* Store final concentration */
  cout << "Writing final concentration...";
  m.WriteConcToFile();
  cout << " done." << endl;

  /* Show metrics */
  PRINT_METRICS();

  /* Show wall clock time */
  wall_clock.stop();
  cout << "Wall clock: " << wall_clock.value()*1e-6 << endl;

  /* Cleanup and exit */
  return 0;
}
