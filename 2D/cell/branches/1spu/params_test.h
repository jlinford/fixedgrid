#ifndef __PARAMS_H__
#define __PARAMS_H__

/* Run identifier */
#define RUN_ID  	100

/* Output directory */
#define OUTPUT_DIR  "Output"

/* Time */
#define START_YEAR  2000
#define START_DOY   100
#define START_HOUR  0
#define START_MIN   0

#define END_YEAR    2000
#define END_DOY     100
#define END_HOUR    1
#define END_MIN     00

/* Timestep size (sec) */
#define STEP_SIZE   50

/* Initial concentration (mol/cm^3 ) */
#define O3_INIT	    8.61E+09

/* Initial wind vector (m/s) */
#define WIND_U_INIT 5.0
#define WIND_V_INIT 0.0

/* Initial diffusion (m^2/s) */
#define DIFF_INIT   100.0

/* Matrix dimensions */
#define NROWS       96
#define NCOLS       170

/* Grid dimensions */
#define NX          NCOLS
#define NY          NROWS

/* Cell dimensions */
#define DX          1250.0
#define DY          1250.0
#define CELL_SIZE	1250.0

/* Emission source point coordinates */
#define SOURCE_X	50
#define SOURCE_Y	50

/* Emission source emission rate (mol/m^2/s) */
#define SOURCE_RATE	4.67E+23

#endif
