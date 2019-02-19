/*
 * Copyright 2019 Holger Hoffmann, Ted Moldenhawer, Thomas Münch, Jürgen Schmidt, Michael Seiler, Frank Spahn
 * 
 * This file is part of DDX.
 * 
 * DDX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * DDX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with DDX.  If not, see <https://www.gnu.org/licenses/>.
 * 
 */


#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gsl/gsl_rng.h"
#include "psvindex.h"
#include "impact.h"
#include "moons.h"
#include "initialvalues.h"

/******************************************************************************
/NUMBERS **********************************************************************
******************************************************************************/

#define val_pi (3.14159265358979323846264338328)
#define val_2pi (6.28318530717958647692528676656)
#define val_pi2 (1.57079632679489661923132169164)
#define val_sqrtpi (1.77245385090551602729816748334)
#define val_euler (2.71828182845904523536028747135)

#define val_dim1 (135)
#define val_dim2 (120)
#define val_dim3 (124)
#define val_dim4 (156)

#define val_dim_phi (310)


/******************************************************************************
/INTERPOLATION IDENTIFIERS ****************************************************
******************************************************************************/

#define POLYNOMIAL (1)
#define CUBIC_SPLINE (2)

/******************************************************************************
/TIMESCALE IDENTIFIERS*********************************************************
******************************************************************************/

#define TIMESCALE_MOON (0)
#define TIMESCALE_PLANET (1)
#define TIMESCALE_SECONDS (2)

/******************************************************************************
/OUTPUT IDENTIFIERS ***********************************************************
******************************************************************************/

#define OUTPUT_ASCII (1)
#define OUTPUT_BINARY (2)
#define FILE_OPEN  (1)
#define FILE_CLOSED  (0)

/******************************************************************************
/PLASMA MODEL IDENTIFIERS *****************************************************
******************************************************************************/

#define MUENCH2012_a (1)
#define MUENCH2012_b (2)

/******************************************************************************
/MACROS ***********************************************************************
******************************************************************************/

#define MACRO_square(x) ((x) * (x))

/******************************************************************************
/ERROR HANDLING ***************************************************************
*******************************************************************************/

#ifndef MACRO_ERROR
#define MACRO_ERROR(CONDITION,MESSAGE)\
do\
  {\
    if(CONDITION)\
      {\
	printf("ERROR: (%s:%s:%d:%s)",\
	       __FILE__,__FUNCTION__,__LINE__,#CONDITION);\
	printf  MESSAGE;\
	exit( EXIT_FAILURE );\
      }\
  }\
while(0)
#endif

#ifndef __GNUC__
#define __FUNCTION__       ""             /* only gcc knows this... */
#endif

/******************************************************************************
/MESSAGES *********************************************************************
******************************************************************************/
#ifndef MACRO_MESSAGE
#define MACRO_MESSAGE(MESSAGE)\
  do\
    {\
      printf("MESSAGE (%s:%s:%d): ",\
	     __FILE__,__FUNCTION__,__LINE__);\
      printf  MESSAGE;\
    }\
while(0)
#endif
#ifndef __GNUC__
#define __FUNCTION__       ""             /* only gcc knows this... */
#endif

/******************************************************************************
/TYPES ************************************************************************
******************************************************************************/

typedef struct {

    /*PHYSICAL CONSTANTS*/

    double G;    /*gravity constant*/
    double c;    /*speed of light*/
    double eps0; /*vacuum permittivity*/
    double qe;   /*elementary charge*/
    double me;   /*electron mass*/
    double qm;   /*charge per mass ratio of electron*/
    double mp;   /*proton mass*/
    double au;   /*astronomical unit*/
    double amu;  /*atomic mass unit*/

}
type_constants;

/*****************************************************************************/
#ifndef _dyne_type_elements
#define _dyne_type_elements
/*
 *    ORBIT ELEMENTS
 */
typedef struct {
    double a;                               // semi-major axis
    double e;                               // eccentricity
    double i;                               // inclination
    double aop;                             // argument of pericentre
    double lan;                             // longitude of ascending node
    double Tper;                            // time of pericentre passage
    double om;                              // mean motion
    double lop;                             // longitude of pericentre
    double pdot;                            // precession rate of pericentre node
    double odot;                            // precession rate of ascending node
    double lambda;                          // mean longitude 
}
type_elements;
#endif

/*****************************************************************************/
#ifndef _dyne_type_coordinates
#define _dyne_type_coordinates
typedef struct {

    /*PHASE SPACE COORDINATES*/

    double t;   /*time*/
    double x;   /*coordinates*/
    double y;
    double z;
    double u;   /*velocities*/
    double v;
    double w;

}
type_coordinates;
#endif

/*****************************************************************************/
typedef struct {

    /*CENTRAL PLANET*/

    char name[50];
    double r;               /*radius*/
    type_elements elem;     /*orbit elements*/
//    double m;             /*mass*/
    double GM;
    double om;              /*spin rate*/
    double obliquity;       /*obliquity of planetary spin axis*/
    double cosobl, sinobl;
    double j2;              /*gravity quadrupole: oblateness*/
    double j4;
    double j6;
    double g10;             /*magnetic dipole*/
    double year;            /*planetary year*/

    double f;               /*fraction of corotation of magnetosphere*/
    double n0;              /*background ion density at Enceladus*/
    double fpd;             // constant azimuthal component of the (specific) plasma drag force
}
type_planet;

/*****************************************************************************/
typedef struct
{
    int TIMESCALE;          // time scaling to source moon years or planet years
    double ti;              // starting time
    double tf;	            // end time
    double dt;	            // time step for output
    double torbit;          // orbit period
    double steps_per_orbit; // timesteps per orbit
    double num_orbits;      // number of orbits to integrate
    long nt;	            // number of time-steps

    int box_integration;    // flag to signal if integration shall be stopped when
			                // grain leaves predefined rectangular box
}
type_time_frame;

/*****************************************************************************/
typedef struct {

    double rad2deg;
    double SecToPlanetYears;
    double PlanetYearsToSec;
    double MeterToPlanetRadii;
    double SecToMoonYears;
    double MoonYearsToSec;
    double SecToEarthYears;

}
type_conversion;

/*****************************************************************************/
typedef struct {

    /*SUN*/

//    double m;   /*mass*/
    double GM;  /*mu of the Sun*/
    double J0;  /*modulus of solar radiation energy flux at 1 AU*/
    double RadPress;   /*acceleration coefficient due to radiation pressure*/

    type_coordinates coords;        // initial coordinates of the sun
                                    // in the planetocentric frame in AU, AU/second
}
type_sun;

/*****************************************************************************/
typedef struct {

    /*GRAIN*/

    double r;	    /*radius*/
    double rho;	    /*density*/
    double mass;    /* grain mass */
    double dm;	    /*secondary electron yield*/
    double em;	    /*optimal energy for secondary emission*/
    double kappa;	    /*efficiency of photoemission*/
    double kTs;	    /*thermal energy of secondary electrons*/
    double kTn; 	    /*thermal energy of photo electrons*/
    double phiconst;  /*constant equilibrium potential*/
    double phi_calc;  /*on-the-fly integrated potential*/
    double Qpr;       /*radiation pressure coefficient*/
    int pos_mode;     /*signal to start grain from other position than from moon*/
}
type_grain;


/*****************************************************************************/
typedef struct {

    int filemod;            // type of output -> ascii=1 or binary=2
    char filemodchar[3];    // w or wb
    int n_out;              // fraction of timesteps for which output is written:
                            // every n_out-th integration timestep
    int PSpaceOut;          // output the phase-space coordinates
    int elemout;            // output osculating orbital elements
    int SatCoordOut;        // output coordinates of the satellites
    int SatElemOut;         // output elements of the satellites
    int combined;           // output positions of grain and reference moon 
    int sun;                // output position and velocity of the sun
    int TestsOut;           // output of physical tests of numerical integration
    int JacobiConstant;     // output Jacobi Integral of the restricted three body problem
    int TrafoCoord;         // output coordinates/velocities in corotating frame
    int RelCoord;           // output relative coordinates dust grain-source moon
    int EARTHTIMESCALE;     // scale output time to Earth years (1: yes, 0: no)
    int EM_INTEGRATION;     // calculate and output energy and momentum integrals of
                            // em. interaction
    int EFFECTIVE_SPEED;    // output effective escape speed of dust grain
    int STATISTICS;         // output statistic data of escapers to E Ring
    int PlasmaParOut;       // output of plasma parameters along grain trajectory
    int LorentzOut;         // output of Lorentz force along grain trajectory
    int No_Output;          // make no other output at all - only to save disk space
}
type_output_flags;

/*****************************************************************************/
typedef struct {

    /* SWITCHES */

    int LORENTZ;                // magnetic dipole
    int LORENTZ2;               // switch on plasma-plume interaction
    int CURRENTS;               // integrate currents
    int CONSTPOT;               // use constant pontential
    int EQUIPOT;                // use equilibrium pontential
    int SINKS;                  // switch on sinks for particles
    int RADPRESS;               // switch on radiation pressure of the Sun
    int INTEGRATOR;             // time-integration routine
    int INTERPOLATION;          // set type of interpolation scheme for plasma data
    int PLASMA_MODEL;           // set which plasma model to use MUENCH2012_a or b
    int DIRECT_DRAG;            // switch on the calculation of direct_drag
    int PLASMA_DRAG;            // switch on the calculation of direct_drag
    int PRINT_INFO;             // switch on/off printing of infos onto command line
    type_output_flags output;   // the output flags
}
type_switches;


/*****************************************************************************/
typedef struct stepsize_t {
    double try;
    double did;
    double next;
}
stepsize_t;


typedef enum {
    IMPACTED_SINK,
    SUCCESSFUL,
    STEPSIZE_TOO_SMALL,
    TOO_MANY_STEPS
}
integ_status_t;

typedef struct {
    integ_status_t kind;
    impact_data_t impact;
} integ_data_t;

/*****************************************************************************/

struct type_force_par;

typedef struct type_inte_par
{
    double *atol;   // array of absolute tolerances for adaptive stepsize control
    double rtol;    // relative tolerance for adaptive stepsize control 
    double wy;      // weight factor for err per step
    double wdydt;   // weight factor for err per unit step
    double h1;      // initial stepsize
    double h;       // momentary stepsize
    double hmax;    // maximum allowed stepsize
    double hmin;    // minimum allowed stepsize */
    int first_step; // true for the trajectorie's first integration step, false otherwise */
    integ_status_t (*step)(double y[],
                           double dydt[],
                           double t,
                           stepsize_t *stepsize,
                           double acc_err[],
                           struct type_inte_par *ipar,
                           struct type_force_par *fpar);
}
type_inte_par;



/*****************************************************************************/
typedef struct {

    /*file type*/

    FILE *ID;               // input stream
    FILE *ODpsc;            // output phase space coordinates
    FILE *ODstv;            // output startvalues
    FILE *ODoscel;          // output osculating elements
    FILE *ODstel;           // output starting orbital elements
    FILE *ODphi;            // output integrated grain potential
    FILE *ODs;              // output sinks
    FILE **ODsatInitCoord;  // output initial satellite coordinates,
    FILE **ODsatCoord;      // output satellite coordinates,
    FILE **ODsatElem;       // output satellite elements; satN-times
    FILE *ODcombined;       // output positions of grain and reference moon 
    FILE *ODsun;            // output position and velocity of the sun
    FILE *ODTestsOut;       // output of physical tests of numerical integration
    FILE *ODJacobiConstant; // output Jacobi integral of restricted 3BP
    FILE *ODTrafoCoord;     // coordinates/velocities in corotating system
    FILE *ODRelCoord;       // relative coordinates dust grain - source moon
    FILE *ODeffsp;          // output effective escape speed once it's reached
    FILE *ODplasma_par;     // output of plasma parameters along grain trajectory
    FILE *ODbfield;         // output of B-field along trajectory
    FILE *ODLorentz;        // output of Lorentz force along grain trajectory
    FILE *ODescapers;       // output statistical data for escape rates
    FILE *ODstartparam;     // output start parameters
#ifdef dyne_debug_stepsize
    FILE *ipar;             // integration parameters for diagnosis
    FILE *stepsize_diag;    // stepsize information for diagnosis
#endif
    char name[256];
    char dir[256];
}
type_file;


/*
 *   structure to pass parameters to the forces subroutine
 */ 
typedef struct  type_force_par
{
    type_sun sun;
    type_planet planet;
    type_switches switches;
    type_moons moons;
    type_conversion conversion;
    type_grain grain;
    type_constants constants;
    type_time_frame time_frame;
    double RadPress;   /*acceleration coefficient due to radiation pressure*/
    long int run_count;
    int integration_error;
    int integration_error_alternativ;
    int lagrange_point;
    type_file files;
    psvidx_t *idx;
    gsl_rng *rng;
    double stepsize_eps;
    int (*forces)(double t, double psv[], double dydt[], struct type_force_par *fpar);
}
type_force_par;


/******************************************************************************
****** FUNCTION DECLARATION ***************************************************
******************************************************************************/


/****** INTEGRATORS ******/

integ_data_t odeint(double psv[], double t1, double t2, type_force_par *fpar, type_inte_par *ipar);
void rkdopr853(double y[], double dydt[], long int neq, double t, double h, double yout[],
                double yerr3[], double yerr5[], type_force_par *fpar);
integ_status_t rkdopr853qs(double y[], double dydt[], double t, stepsize_t *stepsize, double acc_err[],
                            type_inte_par *ipar, type_force_par *fpar);


/****** FROM forces.c ******/

int dyne_forces2(double t, double psv[], double dxx[], type_force_par *force_par);


/****** FROM output.c ******/

void WriteCoordinates(double[], type_elements *, type_file *,
                      type_force_par *);
void WriteConstants(type_sun *, type_planet *, type_moons *, type_constants *,
                    type_grain *, type_file *, type_file *, psvidx_t *);
void WriteSinklog(type_file *output_file, type_force_par *fpar, integ_data_t status, double t);


/****** FROM input.c ******/

void read_jobfile(int, char **, type_file *, type_file *, type_switches *,
                  type_grain *, type_sun *, type_planet *, type_moons *,
                  type_time_frame *, type_inte_par *, psvidx_t *, iv_t *);
void read_constants(type_constants *constants, int verbose);
void precalculation(type_constants *, type_sun *, type_planet *,
                    type_conversion *, type_moons *, type_grain *,
                    type_switches *, psvidx_t *);
void prepare_files(type_file *, type_file *, type_switches *, type_moons *, psvidx_t *);


