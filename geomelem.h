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


#ifndef _dyne_val_pi
#define _dyne_val_pi
#define val_pi (3.14159265358979323846264338328)
#endif

#ifndef _dyne_val_2pi
#define _dyne_val_2pi
#define val_2pi (6.28318530717958647692528676656)
#endif

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


#ifndef _dyne_type_coordinates
#define _dyne_type_coordinates
/*
 *    PHASE SPACE COORDINATES
 */
typedef struct {
    double t;                               // time
    double x;                               // coordinates
    double y;                               //
    double z;                               //
    double u;                               // velocities
    double v;                               //
    double w;                               //
}
type_coordinates;
#endif


#ifndef _dyne_geomelem_h
#define _dyne_geomelem_h

// Input Parameter: double GM_planet, double radius_planet, double J2, double J4, double J6
void geomElemInit(double, double, double, double, double);

// Input Parameter:  double GM_moon, type_elements *orbital_elements_moon
// Output Parameter: double *n, double *kappa, double* nu
void geomElemFreq(double, type_elements*, double*, double*, double*);

void circOrbit_elem2psv(double t, type_elements *elem, double pos[], double vel[]);

// Input Parameter:  double t, double GM_moon, type_elements *orbital_elements_moon
// Output Parameter: type_coordinates *coord_moon
void geomElem2Coord(double, double, type_elements*, type_coordinates*);

// Input Parameter:  double GMm, type_elements *orbital_elements, double nu, double kappa
// Output Parameter: double *eta_sq, double *chi_sq, double *alpha1, double *alpha2, double *alpha_sq
void geomElemHelpFreq( double, type_elements*, double, double, double*, double*, double*, double*, double*);

// Input Parameter:  type_coordinates *coordinates
// Output Parameter: double *L
void Longitude(type_coordinates*, double*);

// Input Parameter:  type_elements *orbital_elements, double r, double rc, double rdot, double rcdot, double kappa, double lambda
// Output Parameter: double *lop
void LongOfPeri(type_elements*, double, double, double, double, double, double, double*);

// Input Parameter:  double z, double zc, double zdot, double zcdot, double nu, double lambda
// Output Parameter: double *lan
void LongAscNode(double z, double zc, double zdot, double zcdot, double nu, double lambda, double *lan );

// Input Parameter:  double t, double GMm, type_coordinates *coordinates
// Output Parameter: type_elements *orbital_elements
void Coord2geomElem( double t, double GMm, type_coordinates *coord, type_elements *elem );

// Input Parameter:  double GMm, type_elements *orbital_elements
// Output Parameter: double *dn, double *dkappa, double *dnu
void DgeomElemFreq(double, type_elements*, double*, double*, double*);

// Input Parameter:  double GMm
// Output Parameter: type_elements *orbital_elements
void SemiMajorAxis(double, type_elements*);

#endif
