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


#ifndef dyne_moons_h
#define dyne_moons_h

#include "geomelem.h"

typedef struct
{
    char name[50];
    double r;                           // radius
//    double m;                           // mass
    double GM;                          // mu of the moon
    type_elements Sat_elem;             // orbital elements of satellite
    type_coordinates Sat_coord;         // initial coordinates of satellite
    double h;                           // hill radius
    double vesc;                        // escape velocity
    int sink;                           // active as a particle sink
    double r_mean;                      // orbit averaged radius
    double year;                        // orbital period
    double MeterPerSecToOrbitSpeed;     // conversion from m/s to orbital speed
}
type_satellite;


typedef struct type_moons
{
    int satN;                   // number of satellites given in jobfile
    int sourceSat;              // which satellite is source (0 = NONE)
    int ref_moon;               // which satellite is the reference for
                                // the corotating coordinate system
    type_satellite *satellite;
    void (*psv)(double t, type_elements *elem, double pos[], double vel[]);
}
type_moons;

#endif

