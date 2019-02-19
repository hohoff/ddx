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
#include <math.h>
#include "rotations.h"

void rotate_vec(double axis[], double angle, double vec[], double rvec[])
{
    double c = cos(angle);
    double s = sin(angle);
    double omc = 1 - c;
    
    double nx = axis[0];
    double ny = axis[1];
    double nz = axis[2];

    double R00 = c + nx*nx * omc;
    double R01 =     nx*ny * omc - nz * s;
    double R02 =     nx*nz * omc + ny * s;

    double R10 =     ny*nx * omc + nz * s;
    double R11 = c + ny*ny * omc;
    double R12 =     ny*nz * omc - nx * s;

    double R20 =     nz*nx * omc - ny * s;
    double R21 =     nz*ny * omc + nx * s;
    double R22 = c + nz*nz * omc;

    rvec[0] = R00 * vec[0] + R01 * vec[1] + R02 * vec[2];
    rvec[1] = R10 * vec[0] + R11 * vec[1] + R12 * vec[2];
    rvec[2] = R20 * vec[0] + R21 * vec[1] + R22 * vec[2];
}


void rotate_on_sphere(double phi, double theta, double r[], double rvec[])
{
    double R = sqrt(r[0]*r[0] + r[1]*r[1]);

    double ex[3] = {0., 0., 0.};
    ex[0] = r[0]/R;
    ex[1] = r[1]/R;

    double ey[3] = {0., 0., 0.};
    ey[0] = -r[1]/R;
    ey[1] = r[0]/R;

    double ez[3] = {0., 0., 1.};

    double aux[3] = {0., 0., 0.};
    rotate_vec(ey, theta, ez, aux);
    rotate_vec(ez, phi, aux, rvec);
}


