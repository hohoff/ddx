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
#include <stdbool.h>
#include <string.h>
#include "initialvalues.h"
#include "moons.h"
#include "rotations.h"
#include "psvindex.h"


static double deg2rad = 0.0174532925199432957692369076849;


bool iv_init(psvidx_t *idx, type_moons *moons, iv_t *iv)
{
    switch(iv->pos_mode)
    {
        case 6:
            iv->moon = moons->ref_moon;
            iv->data.psv.dummy = 0.0;
            iv->iv2psv = iv2psv_from_psv;
            iv->read_vals = iv_read_vals_from_psv;
            iv->write_vals = iv_write_vals_from_psv;
            break;
        case 9:
            iv->moon = moons->sourceSat;
            iv->data.mn.R = moons->satellite[moon_array_idx(idx, iv->moon)].r;
            iv->iv2psv = iv2psv_from_moon_normal;
            iv->read_vals = iv_read_vals_from_moon_normal;
            iv->write_vals = iv_write_vals_from_moon_normal;
            break;
        default:
            return false;           // unknown pos_mode
    }
    return true;
}


static bool is_comment(char *line)
{
    if ( strncmp(line, "#", 1) == 0 )
    {
        return true;
    }
    return false;
}

void iv2psv_from_moon_normal(double t, struct iv_t *iv, double rm[], double vm[], double r[], double v[])
{
    double ev[3]  = {0., 0., 0.};

    double phi = deg2rad * iv->values.mn.phi;
    double theta = deg2rad * iv->values.mn.theta;

    rotate_on_sphere(phi, theta, rm, ev);

    double Rmoon = iv->data.mn.R;
    r[0] = rm[0] + Rmoon * ev[0];
    r[1] = rm[1] + Rmoon * ev[1];
    r[2] = rm[2] + Rmoon * ev[2];

    double vstart = iv->values.mn.v;
    v[0] = vm[0] + vstart * ev[0];
    v[1] = vm[1] + vstart * ev[1];
    v[2] = vm[2] + vstart * ev[2];
}


bool iv_read_vals_from_moon_normal(struct iv_t *iv)
{
    char line[LINE_LENGTH + 1];

    do
    {   
        if (fgets(line, LINE_LENGTH, iv->file) == NULL)
        {
            return false;
        }
    }   
    while (is_comment(line));

    int n = sscanf(line, "%lf %lf %lf", &iv->values.mn.phi,
                                        &iv->values.mn.theta,
                                        &iv->values.mn.v);
    if (n != 3)
    {
        printf("Wrong format of input line: %s\n", line);
        return false;
    }
    return true;

}


void iv_write_vals_from_moon_normal(FILE *file, long int run_count, struct iv_t *iv)
{
    fprintf(file, "%ld %.15e %.15e %.15e\n", run_count,
                                             iv->values.mn.phi,
                                             iv->values.mn.theta,
                                             iv->values.mn.v);
}

void iv2psv_from_psv(double t, struct iv_t *iv, double rm[], double vm[], double r[], double v[])
{
    r[0] = iv->values.psv.x;
    r[1] = iv->values.psv.y;
    r[2] = iv->values.psv.z;

    v[0] = iv->values.psv.u;
    v[1] = iv->values.psv.v;
    v[2] = iv->values.psv.w;
}


bool iv_read_vals_from_psv(struct iv_t *iv)
{
    char line[LINE_LENGTH + 1];

    do
    {   
        if (fgets(line, LINE_LENGTH, iv->file) == NULL)
        {
            return false;
        }
    }   
    while (is_comment(line));

    int n = sscanf(line, "%lf %lf %lf %lf %lf %lf", &iv->values.psv.x,
                                                    &iv->values.psv.y,
                                                    &iv->values.psv.z,
                                                    &iv->values.psv.u,
                                                    &iv->values.psv.v,
                                                    &iv->values.psv.w);
    if (n != 6)
    {
        printf("Wrong format of input line: %s\n", line);
        return false;
    }
    return true;
}


void iv_write_vals_from_psv(FILE *file, long int run_count, struct iv_t *iv)
{
    fprintf(file, "%ld %.15e %.15e %.15e %.15e %.15e %.15e\n", run_count,
                                                               iv->values.psv.x,
                                                               iv->values.psv.y,
                                                               iv->values.psv.z,
                                                               iv->values.psv.u,
                                                               iv->values.psv.v,
                                                               iv->values.psv.w);
}



