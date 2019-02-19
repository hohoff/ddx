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


#ifndef dyne_initialvalues_h
#define dyne_initialvalues_h

#ifndef LINE_LENGTH
#define LINE_LENGTH 256
#endif

#include <stdbool.h>
#include "moons.h"
#include "rotations.h"
#include "psvindex.h"

typedef struct
{
    double x;
    double y;
    double z;
    double u;
    double v;
    double w;
}
iv_values_psv_t;


typedef struct
{
    double dummy;
}
iv_data_psv_t;


typedef struct
{
    double phi;
    double theta;
    double v;
}
iv_values_moon_normal_t;

typedef struct
{
    double R;               // moon radius
}
iv_data_moon_normal_t;


typedef union iv_kind_t
{
    iv_values_psv_t psv;
    iv_values_moon_normal_t mn;
    // iv_values_elem_t elem;
}
iv_kind_t;


typedef union iv_data_kind_t
{
    iv_data_psv_t psv;
    iv_data_moon_normal_t mn;
}
iv_data_kind_t;


typedef struct iv_t
{
    int pos_mode;
    FILE *file;
    iv_kind_t values;
    iv_data_kind_t data;
    int moon;
    void (*iv2psv)(double t, struct iv_t *iv, double rm[], double vm[], double r[], double v[]);
    bool (*read_vals)(struct iv_t *iv);
    void (*write_vals)(FILE *file, long int run_count, struct iv_t *iv);
}
iv_t;

bool iv_init(psvidx_t *idx, type_moons *moons, iv_t *iv);

void iv2psv_from_psv(double t, struct iv_t *iv, double rm[], double vm[], double r[], double v[]);
bool iv_read_vals_from_psv(struct iv_t *iv);
void iv_write_vals_from_psv(FILE *file, long int run_count, struct iv_t *iv);

void iv2psv_from_moon_normal(double t, struct iv_t *iv, double rm[], double vm[], double r[], double v[]);
bool iv_read_vals_from_moon_normal(struct iv_t *iv);
void iv_write_vals_from_moon_normal(FILE *file, long int run_count, struct iv_t *iv);
#endif
