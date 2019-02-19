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


// Module for calculations to decide if grain has impacted a sink
// 
//   see
// 
//     X. Liu, M. Sachse, F. Spahn, and J. Schmidt (2016),
//     Dynamics and distribution of Jovian dust ejected from the Galilean satellites,
//     J. Geophys. Res. Planets, 121, 1141–1173,
//     doi:10.1002/2016JE004999.
// 
//   and
// 
//     J. E. Chambers (1999),
//     A hybrid symplectic integrator that permits close encounters between massive bodies,
//     Mon. Not. R. Astron. Soc. 304, 793-799,
//     doi:10.1046/j.1365-8711.1999.02379.x
// 
//   for details.


#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "impact.h"

static double val_eps = 1.e-50;

static inline double norm_sq(double x, double y, double z)
{
    return x*x + y*y + z*z;
}

static double min_dist2saturn(double GM, psvidx_t *idx, double psv[])
{
    // distance to Saturn
    double r_sq = norm_sq(psv[grain_pos_idx(idx, 0)],
                          psv[grain_pos_idx(idx, 1)],
                          psv[grain_pos_idx(idx, 2)]);

    // speed (squared) relative to Saturn
    double v_sq = norm_sq(psv[grain_vel_idx(idx, 0)],
                          psv[grain_vel_idx(idx, 1)],
                          psv[grain_vel_idx(idx, 2)]);

    // \vec{r} \dot \vec{v}
    double rv =   psv[grain_pos_idx(idx, 0)]*psv[grain_vel_idx(idx, 0)]
                + psv[grain_pos_idx(idx, 1)]*psv[grain_vel_idx(idx, 1)]
                + psv[grain_pos_idx(idx, 2)]*psv[grain_vel_idx(idx, 2)];

    // return semi latus rectum p = h^2/mu = a*(1-e^2)
    return (r_sq*v_sq - rv*rv) / GM;
}


// This is a utility structure which should not be visible
// outside of this file.
typedef struct {
    double val;
    double deriv;
} distance_t;


static inline distance_t squared_distance2moon(psvidx_t *idx, double psv[], int j)
{
    double dx = psv[grain_pos_idx(idx, 0)] - psv[moon_pos_idx(idx, j, 0)];
    double dy = psv[grain_pos_idx(idx, 1)] - psv[moon_pos_idx(idx, j, 1)];
    double dz = psv[grain_pos_idx(idx, 2)] - psv[moon_pos_idx(idx, j, 2)];
    double du = psv[grain_vel_idx(idx, 0)] - psv[moon_vel_idx(idx, j, 0)];
    double dv = psv[grain_vel_idx(idx, 1)] - psv[moon_vel_idx(idx, j, 1)];
    double dw = psv[grain_vel_idx(idx, 2)] - psv[moon_vel_idx(idx, j, 2)];

    distance_t dist;
    dist.val   = dx*dx + dy*dy + dz*dz;
    dist.deriv = 2.0 * (dx*du + dy*dv + dz*dw);
    return dist;
}


static int inrange(double t0, double t1, double t)
{
    return (t0 <= t && t <= t1);
}


static double impact_time(double t0, double t1, double d0, double d1, double rmoon_sq)
{
    double m = (d1 - d0) / (t1 - t0);
    double n = (t1*d0 - t0*d1) / (t1 - t0);

    double t_imp = (rmoon_sq - n) / m;

    if ( !inrange(t0, t1, t_imp) )
    {
        t_imp = -1.;
    }
    
    return t_imp;
}


#ifndef val_2pi
#define val_2pi (6.28318530717958647692528676656)
#endif


// We interpolate the squared distance with a third order polynomial and
// assume that there is only one minimal value in the interval [told, tnew]
static double min_sq_dist(double h, distance_t d0, distance_t d1)
{
    double d2min = 0.0;

    // Condition for a minimum in (told, tnew) is d0.deriv*h < 0 && d1.deriv*h > 0.
    if (d0.deriv * h > 0. || d1.deriv * h < 0.)
    {
        // minimum at the boundaries
        if (d0.val <= d1.val)
        {
            // d0 is minimal distance to j'th moon in [told, tnew]
            d2min = d0.val;
        }
        else
        {
            // d1 is minimal distance to j'th moon
            d2min = d1.val;
        }
    }
    else
    {
        // minimum in (told, tnew)
        double temp = 6. * (d0.val - d1.val);
        double A = temp + 3. * h * (d0.deriv + d1.deriv);
        double B = temp + 2. * h * (d0.deriv + 2.*d1.deriv);
        double C = h * d1.deriv;

        temp = -0.5 * (B + sqrt(fmax(B*B - 4.*A*C, 0.)));
        double tau = 0.0;
        if (fabs(temp) > val_eps)
        {
            tau = C/temp;
        }

        tau = fmin(tau,  0.);
        tau = fmax(tau, -1.);

        temp = 1. + tau;
        d2min = tau*tau * ((3.+2.*tau)*d0.val + temp*h*d0.deriv) + temp*temp * ((1.-2.*tau)*d1.val + tau*h*d1.deriv);
        d2min = fmax(d2min, 0.);
    }
    return d2min;
}


// see Spitale and Porco (2009), AJ 138:1520-1528.
// see also Moutamid et al. (2016), Icarus 279:125-140
#ifndef RING_RADIUS
#define RING_RADIUS (136768.9e3)
#endif
void check4impacts(double GMplanet, type_moons *moons, psvidx_t *idx,
                   double oldpsv[], double newpsv[], impact_data_t *impact)
{
    /**** check for impacts to moons ****/
    distance_t d0 = {0.0, 0.0};        // squared distance and deriv. of distance to sink at told
    distance_t d1 = {0.0, 0.0};        // squared distance and deriv. of distance to sind at tnew

    double told = oldpsv[time_idx(idx)];
    double tnew = newpsv[time_idx(idx)];
    double h = tnew - told;

    for(int j = moon_start_idx(); exist_moon(idx, j); j++)
    {
        if (moons->satellite[moon_array_idx(idx, j)].sink)
        {
            double rmoon = moons->satellite[moon_array_idx(idx, j)].r;
            double rmoon_sq = rmoon*rmoon;

            // old and new squared distance and deriv. of distance to moon
            d0 = squared_distance2moon(idx, oldpsv, j);
            d1 = squared_distance2moon(idx, newpsv, j);

            // calculate minimum squared distance to moon in [told, tnew]
            double d2min = min_sq_dist(h, d0, d1);

            if (d2min <= rmoon_sq)
            {
                impact->status = IMPACT;
                strncpy(impact->sink_name, moons->satellite[moon_array_idx(idx, j)].name, MOON_NAME_LENGTH-1);

                impact->time = impact_time(told, tnew, d0.val, d1.val, rmoon_sq);
                if (impact->time < 0.)
                {
                    // return of impact_time not in [told, tnew]
                    impact->time = 0.5 * (tnew + told);
                }
                return;
            }
        }
    }

    /**** check for impacts to A ring ****/
    if (RING_RADIUS >= min_dist2saturn(GMplanet, idx, newpsv))
    {
        impact->status = IMPACT;
        strncpy(impact->sink_name, "A ring", MOON_NAME_LENGTH-1);
        impact->time = 0.5 * (tnew - told);
    }
    return;
}

