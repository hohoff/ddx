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


#include "dd.h"
#include "impact.h"

#define MAXSTP 100000
#define TINY 1.0e-30



integ_data_t odeint(double psv[], double t1, double t2, type_force_par *fpar, type_inte_par *ipar)
{
    // Return values
    integ_data_t integ_data;
    integ_status_t status;
    impact_data_t impact;

    // Phase space vectors and index structure
    //
    psvidx_t *idx = fpar->idx;
    double oldpsv[idx->nvars];
    double newpsv[idx->nvars];

    // Integration vectors (possibly different length than phase space vectors)
    //
    double    dydt[idx->neq];      // vector of forces (rh side derivatives)
    double acc_err[idx->neq];      // scaling vector to monitor accuracy

    stepsize_t stepsize = {0.0, 0.0, 0.0};
    stepsize.try = ipar->h;

    double t = t1;                  // start time for current step

    // copy starting values to local phase space vectors
    for (int i = 0; i < idx->nvars; i++)
    {
        oldpsv[i] = psv[i];
        newpsv[i] = psv[i];
    }


    for (int nsteps = 1; nsteps <= MAXSTP; nsteps++)
    {
        // Forces at time t, newpsv is uptodate
        (fpar->forces)(t, newpsv, dydt, fpar);

        // If stepsize can overshoot, decrease
        if (t+stepsize.try > t2)
        {
            stepsize.try = t2 - t;
        }

        // Desired error level vector
        for (int i = 0; i < idx->neq; i++)
            acc_err[i] = ipar->atol[i]
                + ipar->rtol * (ipar->wy*fabs(newpsv[i]) + ipar->wdydt*fabs(dydt[i]*stepsize.try))
                + TINY;

        status = (ipar->step)(newpsv, dydt, t, &stepsize, acc_err, ipar, fpar);

        if (status == STEPSIZE_TOO_SMALL)
        {
            integ_data.kind = status;
            return integ_data;
        }

        t += stepsize.did;
        newpsv[time_idx(idx)] = t;

        if (fpar->switches.SINKS)
        {
            if (!ipar->first_step)
            {
                check4impacts(fpar->planet.GM, &fpar->moons, fpar->idx, oldpsv, newpsv, &impact);

                if (impact.status == IMPACT)
                {
                    integ_data.kind = IMPACTED_SINK;
                    integ_data.impact = impact;
                    return integ_data;
                }
            }

            for (int i = 0; i < idx->nvars; i++)
                oldpsv[i] = newpsv[i];
        }

        if (ipar->first_step)
        {
            ipar->first_step = 0;    // first_step = false
        }

        // Check if integration step is done
        if (t >= t2)
        {
            for (int i = 0; i < idx->nvars; i++)
                psv[i] = newpsv[i];

            ipar->h = stepsize.next;

            integ_data.kind = SUCCESSFUL;
            return integ_data;
        }
        
        stepsize.try = stepsize.next;
    }

    integ_data.kind = TOO_MANY_STEPS;
    return integ_data;
}

double debug_distance2moon(psvidx_t *idx, double psv[], int j)
{
    double dx = psv[grain_pos_idx(idx, 0)] - psv[moon_pos_idx(idx, j, 0)];
    double dy = psv[grain_pos_idx(idx, 1)] - psv[moon_pos_idx(idx, j, 1)];
    double dz = psv[grain_pos_idx(idx, 2)] - psv[moon_pos_idx(idx, j, 2)];
    return sqrt(dx*dx + dy*dy + dz*dz);
}

#undef MAXSTP
#undef TINY

