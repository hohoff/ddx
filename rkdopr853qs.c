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
#define SAFETY 0.9
#define ALPHA 0.125
#define MINSCALE 0.333
#define MAXSCALE 6.0

// Function prototype of the step function
void rkdopr853(double y[], double dydt[], long int neq, double t, double h,
              double yout[], double yerr3[], double yerr5[], type_force_par *fpar);



integ_status_t rkdopr853qs(double y[],
                           double dydt[],
                           double t,
                           stepsize_t *stepsize,
                           double acc_err[],
                           type_inte_par *ipar,
                           type_force_par *fpar)
{
    double err = 0.0, h = 0.0, scale = 0.0;
    int neq = fpar->idx->neq;
    int nvars = fpar->idx->nvars;

    double yerr3[neq];    // error estimate, will be calculated by rkdopr853
    double yerr5[neq];    // error estimate, will be calculated by rkdopr853
    double ytemp[nvars];  // new phase space data after successful Runge-Kutta step

    h = stepsize->try;

    /* Runge-Kutta step */
    for ( ; ; )
    {
        /* Take a step */
        rkdopr853(y, dydt, neq, t, h, ytemp, yerr3, yerr5, fpar);

        /* Evaluate accuracy */
        double err3 = 0.0, err5 = 0.0;
        for(int i=0; i < neq; i++)
        {
            double acc = acc_err[i]*acc_err[i];
            err3 += yerr3[i]*yerr3[i]/acc;    
            err5 += yerr5[i]*yerr5[i]/acc;    
        }
        double deno = err5 + 0.01*err3;
        if (deno < 0.0)
        {
            deno = 1.0;
        }
        err = h * err5 * sqrt(1.0/(deno * (double) neq));

        // If step succeeded, compute size of next step and break
        if (err <= 1.0)
        {
            if (err == 0.0)
            {
                scale = MAXSCALE;
            }
            else
            {
                scale = SAFETY * pow(err, -ALPHA);
                if (scale < MINSCALE) scale=MINSCALE;
                if (scale > MAXSCALE) scale=MAXSCALE;
            }
            stepsize->next = scale * h;

            // restrict stepsize to be smaller than hmax
            if (stepsize->next > ipar->hmax) stepsize->next = ipar->hmax;

            stepsize->did = h;

            for(int i = 0; i < neq; i++)
                y[i] = ytemp[i];

            return SUCCESSFUL;
        }

        // Truncation error too large, reduce stepsize
        scale = fmax(SAFETY * pow(err, -ALPHA), MINSCALE);
        h *= scale;

        // Error if stepsize is too small
        if (fabs(h) <= ipar->hmin)
        {
            break;
        }

        // Error if stepsize does not change t anymore
        if (t+h == t)
        {
           break;
        }
    }
    return STEPSIZE_TOO_SMALL;
}

