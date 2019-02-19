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
 * along with DDX.  If not, see <https://www.gnu.org/licenses/>.\
 * 
 */


#include "dd.h"
#include "math.h"
// #include "gsl/gsl_integration.h"
// #include "gsl/gsl_sf_erf.h"

/***************************************************************************************/
void gravity_spherical_planet(double mu, double pos[], double acc[])
{
    double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    double r3 = r*r*r;

    acc[0] += - mu * pos[0] / r3;
    acc[1] += - mu * pos[1] / r3;
    acc[2] += - mu * pos[2] / r3;
}


/***************************************************************************************/
void orbit_sun(double mu, double pos[], double acc[])
{
    // Distances are measured in astronomical units and
    // velocities in AU/second
    double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    double r3 = r*r*r;

    acc[0] += - mu * pos[0] / r3;
    acc[1] += - mu * pos[1] / r3;
    acc[2] += - mu * pos[2] / r3;
}


/***************************************************************************************/
void gravity_oblate_planet(type_planet *planet, double pos[], double acc[])
{
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double r = sqrt(x*x + y*y + z*z);

    double rho   = planet->r/r;
    double zeta  = z/r;
    double gamma = planet->GM/(r*r*r);

    /* J2-oblateness prefactor */
    double fj2_fac = 3.* gamma* planet->j2* rho* rho* 0.5;             /*   1/2 = 0.5 */ 
    /* J4-oblateness prefactor */
    double fj4_fac = 5.* gamma* planet->j4* rho* rho* rho* rho* 0.125; /*   1/8 = 0.125 */
    /* J6-oblateness prefactor */
    double fj6_fac = 7.* gamma* planet->j6* pow(rho,6)* 0.0625;        /*   1/16= 0.0625 */

    /* x-component of Gravity due to oblateness on Dust grain */
    double fj2 = fj2_fac* (5.* zeta* zeta - 1);
    double fj4 = fj4_fac* (63*pow(zeta,4) - 42*zeta*zeta + 3);
    double fj6 = fj6_fac* (429*pow(zeta,6) - 495*pow(zeta,4) + 135*zeta*zeta - 5);
    acc[0] += (fj2 + fj4 + fj6)*x;

    /* y-component of Gravity due to oblateness on Dust grain */
    acc[1] += (fj2 + fj4 + fj6)*y;

    /* z-component of Gravity due to oblateness on Dust grain */
    fj2 = fj2_fac* (5.* zeta* zeta - 3);
    fj4 = fj4_fac* (63*pow(zeta,4) - 70*zeta*zeta + 15);
    fj6 = fj6_fac* (429*pow(zeta,6) - 693*pow(zeta,4) + 315*zeta*zeta - 35);
    acc[2] += (fj2 + fj4 + fj6)*z;
}



/***************************************************************************************/
void gravity_spherical_moon(double GMm, double grain_pos[], double moon_pos[], double acc[])
{
    double dx = grain_pos[0] - moon_pos[0];
    double dy = grain_pos[1] - moon_pos[1];
    double dz = grain_pos[2] - moon_pos[2];
    double d2 = dx*dx + dy*dy + dz*dz;
    double d3 = d2 * sqrt(d2);

    double xm = moon_pos[0];
    double ym = moon_pos[1];
    double zm = moon_pos[2];
    double rm2 = xm*xm + ym*ym + zm*zm;
    double rm3 = rm2 * sqrt(rm2);

    // keep in mind, this is for a planetocentric coordinate
    // frame, a term -ddot{r}_p appears in this case
    acc[0] += -GMm * dx / d3 - GMm * xm / rm3;
    acc[1] += -GMm * dy / d3 - GMm * ym / rm3;
    acc[2] += -GMm * dz / d3 - GMm * zm / rm3;
    //        -----|--------   ------|-------
    //             v                 v
    //           ddot{r}          ddot{rp}
}



/*****************************************************************************/
void lorentz_force(type_planet *planet,
                   type_grain *grain,
                   double eps0,
                   double phi,
                   double grain_pos[],
                   double grain_vel[],
                   double acc[])
{
    /* grain position and distance to planet */
    double x = grain_pos[0];
    double y = grain_pos[1];
    double z = grain_pos[2];
    double r = sqrt(x*x + y*y + z*z);
    double r3 = r*r*r;

    /* radius of planet */
    double Rp3 = (planet->r * planet->r * planet->r);

    /* prefactors */
    double LfacB = planet->g10 * Rp3 / r3;      // g10 of central planet
    double ct = z/r;                            // ct = cos(polar angle)
    double DynfacB = 3./r*ct;

    /* grain parameters */
    double rho = grain->rho;                    // density of dust grain
    double rg2 = grain->r * grain->r;           // squared radius of dust grain

    /* relative velocity to corotation */
    double u_rel = grain_vel[0] + planet->f * y * planet->om;  // f: fraction of corotation
    double v_rel = grain_vel[1] - planet->f * x * planet->om;
    double w_rel = grain_vel[2];

    /* auxiliary quantities for calculating Lorentz force */
    double cross_x = v_rel * z - w_rel * y;
    double cross_y = w_rel * x - u_rel * z;
    double cross_z = u_rel * y - v_rel * x;

    double Qm = 3.*eps0*phi/(rho*rg2);     // current grain charge-to-mass ratio*/

    /* Lorentz forces due to planetary magnetic field including
       corotational part */
    acc[0] += Qm * LfacB * (DynfacB*cross_x - v_rel);
    acc[1] += Qm * LfacB * (DynfacB*cross_y + u_rel);
    acc[2] += Qm * LfacB *  DynfacB*cross_z;
}



void radiation_pressure(type_constants *constants, double Rp, double konst, double t, double r[], double v[], double rs[], double vs[], double acc[])
{
    /*
     * For details see
     *
     *     X. Liu, M. Sachse, F. Spahn, and J. Schmidt (2016),
     *     Dynamics and distribution of Jovian dust ejected from the Galilean satellites,
     *     J. Geophys. Res. Planets, 121, 1141–1173,
     *     doi:10.1002/2016JE004999.
     *
     */

    double dot = r[0]*rs[0] + r[1]*rs[1] + r[2]*rs[2];
    if (dot < 0.0)      // alpha > 90 deg, planetary shadow possible
    {
        double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        double rs2 = rs[0]*rs[0] + rs[1]*rs[1] + rs[2]*rs[2];
        double Rp2 = Rp*Rp;
        if (r2 - dot*dot/rs2 < Rp2)     // True, if dust grain is inside the planetary shadow
        {
            return;
        }
    }

    double c = constants->c;
    
    double xi[3] = {0.0, 0.0, 0.0};
    xi[0] = (r[0] - rs[0]) / constants->au;
    xi[1] = (r[1] - rs[1]) / constants->au;
    xi[2] = (r[2] - rs[2]) / constants->au;
    
    double xi_sq  = xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2];
    double xi_abs = sqrt(xi_sq);
    
    double e_xi[3] = {0.0, 0.0, 0.0};
    e_xi[0] = xi[0] / xi_abs;
    e_xi[1] = xi[1] / xi_abs;
    e_xi[2] = xi[2] / xi_abs;
    
    double xi_dot[3] = {0.0, 0.0, 0.0};
    xi_dot[0] = v[0] - vs[0];
    xi_dot[1] = v[1] - vs[1];
    xi_dot[2] = v[2] - vs[2];
    
    double dotp = xi_dot[0]*e_xi[0] + xi_dot[1]*e_xi[1] + xi_dot[2]*e_xi[2];
    double term1 = konst/xi_sq;
    double term2 = 1.0 - dotp/c;
    
    acc[0] += term1 * ( term2 * e_xi[0] + xi_dot[0]/c );
    acc[1] += term1 * ( term2 * e_xi[1] + xi_dot[1]/c );
    acc[2] += term1 * ( term2 * e_xi[2] + xi_dot[2]/c );

    return;
}



void zero(double array[], int len);
void zero(double array[], int len)
{
    for(int j = 0; j < len; j++)
    {
        array[j] = 0.0;
    }
}



/*****************************************************************************/
int dyne_forces2(double t, double psv[], double dxx[], type_force_par *sim)
{
    /*
     * dyne_forces provides all right-hand sides for the numerical integration
     * 
     * Input:  *t:   time at which forces are evaluated;
     *               i.e., time at beginning of step
     *       psv[]:  vector of integrated phase space values
     * 
     * Output: dxx[]: vector containing the right-hand sides
     *                of the integrated phase space values
     */

    /* index structue for access to the phasespacevector (psv) */
    psvidx_t *idx = sim->idx;

    
    zero(dxx, idx->neq);            // initialise vector of right hand sides to zero,
                                    // this is *important* because we only add to this
                                    // vector and it is never overwritten
   
/*-------------------------------------------------------------------------------------*/
/*          initialise phase space data                                                */
/*-------------------------------------------------------------------------------------*/

    // we can't be sure that psv is properly initialised for non-integrated quantities
    psv[time_idx(idx)] = t;
    psv[grain_pot_idx(idx)] = sim->grain.phiconst;


/*-------------------------------------------------------------------------------------*/
/*          moon accelerations                                                         */
/*-------------------------------------------------------------------------------------*/
    
    for (int j = moon_start_idx(); exist_moon(idx, j); j++)
    {
        dxx[moon_pos_idx(idx, j, 0)] = psv[moon_vel_idx(idx, j, 0)];     // change of moon position is 
        dxx[moon_pos_idx(idx, j, 1)] = psv[moon_vel_idx(idx, j, 1)];     // given by moon velocity times
        dxx[moon_pos_idx(idx, j, 2)] = psv[moon_vel_idx(idx, j, 2)];     // timestep
    
        double mu = sim->planet.GM + sim->moons.satellite[j-1].GM;
        gravity_spherical_planet(mu, &psv[moon_pos_idx(idx, j, 0)],
                                     &dxx[moon_vel_idx(idx, j, 0)]);
    
        // moon acceleration due to gravity of oblate planet
        gravity_oblate_planet(&sim->planet, &psv[moon_pos_idx(idx, j, 0)],
                                            &dxx[moon_vel_idx(idx, j, 0)]);
    
        // moon acceleration due to gravity of the other moons
        for (int k = moon_start_idx(); exist_moon(idx, k); k++)
        { 
            if ( k != j )
            {  
              gravity_spherical_moon(sim->moons.satellite[moon_array_idx(idx, k)].GM,
                                   &psv[moon_pos_idx(idx, j, 0)],
                                   &psv[moon_pos_idx(idx, k, 0)],
                                   &dxx[moon_vel_idx(idx, j, 0)]);
            }
        }
    }
     

/*-------------------------------------------------------------------------------------*/
/*          Orbit of the sun relative to Saturn                                        */
/*-------------------------------------------------------------------------------------*/
    
    dxx[sun_pos_idx(idx, 0)] = psv[sun_vel_idx(idx, 0)];    // change of sun position is 
    dxx[sun_pos_idx(idx, 1)] = psv[sun_vel_idx(idx, 1)];    // given by sun velocity times
    dxx[sun_pos_idx(idx, 2)] = psv[sun_vel_idx(idx, 2)];    // timestep

    orbit_sun(sim->sun.GM + sim->planet.GM, &psv[sun_pos_idx(idx, 0)],  &dxx[sun_vel_idx(idx, 0)]);


/*-------------------------------------------------------------------------------------*/
/*          dust grain accelerations                                                   */
/*-------------------------------------------------------------------------------------*/
    
    dxx[grain_pos_idx(idx, 0)] = psv[grain_vel_idx(idx, 0)];    // change of dust grain position
    dxx[grain_pos_idx(idx, 1)] = psv[grain_vel_idx(idx, 1)];    // is given by grain velocity times
    dxx[grain_pos_idx(idx, 2)] = psv[grain_vel_idx(idx, 2)];    // timestep


/*-------------------------------------------------------------------------------------*/
/*         Gravity of spherical planet                                                 */
/*-------------------------------------------------------------------------------------*/
    
    // dust grain acceleration due to gravity of spherical planet
    gravity_spherical_planet(sim->planet.GM, &psv[grain_pos_idx(idx, 0)],  &dxx[grain_vel_idx(idx, 0)]);

    // dust grain acceleration due to additional gravity term of oblate planet
    gravity_oblate_planet(&sim->planet, &psv[grain_pos_idx(idx, 0)],
                                        &dxx[grain_vel_idx(idx, 0)]);
     
    // dust grain acceleration due to gravity of the moons
    for (int j = moon_start_idx(); exist_moon(idx, j); j++)
    {
        gravity_spherical_moon(sim->moons.satellite[moon_array_idx(idx, j)].GM,
                               &psv[grain_pos_idx(idx, 0)],
                               &psv[moon_pos_idx(idx, j, 0)],
                               &dxx[grain_vel_idx(idx, 0)]);
    }


/*-------------------------------------------------------------------------------------*/
/*          Lorentz force onto dust grain                                              */
/*-------------------------------------------------------------------------------------*/

    if (sim->switches.LORENTZ)
    {
        //dust grain acceleration due to lorentz force
        lorentz_force(&sim->planet,
                      &sim->grain,
                      sim->constants.eps0,
                      psv[grain_pot_idx(idx)],
                      &psv[grain_pos_idx(idx, 0)],
                      &psv[grain_vel_idx(idx, 0)],
                      &dxx[grain_vel_idx(idx, 0)]);
    }


/*-------------------------------------------------------------------------------------*/
/*          acceleration due to the radiation pressure                                 */
/*-------------------------------------------------------------------------------------*/

    if (sim->switches.RADPRESS)
    {
        radiation_pressure(&sim->constants,             
                           sim->planet.r,                           // Radius of Saturn
                           sim->RadPress,                           // Constant term for the calculation of the radiation pressure
                           t,                                       // Time
                           &psv[grain_pos_idx(idx, 0)],             // Position of the dust grain
                           &psv[grain_vel_idx(idx, 0)],             // Velocitu of the dust grain
                           &psv[sun_pos_idx(idx, 0)],               // Position of the sun
                           &psv[sun_pos_idx(idx, 0)],               // Velocity of the sun
                           &dxx[grain_vel_idx(idx, 0)]);            // Acceleration of the dust grain due to radiation pressure
    }


    return 0;   
}

