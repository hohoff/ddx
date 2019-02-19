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
#include "geomelem.h"
#include "psvindex.h"
#include "moons.h"
#include "time.h"
#include "initialvalues.h"

FILE* open_file(char *filename, char *mode);

void initialise_moons(psvidx_t *idx, type_satellite *satellites, double GMplanet, double t, double psv[])
{
    for (int j = moon_start_idx(); exist_moon(idx, j); j++)
    {
        psv[moon_pos_idx(idx, j, 0)] = 
            satellites[moon_array_idx(idx, j)].Sat_coord.x;
        psv[moon_pos_idx(idx, j, 1)] =
            satellites[moon_array_idx(idx, j)].Sat_coord.y;
        psv[moon_pos_idx(idx, j, 2)] =
            satellites[moon_array_idx(idx, j)].Sat_coord.z;
    
        psv[moon_vel_idx(idx, j, 0)] =
            satellites[moon_array_idx(idx, j)].Sat_coord.u;
        psv[moon_vel_idx(idx, j, 1)] =
            satellites[moon_array_idx(idx, j)].Sat_coord.v;
        psv[moon_vel_idx(idx, j, 2)] =
            satellites[moon_array_idx(idx, j)].Sat_coord.w;
    }
}

void initialise_sun(psvidx_t *idx, type_sun *sun, double psv[])
{
    psv[sun_pos_idx(idx, 0)] = sun->coords.x;
    psv[sun_pos_idx(idx, 1)] = sun->coords.y;
    psv[sun_pos_idx(idx, 2)] = sun->coords.z;
    
    psv[sun_vel_idx(idx, 0)] = sun->coords.u;
    psv[sun_vel_idx(idx, 1)] = sun->coords.v;
    psv[sun_vel_idx(idx, 2)] = sun->coords.w;

    printf("Sun coords.x = %.15e \n", sun->coords.x);
    printf("Sun coords.y = %.15e \n", sun->coords.y);
    printf("Sun coords.z = %.15e \n", sun->coords.z);
}


/*--------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------*/
/*--------- OneParticle ----------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------*/
void OneParticle(type_time_frame *time_frame, type_inte_par *inte_par,
                 type_conversion *conversion, type_force_par *force_par,
                 type_elements *elem, type_file *output_file,
                 type_file *input_file, psvidx_t *idx, iv_t *iv)
{

    double t = time_frame->ti;                 /* time initialised to the starting time */
    double psv[idx->nvars];                    /* vector of phase space values */
    double dt = 0.0, tend = 0.0;

    integ_data_t status;

    int n_out = force_par->switches.output.n_out;     // get desired fraction of output
                                                      // steps to integration steps

    long int n;
   

    // initial conditions of the sun
    initialise_sun(idx, &force_par->sun, psv);

    // initial conditions of the mooons
    initialise_moons(idx, force_par->moons.satellite,
                          force_par->planet.GM,
                          t, psv);

    // initial conditions of the grain
    iv->iv2psv(t, iv,                                   // initial value structure
               &psv[moon_pos_idx(idx, iv->moon, 0)],    // rm[]:    moon position
               &psv[moon_vel_idx(idx, iv->moon, 0)],    // vm[]:    moon velocity
               &psv[grain_pos_idx(idx, 0)],             // r[]:     grain position
               &psv[grain_vel_idx(idx, 0)]);            // v[]:     grain velocity
    iv->write_vals(output_file->ODstartparam, force_par->run_count, iv);


    /* ---------- Grain potential ---------- */
    if (force_par->switches.CONSTPOT)
    {
        psv[grain_pot_idx(idx)] = force_par->grain.phiconst;
    }
    else
    {
        psv[grain_pot_idx(idx)] = 0.;
    }

    /* --------------- time ---------------- */
    psv[time_idx(idx)] = t;

    /* ------------ first step ------------- */
    inte_par->first_step = 1;       // true for the trajectory's first
                                    // integration step, false otherwise

    // write initial values to output files
    WriteCoordinates(psv, elem, output_file, force_par);


    /*-------------------------------------------------------------------------*/
    /*MAIN INTEGRATION LOOP*/
    /*-------------------------------------------------------------------------*/

    dt = time_frame->dt; /*integration time step*/

    // first step might be choosen at random to avoid artifacts in the density grid
    #ifdef dyne_random_first_step
        double dt1 = dt * gsl_rng_uniform(force_par->rng);
    #else
        double dt1 = dt;
    #endif

    if (force_par->switches.PRINT_INFO)
    {
        printf("%ld: First stepsize:  dt1 = %.15e s", force_par->run_count, dt1);
        printf("     Normal stepsize: dt  = %.15e s\n", dt);
    }

    n = 0;
    force_par->integration_error = 0;

    // reset stepsize to initial stepsize for every particle
    inte_par->h = inte_par->h1;

    while (t <= time_frame->tf)
    {
        if (n != 0)
        {
            tend = t + dt;          // subsequent step: use normal stepsize
        }
        else
        {
            tend = t + dt1;         // first step: use random stepsize if
                                    // compiled with -Ddyne_random_first_step
        }
        
        status = odeint(psv, t, tend, force_par, inte_par);

        if (status.kind == IMPACTED_SINK)
        {
            break;
        }
        else if  (status.kind == STEPSIZE_TOO_SMALL)
        {
            printf("Stepsize underflow in rkqs.\n");
            break;
        }
        else if (status.kind == TOO_MANY_STEPS)
        {
            printf("Too many steps in routine odeint.\n");
            break;
        }

        if (n != 0)
        {
            t += dt;                // subsequent step: use normal stepsize
        }
        else
        {
            t += dt1;               // first step: use random stepsize
        }

        n++;

        if (n % n_out == 0)
        {
            WriteCoordinates(psv, elem, output_file, force_par);
        }
    }

    if (!force_par->switches.output.No_Output)
        WriteSinklog(output_file, force_par, status, t);
}


unsigned long int seed_value(void)
{
    FILE* fh;
    unsigned long int r;
    unsigned int r_size = sizeof(r);

    fh = open_file("/dev/urandom", "rb");

    size_t n = fread(&r, r_size, 1, fh);
    if (n != 1)
    {
        printf("Read %lu seeds from /dev/urandom\n", n);
        exit(EXIT_FAILURE);
    }
    fclose(fh);
    return r;
}


/*-----------------------------------------------------------------------------
 --------- MAIN ROUTINE -------------------------------------------------------
 ----------------------------------------------------------------------------*/
/*****************************************************************************/
int main(int in_argc, char **in_argv)
{

    double runtime, sec;
    int min = 0;

    /*declaration of data structures*/

    /*
      HINT:
      Most of them are, after having been initialised, copied to the main
      data structure force_par. This is done in order to reduce the amount of
      structures that are passed over between the subroutines.
    */

    /* for the index functions for the phase space vector (psv) */
    psvidx_t idx;

    // Initial values structure
    iv_t iv;

    type_elements elem;
    type_inte_par inte_par;
    type_planet planet;
    type_sun sun;
    type_constants constants;
    type_time_frame time_frame;
    type_conversion conversion;
    type_switches switches;
    type_moons moons;
    type_force_par force_par;
    type_file output_file;
    type_file input_file;
    type_grain grain;


    /*-------------------------------------------------------------------------*/
    /* INITIALISE DATA STRUCTURES*/
    /*-------------------------------------------------------------------------*/

    read_jobfile(in_argc,in_argv, &input_file, &output_file, &switches,
                 &grain, &sun, &planet, &moons, &time_frame, &inte_par, &idx, &iv);

    read_constants(&constants, switches.PRINT_INFO);

    // Initialise geometric orbital elements module
    geomElemInit(planet.GM, planet.r, planet.j2, planet.j4, planet.j6);

    precalculation(&constants, &sun, &planet, &conversion, &moons, &grain,
                   &switches, &idx);

    grain.mass = 4.*val_pi * (grain.r*grain.r*grain.r) / 3. * grain.rho;

    prepare_files(&input_file, &output_file, &switches, &moons, &idx);


    // Initialise initial values structure
    if ( !iv_init(&idx, &moons, &iv) )
    {
        printf("This kind of initial values is not supported yet: pos_mode = %d\n", iv.pos_mode);
        exit(1);
    }

    // Initialise random number generator
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);      // Mersenne Twister: MT19937
    unsigned long int seed = seed_value();              // seed the random number generator
    if (switches.PRINT_INFO)
    {
        printf("RNG seed value = %lu\n", seed);
    }
    gsl_rng_set(rng, seed);
    for(int i=0; i < 1000; i++)                         // warm the rng up
    {
        gsl_rng_uniform(rng);
    }

    printf("Sun coords.x = %.15e m\n", sun.coords.x);

    /*-------------------------------------------------------------------------*/
    /*print some messages*/

    if (switches.PRINT_INFO)
    {
        MACRO_MESSAGE(("Using Runge Kutta (DOPR853) integrator\n"));

        if (time_frame.TIMESCALE == TIMESCALE_MOON)
        {
            MACRO_MESSAGE(("SCALING TIME TO SOURCE MOON YEARS\n"));
        }
        else
        {
            if (time_frame.TIMESCALE == TIMESCALE_PLANET)
                MACRO_MESSAGE(("SCALING TIME TO PLANET YEARS\n"));
            else
                MACRO_MESSAGE(("SCALING TIME TO SECONDS\n"));
        }
    }


    /*-------------------------------------------------------------------------*/
    /*time frame*/

    /*get ti and tf, i.e. starting and end integration time, in seconds*/
    if (time_frame.TIMESCALE == TIMESCALE_MOON)
    {
        /*time scale is Enceladus orbital periods*/
        time_frame.ti *= conversion.MoonYearsToSec;
        time_frame.tf *= conversion.MoonYearsToSec;
    }
    else
    {
        if (time_frame.TIMESCALE == TIMESCALE_PLANET)
        {
            /*time scale is Earth years*/
            time_frame.ti *= conversion.PlanetYearsToSec;
            time_frame.tf *= conversion.PlanetYearsToSec;
        }
        else
        {
            // time scale is seconds; this code is only needed for readability
            time_frame.ti *= 1.0;
            time_frame.tf *= 1.0;
        }
    }

    /*integration time step*/
    time_frame.dt = (time_frame.tf-time_frame.ti)/(double)time_frame.nt;

    if (switches.PRINT_INFO)
    {
        printf("ti = %f\n", time_frame.ti);
        printf("tf = %f\n", time_frame.tf);
        printf("dt = %f\n", time_frame.dt);
        /*
            Torbit, dt, timesteps/orbit
        */
    }


    /*-------------------------------------------------------------------------*/
    /*parameters for odeint integrator*/

    // calculate hmax and hmin
    inte_par.h1 *= val_2pi/moons.satellite[0].Sat_elem.om;
    inte_par.hmax *= val_2pi/moons.satellite[0].Sat_elem.om;
    inte_par.hmin *= val_2pi/moons.satellite[0].Sat_elem.om;

    // The guarantee is that the stepsize is larger than hmin.
    // To test if times differ, we need a constant by which we
    // judge if t1 == t2. This constant is stepsize_eps and we
    // choose the value here to be hmin/100.
    force_par.stepsize_eps = 0.01 * inte_par.hmin;

    // check if initial stepsize is in the range [hmin, hmax]
    if (inte_par.h1 < inte_par.hmin)
    {
        inte_par.h1 = inte_par.hmin;
    }
    else if (inte_par.h1 > inte_par.hmax)
    {
        inte_par.h1 = inte_par.hmax;
    }

    printf("h1   = %.15e\n", inte_par.h1);
    printf("hmin = %.15e\n", inte_par.hmin);
    printf("hmax = %.15e\n", inte_par.hmax);


#ifdef dyne_debug_stepsize
    fprintf(output_file.ipar, "%.15e   ", inte_par.rtol);
    fprintf(output_file.ipar, "%.15e   ", inte_par.wy);
    fprintf(output_file.ipar, "%.15e   ", inte_par.wdydt);
    fprintf(output_file.ipar, "%.15e   ", inte_par.hmin);
    fprintf(output_file.ipar, "%.15e   ", inte_par.hmax);
    fprintf(output_file.ipar, "%.15e   ", val_2pi/moons.satellite[0].Sat_elem.om);
    fprintf(output_file.ipar, "\n");
    for (int i = 0; i < idx.neq; i++)
    {
        fprintf(output_file.ipar, "%.15e   ", inte_par.atol[i]);
    }
    fprintf(output_file.ipar, "\n");
#endif


    /*-------------------------------------------------------------------------*/
    /*pass structures and parameters to the force_par structure*/

    force_par.planet = planet;
    force_par.sun = sun;
    force_par.switches = switches;
    force_par.moons = moons;
    force_par.conversion = conversion;
    force_par.constants = constants;
    force_par.grain = grain;
    force_par.time_frame = time_frame;
    force_par.RadPress = sun.RadPress;
    force_par.run_count = 1;
    force_par.integration_error = 0;
    force_par.integration_error_alternativ = 0;
    force_par.files = output_file;
    force_par.idx = &idx;
    force_par.rng = rng;

    // function to calculate the right hand derivatives (func. pointer)
    // force_par.forces = odeint_forces2;
    force_par.forces = dyne_forces2;

    // integration step routine (func. pointer)
    switch (switches.INTEGRATOR)
    {
        case 2:     inte_par.step = rkdopr853qs;
                    break;
        default:    printf("Unknown integrator type: %d", switches.INTEGRATOR);
                    exit(1);
    }

    // Choose the right function for the analytic calculation of the moon
    // phase space values.
    force_par.moons.psv = circOrbit_elem2psv;


    /*-------------------------------------------------------------------------*/
    /*CALL OF MAIN SUBROUTINE AND ADDITIONAL OUTPUT*/
    /*-------------------------------------------------------------------------*/

    /*-------------------------------------------------------------------------*/
    /*write some important physical constants and parameters to file*/

    if (switches.output.No_Output == 0)
        WriteConstants(&sun, &planet, &moons, &constants, &grain, &output_file,
                       &input_file, &idx);


    while (iv.read_vals(&iv))
    {
        OneParticle(&time_frame, &inte_par, &conversion, &force_par, &elem, &output_file, &input_file, &idx, &iv);
        force_par.run_count = force_par.run_count + 1;
    }



    /*-------------------------------------------------------------------------*/
    /*CLOSE FILE STREAMS, FREE MEMORY  AND EXIT*/
    /*-------------------------------------------------------------------------*/

    /*-------------------------------------------------------------------------*/
    /*close input/output streams*/

    fclose(input_file.ID);          // jobfile
    fclose(iv.file);                // initial values

    if (switches.output.PSpaceOut)
        fclose(output_file.ODpsc);

    if (switches.output.elemout)
        fclose(output_file.ODoscel);

    if (switches.CURRENTS == 1)
    {
        if (switches.output.No_Output == 0)
            fclose(output_file.ODphi);

        if (switches.output.PlasmaParOut)
        {
            fclose(output_file.ODplasma_par);
            fclose(output_file.ODbfield);
        }
    }

    if (switches.output.LorentzOut)
        fclose(output_file.ODLorentz);

    if (switches.output.SatCoordOut)
        for (int j = moon_start_idx(); exist_moon(&idx, j); j++)
            fclose(output_file.ODsatCoord[moon_array_idx(&idx, j)]);

    if (switches.output.SatElemOut)
        for (int j = moon_start_idx(); exist_moon(&idx, j); j++)
            fclose(output_file.ODsatElem[moon_array_idx(&idx, j)]);

    if (switches.output.combined)
        fclose(output_file.ODcombined);

    if (switches.output.sun)
        fclose(output_file.ODsun);

    if (switches.output.TestsOut)
        fclose(output_file.ODTestsOut);

    if (switches.output.JacobiConstant)
        fclose(output_file.ODJacobiConstant);

    if (switches.output.TrafoCoord)
        fclose(output_file.ODTrafoCoord);

    if (switches.output.RelCoord)
        fclose(output_file.ODRelCoord);

    if (!switches.output.No_Output)
        fclose(output_file.ODs);

    fclose(output_file.ODstartparam);

#ifdef dyne_debug_stepsize
    fclose(output_file.ipar);
    fclose(output_file.stepsize_diag);
#endif



    /*-------------------------------------------------------------------------*/
    /*free allocated memory*/

    free(moons.satellite);
    free(inte_par.atol);


    if (switches.output.SatCoordOut)
        free(output_file.ODsatCoord);
    if (switches.output.SatElemOut)
        free(output_file.ODsatElem);


    // free all memory used for the random number generator
    gsl_rng_free(rng);


    /*-------------------------------------------------------------------------*/
    /*calculate runtime of the programme*/

    if (switches.PRINT_INFO)
    {
        runtime = clock()/(60.*CLOCKS_PER_SEC);
        while(min < runtime-1)
            min++;
        sec = (runtime - min)*60;

        printf("\nIntegration runtime: %d min : %2.0lf s\n\n",min,sec);
    }


    return 0;

}
