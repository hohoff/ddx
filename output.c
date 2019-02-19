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

/*-----------------------------------------------------------------------------
 --------- WriteCoordinates ---------------------------------------------------
 ----------------------------------------------------------------------------*/
/*****************************************************************************/
void WriteCoordinates(double psv[], type_elements *elem, type_file *output_file, type_force_par *force_par)
{

    /*main output routine*/

    long arraylength;
    size_t size;

    double t,x,y,z,u,v,w;
    double tscale;

    psvidx_t *idx = force_par->idx;

    type_coordinates coord = {0., 0., 0., 0., 0., 0., 0.};

    /*dust grain coordinates/velocity*/
    coord.t = psv[time_idx(idx)];
    t = coord.t;

    coord.x = psv[grain_pos_idx(idx, 0)];
    coord.y = psv[grain_pos_idx(idx, 1)];
    coord.z = psv[grain_pos_idx(idx, 2)];
    coord.u = psv[grain_vel_idx(idx, 0)];
    coord.v = psv[grain_vel_idx(idx, 1)];
    coord.w = psv[grain_vel_idx(idx, 2)];
    
    
    /*dust grain potential*/
    if (force_par->switches.CURRENTS)
        force_par->grain.phi_calc = psv[grain_pot_idx(idx)];

    /*scaling output time*/

    if (force_par->switches.output.EARTHTIMESCALE)
    {
        tscale = t*force_par->conversion.SecToEarthYears;       // Earth years
    }
    else
    {
        if (force_par->time_frame.TIMESCALE == TIMESCALE_MOON)
        {
            tscale = t*force_par->conversion.SecToMoonYears;        // Moon orbital periods
        }
        else
        {
            if (force_par->time_frame.TIMESCALE == TIMESCALE_PLANET)
            {
                tscale = t*force_par->conversion.SecToPlanetYears;  // Saturnian years
            }
            else
            {
                tscale = t;     // don't scale, just seconds
            }
        }
    }


    /*-------------------------------------------------------------------------*/
    /*output phase space data of dust grain*/

    if (force_par->switches.output.PSpaceOut)
    {
        x = coord.x;
        y = coord.y;
        z = coord.z;
        u = coord.u;
        v = coord.v;
        w = coord.w;

        if (force_par->switches.output.filemod == OUTPUT_BINARY)
        {

            /*Binary output*/

            arraylength = 1;
            size = sizeof(double);
            fwrite(&t,size,arraylength,output_file->ODpsc);
            fwrite(&tscale,size,arraylength,output_file->ODpsc);
            fwrite(&x,size,arraylength,output_file->ODpsc);
            fwrite(&y,size,arraylength,output_file->ODpsc);
            fwrite(&z,size,arraylength,output_file->ODpsc);
            fwrite(&u,size,arraylength,output_file->ODpsc);
            fwrite(&v,size,arraylength,output_file->ODpsc);
            fwrite(&w,size,arraylength,output_file->ODpsc);
        }

        else
        {

            /*ASCII output*/

            fprintf(output_file->ODpsc,
                    "%21.15e %21.15e %21.15e %21.15e %21.15e %21.15e %21.15e "
                    "%21.15e\n",t,tscale,x,y,z,u,v,w);
        }
    }




    /*-------------------------------------------------------------------------*/
    /*output phase space data and orbital elements of included satellites*/

    if (force_par->switches.output.SatCoordOut)
    {
        for ( int j = moon_start_idx(); exist_moon(idx, j); j++)
        {
            x = psv[moon_pos_idx(idx, j, 0)];
            y = psv[moon_pos_idx(idx, j, 1)];
            z = psv[moon_pos_idx(idx, j, 2)];
            u = psv[moon_vel_idx(idx, j, 0)];
            v = psv[moon_vel_idx(idx, j, 1)];
            w = psv[moon_vel_idx(idx, j, 2)];
    
            if (force_par->switches.output.filemod == OUTPUT_BINARY)
            {
                arraylength = 1;
                size = sizeof(double);
                fwrite(&t,size,arraylength,output_file->ODsatCoord[moon_array_idx(idx, j)]);
                fwrite(&tscale,size,arraylength,output_file->ODsatCoord[moon_array_idx(idx, j)]);
                fwrite(&x,size,arraylength,output_file->ODsatCoord[moon_array_idx(idx, j)]);
                fwrite(&y,size,arraylength,output_file->ODsatCoord[moon_array_idx(idx, j)]);
                fwrite(&z,size,arraylength,output_file->ODsatCoord[moon_array_idx(idx, j)]);
                fwrite(&u,size,arraylength,output_file->ODsatCoord[moon_array_idx(idx, j)]);
                fwrite(&v,size,arraylength,output_file->ODsatCoord[moon_array_idx(idx, j)]);
                fwrite(&w,size,arraylength,output_file->ODsatCoord[moon_array_idx(idx, j)]);
            }
    
            else
            {
                fprintf(output_file->ODsatCoord[moon_array_idx(idx, j)],"%21.15e %21.15e %21.15e "
                        "%21.15e %21.15e %21.15e %21.15e %21.15e\n",
                        t,tscale,x,y,z,u,v,w);
            }
        }
    }
    

    if (force_par->switches.output.combined)
    {
        int ref_moon = force_par->moons.ref_moon;
        fprintf(output_file->ODcombined,
                   "%.15e %.15e   %.15e %.15e %.15e   %.15e %.15e %.15e\n",
                      t, tscale,
                      psv[grain_pos_idx(idx, 0)],
                      psv[grain_pos_idx(idx, 1)],
                      psv[grain_pos_idx(idx, 2)],
                      psv[moon_pos_idx(idx, ref_moon, 0)],
                      psv[moon_pos_idx(idx, ref_moon, 1)],
                      psv[moon_pos_idx(idx, ref_moon, 2)]);
    }


    if (force_par->switches.output.sun)
    {
        fprintf(output_file->ODsun,
                   "%.15e %.15e   %.15e %.15e %.15e   %.15e %.15e %.15e\n",
                      t, tscale,
                      psv[sun_pos_idx(idx, 0)],
                      psv[sun_pos_idx(idx, 1)],
                      psv[sun_pos_idx(idx, 2)],
                      psv[sun_vel_idx(idx, 0)],
                      psv[sun_vel_idx(idx, 1)],
                      psv[sun_vel_idx(idx, 2)]);
    }

}


/*-----------------------------------------------------------------------------
 --------- WriteConstants -----------------------------------------------------
 ----------------------------------------------------------------------------*/
/*****************************************************************************/
void WriteConstants(type_sun *sun, type_planet *planet, type_moons *moons,
                    type_constants *constants, type_grain *grain,
                    type_file *output_file, type_file *input_file, 
                    psvidx_t *idx)
{

    double GM_Sun, GM_Sat, J0, dS, RadPress, R_Sat, J2, g10, om_Sat, Q_m_perPhi;
    double phi;
    double GM_M, n_M, T_M, h_M, v_esc_M, v_orbit_M, pom_J2_M, pom_EM_M, gamma;

    FILE *ConstOut;
    char OUTFILE[181];
    char jobname[181];
    char aux1[181];
    char aux2[11];
    int i, n;


    /*open output file stream "constants.txt"*/

    i = strlen(input_file->name);
    strncpy(jobname,input_file->name, (i-4));
    jobname[i-4] = '\0';

    strcpy(OUTFILE,output_file->dir);
    strcat(OUTFILE,"/");
    strcat(OUTFILE,jobname);
    strcat(OUTFILE,"/constants.txt");

    ConstOut = fopen(OUTFILE,"w");
    MACRO_ERROR(ConstOut == NULL,
                ("CANNOT OPEN OUTPUT FILE (.txt) FOR WRITING\n"));


    /*constants and parameters for output*/

    GM_Sun = sun->GM;        /*mu of the Sun*/
    GM_Sat = planet->GM;     /*mu of Saturn*/
    J0 = sun->J0;            /*solar constant at 1 AU*/

    /*orbit averaged distance of Saturn from the Sun in AU*/
    dS = planet->elem.a/constants->au*(1 + 0.5*pow(planet->elem.e,2));

    RadPress = sun->RadPress;/*acceleration due to radiation pressure*/
    R_Sat = planet->r;       /*radius of Saturn*/
    J2 = planet->j2;         /*J2 of Saturn*/
    g10 = planet->g10;       /*g10 of Saturn's magnetic field*/
    om_Sat = planet->om;     /*spin rate of Saturn*/
    phi = grain->phiconst;   /*constant grain potential*/

    /*Q/m per potential*/
    Q_m_perPhi = 3*constants->eps0/(grain->rho*grain->r*grain->r);


    /*write output*/

    fprintf(ConstOut,"Constants for jobfile: %s\n\n",jobname);

    fprintf(ConstOut,
            "mu of the Sun:\t\t\t%21.15e\nmu of Saturn:\t\t\t%21.15e\n"
            "Solar constant:\t\t\t%21.15e\nDistance Sun-Saturn [AU]:\t%21.15e\n"
            "Radiation pressure acc.:\t%21.15e\nSaturn radius:\t\t\t%21.15e\n"
            "J2 of Saturn:\t\t\t%21.15e\ng10 of B-field:\t\t\t%21.15e\n"
            "Saturn spin rate:\t\t%21.15e\nQ/m per potential:\t\t%21.15e\n"
            "Const. grain potential:\t\t%4.1lf\n\n",
            GM_Sun,GM_Sat,J0,dS,RadPress,R_Sat,J2,g10,om_Sat,Q_m_perPhi,phi);


    /*write properties of included moons*/

    n = 0;

    for (int j = moon_start_idx(); exist_moon(idx, j); j++)
    {
        n++;

        /*constants of jth moon*/

        GM_M = moons->satellite[moon_array_idx(idx, j)].GM;         /*mu of moon*/
        n_M = moons->satellite[moon_array_idx(idx, j)].Sat_elem.om; /*Kepler frequency*/
        T_M = moons->satellite[moon_array_idx(idx, j)].year;        /*orbital period*/
        h_M = moons->satellite[moon_array_idx(idx, j)].h;           /*Hill radius*/
        v_esc_M = moons->satellite[moon_array_idx(idx, j)].vesc;    /*escape speed*/

        /*orbital speed*/
        v_orbit_M = 1./(moons->satellite[moon_array_idx(idx, j)].MeterPerSecToOrbitSpeed);

        /*precession rate due to J2*/
        pom_J2_M = 1.5*n_M*planet->j2*pow(planet->r,2)/
                   pow(moons->satellite[moon_array_idx(idx, j)].Sat_elem.a,2);

        /*precession rate due to EM*/
        pom_EM_M = 2*Q_m_perPhi*grain->phiconst*g10*
                   pow(R_Sat/moons->satellite[moon_array_idx(idx, j)].Sat_elem.a,3);

        /*total precession rate*/
        gamma = pom_J2_M + pom_EM_M;

        strcpy(aux1,"Satellite no. ");
        sprintf(aux2,"%d",n);
        strcat(aux1,aux2);
        strcat(aux1,":\n\n");

        fprintf(ConstOut,"%s",aux1);
        fprintf(ConstOut,
                "mu of the moon:\t\t\t%21.15e\nKepler frequency:\t\t%21.15e\n"
                "orbital period:\t\t\t%21.15e\nHill radius:\t\t\t%21.15e\n"
                "escape speed:\t\t\t%21.15e\norbital speed:\t\t\t%21.15e\n"
                "J2 precession rate:\t\t%21.15e\n"
                "EM precession rate* (%4.1lf V):\t%21.15e\n"
                "total precession rate gamma*:\t%21.15e\n"
                "*: valid for particle located at position of this moon!\n\n",
                GM_M,n_M,T_M,h_M,v_esc_M,v_orbit_M,pom_J2_M,phi,pom_EM_M,gamma);
    }


    /*close file*/

    fclose(ConstOut);

}


void WriteSinklog(type_file *output_file, type_force_par *fpar, integ_data_t status, double t)
{
    char sinkstatus[181];

    if (status.kind == IMPACTED_SINK)
    {
        strcpy(sinkstatus,"Dust grain impacted ");
        strcat(sinkstatus, status.impact.sink_name);
        strcat(sinkstatus,".");

        t = status.impact.time;
    }
    else if  (status.kind == STEPSIZE_TOO_SMALL)
    {
        strcpy(sinkstatus,"Stepsize underflow");
    }
    else if (status.kind == TOO_MANY_STEPS)
    {
        strcpy(sinkstatus,"Too many steps in routine odeint.");
    }
    else
    {
        strcpy(sinkstatus,"Dust grain survived.");
    }

    double tscale = -1.;

    if (fpar->switches.output.EARTHTIMESCALE)
    {
        tscale = t*fpar->conversion.SecToEarthYears;       // Earth years
    }
    else
    {
        if (fpar->time_frame.TIMESCALE == TIMESCALE_MOON)
        {
            tscale = t*fpar->conversion.SecToMoonYears;        // Moon orbital periods
        }
        else
        {
            if (fpar->time_frame.TIMESCALE == TIMESCALE_PLANET)
            {
                tscale = t*fpar->conversion.SecToPlanetYears;  // Saturnian years
            }
            else
            {
                tscale = t;     // don't scale, just seconds
            }
        }
    }

    fprintf(output_file->ODs, "%ld     %21.15e     %21.15e     %s\n",
                                    fpar->run_count, t, tscale, sinkstatus);
}


