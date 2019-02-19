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
#include "jobfile.h"
#include "psvindex.h"
#include "initialvalues.h"

static const long int expected_format_version = 27022017;

/*-----------------------------------------------------------------------------
 --------- open files with proper error messages ------------------------------
 ----------------------------------------------------------------------------*/
/*****************************************************************************/
FILE* open_file(char *filename, char *mode)
{
    errno = 0;
    FILE* fh;
    fh = fopen(filename, mode);
    if (fh == NULL)
    {
        printf("Could'nt open file '%s': %s\n", filename, strerror(errno)); 
        exit(EXIT_FAILURE);
    }
    return fh;
}


/*-----------------------------------------------------------------------------
 --------- STRUCTURES FOR PHYSICAL CONSTANTS AND DATA -------------------------
 ----------------------------------------------------------------------------*/
/*****************************************************************************/
void read_constants(type_constants *constants, int verbose)
{
    char apu[181];
    char varname[181];
    long int formatversion;

    FILE* const_fh = open_file("constants.conf", "r");

    fgets(apu,80,const_fh);
    sscanf(apu,"%s %ld",varname,&formatversion);
    if (formatversion != expected_format_version)
    {
        printf("Unexpected format of jobfile (%s:%d).\n", __FILE__, __LINE__);
        printf("   expected format version: %ld\n", expected_format_version); 
        printf("    jobfile format version: %ld\n", formatversion); 
        exit(EXIT_FAILURE);
    }

    read_option(const_fh, REAL, "G", &constants->G);        // gravitational constant [m3 kg-1 s-2]
    read_option(const_fh, REAL, "c", &constants->c);        // speed of light [m/s]
    read_option(const_fh, REAL, "eps0", &constants->eps0);  // vacuum permittivity [A s V-1 m-1]
    read_option(const_fh, REAL, "qe", &constants->qe);      // elementary charge [C]
    read_option(const_fh, REAL, "me", &constants->me);      // electron mass [kg]
    constants->qm = constants->qe/constants->me;
    read_option(const_fh, REAL, "mp", &constants->mp);      // proton mass [kg]
    read_option(const_fh, REAL, "au", &constants->au);      //astronomical unit [m]
    read_option(const_fh, REAL, "amu", &constants->amu);    // atomic mass unit [kg]

    if(verbose)
    {
        printf("G = %.15e m^3 kg^-1 s^-2\n", constants->G);
        printf("c = %.15e m s^-1\n", constants->c);
        printf("eps0 = %.15e A s V^-1 m^-1 \n", constants->eps0);
        printf("qe = %.15e C\n", constants->qe);
        printf("me = %.15e kg\n", constants->me);
        printf("qm = %.15e C kg^-1\n", constants->qm);
        printf("mp = %.15e kg\n", constants->mp);
        printf("au = %.15e m\n", constants->au);
        printf("amu = %.15e kg\n", constants->amu);
    }

    fclose(const_fh);
}


/*****************************************************************************/
void precalculation(type_constants *constants, type_sun *sun,
                    type_planet *planet, type_conversion *conversion,
                    type_moons *moons, type_grain *grain,
                    type_switches *switches, psvidx_t *idx)
{

    /*fill Solar System data structures*/

    // double dS;
    double deg2rad;
    

    /*---ANGLE CONVERSION---*/

    /*convert all input angles in JOBFILE from degree to radian*/

    conversion->rad2deg = 180./val_pi;
    deg2rad = 1./(conversion->rad2deg);

    /*central planet (Saturn)*/
    planet->elem.i *= deg2rad;
    planet->elem.aop *= deg2rad;
    planet->elem.lan *= deg2rad;
    planet->obliquity *= deg2rad;

    /*satellites*/
    for (int j = moon_start_idx(); exist_moon(idx, j); j++)
    {
        moons->satellite[moon_array_idx(idx, j)].Sat_elem.i *= deg2rad;
        moons->satellite[moon_array_idx(idx, j)].Sat_elem.aop *= deg2rad;
        moons->satellite[moon_array_idx(idx, j)].Sat_elem.lan *= deg2rad;
    }

    /*grain launching location*/
    // grain->phi *= deg2rad;
    // grain->theta *= deg2rad;
    // grain->alpha_min *= deg2rad;
    // grain->alpha_max *= deg2rad;
    // grain->beta_min *= deg2rad;
    // grain->beta_max *= deg2rad;


    /*---SUN---*/

    /*orbit averaged distance of Saturn to the Sun in AU*/
    // dS = planet->elem.a/constants->au*(1 + 0.5*pow(planet->elem.e,2.));
    /*coefficient of radiation pressure force*/
    // printf("# Old acceleration due to radiation pressure: %.15e\n", 0.75*sun->J0*grain->Qpr/(grain->rho*grain->r*constants->c*pow(dS,2.)));

    // coefficient of radiation pressure force
    sun->RadPress = 0.75 * sun->J0 * grain->Qpr / (grain->rho * grain->r * constants->c);


    /*---PLANET -> SATURN---*/

    planet->elem.om = sqrt((sun->GM+planet->GM)/pow(planet->elem.a,3));
    /*Kepler frequency of planet*/
    planet->year = val_2pi/planet->elem.om; /*planetary year in seconds*/
    planet->cosobl = cos(planet->obliquity);
    planet->sinobl = sin(planet->obliquity);


    /*---MOONS---*/

    for (int j = moon_start_idx(); exist_moon(idx, j); j++)
    {
        // geometric orbital elements of the moon
        Coord2geomElem(0.0,
                       moons->satellite[moon_array_idx(idx, j)].GM,
                       &moons->satellite[moon_array_idx(idx, j)].Sat_coord,
                       &moons->satellite[moon_array_idx(idx, j)].Sat_elem);

        // mean motion of the moon
        double n = 0.0;
        double kappa = 0.0;
        double nu = 0.0;
        geomElemFreq(moons->satellite[moon_array_idx(idx, j)].GM, &moons->satellite[moon_array_idx(idx, j)].Sat_elem, &n, &kappa, &nu);
        moons->satellite[moon_array_idx(idx, j)].Sat_elem.om = n;

        // orbital period of the moon
        moons->satellite[moon_array_idx(idx, j)].year = val_2pi/moons->satellite[moon_array_idx(idx, j)].Sat_elem.om;

        // Hill radius of the moon
        moons->satellite[moon_array_idx(idx, j)].h =
            moons->satellite[moon_array_idx(idx, j)].Sat_elem.a*
            pow(moons->satellite[moon_array_idx(idx, j)].GM/planet->GM/3.,1./3.);

        // two-body escape velocity from moon
        moons->satellite[moon_array_idx(idx, j)].vesc =
            sqrt(2.*moons->satellite[moon_array_idx(idx, j)].GM/moons->satellite[moon_array_idx(idx, j)].r);

        // orbit averaged distance moon to planet
        moons->satellite[moon_array_idx(idx, j)].r_mean = moons->satellite[moon_array_idx(idx, j)].Sat_elem.a*
                                     (1. + 0.5*pow(moons->satellite[moon_array_idx(idx, j)].Sat_elem.e,2.));

        // conversion from m/s to orbital speed
        moons->satellite[moon_array_idx(idx, j)].MeterPerSecToOrbitSpeed =
            1./(moons->satellite[moon_array_idx(idx, j)].Sat_elem.om*moons->satellite[moon_array_idx(idx, j)].r_mean);
    }


    /*---OTHER CONVERSIONS - in the following order:---*/

    /*Seconds to Saturn years (orbital period)
      Saturn years to seconds
      Meter to Saturnian radii
      Seconds to Enceladus orbital periods
      Enceladus orbital periods to seconds
      Seconds to Earth years*/

    conversion->SecToPlanetYears = 1./planet->year;
    conversion->PlanetYearsToSec = planet->year;
    conversion->MeterToPlanetRadii = 1./planet->r;
    conversion->SecToMoonYears = 1./moons->satellite[moons->sourceSat-1].year;
    conversion->MoonYearsToSec = moons->satellite[moons->sourceSat-1].year;
    conversion->SecToEarthYears = 1./31558118.4;


}


/*-----------------------------------------------------------------------------
 --------- READ THE JOBILE ----------------------------------------------------
 ----------------------------------------------------------------------------*/
/*****************************************************************************/
void read_jobfile(int argc, char **argv, type_file *input_file,
                  type_file *output_file, type_switches *switches,
                  type_grain *grain, type_sun *sun, type_planet *planet,
                  type_moons *ptr_moons,type_time_frame *time_frame,
                  type_inte_par *ipar, psvidx_t *idx, iv_t *iv)
{

    /*----------------------------------------------------*/
    /*read input parameters from jobfile*/
    /*
      Input:
      argc: number of command line arguments
      argv: string vector containing command line arguments
            argv[0]: "dd.x"
      argv[1]: name of JOBFILE
      + all necessary data strucures
    */
    /*----------------------------------------------------*/
    // printf("expected_format_version: %ld\n", expected_format_version);

    char apu[181];
    char varname[181];

    long int formatversion;
    int prnt;

    type_satellite *satellites;


    /*=================================*/
    /*-----FILES AND DIRECTOIRES:------*/
    /*=================================*/

    /*error checking*/

    if (argc < 2)
    {
        printf("USAGE: dd.x JOBFILE \n");
        exit(-1);
    }

    strcpy(input_file->name,argv[1]);

    /*if missing, add a ".job" extension to read the file*/

    if (strstr(input_file->name,".job") == NULL)
    {
        strcat(input_file->name,".job");
    }

    /*open jobfile*/

    input_file->ID=fopen(input_file->name,"r");
    MACRO_ERROR(input_file->ID==NULL,("\nCANNOT OPEN .job FILE FOR READING\n"));


    /*=================================*/
    /*-----READ DATA FROM JOBFILE:-----*/
    /*=================================*/


    /***************************************************************************/
    /******************** FORMATVERSION **************************************/
    /***************************************************************************/

    fgets(apu,80,input_file->ID);
    sscanf(apu,"%s %ld",varname,&formatversion);

    if (formatversion != expected_format_version)
    {
        printf("Unexpected format of jobfile (%s:%d).\n", __FILE__, __LINE__);
        printf("   expected format version: %ld\n", expected_format_version); 
        printf("    jobfile format version: %ld\n", formatversion); 
        exit(EXIT_FAILURE);
    }

    /***************************************************************************/
    /******************** SWITCHES *********************************************/
    /***************************************************************************/

    read_option(input_file->ID, INTEGER, "LORENTZ", &switches->LORENTZ);
    read_option(input_file->ID, INTEGER, "LORENTZ2", &switches->LORENTZ2);
    read_option(input_file->ID, INTEGER, "CURRENTS", &switches->CURRENTS);
    read_option(input_file->ID, INTEGER, "CONSTPOT", &switches->CONSTPOT);
    read_option(input_file->ID, INTEGER, "EQUIPOT", &switches->EQUIPOT);
    read_option(input_file->ID, INTEGER, "SINKS", &switches->SINKS);
    read_option(input_file->ID, INTEGER, "RADPRESS", &switches->RADPRESS);
    read_option(input_file->ID, INTEGER, "INTEGRATOR", &switches->INTEGRATOR);
    read_option(input_file->ID, INTEGER, "INTERPOLATION", &switches->INTERPOLATION);


    if (!switches->CURRENTS && !switches->CONSTPOT && !switches->EQUIPOT)
    {
        printf("Please choose exactly one kind of dustgrain potential calculation\n"
               "in the jobfile (CURRENTS, CONSTPOT or EQUIPOT).\n");
        exit(EXIT_FAILURE);
    }
    else if (( switches->CURRENTS &&  switches->CONSTPOT &&  switches->EQUIPOT) ||
             ( switches->CURRENTS &&  switches->CONSTPOT && !switches->EQUIPOT) ||
             ( switches->CURRENTS && !switches->CONSTPOT &&  switches->EQUIPOT) ||
             (!switches->CURRENTS &&  switches->CONSTPOT &&  switches->EQUIPOT)   )
    {
        printf("Please choose only one kind of dustgrain potential calculation\n"
               "in the jobfile (CURRENTS, CONSTPOT or EQUIPOT).\n");
        exit(EXIT_FAILURE);
    }

    MACRO_ERROR(((switches->INTERPOLATION != POLYNOMIAL) &&
                 (switches->INTERPOLATION != CUBIC_SPLINE)),
                ("\nUnknown choice of interpolation scheme "
                 "(choose 1 or 2 only)\n"));

    MACRO_ERROR((switches->INTERPOLATION == POLYNOMIAL),
                ("\nProgramme is no longer intended to be used with polynomial "
                 "interpolation! Use CUBIC SPLINE interpolation!\n"));


    read_option(input_file->ID, INTEGER, "PLASMA_MODEL", &switches->PLASMA_MODEL);
    read_option(input_file->ID, INTEGER, "DIRECT_DRAG", &switches->DIRECT_DRAG);
    // read_option(input_file->ID, INTEGER, "PLASMA_DRAG", &switches->PLASMA_DRAG);
    read_option(input_file->ID, INTEGER, "PRINT_INFO", &switches->PRINT_INFO);

    prnt = switches->PRINT_INFO;

    if(prnt)
    {
        printf("input_file->name = %s\n",input_file->name);
        printf("formatversion = %ld\n",formatversion);

        printf("LORENTZ = %d\n",switches->LORENTZ);
        printf("LORENTZ2 = %d\n",switches->LORENTZ2);
        printf("CURRENTS = %d\n",switches->CURRENTS);
        printf("CONSTPOT = %d\n",switches->CONSTPOT);
        printf("EQUIPOT = %d\n",switches->EQUIPOT);
        printf("SINKS = %d\n",switches->SINKS);
        printf("RADPRESS = %d\n",switches->RADPRESS);
        printf("INTEGRATOR = %d\n",switches->INTEGRATOR);
        printf("INTERPOLATION = %d\n",switches->INTERPOLATION);
        printf("PLASMA MODEL = %d\n",switches->PLASMA_MODEL);
        printf("DIRECT DRAG = %d\n",switches->DIRECT_DRAG);
        // printf("PLASMA DRAG = %d\n",switches->PLASMA_DRAG);
        printf("PRINT INFOS = %d\n",switches->PRINT_INFO);
    }

    /***************************************************************************/
    /******************** GRAIN PROPERTIES *************************************/
    /***************************************************************************/

    read_option(input_file->ID, REAL, "r", &grain->r);
    read_option(input_file->ID, REAL, "rho", &grain->rho);
    read_option(input_file->ID, REAL, "dm", &grain->dm);
    read_option(input_file->ID, REAL, "em", &grain->em);
    read_option(input_file->ID, REAL, "kappa", &grain->kappa);
    read_option(input_file->ID, REAL, "kTs", &grain->kTs);
    read_option(input_file->ID, REAL, "kTn", &grain->kTn);
    read_option(input_file->ID, REAL, "phiconst", &grain->phiconst);
    read_option(input_file->ID, REAL, "Qpr", &grain->Qpr);

    if(prnt)
    {
        printf("r = %le\n", grain->r);
        printf("rho = %le\n", grain->rho);
        printf("dm = %le\n", grain->dm);
        printf("em = %le\n", grain->em);
        printf("kappa = %le\n", grain->kappa);
        printf("kTs = %le\n", grain->kTs);
        printf("kTn = %le\n", grain->kTn);
        printf("phiconst = %le\n", grain->phiconst);
        printf("Qpr = %le\n", grain->Qpr);
    }

    /***************************************************************************/
    /******************** SUN **************************************************/
    /***************************************************************************/

    read_option(input_file->ID, REAL, "GM", &sun->GM);
    read_option(input_file->ID, REAL, "J0", &sun->J0);
    read_option(input_file->ID, REAL, "x", &sun->coords.x);
    read_option(input_file->ID, REAL, "y", &sun->coords.y);
    read_option(input_file->ID, REAL, "z", &sun->coords.z);
    read_option(input_file->ID, REAL, "u", &sun->coords.u);
    read_option(input_file->ID, REAL, "v", &sun->coords.v);
    read_option(input_file->ID, REAL, "w", &sun->coords.w);

    if(prnt)
    {
        printf("Sun GM = %.15e\n", sun->GM);
        printf("Sun J0 = %.15e\n", sun->J0);
        printf("Sun x0 = %.15e\n", sun->coords.x);
        printf("Sun y0 = %.15e\n", sun->coords.y);
        printf("Sun z0 = %.15e\n", sun->coords.z);
        printf("Sun u0 = %.15e\n", sun->coords.u);
        printf("Sun v0 = %.15e\n", sun->coords.v);
        printf("Sun w0 = %.15e\n", sun->coords.w);
    }

    /***************************************************************************/
    /******************** PLANET ************************************************/
    /***************************************************************************/

    read_option(input_file->ID, STRING, "name", planet->name);
    read_option(input_file->ID, REAL, "r", &planet->r);
    read_option(input_file->ID, REAL, "a", &planet->elem.a);
    read_option(input_file->ID, REAL, "e", &planet->elem.e);
    read_option(input_file->ID, REAL, "i", &planet->elem.i);
    read_option(input_file->ID, REAL, "aop", &planet->elem.aop);
    read_option(input_file->ID, REAL, "lan", &planet->elem.lan);
    read_option(input_file->ID, REAL, "Tper", &planet->elem.Tper);
    read_option(input_file->ID, REAL, "om", &planet->om);
    read_option(input_file->ID, REAL, "GM", &planet->GM);
    read_option(input_file->ID, REAL, "obliquity", &planet->obliquity);
    read_option(input_file->ID, REAL, "j2", &planet->j2);
    read_option(input_file->ID, REAL, "j4", &planet->j4);
    read_option(input_file->ID, REAL, "j6", &planet->j6);
    read_option(input_file->ID, REAL, "g10", &planet->g10);
    read_option(input_file->ID, REAL, "f", &planet->f);
    read_option(input_file->ID, REAL, "n0", &planet->n0);
    // read_option(input_file->ID, REAL, "fpd", &planet->fpd);

    if(prnt)
    {
        printf("name = %s\n",  planet->name);
        printf("r = %le\n", planet->r);
        printf("a = %le\n", planet->elem.a);
        printf("e = %le\n", planet->elem.e);
        printf("i = %le\n", planet->elem.i);
        printf("aop = %le\n", planet->elem.aop);
        printf("lan = %le\n", planet->elem.lan);
        printf("Tper = %le\n", planet->elem.Tper);
        printf("om = %le\n", planet->om);
        printf("GM = %le\n", planet->GM);
        printf("obliquity = %le\n", planet->obliquity);
        printf("J2 = %le\n", planet->j2);
        printf("J4 = %le\n", planet->j4);
        printf("J6 = %le\n", planet->j6);
        printf("g10 = %le\n", planet->g10);
        printf("f = %le\n", planet->f);
        printf("n0 = %le\n", planet->n0);
        // printf("fpd = %le\n", planet->fpd);
    }

    /***************************************************************************/
    /******************** SATELLITES *******************************************/
    /***************************************************************************/

    read_option(input_file->ID, INTEGER, "satN", &ptr_moons->satN);
    read_option(input_file->ID, INTEGER, "sourceSat", &ptr_moons->sourceSat);
    read_option(input_file->ID, INTEGER, "ref_moon", &ptr_moons->ref_moon);

    /* error checking */
    MACRO_ERROR((ptr_moons->sourceSat == 0),
                ("\nSOURCE SATELLITE NEEDED FOR CALCULATING "
                 "INITIAL DUST GRAIN POSITION!\n"));

    if (ptr_moons->sourceSat < 1 || ptr_moons->sourceSat > ptr_moons->satN)
    {
        printf("The source moon number has to be between 1 and %d. Source moon number is %d\n", ptr_moons->satN, ptr_moons->sourceSat); 
        exit(EXIT_FAILURE);
    }

    if (ptr_moons->ref_moon < 1 || ptr_moons->ref_moon > ptr_moons->satN)
    {
        printf("The reference moon number has to be between 1 and %d. Reference moon number is %d\n", ptr_moons->satN, ptr_moons->ref_moon); 
        exit(EXIT_FAILURE);
    }

    /* allocate memory for the moons structure */
    satellites =
        (type_satellite*)malloc((ptr_moons->satN+1)*sizeof(type_satellite));
    MACRO_ERROR(satellites==NULL,("\nCANNOT ALLOCATE MEMORY\n"));

    if(prnt)
    {
        printf("satN = %d\n", ptr_moons->satN);
        printf("sourceSat = %d\n", ptr_moons->sourceSat);
        printf("ref_moon = %d\n", ptr_moons->ref_moon);
    }

    /***************************************************************************/
    /******************** j-th SATELLITE ***************************************/
    /***************************************************************************/

    idx_initialise(idx, ptr_moons->satN, switches->CURRENTS);
    /*loop over the satN moons to fill each satellite structure with
      the respective moon data*/
    for (int j = moon_start_idx(); exist_moon(idx, j); j++)
    {
        read_option(input_file->ID, STRING, "name", &satellites[moon_array_idx(idx, j)].name);
        read_option(input_file->ID, REAL, "r", &satellites[moon_array_idx(idx, j)].r);
        read_option(input_file->ID, REAL, "GM", &satellites[moon_array_idx(idx, j)].GM);
        read_option(input_file->ID, INTEGER, "sink", &satellites[moon_array_idx(idx, j)].sink);
        read_option(input_file->ID, REAL, "x", &satellites[moon_array_idx(idx, j)].Sat_coord.x);
        read_option(input_file->ID, REAL, "y", &satellites[moon_array_idx(idx, j)].Sat_coord.y);
        read_option(input_file->ID, REAL, "z", &satellites[moon_array_idx(idx, j)].Sat_coord.z);
        read_option(input_file->ID, REAL, "u", &satellites[moon_array_idx(idx, j)].Sat_coord.u);
        read_option(input_file->ID, REAL, "v", &satellites[moon_array_idx(idx, j)].Sat_coord.v);
        read_option(input_file->ID, REAL, "w", &satellites[moon_array_idx(idx, j)].Sat_coord.w);
    }


    /*copy satellite structure to structure given by calling routine*/
    ptr_moons->satellite = satellites;

    if(prnt)
    {
        for (int j = moon_start_idx(); exist_moon(idx, j); j++)
        {
            printf("name = %s\n", ptr_moons->satellite[moon_array_idx(idx, j)].name);
            printf("r = %le\n",ptr_moons->satellite[moon_array_idx(idx, j)].r);
            printf("GM = %le\n",ptr_moons->satellite[moon_array_idx(idx, j)].GM);
            printf("sink = %d\n", ptr_moons->satellite[moon_array_idx(idx, j)].sink);
            printf("x = %.15e\n",ptr_moons->satellite[moon_array_idx(idx, j)].Sat_coord.x);
            printf("y = %.15e\n",ptr_moons->satellite[moon_array_idx(idx, j)].Sat_coord.y);
            printf("z = %.15e\n",ptr_moons->satellite[moon_array_idx(idx, j)].Sat_coord.z);
            printf("u = %.15e\n",ptr_moons->satellite[moon_array_idx(idx, j)].Sat_coord.u);
            printf("v = %.15e\n",ptr_moons->satellite[moon_array_idx(idx, j)].Sat_coord.v);
            printf("w = %.15e\n",ptr_moons->satellite[moon_array_idx(idx, j)].Sat_coord.w);
        }
    }

    /***************************************************************************/
    /******************** INTEGRATION TIME FRAME *******************************/
    /***************************************************************************/

    read_option(input_file->ID, INTEGER, "TIMESCALE", &time_frame->TIMESCALE);
    read_option(input_file->ID, REAL, "ti", &time_frame->ti);
    read_option(input_file->ID, REAL, "tf", &time_frame->tf);
    read_option(input_file->ID, LONG_INTEGER, "nt", &time_frame->nt);
    read_option(input_file->ID, INTEGER, "box_integration", &time_frame->box_integration);

    /*error checking*/
    MACRO_ERROR(((time_frame->TIMESCALE != TIMESCALE_MOON) &&
                 (time_frame->TIMESCALE != TIMESCALE_PLANET) &&
                 (time_frame->TIMESCALE != TIMESCALE_SECONDS)),
                ("\nUnknown choice of timescale (choose 0 or 1 only)\n"));

    if(prnt)
    {
        printf("TIMESCALE = %d\n",time_frame->TIMESCALE);
        printf("ti = %lf\n",time_frame->ti);
        printf("tf = %lf\n",time_frame->tf);
        printf("nt = %ld\n",time_frame->nt);
        printf("box_integration = %d\n",time_frame->box_integration);
    }

    /***************************************************************************/
    /******************** PARTICLE INITIAL COORDINATES *************************/
    /***************************************************************************/

    read_option(input_file->ID, INTEGER, "pos_mode", &grain->pos_mode);
    iv->pos_mode = grain->pos_mode;

    char iv_filename[256];
    read_option(input_file->ID, STRING, "iv_filename", iv_filename);
    iv->file = open_file(iv_filename, "r");


    if(prnt)
    {
        printf("start position mode = %d\n",grain->pos_mode);
        printf("initial values are read from: %s\n", iv_filename);
    }


    /***************************************************************************/
    /******************** OUTPUT ***********************************************/
    /***************************************************************************/

    read_option(input_file->ID, INTEGER, "filemod", &switches->output.filemod);

    /*error checking*/
    MACRO_ERROR(((switches->output.filemod != OUTPUT_ASCII) &&
                 (switches->output.filemod != OUTPUT_BINARY)),
                ("Unknown choice of output (choose 1 or 2 ONLY)\n"));

    if (switches->output.filemod == OUTPUT_ASCII)
    {
        strcpy(switches->output.filemodchar,"w");
    }
    if (switches->output.filemod == OUTPUT_BINARY)
    {
        strcpy(switches->output.filemodchar,"wb");
    }

    read_option(input_file->ID, INTEGER, "OUTPUTTIMESCALE", &switches->output.EARTHTIMESCALE);
    read_option(input_file->ID, INTEGER, "n_out", &switches->output.n_out);
    read_option(input_file->ID, STRING, "outputdir", output_file->dir);
    read_option(input_file->ID, INTEGER, "PSpaceOut", &switches->output.PSpaceOut);
    read_option(input_file->ID, INTEGER, "elemout", &switches->output.elemout);
    read_option(input_file->ID, INTEGER, "SatCoordOut", &switches->output.SatCoordOut);
    read_option(input_file->ID, INTEGER, "SatElemOut", &switches->output.SatElemOut);
    read_option(input_file->ID, INTEGER, "combined", &switches->output.combined);
    read_option(input_file->ID, INTEGER, "Sun", &switches->output.sun);
    read_option(input_file->ID, INTEGER, "TestsOut", &switches->output.TestsOut);
    read_option(input_file->ID, INTEGER, "JacobiConstant", &switches->output.JacobiConstant);
    read_option(input_file->ID, INTEGER, "TrafoCoord", &switches->output.TrafoCoord);
    read_option(input_file->ID, INTEGER, "RelCoord", &switches->output.RelCoord);
    read_option(input_file->ID, INTEGER, "EM_INTEGRATION", &switches->output.EM_INTEGRATION);
    read_option(input_file->ID, INTEGER, "EFFECTIVE_SPEED", &switches->output.EFFECTIVE_SPEED);
    read_option(input_file->ID, INTEGER, "STATISTICS", &switches->output.STATISTICS);
    read_option(input_file->ID, INTEGER, "PlasmaPar", &switches->output.PlasmaParOut);
    read_option(input_file->ID, INTEGER, "LorentzOut", &switches->output.LorentzOut);
    read_option(input_file->ID, INTEGER, "No_Output", &switches->output.No_Output);

    if(prnt)
    {
        printf("filemodus = %d\n",switches->output.filemod);
        printf("OutputTimeScale = %d\n",switches->output.EARTHTIMESCALE);
        printf("OutputFraction = %d\n",switches->output.n_out);
        printf("outputdir = %s\n",output_file->dir);
        printf("PSpaceOut = %d\n",switches->output.PSpaceOut);
        printf("elemout = %d\n",switches->output.elemout);
        printf("SatCoordOut = %d\n",switches->output.SatCoordOut);
        printf("SatElemOut = %d\n",switches->output.SatElemOut);
        printf("combined = %d\n",switches->output.combined);
        printf("TestsOut = %d\n",switches->output.TestsOut);
        printf("TrafoCoordOut = %d\n",switches->output.TrafoCoord);
        printf("RelCoordOut = %d\n",switches->output.RelCoord);
        printf("EM_INTEGRATION = %d\n",switches->output.EM_INTEGRATION);
        printf("EFFECTIVE_SPEED = %d\n",switches->output.EFFECTIVE_SPEED);
        printf("STATISTICS = %d\n",switches->output.STATISTICS);
        printf("PlasmaParOut = %d\n",switches->output.PlasmaParOut);
        printf("LorentzOut = %d\n",switches->output.LorentzOut);
        printf("No_Output_at_all = %d\n",switches->output.No_Output);
    }

    /*****************************************************************************/
    /******************** ADAPTIVE STEPSIZE CONTROL ******************************/
    /*****************************************************************************/


    double *atol = (double *) malloc(idx->neq*sizeof(double));
    if (atol == NULL)
    {
        printf("%s:%d: Memory allocation failed.n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    ipar->atol = atol;

    read_option(input_file->ID, REAL, "rtol", &ipar->rtol);
    read_option(input_file->ID, REAL, "wy", &ipar->wy);
    read_option(input_file->ID, REAL, "wdydt", &ipar->wdydt);
    read_option(input_file->ID, REAL, "ini_stepsize", &ipar->h1);
    read_option(input_file->ID, REAL, "min_stepsize", &ipar->hmin);
    read_option(input_file->ID, REAL, "max_stepsize", &ipar->hmax);

    read_option(input_file->ID, REAL, "grain.x.atol",  &ipar->atol[grain_pos_idx(idx, 0)]);
    read_option(input_file->ID, REAL, "grain.y.atol",  &ipar->atol[grain_pos_idx(idx, 1)]);
    read_option(input_file->ID, REAL, "grain.z.atol",  &ipar->atol[grain_pos_idx(idx, 2)]);
    read_option(input_file->ID, REAL, "grain.vx.atol", &ipar->atol[grain_vel_idx(idx, 0)]);
    read_option(input_file->ID, REAL, "grain.vy.atol", &ipar->atol[grain_vel_idx(idx, 1)]);
    read_option(input_file->ID, REAL, "grain.vz.atol", &ipar->atol[grain_vel_idx(idx, 2)]);

    if (switches->CURRENTS)
    {
        read_option(input_file->ID, REAL, "grain.phi.atol", &ipar->atol[grain_pot_idx(idx)]);
    }

    for (int j = moon_start_idx(); exist_moon(idx, j); j++)
    {
        read_option(input_file->ID, REAL, "moon.x.atol",  &ipar->atol[moon_pos_idx(idx, j, 0)]);
        read_option(input_file->ID, REAL, "moon.y.atol",  &ipar->atol[moon_pos_idx(idx, j, 1)]);
        read_option(input_file->ID, REAL, "moon.z.atol",  &ipar->atol[moon_pos_idx(idx, j, 2)]);
        read_option(input_file->ID, REAL, "moon.vx.atol", &ipar->atol[moon_vel_idx(idx, j, 0)]);
        read_option(input_file->ID, REAL, "moon.vy.atol", &ipar->atol[moon_vel_idx(idx, j, 1)]);
        read_option(input_file->ID, REAL, "moon.vz.atol", &ipar->atol[moon_vel_idx(idx, j, 2)]);
    }

    read_option(input_file->ID, REAL, "sun.x.atol",  &ipar->atol[sun_pos_idx(idx, 0)]);
    read_option(input_file->ID, REAL, "sun.y.atol",  &ipar->atol[sun_pos_idx(idx, 1)]);
    read_option(input_file->ID, REAL, "sun.z.atol",  &ipar->atol[sun_pos_idx(idx, 2)]);
    read_option(input_file->ID, REAL, "sun.vx.atol", &ipar->atol[sun_vel_idx(idx, 0)]);
    read_option(input_file->ID, REAL, "sun.vy.atol", &ipar->atol[sun_vel_idx(idx, 1)]);
    read_option(input_file->ID, REAL, "sun.vz.atol", &ipar->atol[sun_vel_idx(idx, 2)]);

    if(prnt)
    {
        printf("--- Adaptive step size control ---\n");
        printf("rtol = %le\n",ipar->rtol);
        printf("wy = %le\n",ipar->wy);
        printf("wdydt = %le\n",ipar->wdydt);
        printf("min_stepsize = %le\n",ipar->hmin);
        printf("max_stepsize = %le\n",ipar->hmax);

        // --- Grain --- Header
        printf("--- Grain: abs. tolerance  ---\n");
        printf("  x.atol = %le\n", ipar->atol[grain_pos_idx(idx, 0)]);
        printf("  y.atol = %le\n", ipar->atol[grain_pos_idx(idx, 1)]);
        printf("  z.atol = %le\n", ipar->atol[grain_pos_idx(idx, 2)]);
        printf(" vx.atol = %le\n", ipar->atol[grain_vel_idx(idx, 0)]);
        printf(" vy.atol = %le\n", ipar->atol[grain_vel_idx(idx, 1)]);
        printf(" vz.atol = %le\n", ipar->atol[grain_vel_idx(idx, 2)]);

        if (switches->CURRENTS)
        {
            printf("phi.atol = %le\n", ipar->atol[grain_pot_idx(idx)]);
        }

        for (int j = moon_start_idx(); exist_moon(idx, j); j++)
        {
            // --- Moon j --- Header
            printf("--- Moon %d: abs. tolerance  ---\n", j);
            printf("  x.atol = %le\n", ipar->atol[moon_pos_idx(idx, j, 0)]);
            printf("  y.atol = %le\n", ipar->atol[moon_pos_idx(idx, j, 1)]);
            printf("  z.atol = %le\n", ipar->atol[moon_pos_idx(idx, j, 2)]);
            printf(" vx.atol = %le\n", ipar->atol[moon_vel_idx(idx, j, 0)]);
            printf(" vy.atol = %le\n", ipar->atol[moon_vel_idx(idx, j, 1)]);
            printf(" vz.atol = %le\n", ipar->atol[moon_vel_idx(idx, j, 2)]);
        }

        // --- Sun --- Header
        printf("--- Sun: abs. tolerance  ---\n");
        printf("  x.atol = %le\n", ipar->atol[sun_pos_idx(idx, 0)]);
        printf("  y.atol = %le\n", ipar->atol[sun_pos_idx(idx, 1)]);
        printf("  z.atol = %le\n", ipar->atol[sun_pos_idx(idx, 2)]);
        printf(" vx.atol = %le\n", ipar->atol[sun_vel_idx(idx, 0)]);
        printf(" vy.atol = %le\n", ipar->atol[sun_vel_idx(idx, 1)]);
        printf(" vz.atol = %le\n", ipar->atol[sun_vel_idx(idx, 2)]);
    }
}


static FILE* create_file(char* dirname, char* filename, char* ext, int moon_num)
{
    char two_digits[3];             // string with the capacity for two digits and \0
    char buffer[256];

    if  (moon_num > 99)
    {
        printf("Too many satellites for two_digit string. Satellite number: %d\n", moon_num); 
        exit(EXIT_FAILURE);
    }

    strcpy(buffer, dirname);
    strcat(buffer, filename);
    strcat(buffer, ext);

    if (moon_num >= 0 && moon_num <= 99)
    {
        sprintf(two_digits,"%02i", moon_num);
        strcat(buffer, two_digits);
    }

    return open_file(buffer, "w");
}



/*---------------------------------------------------------------------------*/
/*-------------- SET UP OUTPUT DIRECTORY AND FILES --------------------------*/
/*---------------------------------------------------------------------------*/
/*****************************************************************************/
void prepare_files(type_file *input_file, type_file *output_file,
                   type_switches *switches, type_moons *moons, psvidx_t *idx)
{

    char jobname[256];
    char syscmd[256];
    char dirname[256];

    // subtract the .job extension of the JOBFILE for the directory name
    if (strstr(input_file->name,".job") != NULL)
    {
        int i;
        i = strlen(input_file->name);
        strncpy(jobname, input_file->name, (i-4));
        jobname[i-4] = '\0';    // It is really essential to give the ending sign!
    }

    // strcat(output_file->dir,"/");
    strcpy(dirname, output_file->dir);
    strcat(dirname,"/");
    strcat(dirname, jobname);
    strcat(dirname, "/");

    // open a directory for the run and put a copy of the jobfile there
    strcpy(syscmd,"mkdir -p -m=rwxrwxrwx ");
    strcat(syscmd, dirname);
    system(syscmd);
    strcpy(syscmd,"cp ");
    strcat(syscmd,input_file->name);
    strcat(syscmd," ");
    strcat(syscmd, dirname);

    // put copy of jobfile only if output is switched on*/
    if (switches->output.No_Output == 0)
        system(syscmd);


    // allocate memory for the satellite's output streams*/
    if (switches->output.SatCoordOut)
    {
        output_file->ODsatCoord = (FILE **)malloc((size_t)(moons->satN)*sizeof(FILE *));
        MACRO_ERROR((output_file->ODsatCoord == NULL), ("CANNOT ALLOCATE MEMORY FOR SATELLITE OUTPUT STREAM"));
    }

    if (switches->output.SatElemOut)
    {
        output_file->ODsatElem = (FILE **)malloc((size_t)(moons->satN)*sizeof(FILE *));
        MACRO_ERROR((output_file->ODsatElem == NULL), ("CANNOT ALLOCATE MEMORY FOR SATELLITE OUTPUT STREAM"));
    }


    // sinklog
    if (switches->output.No_Output == 0)
    {
        output_file->ODs = create_file(dirname, "sinklog", ".txt", -1);
    }

    // phase space coordinates
    if (switches->output.PSpaceOut)
    {
        output_file->ODpsc = create_file(dirname, jobname, ".dat", -1);
    }

    // particle's osculating elements
    if (switches->output.elemout)
    {
        output_file->ODoscel = create_file(dirname, jobname, ".elem", -1);
    }

    // files to store integrated grain potential and plasma-parameters
    if (switches->CURRENTS == 1)
    {
        // integrated grain potential
        if (switches->output.No_Output == 0)
        {
            output_file->ODphi = create_file(dirname, jobname, ".phi", -1);
        }

        // output of plasma parameters and B-field along grain trajectory
        if (switches->output.PlasmaParOut)
        {
            output_file->ODplasma_par = create_file(dirname, jobname, ".plasma_par", -1);
            output_file->ODbfield = create_file(dirname, jobname, ".B_field", -1);
        }
    }

    // file for output of Lorentz force and total relative force acting on the dust grain along its trajectory
    if (switches->output.LorentzOut)
    {
        output_file->ODLorentz = create_file(dirname, jobname, ".lorentz_force", -1);
    }

    // files for coordinates of the satellites
    if (switches->output.SatCoordOut)
    {
        for (int j = moon_start_idx(); exist_moon(idx, j); j++)
        {
            output_file->ODsatCoord[moon_array_idx(idx, j)] = create_file(dirname, jobname, ".SatCoord", j);
        }
    }

    // files for orbital elements of the satellites
    if (switches->output.SatElemOut)
    {
        for (int j = moon_start_idx(); exist_moon(idx, j); j++)
        {
            output_file->ODsatElem[moon_array_idx(idx, j)] = create_file(dirname, jobname, ".SatElem", j);
        }
    }

    // positions of the particle and the reference moon
    if (switches->output.combined)
    {
        output_file->ODcombined = create_file(dirname, jobname, ".combined", -1);
    }

    // position and velocity of the sun 
    if (switches->output.sun)
    {
        output_file->ODsun = create_file(dirname, jobname, ".sun", -1);
    }

    // file to store testing of the numerical integration results
    if (switches->output.TestsOut)
    {
        output_file->ODTestsOut = create_file(dirname, jobname, ".Tests", -1);
    }

    // file to store the Jacobi integral of the restricted three body problem,
    // this is for testing purposes. Make sure that only one moon is present
    if (switches->output.JacobiConstant)
    {
        output_file->ODJacobiConstant = create_file(dirname, jobname, ".jacobi", -1);
    }

    // file to store the coordinates of the dust grain and the main satellite
    // (Enceladus) in the corotating reference frame
    if(switches->output.TrafoCoord)
    {
        output_file->ODTrafoCoord = create_file(dirname, jobname, ".TrafoCoord", -1);
    }

    // file to store the relative coordinates between dust grain and Enceladus
    if(switches->output.RelCoord)
    {
        output_file->ODRelCoord = create_file(dirname, jobname, ".RelCoord", -1);
    }

    // open file to store start parameters
    output_file->ODstartparam = create_file(dirname, "startparam", ".txt", -1);

#ifdef dyne_debug_stepsize
    output_file->ipar = create_file(dirname, "ipar", ".dat", -1);
    output_file->stepsize_diag = create_file(dirname, "stepsize", ".dat", -1);
#endif

#ifdef dyne_debug_accelerations
    output_file->acc = create_file(dirname, "accelerations", ".dat", -1);
#endif

}

