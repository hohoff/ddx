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






/*
 * Module for calculations to and from geometrical orbital elements
 * after:
 *
 *     Renner and Sicardy (2006)
 *     Use of geometric elements in numerical simulations
 *     Celestial Mechanics and Dynamical Astronomy 94:237–248
 *
 */

#include <stdio.h> 
#include <math.h>
#include "geomelem.h"



/*
 *   Module wide variables, describing the planet. 
 */

static double J2 = 16298.0e-6;         // Zonal harmonic cooefficients. Default values
static double J4 = -915.0e-6;          // are for Saturn taken from Campbell and Anderson
static double J6 = 103.0e-6;           // (1989).

static double Rp = 60330.0e3;          // Equatorial radius of the planet. Default values
                                       // are for Saturn taken from Campbell and Anderson
                                       // (1989). Unit: meter

static double GM = 3.7931272e16;       // Value from Campbell and Anderson (1989). Unit: m^3 s^-2


/*
 *   Module functions.
 */

void geomElemInit(double GMplanet, double rplanet, double jtwo, double jfour, double jsix)
{
    GM = GMplanet;
    Rp = rplanet;
    J2 = jtwo;
    J4 = jfour;
    J6 = jsix;
}


void geomElemFreq(double GMm, type_elements *elem, double *n, double *kappa, double* nu)
{
    double zeta = Rp/(elem->a);
    double zeta2 = pow(zeta, 2);
    double zeta4 = pow(zeta, 4);
    double zeta6 = pow(zeta, 6);

    /* Kepler frequency */
    double n0_sq = GM / pow(elem->a, 3);

    /* mean motion n, Eq. (xxx) Murray/Dermott (200x) */
    *n  = 1. + GMm/GM;
    *n += 1.5 * J2 * zeta2;                                //    3./2.   = 1.5
    *n += -1.875 * J4 * zeta4;                             //   15./8.  = 1.875
    *n += 2.1875 * J6 * zeta6;                             //   35./16.  = 2.1875
    *n *= n0_sq;
    *n  = sqrt(*n);

    /* epicyclic frequency kappa, Eq. (15) in Renner and Sicardy (2006) */
    *kappa  = 1. + GMm/GM;
    *kappa += -1.5 * J2 * zeta2;                            //    3./2.   = 1.5
    *kappa += 5.625 * J4 * zeta4;                           //   45./8.   = 5.625
    *kappa += -10.9375 * J6 * zeta6;                        //  175./16.  = 10.9375
    *kappa *= n0_sq;
    *kappa  = sqrt(*kappa);

    /* vertical frequency nu, Eq. (16) in Renner and Sicardy (2006) */
    *nu  = 1. + GMm/GM;
    *nu += 4.5 * J2 * zeta2;                               //    9./2.   = 4.5
    *nu += -9.375 * J4 * zeta4;                            //   75./8.   = 9.375
    *nu += 15.3125 * J6 * zeta6;                           //  245./16.  = 15.3125
    *nu *= n0_sq;
    *nu  = sqrt(*nu);
}


void DgeomElemFreq(double GMm, type_elements *elem, double *dn, double *dkappa, double* dnu)
{
    double zeta = Rp/(elem->a);
    double zeta2 = pow(zeta, 2);
    double zeta4 = pow(zeta, 4);
    double zeta6 = pow(zeta, 6);

    double n0_sq = GM / pow(elem->a, 5);

    /* change in mean motion dn */
    double nbot = 0.0;
    nbot  = 16.*(1.+ GMm/GM);
    nbot += 24.* J2* zeta2;
    nbot += -30.* J4* zeta4;
    nbot += 35.* J6* zeta6;
    nbot  = sqrt(nbot);

    double ntop = 0.0;
    ntop  = 16.*(1.+ GMm/GM);
    ntop += 40.* J2* zeta2;
    ntop += -70.* J4* zeta4;
    ntop += 105.* J6* zeta6;
    ntop *= sqrt(n0_sq);
    ntop *= -0.375;                                           // 3./8.    = 0.375

    *dn = ntop/nbot;

    /* change in epicyclic frequency dkappa */
    double kapbot = 0.0;
    kapbot  = 16.*(1.+ GMm/GM);
    kapbot += -24.* J2* zeta2;
    kapbot += 90.* J4* zeta4;
    kapbot += -175.* J6* zeta6;
    kapbot  = sqrt(kapbot);

    double kaptop = 0.0;
    kaptop  = 16.*(1.+ GMm/GM);
    kaptop += -40.* J2* zeta2;
    kaptop += 210.* J4* zeta4;
    kaptop += -525.* J6* zeta6;
    kaptop *= sqrt(n0_sq);
    kaptop *= -0.375;                                           // 3./8.    = 0.375

    *dkappa = kaptop/kapbot;

    /* change in vertical frequency dnu */
    double nubot = 0.0;
    nubot  = 16.*(1.+ GMm/GM);
    nubot += 72.* J2* zeta2;
    nubot += -150.* J4* zeta4;
    nubot += 245.* J6* zeta6;
    nubot  = sqrt(nubot);

    double nutop = 0.0;
    nutop  = 16.*(1.+ GMm/GM);
    nutop += 120.* J2* zeta2;
    nutop += -350.* J4* zeta4;
    nutop += 735.* J6* zeta6;
    nutop *= sqrt(n0_sq);
    nutop *= -0.375;                                           // 3./8.    = 0.375

    *dnu = nutop/nubot;
}


void SemiMajorAxis(double GMm, type_elements *elem)
{
  double a = elem->a;

  double n = 0.0;
  double kappa = 0.0;
  double nu = 0.0;
  geomElemFreq(GMm, elem, &n, &kappa, &nu);
  
  double n_sat = elem->om;
  double accuracy = 1e-15;

  if ( fabs(n-n_sat)/n_sat<= accuracy )
  {  elem->a = a;  }

  else
  {  double dn = 0.0;
     double dkappa = 0.0;
     double dnu = 0.0;
     DgeomElemFreq(GMm, elem, &dn, &dkappa, &dnu);
     elem->a= a- ((n-n_sat)/dn);

     geomElemFreq(GMm, elem, &n, &kappa, &nu);
     DgeomElemFreq(GMm, elem, &dn, &dkappa, &dnu);
     long int steps=1;
     long int max_steps = 1000;
   
     while (fabs(n-n_sat)/n_sat > accuracy)
     {   
       a=elem->a;
       elem->a = a- ((n-n_sat)/dn);
       geomElemFreq(GMm, elem, &n, &kappa, &nu);
       DgeomElemFreq(GMm, elem, &dn, &dkappa, &dnu);
       steps+=1;
       
       if (steps > max_steps)
       {  elem->a = a;
	  printf("Maximum Number of Steps reached!"); 
          break;  }
//       if (steps == 5)
//       {  break; }
     }
  }
}


void circOrbit_elem2psv(double t, type_elements *elem, double pos[], double vel[])
{
    double a = elem->a;
    double lop = elem->aop + elem->lan;
    double lambda = elem->om * (t - elem->Tper) + lop;    // mean longitude

    pos[0] = a * cos(lambda);
    pos[1] = a * sin(lambda);
    pos[2] = 0.0;

    vel[0] = - a * elem->om * sin(lambda);
    vel[1] = a * elem->om * cos(lambda);
    vel[2] = 0.0;
}


void geomElem2Coord(double t, double GMm, type_elements *elem, type_coordinates *coord)
{
    /* mean motion and the epicyclic and vertical frequencies */
    double n = 0.0;
    double kappa = 0.0;
    double nu = 0.0;
    geomElemFreq(GMm, elem, &n, &kappa, &nu);
    double n_sq = pow(n, 2);
    double kappa_sq = pow(kappa, 2);

    double zeta = Rp/(elem->a);
    double zeta2 = pow(zeta, 2);
    double zeta4 = pow(zeta, 4);
    double zeta6 = pow(zeta, 6);

    double n0_sq = GM / pow(elem->a, 3);

    /* Eq. (17) in Renner and Sicardy (2006) */
    double eta_sq = 0.0;
    eta_sq  = 1.;
    eta_sq += -2. * J2 * zeta2;
    eta_sq += 9.375 * J4 * zeta4;                           //  75./8. = 9.375
    eta_sq += -21.875 * J6 * zeta6;                         // 175./8. = 21.875
    eta_sq *= n0_sq;

    /* Eq. (18) in Renner and Sicardy (2006) */
    double chi_sq = 0.0;
    chi_sq  = 1.;
    chi_sq += 7.5 * J2 * zeta2;                             //  15./2.  = 7.5
    chi_sq += -21.875 * J4 * zeta4;                         // 175./8.  = 21.875
    chi_sq += 45.9375 * J6 * zeta6;                         // 735./16. = 45.9375
    chi_sq *= n0_sq;

    /* Eq. (19) in Renner and Sicardy (2006) */
    double alpha1 = 0.0;
    alpha1 = (2.*nu + kappa)/3.;

    /* Eq. (20) in Renner and Sicardy (2006) */
    double alpha2 = 0.0;
    alpha2 = 2.*nu - kappa;

    /* Eq. (21) in Renner and Sicardy (2006) */
    double alpha_sq = 0.0;
    alpha_sq = alpha1 * alpha2;

    double a = elem->a;
    double e = elem->e;
    double i = elem->i;
    double e_sq = pow(elem->e, 2);
    double i_sq = pow(elem->i, 2);
    double aop = elem->aop;
    double lan = elem->lan;

    /* mean longitude */
    double lop = elem->aop + elem->lan;
    double lambda = n * (t - elem->Tper) + lop;

    /* Eq. (2) in Renner and Sicardy (2006) */
    double r = 0.0;
    r  = 1.;
    r += -e * cos(lambda - lop);
    r += e_sq * (1.5*eta_sq/kappa_sq - 1.  - 0.5*eta_sq/kappa_sq*cos(2.*(lambda - lop)));
    r += i_sq * (0.75*chi_sq/kappa_sq - 1. + 0.25*chi_sq/alpha_sq*cos(2.*(lambda - lan)));
    r *= a;

    /* Eq. (3) in Renner and Sicardy (2006) */
    double L = 0.0;
    L  = lambda;
    L += 2.*e*n*sin(lambda - lop)/kappa;
    L += e_sq * (0.75 + 0.5*eta_sq/kappa_sq) * n/kappa * sin(2.*(lambda - lop));
    L += -i_sq * 0.25*chi_sq/alpha_sq * n/nu * sin(2.*(lambda - lan));

    /* Eq. (4) in Renner and Sicardy (2006) */
    double z = 0.0;
    z  = sin(lambda - lan);
    z += e * 0.5*chi_sq/(kappa*alpha1) * sin(2.*lambda - lop -lan);
    z += -1.5 * e * chi_sq/(kappa*alpha2) * sin(lop - lan);
    z *= a * i;

    /* Eq. (5) in Renner and Sicardy (2006) */
    double rdot = 0.0;
    rdot  = e*sin(lambda - lop);
    rdot += e_sq*eta_sq/kappa_sq * sin(2.*(lambda - lop));
    rdot += i_sq*0.5*chi_sq/alpha_sq * nu/kappa * sin(2.*(lambda - lan));
    rdot *= a * kappa;

    /* Eq. (6) in Renner and Sicardy (2006) */
    double Ldot = 0.0;
    Ldot  = 1.;
    Ldot += 2.*e*cos(lambda - lop);
    Ldot += e_sq * (3.5 - 3.*eta_sq/kappa_sq - 0.5*kappa_sq/n_sq + (1.5 + eta_sq/kappa_sq)*cos(2.*(lambda - lop)));
    Ldot += i_sq * (2. - 0.5*kappa_sq/n_sq - 1.5*chi_sq/kappa_sq - 0.5*chi_sq/alpha_sq * cos(2.*(lambda - lan)));
    Ldot *= n;

    /* Eq. (7) in Renner and Sicardy (2006) */
    double zdot = 0.0;
    zdot  = cos(lambda - lan);
    zdot += 0.5*e*chi_sq*(kappa+nu)/(kappa*alpha1*nu) * cos(2.*lambda - lop - lan);
    zdot += 1.5*e*chi_sq*(kappa-nu)/(kappa*alpha2*nu) * cos(lop - lan);
    zdot *= a * i * nu;

    /* Eqs. (8) -- (13) in Renner and Sicardy (2006) */
    coord->t = t;
    coord->x = r*cos(L);
    coord->y = r*sin(L);
    coord->z = z;
    coord->u = rdot*cos(L) - r*Ldot*sin(L);
    coord->v = rdot*sin(L) + r*Ldot*cos(L);
    coord->w = zdot;
}


/* Function to Calculate Help-Frequencies after Renner and Sicardy 2006 */
/* -------------------------------------------------------------------- */
void geomElemHelpFreq( double GMm, type_elements *elem, double nu, double kappa, double *eta_sq, double *chi_sq, double *alpha1, double *alpha2, double *alpha_sq)
{
  double zeta  = Rp/elem->a;
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double zeta6 = zeta2*zeta4;

  double n0_sq = GM/pow(elem->a, 3);

  /* Eq. (17) in Renner and Sicardy (2006) */
  *eta_sq  = 1.;
  *eta_sq += -2.*J2*zeta2;
  *eta_sq += 9.375*J4*zeta4;                           //  75./8. = 9.375
  *eta_sq += -21.875*J6*zeta6;                         // 175./8. = 21.875
  *eta_sq *= n0_sq;

  /* Eq. (18) in Renner and Sicardy (2006) */
  *chi_sq  = 1.;
  *chi_sq += 7.5*J2*zeta2;                             //  15./2.  = 7.5
  *chi_sq += -21.875*J4*zeta4;                         // 175./8.  = 21.875
  *chi_sq += 45.9375*J6*zeta6;                         // 735./16. = 45.9375
  *chi_sq *= n0_sq;

  /* Eq. (19) in Renner and Sicardy (2006) */
  *alpha1 = (2.*nu + kappa)/3.;


  /* Eq. (20) in Renner and Sicardy (2006) */
  *alpha2 = 2.*nu - kappa;


  /* Eq. (21) in Renner and Sicardy (2006) */
  *alpha_sq = *alpha1 * *alpha2;
}


/* Function to calculate the Longitude (exact value) */
/* ------------------------------------------------- */
void Longitude(type_coordinates *coord, double *L)
{  
  /* Equation 26 - 29 in Renner and Sicardy (2006) */

  if (coord->x < 0.)
  {  *L = atan( coord->y/coord->x ) + val_pi; }

  else if (coord->x == 0. && coord->y > 0.)
  {  *L = 0.5*val_pi; }

  else if (coord->x == 0. && coord->y < 0.)
  {  *L = 1.5*val_pi; }

  else if (atan(coord->y/coord->x) < 0.)
  {  *L = atan( coord->y/coord->x) + val_2pi; }

  else
  {  *L = atan( coord->y/coord->x ); }
}


/* Function to calculate the Argument of Pericentre */
/* ------------------------------------------------ */
void LongOfPeri(type_elements *elem, double r, double rc, double rdot, double rcdot, double kappa, double lambda, double *lop )
{
  double a = elem->a;

  double numerator   = 0.;
  double denumerator = 0.;
  double LamMinLop   = 0.;

  if (rc == 0. && rcdot == 0.)
  {  numerator   = rdot;
     denumerator = a*kappa*(1.-r/a); }

  else
  {  numerator   = rdot - rcdot;
     denumerator = a*kappa*( 1.-(r-rc)/a ); }
  
  if (denumerator == 0. && numerator >0.)
  {  LamMinLop = 0.5*val_pi; }

  else if (denumerator == 0. && numerator < 0.)
  {  LamMinLop = 1.5*val_pi; }

  else if (denumerator < 0.)
  {  LamMinLop = atan( numerator/denumerator ) + val_pi; }

  else if (denumerator == 0. && numerator == 0.)
  {  LamMinLop = lambda; }

  else if (atan(numerator/denumerator) < 0.)
  {  LamMinLop = atan( numerator/denumerator) + val_2pi; }

  else
  {  LamMinLop = atan( numerator/denumerator ); }
  
  *lop = lambda - LamMinLop;

  if (*lop < 0.)
  {  *lop += val_2pi; }
}


/* Function to calculate the Longitude of Ascending Node */
/* ----------------------------------------------------- */
void LongAscNode(double z, double zc, double zdot, double zcdot, double nu, double lambda, double *lan )
{  
  double numerator   = 0.0;
  double denumerator = 0.0;
  double LamMinLan   = 0.0;

  if (zc == 0. && zcdot == 0.)
  {  numerator   = nu*z;
     denumerator = zdot; }

  else
  {  numerator   = nu*(z-zc);
     denumerator = zdot-zcdot; }

  if (denumerator == 0. && numerator >0.)
  {  LamMinLan = 0.5*val_pi; }

  else if (denumerator == 0. && numerator < 0.)
  {  LamMinLan =  1.5*val_pi; }

  else if (denumerator < 0.)
  {  LamMinLan = atan( numerator/denumerator ) + val_pi; }

  else if (denumerator == 0. && numerator == 0.)
  {  LamMinLan = lambda; }

  else if (atan(numerator/denumerator) < 0.)
  {  LamMinLan = atan( numerator/denumerator) + val_2pi; }

  else
  {  LamMinLan = atan( numerator/denumerator ); }
  
  *lan = lambda - LamMinLan;

  if (*lan <0.)
  {  *lan += val_2pi; }
}


/* Transformation of Coordinates into geometrical Elements */
/* ------------------------------------------------------- */
void Coord2geomElem( double t, double GMm, type_coordinates *coord, type_elements *elem )
{

  double r    = 0.0;
  double L    = 0.0;
  double z    = 0.0;
  double rdot = 0.0;
  double Ldot = 0.0;
  double zdot = 0.0;

  double rc    = 0.0;
  double Lc    = 0.0;
  double zc    = 0.0;
  double rcdot = 0.0;
  double Lcdot = 0.0;
  double zcdot = 0.0;

  /* Equation 22 - 25 in Renner and Sicardy (2006) */
  r    = sqrt( coord->x*coord->x + coord->y*coord->y );
  Longitude(coord, &L);
  rdot = coord->u*cos(L) + coord->v*sin(L);
  Ldot = ( coord->v*cos(L) - coord->u*sin(L) )/r;
  z    = coord->z;
  zdot = coord->w;

  elem->a      = r;
  elem->e      = 0.;
  elem->i      = 0.;
  elem->lan    = 0.;
  elem->aop    = 0.;
  elem->lambda = 0.;
  elem->lop    = 0.;
  elem->pdot   = 0.;
  elem->odot   = 0.;
 
  double n     = 0.0;
  double kappa = 0.0;
  double nu    = 0.0;
  double pdot  = 0.0;
  double odot  = 0.0;
  geomElemFreq(GMm,elem,&n,&kappa,&nu);
  pdot = n - kappa;
  odot = n - nu;

  double eta_sq   = 0.0;
  double chi_sq   = 0.0;
  double alpha1   = 0.0;
  double alpha2   = 0.0;
  double alpha_sq = 0.0;
  geomElemHelpFreq(GMm, elem, nu, kappa, &eta_sq, &chi_sq, &alpha1, &alpha2, &alpha_sq);
 
  /* Equation 42 - 47 in Renner and Sikardy (2006) */
  double a = 0.0;
  a  = r;
  a /= 1. - 0.5*(Ldot-n)/n;
  
  double e = 0.0;
  e  = 0.25*pow((Ldot-n)/n,2);
  e += pow(rdot/(elem->a*kappa),2);
  e  = sqrt(e);
  
  double i = 0.0;
  i  = pow(z/elem->a,2);
  i += pow(zdot/(elem->a*nu),2);
  i  = sqrt(i);
  
  double lambda = 0.0;
  lambda = L-2.*n/kappa*rdot/(elem->a*kappa);

  double lop = 0.0;
  double lan = 0.0;
  double aop = 0.0;
  LongOfPeri(elem, r, rc, rdot, rcdot, kappa, lambda, &lop);
  LongAscNode(z, zc, zdot, zcdot, nu, lambda, &lan);
  aop = lop-lan;

  elem->om     = n;
  elem->pdot   = pdot;
  elem->odot   = odot;
  elem->e      =  e;
  elem->i      =  i;
  elem->lambda =  lambda;
  elem->lan    = lan;
  elem->aop    = aop;
  elem->lop    = lop;

  double eps = 1e-15;
  long int steps = 1;
  long int maxSteps = 1000;

  while (fabs(elem->a-a)/elem->a > eps)
  {
    elem->a   = a;
    elem->e   = e;
    elem->i   = i;
    elem->aop = aop;
    elem->lan = lan;
    elem->lop = lop;
    elem->lambda = lambda;
    elem->pdot = pdot;
    elem->odot = odot;
    elem->om   = n;

    geomElemFreq(GMm, elem, &n, &kappa, &nu);
    geomElemHelpFreq(GMm, elem, nu, kappa, &eta_sq, &chi_sq, &alpha1, &alpha2, &alpha_sq);
    pdot = n-kappa;
    odot = n-nu;

    double e_sq     = pow(e,2);
    double i_sq     = pow(i,2);
    double n_sq     = pow(n,2);
    double kappa_sq = pow(kappa,2);
  
    /* Equation 36 - 41 in Renner and Sicardy (2006) */
    rc  = a*e_sq*( 1.5*eta_sq/kappa_sq - 1. - 0.5*eta_sq/kappa_sq*cos( 2.*(lambda-lop)) );
    rc += a*i_sq*( 0.75*chi_sq/kappa_sq - 1. + 0.25*chi_sq/alpha_sq*cos( 2.*(lambda-lan)) );
  
    Lc  = e_sq*( 0.75 + 0.5*eta_sq/kappa_sq )*n/kappa*sin( 2.*(lambda-lop) );
    Lc += -i_sq*0.25*(chi_sq/alpha_sq)*(n/nu)*sin( 2.*(lambda-lan) );
  
    zc = a*i*e*( 0.5*chi_sq/(kappa*alpha1)*sin( 2.*lambda-lop-lan ) - 1.5*chi_sq/(alpha2*kappa)*sin( lop-lan ) );
  
    rcdot  = a*e_sq*eta_sq/kappa*sin( 2.*(lambda-lop) );
    rcdot += -a*i_sq*0.5*chi_sq/alpha_sq*nu*sin( 2.*(lambda-lan) );
  
    Lcdot  = e_sq*n*( 3.5 - 3.*eta_sq/kappa_sq - 0.5*kappa_sq/n_sq + (1.5+eta_sq/kappa_sq)*cos( 2.*(lambda-lop) ) );
    Lcdot += i_sq*n*( 2. - 0.5*kappa_sq/n_sq - 1.5*chi_sq/kappa_sq - 0.5*chi_sq/alpha_sq*cos( 2.*(lambda-lan) ) );
  
    zcdot = a*i*e*( 0.5*chi_sq*(kappa+nu)/(kappa*alpha1)*cos( 2.*lambda - lop - lan ) + 1.5*chi_sq*(kappa-nu)/(kappa*alpha2)*cos( lop-lan ) );

    a  = r - rc;
    a /= 1. - 0.5*(Ldot-Lcdot-n)/n;
  
    e  = 0.25*pow((Ldot-Lcdot-n)/n,2);
    e += pow((rdot-rcdot)/(elem->a*kappa),2);
    e  = sqrt(e);
  
    i  = pow((z-zc)/elem->a,2);
    i += pow((zdot-zcdot)/(elem->a*nu),2);
    i  = sqrt(i);
  
    lambda = L-2.*n/kappa*(rdot-rcdot)/(elem->a*kappa);
  
    LongOfPeri(elem, r, rc, rdot, rcdot, kappa, lambda, &lop);
    LongAscNode(z, zc, zdot, zcdot, nu, lambda, &lan);
    aop = lop - lan;

    steps += 1;

    if (steps > maxSteps)
    {  printf("Maximum Steps reached!"); 
      break; }
  }
}    
