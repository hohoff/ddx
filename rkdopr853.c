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

/*
 * Runge Kutta Dormand Prince eigth order
 * from "Hairer, Norset, Wanner: Solving Ordinary Differential Equations 1"
 *
 */


static double c2  = 0.526001519587677318785587544488E-01;
static double c3  = 0.789002279381515978178381316732E-01;
static double c4  = 0.118350341907227396726757197510E+00;
static double c5  = 0.281649658092772603273242802490E+00;
static double c6  = 0.333333333333333333333333333333E+00;
static double c7  = 0.25E+00;
static double c8  = 0.307692307692307692307692307692E+00;
static double c9  = 0.651282051282051282051282051282E+00;
static double c10 = 0.6E+00;
static double c11 = 0.857142857142857142857142857142E+00;

static double b1 =   5.42937341165687622380535766363E-2;
static double b6 =   4.45031289275240888144113950566E0;
static double b7 =   1.89151789931450038304281599044E0;
static double b8 =  -5.8012039600105847814672114227E0;
static double b9 =   3.1116436695781989440891606237E-1;
static double b10 = -1.52160949662516078556178806805E-1;
static double b11 =  2.01365400804030348374776537501E-1;
static double b12 =  4.47106157277725905176885569043E-2;

static double bhh1 = 0.244094488188976377952755905512E+00;
static double bhh2 = 0.733846688281611857341361741547E+00;
static double bhh3 = 0.220588235294117647058823529412E-01;

static double er1  =  0.1312004499419488073250102996E-01;
static double er6  = -0.1225156446376204440720569753E+01;
static double er7  = -0.4957589496572501915214079952E+00;
static double er8  =  0.1664377182454986536961530415E+01;
static double er9  = -0.3503288487499736816886487290E+00;
static double er10 =  0.3341791187130174790297318841E+00;
static double er11 =  0.8192320648511571246570742613E-01;
static double er12 = -0.2235530786388629525884427845E-01;

static double a21 =    5.26001519587677318785587544488E-2;
static double a31 =    1.97250569845378994544595329183E-2;
static double a32 =    5.91751709536136983633785987549E-2;
static double a41 =    2.95875854768068491816892993775E-2;
static double a43 =    8.87627564304205475450678981324E-2;
static double a51 =    2.41365134159266685502369798665E-1;
static double a53 =   -8.84549479328286085344864962717E-1;
static double a54 =    9.24834003261792003115737966543E-1;
static double a61 =    3.7037037037037037037037037037E-2;
static double a64 =    1.70828608729473871279604482173E-1;
static double a65 =    1.25467687566822425016691814123E-1;
static double a71 =    3.7109375E-2;
static double a74 =    1.70252211019544039314978060272E-1;
static double a75 =    6.02165389804559606850219397283E-2;
static double a76 =   -1.7578125E-2;

static double a81 =    3.70920001185047927108779319836E-2;
static double a84 =    1.70383925712239993810214054705E-1;
static double a85 =    1.07262030446373284651809199168E-1;
static double a86 =   -1.53194377486244017527936158236E-2;
static double a87 =    8.27378916381402288758473766002E-3;
static double a91 =    6.24110958716075717114429577812E-1;
static double a94 =   -3.36089262944694129406857109825E0;
static double a95 =   -8.68219346841726006818189891453E-1;
static double a96 =    2.75920996994467083049415600797E1;
static double a97 =    2.01540675504778934086186788979E1;
static double a98 =   -4.34898841810699588477366255144E1;
static double a101 =   4.77662536438264365890433908527E-1;
static double a104 =  -2.48811461997166764192642586468E0;
static double a105 =  -5.90290826836842996371446475743E-1;
static double a106 =   2.12300514481811942347288949897E1;
static double a107 =   1.52792336328824235832596922938E1;
static double a108 =  -3.32882109689848629194453265587E1;
static double a109 =  -2.03312017085086261358222928593E-2;

static double a111 =  -9.3714243008598732571704021658E-1;
static double a114 =   5.18637242884406370830023853209E0;
static double a115 =   1.09143734899672957818500254654E0;
static double a116 =  -8.14978701074692612513997267357E0;
static double a117 =  -1.85200656599969598641566180701E1;
static double a118 =   2.27394870993505042818970056734E1;
static double a119 =   2.49360555267965238987089396762E0;
static double a1110 = -3.0467644718982195003823669022E0;
static double a121 =   2.27331014751653820792359768449E0;
static double a124 =  -1.05344954667372501984066689879E1;
static double a125 =  -2.00087205822486249909675718444E0;
static double a126 =  -1.79589318631187989172765950534E1;
static double a127 =   2.79488845294199600508499808837E1;
static double a128 =  -2.85899827713502369474065508674E0;
static double a129 =  -8.87285693353062954433549289258E0;
static double a1210 =  1.23605671757943030647266201528E1;
static double a1211 =  6.43392746015763530355970484046E-1;



/*
 * Input:
 *     t       time
 *     h       stepsize
 *     neq     number of equations to integrate
 *     
 *     y[nvars]     values at time t
 *     dydt[neq]    derivatives at time t
 * 
 *     fpar    all values which are needed for the calculation of the
 *             derivatives at time t'
 * 
 * Output:
 *     yout[nvars]  values at time t+h
 *     yerr[neq]    estimate of the local truncation error
 * 
 *
 * The twelve stages:
 * 
 *    0   |
 *    c2  |  a21
 *    c3  |  a31  a32h
 *    .   |  .
 *    .   |  .
 *    .   |  .
 *    c12 |  a121 a122 a123 ...  a1211
 *    ----+---------------------------------
 *        |  b1   b2   b3        b11    b12
 *    ----+---------------------------------
 *        |  b'1  b'2  b'3  ...  b'11   b'12
 * 
 * 
 */




void rkdopr853(double y[], double dydt[], long int neq, double t, double h,
              double yout[], double yerr3[], double yerr5[], type_force_par *fpar)
{
    psvidx_t *idx = fpar->idx;

    double k2[neq];
    double k3[neq];
    double k4[neq];
    double k5[neq];
    double k6[neq];
    double k7[neq];
    double k8[neq];
    double k9[neq];
    double k10[neq];
    double k11[neq];
    double k12[neq];
    double k13[neq];
    double ytemp[idx->nvars];

    for (int i = 0; i < idx->nvars; i++)
        ytemp[i] = 0.0;
    ytemp[time_idx(idx)] = -1.e30;

    for (int i = 0; i < neq; i++)
        ytemp[i] = y[i] + a21*h*dydt[i];

    (fpar->forces)(t+c2*h, ytemp, k2, fpar);
    for (int i = 0; i < neq; i++)
        ytemp[i] = y[i] + h*(a31*dydt[i]+a32*k2[i]);

    (fpar->forces)(t+c3*h, ytemp, k3, fpar);
    for (int i = 0; i < neq ; i++)
        ytemp[i] = y[i] + h*(a41*dydt[i]+a43*k3[i]);

    (fpar->forces)(t+c4*h, ytemp, k4, fpar);
    for (int i = 0; i < neq; i++)
        ytemp[i] = y[i] + h*(a51*dydt[i]+a53*k3[i]+a54*k4[i]);

    (fpar->forces)(t+c5*h, ytemp, k5, fpar);
    for (int i = 0; i < neq; i++)
        ytemp[i] = y[i] + h*(a61*dydt[i]+a64*k4[i]+a65*k5[i]);

    (fpar->forces)(t+c6*h, ytemp, k6, fpar);
    for (int i = 0; i < neq; i++)
        ytemp[i] = y[i] + h*(a71*dydt[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);

    (fpar->forces)(t+c7*h, ytemp, k7, fpar);
    for (int i = 0; i < neq; i++)
        ytemp[i] = y[i] + h*(a81*dydt[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]);

    (fpar->forces)(t+c8*h, ytemp, k8, fpar);
    for (int i = 0; i < neq; i++)
        ytemp[i] = y[i] + h*(a91*dydt[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i]);

    (fpar->forces)(t+c9*h, ytemp, k9, fpar);
    for (int i = 0; i < neq; i++)
        ytemp[i] = y[i] + h*(a101*dydt[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i]+a108*k8[i]+a109*k9[i]);

    (fpar->forces)(t+c10*h, ytemp, k10, fpar);
    for (int i = 0; i < neq; i++)
        ytemp[i] = y[i] + h*(a111*dydt[i]+a114*k4[i]+a115*k5[i]+a116*k6[i]+a117*k7[i]+a118*k8[i]+a119*k9[i]+a1110*k10[i]);

    (fpar->forces)(t+c11*h, ytemp, k11, fpar);
    for (int i = 0; i < neq; i++)
        ytemp[i] = y[i] + h*(a121*dydt[i]+a124*k4[i]+a125*k5[i]+a126*k6[i]+a127*k7[i]+a128*k8[i]+a129*k9[i]+a1210*k10[i]+a1211*k11[i]);

    (fpar->forces)(t+h, ytemp, k12, fpar);
    for (int i = 0 ; i < neq; i++)
    {
        k13[i] = b1*dydt[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k11[i]+b12*k12[i];
        yout[i] = y[i] + h*k13[i];
    }

    for (int i = 0; i < neq; i++)
    {
        yerr5[i] = er1*dydt[i]+er6*k6[i]+er7*k7[i]+er8*k8[i]+er9*k9[i]+er10*k10[i]+er11*k11[i]+er12*k12[i];
        yerr3[i] = k13[i] - bhh1*dydt[i] - bhh2*k9[i] - bhh3*k12[i];
    }
}

