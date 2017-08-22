/** \file threeBodydecays.h
   - Project:     SOFTSUSY 
   - Author:      Tom Cridge, Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Code calculates three body decay modes 
*/

#ifndef THREEBODYDECAYS_H
#define THREEBODYDECAYS_H

#include "nmssmsoftsusy.h"
#include "decays.h"
#include "softsusy.h"
#include "physpars.h"
#include "lowe.h"
#include "softpars.h"
#include "softsusy.h"
#include "flavoursoft.h"
#include "susy.h"
#include "particle.h"
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include <complex>

using namespace std;

//double gluinoamplitudedecaydgausscharginoqqpbarfirsttwogen (double mgluino, double mchargino, double mquark, double mquarkp, double msqL, double msqpL, double g, double thetaL, double thetaR, double alphas, int charg, bool onetothree);
const double fpi = 0.13041 / sqrt(2.0), mpiplus = 0.13957018, mpi0 = 0.1349766;
double charginoToNeutralino1pion(const MssmSoftsusy * m);

double gluinoamplitudedecay (double m1, double m2, double m3,
			     double alphastrong); 
double gluinoamplitudedecaymix (double m1, double m2, double m3,
				double alphastrong, double squarkmix,
				double theta);
double squarkamplitudedecaygluino (double m1, double m2, double m3,
				   double alphastrong);
double squarkamplitudedecaygluinomix (double m1, double m2, double m3,
				      double alphastrong, double squarkmix,
				      double theta);
double squarkamplitudedecaycharginoW1 (double m1, double m2, double m3,
				       double g, double gamma);
double squarkamplitudedecaycharginoW2 (double m1, double m2, double m3,
				       double g, double gamma);
double squark1amplitudedecaycharginoW1mix
(double m1, double m2, double m3, double g, double gammaL, double gammaR,
 double theta, double beta, double mWboson, double runmt, double runmb,
 double torb); 
double squark1amplitudedecaycharginoW2mix
(double m1, double m2, double m3, double g, double gammaL, double gammaR,
 double theta, double beta, double mWboson, double runmt, double runmb,
 double torb);
double squark2amplitudedecaycharginoW1mix
(double m1, double m2, double m3, double g, double gammaL, double gammaR,
 double theta, double beta, double mWboson, double mup, double mdown,
 double torb);
double squark2amplitudedecaycharginoW2mix
(double m1, double m2, double m3, double g, double gammaL, double gammaR,
 double theta, double beta, double mWboson, double mup, double mdown,
 double torb);
double squarkLamplitudedecayneutralino
(double m1, double m2, double m3, double g, double gprime,
 DoubleMatrix & mixNeut, int neutralino, int uord ); 
double squarkRamplitudedecayneutralino
(double m1, double m2, double m3, double g, double gprime,
 DoubleMatrix & mixNeut, int neutralino, int uord );
double squark3amplitudedecayneutralino
(double m1, double m2, double m3, double mWboson, double theta, double beta,
 DoubleMatrix & mixNeut, double g, double gp, double runmt, int squark,
 int oneortwo,  int neutralino);
double squark3amplitudedecaysquark3Wboson
(double m1, double m2, double m3, double g, double thetat, double thetab,
 int m1torb, int m1oneortwo, int m3torb, int m3oneortwo);
double squark3amplitudedecaychargedHiggssquark3
(double m1, double m2, double m3, double g, double mWboson, double beta,
 double thetat, double thetab, double greekmu, double At, double Ab,
 double mt, double mb, int t1or2, int b1or2);
#endif
