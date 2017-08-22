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
#endif
