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
const double fpi = 0.13041, mpiplus = 0.13957018, mpi0 = 0.1349766;
double charginoToNeutralino1pion(const MssmSoftsusy * m);

#endif
