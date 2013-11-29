
/** \file rgeroutines.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: edges etc, NOT for public release!

*/

#ifndef RGEROUTINES_H
#define RGEROUTINES_H

#include <iostream>
#include "utils.h"
#include "softsusy.h"
#include "softpars.h"
#include "susy.h"
#include "lowe.h"
#include "linalg.h"
#include "def.h"
#include <string.h>
#include "rpvsusypars.h"
#include "rpvsoft.h"

namespace softsusy {

const double factor = 1.;

/// mll'(max). Input: BR(chi20 -> l_R -> chi_1^0) and 
/// BR(chi20 -> l_L -> chi_1^0)
double mllMax(const MssmSoftsusy & r, double brchi20tolRchi10, 
	      double brchi20tolLchi10);
double mllMax(double msqL, double mchi20, double mer, double chi10);

double mlqspc(const MssmSoftsusy & r);
double mlqspc(double msqL, double mchi20, double mer, double mchi10);
double findChiSqSugra(double m0, double m12, double a0, double tanb, 
		      int sgnMu, int numRoutine);
double findChiSq(DoubleVector & pars, double tanb, int sgnMu, 
		 int numRoutine, void (*boundaryCondition)
		 (MssmSoftsusy &, const DoubleVector &), bool unified, 
		 MssmSoftsusy &);
double mhqMax(const MssmSoftsusy & r);


double mtautauMax(const MssmSoftsusy & r);
double mlqMax(const MssmSoftsusy & r);
double mlqMax(double msqL, double mchi20, double mer, double chi10);
double edgeRatio(const MssmSoftsusy & r);
double mhqMin(const MssmSoftsusy & r);
double mllqMin(const MssmSoftsusy & r);
double mllqMin(double msqL, double mchi20, double mer, double chi10);
double mttMax(const MssmSoftsusy & r);
double mll4max(const MssmSoftsusy & r);
double mlqMaxFar(const MssmSoftsusy & r);
double mlqMaxFar(double msqL, double mchi20, double mer, double chi10);
double mllqMax(const MssmSoftsusy & r);
double mllqMax(double msqL, double mchi20, double mer, double mchi10);
double mllbMin(const MssmSoftsusy & r);

// edge variables at SU3
double chiSqSU3(const MssmSoftsusy & r, double & chiMllMax, 
		double & chiMllqMax, double & chiMllqMin, 
		double & chiMlqMin, double & chiMlqMax, 
		double & chiMtautauMax, double msl, double msq);

} // namespace softsusy

#endif
