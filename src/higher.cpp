
/** 
   Project:     SOFTSUSY 
   File:        main.cpp
   Author:      Ben Allanach 
   Manual:      B.C. Allanach,hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
                B.C. Allanach, M.A. Bernhardt, arXiv:0903.1805, Comp. Phys. 
		Commun. 181 (2010) 232-245
   Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   Description: main calling program example: performs a scan of tan beta 
   (starting at CMSSM10.1.1) and prints out Higgs masses as a result in
   the format
       <tan beta>      <mh0>     <mA0>    <mH0>     <mH+/->
*/

#include <iostream>
#include "mycomplex.h"
#include "def.h"
#include "linalg.h"
#include "lowe.h"
#include "rge.h"
#include "softsusy.h"
#include "softpars.h"
#include "susy.h"
#include "utils.h"
#include "numerics.h"

int main(int argc, char *argv[]) {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  try {
    /// Sets format of output: 6 decimal places
    outputCharacteristics(6);

    double m12 = 1000., m0 = 0., m0Overm12 = 0., a0 = 0., tanb = 10.;

    int sgnMu = 1;      ///< sign of mu parameter 
    
    QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

    /// most important Standard Model inputs: you may change these and recompile
    double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
    double lowRatio = 0.1, highRatio = 4.1;
    oneset.setAlpha(ALPHAS, alphasMZ);
    oneset.setPoleMt(mtop);
    oneset.setMbMb(mbmb);    
    oneset.toMz();      ///< Runs SM fermion masses to MZ
    
    MssmSoftsusy r, ho; 

    DoubleVector pars(3); 
    bool uni = true, ewsbBCscale = false; double mGutGuess = 1.e16;
    int numPoints = 50;
    int i; for (i=0; i<=numPoints; i++) {
      r = MssmSoftsusy(); ho = MssmSoftsusy();
      
      m0Overm12 = lowRatio + 
	(highRatio - lowRatio) / double(numPoints) * double(i);
      m0 = m0Overm12 * m12;
      pars(1) = m0; pars(2) = m12; pars(3) = a0; 
      
      /// Calculate the spectrum
      USE_TWO_LOOP_SPARTICLE_MASS = false;
      r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni, 
	       ewsbBCscale);

      USE_TWO_LOOP_SPARTICLE_MASS = true;
      ho.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni, 
		ewsbBCscale);
      
      if (r.displayProblem().test()) cout << "# ";
      cout << m0Overm12 << " " << r.displayPhys().mGluino << " " 
	   << ho.displayPhys().mGluino << " " 
	   << r.displayPhys().mu(1, 1) << " " 
	   << ho.displayPhys().mu(1, 1) << " " 
	   << r.displayPhys().mu(1, 3) << " " 
	   << r.displayPhys().mu(2, 3) << endl;
      if (r.displayProblem().test()) cout << " " << r.displayProblem();      
    }
  }
  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}
