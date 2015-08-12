
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

  TOLERANCE = 1.0e-4;

  try {
    /// Sets format of output: 6 decimal places
    outputCharacteristics(9);

    double m12 = 1000., m0 = 0., m0Overm12 = 0., a0 = 0., tanb = 10.;

    int sgnMu = 1;      ///< sign of mu parameter 
    
    QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

    /// most important Standard Model inputs: you may change these and recompile
    double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
    double lowRatio = 0.1, highRatio = 5.0;
    // double lowRatio = 0.7817552, highRatio = 0.7817581;
    oneset.setAlpha(ALPHAS, alphasMZ);
    oneset.setPoleMt(mtop);
    oneset.setMbMb(mbmb);    
    oneset.toMz();      ///< Runs SM fermion masses to MZ
    

    DoubleVector pars(3); 
    bool uni = true, ewsbBCscale = false; double mGutGuess = 1.e16;
    int numPoints = 100;
    int i; for (i=0; i<=numPoints; i++) {
      USE_TWO_LOOP_SPARTICLE_MASS = false;
      USE_TWO_LOOP_GAUGE_YUKAWA = true;
      USE_THREE_LOOP_RGE = true;

      MssmSoftsusy r; 
      
      m0Overm12 = lowRatio + 
	(highRatio - lowRatio) / double(numPoints) * double(i);
      m0 = m0Overm12 * m12;
      pars(1) = m0; pars(2) = m12; pars(3) = a0; 
      /// Calculate the 
      r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni, 
	       ewsbBCscale);
      //      const char* modelIdent = "sugra"; 
      //double qMax = 0.;
      //      r.lesHouchesAccordOutput(cout, modelIdent, pars, sgnMu, tanb, qMax, 
      //			       0, ewsbBCscale);
      
      MssmSoftsusy ho;
      USE_TWO_LOOP_SPARTICLE_MASS = true;
      ho.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni, 
		ewsbBCscale);
      //      ho.lesHouchesAccordOutput(cout, modelIdent, pars, sgnMu, tanb, qMax, 
      //			       0, ewsbBCscale);

      if (r.displayProblem().test()) cout << "# ";
      cout << m0Overm12 << " "                 // 1
	   << r.displayPhys().mGluino << " "   // 2
	   << ho.displayPhys().mGluino << " "  // 3
	   << r.displayPhys().mu(1, 1) << " "  // 4
	   << ho.displayPhys().mu(1, 1) << " " // 5
	   << r.displayPhys().mu(1, 3) << " "  // 6
	   << ho.displayPhys().mu(1, 3) << " " // 7
	   << r.displayPhys().mu(2, 3) << " "  // 8
	   << ho.displayPhys().mu(2, 3) << " " // 9
	   << r.displayPhys().md(1, 1) << " "  // 10
	   << ho.displayPhys().md(1, 1) << " " // 11
	   << r.displayPhys().md(1, 3) << " "  // 12
	   << ho.displayPhys().md(1, 3) << " " // 13
	   << r.displayPhys().mu(2, 1) << " "  // 14
	   << ho.displayPhys().mu(2, 1) << " " // 15
	   << r.displayPhys().md(2, 1) << " "  // 16
	   << ho.displayPhys().md(2, 1) << " " // 17
	   << r.displayPhys().md(2, 3) << " "  // 18
	   << ho.displayPhys().md(2, 3) << " ";// 19

      r.runto(r.displayMsusy());
      r.calcDrBarPars();
      cout << r.displayDrBarPars().mGluino  << " " // 20
	   << r.displayDrBarPars().mu(1, 1) << " " // 21
	   << r.displayDrBarPars().mu(1, 3) << " " // 22
	   << r.displayDrBarPars().mu(2, 3) << " " // 23
	   << r.displayDrBarPars().md(1, 1) << " " // 24
	   << r.displayDrBarPars().md(1, 3) << " " // 25
	   << r.displayDrBarPars().mu(2, 1) << " " // 26
	   << r.displayDrBarPars().md(2, 1) << " " // 27
	   << r.displayDrBarPars().md(2, 3) << " " // 28
	   << r.displayDrBarPars().mt       << " " // 29
	   << r.displayDrBarPars().mneut(1) << " " // 30
	   << r.displayDrBarPars().mneut(2) << " " // 31
	   << r.displayDrBarPars().mneut(3) << " " // 32
	   << r.displayDrBarPars().mneut(4) << " ";// 33

      if (r.displayProblem().test()) cout << " " << r.displayProblem();      
      cout << endl;
      //      cout << r.displayDrBarPars(); 
      //      exit(0);
    }
  }
  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}
