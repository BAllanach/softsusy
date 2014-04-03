
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

int main() {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  try {
 /// Sets format of output: 6 decimal places
  outputCharacteristics(6);

  cerr << "SOFTSUSY" << SOFTSUSY_VERSION 
       << " test program, Ben Allanach 2002\n";
  cerr << "If you use SOFTSUSY, please refer to B.C. Allanach,\n";
  cerr << "Comput. Phys. Commun. 143 (2002) 305, hep-ph/0104145\n";

  /// Parameters used: CMSSM parameters
  double m12 = 500., a0 = 0., mGutGuess = 2.0e16, tanb = 30.0, m0 = 125.;
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 0; ///< number of scan points

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Print out the SM data being used, as well as quark mixing assumption and
  /// the numerical accuracy of the solution
  cout << "# Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
       << TOLERANCE << endl << oneset << endl;

  /// Print out header line
  cout << " #        m0          m12           a0         tanb           mh           Dmh           mA           dmA          mH           DmH           mH+          DmH+\n";

  int i, j; 
  /// Set limits of tan beta scan
  double startM0 = 100.0, endM0 = 6000.0;
  double startM12 = 250., endM12 = 1000.;
  /// Cycle through different points in the scan
  for (i = 0; i<=numPoints; i++) {

    m0 = (endM0 - startM0) / double(numPoints) * double(i) +
      startM0;
    a0 = -2.0 * m0;

    for (j = 0; j<=numPoints; j++) {

    m12 = (endM12 - startM12) / double(numPoints) * double(j) +
      startM12; 

    /// Preparation for calculation: set up object and input parameters
    MssmSoftsusy r; 
    DoubleVector pars(3); 
    pars(1) = m0; pars(2) = m12; pars(3) = a0;
    bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
    
    /// Switch off 3-loop RGEs etc
    SOFTSUSY_THREELOOP_RGE = false;
    SOFTSUSY_TWOLOOP = false;
    SOFTSUSY_TWOLOOP_TQUARK_STRONG = false;
    SOFTSUSY_TWOLOOP_BQUARK_STRONG = false;
    SOFTSUSY_TWOLOOP_BQUARK_YUKAWA = false;
    SOFTSUSY_TWOLOOP_TAU_YUKAWA = false;
    SOFTSUSY_TWOLOOP_GS = false;
    /// Calculate the spectrum
    r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);

    MssmSoftsusy s; mGutGuess = 2.0e16;
    SOFTSUSY_THREELOOP_RGE = true;
    SOFTSUSY_TWOLOOP = true;
    SOFTSUSY_TWOLOOP_TQUARK_STRONG = true;
    SOFTSUSY_TWOLOOP_BQUARK_STRONG = true;
    SOFTSUSY_TWOLOOP_BQUARK_YUKAWA = true;
    SOFTSUSY_TWOLOOP_TAU_YUKAWA = true;
    SOFTSUSY_TWOLOOP_GS = true;
    s.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);

    /// check the point in question is problem free: if so print the output
    if (r.displayProblem().test() || s.displayProblem().test()) cout << "# ";
    cout << m0 << " " << m12 << " " << a0 << " " << tanb 
	 << " " << r.displayPhys().mh0(1)
	 << " " << (1. - r.displayPhys().mh0(1) / s.displayPhys().mh0(1))
	 << " " << r.displayPhys().mA0(1)
	 << " " << (1. - r.displayPhys().mA0(1) / s.displayPhys().mA0(1))
	 << " " << r.displayPhys().mh0(2)
	 << " " << (1. - r.displayPhys().mh0(2) / s.displayPhys().mh0(2))
	 << " " << r.displayPhys().mHpm 
	 << " " << (1. - r.displayPhys().mHpm / s.displayPhys().mHpm);
    if (r.displayProblem().test()) cout << " 2-loop problem: " 
					 << r.displayProblem();
    if (s.displayProblem().test()) cout << " 3-loop problem " 
					 << s.displayProblem(); 
    cout << endl;
  }
  }
  }
  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}
