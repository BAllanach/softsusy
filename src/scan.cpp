
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
  double alphasMZ = 0.1187, mtop = 173.2, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Print out the SM data being used, as well as quark mixing assumption and
  /// the numerical accuracy of the solution
  cout << "# Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
       << TOLERANCE << endl << oneset << endl;

  /// Print out header line
  cout << "#         m0          m12           a0         tanb           mh           Dmh           mA           dmA          mH           DmH           mH+          DmH+ ";
  cout << "         Dmg " << "          mg " << "        Dmsq " 
       << "         msq " << "        DmeL " << "         meL " 
       << "        DmeR " << "         meR " << "     Dmneut1 " 
       << "      mneut1 " << "     Dmneut2 " << "     mneut2  "
       << "     Dmneut3 " << "      mneut3 " << "    Dmneut4  " 
       << "      mneut4 " << "        dmtL " << "         mtL " 
       << "       DmtR  " << "         mtR " << "       DmbL  " 
       << "        mbL  " << "    DmtauL   " << "      mtauL  " 
       << "        dht  " << "          ht " << "       dhb   " 
       << "         hb  " << "      dhtau  " << "       htau  "
       << "    Dmchi+1  " << "      mchi+1 " << "     Dmchi+2 " 
       << "      mchi+2 " << endl;

  int i, j; 
  /// Set limits of tan beta scan
  double startM0 = 100.0, endM0 = 6000.0;
  double startM12 = 250., endM12 = 1000.;
  /// Cycle through different points in the scan
  for (i = 0; i<=numPoints; i++) {

    m12 = (endM12 - startM12) / double(numPoints) * double(i) +
      startM12; 

    for (j = 0; j<=numPoints; j++) {

      m0 = (endM0 - startM0) / double(numPoints) * double(j) +
	startM0;
      a0 = -2.0 * m0;
      
      if (numPoints == 0) { m0 = startM0; m12 = startM12; a0 = -2.0 * m0; }
      
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
      
      double msq2loop = (r.displayPhys().mu(2, 1) + r.displayPhys().mu(1, 1) +
			 r.displayPhys().md(2, 1) + r.displayPhys().md(1, 1)) * 
	0.25;
      
      MssmSoftsusy s; mGutGuess = 2.0e16;
      SOFTSUSY_THREELOOP_RGE = true;
      SOFTSUSY_TWOLOOP = true;
      SOFTSUSY_TWOLOOP_TQUARK_STRONG = true;
      SOFTSUSY_TWOLOOP_BQUARK_STRONG = true;
      SOFTSUSY_TWOLOOP_BQUARK_YUKAWA = true;
      SOFTSUSY_TWOLOOP_TAU_YUKAWA = true;
      SOFTSUSY_TWOLOOP_GS = true;
      s.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
      
      double msq3loop = (s.displayPhys().mu(2, 1) + s.displayPhys().mu(1, 1) +
			 s.displayPhys().md(2, 1) + s.displayPhys().md(1, 1)) * 
	0.25;
      
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
	   << " " << (1. - r.displayPhys().mHpm / s.displayPhys().mHpm)
	   << " " << r.displayPhys().mGluino
	   << " " << (1. - r.displayPhys().mGluino / s.displayPhys().mGluino)
	   << " " << msq2loop
	   << " " << (1. - msq2loop / msq3loop)
	   << " " << r.displayPhys().me(1, 1)
	   << " " << (1. - r.displayPhys().me(1, 1) / s.displayPhys().me(1, 1))
	   << " " << r.displayPhys().me(2, 1)
	   << " " << (1. - r.displayPhys().me(2, 1) / s.displayPhys().me(2, 1))
	   << " " << fabs(r.displayPhys().mneut(1))
	   << " " << (1. - fabs(r.displayPhys().mneut(1)) / fabs(s.displayPhys().mneut(1)))
	   << " " << fabs(r.displayPhys().mneut(2))
	   << " " << (1. - fabs(r.displayPhys().mneut(2)) / fabs(s.displayPhys().mneut(2)))
	   << " " << fabs(r.displayPhys().mneut(3))
	   << " " << (1. - fabs(r.displayPhys().mneut(3)) / fabs(s.displayPhys().mneut(3)))
	   << " " << fabs(r.displayPhys().mneut(4))
	   << " " << (1. - fabs(r.displayPhys().mneut(4)) / fabs(s.displayPhys().mneut(4)))
	   << " " << fabs(r.displayPhys().mu(1, 3))
	   << " " << (1. - fabs(r.displayPhys().mu(1, 3)) / fabs(s.displayPhys().mu(1, 3)))
	   << " " << r.displayPhys().mu(2, 3)
	   << " " << (1. - r.displayPhys().mu(2, 3) / s.displayPhys().mu(2, 3))
	   << " " << r.displayPhys().md(1, 3)
	   << " " << (1. - r.displayPhys().md(1, 3) / s.displayPhys().md(1, 3))
	   << " " << r.displayPhys().md(2, 3)
	   << " " << (1. - r.displayPhys().md(2, 3) / s.displayPhys().md(2, 3))
	   << " " << r.displayPhys().me(1, 3)
	   << " " << (1. - r.displayPhys().me(1, 3) / s.displayPhys().me(1, 3))
	   << " " << r.displayPhys().me(2, 3)
	   << " " << (1. - r.displayPhys().me(2, 3) / s.displayPhys().me(2, 3))
	   << " " << r.displayYukawaElement(YU, 3, 3)
	   << " " << (1. - r.displayYukawaElement(YU, 3, 3) / 
		      s.displayYukawaElement(YU, 3, 3))
	   << " " << r.displayYukawaElement(YD, 3, 3)
	   << " " << (1. - r.displayYukawaElement(YD, 3, 3) / 
		      s.displayYukawaElement(YD, 3, 3))
	   << " " << r.displayYukawaElement(YE, 3, 3)
	   << " " << (1. - r.displayYukawaElement(YE, 3, 3) / 
		      s.displayYukawaElement(YE, 3, 3))
	   << " " << r.displaySusyMu() 
	   << " " << (1. - r.displaySusyMu() / s.displaySusyMu())
	   << " " << r.displayPhys().mch(1)
	   << " " << (1. - r.displayPhys().mch(1) / s.displayPhys().mch(1))
	   << " " << r.displayPhys().mch(2)
	   << " " << (1. - r.displayPhys().mch(2) / s.displayPhys().mch(2));
      
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
