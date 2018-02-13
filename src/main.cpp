
/** 
   Project:     SOFTSUSY 
   File:        main.cpp
   Author:      Ben Allanach 
   Manual:      B.C. Allanach, hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   Webpage:     http://softsusy.hepforge.org/
   Description: main calling program example: performs a scan of tan beta 
   (starting at CMSSM10.1.1) and prints out Higgs masses as a result in
   the format
       <tan beta>      <mh0>     <mA0>    <mH0>     <mH+/->
*/

#include "main.h"
#include "nmssmsoftsusy.h"
#include "decays.h"
using namespace softsusy;

int main() {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  try {
    /// Sets format of output: 6 decimal places
    outputCharacteristics(6);

    cout << "# SOFTSUSY" << PACKAGE_VERSION 
	 << " test program, Ben Allanach 2002\n";
    cout << "# If you use SOFTSUSY, please refer to B.C. Allanach,\n";
    cout << "# Comput. Phys. Commun. 143 (2002) 305, hep-ph/0104145\n";
    
    /// Parameters used: CMSSM parameters to get a high mh
    double m12 = 5000., a0 = -14.e3, mGutGuess = 2.0e16,
      tanb = 20.0, m0 = 5000.;
    int sgnMu = 1;       ///< sign of mu parameter 
    int numPoints = 500;  ///< number of scan points
    
    QedQcd oneset; ///< See "lowe.h" for default definitions parameters
    
    /// most important Standard Model inputs: you may change these and recompile
    double alphasMZ = 0.1181, mtop = 173.34, mbmb = 4.18;
    oneset.setAlpha(ALPHAS, alphasMZ);
    oneset.setPoleMt(mtop);
    oneset.setMbMb(mbmb);
    
    oneset.toMz();       ///< Runs SM fermion masses to MZ
    
    /// Print out the SM data being used, as well as quark mixing assumption and
    /// the numerical accuracy of the solution
    cout << "# Low energy data in SOFTSUSY: mixing=0" << " TOLERANCE=" 
	 << TOLERANCE << endl;
    
    /// Print out header line
    cout << "# mQEDxQCD   mh           mA           mH0          mH+-         mzrun        mwrun        mtrun        mA0run       mh0run\n";
    
    int i;
    /// Set limits of tan beta scan
    //    double startlnM = 1.53164e+02, endlnM = 1.547340e+02;
    double startlnM = MZ, endlnM = 6.0e3;    
    numHiggsMassLoops = 3;
    
    /// Cycle through different points in the scan
    for (i = 0; i <= numPoints; i++) {
      double mScale =
			  (endlnM - startlnM) / static_cast<double>(numPoints) *
			  static_cast<double>(i) + startlnM
			  ; /// set mScale ready for the scan.

      QedQcd twoset(oneset);
      twoset.runto(mScale);

      /// Preparation for calculation: set up object and input parameters
      USE_TWO_LOOP_GAUGE_YUKAWA = true;      
      MssmSoftsusy r;
      r.setLoops(3);     
      r.included_thresholds = 31; ///< all thresholds included
      
      DoubleVector pars(3);
      pars(1) = m0; pars(2) = m12; pars(3) = a0;
      bool uni = true;  ///< MGUT defined by g1(MGUT)=g2(MGUT)
      
      /// Calculate the spectrum
      r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, twoset, uni);

      r.runto(mScale);
      r.calcDrBarPars();

      double mtrun = r.displayDrBarPars().mt;


      
      /// check the point in question is problem free: if so print the output
      if (!r.displayProblem().test()) 
	cout << mScale << " " << r.displayPhys().mh0(1) << " " 
	     << r.displayPhys().mA0(1) << " " 
	     << r.displayPhys().mh0(2) << " " 
	     << r.displayPhys().mHpm << " "
	     << r.displayMzRun() << " "
	     << r.displayMwRun() << " "
	     << mtrun << " " 
	     << r.printLongDrbar()
	     << endl;
      else
	/// print out what the problem(s) is(are)
	cout << "#" << mScale << " " << r.displayProblem() << endl;
  }
  }
  catch(const string & a) { cerr << a; return -1; }
  catch(const char * a) { cerr << a; return -1; }
  catch(...) { cerr << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}
