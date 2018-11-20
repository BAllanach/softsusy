
/** 
   Project:     SOFTSUSY 
   File:        main.cpp
   Author:      Ben Allanach 
   Manual:      B.C. Allanach, hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
                B.C. Allanach, M.A. Bernhardt, arXiv:0903.1805, Comp. Phys. 
		Commun. 181 (2010) 232-245
   Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   Description: main calling program example: performs a scan of tan beta 
   (starting at CMSSM10.1.1) and prints out Higgs masses as a result in
   the format
       <tan beta>      <mh0>     <mA0>    <mH0>     <mH+/->
*/

#include "main.h"
#include "nmssmsoftsusy.h"
#include "decays.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

using namespace softsusy;

int main() {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  try {
    cout << "# SOFTSUSY" << PACKAGE_VERSION
	 << " test program, Ben Allanach 2002\n";
    cout << "# If you use SOFTSUSY, please refer to B.C. Allanach,\n";
    cout << "# Comput. Phys. Commun. 143 (2002) 305, hep-ph/0104145\n";
    
    /// Parameters used: CMSSM parameters
    double m12 = 500., a0 = 0., mGutGuess = 2.0e16, tanb = 20.0, m0 = 500.;
    int sgnMu = 1;       ///< sign of mu parameter 
    int numPoints = 2000;  ///< number of scan points
    
    QedQcd oneset;       ///< See "lowe.h" for default definitions parameters
    
    /// most important Standard Model inputs: you may change these and recompile
    double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
    oneset.setAlpha(ALPHAS, alphasMZ);
    oneset.setPoleMt(mtop);
    oneset.setMbMb(mbmb);
    
    oneset.toMz();       ///< Runs SM fermion masses to MZ
    oneset.runto(oneset.displayPoleMt());   ///< Runs SM fermion masses to mt
    
    /// Print out the SM data being used, as well as quark mixing assumption and
    /// the numerical accuracy of the solution
    cout << "# Low energy data in SOFTSUSY: mixing=0" << " TOLERANCE=" 
	 << TOLERANCE << endl;
    
    /// Print out header line
    cout << "#         M2       mchi1+      delta M     BR(jets)      BR(pi+)  BR(pi+ pi0) BR(nu e+)    BR(nu mu+)   BR(nu tau+)         Gamma\n";
    
    int i;
    /// Set limits of M2 scan
    double startM2 = 255, endM2 = 280.; ///< DEBUG
    
    /// Cycle through different points in the scan
    for (i = 0; i <= numPoints; i++) {
      double m2 = (endM2 - startM2) / static_cast<double>(numPoints) *
	static_cast<double>(i) + startM2;  /// set tan beta ready for the scan.

    /// Preparation for calculation: set up object and input parameters
    MssmSoftsusy r;

    DoubleVector pars(5);
    pars(1) = m0; pars(2) = m12; pars(3) = a0; pars(4) = m2;
    pars(5) = m12;
    bool uni = true;  ///< MGUT defined by g1(MGUT)=g2(MGUT)
    
    /// Calculate the spectrum
    const char * modelIdent = "nonUniversal";
    double qMax = 0.;
    r.lowOrg(nonUniGauginos, mGutGuess, pars, sgnMu, tanb, oneset, uni);

    
    NmssmSoftsusy nmssm;
    vector<Particle> decayTable;
    int err = calculateDecays(cout, &r, decayTable, nmssm, false);

    r.lesHouchesAccordOutput(cout, modelIdent, pars, sgnMu, tanb, qMax, 0,
    			     false);
    slhaDecays(cout, decayTable, false);
    exit(0);
    
    /// check the point in question is problem free: if so print the output
    if (!r.displayProblem().test()) {
      /// Sets format of output: 6 decimal places
      outputCharacteristics(6);

      
      if (decayTable.size() > 0 && err == 0) 
      cout << m2 << " " << r.displayPhys().mch(1) << " "
	   << r.displayPhys().mch(1) - r.displayPhys().mneut(1) << " "
	   << decayTable[chargino1].Array_Decays[45][5] +
	decayTable[chargino1].Array_Decays[26][5] << " " 
	   << decayTable[chargino1].Array_Decays[23][5] << " "
	   << decayTable[chargino1].Array_Decays[46][5] << " "
	   << decayTable[chargino1].Array_Decays[27][5] << " "
	   << decayTable[chargino1].Array_Decays[28][5] << " "
	   << decayTable[chargino1].Array_Decays[29][5] << " "
	   << decayTable[chargino1].total_width << endl;
    }
    else
      /// print out what the problem(s) is(are)
      cout << m2 << " " << r.displayProblem() << endl;
    }
    exit(0);
  }
  catch(const string & a) { cerr << a; return -1; }
  catch(const char * a) { cerr << a; return -1; }
  catch(...) { cerr << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}
