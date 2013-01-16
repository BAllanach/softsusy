
/**
   Project:     SOFTSUSY 3.0
   File:        rpvmain.cpp
   Authors:     B.C. Allanach, Markus Bernhardt 
   Manual:      B.C. Allanach and M.A. Bernhardt, arXiv:0903.1805 [hep-ph] and
                B.C. Allanach, M. Hanussek and S. Kom, arXiv:1109.3735
   Webpage:     http://projects.hepforge.org/softsusy/
   Description: main calling program example:
		- scanning SPS1a with one varying RPV coupling
*/

#include <rpvmain.h>

int main() {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  bool gaugeUnification = true, ewsbBCscale = false;

  /// Do we include 2-loop RGEs of *all* scalar masses and A-terms, or only the
  /// scalar mass Higgs parameters? (Other quantities all 2-loop anyway): the
  /// default in SOFTSUSY 3 is to include all 2-loop terms, except for RPV,
  /// which is already slow and calculated to less accuracy than the R-parity
  /// conserving version
  //  bool INCLUDE_2_LOOP_SCALAR_CORRECTIONS = false;

  /// Sets format of output: 6 decimal places
  outputCharacteristics(6);

  /// Header  
  cerr << "SOFTSUSY" << SOFTSUSY_VERSION 
       << " Ben Allanach, Markus Bernhardt 2009\n";
  cerr << "If you use SOFTSUSY, please refer to B.C. Allanach, ";
  cerr << "Comput. Phys. Commun. 143 (2002) 305, hep-ph/0104145;\n";
  cerr << "For RPV aspects, B.C. Allanach and M.A. Bernhardt, ";
  cerr << "Comp. Phys. Commun. 181 (2010) 232, ";
  cerr << "arXiv:0903.1805.\n";
  cerr << "For RPV neutrino masses, B.C. Allanach, M. Hanussek and S. Kom,\n";
  cerr << "arXiv:1109.3735 [hep-ph]\n";

  /// "try" catches errors in main program and prints them out
  try {

    /// Contains default quark and lepton masses and gauge coupling
    /// information 
    QedQcd oneset;      ///< See "lowe.h" for default parameter definitions 
    oneset.toMz();      ///< Runs SM fermion masses to MZ
    
    /// Print out the Standard Model data being used, as well as quark mixing
    /// assumption and the numerical accuracy of the solution
    cerr << "Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
	 << TOLERANCE << endl << oneset << endl;
    /// set parameters
    double tanb = 10.;
    int sgnMu = 1;
    double mgutGuess = 2.e16; 
    double a0 = 0.0, m12 = 250.0, m0 = 100.0; 
    
    /// number of points for scan
    const int numPoints = 20; 
    
    /// parameter region
    double Start = 0. , End = 0.144;
    
    DoubleVector pars(3);
    /// set basic entries in pars
    pars(1) = m0; pars(2) = m12; pars(3) = a0;
      
    cout << "l'_{331}(M_Z) m_snutau     # Problem flag" << endl;
    /// loop over parameter space region
    int ii; for (ii=0; ii<=numPoints; ii++){
      double lambda = Start + ((End - Start) / double(numPoints) * double(ii));
      
      lambda = 0.01;

      /// define rpvSoftsusy object
      RpvSoftsusy kw;
      
      /// Check for MSUSY BCs on SUSY RPV couplings
      susyRpvBCatMSUSY = true;
      TOLERANCE = 1.0e-5;

      /// set lambda coupling at MZ
      kw.setLam(1, 2, 1, lambda); 
      
      /// output parameters into double vector pars used by lowOrg
      kw.rpvDisplay(pars);
      
      /// generate spectrum in RpvSoftsusy object kw
      double mgut = kw.lowOrg(rpvSugraBcs, mgutGuess, pars, sgnMu,
		tanb, oneset, gaugeUnification, ewsbBCscale);
      
      /// outputs for this scan
      cout << lambda << "  " << kw.displayPhys().msnu.display(3) << " # " 
	   << kw.displayProblem() << endl;

      cout << kw;
      kw.runto(mgut);
      cout << kw;
      exit(0);
    }
  }
  catch(const string & a) {
    cout << a; exit(-1);
  }
  catch(const char *a) {
    printf("%s", a); exit(-1);
  }
  
}
