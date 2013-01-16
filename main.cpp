
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
void line2() {
  /// Parameters used: CMSSM parameters
  double m12 = 500., a0 = 0., mGutGuess = 2.0e16, tanb = 10.0, m0 = 125.;
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 50; ///< number of scan points

  QedQcd oneset, twoset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(171.5);
  oneset.setMbMb(mbmb);
  oneset.toMz();      ///< Runs SM fermion masses to MZ
  twoset.setAlpha(ALPHAS, alphasMZ);
  twoset.setPoleMt(175.5);
  twoset.setMbMb(mbmb);
  twoset.toMz();      ///< Runs SM fermion masses to MZ

  /// Print out the SM data being used, as well as quark mixing assumption and
  /// the numerical accuracy of the solution
  /// Print out header line
  cout << "#     m0     m12   mchi10       mg     mstopR      mstopL      mchi1+   mchi2+        mA0\n";

  int i, j;  tanb = 10.;
  /// Set limits of tan beta scan
  double startM0 = 3000.0, endM0 = 3500.0;
  double startM12 = 300., endM12 = 450.;
  /// Cycle through different points in the scan
  //  for (i = 0; i<=numPoints; i++) {
    for (j = 0; j<=numPoints; j++) {

    m0 = (endM0 - startM0) / double(numPoints) * double(j) +
      startM0; // set tan beta ready for the scan.
    m12 = -m0 * 0.28 +1280.;

    /// Preparation for calculation: set up object and input parameters
    MssmSoftsusy r, s; 
    DoubleVector pars(3); 
    pars(1) = m0; pars(2) = m12; pars(3) = a0;
    bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
    
    cout << m0 << " " << m12 << " "; 

    /// Calculate the spectrum
    
    r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
    s.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, twoset, uni);

    cout << r.displayPhys().mneut(1) << " " << r.displayPhys().mGluino << " " 
	 << r.displayPhys().mu(1, 1) << " ";
	 
    cout << s.displayPhys().mneut(1) << " " << s.displayPhys().mGluino << " " 
	 << s.displayPhys().mu(1, 1) << " # ";
    cout << r.displayProblem() << endl;
  }
}

void point() {
  printMuState = true;
  TOLERANCE=1.0e-3;
  int sgnMu = 1;      ///< sign of mu parameter 

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  cout << oneset;

  /// Calculate the spectrum: convergent point
  double m12 = 337.5, a0 = 0., mGutGuess = 2.0e16, tanb = 10.0, m0 = 3000.;
  MssmSoftsusy r; 
  DoubleVector pars(3); 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
  
  cout << "# convergent point\n# iteration mu(MSUSY) m3^2(MSUSY)\n";
  r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
  cout << endl << endl;

  cout << "# Non-convergent point\n# iteration mu(MSUSY) m3^2(MSUSY)\n";
  m0 = 3400.; pars(1) = m0; 
  r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);  
  cout << endl << endl;
  
  cout << "# no EWSB point\n# iteration mu(MSUSY) m3^2(MSUSY)\n";
  m0 = 3500.; pars(1) = m0; 
  r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);  
}

void cmssmMap(double mtop, double tanb) {
  /// Parameters used: CMSSM parameters
  double m12 = 500., a0 = 0., mGutGuess = 2.0e16, m0 = 125.;
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 20; ///< number of scan points

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Print out header line
  cout << "# mt=" << mtop << " tanb=" << tanb << endl;
  cout << "#     m0     m12   mu    m3sq   fracDiff   \n";

  int i, j;  
  /// Set limits of tan beta scan
  double startM0 = 60., endM0 = 4000.0;
  double startM12 = 60., endM12 = 1200.;
  /// Cycle through different points in the scan
  for (i = 0; i<=numPoints; i++) {
    for (j = 0; j<=numPoints; j++) {
      
      m0 = (endM0 - startM0) / double(numPoints) * double(j) +
	startM0; // set tan beta ready for the scan.
      m12 = (endM12 - startM12) / double(numPoints) * double(i) +
	startM12;
      
      /// Preparation for calculation: set up object and input parameters
      MssmSoftsusy r; 
      DoubleVector pars(3); 
      pars(1) = m0; pars(2) = m12; pars(3) = a0;
      bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)     

      /// Calculate the spectrum
      TOLERANCE = 1.0e-3;
      r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
      
      cout << m0 << " " << m12 << " ";
      if (r.displayProblem().muSqWrongSign) cout << -r.displaySusyMu();
      else cout << r.displaySusyMu();
      cout << " " << r.displayM3Squared() << " " << r.displayFracDiff() 
	   << " " << fabs(r.displayPhys().mch(1)) << " " 
	   << minimum(r.displayPhys().me(1, 3), r.displayPhys().me(2, 3)) - 
	fabs(r.displayPhys().mneut(1)) << " " << fabs(r.displayPhys().mneut(1))
	   << " " << r.displayPhys().mGluino << " " 
	   << (r.displayPhys().mu(1, 1) + r.displayPhys().mu(2, 1) +
	       r.displayPhys().md(1, 1) + r.displayPhys().md(2, 1)) * 0.25
	   << " ";
      cout << " # " << r.displayProblem() << endl;
    }
  }
}

void zoom() {
  /// Parameters used: CMSSM parameters
  double m12 = 500., a0 = 0., mGutGuess = 2.0e16, tanb = 10.0, m0 = 125.;
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 20; ///< number of scan points

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Print out header line
  cout << "#     m0     m12   mchi10       mg     mstopR      mstopL      mchi1+   mchi2+        mA0\n";

  int i, j;  tanb = 10.;
  /// Set limits of tan beta scan
  double startM0 = 3000.0, endM0 = 3500.0;
  double startM12 = 300., endM12 = 450.;
  /// Cycle through different points in the scan
  for (i = 0; i<=numPoints; i++) {
    for (j = 0; j<=numPoints; j++) {
      
      m0 = (endM0 - startM0) / double(numPoints) * double(j) +
	startM0; // set tan beta ready for the scan.
      m12 = (endM12 - startM12) / double(numPoints) * double(i) +
	startM12;
      
      /// Preparation for calculation: set up object and input parameters
      MssmSoftsusy r; 
      DoubleVector pars(3); 
      pars(1) = m0; pars(2) = m12; pars(3) = a0;
      bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
      
      cout << m0 << " " << m12 << " "; 
      
      /// Calculate the spectrum
      TOLERANCE = 1.0e-3; //PRINTOUT = 1;
      r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
      
      /// check the point in question is problem free: if so print the output
      if (!r.displayProblem().noConvergence) {
	if (r.displayProblem().muSqWrongSign) cout << 0.5 << " "; 
	else cout << 0. << " ";
      } else {
	if (r.displayProblem().muSqWrongSign) cout << 0.8 << " "; 
	else cout << 1. << " ";      
      }
      sProblem s(r.displayProblem());
      
      double tb = 0.;
      cout << sqrt(r.displayPredMzSq()) << " ";
      cout << fabs(r.displayPhys().mch(1));
      cout << " # " << s << flush << endl;
    }
  }

  cout << "\n\n3000 450\n3500 300\n\n\n343 20\n343 100\n";
}


void test() {
  /// Parameters used: CMSSM parameters
  double m12 = 500., a0 = 0., mGutGuess = 2.0e16, tanb = 10.0, m0 = 125.;
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 20; ///< number of scan points

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Print out header line
  cout << "#     m0     m12   mchi10       mg     mstopR      mstopL      mchi1+   mchi2+        mA0\n";

  int numNoCons = 0, numBadMz = 0;
  int i, j;  tanb = 10.;
  /// Set limits of tan beta scan
  double startM0 = 3000.0, endM0 = 3500.0;
  double startM12 = 300., endM12 = 450.;
  /// Cycle through different points in the scan
  int count = 0;
  for (i = 0; i<=numPoints; i++) {
    for (j = 0; j<=numPoints; j++) {
      
      m0 = (endM0 - startM0) / double(numPoints) * double(j) +
	startM0; // set tan beta ready for the scan.
      m12 = (endM12 - startM12) / double(numPoints) * double(i) +
	startM12;
      
      /// Preparation for calculation: set up object and input parameters
      MssmSoftsusy r; 
      DoubleVector pars(3); 
      pars(1) = m0; pars(2) = m12; pars(3) = a0;
      bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
      
      cout << m0 << " " << m12 << " "; 
      
      /// Calculate the spectrum
      TOLERANCE = 1.0e-3; //PRINTOUT = 1;
      r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
      
      double tbIn = 0., tbOut; double sTin = r.displayPredMzSq();
      cout << sqrt(sTin) << " ";
      double sTout = sqr(MZ);
      if (r.displayProblem().noConvergence) numNoCons++;

      
      cout << r.displayFracDiff() << " " << " # " 
	   << r.displayProblem() << flush << endl;
    }
  }

  cout << "number of non convergences=" << numNoCons << endl;  
}


void line() {
  /// Parameters used: CMSSM parameters
  double m12 = 500., a0 = 0., mGutGuess = 2.0e16, tanb = 10.0, m0 = 125.;
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 50; ///< number of scan points

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Print out the SM data being used, as well as quark mixing assumption and
  /// the numerical accuracy of the solution
  /// Print out header line
  cout << "#     m0     m12   mchi10       mg     mstopR      mstopL      mchi1+   mchi2+        mA0\n";

  int i, j;  tanb = 10.;
  /// Set limits of tan beta scan
  double startM0 = 3000.0, endM0 = 3500.0;
  double startM12 = 300., endM12 = 450.;
  /// Cycle through different points in the scan
  //  for (i = 0; i<=numPoints; i++) {
    for (j = 0; j<=numPoints; j++) {

    m0 = (endM0 - startM0) / double(numPoints) * double(j) +
      startM0; // set tan beta ready for the scan.
    m12 = -m0 * 0.28 +1280.;

    /// Preparation for calculation: set up object and input parameters
    MssmSoftsusy r; 
    DoubleVector pars(3); 
    pars(1) = m0; pars(2) = m12; pars(3) = a0;
    bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
    
    cout << m0 << " " << m12 << " "; 

    /// Calculate the spectrum
    TOLERANCE = 1.0e-3;
    r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);

    cout << r.displayPhys().mneut(1) << " " << r.displayPhys().mGluino << " " 
	 << r.displayPhys().mu(1, 3) << " " << r.displayPhys().mu(2, 3) << " "
	 << r.displayPhys().mch(1) << " " << r.displayPhys().mch(2) << " " 
	 << r.displayPhys().mA0 << " ";
    cout <<  r.displayProblem().noConvergence << " ";
    
    sProblem s(r.displayProblem());

    TOLERANCE = 1.0e-2;
    r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
    cout << r.displayPhys().mneut(1) << " " << r.displayPhys().mGluino << " " 
	 << r.displayPhys().mu(1, 3) << " " << r.displayPhys().mu(2, 3) << " "
	 << r.displayPhys().mch(1) << " " << r.displayPhys().mch(2) << " " 
	 << r.displayPhys().mA0 << " ";
    cout <<  r.displayProblem().noConvergence 
	 << " # " << s << endl;
  }
}

int main() {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

 /// Sets format of output: 6 decimal places
  outputCharacteristics(6);

  zoom();
    //line();
  exit(0);

  cmssmMap(173.5, 10.);
  cout << endl << endl;

  cmssmMap(171.5, 10.);
  cout << endl << endl;

  cmssmMap(175.5, 10.);
  cout << endl << endl;
  
  cmssmMap(173.5, 50.);
  cout << endl << endl;

  cmssmMap(171.5, 50.);
  cout << endl << endl;

  cmssmMap(175.5, 50.);
  cout << endl << endl;

  exit(0);
}
