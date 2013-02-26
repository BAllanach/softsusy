
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

/// Scans through mu for a given value of m0
void muScan(double m0, double mtop, double alphasMZ, double mbmb, double m12, 
	  double a0, double tanb, double start, double end) {
  
  MssmSoftsusy r;
  double mGutGuess = 2.0e16;
  /// Parameters used: CMSSM parameters
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 50; ///< number of scan points
  
  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters
  
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);
  oneset.toMz();      ///< Runs SM fermion masses to MZ

  DoubleVector pars(3); 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)

  TOLERANCE = 1.0E-4;

  cout << "# m0=" << m0 << " mt=" << mtop << " a_s(M_Z)=" << alphasMZ 
       << " mb(mb)=" << mbmb << "\n# m12=" << m12 << " a0=" << a0 
       << " tanb=" << tanb << endl;
  cout << "# mu(MSUSY)    MZ          (MZ:P/E^2)   MW(MW)      "
       << "mch(1)         " 
       << "mneut(1)      mneut(2)     PIZZT(MZ)    PIWWT(0)     "
       << " PIWWT(MW)    g1(MZ)       g2(MZ)       g3(MZ)       "
       << "t1/v1        t2/v2\n";

  /// Calculate the spectrum
  PRINTOUT = 0;
  for (int i = 0; i <= numPoints; i++) {
    double mu = (end - start) / double(numPoints) * double(i) + start;
    trialMuSq = sqr(mu);
    double mx = r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);

    r.calcDrBarPars();
    drBarPars s(r.displayDrBarPars());

    if (r.displayProblem().test()) cout << "# ";
    /// check the point in question is problem free: if so print the output
    cout << sqrt(trialMuSq) << " " 
	 << r.displayMu() << " " 
	 << r.displayPredMzSq() / sqr(MZ) << " " 
	 << s.mw << " " << " " << s.mch(1) << " " 
	 << " " << s.mneut(1) << " " << s.mneut(2) << " " 
	 << r.piZZT(r.displayMz(), r.displayMu(), true) << " "
	 << r.piWWT(0., r.displayMu(), true) << " " 
	 << r.piWWT(r.displayMw(), r.displayMu(), true) << " ";
    cout  << r.displayGaugeCoupling(1) << " " 
	  << r.displayGaugeCoupling(2) << " " 
	  << r.displayGaugeCoupling(3) << " " 
	  << r.displayTadpole1Ms() << " "
	  << r.displayTadpole2Ms() << " "
	  << mx << " ";
    r.runto(r.displayMsusy());
    r.calcDrBarPars();
    s = r.displayDrBarPars();      
    cout << s.mw << " " << " " << s.mch(1) << " " 
	 << " " << s.mneut(1) << " " << s.mneut(2) << " " 
	 << r.piZZT(r.displayMz(), r.displayMu(), true) << " "
	 << r.piWWT(0., r.displayMu(), true) << " " 
	 << r.piWWT(r.displayMw(), r.displayMu(), true) << " " << endl;
    if (r.displayProblem().test()) cout << r.displayProblem() << endl;
  }
  cout << endl;
}	       

/// Scans over m0, providing a scan of mu bet
void m0Scan(double mtop, double alphasMZ, double mbmb, double m12, double a0, 
	  double tanb, double muStart, double muEnd) {
  double mGutGuess = 2.0e16;
  /// Parameters used: CMSSM parameters
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 10; ///< number of scan points
  
  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters
  
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  TOLERANCE = 1.0E-4;

  /// Print out the SM data being used, as well as quark mixing assumption and
  /// the numerical accuracy of the solution

  for (int k=0; k <=numPoints; k++) {

    double m0 = 600. + k * 300.;
    muScan(m0, mtop, alphasMZ, mbmb, m12, a0, tanb, muStart, muEnd);
  }
}

int main() {
  printDEBUG.push_back(false);
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  try {
    /// Sets format of output: 6 decimal places
    outputCharacteristics(6);
    
    /// most important Standard Model inputs: you may change these and recompile
    double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
    double m12 = 300., a0 = 0., tanb = 10.0, m0 = 3050.;
    double start = 0.1, end = 200.;
    /*    muScan(m0, mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl << endl;

    //    start = 58., end = 59;
    //    muScan(m0, mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl << endl;
    exit(0);
    */
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;

    mtop = 172.5; cout << "# mt=172.5\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;

    mtop = 174.5; cout << "# mt=174.5\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;
    mtop = 173.5;

    alphasMZ = 0.1177; cout << "# as=0.1177\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;

    alphasMZ = 0.1197; cout << "# as=0.1197\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;
    alphasMZ = 0.1187;

    mbmb = 4.0; cout << "# mb=4.0\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;

    mbmb = 4.4; cout << "# mb=4.4\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;
    mbmb = 4.2;

    a0 = 500.; cout << "# a0=500\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;

    a0 = -500.; cout << "# a0=-500\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;
    a0 = 0.;

    tanb = 40.; cout << "# tb=40\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;
    tanb = 10.;

    m12 = 200.; cout << "# m12=200\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;

    m12 = 400.; cout << "# m12=400\n";
    m0Scan(mtop, alphasMZ, mbmb, m12, a0, tanb, start, end); cout << endl;
    m12 = 300.;
  }
  catch(const string & a) { cout << a; }
  catch(const char * a) { cout << a; }
  catch(...) { cout << "Unknown type of exception caught.\n"; }

  exit(0);
}
