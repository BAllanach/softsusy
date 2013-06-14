
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
  /// Sets format of output: 6 decimal places
  outputCharacteristics(6);

  /*  try{
  int check = 0;
  DoubleVector x(2); x(1) = 0.1; x(2) = 0.5; 
  bool err = newt(x, shoot);
  cout << "Finished. x=" << x << " err=" << err << endl;  
  }
  catch(const string & a) { cout << a; }
  exit(0);*/

  try {

  cerr << "SOFTSUSY" << SOFTSUSY_VERSION 
       << " test program, Ben Allanach 2002\n";
  cerr << "If you use SOFTSUSY, please refer to B.C. Allanach,\n";
  cerr << "Comput. Phys. Commun. 143 (2002) 305, hep-ph/0104145\n";

  /// Parameters used: CMSSM parameters
  double m12 = 660., a0 = 0., mGutGuess = 2.0e16, tanb = 40.0, m0 = 2800.;
  int sgnMu = -1;      ///< sign of mu parameter 
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
  TOLERANCE = 1.0e-4; MIXING=-1; 
  cout << "# Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
       << TOLERANCE << endl << oneset << endl;

  /// Print out header line
  cout << "# tan beta   mh           mA           mH0          mH+-\n";

  int i; 
  /// Set limits of tan beta scan
  double startM0 = 3000., endM0 = 6000.;
  /// Cycle through different points in the scan
  for (i = 0; i<=numPoints; i++) {

    m0 = (endM0 - startM0) / double(numPoints) * double(i) +
      startM0; // set tan beta ready for the scan.

    /// Preparation for calculation: set up object and input parameters
    MssmSoftsusy r, newtM; 
    DoubleVector pars(3); 
    pars(1) = m0; pars(2) = m12; pars(3) = a0;
    bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
    
    newtonMethod = false;
    /// Calculate the spectrum
    //  double mx = r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
    //  cout << m0 << " "     << r.displaySusyMu() << " ";

    vector<MssmSoftsusy> solutions;
    newtonMethod = true; 
    for (int j=1; j<=100; j++) { 
      numTry++;
      double mx2 = newtM.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, 
			      uni);

      if (!newtM.displayProblem().test()) {
	if (solutions.end() == solutions.begin()) {
	  solutions.push_back(newtM);
	}
	else {
	  bool differentSolution = true;
	  /// check through solutions to see if they are the same or not
	  for (vector<MssmSoftsusy>::iterator it = solutions.begin(); 
	       it!=solutions.end(); it++) 
	    if (sumTol(newtM, *it, 0) < 1.e-3) differentSolution = false;
	  if (differentSolution) {
	    solutions.push_back(newtM);
	  }
	}
      }
    }
    
    cout << m0 << " " << m12 << " " << a0 << " " << tanb << " ";
    cout << solutions.size() << " ";
    for (vector<MssmSoftsusy>::iterator it = solutions.begin(); 
	 it!=solutions.end(); it++) cout << it->displaySusyMu() << " ";
    cout << endl;
  }
  }
  catch(const string & a) { cout << a; }
  catch(const char * a) { cout << a; }
  catch(...) { cout << "Unknown type of exception caught.\n"; }

  exit(0);
}
