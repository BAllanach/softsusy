
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
    outputCharacteristics(16);
    ltini();
    
    //    cout << b0(2.406809e+03, 2.404110e+00, 2.404110e+00, 8.829185e+02) << endl; exit(0);
    //  cout << b1(6.554676e+2, 2.302428e+00, 2.534452e+03, 2.176440e+03) <<
    //  endl;
    //    m1 = 2.423902e+00, m2 = 3.858506e+03, q = 9.118760e+01;
    //cout << b1(p, m1, m2, q) << endl;
    //  cout << b1(2.849934e+00, 2.423902e+00, 3.858506e+03, 9.118760e+01) << endl;
    //    exit(0);

    //  setmudim(q * q);
  /*  cout << b0(p, m1, m2, q) << endl; 
  cout << B0(p*p, m1*m1, m2*m2) << endl;
  ltexi(); exit(0);
  */
  printDEBUG.push_back(false);
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  try {
    /// Sets format of output: 6 decimal places
    outputCharacteristics(6);
    
    /// most important Standard Model inputs: you may change these and recompile
    double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
    double m12 = 300., a0 = 0., tanb = 10.0, m0 = 3100.;
    double start = 100., end = 2000.;
    int numPoints = 5.;
    MssmSoftsusy r;
    QedQcd oneset;     oneset.toMz(); 
    int sgnMu = 1;
    for (int i=0; i<=numPoints; i++) {
      for (int j = 0; j<=numPoints; j++) {
	for (int k = 0; k<=numPoints; k++) {
	  for (int l = 0; l<=numPoints; l++) {
	    m0 = (end - start) * double(i) / double(numPoints) + start;
	    m12 = (end - start) * double(j) / double(numPoints) + start;	
	    a0 = -2. * m0 * double(k) / double(numPoints) + m0;
	    tanb = 58. * double(l) / double(numPoints) + 2.;

	    double mgutGuess = 2.e16;
	    DoubleVector pars(3);
	    pars(1) = m0; pars(2) = m12; pars(3) = a0;
	    void (*boundaryCondition)(MssmSoftsusy &, 
				      const DoubleVector &)=sugraBcs;	  
	    cout << "m0=" << m0 << " m12=" << m12 << " a0=" << a0 << " tb=" << tanb << " " << r.displayProblem() << endl;
	    double mgut = r.lowOrg(boundaryCondition, mgutGuess, pars, sgnMu,
				   tanb, oneset, true);
	  }
	}
      }
    }
  }
  catch(const string & a) { cout << a; }
  catch(const char * a) { cout << a; }
  catch(...) { cout << "Unknown type of exception caught.\n"; }
  
  exit(0);
}
