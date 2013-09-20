
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

/// Returns an object with the relevant SM inputs
MssmSoftsusy doStuff(double alphasMZ, double alphaMZinv, double mtop, double mbmb, double mz) {

  MZ = mz;
  TOLERANCE = 1.e-5;
 /// Parameters used: CMSSM 40.2.5 parameters
  double m12 = 650., a0 = -500., mGutGuess = 2.0e16, tanb = 40.0, m0 = 750.;
  int sgnMu = 1;      ///< sign of mu parameter 

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setAlpha(ALPHA, 1. / alphaMZinv);
  oneset.setPoleMt(mtop);
  oneset.setMass(mBottom, mbmb);
  oneset.setMu(MZ);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Preparation for calculation: set up object and input parameters
  MssmSoftsusy r; 
  r.setMu(MZ);
  DoubleVector pars(3); 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
  
  /// Calculate the spectrum
  r.lowOrg(cmssmBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);

  cerr << r.displayProblem();

  return r;
}

int main() {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

 /// Sets format of output: 6 decimal places
  outputCharacteristics(6);

  ///< Central values and errors on SM inputs
  double alphasMZ = 0.1184, sigAs = 0.0007, sigAs2 = 0.1184 * 1.1 / 1.e3;  
  double alphaMZinv = 127.916, sigA = 0.015;   
  double mt = 173.18, sigMt = 0.6;
  double mb = 4.2, sigMb = 0.2;
  double mz = 91.1876, sigMz = 0.0021;

  /// To get high alphas
  MssmSoftsusy highAs(doStuff(alphasMZ + 2.0 * sigAs, alphaMZinv, 
			      mt + 2.0 * sigMt, mb + 2.0 * 
			      sigMb, mz));
  /// To get low alphas
  MssmSoftsusy lowAs(doStuff(alphasMZ - 2.0 * sigAs, alphaMZinv, 
			     mt - 2.0 * sigMt, mb - 2.0 * sigMb, mz));
  /// To get high alpha1
  MssmSoftsusy highA1(doStuff(alphasMZ - 2.0 * sigAs, alphaMZinv + 
			      2.0 * sigA, mt + 2.0 * sigMt,
			      mb + 2.0 * sigMb, mz + 2.0 * sigMz));
  /// To get low alpha1
  MssmSoftsusy lowA1(doStuff(alphasMZ + 2.0 * sigAs, alphaMZinv - 2.0 * sigA, 
			     mt - 2.0 * sigMt,
			     mb - 2.0 * sigMb, mz - 2.0 * sigMz)); 
  /// To get high alpha2
  MssmSoftsusy highA2(doStuff(alphasMZ - 2.0 * sigAs, alphaMZinv + 
			      2.0 * sigA, mt + 2.0 * sigMt,
			      mb + 2.0 * sigMb, mz - 2.0 * sigMz));
  /// To get low alpha2
  MssmSoftsusy lowA2(doStuff(alphasMZ + 2.0 * sigAs, alphaMZinv - 2.0 * sigA, 
			     mt - 2.0 * sigMt,
			     mb - 2.0 * sigMb, mz + 2.0 * sigMz));
  /// To get high alphas with LHeC errors
  MssmSoftsusy lhecHigh(doStuff(alphasMZ + 2.0 * sigAs2, alphaMZinv, 
			       mt + 2.0 * sigMt, mb + 2.0 *
			       sigMb, mz));
  /// To get low alphas
  MssmSoftsusy lhecLow(doStuff(alphasMZ - 2.0 * sigAs2, alphaMZinv, 
			       mt - 2.0 * sigMt, mb - 2.0 * 
			       sigMb, mz));

  const int numPoints = 50;
  /// Set limits of tan beta scan
  double startLogMx = log(MZ), endLogMx = log(3.e16);
  /// Cycle through different points in the scan
  for (int i = 0; i<=numPoints; i++) {

    double lnq = (endLogMx - startLogMx) / double(numPoints) * double(i) +
      startLogMx; // set tan beta ready for the scan.

    lowA1.runto(exp(lnq));
    highA1.runto(exp(lnq));
    lowA2.runto(exp(lnq));
    highA2.runto(exp(lnq));
    lowAs.runto(exp(lnq));
    highAs.runto(exp(lnq));
    lhecLow.runto(exp(lnq));
    lhecHigh.runto(exp(lnq));

    cout << exp(lnq) << " " 
	 << 4. * PI / sqr(lowA1.displayGaugeCoupling(1)) << " " 
	 << 4. * PI / sqr(highA1.displayGaugeCoupling(1)) << " " 
	 << 4. * PI / sqr(lowA2.displayGaugeCoupling(2)) << " "
	 << 4. * PI / sqr(highA2.displayGaugeCoupling(2)) << " "
	 << 4. * PI / sqr(lowAs.displayGaugeCoupling(3)) << " "
	 << 4. * PI / sqr(highAs.displayGaugeCoupling(3)) << " "
	 << 4. * PI / sqr(lhecHigh.displayGaugeCoupling(3)) << " "
	 << 4. * PI / sqr(lhecLow.displayGaugeCoupling(3)) << " # "
	 << endl;

  }

  cout << endl << endl;
  /// Coordinates of the box for gauge unification
  cout << "15 24.6\n15 25.8\n16.3 25.8\n16.3 24.6\n15 24.6\n";
 }
