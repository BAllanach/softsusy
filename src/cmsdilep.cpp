
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

double predMchi10(MssmSoftsusy & m, double mllmax) {
  double ml = m.displayPhys().me(2, 1);
  double m2 = m.displayPhys().mneut(2);

  double mChi10PredictionSquared = 
    (sqr(ml) * sqr(m2) - sqr(sqr(ml)) - sqr(mllmax) * sqr(ml)) / 
    (sqr(m2) - sqr(ml));
  return sqrt(mChi10PredictionSquared);
}

double predMchi20(MssmSoftsusy & m, double mllmax) {
  double ml = m.displayPhys().me(2, 1);
  double m1 = m.displayPhys().mneut(1);

  double mChi20PredictionSquared = 
    (sqr(ml) * sqr(m1) - sqr(sqr(ml)) - sqr(mllmax) * sqr(ml)) / 
    (sqr(m1) - sqr(ml));

  return sqrt(mChi20PredictionSquared);
}


double predMslep(MssmSoftsusy & m, double mllmax) {
  double m1 = m.displayPhys().mneut(1);
  double m2 = m.displayPhys().mneut(2);

  double twomSlepPredictionSquared = 
    sqr(m1) + sqr(m2) - sqr(mllmax) + 
    sqrt(sqr(sqr(mllmax) - sqr(m2) - sqr(m1)) - 4.0 * sqr(m1) * sqr(m2));
  return sqrt(twomSlepPredictionSquared * 0.5);
}

void cmsBcs(MssmSoftsusy & m, const DoubleVector & pars) {
  double m0 = 3500.;
  double m12 = 1600.;
  double a0 = 0.;

  double mh1sq = m.displayMh1Squared();
  double mh2sq = m.displayMh2Squared();

  m.standardSugra(m.displayMssmSusy(), m0, m12, a0);
  /// Makes sure that EWSB isn't screwed up
  m.setMh1Squared(mh1sq); m.setMh2Squared(mh2sq);
  m.setGauginoMass(1, pars.display(1));
  m.setGauginoMass(2, pars.display(2));

  double mslepSq = sqr(pars.display(3));
  double mslepLsq = sqr(pars.display(5));
  m.setSoftMassElement(mLl, 1, 1, mslepLsq);
  m.setSoftMassElement(mLl, 2, 2, mslepLsq);
  m.setSoftMassElement(mEr, 1, 1, mslepSq);
  m.setSoftMassElement(mEr, 2, 2, mslepSq);
  
  double msqSq = sqr(pars.display(4));
  m.setSoftMassElement(mQl, 1, 1, msqSq);
  m.setSoftMassElement(mQl, 2, 2, msqSq);
  m.setSoftMassElement(mUr, 1, 1, msqSq);
  m.setSoftMassElement(mUr, 2, 2, msqSq);
  m.setSoftMassElement(mDr, 1, 1, msqSq);
  m.setSoftMassElement(mDr, 2, 2, msqSq);

  return;
}


int main(int argc, char *argv[]) {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  try {
 /// Sets format of output: 6 decimal places
  outputCharacteristics(6);

  double m1, m2, msq, mslep, tanb, mslepL;
  bool m1SmallerThanM2 = true;  
  if (argc == 5) {
    mslep = atof(argv[1]);
    m2 = atof(argv[2]);
    msq = atof(argv[3]); 
    tanb = atof(argv[4]);
    if (m2 < mslep) m1SmallerThanM2 = false;
    mslepL = 2.0 * mslep;
  } else {
    throw("SOFTSUSY requires 4 arguments: mslepR m2 msq tanb\n");
    exit(0);
  }

  int sgnMu = 1;      ///< sign of mu parameter 

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Initial guess for mass
  if (m1SmallerThanM2) m1 = mslep * 0.5;
  else m1 = mslep * 2.0; 
  int numPoints = 40;
  MssmSoftsusy r; 

  DoubleVector pars(5); 
  bool uni = false, ewsbBCscale = true; double mGutGuess = 3.5e3;
  int i; for (i=0; i<=numPoints; i++) {
    MssmSoftsusy oldOne(r);
    r = MssmSoftsusy();
    r.useAlternativeEwsb();
    r.setMaCond(3.5e3); r.setMuCond(3.5e3);
    
    pars(1) = m1; pars(2) = m2; pars(3) = mslep; pars(4) = msq;
    pars(5) = mslepL;
    
    /// Calculate the spectrum
    r.lowOrg(cmsBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni, ewsbBCscale);
      
    double mchi2 = r.displayPhys().mneut(2), mchi=r.displayPhys().mneut(1), 
      msel = r.displayPhys().me(2, 1);
    double mllmax = sqrt((sqr(mchi2)-sqr(msel)) * 
			 (sqr(msel) - sqr(mchi)) / 
			 sqr(msel)); 
    const double tolerance = 1.0e-4;
    cout << "# diff=" << fabs((mllmax - 78.4) / 78.4);
    cout << " mll(max)=" << mllmax;
    if (m1SmallerThanM2) cout << " mchi10(pred)=" << predMchi10(r, 78.4);
    else cout << " mchi20(pred)=" << predMchi20(r, 78.4);
    if (fabs((mllmax - 78.4) / 78.4) < tolerance) {
      cout << endl; break;
    } else {
      if (m1SmallerThanM2) 
	m1 *= fabs(predMchi10(r, 78.4) / r.displayPhys().mneut(1));
      else m1 *= fabs(predMchi20(r, 78.4) / r.displayPhys().mneut(2));
    }
      cout << endl;
  }

  numPoints = 1;
  /// check the point in question is problem free: if so print the output
      double qMax = 0.;
      const char* modelIdent = "cmsDilep";
      if (!r.displayProblem().test()) 
      r.lesHouchesAccordOutput(cout, modelIdent, pars, sgnMu, tanb, qMax,  
      numPoints, ewsbBCscale);
      else
      /// print out what the problem(s) is(are)
      cout << tanb << " " << r.displayProblem() << endl;
  

  //  cout << r;
  /// Run up to GUT scale
  /*  double start = log(r.displayMu()), end = log(2.e16);
  numPoints = 40;
  for (i=0; i<=numPoints; i++) {
    double lnMu = (end - start) * double(i) / double (numPoints) + start;
    r.runto(exp(lnMu));
    cout << lnMu << " " << r.displayGaugino(3) << " " 
	 << signedSqrt(r.displaySoftMassSquared(mUr, 1, 1)) << endl;
	 }*/
  
  }
  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}
