#include <iostream>
#include "nmssmsoftsusy.h"
#include "mycomplex.h"
#include "def.h"
#include "linalg.h"
#include "lowe.h"
#include "rge.h"
#include "softsusy.h"
#include "nmssmsoftpars.h"
#include "nmssmsusy.h"
#include "softpars.h"
#include "susy.h"
#include "utils.h"
#include "numerics.h"

int main() {

  /// Sets format of output: 6 decimal places
  outputCharacteristics(6);
  
  /// Parameters used: CMSSM parameters
  // double m12 = 500., a0 = -900., mGutGuess = 2.0e16, tanb = 5.0, m0 = 1000.;
    double m12 = 200., a0 = -500., mGutGuess = 2.0e16, tanb = 10.0, m0 = 200.;
  
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 1; ///< number of scan points
  //PA: new nmssm parameters to be fed into lowOrg
   double lambda = 0.1, kappa = 0.1, Al = a0, Ak = a0,  s = 1e3, xiF = 0.0, mupr = 0.0;
  // double lambda = 0.25, kappa = 0.1, Al = a0, Ak = -140,  s = 1e3, xiF = 0.0, mupr = 0.0;
  
  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters
 
  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMass(mBottom, mbmb);
  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Print out the SM data being used, as well as quark mixing assumption and
  /// the numerical accuracy of the solution
  cout << "# Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
       << TOLERANCE << endl << oneset << endl;

  /// Print out header line
  cout << "# tan beta   mh           mA           mH0          mH+-\n";

  // int i; 
  // /// Set limits of tan beta scan
  // double startTanb = 5.0, endTanb = 50.0;
  // /// Cycle through different points in the scan
  // for (i = 1; i<=numPoints; i++) {

  //   tanb = (endTanb - startTanb) / double(numPoints) * double(i) +
  //     startTanb; // set tan beta ready for the scan.
  
  NmssmSoftsusy n;
  //   cout << "n = "  << n << endl;
  // Fillfortest(n);
  //  cout << "n = "  << n << endl;
  // cout << "before calcDrBarPars " << endl;
  // n.calcDrBarPars();
  // cout << "after calcDrBarPars " << endl;
  //  n.physical(0);
  // cout << "after physical " << endl;
  // n.calcRunMtHiggs();

  DoubleVector pars(5); 
  pars(1) = m0; pars(2) = m12; pars(3) = a0, pars(4) = Al, pars(5) = Ak;
  DoubleVector nmpars(5);
  nmpars(1) = lambda; nmpars(2) = kappa; nmpars(3) = s; 
  nmpars(4) = xiF; nmpars(5) = mupr;
  bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
  Z3 = true;
  n.lowOrg(SemiMsugraBcs, mGutGuess, pars, nmpars, sgnMu, tanb, oneset, uni);

   /// check the point in question is problem free: if so print the output
    if (!n.displayProblem().test()) {

      cout << "  tanb =  "  << tanb << endl;
      cout << "mh1 = " << n.displayPhys().mh0(1) << endl;
      cout << "mh2 = " << n.displayPhys().mh0(2) << endl;
      cout << "mh3 = " << n.displayPhys().mh0(3) << endl;
       cout << "mA1 = " << n.displayPhys().mA0(1) << endl;
      cout << "mA2 = " << n.displayPhys().mA0(2) << endl;
      cout << "mHpm = " << n.displayPhys().mHpm << endl;
      cout << "msnu = " <<  n.displayPhys().msnu << endl;
      cout << "mch = " <<  n.displayPhys().mch << endl;
      cout << "mneut = " <<   n.displayPhys().mneut << endl;
      cout << "mGluino = " <<   n.displayPhys().mGluino << endl;
      cout << "mu = " <<   n.displayPhys().mu << endl;
      cout << "md = " <<   n.displayPhys().md << endl;
      cout << "me = " <<   n.displayPhys().me << endl;
      

      cout << "lambda = " << n.displayLambda() << endl;
      cout << "kappa = " << n.displayKappa() << endl;
      n.printall();

      cout << " At = "  << n.displaySoftA(UA, 3, 3) << endl;
      cout << " Ab = "  << n.displaySoftA(DA, 3, 3) << endl;
      cout << " Atau = "  << n.displaySoftA(EA, 3, 3) << endl;
      cout << " mQl = "  <<n.displaySoftMassSquared(mQl) << endl;
      cout << " mUr = "  <<n.displaySoftMassSquared(mUr) << endl;
      cout << " mDr = "  <<n.displaySoftMassSquared(mDr) << endl;
      cout << " mEr = "  <<n.displaySoftMassSquared(mEr) << endl;
      cout << " mLl = "  <<n.displaySoftMassSquared(mLl) << endl;
      
      cout << "mH1sq = " << n.displayMh1Squared() << endl;
      cout << "mH2sq = " << n.displayMh2Squared() << endl;
      cout << "mu = "  << n.displaySusyMu() << endl;
      cout << "m3sq = " << n.displayM3Squared() << endl;
      
       cout << "M1 = "  << n.displayGaugino(1) << endl;
       cout << "M2 = "  << n.displayGaugino(2) << endl;
       cout << "M3 = "  << n.displayGaugino(3) << endl;
    }



      // cout << tanb << " " << n.displayPhys().mh0(1) << " " 
      // 	   << n.displayPhys().mA0(1) << " " 
      // 	   << n.displayPhys().mh0(2) << " " 
      // 	   << n.displayPhys().mHpm << endl;
    else {
      cout << "Problem point! " << n.displayProblem() << endl;
       cout << "  tanb =  "  << tanb << endl;
       cout << "mh0 = " << n.displayPhys().mh0(1) << endl;
       cout << "mH0 = " << n.displayPhys().mh0(2) << endl;
       cout << "mH2 = " << n.displayPhys().mh0(3) << endl;
       cout << "mA1 = " << n.displayPhys().mA0(1) << endl;
       cout << "mA2 = " << n.displayPhys().mA0(2) << endl;
      cout << "mHpm = " << n.displayPhys().mHpm << endl;
      cout << "msnu = " <<  n.displayPhys().msnu << endl;
      cout << "mch = " <<  n.displayPhys().mch << endl;
      cout << "mneut = " <<   n.displayPhys().mneut << endl;
      cout << "mGluino = " <<   n.displayPhys().mGluino << endl;
      cout << "mu = " <<   n.displayPhys().mu << endl;
      cout << "md = " <<   n.displayPhys().md << endl;
      cout << "me = " <<   n.displayPhys().me << endl;
      n.printall();
      cout << " At = "  << n.displaySoftA(UA, 3, 3) << endl;
      cout << " Ab = "  << n.displaySoftA(DA, 3, 3) << endl;
      cout << " Atau = "  << n.displaySoftA(EA, 3, 3) << endl;
      cout << " mQl = "  <<n.displaySoftMassSquared(mQl) << endl;
      cout << " mUr = "  <<n.displaySoftMassSquared(mUr) << endl;
      cout << " mDr = "  <<n.displaySoftMassSquared(mDr) << endl;
      cout << " mEr = "  <<n.displaySoftMassSquared(mEr) << endl;
      cout << " mLl = "  <<n.displaySoftMassSquared(mLl) << endl;
      
      cout << "mH1sq = " << n.displayMh1Squared() << endl;
      cout << "mH2sq = " << n.displayMh2Squared() << endl;
         cout << "mu = "  << n.displaySusyMu() << endl;
      cout << "m3sq = " << n.displayM3Squared() << endl;
      
       cout << "M1 = "  << n.displayGaugino(1) << endl;
       cout << "M2 = "  << n.displayGaugino(2) << endl;
       cout << "M3 = "  << n.displayGaugino(3) << endl;


    }
 

  //  cout << "after lowOrg() " << endl;
  // n.printall();
  cout << "End of mainDev " << endl;
  return 1;
}
