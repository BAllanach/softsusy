/** \file main-nmssm.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach    
   - Manual:      B.C. Allanach, P. Athron, L. Tunstall, A. Voigt and 
   A. Williams, Comput. Phys. Comm. 185 (2014) 2322, arXiv:1311.7659;
   B.C. Allanach, hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305; 

   \brief A main C++ program to calculate Higgs masses as a function of tan
   beta in the NMSSM. Other than tan beta, this mimics the Z3 violating MSSM
   example SLHA input point. It provides an example of a Bayesian naturalness
   calculation (`BN' in the output).
*/

#include <iostream>
#include "def.h"
#include "linalg.h"
#include "lowe.h"
#include "nmssmsoftsusy.h"

using namespace softsusy;

int main() {
  /// Sets format of output: 6 decimal places
  outputCharacteristics(6);
  softsusy::PRINTOUT = 0;

  /// Parameters used: CMSSM parameters
  double m12 = 350., a0 = -300., mGutGuess = 2.0e16, tanb = 10.0, m0 = 400.;
  int sgnMu = 1;      ///< sign of mu parameter
  int numPoints = 10; ///< number of scan points
  double lambda = 0.1,
    kappa = 0.1,
    s = 2000.0,
    xiF = 100.0,
    mupr = 0.0;

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1181, mtop = 170.9, mbmb = 4.20;
  double aInv = 127.918;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setAlphaMz(ALPHAS, alphasMZ);
  oneset.setAlphaMz(ALPHA, 1.0 / aInv);
  oneset.setPoleMt(mtop);
  oneset.setMass(mBottom, mbmb);
  oneset.toMz();      ///< Runs SM fermion masses to MZ
  oneset.runto(oneset.displayPoleMt());      ///< Runs SM fermion masses to mt

  /// Print out the SM data being used, as well as quark mixing assumption and
  /// the numerical accuracy of the solution
  cout << "# Data in SOFTSUSY: mixing=0" << " TOLERANCE="
       << TOLERANCE << endl;

  /// Print out header line
  cout << "# tan beta   mh(1)        mh(2)        mA(1)        mA(2)"
       << "        mH+-         BN   \n";

  /// Set limits of tan beta scan
  double startTanb = 5.0, endTanb = 55.0;

  DoubleVector pars(5);
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  pars(4) = a0, pars(5) = a0;
  DoubleVector nmpars(5);
  nmpars(1) = lambda; nmpars(2) = kappa; nmpars(3) = s;
  nmpars(4) = xiF; nmpars(5) = mupr;
  bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)

  for (int i = 0; i < numPoints; i++) {
     tanb = (endTanb - startTanb) / double(numPoints) * double(i) +
        startTanb; // set tan beta ready for the scan.

     NmssmSoftsusy n;
     n.setZ3(false);
     NMSSM_input nmssm_input; // NMSSM input parameters     
     nmssm_input.set(NMSSM_input::lambda, lambda);
     nmssm_input.set(NMSSM_input::kappa, kappa);
     nmssm_input.set(NMSSM_input::lambdaS, lambda * s);
     nmssm_input.set(NMSSM_input::xiF, xiF);
     nmssm_input.check_setup();
     DoubleVector nmpars(nmssm_input.get_nmpars());

     n.setGUTlambda(true);
     n.setGUTkappa(true);
     n.setGUTmuPrime(true);
     n.setGUTxiF(true);
     n.setGUTsVev(false);
     n.setMixing(2);
     
     try {
       n.NmssmSoftsusy::lowOrg(NmssmMsugraBcs, mGutGuess, pars,
			       nmpars, sgnMu, tanb, oneset, uni);
     } catch (const std::string& error) {
       n.flagProblemThrown(true);
     } catch (const char* error) {
       n.flagProblemThrown(true);
     }

     /// check the point in question is problem free: if so print the output
     if (!n.displayProblem().test()) {
        cout << tanb << ' '
             << n.displayPhys().mh0(1) << ' '
             << n.displayPhys().mh0(2) << ' '
             << n.displayPhys().mA0(1) << ' '
             << n.displayPhys().mA0(2) << ' '
             << n.displayPhys().mHpm   << ' '
	     << n.calcBayesianNaturalness() << '\n';
     } else {
        cout << tanb << ' ' << n.displayProblem() << '\n';
     }
  }

  return 0;
}
