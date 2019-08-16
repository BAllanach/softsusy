/** \file main-nmssm.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach    
   - Manual:      B.C. Allanach, P. Athron, L. Tunstall, A. Voigt and 
   A. Williams, Comput. Phys. Comm. 185 (2014) 2322, arXiv:1311.7659;
   B.C. Allanach, hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305; 

   \brief a main C++ program to calculate Higgs masses as a function of tan
   beta in the NMSSM
*/

#include <iostream>
#include <cmath>
#include <cfloat>
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
  // double m12 = 300., a0 = -300., mGutGuess = 2.0e16, tanb = 10.0, m0 = 500.;
  int sgnMu = 1;      ///< sign of mu parameter
  int numPoints = 10; ///< number of scan points
  // double lambda = 0.1, kappa = 0.1, s = 0.0, xiF = 0.0, mupr = 0.0;

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 170.9, mbmb = 4.18;
  oneset.setAlphaMz(ALPHAS, alphasMZ);
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

  double LAMBDA3, LAMBDAX, mMess, cgrav;
  double mGutGuess;
  LAMBDA3 = 2.e5;
  LAMBDAX = 2.e5;
  mMess = 5.e5;
  cgrav = 1.0;
  mGutGuess = mMess;
    
  DoubleVector pars(4);
  pars(1) = LAMBDA3; pars(2) = LAMBDAX; pars(3) = mMess; pars(4) = cgrav;

  if (mMess < 1.0e3) {
  	  ostringstream ii; 
  	  ii << " mMess=" << mMess
  	     << " in SUGRA input (too low). The point will not yield a sensible answer\n";
  	  throw ii.str();
  	}
	
  	// r = &m;
  	if (LAMBDA3 > mMess) {
  	  ostringstream ii;
  	  ii << "Input LAMBDA3=" << LAMBDA3 << " should be less than mMess="
  	     << mMess << endl;
  	  throw ii.str();
  	}

  	if (LAMBDAX > mMess) {
  	  ostringstream ii;
  	    	  ii << "Input LAMBDAX=" << LAMBDAX << " should be less than mMess="
  	     << mMess << endl;
  	  throw ii.str();
  	}

  	if (cgrav > 1.0) {
  	  ostringstream ii;
  	  ii << "Input cgrav=" << cgrav << " a real number bigger than or "
  	     << " equal to 1 (you can use 1 as a default value).\n";
  	  throw ii.str();
  	}

	double lambda, kappa, mupr, xiF, s, xiS;	
  //point I:
   lambda = 0.8;
   kappa  = 0.1;
   mupr =  2.2e3;     //order 1e3;
   xiF  = -2.0e3;    //is it positive or negative ??
   xiS  = -8.0e7;
   s    =  730.;

  //point III:
  // lambda = 0.8;
  // kappa = 0.0;
  // mupr = 2.2e3;     //order 1e3;
  // xiF = 4.84e6;    //order 1e6;
  // s   = 730.;

   /*  DoubleVector nmpars(5);
  nmpars(1) = lambda; nmpars(2) = kappa; nmpars(3) = s;
  nmpars(4) = xiF; nmpars(5) = mupr;*/
  
  bool uni = false; // MGUT defined by g1(MGUT)=g2(MGUT)  // the final point of running is ***not*** MGUT
 
  double tanb = 4.;  
  //  for (int i = 0; i < numPoints; i++) {
  //     tanb = (endTanb - startTanb) / double(numPoints) * double(i) +
  //        startTanb; // set tan beta ready for the scan.


     NmssmSoftsusy n;
     //this was originally true - i.e. it was for Z_3 symmetric models     
     n.setZ3(false);   
     NMSSM_input nmssm_input; // NMSSM input parameters     
     try {     
     nmssm_input.set(NMSSM_input::lambda, lambda);
     nmssm_input.set(NMSSM_input::kappa, kappa);
     //     nmssm_input.set(NMSSM_input::lambdaS, lambda * s);
     nmssm_input.set(NMSSM_input::xiF, xiF);
     nmssm_input.set(NMSSM_input::xiS, xiS);
     nmssm_input.check_setup();
     DoubleVector nmpars(nmssm_input.get_nmpars());

     n.setGUTlambda(true);
     n.setGUTkappa(true);
     n.setGUTmuPrime(true);
     n.setGUTxiF(true);
     n.setGUTxiS(true);     
     //     n.setGUTsVev(true);
     
       n.NmssmSoftsusy::lowOrg(focusgmsb, mMess, pars, nmpars, 
  			       sgnMu, tanb, oneset, uni);
     } catch (const std::string& error) {
       cout << error;
       n.flagProblemThrown(true);
     } catch (const char* error) {
       cout << error;
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
     //  }

  return 0;
}
