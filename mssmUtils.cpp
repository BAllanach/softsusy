
#include "mssmUtils.h"
#include "softsusy.h"

#include <iostream>

/// Returns true if a point passes the Higgs constraint from LEP2, false
/// otherwise.  Error is the amount of uncertainty on SOFTSUSY's mh prediction
bool testLEPHiggs(const MssmSoftsusy & r, double error) {
  double Mh = r.displayPhys().mh0(1);
  Mh = Mh + error;
  double sinba2 = sqr(sin(atan(r.displayTanb()) - r.displayPhys().thetaH));

  ///  cout << "sinba2=" << sinba2 << endl;

  if (Mh < 90.0) return false;
  else if (90.0 <= Mh &&  Mh < 99.0) {
      if (sinba2 < -6.1979 + 0.12313 * Mh - 0.00058411 * sqr(Mh)) return true;
      else return false;
    }
  else if (99.0 <= Mh &&  Mh < 104.0) {
      if (sinba2 < 35.73 - 0.69747 * Mh + 0.0034266 * sqr(Mh)) return true;
      else return false;
    }
  else if (104.0 <= Mh &&  Mh < 109.5) {
    if (sinba2 < 21.379 - 0.403 * Mh + 0.0019211 * sqr(Mh)) return true;
    else return false;
  }
  else if (109.5 <= Mh &&  Mh < 114.4) {
    if (sinba2 <  1/(60.081 - 0.51624 * Mh)) return true;
    else return false;
  }
  return true;
}

/// from hep-ph/9507294 -- debugged 19/11/04
double ufb3fn(double mu, double htau, double h2, int family, const MssmSoftsusy
	      & temp) { 
  double vufb3 = 0.0;
  /// potential value for these VEVs
  if (fabs(h2) > 
      sqrt(sqr(mu) / (4.0 * sqr(htau)) + 
	   4.0 * temp.displaySoftMassSquared(mLl, family, family) /  
	   (0.6 * sqr(temp.displayGaugeCoupling(1)) +
	    sqr(temp.displayGaugeCoupling(2)))) - fabs(mu) / 
      temp.displayYukawaElement(YE, 3, 3) * 0.5)
    vufb3 = 
      sqr(h2) * (temp.displayMh2Squared() +
		 temp.displaySoftMassSquared(mLl, family, family)) + 
      fabs(mu * h2) / htau * 
      (temp.displaySoftMassSquared(mLl, 3, 3) +
       temp.displaySoftMassSquared(mEr, 3, 3) 
       + temp.displaySoftMassSquared(mLl, family, family)) -
      2.0 * sqr(temp.displaySoftMassSquared(mLl, family, family)) / 
      (0.6 * sqr(temp.displayGaugeCoupling(1)) +
       sqr(temp.displayGaugeCoupling(2)));
  else
    vufb3 = 
      sqr(h2) * temp.displayMh2Squared() + 
      fabs(mu * h2) / htau * 
      (temp.displaySoftMassSquared(mLl, 3, 3) +
       temp.displaySoftMassSquared(mEr, 3, 3)) +  
      1.0 / 8.0 * (0.6 * sqr(temp.displayGaugeCoupling(1)) +
		   sqr(temp.displayGaugeCoupling(2))) * 
      sqr(sqr(h2) + fabs(mu * h2) / htau);
  
  if (PRINTOUT > 1) cout << vufb3 << endl;
  return vufb3;
}

/// For ufb3direction, returns scale at which one-loop corrections are smallest
double getQhat(double inminTol,double eR, double h2, double Lisq, double mx,
		MssmSoftsusy & temp) {
  double oldQhat = -1.0e16;
  int maxNum = 40;
  
  int d; for (d = 1; d <= maxNum; d++)     {
    double qhat = 
      maximum(maximum(maximum(temp.displayGaugeCoupling(2) * eR, 
		     temp.displayGaugeCoupling(2) * fabs(h2)), 
		temp.displayGaugeCoupling(2) * sqrt(fabs(Lisq))),
	   temp.displayYukawaElement(YU, 3, 3) * fabs(h2));
    /// Run all paramaters to that scale
    if (qhat < mx) temp.runto(qhat);
    else temp.runto(mx); 
    if (PRINTOUT > 1) cout << qhat << " ";
    
    if (fabs((qhat - oldQhat) / qhat) < inminTol) return qhat;
    oldQhat = qhat;
  }
  /// Return NOB if no convergence on qhat
  return -6.66e66;
}

/// Difference between two SOFTSUSY objects in and out: EWSB terms only
double sumTol(const MssmSoftsusy & in, const MssmSoftsusy & out, int numTries) {

  drBarPars inforLoops(in.displayDrBarPars()), 
    outforLoops(out.displayDrBarPars());  

  DoubleVector sT(34);
  int k = 1;

  double sTin  = fabs(inforLoops.mh0(1)); double sTout = fabs(outforLoops.mh0(1));
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mA0(1)); sTout = fabs(outforLoops.mA0(1));
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mh0(2)); sTout = fabs(outforLoops.mh0(2));
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mHpm); sTout = fabs(outforLoops.mHpm);
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  int i; for (i=1; i<=3; i++) {
    sTin  = fabs(inforLoops.msnu(i));
    sTout = fabs(outforLoops.msnu(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k++;
  }
  for (i=1; i<=2; i++) {
    sTin = fabs(inforLoops.mch(i));
    sTout = fabs(outforLoops.mch(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k++;
  }
  for (i=1; i<=4; i++) {
    sTin = fabs(inforLoops.mneut(i));
    sTout = fabs(outforLoops.mneut(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k++;
  }
  sTin = fabs(inforLoops.mGluino);
  sTout = fabs(outforLoops.mGluino);
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
  k++;
  int j; for (j=1; j<=3; j++)
    for(i=1; i<=2; i++) {
      sTin = fabs(inforLoops.mu(i, j));
      sTout = fabs(outforLoops.mu(i, j));
      sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
      k++;
      sTin = fabs(inforLoops.md(i, j));
      sTout = fabs(outforLoops.md(i, j));
      sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
      k++;
      sTin = fabs(inforLoops.me(i, j));
      sTout = fabs(outforLoops.me(i, j));
      sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
      k++;
    }
  /// The predicted value of MZ^2 is an absolute measure of how close to a
  /// true solution we are:
  double tbPred = 0.;
  double predictedMzSq = in.displayPredMzSq();
  /// We allow an extra factor of 10 for the precision in the predicted value
  /// of MZ compared to TOLERANCE if the program is struggling and gone beyond
  /// 10 tries - an extra 2 comes from MZ v MZ^2
  if (!in.displayProblem().testSeriousProblem()) {
    sT(k) = 0.5 * 
      fabs(1. - minimum(predictedMzSq, sqr(MZ)) / 
	   maximum(sqr(MZ), predictedMzSq));
    if (numTries > 10) sT(k) *= 0.1;
  }

  return sT.max();
}

/// Prints out what the lsp is
string recogLsp(int temp, int posj) {
  string out;
  switch(temp) {
  case -1: out = "gravitino"; break;
  case 0: out = "neutralino"; break;
  case 1: 
    switch(posj) {
      case 3: out = "stop"; break;
      case 2: out = "scharm"; break;
      case 1: out = "sup"; break;
      } break;
  case 2:
    switch(posj) {
      case 3: out = "sbottom"; break;
      case 2: out = "sstange"; break;
      case 1: out = "sdown"; break;
      } break;
  case 3:
    switch(posj) {
      case 3: out = "stau"; break;
      case 2: out = "smu"; break;
      case 1: out = "selectron"; break;
      } break;
  case 4: out = "chargino"; break;
  case 5: out = "sneutrino"; break;
  case 6: out = "gluino"; break;
  default:
    ostringstream ii;
    ii << "Wrong input to lsp printing routine\n";
    throw ii.str(); break;
  }
  return out;
}

ostream & operator <<(ostream &left, const MssmSoftsusy &s) {
  left << HR << endl;
  left << "Gravitino mass M3/2: " << s.displayGravitino() << endl;
  left << "Msusy: " << s.displayMsusy() << " MW: " << s.displayMw() 
       << " Predicted MZ: " << sqrt(s.displayPredMzSq()) << endl;  
  left << "Data set:\n" << s.displayDataSet();
  left << HR << endl;
  left << s.displaySoftPars();
  left << "t1/v1(MS)=" << s.displayTadpole1Ms() 
       << " t2/v2(MS)=" << s.displayTadpole2Ms() << endl;
  left << HR << "\nPhysical MSSM parameters:\n";
  left << s.displayPhys();
  double mass; int posi, posj, id;
  id = s.lsp(mass, posi, posj);

  /// If the gravitino mass is non-zero, and if it is smaller than the visible
  /// sector LSP mass, make it clear that the particle is the NLSP
  left << "lsp is " << recogLsp(id, posj);
  left << " of mass " << mass << " GeV\n";
  if (s.displayProblem().test()) left << "***** PROBLEM *****" <<
				   s.displayProblem() << " *****" << endl;
  left << HR << endl;

  if (s.displaySetTbAtMX()) left << "Tan beta is set at user defined scale\n";
  if (s.displayAltEwsb()) left << "Alternative EWSB conditions: mu=" 
			       << s.displayMuCond() 
			       << " mA=" << s.displayMaCond() << endl;

  return left;
}
