
/** \file physpars.cpp
   - Project:     SOFTSUSY 
   - File:        physpars.cpp
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include <physpars.h>

const drBarPars & drBarPars::operator=(const drBarPars &s) {
  if (this == &s) return *this;
  mz = s.mz; mw = s.mw;
  mt = s.mt; mb = s.mb; mtau = s.mtau;
  ht = s.ht; hb = s.hb; htau = s.htau;
  ut = s.ut; ub = s.ub; utau = s.utau;
  nBpmz = s.nBpmz; uBpmz = s.uBpmz; vBpmz = s.vBpmz; 
  mnBpmz = s.mnBpmz; mchBpmz = s.mchBpmz;
  setsPhysical(s.displaysPhysical());
  
  return *this;
}

const sPhysical & sPhysical::operator=(const sPhysical &s) {
  if (this == &s) return *this;
  mh0 = s.mh0; mA0 = s.mA0; mH0 = s.mH0; mHpm = s.mHpm;
  msnu = s.msnu; 
  mch = s.mch; mneut = s.mneut; mixNeut = s.mixNeut;
  thetaL = s.thetaL; thetaR = s.thetaR; mGluino = s.mGluino;
  thetat = s.thetat; thetab = s.thetab; thetatau = s.thetatau;
  mu = s.mu; md = s.md; me = s.me; thetaH = s.thetaH;
  return *this;
}

// a should be in C convention ie start from index zero
void sPhysical::display(double *a) const {
  a[0] = mh0; a[1] = mA0;
  a[2] = mH0; a[3] = mHpm;

  a[4] = msnu.display(1); a[5] = msnu.display(2); a[6] = msnu.display(3);

  a[7] = mch.display(1); a[8] = mch.display(2);

  a[9] = mneut.display(1); a[10] = mneut.display(2); 
  a[11] = mneut.display(3); a[12] = mneut.display(4);

  a[13] = mGluino;

  int i, j, k = 13; 
  for (i=1; i<=4; i++)
    for (j=1; j<=4; j++) {
      k++;
      a[k] = mixNeut.display(i, j);  
    }

  a[30] = thetaL; a[31] = thetaR; 
  a[32] = thetat; a[33] = thetab; a[34] = thetatau;

  k = 34;
  for (i=1; i<=2; i++)
    for (j=1; j<=3; j++) {
      k++;
      a[k] = mu.display(i, j);
      a[k+6] = md.display(i, j);
      a[k+12] = me.display(i, j);
    }  

  a[53] = thetaH;
}

#define HR "---------------------------------------------------------------\n"

std::ostream & operator <<(std::ostream & left, const drBarPars &s) {
  left << s.displaysPhysical();
  left << "BPMZ conventions, N" << s.nBpmz << "U" << s.uBpmz << "V" 
       << s.vBpmz;
  left << "mt: "  << s.mt << " mb: " << s.mb << " mtau: " << s.mtau << endl;
  left << "ht: "  << s.ht << " hb: " << s.hb << " htau: " << s.htau << endl;
  left << "Ut: "  << s.ut << " Ub: " << s.ub << " Utau: " << s.utau << endl;
  left << "mz: "  << s.mz << " mw: " << s.mw << endl;

  return left;
}

std::ostream & operator <<(std::ostream & left, const sPhysical &s) {
  left << "mh^0: " << s.mh0 << " mA^0: " << s.mA0
       << " mH^0: " << 
    s.mH0 << " mH^+-: " << s.mHpm << "\n";
  left << "alpha: " << s.thetaH << "\n";
  left << "sneutrinos" << s.msnu; 
  left << "mU~" << s.mu << "mD~" << s.md << "mE~" << s.me;
  left << "thetat: " << s.thetat << " thetab: " << s.thetab << 
    " thetatau: " << s.thetatau << "\n";
  left << "mGluino:  " << s.mGluino << "\n";
  left << "charginos" << s.mch;
  left << "thetaL: " << s.thetaL << " thetaR: " << s.thetaR << "\n";
  left << "neutralinos" << s.mneut;
  left << "neutralino mixing matrix " << s.mixNeut;
  return left;
}

std::istream & operator >>(std::istream & left, sPhysical &s) {
  char c[70];
  left >> c >> c >> c >> c;
  left >> c >> s.mh0 >> c >> s.mA0
       >> c >> s.mH0 >> c >> s.mHpm;
  left >> c >> s.thetaH;
  left >> s.msnu; 
  left >> c >> s.mu >> c >> s.md >> c >> s.me;
  left >> c >> s.thetat >> c >> s.thetab >> 
    c >> s.thetatau;
  left >> c >> s.mGluino;
  left >> s.mch;
  left >> c >> s.thetaL >> c >> s.thetaR;
  left >> s.mneut;
  left >> c >> c >> c >> c >> s.mixNeut;
  return left;
}

#undef HR

ostream & operator <<(ostream &st, const sProblem & p) {
  if (!p.test()) return st;
  st << "[ ";
  if (p.mgutOutOfBounds) st << "GUT scale too high or too low ";
  if (p.badConvergence) st << "No acceptable solution found ";
  if (p.irqfp) st << "Quasi-fixed point breached ";
  if (p.noMuConvergence) st << "No mu convergence ";
  if (p.noRhoConvergence) st << "No rho convergence ";
  if (p.nonperturbative) st << "Non-perturbative ";
  if (p.noConvergence) st << "No convergence ";
  if (p.tachyon) {
    /*    switch (p.tachyonType) {
          case selectron:
    case smuon:
      if (p.tachyonType == )
      }*/
    st << tachyonNames[p.tachyon] << " tachyon ";
  }
  if (p.muSqWrongSign) st << "MuSqWrongsign ";
  if (p.m3sq) st << "m3sq-problem ";
  if (p.higgsUfb) st << "Higgs potential ufb ";
  if (p.inaccurateHiggsMass) st << "Inaccurate Higgs mass ";
  if (p.problemThrown) st << "Numerical problemThrown ";
  st << "]";
  return st;
}

const sProblem & sProblem::operator=(const sProblem &s) {
  if (this == &s) return *this;
  mgutOutOfBounds = s.mgutOutOfBounds;
  irqfp = s.irqfp;
  badConvergence = s.badConvergence;
  noMuConvergence = s.noMuConvergence;
  noRhoConvergence = s.noRhoConvergence;
  nonperturbative = s.nonperturbative;
  noConvergence = s.noConvergence;
  tachyon = s.tachyon;
  muSqWrongSign = s.muSqWrongSign;
  higgsUfb = s.higgsUfb;
  m3sq = s.m3sq;
  problemThrown = s.problemThrown;
  return *this;
}

// Returns mixing matrix o and neutralino masses mn in the MPZ convention
// (hep-ph/9606211), n is 4 by 4 and mneut is 1->4.
void drBarPars::mpzNeutralinos() { 
  // We want to change the PHASES of the neutralino mixing matrix in order to
  // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang

  DoubleVector temp(mneut);
  
  ComplexMatrix K(4, 4);
  int i; for (i=1; i<=4; i++) 
    if (mneut.display(i) < 0.0) K(i, i) = Complex(0.0, 1.0);
    else
      K(i, i) = Complex(1.0, 0.0);
  
  mnBpmz = temp.apply(fabs);
  nBpmz = K.hermitianConjugate() * mixNeut.transpose();
}

// Returns mixing matrices u,v and neutralino masses mneut in the MPZ
// convention (hep-ph/9606211),  u+v are (2,2) and mch is 1->2.
void drBarPars::mpzCharginos() {
  // We want to change the PHASES of the neutralino mixing matrix in order to
  // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
  ComplexMatrix u(2, 2), v(2, 2);
  positivise(thetaL, thetaR, mch, u, v);
  uBpmz = u; vBpmz = v;
  mchBpmz = mch.apply(fabs); 
}
