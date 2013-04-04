
/** \file physpars.cpp
   - Project:     SOFTSUSY 
   - File:        physpars.cpp
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include <nmssmphyspars.h>

const NMdrBarPars & NMdrBarPars::operator=(const NMdrBarPars &s) {
  if (this == &s) return *this;
  mz = s.mz; mw = s.mw;
  mt = s.mt; mb = s.mb; mtau = s.mtau;
  ht = s.ht; hb = s.hb; htau = s.htau;
  ut = s.ut; ub = s.ub; utau = s.utau;
  nBpmz = s.nBpmz; uBpmz = s.uBpmz; vBpmz = s.vBpmz; 
  mnBpmz = s.mnBpmz; mchBpmz = s.mchBpmz;
  setnmsPhysical(s.displaynmsPhysical());
  
  return *this;
}

const nmsPhysical & nmsPhysical::operator=(const nmsPhysical &s) {
  if (this == &s) return *this;
  mh0 = s.mh0; mA0 = s.mA0; mHpm = s.mHpm;
  mH1 = s.mH1; mH2 = s.mH2; mA1 = s.mA1; mA2 = s.mA2;
  msnu = s.msnu; 
  mch = s.mch; mneut = s.mneut; mixNeut = s.mixNeut;
  thetaL = s.thetaL; thetaR = s.thetaR; mGluino = s.mGluino;
  thetat = s.thetat; thetab = s.thetab; thetatau = s.thetatau;
  mu = s.mu; md = s.md; me = s.me; thetaH = s.thetaH;
  t1OV1Ms = s.t1OV1Ms; t2OV2Ms = s.t2OV2Ms;
  t1OV1Ms1loop = s.t1OV1Ms1loop; t2OV2Ms1loop = s.t2OV2Ms1loop;
  return *this;
}

// a should be in C convention ie start from index zero
void nmsPhysical::display(double *a) const {
  std::size_t k = 0;
  mh0.fillArray(a,k);   k += mh0.size();
  mA0.fillArray(a,k);   k += mA0.size();
  a[k++] = mHpm;
  msnu.fillArray(a,k);  k += msnu.size();
  mch.fillArray(a,k);   k += mch.size();
  mneut.fillArray(a,k); k += mneut.size();
  a[k++] = mGluino;
  mixNeut.fillArray(a,k); k += mixNeut.size();
  a[k++] = thetaL;
  a[k++] = thetaR;
  a[k++] = thetat;
  a[k++] = thetab;
  a[k++] = thetatau;
  mu.fillArray(a,k); k += mu.size();
  md.fillArray(a,k); k += md.size();
  me.fillArray(a,k); k += me.size();
  a[k++] = thetaH;
  a[k++] = mH1;
  a[k++] = mH2;
  a[k++] = mA1;
  a[k++] = mA2;
  a[k++] = thetaHP;
}
#define HR "---------------------------------------------------------------\n"

std::ostream & operator <<(std::ostream & left, const NMdrBarPars &s) {
  left << s.displaynmsPhysical();
  left << "BPMZ conventions, N" << s.nBpmz << "U" << s.uBpmz << "V" 
       << s.vBpmz;
  left << "mt: "  << s.mt << " mb: " << s.mb << " mtau: " << s.mtau << endl;
  left << "ht: "  << s.ht << " hb: " << s.hb << " htau: " << s.htau << endl;
  left << "Ut: "  << s.ut << " Ub: " << s.ub << " Utau: " << s.utau << endl;
  left << "mz: "  << s.mz << " mw: " << s.mw << endl;

  return left;
}

std::ostream & operator <<(std::ostream & left, const nmsPhysical &s) {
  left << "mh^0: " << s.mh0(1) << " mA^0: " << s.mA0(1)
       << " mH^0: " << 
    s.mh0(2) << " mH^+-: " << s.mHpm << "\n";
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
  left << "mH1 = " << s.mH1 << "mH2 = "  << s.mH2 << endl;
  left << "mA1 = " << s.mA1 << "mA2 = "  << s.mA2 << endl;
  
 return left;
}

std::istream & operator >>(std::istream & left, nmsPhysical &s) {
  char c[70];
  left >> c >> c >> c >> c;
  left >> c >> s.mh0(1) >> c >> s.mA0(1)
       >> c >> s.mh0(2) >> c >> s.mHpm;
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

ostream & operator <<(ostream &st, const NMsProblem& p) {
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

const NMsProblem& NMsProblem::operator=(const NMsProblem&s) {
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
void NMdrBarPars::mpzNeutralinos() { 
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
void NMdrBarPars::mpzCharginos() {
  // We want to change the PHASES of the neutralino mixing matrix in order to
  // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
  ComplexMatrix u(2, 2), v(2, 2);
  positivise(thetaL, thetaR, mch, u, v);
  uBpmz = u; vBpmz = v;
  mchBpmz = mch.apply(fabs); 
}
