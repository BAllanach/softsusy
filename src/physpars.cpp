
/** \file physpars.cpp
   - Project:     SOFTSUSY 
   - File:        physpars.cpp
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "physpars.h"

namespace softsusy {

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
  mh0 = s.mh0; mA0 = s.mA0; mHpm = s.mHpm;
  msnu = s.msnu; 
  mch = s.mch; mneut = s.mneut; mixNeut = s.mixNeut;
  thetaL = s.thetaL; thetaR = s.thetaR; mGluino = s.mGluino;
  thetat = s.thetat; thetab = s.thetab; thetatau = s.thetatau;
  thetamu = s.thetamu, mu = s.mu; md = s.md; me = s.me; thetaH = s.thetaH;  
  thetaA0 = s.thetaA0; mixh0 = s.mixh0; 
  return *this;
}

// a should be in C convention ie start from index zero
void sPhysical::display(double *a) const {
   const DoubleVector v = display();
   for (std::size_t i = 0; i < v.size(); i++)
      a[i] = v(i+1);
}

int sPhysical::size() const
{
   return mh0.size() + mA0.size() + 1 // Higgses
        + msnu.size() + mch.size() + mneut.size()
        + 1                                 // gluino
        + mixNeut.size() + 8 + mixh0.size() // mixing angles
        + mu.size() + md.size() + me.size();
}

DoubleVector sPhysical::display() const
{
   const int len = size();
   DoubleVector v(len);
   int k = 1;

   for (std::size_t i = 1; i <= mh0.size()    ; i++) v(k++) = mh0(i);
   for (std::size_t i = 1; i <= mA0.size()    ; i++) v(k++) = mA0(i);
   v(k++) = mHpm;
   for (std::size_t i = 1; i <= msnu.size()   ; i++) v(k++) = msnu(i);
   for (std::size_t i = 1; i <= mch.size()    ; i++) v(k++) = mch(i);
   for (std::size_t i = 1; i <= mneut.size()  ; i++) v(k++) = mneut(i);
   v(k++) = mGluino;
   for (int i = 1; i <= mixNeut.displayRows(); i++)
      for (int j = 1; j <= mixNeut.displayCols(); j++)
         v(k++) = mixNeut(i,j);
   v(k++) = thetaL;
   v(k++) = thetaR;
   v(k++) = thetat;
   v(k++) = thetab;
   v(k++) = thetatau;
   v(k++) = thetamu;
   for (int i = 1; i <= mu.displayRows(); i++)
      for (int j = 1; j <= mu.displayCols(); j++)
         v(k++) = mu(i,j);
   for (int i = 1; i <= md.displayRows(); i++)
      for (int j = 1; j <= md.displayCols(); j++)
         v(k++) = md(i,j);
   for (int i = 1; i <= me.displayRows(); i++)
      for (int j = 1; j <= me.displayCols(); j++)
         v(k++) = me(i,j);
   v(k++) = thetaH;
   for (int i = 1; i <= mixh0.displayRows(); i++)
      for (int j = 1; j <= mixh0.displayCols(); j++)
         v(k++) = mixh0(i,j);
   v(k++) = thetaA0;

   if (len+1 != k)
      throw std::string("Bug: access outside array boundaries");

   return v;
}

void sPhysical::set(const DoubleVector& v)
{
   if (size() != v.size())
      throw std::string("Bug: vector is too small!");

   std::size_t k = 1;

   for (std::size_t i = 1; i <= mh0.size()    ; i++) mh0(i) = v(k++);
   for (std::size_t i = 1; i <= mA0.size()    ; i++) mA0(i) = v(k++);
   mHpm = v(k++);
   for (std::size_t i = 1; i <= msnu.size()   ; i++) msnu(i) = v(k++);
   for (std::size_t i = 1; i <= mch.size()    ; i++) mch(i) = v(k++);
   for (std::size_t i = 1; i <= mneut.size()  ; i++) mneut(i) = v(k++);
   mGluino = v(k++);
   for (int i = 1; i <= mixNeut.displayRows(); i++)
      for (int j = 1; j <= mixNeut.displayCols(); j++)
         mixNeut(i,j) = v(k++);
   thetaL = v(k++);
   thetaR = v(k++);
   thetat = v(k++);
   thetab = v(k++);
   thetatau = v(k++);
   thetamu = v(k++);
   for (int i = 1; i <= mu.displayRows(); i++)
      for (int j = 1; j <= mu.displayCols(); j++)
         mu(i,j) = v(k++);
   for (int i = 1; i <= md.displayRows(); i++)
      for (int j = 1; j <= md.displayCols(); j++)
         md(i,j) = v(k++);
   for (int i = 1; i <= me.displayRows(); i++)
      for (int j = 1; j <= me.displayCols(); j++)
         me(i,j) = v(k++);
   thetaH = v(k++);
   for (int i = 1; i <= mixh0.displayRows(); i++)
      for (int j = 1; j <= mixh0.displayCols(); j++)
         mixh0(i,j) = v(k++);
   thetaA0 = v(k++);
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
  left << "mh^0: " << s.mh0 << "mA^0: " << s.mA0
       << "mH^+-: " << s.mHpm << "\n";
  left << "alpha: " << s.thetaH << "\n";
  left << "sneutrinos" << s.msnu; 
  left << "mU~" << s.mu << "mD~" << s.md << "mE~" << s.me;
  left << "thetat: " << s.thetat << " thetab: " << s.thetab << 
     " thetatau: " << s.thetatau <<"\nthetamu: "  << s.thetamu << "\n";
  left << "mGluino:  " << s.mGluino << "\n";
  left << "charginos" << s.mch;
  left << "thetaL: " << s.thetaL << " thetaR: " << s.thetaR << "\n";
  left << "neutralinos" << s.mneut;
  left << "neutralino mixing matrix " << s.mixNeut;
  left << "CP even Higgs mixing matrix: " << s.mixh0;
  left << "CP odd mixing angle: " << s.thetaA0;
  return left;
}

std::istream & operator >>(std::istream & left, sPhysical &s) {
  string c;
  left >> c >> c >> c >> c;
  left >> c >> s.mh0 >> c >> s.mA0
       >> c >> s.mHpm;
  left >> c >> s.thetaH;
  left >> s.msnu; 
  left >> c >> s.mu >> c >> s.md >> c >> s.me;
  left >> c >> s.thetat >> c >> s.thetab >> 
     c >> s.thetatau >> c >> s.thetamu;
  left >> c >> s.mGluino;
  left >> s.mch;
  left >> c >> s.thetaL >> c >> s.thetaR;
  left >> s.mneut;
  left >> c >> c >> c >> c >> s.mixNeut;
  left >> c >> c >> c >> c >> s.mixh0;
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
  if (p.tachyon) st << tachyonNames[p.tachyon] << " tachyon ";
  if (p.tachyonWarning) st << tachyonNames[p.tachyonWarning] 
			   << " is tree-level tachyon at MZ ";
  if (p.muSqWrongSign) st << "MuSqWrongsign ";
  if (p.m3sq) st << "m3sq-problem ";
  if (p.higgsUfb) st << "Higgs potential ufb ";
  if (p.notGlobalMin) st << "Not in global min of Higgs potential " ;
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
  tachyonWarning = s.tachyonWarning;
  muSqWrongSign = s.muSqWrongSign;
  higgsUfb = s.higgsUfb;
  notGlobalMin = s.notGlobalMin;
  m3sq = s.m3sq;
  problemThrown = s.problemThrown;
  return *this;
}

// Returns mixing matrix o and neutralino masses mn in the MPZ convention
// (hep-ph/9606211), n is 4 by 4 and mneut is 1->4.
void drBarPars::mpzNeutralinos() { 
  // We want to change the PHASES of the neutralino mixing matrix in order to
  // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
  const int rank = mneut.displayEnd();
  DoubleVector temp(mneut);
  
  ComplexMatrix K(rank, rank);
  int i; for (i=1; i<=rank; i++) 
    if (mneut.display(i) < 0.0) K(i, i) = Complex(0.0, -1.0);
    else
      K(i, i) = Complex(1.0, 0.0);
  
  mnBpmz = temp.apply(fabs);
  nBpmz = K * mixNeut.transpose();
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

} // namespace softsusy

double sTfn(double sTins, double sTouts) {
  double sTin  = fabs(sTins);
  double sTout = fabs(sTouts);
  if (sTin < 1. && sTout < 1.) return fabs(sTin - sTout);
  else return fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
}

/// LCT: Difference between two drBarPars objects
void sumTol(const softsusy::drBarPars & a, const softsusy::drBarPars & b, DoubleVector & sT) {
  int k = 1;

  sT(k) = sTfn(a.mGluino, b.mGluino); k++;
  int i; for (i=1; i<=a.mh0.displayEnd(); i++) {
    sT(k) = sTfn(a.mh0(i), b.mh0(i)); k++;
  }
  for (i=1; i<=a.mA0.displayEnd(); i++) {
    sT(k) = sTfn(a.mA0(i), b.mA0(i)); k++;
  }
  sT(k) = sTfn(a.mHpm, b.mHpm); k++;
  for (i=1; i<=3; i++) {
    sT(k) = sTfn(a.msnu(i), b.msnu(i)); k++;
  }
  for (i=1; i<=2; i++) {
    sT(k) = sTfn(a.mch(i), b.mch(i)); k++;
  }
  for (i=1; i<=a.mneut.displayEnd(); i++) {
    sT(k) = sTfn(a.mneut(i), b.mneut(i)); k++;
  }
  int j; for (j=1; j<=3; j++)
    for(i=1; i<=2; i++) {
      sT(k) = sTfn(a.mu(i, j), b.mu(i, j)); k++;
      sT(k) = sTfn(a.md(i, j), b.md(i, j)); k++;
      sT(k) = sTfn(a.me(i, j), b.me(i, j)); k++;
    }
}
