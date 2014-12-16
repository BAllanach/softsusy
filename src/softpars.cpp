
/** \file softpars.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: includes all MSSM parameters in the Lagrangian

*/

#include "softpars.h"

const MssmSoftPars & MssmSoftPars::operator=(const MssmSoftPars & s) {
  if (this == &s) return *this;
  mGaugino = s.mGaugino;
  ua = s.ua;
  da = s.da;
  ea = s.ea;
  m32 = s.m32;  
  mQLsq = s.mQLsq;
  mURsq = s.mURsq;
  mDRsq = s.mDRsq;
  mLLsq = s.mLLsq;
  mSEsq = s.mSEsq;
  m3sq = s.m3sq;
  mH1sq = s.mH1sq;
  mH2sq = s.mH2sq;
  return *this;
}

const DoubleMatrix & MssmSoftPars::displayTrilinear(trilinears k) const {
  switch(k) {
  case UA: return ua; break;
  case DA: return da; break;
  case EA: return ea; break;
  default: 
    ostringstream ii;
    ii << "In SoftPars::displayTrilinear, called with illegal argument";
    ii << " " << k << endl;
    throw ii.str();
    break;
  }
}

double MssmSoftPars::displayTrilinear(trilinears k, int i, int j) 
  const {
  switch(k) {
  case UA: return ua.display(i, j); break;
  case DA: return da.display(i, j); break;
  case EA: return ea.display(i, j); break;
  default: 
    ostringstream ii;
    ii << "In SoftPars::displayTrilinear, called with illegal argument";
    ii << " " << k << endl;
    throw ii.str();
    break;
  }
}

const DoubleMatrix & MssmSoftPars::displaySoftMassSquared(softMasses k) const {
  switch(k) {
  case mQl: return mQLsq; break;
  case mUr: return mURsq; break;
  case mDr: return mDRsq; break;
  case mLl: return mLLsq; break;
  case mEr: return mSEsq; break;
  default: 
    ostringstream ii;
    ii << "SoftPars::displaySoftMassSquared with illegal argument ";
    ii << k << endl;
    throw ii.str();
    break;
  }
}

double MssmSoftPars::displaySoftMassSquared(softMasses k, int i, int j) 
  const {
  switch(k) {
  case mQl: return mQLsq.display(i, j); break;
  case mUr: return mURsq.display(i, j); break;
  case mDr: return mDRsq.display(i, j); break;
  case mLl: return mLLsq.display(i, j); break;
  case mEr: return mSEsq.display(i, j); break;
  default: 
    ostringstream ii;
    ii << "SoftPars::displaySoftMassSquared with illegal argument ";
    ii << k << endl;
    throw ii.str();
    break;
  }
}

void MssmSoftPars::setSoftMassElement(softMasses k, int i, int j, 
					  double f) { 
  switch(k) {
  case mQl: mQLsq(i, j) = f; break;
  case mUr: mURsq(i, j) = f; break;
  case mDr: mDRsq(i, j) = f; break;
  case mLl: mLLsq(i, j) = f; break;
  case mEr: mSEsq(i, j) = f; break;
  }
}

void MssmSoftPars::setSoftMassMatrix(softMasses k, const DoubleMatrix & m) { 
  switch(k) {
  case mQl: mQLsq = m; break;
  case mUr: mURsq = m; break;
  case mDr: mDRsq = m; break;
  case mLl: mLLsq = m; break;
  case mEr: mSEsq = m; break;
  }
}

void MssmSoftPars::setTrilinearMatrix(trilinears k, const DoubleMatrix & m) { 
  switch(k) {
  case UA: ua = m; break;
  case DA: da = m; break;
  case EA: ea = m; break;
  }
}

void MssmSoftPars::setTrilinearElement(trilinears k, int i, int j, 
					  double m) { 
  switch(k) {
  case UA: ua(i, j) = m; break;
  case DA: da(i, j) = m; break;
  case EA: ea(i, j) = m; break;
  }
}

void MssmSoftPars::setAllGauginos(const DoubleVector & v) { 
  if (v.displayStart() != 1 || v.displayEnd() !=3) {
    ostringstream ii;
    ii << "Initialising SoftPars::setAllGauginos with vector"
	 << v;
    throw ii.str();
  }
  mGaugino = v; 
}




/* 
   Give it a SUSY object and a value of M3/2, and it will return a soft
   object with AMSB soft breaking terms. Note that the sleptons will be
   tachyonic, ie nothing has been done to fix that problem.
   Note that in the following, we are neglecting all Yukawa couplings except
   that of the third family.
   
   THE CURRENT STATE OF PLAY:
   Two loop additions are possible, but a pain.
   */
/*void MssmSoftPars::addAmsb(double maux) {
  MssmSusy run(displayMssmSusy());
  const double ONEO16pisq = 1.0 / (16. * sqr(PI));
  const double ONEO16pif = sqr(ONEO16pisq);
  double     g1   = run.displayGaugeCoupling(1), 
    g2   = run.displayGaugeCoupling(2), 
    g3   = run.displayGaugeCoupling(3); 

  // For calculational brevity: 
  static sBrevity a;
  static MssmSusy dsb;
  
  // calculate derivatives for full SUSY spectrum. Brevity calculations come
  // out encoded in a
  dsb = run.beta(a);

  /// three family version -- will be out by some two-loop terms. See
  /// hep-ph/9904378 
  mQLsq = mQLsq + sqr(maux) * 
    (ONEO16pif * (-11. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2)) + 
		  8.0 * sqr(sqr(g3))) +
     ONEO16pisq * 0.5 *  
     (dsb.displayYukawaMatrix(YD) * displayYukawaMatrix(YD).transpose() +
      displayYukawaMatrix(YD) * dsb.displayYukawaMatrix(YD).transpose() +
      displayYukawaMatrix(YU) * dsb.displayYukawaMatrix(YU).transpose() +
      dsb.displayYukawaMatrix(YU) * displayYukawaMatrix(YU).transpose()));
  mLLsq = mLLsq + sqr(maux) * 
    (ONEO16pif * (-99. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2))) + 
     ONEO16pisq * 0.5 * 
     (dsb.displayYukawaMatrix(YE) * displayYukawaMatrix(YE).transpose() +
      displayYukawaMatrix(YE) * dsb.displayYukawaMatrix(YE).transpose())); 
  mURsq = mURsq + sqr(maux) * 
    (ONEO16pif * (-88. / 25. * sqr(sqr(g1)) + 8 * sqr(sqr(g3))) + 
     ONEO16pisq * 
     (dsb.displayYukawaMatrix(YU).transpose() * displayYukawaMatrix(YU) +
      displayYukawaMatrix(YU).transpose() * dsb.displayYukawaMatrix(YU)));
  mDRsq = mDRsq + sqr(maux) * 
    (ONEO16pif * (-22. / 25. * sqr(sqr(g1)) + 8 * sqr(sqr(g3))) + 
     ONEO16pisq * 
     (dsb.displayYukawaMatrix(YD).transpose() * displayYukawaMatrix(YD) +
      displayYukawaMatrix(YD).transpose() * dsb.displayYukawaMatrix(YD)));
  mSEsq = mSEsq + sqr(maux) * 
    (ONEO16pif * (-198. / 25. * sqr(sqr(g1))) + 
     ONEO16pisq * 
     (dsb.displayYukawaMatrix(YE).transpose() * displayYukawaMatrix(YE) +
      displayYukawaMatrix(YE).transpose() * dsb.displayYukawaMatrix(YE)));
  mH1sq = mH1sq + sqr(maux) * 
    (ONEO16pif * (-99. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2))) + 
     ONEO16pisq * (3.0 * (dsb.displayYukawaMatrix(YD) *
			  displayYukawaMatrix(YD).transpose()).trace() 
		   + (dsb.displayYukawaMatrix(YE) * 
		      displayYukawaMatrix(YE).transpose()).trace()));
  mH2sq = mH2sq + sqr(maux) * 
    (ONEO16pif * (-99. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2))) + 
     ONEO16pisq * (3.0 * (dsb.displayYukawaMatrix(YU) * 
			  displayYukawaMatrix(YU).transpose()).trace()));

  // AMSB spectrum
  DoubleVector amsbGaugino(3);
  DoubleMatrix temp(3, 3);

  int i; for (i=1; i<=3; i++)
    amsbGaugino(i) = dsb.displayGaugeCoupling(i) * maux / 
      run.displayGaugeCoupling(i);
  mGaugino = mGaugino + amsbGaugino;
  
  ua = ua - dsb.displayYukawaMatrix(YU) * maux;
  da = da - dsb.displayYukawaMatrix(YD) * maux;
  ea = ea - dsb.displayYukawaMatrix(YE) * maux;
  
  m32 = maux;
  
  //  u1R_PQflip();
  return;
}


void MssmSoftPars::u1R_PQflip() {
  setSusyMu(-displaySusyMu());
  mGaugino = -1. * mGaugino;
  ua = -1. * ua;
  da = -1. * da;
  ea = -1. * ea;
  }*/

// Reads in universal boundary conditions at the current scale:
// m0, M1/2, A0, B and sign of mu
/*void MssmSoftPars::universal(double m0,  double m12,  double a0,  double mu,
			      double m3sq) {
  standardSugra(m0, m12, a0);  
  setSusyMu(mu);
  setM3Squared(m3sq);
  }*/

void MssmSoftPars::universalScalars(double m0) {
  // scalar masses
  DoubleMatrix ID(3, 3), mm0(3, 3);
  int i; for (i=1; i<=3; i++) ID(i, i) = 1.0;
  mm0 = ID * sqr(m0);
  setSoftMassMatrix(mQl, mm0); setSoftMassMatrix(mUr, mm0);
  setSoftMassMatrix(mDr, mm0); setSoftMassMatrix(mLl, mm0);
  setSoftMassMatrix(mEr, mm0);
  setMh1Squared(sqr(m0)); setMh2Squared(sqr(m0));
}

void MssmSoftPars::universalGauginos(double m12) {  
  // gaugino masses
  int i; for (i=1; i<=3; i++) setGauginoMass(i, m12);
}

/*
void MssmSoftPars::universalTrilinears(double a0)  {  
  // trilinears
  setTrilinearMatrix(UA, a0 * displayYukawaMatrix(YU)); 
  setTrilinearMatrix(DA, a0 * displayYukawaMatrix(YD));
  setTrilinearMatrix(EA, a0 * displayYukawaMatrix(YE));
} 

// Input m0, NOT m0 squared.
void MssmSoftPars::standardSugra(double m0,  double m12, double a0) {
  universalScalars(m0);
  universalGauginos(m12);
  universalTrilinears(a0);
}
*/

#define HR "---------------------------------------------------------------\n"

ostream & operator <<(ostream &left, const MssmSoftPars &s) {
  left << " UA" << s.displayTrilinear(UA) 
       << " UD" << s.displayTrilinear(DA) 
       << " UE" << s.displayTrilinear(EA); 
  left << " mQLsq" << s.displaySoftMassSquared(mQl) 
       << " mURsq" << s.displaySoftMassSquared(mUr) 
       << " mDRsq" << s.displaySoftMassSquared(mDr) 
       << " mLLsq" << s.displaySoftMassSquared(mLl) 
       << " mSEsq" << s.displaySoftMassSquared(mEr);
  left << "m3sq: " << s.displayM3Squared() << " mH1sq: " <<
    s.displayMh1Squared() << " mH2sq: " << s.displayMh2Squared() << '\n';
  left << "Gaugino masses" << s.displayGaugino();
  return left;
}

#undef HR

void MssmSoftPars::inputSoftParsOnly() {
  string c;

  cin >> c >> c >> c >> c >> c;
  cin >> c >> ua
       >> c >> da
       >> c >> ea;
  cin >> c >> mQLsq
       >> c >> mURsq
       >> c >> mDRsq
       >> c >> mLLsq
       >> c >> mSEsq;
  cin >> c >> m3sq >> c >> mH1sq >> c >> mH2sq;
  cin >> c >> mGaugino; 
}

istream & operator >>(istream &left, MssmSoftPars &s) {
  string c;

  left >> c >> c >> c >> c >> c >> c >> c;
  DoubleMatrix ua(3, 3), da(3, 3), ea(3, 3);
  left >> c >> ua
       >> c >> da
       >> c >> ea;
  s.setTrilinearMatrix(UA, ua);
  s.setTrilinearMatrix(DA, da);
  s.setTrilinearMatrix(EA, ea);
  DoubleMatrix mqlsq(3, 3), mursq(3, 3), mdrsq(3, 3), mllsq(3, 3), mersq(3, 3);
  left >> c >> mqlsq
       >> c >> mursq
       >> c >> mdrsq
       >> c >> mllsq
       >> c >> mersq;
  s.setSoftMassMatrix(mQl, mqlsq); 
  s.setSoftMassMatrix(mUr, mursq); 
  s.setSoftMassMatrix(mDr, mdrsq); 
  s.setSoftMassMatrix(mLl, mllsq); 
  s.setSoftMassMatrix(mEr, mersq); 
  double m3sq, mh1sq, mh2sq;
  left >> c >> m3sq >> c >> mh1sq >> c >> mh2sq;
  s.setM3Squared(m3sq); s.setMh1Squared(mh1sq); s.setMh2Squared(mh2sq);
  DoubleVector mg(3);
  left >> c >> mg; 
  s.setAllGauginos(mg);
  return left;
}

// Boundary conditions to be applied at messenger scale for Gauge mediated
// SUSY breaking (see hep-ph/9703211 for example)
/*
void MssmSoftPars::minimalGmsb(int n5, double LAMBDA, double mMess, 
			       double cgrav) {

// Modified thresholds by JEL 1-26-04 to accomodate numerical infinities

  const double epstol = 1.0e-4;
  double x = LAMBDA / mMess;

  double f, g;

  if(fabs(x) < epstol) { /// hep-ph/9801271
    g = 1.0 + x*x/6.0 + sqr(x*x)/15.0;
    f = 1.0 + x*x/36.0 - 11.0*sqr(x*x)/450.0;
  }
  else if(fabs(x-1.0) < 0.0001) {
    g  =  log(4.0);
    f  = -sqr(PI)/6.0 + log(4.0) + 0.5*sqr(log(4.0));
    g -=  0.0008132638905771205626;
    f -= -0.0049563838821509165200;
  }
  else {
    g = 1.0 / sqr(x) * 
      ((1.0 + x) * log(1.0 + x) + (1.0 - x) * log(1.0 - x));
    f = (1.0 + x) / sqr(x) * 
    (log(1.0 + x) - 2.0 * dilog(x / (1.0 + x)) + 0.5 * 
     dilog(2.0 * x / (1.0 + x))) + 
     (1.0 - x) / sqr(x) * (log(1.0 - x) - 2.0 * dilog(-x / (1.0 - x)) +
			 0.5 * dilog(-2.0 * x / (1.0 - x)));
  }

  double n5d = double(n5);

  /// There is a relative minus in the mGMSB conditions for gaugino masses,
  /// since these equations are for L=-M/2 gaugino gaugino. See hep-ph/9801271:
  /// BCA 27/7/12
  double m1, m2, m3;
  m1 = n5d * sqr(displayGaugeCoupling(1)) / (16.0 * sqr(PI)) * LAMBDA * g; 
  m2 = n5d * sqr(displayGaugeCoupling(2)) / (16.0 * sqr(PI)) * LAMBDA * g; 
  m3 = n5d * sqr(displayGaugeCoupling(3)) / (16.0 * sqr(PI)) * LAMBDA * g; 
  setGauginoMass(1, m1);   setGauginoMass(2, m2);   setGauginoMass(3, m3);

  setM32(2.37e-19 * LAMBDA * mMess * cgrav);

  double g1f = sqr(sqr(displayGaugeCoupling(1)));
  double g2f = sqr(sqr(displayGaugeCoupling(2)));
  double g3f = sqr(sqr(displayGaugeCoupling(3)));

  double mursq, mdrsq, mersq, mqlsq, mllsq;
  mursq = 2.0 * f * sqr(LAMBDA) * n5d * 
    (4.0 / 3.0 * g3f + 0.6 * 4.0 / 9.0 * g1f) 
    / sqr(16.0 * sqr(PI));
  mdrsq = 2.0 * f * sqr(LAMBDA) * n5d * 
    (4.0 / 3.0 * g3f + 0.6 * 1.0 / 9.0 * g1f) 
    / sqr(16.0 * sqr(PI));
  mersq = 2.0 * f * sqr(LAMBDA) * n5d * 
    (0.6 * g1f) 
    / sqr(16.0 * sqr(PI));
  mqlsq = 2.0 * f * sqr(LAMBDA) * n5d * 
    (4.0 / 3.0 * g3f + 0.75 * g2f + 0.6 * g1f / 36.0) 
    / sqr(16.0 * sqr(PI));
  mllsq = 2.0 * f * sqr(LAMBDA) * n5d * 
    (                  0.75 * g2f + 0.6 * 0.25 * g1f) 
    / sqr(16.0 * sqr(PI));

  // You need Higgs masses too!

  DoubleMatrix id(3, 3);
  id(1, 1) = 1.0; id(2, 2) = 1.0; id(3, 3) = 1.0;

  setSoftMassMatrix(mQl, mqlsq * id);
  setSoftMassMatrix(mUr, mursq * id);
  setSoftMassMatrix(mDr, mdrsq * id);
  setSoftMassMatrix(mLl, mllsq * id);  
  setMh1Squared(mllsq);
  setMh2Squared(mllsq);
  setSoftMassMatrix(mEr, mersq * id);

  universalTrilinears(0.0);
}
*/
