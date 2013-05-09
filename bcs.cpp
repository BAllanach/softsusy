
#include "bcs.h"
#include "softsusy.h"

void generalBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  MssmSusy s; SoftParsMssm r;
  double m3sq = m.displayM3Squared();
  s = m.displaySusy();
  r.set(inputParameters);
  r.setM3Squared(m3sq);
  m.setSoftPars(r);
  m.setSusy(s);

  return;
}

/// This one doesn't overwrite mh1sq or mh2sq at the high scale
void generalBcs2(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  MssmSusy s; SoftParsMssm r;
  double mh1sq = m.displayMh1Squared(); 
  double mh2sq = m.displayMh2Squared();
  double m3sq = m.displayM3Squared();
  s = m.displaySusy();
  r.set(inputParameters);
  r.setMh1Squared(mh1sq);
  r.setMh2Squared(mh2sq);
  r.setM3Squared(m3sq);
  m.setSoftPars(r);
  m.setSusy(s);

  return;
}

void extendedSugraBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  int i;
  for (i=1; i<=3; i++) m.setGauginoMass(i, inputParameters.display(i));
  if (inputParameters.display(25) > 1. && m.displaySetTbAtMX()) 
    m.setTanb(inputParameters.display(25));
  m.setTrilinearElement(UA, 1, 1, m.displayYukawaElement(YU, 1, 1) * 
			inputParameters.display(11));
  m.setTrilinearElement(UA, 2, 2, m.displayYukawaElement(YU, 2, 2) * 
			inputParameters.display(11));
  m.setTrilinearElement(UA, 3, 3, m.displayYukawaElement(YU, 3, 3) * 
			inputParameters.display(11));
  m.setTrilinearElement(DA, 1, 1, m.displayYukawaElement(YD, 1, 1) * 
			inputParameters.display(12));
  m.setTrilinearElement(DA, 2, 2, m.displayYukawaElement(YD, 2, 2) * 
			inputParameters.display(12));
  m.setTrilinearElement(DA, 3, 3, m.displayYukawaElement(YD, 3, 3) * 
			inputParameters.display(12));
  m.setTrilinearElement(EA, 1, 1, m.displayYukawaElement(YE, 1, 1) * 
			inputParameters.display(13));
  m.setTrilinearElement(EA, 2, 2, m.displayYukawaElement(YE, 2, 2) * 
			inputParameters.display(13));
  m.setTrilinearElement(EA, 3, 3, m.displayYukawaElement(YE, 3, 3) * 
			inputParameters.display(13));
  m.setSoftMassElement(mLl, 1, 1, signedSqr(inputParameters.display(31)));
  m.setSoftMassElement(mLl, 2, 2, signedSqr(inputParameters.display(32)));
  m.setSoftMassElement(mLl, 3, 3, signedSqr(inputParameters.display(33)));
  m.setSoftMassElement(mEr, 1, 1, signedSqr(inputParameters.display(34)));
  m.setSoftMassElement(mEr, 2, 2, signedSqr(inputParameters.display(35)));
  m.setSoftMassElement(mEr, 3, 3, signedSqr(inputParameters.display(36)));
  m.setSoftMassElement(mQl, 1, 1, signedSqr(inputParameters.display(41)));
  m.setSoftMassElement(mQl, 2, 2, signedSqr(inputParameters.display(42)));
  m.setSoftMassElement(mQl, 3, 3, signedSqr(inputParameters.display(43)));
  m.setSoftMassElement(mUr, 1, 1, signedSqr(inputParameters.display(44)));
  m.setSoftMassElement(mUr, 2, 2, signedSqr(inputParameters.display(45)));
  m.setSoftMassElement(mUr, 3, 3, signedSqr(inputParameters.display(46)));
  m.setSoftMassElement(mDr, 1, 1, signedSqr(inputParameters.display(47)));
  m.setSoftMassElement(mDr, 2, 2, signedSqr(inputParameters.display(48)));
  m.setSoftMassElement(mDr, 3, 3, signedSqr(inputParameters.display(49)));

  if (!m.displayAltEwsb()) {
    m.setMh1Squared(inputParameters.display(21));
    m.setMh2Squared(inputParameters.display(22));
  }
}

/// universal mSUGRA boundary conditions
void sugraBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);

  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be a0
  ///  if (m0 < 0.0) m.flagTachyon(true); Deleted on request from A Pukhov
  m.standardSugra(m0, m12, a0);
    
  return;
}

void nuhmI(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double mH  = inputParameters.display(3);
  double a0 = inputParameters.display(4);

  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be a0
  ///  if (m0 < 0.0) m.flagTachyon(true); Deleted on request from A Pukhov
  m.standardSugra(m0, m12, a0);
  m.setMh1Squared(mH * mH); m.setMh2Squared(mH * mH);
    
  return;
}

void nuhmII(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double mH1  = inputParameters.display(3);
  double mH2  = inputParameters.display(4);
  double a0 = inputParameters.display(5);

  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be a0
  ///  if (m0 < 0.0) m.flagTachyon(true); Deleted on request from A Pukhov
  m.standardSugra(m0, m12, a0);
  m.setMh1Squared(mH1 * mH1); m.setMh2Squared(mH2 * mH2);
    
  return;
}

/// Other types of boundary condition
void amsbBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m32 = inputParameters.display(1);
  double m0 = inputParameters.display(2);

  m.standardSugra(m0, 0., 0.);
  m.addAmsb(m32);
  return;
}

void lvsBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0  = inputParameters.display(1);
  double m12 = inputParameters.display(1) * sqrt(3.);
  double a0  = -inputParameters.display(1) * sqrt(3.);

  m.standardSugra(m0, m12, a0);

  return;
}

void gmsbBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  int n5     = int(inputParameters.display(1));
  double mMess  = inputParameters.display(2);
  double lambda = inputParameters.display(3);
  double cgrav = inputParameters.display(4);

  m.minimalGmsb(n5, lambda, mMess, cgrav);
    
  return;
}

void userDefinedBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  m.methodBoundaryCondition(inputParameters);
  sugraBcs(m, inputParameters);
}

void nonUniGauginos(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);

  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be a0
  ///  if (m0 < 0.0) m.flagTachyon(true); Deleted on request from A Pukhov
  m.standardSugra(m0, m12, a0);
    
  m.setGauginoMass(2, inputParameters.display(4));
  m.setGauginoMass(3, inputParameters.display(5));

  return;
}

// Boundary conditions of split gauge mediated SUSY breaking (see
// http://www.physics.rutgers.edu/~somalwar/conlsp/slepton-coNLSP.pdf 
// for example). Note that here, mu is set at mMess instead of at the
// electroweak scale.
void splitGmsb(MssmSoftsusy & m, const DoubleVector & inputParameters) {

  double n5 = inputParameters(1);
  double lambdaL = inputParameters(2);
  double lambdaD = inputParameters(3); 
  double mMess = inputParameters(4);
  double muOm2 = inputParameters(5);
  double mAOm2 = inputParameters(6);
  double cgrav = inputParameters(7);

  double lambda1 = n5 * (0.6 * lambdaL + 0.4 * lambdaD);
  double lambda2 = n5 * lambdaL;
  double lambda3 = n5 * lambdaD;

  double m1, m2, m3;
  m1 = sqr(m.displayGaugeCoupling(1)) / (16.0 * sqr(PI)) * lambda1; 
  m2 = sqr(m.displayGaugeCoupling(2)) / (16.0 * sqr(PI)) * lambda2; 
  m3 = sqr(m.displayGaugeCoupling(3)) / (16.0 * sqr(PI)) * lambda3; 
  m.setGauginoMass(1, m1);   
  m.setGauginoMass(2, m2);   
  m.setGauginoMass(3, m3);

  m.setM32(2.37e-19 * sqrt((sqr(lambdaL) + sqr(lambdaD)) * 0.5) * 
	   mMess * cgrav);

  m.setM32(2.37e-19 * sqrt((sqr(lambdaL) + sqr(lambdaD)) * 0.5) * 
	   mMess * cgrav);

  double g1f = sqr(sqr(m.displayGaugeCoupling(1)));
  double g2f = sqr(sqr(m.displayGaugeCoupling(2)));
  double g3f = sqr(sqr(m.displayGaugeCoupling(3)));

  double lambdaP1sq = n5 * (0.6 * sqr(lambdaL) + 0.4 * sqr(lambdaD));
  double lambdaP2sq = n5 * sqr(lambdaL);
  double lambdaP3sq = n5 * sqr(lambdaD);

  double mursq, mdrsq, mersq, mqlsq, mllsq;
  mursq = 2.0 * 
    (4.0 / 3.0 * g3f * lambdaP3sq + 0.6 * 4.0 / 9.0 * g1f * lambdaP1sq) 
    / sqr(16.0 * sqr(PI));
  mdrsq = 2.0 * 
    (4.0 / 3.0 * g3f * lambdaP3sq + 0.6 * 1.0 / 9.0 * g1f * lambdaP1sq) 
    / sqr(16.0 * sqr(PI));
  mersq = 2.0 * 
    (0.6 * g1f * lambdaP1sq) 
    / sqr(16.0 * sqr(PI));
  mqlsq = 2.0 * 
    (4.0 / 3.0 * g3f * lambdaP3sq + 0.75 * g2f * lambdaP2sq + 
     0.6 * g1f / 36.0 * lambdaP1sq) 
    / sqr(16.0 * sqr(PI));
  mllsq = 2.0 * 
    (0.75 * g2f * lambdaP2sq + 0.6 * 0.25 * g1f * lambdaP1sq) 
    / sqr(16.0 * sqr(PI));

  // You need Higgs masses too!

  DoubleMatrix id(3, 3);
  id(1, 1) = 1.0; id(2, 2) = 1.0; id(3, 3) = 1.0;

  m.setSoftMassMatrix(mQl, mqlsq * id);
  m.setSoftMassMatrix(mUr, mursq * id);
  m.setSoftMassMatrix(mDr, mdrsq * id);
  m.setSoftMassMatrix(mLl, mllsq * id);  
  m.setSoftMassMatrix(mEr, mersq * id);

  m.universalTrilinears(0.0);
  DoubleVector pars(2); ///< encodes EWSB BC
  pars(1) = muOm2 * m2; 
  pars(2) = mAOm2 * m2;

  /// Save the two parameters
  m.setEwsbConditions(pars);
}
