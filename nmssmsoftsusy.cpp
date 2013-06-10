/** \file nmssmsoftsusy.cpp
    Project: NMSSMSOFTSUSY
    Author: Ben Allanach, Peter Athron, Lewis Tunstall, Alexander Voigt
    Manual: TBW
    Webpage:  https://github.com/Expander/softsusy.git
*/

#include "nmssmsoftsusy.h"

#ifdef NMSSMSOFTSUSY_H

extern double sw2, gnuL, guL, gdL, geL, guR, gdR, geR, yuL, yuR, ydL,
  ydR, yeL, yeR, ynuL;

 //PA: A print method used in development.  I find it useful and easier to read than couting the normal display function or calling printlong etc.
void NmssmSoftsusy::printall(){

  cout << "mH1 = " << displayPhys().mh0(1) << endl;
  cout << "mH2 = " << displayPhys().mh0(2) << endl;
  cout << "mH3 = " << displayPhys().mh0(3) << endl;

  cout << "mA1 = " << displayPhys().mA0(1) << endl;
  cout << "mA2 = " << displayPhys().mA0(2) << endl;
  cout << "mHpm =" << displayPhys().mHpm << endl;


  cout << "mch = " <<  displayPhys().mch << endl;
  cout << "mneut = " <<   displayPhys().mneut << endl;
  cout << "mGluino = " <<   displayPhys().mGluino << endl;
  cout << "mu = " <<   displayPhys().mu << endl;
  cout << "md = " <<   displayPhys().md << endl;
  cout << "me = " <<   displayPhys().me << endl;
  cout << "msnu = " <<  displayPhys().msnu << endl;

    for(int i = 1; i <= 3; i++){
      for(int j = 1; j <= 3; j++){
	cout << " Au(" << i << "," << j << ") = "  << displaySoftA(UA, i, j) << endl;
  cout << " Ad(" << i << "," << j << ") = "  << displaySoftA(DA, i, j) << endl;
  cout << " Ae(" << i << "," << j << ") = "  << displaySoftA(EA, i, j) << endl;
      }
    }
  cout << " Alam = "  <<displaySoftAlambda() << endl;
  cout << " Akap = "  << displaySoftAkappa() << endl;

  cout << " mQl = "  <<displaySoftMassSquared(mQl) << endl;
  cout << " mUr = "  <<displaySoftMassSquared(mUr) << endl;
  cout << " mDr = "  <<displaySoftMassSquared(mDr) << endl;
  cout << " mEr = "  <<displaySoftMassSquared(mEr) << endl;
  cout << " mLl = "  <<displaySoftMassSquared(mLl) << endl;

  cout << "mH1sq = " << displayMh1Squared() << endl;
  cout << "mH2sq = " << displayMh2Squared() << endl;
  cout << "m3sq = " << displayM3Squared() << endl;
  cout << "mSsq = " << displayMsSquared() << endl;
  cout << "mSpsq = " << displayMspSquared() << endl;
  cout << "xiS = " << displayXiS() << endl;

  cout << "M1 = "  << displayGaugino(1) << endl;
  cout << "M2 = "  << displayGaugino(2) << endl;
  cout << "M3 = "  << displayGaugino(3) << endl;

  cout << " lam = "  <<displayLambda() << endl;
  cout << " kap = "  <<displayKappa() << endl;
  cout << " Svev = "  <<displaySvev() << endl;
  cout <<"mupr = " <<displayMupr() << endl;
  cout << "xiF = "  <<displayXiF() << endl;

  double Beff =displaySoftAlambda() +displayKappa() *displaySvev() / (sqrt(2.0));
  double  m3hatsq =displayM3Squared() + displayLambda() * (displayMupr() *displaySvev() +displayXiF() );
  double  mueff =displayLambda() *displaySvev() / (sqrt(2.0)) ;
  cout << "mueff = "  << displayLambda() *displaySvev() / (sqrt(2.0)) << endl;
  cout << " Beff = "  << Beff << endl;
  cout << "m3hatsq = " << m3hatsq << endl;
  cout << "m3effssq = "  << mueff * Beff + m3hatsq << endl;
  cout << "normal mu = "  << displaySusyMu() << endl;

  cout << "mH1sq = " << displayMh1Squared() << endl;
  cout << "mH2sq = " << displayMh2Squared() << endl;
  cout << "m3sq = " << displayM3Squared() << endl;
  cout << "mSsq = " << displayMsSquared() << endl;
  cout << "mSpsq = " << displayMspSquared() << endl;
  cout << "xiS = " << displayXiS() << endl;
}
//PA: obtains NMSSM H1-sfermion-sfermion couplings
  //for 3rd generation sfermions
void NmssmSoftsusy::H1SfSfCouplings(DoubleMatrix & lTS1Lr, DoubleMatrix & lBS1Lr, DoubleMatrix  & lTauS1Lr, double gmzOcthW, double mu,  double cosb, double v1){
  //PA: NMSSM parameters required.
  double s = displaySvev();
  double lam = displayLambda();
  //PA: fill with parts from MSSM
   //PA: Add extra NMSSM coupling.
  //PA: minus sign since mu = - displaySusyMu for BPMZ conventions here.
  double    mueff   = mu - lam * s / root2;
  Softsusy<SoftParsNmssm>::H1SfSfCouplings(lTS1Lr, lBS1Lr, lTauS1Lr, gmzOcthW, mueff, cosb, v1);
}
//PA: obtains NMSSM H2-sfermion-sfermion couplings
  //for 3rd generation sfermions
void NmssmSoftsusy::H2SfSfCouplings(DoubleMatrix & lTS2Lr, DoubleMatrix & lBS2Lr, DoubleMatrix  & lTauS2Lr, double gmzOcthW, double mu,  double sinb){
  //PA: NMSSM parameters required.
  double s = displaySvev();
  double lam = displayLambda();
  //PA: fill with parts from MSSM
   //PA: Add extra NMSSM coupling.
  //PA: minus sign since mu = - displaySusyMu for BPMZ conventions here.
  double    mueff   = mu - lam * s / root2;
  Softsusy<SoftParsNmssm>::H2SfSfCouplings(lTS2Lr, lBS2Lr, lTauS2Lr, gmzOcthW, mueff, sinb);
}
//PA: obtains NMSSM S-sfermion-sfermion couplings
  //for 3rd generation sfermions
void NmssmSoftsusy::SSfSfCouplings(DoubleMatrix & lTS3Lr, DoubleMatrix & lBS3Lr, DoubleMatrix  & lTauS3Lr,  double lam){
  double v1 = displayHvev() * cos(atan(displayTanb()));
  double v2 = displayHvev() * sin(atan(displayTanb()));
  double ht = displayDrBarPars().ht,  hb = displayDrBarPars().hb;
  double htau = displayDrBarPars().htau;
  lTS3Lr(1, 1) = 0.0;
  lTS3Lr(1, 2) = - 0.5 * ht * lam * v1;
  lTS3Lr(2, 1) = lTS3Lr(1, 2);
  lTS3Lr(2, 2) = 0.0;

  lBS3Lr(1, 1) = 0.0;
  lBS3Lr(1, 2) = - 0.5 * hb * lam * v2;
  lBS3Lr(2, 1) = lBS3Lr(1, 2);
  lBS3Lr(2, 2) = 0.0;

  lTauS3Lr(1, 1) = 0.0;
  lTauS3Lr(1, 2) = - 0.5 * htau * lam * v2;
  lTauS3Lr(2, 1) = lTauS3Lr(1, 2);
  lTauS3Lr(2, 2) = 0.0;
}

double NmssmSoftsusy::doCalcTadSSfermions(DoubleMatrix lTS3Lr, DoubleMatrix lBS3Lr, DoubleMatrix lTauS3Lr, double q, double s){
  drBarPars fL(displayDrBarPars());
  DoubleMatrix lTS312(2, 2), lBS312(2, 2), lTauS312(2, 2), rotate(2, 2);
  rotate = rot2d(fL.thetat);
  lTS312 = rotate * lTS3Lr * rotate.transpose();
  rotate = rot2d(fL.thetab);
  lBS312 = rotate * lBS3Lr * rotate.transpose();
  rotate = rot2d(fL.thetatau);
  lTauS312 = rotate * lTauS3Lr * rotate.transpose();
  /// third generation sfermions
  double sfermions = 3.0 * lTS312(1, 1) / s * a0(fL.mu(1, 3), q);
  sfermions = sfermions + 3.0 * lTS312(2, 2) / s * a0(fL.mu(2, 3), q);
  sfermions = sfermions + 3.0 * lBS312(1, 1) / s * a0(fL.md(1, 3), q);
  sfermions = sfermions + 3.0 * lBS312(2, 2) / s * a0(fL.md(2, 3), q);
  sfermions = sfermions + lTauS312(1, 1) / s * a0(fL.me(1, 3), q);
  sfermions = sfermions + lTauS312(2, 2) / s * a0(fL.me(2, 3), q);

  return sfermions;

}

//PA: for loop corrections, helps adding Higgs corrections in a tidy way
void NmssmSoftsusy::assignHiggs(DoubleVector & higgsm, DoubleVector & higgsa,
                                  DoubleVector & higgsc) const {
  drBarPars f(displayDrBarPars());

  higgsm(1) = f.mh0(1);
  higgsm(2) = f.mh0(2);
  higgsm(3) = f.mh0(3);
  higgsa(1) = displayMzRun();
  higgsa(2) = f.mA0(1);
  higgsa(3) = f.mA0(2);
  higgsc(1) = displayMwRun();
  higgsc(2) = f.mHpm;
}

//PA: NMSSM routine to obtain Higgs loop parts of (16 \pi^2) t1/v1
//Includes goldstone bosons.
double NmssmSoftsusy::doCalcTad1Higgs(double q, double costhDRbar,
                                           double g, double tanb){
/// LCT: NMSSM parameters
  double lam = displayLambda(), lsq = sqr(lam);
  double kap = displayKappa();
  double al  = displayTrialambda();
  double s   = displaySvev();
  double cb = cos(atan(tanb)), sb = sin(atan(tanb)), sb2 = sqr(sb);
  double v1 = displayHvev() * cb, gsq = sqr(g);
  double c2b = cos(2.0 * atan(tanb));
  double costhDRbar2 = sqr(costhDRbar);
  double mupr  = displayMupr();

  /// LCT: new variables for Higgs mixing matrices
  double alphaP	= displayDrBarPars().thetaH;
  DoubleMatrix Ppr(2, 2), P(3, 3), S(3, 3), C(2, 2);
  Ppr(1, 1) = cos(alphaP);   Ppr(1, 2) = sin(alphaP);
  Ppr(2, 1) = - Ppr(1, 2);   Ppr(2, 2) = cos(alphaP);

  /// LCT: Needs transpose to match Slavich's definition
  S = displayDrBarPars().mixh0.transpose();

  P(1, 1) = - cb;           P(1, 2) = sb;             P(1, 3) = 0.0;
  P(2, 1) = sb * Ppr(1, 1); P(2, 2) = cb * Ppr(1, 1); P(2, 3) = Ppr(1, 2);
  P(3, 1) = sb * Ppr(2, 1); P(3, 2) = cb * Ppr(2, 1); P(3, 3) = Ppr(2, 2);

  C(1, 1) = - cb;  C(1, 2) = sb;
  C(2, 1) = C(1, 2); C(2, 2) = cb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Trilinear CP-even Higgs couplings to s1 in basis HdR HuR SR
  /// Have divided Slavich's expressions by v1
  DoubleMatrix sss1(3, 3);
  sss1(1, 1) = 0.125 * gsq / costhDRbar2;
  sss1(2, 2) = 1.0 / 12.0 * (2.0 * lsq - 0.5 * gsq / costhDRbar2);
  sss1(3, 3) = lam / 6.0 * (lam - kap * tanb);
  sss1(1, 2) = tanb / 12.0 * (2.0 * lsq - 0.5 * gsq / costhDRbar2);
  sss1(1, 3) = lsq * s / (6.0 * v1);
  sss1(2, 3) = - (al / root2 + lam * kap * s + lam * mupr / root2) / (6.0 * v1);
  sss1.symmetrise();

  /// LCT: Trilinear CP-odd Higgs couplings to s1 in basis HdI HuI SI
  DoubleMatrix pps1(3, 3);
  pps1(1, 1) = 0.125 * gsq / costhDRbar2;
  pps1(2, 2) = 0.25 * (2.0 * lsq - 0.5 * gsq / costhDRbar2);
  pps1(3, 3) = 0.5 * lam * (lam + kap * tanb);
  pps1(2, 3) = 0.5 * (al / root2 - lam * kap * s - lam * mupr / root2) / v1;
  pps1.symmetrise();

  /// LCT: Rotate to mass basis s1 Hi Hj
  DoubleMatrix hhs1(3, 3), aas1(3, 3);
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <=3; j++) {
        for (int a = 1; a <= 3; a++) {
           for (int b = 1; b <= 3; b++) {
              hhs1(i, j) = hhs1(i, j) + 3.0 * S(i, a) * S(j, b) * sss1(a, b);
              aas1(i, j) = aas1(i, j) + P(i, a) * P(j, b) * pps1(a, b);
           }
        }
     }
  }

  /// LCT: Trilinear with charged Higgs'. Only diagonal entries required
  /// (Amended for NMSSM) Basis (G+ G- H+ H-)
  DoubleMatrix hphps1(2, 2);
  hphps1(1, 1) = 0.25 * gsq / costhDRbar2 * c2b + lsq * sb2;
  hphps1(2, 2) = 0.25 * gsq * (2.0 - c2b / sqr(costhDRbar)) - lsq * sb2;
  double higgs = 0.0;
  //PA: add cpeven and cpodd (inc. neutral goldstone) contriobutions to higgs
 for (int i=1; i <=3; i++) {
  higgs = higgs + hhs1(i, i) * a0(higgsm(i), q) + aas1(i, i) * a0(higgsa(i), q);
  }
  /// LCT: Charged Higgs + Goldstone
  for (int i=1; i <= 2; i++) {
     higgs = higgs + hphps1(i, i) * a0(higgsc(i), q);
  }

 return higgs;
}


//PA: NMSSM routine to obtain Higgs loop parts of (16 \pi^2) t2/v2
//Includes goldstone bosons.
double NmssmSoftsusy::doCalcTad2Higgs(double q, double costhDRbar,
                                           double g, double tanb){
/// LCT: NMSSM parameters
  double lam = displayLambda(), lsq = sqr(lam);
  double kap = displayKappa();
  double al  = displayTrialambda();
  double s   = displaySvev();
  double cb = cos(atan(tanb)), sb = sin(atan(tanb)), cb2 = sqr(cb);
  double v2 = displayHvev() * sb, gsq = sqr(g);
  double c2b = cos(2.0 * atan(tanb));
  double costhDRbar2 = sqr(costhDRbar);
  /// mueff defined to be negative to remain consistent with SOFTSUSY convention
  // double mueff = mu - lam * s / root2;
  double mupr  = displayMupr();

  /// LCT: new variables for Higgs mixing matrices
  double alphaP	= displayDrBarPars().thetaH;
  DoubleMatrix Ppr(2, 2), P(3, 3), S(3, 3), C(2, 2);
  Ppr(1, 1) = cos(alphaP);   Ppr(1, 2) = sin(alphaP);
  Ppr(2, 1) = - Ppr(1, 2);   Ppr(2, 2) = cos(alphaP);

  /// LCT: Needs transpose to match Slavich's definition
  S = displayDrBarPars().mixh0.transpose();

  P(1, 1) = - cb;           P(1, 2) = sb;             P(1, 3) = 0.0;
  P(2, 1) = sb * Ppr(1, 1); P(2, 2) = cb * Ppr(1, 1); P(2, 3) = Ppr(1, 2);
  P(3, 1) = sb * Ppr(2, 1); P(3, 2) = cb * Ppr(2, 1); P(3, 3) = Ppr(2, 2);

  C(1, 1) = - cb;  C(1, 2) = sb;
  C(2, 1) = C(1, 2); C(2, 2) = cb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);
 /// LCT: Trilinear CP-even Higgs couplings to s2 in basis HdR HuR SR
  /// Have divided Slavich's expressions by v2
  DoubleMatrix sss2(3, 3);
  sss2(1, 1) = (2.0 * lsq - 0.5 * gsq / costhDRbar2) / 12.0;
  sss2(2, 2) = 0.125 * gsq / costhDRbar2;
  sss2(3, 3) = lam / 6.0 * (lam - kap / tanb);
  sss2(1, 2) =  (2.0 * lsq - 0.5 * gsq / costhDRbar2) / (12.0 * tanb);
  sss2(1, 3) = - (al / root2 + lam * kap * s + lam * mupr / root2) / (6.0 * v2);
  sss2(2, 3) = lsq * s / (6.0 * v2);
  sss2.symmetrise();

  /// LCT: Trilinear CP-odd Higgs couplings to s2 in basis HdI HuI SI
  DoubleMatrix pps2(3, 3);
  pps2(1, 1) = 0.25 * (2.0 * lsq - 0.5 * gsq / costhDRbar2);
  pps2(2, 2) = 0.125 * gsq / costhDRbar2;
  pps2(3, 3) = 0.5 * lam * (lam + kap / tanb);
  pps2(1, 3) = 0.5 * (al / root2 - lam * kap * s - lam * mupr / root2) / v2;
  pps2.symmetrise();

  /// LCT: Rotate to mass basis s2 Hi Hj
  DoubleMatrix hhs2(3, 3), aas2(3, 3);
  for (int i=1; i <= 3; i++) {
    for (int j=1; j <=3; j++) {
      for (int a = 1; a <= 3; a++) {
	for (int b = 1; b <= 3; b++) {
	  hhs2(i, j) = hhs2(i, j) + 3.0 * S(i, a) * S(j, b) * sss2(a, b);
	  aas2(i, j) = aas2(i, j) + P(i, a) * P(j, b) * pps2(a, b);
	}
      }
    }
  }

  /// LCT: Trilinear with charged Higgs'. Only need diagonal entries for tadpole
  /// (Amended for NMSSM) Basis (G+ G- H+ H-)
  DoubleMatrix hphps2(2, 2);
  hphps2(1, 1) = - 0.25 * gsq * c2b / costhDRbar2 + lsq * cb2;
  hphps2(2, 2) = 0.25 * (gsq * (2.0 + c2b / costhDRbar2)) - lsq * cb2;

  double higgs = 0.0;
  /// CP-even/-odd Higgs contributions
  for (int i=1; i <=3; i++) {
    higgs = higgs + hhs2(i, i) * a0(higgsm(i), q)
                  + aas2(i, i) * a0(higgsa(i), q);
  }
  /// Charged Higgs
  for (int i=1; i <=2; i++) {
    higgs = higgs + hphps2(i, i) * a0(higgsc(i), q);
  }

  return higgs;
}

//PA: NMSSM routine to obtain Higgs loop parts of (16 \pi^2) t1/v1
//Includes goldstone bosons.
double NmssmSoftsusy::doCalcTadSHiggs(double q, double tb){
  double lam = displayLambda(), lsq = sqr(lam);
  double kap = displayKappa(), ksq = sqr(kap);
  double s = displaySvev();
  double ak = displayTriakappa();
  double al = displayTrialambda();
  double mupr = displayMupr();
  double beta = atan(tb);
  double cb = cos(beta), sb = sin(beta);
  double cos2b = cos(2.0 * beta), sin2b = sin(2.0 * beta);
  double v1 = displayHvev() * cos(beta);
  double v2 = displayHvev() * sin(beta);

  /// LCT: new variables for Higgs mixing matrices
  double alphaP	= displayDrBarPars().thetaH;
  DoubleMatrix Ppr(2, 2), P(3, 3), S(3, 3), C(2, 2);
  Ppr(1, 1) = cos(alphaP);   Ppr(1, 2) = sin(alphaP);
  Ppr(2, 1) = - Ppr(1, 2);   Ppr(2, 2) = cos(alphaP);

  /// LCT: Needs transpose to match Slavich's definition
  S = displayDrBarPars().mixh0.transpose();

  P(1, 1) = - cb;           P(1, 2) = sb;             P(1, 3) = 0.0;
  P(2, 1) = sb * Ppr(1, 1); P(2, 2) = cb * Ppr(1, 1); P(2, 3) = Ppr(1, 2);
  P(3, 1) = sb * Ppr(2, 1); P(3, 2) = cb * Ppr(2, 1); P(3, 3) = Ppr(2, 2);

  C(1, 1) = - cb;  C(1, 2) = sb;
  C(2, 1) = C(1, 2); C(2, 2) = cb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// All of Slavich's couplings divided by svev to get TSOVS
  /// LCT: Trilinear CP-even Higgs couplings to s3 in basis HdR HuR SR
  DoubleMatrix sss3(3, 3);
  sss3(1, 1) = lsq / 6.0;
  sss3(2, 2) = lsq / 6.0;
  sss3(3, 3) = (ak / 3.0 + root2 * ksq * s + kap * mupr) / (s * root2); ///<< Bug fixed.
  sss3(1, 2) = -(al / root2 + lam * kap * s + lam * mupr / root2) / (6.0 * s);
  sss3(1, 3) = lam * (lam * v1 - kap * v2) / (6.0 * s);
  sss3(2, 3) = lam * (lam * v2 - kap * v1) / (6.0 * s);
  sss3.symmetrise();

  /// LCT: Trilinear CP-odd Higgs couplings to s3 in basis HdI HuI SI
  DoubleMatrix pps3(3, 3);
  pps3(1, 1) = 0.5 * lsq;
  pps3(2, 2) = pps3(1, 1);
  pps3(3, 3) = (-ak + root2 * ksq * s + kap * mupr) / (s * root2);
  pps3(1, 2) = 0.5 * (al / root2 + lam * kap * s + lam * mupr / root2) / s;
  pps3(1, 3) = -0.5 * lam * kap * v2 / s;
  pps3(2, 3) = -0.5 * lam * kap * v1 / s;
  pps3.symmetrise();

  /// LCT: Rotate to mass basis s3 Hi Hj, s3 Ai Aj
  DoubleMatrix hhs3(3, 3), aas3(3, 3);
  for (int i=1; i <= 3; i++) {
    for (int j=1; j <= 3; j++) {
      for (int a = 1; a <= 3; a++) {
	for (int b = 1; b <= 3; b++) {
	  hhs3(i, j) = hhs3(i, j) + 3.0 * S(i, a) * S(j, b) * sss3(a, b);
	  aas3(i, j) = aas3(i, j) + P(i, a) * P(j, b) * pps3(a, b);
	}
      }
    }
  }

  /// LCT: Trilinear with charged Higgs. Basis (G+ G- H+ H-)
  DoubleMatrix hphps3(2, 2);
  hphps3(1, 1) = 0.5 * (2.0 * lsq * s - (root2 * al + 2.0 * lam * kap * s
					 + root2 * lam * mupr) * sin2b) / s;
  hphps3(2, 2) = 0.5 * (2.0 * lsq * s + (root2 * al + 2.0 * lam * kap * s
					 + root2 * lam * mupr) * sin2b) / s;

  double higgs = 0.0;
  /// CP-even/-odd Higgs
  for (int i=1; i <=3; i++) {
    higgs = higgs + hhs3(i, i) * a0(higgsm(i), q) + aas3(i, i) * a0(higgsa(i), q);
  }
  /// Charged Higgs
  for (int i=1; i <= 2; i++) {
    higgs = higgs + hphps3(i, i) * a0(higgsc(i), q);
  }
  return higgs;

}

double NmssmSoftsusy::doCalcTad1Neutralinos(double q, double costhDRbar,
                                           double g, double cosb){
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  double neutralinos = 0.0;
  double lam = displayLambda();
  double tanthDRbar = tan(acos(costhDRbar));
  for (int family = 1; family <= 5; family++){
    neutralinos = neutralinos -
      sqr(g) * mneut(family) / (displayMwRun() * cosb) *
       (n(family, 3) * (n(family, 2) - n(family, 1) * tanthDRbar)).real() *
       a0(mneut(family), q)
      + 2 * root2 * lam  * mneut(family) / (displayHvev() * cosb) *
       (n(family, 4) * n(family, 5)).real() * a0(mneut(family), q);
  }
    return neutralinos;
}
double NmssmSoftsusy::doCalcTad2Neutralinos(double q, double costhDRbar,
                                           double g, double sinb){
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  double lam = displayLambda();
  double tanthDRbar = tan(acos(costhDRbar));
  double neutralinos = 0.0;
  for (int family = 1; family <= 5; family++)
    neutralinos = neutralinos + sqr(g) * mneut(family) /
      (displayMwRun() * sinb) *
      (n(family, 4) * (n(family, 2) - n(family, 1) * tanthDRbar)).real() *
      a0(mneut(family), q)
      + 2.0 * root2 * lam  * mneut(family) / (displayHvev() * sinb) *
      (n(family, 3) * n(family, 5)).real() * a0(mneut(family), q); ///<< Extra NMSSM piece
  return neutralinos;
}
double NmssmSoftsusy::doCalcTadSNeutralinos(double q,  double lam, double kap){
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  double s = displaySvev();
  double neutralinos = 0.0;
  for (int family = 1; family <= 5; family++)
    neutralinos = neutralinos
    - 2.0 * root2 * kap * mneut(family) / s
    * (n(family, 5) * n(family, 5)).real() * a0(mneut(family), q)
		+ 2.0 * root2 * lam  * mneut(family) / s *
    (n(family, 3) * n(family, 4)).real() * a0(mneut(family), q);
 return neutralinos;
}
double NmssmSoftsusy::doCalcTadSCharginos(double q,  double lam){
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz);
  DoubleVector mch(displayDrBarPars().mchBpmz);
  double s = displaySvev();
  double charginos = 0.0;
  for (int family=1; family<=2; family++)
    charginos = charginos - root2 * 2.0 * lam / s
		* mch(family) * (v(family, 2) * u(family, 2)).real()
		* a0(mch(family), q);
  return charginos;
}

double NmssmSoftsusy::doCalcTadpole1oneLoop(double mt, double sinthDRbar) {
 if (displayDrBarPars().mu(1, 3) == 0.0 || displayDrBarPars().mu(2, 3) == 0.0) {
   if (PRINTOUT > 1)
    cout << "Trying to calculate tadpole without having first calculated"
	 << " the DRbar masses.\n";
   return 0.0;
  }
 double g = displayGaugeCoupling(2), costhDRbar = cos(asin(sinthDRbar)),
     tanb = displayTanb(), cosb = cos(atan(tanb)),
    mu = -displaySusyMu(), q = displayMu(),
    mz = displayMzRun();
  double beta = atan(displayTanb());
  double v1 = displayHvev() * cos(beta);

  const double gmzOcthW =  g * mz / costhDRbar;
  double costhDRbar2 = sqr(costhDRbar);
  double fermions = Softsusy<SoftParsNmssm>::doCalcTad1fermions(q, v1);
  /// PA: stop, sbottom, stau, couplings in the left right basis
  // will be stored in these matrices
  DoubleMatrix lTS1Lr(2, 2), lBS1Lr(2, 2), lTauS1Lr(2, 2);
  H1SfSfCouplings(lTS1Lr, lBS1Lr, lTauS1Lr, gmzOcthW, mu, cosb, v1);

  //PA: Now we take these couplings and obtain sfermion contributions
  double sfermions =  Softsusy<SoftParsNmssm>::doCalcTad1Sfermions(lTS1Lr, lBS1Lr, lTauS1Lr,
					  costhDRbar);
  //PA: Higgs contributions, including goldstone bosons
  double higgs = doCalcTad1Higgs(q, costhDRbar, g, tanb);
  /// Neutralinos
  double neutralinos = doCalcTad1Neutralinos(q, costhDRbar, g, cosb);
  /// Charginos PA::Unchanged from MSSM
  double charginos = doCalcTad1Charginos(q, g, cosb);
  //PA:  gauge boson contributions are not changed from the MSSM,
  //     *but* the goldstone parts are now included in the Higgs
  //      contributions so this is *not* the same as the Softsusy version
   double gaugeBosons = 3.0 * sqr(g) / 4.0 *
    (2.0 * a0(displayMwRun(), q) + a0(mz, q) / costhDRbar2);

   double delta = fermions + sfermions + higgs + charginos + neutralinos +
      gaugeBosons;

  return delta / (16.0 * sqr(PI));

}

double NmssmSoftsusy::doCalcTadpole2oneLoop(double mt, double sinthDRbar) {
 if (displayDrBarPars().mu(1, 3) == 0.0 || displayDrBarPars().mu(2, 3) == 0.0) {
   if (PRINTOUT > 1)
    cout << "Trying to calculate tadpole without having first calculated"
	 << " the DRbar masses.\n";
   return 0.0;
  }
  double g = displayGaugeCoupling(2),
    costhDRbar = cos(asin(sinthDRbar)), costhDRbar2 = sqr(costhDRbar),
    tanb = displayTanb(), sinb = sin(atan(tanb)),
    beta = atan(displayTanb()), mu = -displaySusyMu(), q = displayMu(),
    mz = displayMzRun();
  const double gmzOcthW =  g * mz / costhDRbar;
  double fermions = Softsusy<SoftParsNmssm>::doCalcTad2fermions(q);
  /// Sfermion couplings
  DoubleMatrix lTS2Lr(2, 2),  lBS2Lr(2, 2),  lTauS2Lr(2, 2);
  H2SfSfCouplings(lTS2Lr, lBS2Lr, lTauS2Lr, gmzOcthW, mu, sinb);
  //PA: Now we take these couplings and obtain sfermion contributions
  double sfermions =  Softsusy<SoftParsNmssm>::doCalcTad2Sfermions(lTS2Lr, lBS2Lr, lTauS2Lr,
						    costhDRbar);
  //PA: Higgs contributions, including goldstone bosons
  double  higgs = doCalcTad2Higgs(q, costhDRbar, g, tanb);
  double neutralinos = doCalcTad2Neutralinos(q, costhDRbar, g, sinb);
  /// Charginos PA::Unchanged from MSSM
  double charginos = doCalcTad2Charginos(q, g, sinb);
  //PA:  gauge boson contributions are not changed from the MSSM,
  //     *but* the goldstone parts are now included in the Higgs
  //      contributions so this is *not* the same as the Softsusy version
  double gaugeBosons = 3.0 * sqr(g) / 4.0 *
    (2.0 * a0(displayMwRun(), q) + a0(mz, q) / costhDRbar2);

  double delta = fermions + sfermions + higgs + charginos + neutralinos +
    gaugeBosons;

  return delta / (16.0 * sqr(PI));

}

double NmssmSoftsusy::doCalcTadpoleSoneLoop(double mt, double sinthDRbar) {
 if (displayDrBarPars().mu(1, 3) == 0.0 || displayDrBarPars().mu(2, 3) == 0.0) {
   if (PRINTOUT > 1)
    cout << "Trying to calculate tadpole without having first calculated"
	 << " the DRbar masses.\n";
   return 0.0;
  }
 /// PA: NMSSM parameters needed
  double q = displayMu();
  double lam = displayLambda();
  double s = displaySvev();
  double kap = displayKappa();
  double tanb = displayTanb();
  /// LCT: No fermion contributions.
  /// Sfermion couplings
  DoubleMatrix lTS3Lr(2, 2),lBS3Lr(2, 2), lTauS3Lr(2, 2);
  SSfSfCouplings(lTS3Lr, lBS3Lr, lTauS3Lr, lam);
  ///PA: 3rd gen sfermions only.
  double sfermions = doCalcTadSSfermions(lTS3Lr, lBS3Lr, lTauS3Lr, q, s);
  double higgs = doCalcTadSHiggs(q, tanb);
   /// Neutralinos
  double neutralinos = doCalcTadSNeutralinos(q, lam, kap);
  /// Charginos
  double charginos = doCalcTadSCharginos(q, lam);

 double delta = sfermions + higgs + charginos + neutralinos;
 return delta / (16.0 * sqr(PI));
}

void NmssmSoftsusy::calcTadpole1Ms1loop(double mt, double sinthDRbar) {
   double t1OV1 = doCalcTadpole1oneLoop(mt, sinthDRbar);
  if (testNan(t1OV1)) {
    flagNoMuConvergence(true);
    t1OV1 = 0.0;
  }
   setT1OV1Ms(t1OV1);
}

void NmssmSoftsusy::calcTadpole2Ms1loop(double mt, double sinthDRbar) {
   double t2OV2 = doCalcTadpole2oneLoop(mt, sinthDRbar);
  if (testNan(t2OV2)) {
    flagNoMuConvergence(true);
    t2OV2 = 0.0;
  }
  setT2OV2Ms(t2OV2);
}

void NmssmSoftsusy::calcTadpoleSMs1loop(double mt, double sinthDRbar) {
  tSOVSMs1loop = doCalcTadpoleSoneLoop(mt, sinthDRbar);
  if (testNan(tSOVSMs1loop)) {
    flagNoMuConvergence(true);
    tSOVSMs1loop = 0.0;
  }
}

void NmssmSoftsusy::treeUpSquark(DoubleMatrix & mass, double mtrun,
				double pizztMS, double sinthDRbarMS,
				int family) {
  //PA: only modification is to add lambda * s / root to mu
  double lam = displayLambda(), svev = displaySvev(), tanb = displayTanb();
  Softsusy<SoftParsNmssm>::treeUpSquark(mass, mtrun, pizztMS, sinthDRbarMS, family);
  if (family == 3){
     mass(1, 2) = mass(1, 2) -  mtrun * lam * svev / (root2 * tanb);
     mass(2, 1) = mass(1, 2);
  }

}



void NmssmSoftsusy::treeDownSquark(DoubleMatrix & mass, double mbrun,
				double pizztMS, double sinthDRbarMS,
				int family) {
  //PA: only modification is to add lambda * s / root to mu
  double lam = displayLambda(), svev = displaySvev(), tanb = displayTanb();
  Softsusy<SoftParsNmssm>::treeDownSquark(mass, mbrun, pizztMS, sinthDRbarMS, family);
  if (family == 3){
     mass(1, 2) = mass(1, 2) -  mbrun * lam * svev * tanb / (root2);
     mass(2, 1) = mass(1, 2);
   }
}


void NmssmSoftsusy::treeChargedSlepton(DoubleMatrix & mass, double mtaurun,
				double pizztMS, double sinthDRbarMS,
				int family) {
   //PA: only modification is to add lambda * s / root to mu
  double lam = displayLambda(), svev = displaySvev(), tanb = displayTanb();
  Softsusy<SoftParsNmssm>::treeChargedSlepton(mass, mtaurun, pizztMS, sinthDRbarMS, family);
  if (family == 3) {
     mass(1, 2) = mass(1, 2) -  mtaurun * lam * svev * tanb / (root2);
     mass(2, 1) = mass(1, 2);
  }
}

/// LCT: new routine to allocate NMSSM chargino masses
void NmssmSoftsusy::calcDrBarCharginos(DoubleMatrix & mass, double beta, double mw) {
  double lam = displayLambda(), svev = displaySvev();

  Softsusy<SoftParsNmssm>::calcDrBarCharginos(mass, beta, mw);
  mass(2, 2) = mass(2, 2) + lam * svev / root2;
}

/// LCT: new routine for NMSSM neutralino masses
void NmssmSoftsusy::calcDrBarNeutralinos(DoubleMatrix & mass, double beta, double mz, double mw, double sinthDRbar) {
  double lam = displayLambda(), kap = displayKappa();
  double smu = displaySusyMu(), mupr = displayMupr();
  double cosb = cos(beta), sinb = sin(beta);
  double vev = displayHvev(), svev = displaySvev();

  /// Call MSSM 4 x4 neutralino mass matrix
  Softsusy<SoftParsNmssm>::calcDrBarNeutralinos(mass, beta, mz, mw, sinthDRbar);

  /// Fill remaining values
  mass(3, 4) = mass(3, 4) - lam * svev / root2;
  mass(3, 5) = - lam * vev * sinb / root2;
  mass(4, 5) = - lam * vev * cosb / root2;
  mass(5, 5) = root2 * kap * svev + mupr;

  /// symmetrise tree-level
  mass.symmetrise();
}


void NmssmSoftsusy::calcDrBarGauginos(double beta, double mw, double mz, double sinth, drBarPars & eg) {
 DoubleMatrix mCh(2, 2);
  calcDrBarCharginos(mCh, beta, mw);
  eg.mch = mCh.asy2by2(eg.thetaL, eg.thetaR);
  eg.mpzCharginos();

  DoubleMatrix mNeut(5, 5);
  calcDrBarNeutralinos(mNeut, beta, mz, mw, sinth);
  if (mNeut.diagonaliseSym(eg.mixNeut, eg.mneut) > TOLERANCE *
       1.0e-3) {
      ostringstream ii;
      ii << "accuracy bad in neutralino diagonalisation"<< flush;
      throw ii.str();
   }

  eg.mpzNeutralinos();
}
/// LCT: new routine to set tree-level NMSSM CP-even and odd Higgs mass matrices squared
void NmssmSoftsusy::calcDrBarHiggs(double beta, double mz2, double mw2, double sinthDRbar, drBarPars & eg) {
  double tanb = displayTanb();
  double sinb = sin(beta), cosb = cos(beta);
  double sinb2 = sqr(sinb), cosb2 = sqr(cosb);

	/// LCT: NMSSM parameters
  double lam = displayLambda(), kap = displayKappa();
  double mupr = displayMupr(), smu = displaySusyMu();
  double al = displayTrialambda(), ak = displayTriakappa();
  double xiF = displayXiF(), xiS = displayXiS();
  double m3sq = displayM3Squared(), mSpsq = displayMspSquared();

  double vev = displayHvev(), v1 = vev * cosb, v2 = vev * sinb;
  double svev = displaySvev();

  double mueff = lam * svev / root2 + smu;
  double m3effsq = m3sq + lam * (mupr * svev / root2 + xiF);
  DoubleMatrix mS(3, 3), mPpr(3, 3), mP(3, 3), mP2(2, 2);

  /// CP-even Higgs in EHT notation and basis (HdR, HuR, SR)
  mS(1, 1) = mz2 * cosb2 + (0.5 * lam * kap * sqr(svev)
                            + al * svev / root2 + m3effsq) * tanb;
  mS(1, 2) = - sinb * cosb * mz2 + sqr(lam) * sqr(vev) * sinb * cosb
     - (0.5 * lam * kap * sqr(svev) + al * svev / root2 + m3effsq);
  mS(1, 3) = root2 * lam * mueff * v1
     - (al + lam * kap * svev * root2 + lam * mupr) * v2 / root2;
  mS(2, 2) = mz2 * sinb2 + (0.5 * lam * kap * sqr(svev)
                            + al * svev / root2 + m3effsq) / tanb;
  mS(2, 3) = root2 * lam * mueff * v2
     - (al + lam * kap * svev * root2 + lam * mupr) * v1 / root2;
  mS(3, 3) = (al + lam * mupr) * v1 * v2 / (root2 * svev)
     + svev / root2 * (ak + 2.0 * root2 * sqr(kap) * svev
                       + 3.0 * kap * mupr) - root2 * (xiS + xiF * mupr) / svev;

  mS.symmetrise();

  /// CP-odd Higgs in EHT notation and basis (HdI, HuI, SI)
  mPpr(1, 1) = (0.5 * lam * kap * sqr(svev) + al * svev / root2
                + m3effsq) * tanb;
  mPpr(1, 2) = 0.5 * lam * kap * sqr(svev) + al * svev / root2 + m3effsq;
  mPpr(1, 3) = v2 / root2 * (al - root2 * lam * kap * svev - lam * mupr);
  mPpr(2, 2) = (0.5 * lam * kap * sqr(svev) + al * svev / root2
                + m3effsq) / tanb;
  mPpr(2, 3) = v1 / root2 * (al - root2 * lam * kap * svev - lam * mupr);
  mPpr(3, 3) = (al + lam * kap * svev /root2
             + 3.0 * lam * kap * svev / root2
             + lam * mupr) * v1 * v2 / (root2 * svev) - 3.0 * ak * svev / root2
     - 2.0 * mSpsq - kap * mupr * svev / root2 - xiF * (4.0 * kap
                + root2 * mupr / svev) - root2 * xiS / svev;

  mPpr.symmetrise();

	/// LCT: Diagonalise 3 x 3 CP-even Higgs mass matrix.
	/// Mass basis (H1, H2, H3) is ordered in increasing mass mH1 < mH2 < mH3
  DoubleVector mhsq(3);
  if (mS.diagonaliseSym(eg.mixh0, mhsq) > TOLERANCE *
      1.0e-3) { /// LCT: default 1.0e-3
     ostringstream ii;
    ii << "accuracy bad in CP-even Higgs diagonalisation"<< flush;
    throw ii.str();
  }

  if (mhsq(1) < 0. || mhsq(2) < 0. || mhsq(3) < 0.) {
     flagTachyon(h0);
     if (PRINTOUT > 1) cout << " mH1/mH2/mH3 tachyon ";
  }

  DoubleVector mHiggs(mhsq.apply(ccbSqrt));
  eg.mh0 = mHiggs;

  // LCT: Rotate CP-odd mass^2 matrix into (G, A, S_I) basis
  mP = rot3d(beta).transpose() * mPpr * rot3d(beta);
  /// LCT: Drop Goldstone from 3 x 3 CP-odd Higgs mass^2 matrix and
  /// construct 2 x 2 matrix in (A, S_I) basis
  mP2(1, 1) = mP(2, 2);
  mP2(1, 2) = mP(2, 3);
  mP2(2, 1) = mP(3, 2);
  mP2(2, 2) = mP(3, 3);


  /// LCT: Diagonalise
  //PA using thetaH for now since in nmssm this is reundent
  DoubleVector mSq = mP2.sym2by2(eg.thetaH);
  //
  if (mSq(1) < 0. || mSq(2) < 0.) {
     flagTachyon(A0);
     if (PRINTOUT > 1) cout << " mA1/mA2 tachyon";
  }
  DoubleVector temp(mSq.apply(ccbSqrt));
  if (temp(1) > temp(2)) eg.thetaH = eg.thetaH + PI * 0.5;

  int pos;
  eg.mA0(1) = temp.min(pos); eg.mA0(2) = temp.max();

  /// LCT: Have not (as yet) included 1-loop corrections to mA1/mA2
  /// in definition of mHpm
  double mHpmsq = mP2(1, 1) + mw2 - 0.5 * sqr(vev) * sqr(lam);
  if (mHpmsq < 0.) {
     flagTachyon(h0); mHpmsq = fabs(mHpmsq);
  }
  eg.mHpm = sqrt(mHpmsq);

}




void NmssmSoftsusy::calcDrBarPars() {
  drBarPars eg(displayDrBarPars());
  //PA: set up drBarPars object to be suitable for NMSSM
  eg.mh0.setEnd(3);
  eg.mA0.setEnd(2);
  eg.mixh0.resize(3,3);
  eg.mixA0.resize(2,2);
  eg.mneut.setEnd(5);
  eg.mixNeut.resize(5,5);

  eg.mnBpmz.setEnd(5);
  eg.nBpmz.resize(5,5);


  /// First, must define mstop,sbot,stau and mixing angles in DRbar scheme
  double beta = atan(displayTanb()), mzPole = displayMz();
  double sinthDRbar = calcSinthdrbar();
  double mz = displayMzRun(), mz2 = sqr(mz);
  double pizzt = sqr(mz) - sqr(mzPole);
  double vev = displayHvev();
  Softsusy<SoftParsNmssm>::setNeutCurrCouplings(sinthDRbar, sw2, guL, gdL, geL, guR, gdR, geR);
  Softsusy<SoftParsNmssm>::calcDRTrilinears(eg, vev, beta);
  eg.mGluino = displayGaugino(3);
  DoubleVector mSq(2);
  int family; for(family = 1; family <= 3; family++) {

     DoubleMatrix mSquared(2, 2);
     treeUpSquark(mSquared, eg.mt, pizzt, sinthDRbar, family);
     mSq = mSquared.sym2by2(eg.thetat);
     if (mSq(1) < 0. || mSq(2) < 0.) {
        switch(family) {
        case 1: flagTachyon(sup); break;
        case 2: flagTachyon(scharm); break;
        case 3: flagTachyon(stop); break;
        default: throw("Bad family number in calcDrBarPars\n");
        }
        if (PRINTOUT > 2) cout << " tree sup(" << family << ") tachyon ";
     }
    DoubleVector mstopDRbar(mSq.apply(zeroSqrt));

    treeDownSquark(mSquared, eg.mb, pizzt, sinthDRbar, family);
    mSq = mSquared.sym2by2(eg.thetab);
    if (mSq(1) < 0. || mSq(2) < 0.) {
      switch(family) {
      case 1: flagTachyon(sdown); break;
      case 2: flagTachyon(sstrange); break;
      case 3: flagTachyon(sbottom); break;
      default: throw("Bad family number in calcDrBarPars\n");
      }
    if (PRINTOUT > 1) cout << " tree sdown(" << family << ") tachyon ";
    }
    DoubleVector msbotDRbar(mSq.apply(zeroSqrt));

    treeChargedSlepton(mSquared, eg.mtau, pizzt, sinthDRbar, family);
    mSq = mSquared.sym2by2(eg.thetatau);
    if (mSq(1) < 0. || mSq(2) < 0.) {
      switch(family) {
      case 1: flagTachyon(selectron); break;
      case 2: flagTachyon(smuon); break;
      case 3: flagTachyon(stau); break;
      default: throw("Bad family number in calcDrBarPars\n");
      }
    if (PRINTOUT > 1) cout << " tree selectron(" << family << ") tachyon ";
    }
    DoubleVector mstauDRbar(mSq.apply(zeroSqrt));

    int i; for (i=1; i<=2; i++) {
      eg.mu(i, family) = mstopDRbar(i);    eg.md(i, family) = msbotDRbar(i);
      eg.me(i, family) = mstauDRbar(i);
    }
    double mSnuSquared;
    //PA unmodified by NMSSM.
    treeSnu(mSnuSquared, pizzt, family);
    if (mSnuSquared < 0.) {
      switch(family) {
      case 1: flagTachyon(snue); break;
      case 2: flagTachyon(snumu); break;
      case 3: flagTachyon(snutau); break;
      default: throw("Bad family number in calcDrBarPars\n");
      }
    if (PRINTOUT > 1) cout << " tree sneutrino(" << family << ") tachyon@"
			   << displayMu() << " ";
    }
    eg.msnu(family) = zeroSqrt(mSnuSquared);
  }

  double mw = displayMwRun();
  double mw2 = sqr(mw);
  eg.mw = mw;
  eg.mz = mz;
  calcDrBarGauginos(beta, mw, mz, sinthDRbar, eg);
  calcDrBarHiggs(beta, mz2, mw2, sinthDRbar, eg);

  setDrBarPars(eg);

  return;

}

#endif
