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

void NmssmSoftsusy::H1SfSfCouplings(DoubleMatrix & lTS1Lr, DoubleMatrix & lBS1Lr, DoubleMatrix  & lTauS1Lr, double gmzOcthW, double mu,  double cosb, double v1) const {
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

void NmssmSoftsusy::H2SfSfCouplings(DoubleMatrix & lTS2Lr, DoubleMatrix & lBS2Lr, DoubleMatrix  & lTauS2Lr, double gmzOcthW, double mu,  double sinb) const
{
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
void NmssmSoftsusy::SSfSfCouplings(DoubleMatrix & lTS3Lr, DoubleMatrix & lBS3Lr, DoubleMatrix  & lTauS3Lr,  double lam) const
{
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

double NmssmSoftsusy::doCalcTadSSfermions(DoubleMatrix lTS3Lr, DoubleMatrix lBS3Lr, DoubleMatrix lTauS3Lr, double q, double s) const
{
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
                                           double g, double tanb) const
{
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

  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  /// LCT: Needs transpose to match Slavich's definition
  S = displayDrBarPars().mixh0;
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
                                           double g, double tanb) const
{
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
  DoubleMatrix Ppr(2, 2), P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  /// LCT: Needs transpose to match Slavich's definition
  S = displayDrBarPars().mixh0;

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
double NmssmSoftsusy::doCalcTadSHiggs(double q, double tb) const
{
  double lam = displayLambda(), lsq = sqr(lam);
  double kap = displayKappa(), ksq = sqr(kap);
  double s = displaySvev();
  double ak = displayTriakappa();
  double al = displayTrialambda();
  double mupr = displayMupr();
  double beta = atan(tb);
  double cb = cos(beta), sb = sin(beta);
  double sin2b = sin(2.0 * beta);
  double v1 = displayHvev() * cos(beta);
  double v2 = displayHvev() * sin(beta);

  /// LCT: new variables for Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  /// LCT: Needs transpose to match Slavich's definition
  S = displayDrBarPars().mixh0;
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
                                           double g, double cosb) const
{
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
                                           double g, double sinb) const
{
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

double NmssmSoftsusy::doCalcTadSNeutralinos(double q,  double lam, double kap) const
{
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

double NmssmSoftsusy::doCalcTadSCharginos(double q,  double lam) const
{
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

double NmssmSoftsusy::doCalcTadpole1oneLoop(double mt, double sinthDRbar) const
{
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

double NmssmSoftsusy::doCalcTadpole2oneLoop(double mt, double sinthDRbar) const
{
 if (displayDrBarPars().mu(1, 3) == 0.0 || displayDrBarPars().mu(2, 3) == 0.0) {
   if (PRINTOUT > 1)
    cout << "Trying to calculate tadpole without having first calculated"
	 << " the DRbar masses.\n";
   return 0.0;
  }
  double g = displayGaugeCoupling(2),
    costhDRbar = cos(asin(sinthDRbar)), costhDRbar2 = sqr(costhDRbar),
    tanb = displayTanb(), sinb = sin(atan(tanb)),
    mu = -displaySusyMu(), q = displayMu(),
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

double NmssmSoftsusy::doCalcTadpoleSoneLoop(double mt, double sinthDRbar) const
{
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
  double mupr = displayMupr();
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
  DoubleMatrix mixh(3,3);
  if (mS.diagonaliseSym(mixh, mhsq) > TOLERANCE *
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
  //PA:  the diagonaliseSym gives us a mixing matrix A such that
  //A^T h^gauge = h^mass, but we want SLHA convention, so take transpose. 
  eg.mixh0 = mixh.transpose();
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
//PA:: fixes The CP odd mixing matrix with the conventions 
// Degrassi and Slavich arXiv:0907.4682
void NmssmSoftsusy::DegrassiSlavicMix(DoubleMatrix & P) const {
  double alphaP = displayDrBarPars().thetaH;
  DoubleMatrix Ppr(2, 2);
  double cb = cos(atan(displayTanb())), sb = sin(atan(displayTanb()));
  Ppr(1, 1) = cos(alphaP);
  Ppr(1, 2) = sin(alphaP);
  Ppr(2, 1) = - Ppr(1, 2);
  Ppr(2, 2) = cos(alphaP);
  
  P(1, 1) = - cb;           P(1, 2) = sb;             P(1, 3) = 0.0;
  P(2, 1) = sb * Ppr(1, 1); P(2, 2) = cb * Ppr(1, 1); P(2, 3) = Ppr(1, 2);
  P(3, 1) = sb * Ppr(2, 1); P(3, 2) = cb * Ppr(2, 1); P(3, 3) = Ppr(2, 2);
}


double NmssmSoftsusy::piZZTHiggs(double p, double q, 
				      double thetaWDRbar) const {
  //PA: NMSSM extensiion of BPMZ terms 
  //(only need to change mixing for new mass eigenstates, no new couplings)
  double    mz      = displayMzRun();
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
   double    beta    = atan(displayTanb());
  /// PA: 3 x 3 Higgs CP-even, S, and CP-odd, P, mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  	
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);
  double nmHiggs = 0.0;
  for(int i = 1; i <= 3; i++){ 
    nmHiggs = nmHiggs 
      + sqr(mz) * sqr(S(i, 1) * cos(beta) + S(i, 2) * sin(beta)) 
      * b0(p, mz, higgsm(i), q);//Z-Z-H
    	for (int j = 1; j <= 3; j++) {
	  nmHiggs =  nmHiggs - sqr(S(i, 1) * P(j, 1) - S(i, 2) * P(j, 2))
	    * b22bar(p, higgsm(i), higgsa(j), q); //CPodd and neut Goldstone
	}
  }	
  
  nmHiggs =  nmHiggs //charged Higgs
     - sqr(cos(2.0 * thetaWDRbar)) * b22bar(p, displayDrBarPars().mHpm, displayDrBarPars().mHpm, q);

  nmHiggs = nmHiggs
    - 2.0 * sqr(cw2DRbar) * (2 * sqr(p) + sqr(displayMwRun()) - sqr(mz) *
			     sqr(sw2DRbar) / cw2DRbar)
    * b0(p, displayMwRun(), displayMwRun(), q) //charged goldstone  
    - (8.0 * sqr(cw2DRbar) + sqr(cos(2.0 * thetaWDRbar))) * 
    b22bar(p, displayMwRun(), displayMwRun(), q); //charged Higgs
  

 return nmHiggs;
}

double NmssmSoftsusy::piZZTNeutralinos(double p, double q, 
					    double thetaWDRbar) const {
  
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    g       = displayGaugeCoupling(2);
  /// Neutralinos
  //static 
  double neutralinos = 0.0;
  ComplexMatrix aPsi(5, 5), bPsi(5, 5), aChi(5, 5), bChi(5, 5);
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);

  aPsi(3, 3) = g / (2.0 * cos(thetaWDRbar)); aPsi(4, 4) = -1. * aPsi(3, 3);
  bPsi = -1. * aPsi;
  
  aChi = n.complexConjugate() * aPsi * n.transpose();
  bChi = n * bPsi * n.hermitianConjugate();
  
  for (int i=1; i<=5; i++)
    for (int j=1; j<=5; j++) {
      neutralinos = neutralinos + cw2DRbar / (2.0 * sqr(g)) * 
	((sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod())) * 
	 hfn(p, mneut(i), mneut(j), q)
	 + 4.0 * (bChi(i, j).conj() * aChi(i, j)).real() *
	 mneut(i) * mneut(j) * b0(p, mneut(i), mneut(j), q)); 
    }
  
  return neutralinos;
}


double NmssmSoftsusy::piZZT(double p, double q, bool usePoleMt) const {
  
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    g       = displayGaugeCoupling(2);
  double rhs = 0.0;
 
  //PA: obtain Higgs contributions in separate method
  double higgs = piZZTHiggs(p, q, thetaWDRbar);
  //PA: obtain sfermion contributions in separate method
  double sfermions = piZZTsfermions(p, q);
  //PA: obtain fermion contributions in separate method
  double fermions = piZZTfermions(p, q, usePoleMt);
   //PA: obtain neutralino contributions in separate method
  double neutralinos = piZZTNeutralinos(p, q, thetaWDRbar);
   //PA: obtain neutralino contributions in separate method
  double charginos = piZZTCharginos(p, q, thetaWDRbar);
  
  
  rhs = higgs + charginos + neutralinos + fermions + sfermions ;

  double pi = rhs * sqr(g) / (cw2DRbar * 16.0 * sqr(PI));

  return pi;
}

double NmssmSoftsusy::piWWTHiggs(double p, double q, double thetaWDRbar) const {
  double    beta      = atan(displayTanb());
  double    cw2DRbar  = sqr(cos(thetaWDRbar));
  double    sw2DRbar  = 1.0 - cw2DRbar;
  double    mHc = displayDrBarPars().mHpm;
  double    mA = displayDrBarPars().mA0(1);
  double    mz      = displayMzRun();
  double sb = sin(beta), cb = cos(beta);
  /// PA: 3 x 3 Higgs CP-even, S, and CP-odd, P, mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cb;  C(1, 2) = sb; 
  C(2, 1) = C(1, 2); C(2, 2) = cb;	
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  double nmHiggs = 0.0;
  for(int i = 1; i <= 3; i++){ 
    nmHiggs =  nmHiggs  
      + sqr(displayMwRun()) * sqr(S(i, 1) * cb + S(i, 2) * sb)
      * b0(p, displayMwRun(), higgsm(i), q); //W-W-H
    for (int j = 1; j <= 2; j++) {
    nmHiggs = nmHiggs - sqr(S(i, 1) * C(j, 1) - S(i, 2) * C(j, 2)) 
    * b22bar(p, higgsm(i), higgsc(j), q);//includes goldstone
    nmHiggs = nmHiggs - sqr(P(i, 1) * C(j, 1) + P(i, 2) * C(j, 2)) 
    * b22bar(p, higgsa(i), higgsc(j), q);//includes goldstones
    }
  }
  nmHiggs +=  - (8.0 * cw2DRbar) * b22bar(p, mz, displayMwRun(), q)
    - sw2DRbar * (8.0 * b22bar(p, displayMwRun(), 0.0, q) + 4.0 * sqr(p) * 
		  b0(p, displayMwRun(), 0.0, q))  
    - ((4.0 * sqr(p) + sqr(mz) + sqr(displayMwRun())) * cw2DRbar - sqr(mz)  *
       sqr(sw2DRbar)) * b0(p, mz, displayMwRun(), q);
  cout << "Dev: nmHiggs = "  << nmHiggs << endl;
  return nmHiggs;
}

double NmssmSoftsusy::piWWTgauginos(double p, double q, double thetaWDRbar) const {
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    g       = displayGaugeCoupling(2);
  ComplexMatrix aPsi0PsicW(5, 2), bPsi0PsicW(5, 2), aChi0ChicW(5, 2),
    bChi0ChicW(5, 2);
  DoubleMatrix fW(5, 2), gW(5, 2);
  
  aPsi0PsicW(2, 1) = - g;
  bPsi0PsicW(2, 1) = - g;
  aPsi0PsicW(4, 2) = g / root2;		     
  bPsi0PsicW(3, 2) = -g / root2;		     
  
  ComplexMatrix aPsi(5, 5), bPsi(5, 5), aChi(5, 5), bChi(5, 5);
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 

  /// These ought to be in physpars
  aChi0ChicW = n.complexConjugate() * aPsi0PsicW * v.transpose();
  bChi0ChicW = n * bPsi0PsicW * u.hermitianConjugate();

  double gauginos = 0.0;

  for(int i=1;i<=5;i++)
    for(int j=1;j<=2;j++) {
      fW(i, j) = sqr(aChi0ChicW(i, j).mod()) + sqr(bChi0ChicW(i, j).mod());
      gW(i, j) = 2.0 * (bChi0ChicW(i, j).conj() * aChi0ChicW(i, j)).real(); 
      gauginos = gauginos + 
	(fW(i, j) * hfn(p, mneut(i), mch(j), q)
	 + 2.0 * gW(i, j) * mneut(i) * mch(j) * b0(p, mneut(i), mch(j), q)) 
	/ sqr(g);
    }
return gauginos;
}
double NmssmSoftsusy::piWWT(double p, double q, bool usePoleMt) const {

  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    g       = displayGaugeCoupling(2);

  double ans = 0.0;
  double higgs = piWWTHiggs(p, q, thetaWDRbar);
  double fermions = piWWTfermions(p, q, usePoleMt);   
  double sfermions = piWWTsfermions(p, q);   
  double gauginos = piWWTgauginos(p, q, thetaWDRbar);
  ans = higgs + sfermions + fermions + gauginos;
  cout << "Dev: gauginos = "  << gauginos << endl;
  cout << "Dev: fermions = "  << fermions << endl;
  cout << "Dev: sfermions = "  << sfermions << endl;
  cout << "Dev: higgs = "  << higgs << endl;
  double pi = ans * sqr(g) / (16.0 * sqr(PI));

  return pi;
}


//PA: Obtains trilnear couplings of s1-higgs-higgs for use in loop functions
void NmssmSoftsusy::getS1HiggsTriCoup(DoubleMatrix & sss1, DoubleMatrix & pps1,DoubleMatrix & hphps1, double thetaWDRbar) const {
  double lam =  displayLambda(), lsq = sqr(lam);
  double s = displaySvev();
  double kap = displayKappa();
  double al = displayTrialambda();
  double mupr = displayMupr();
  double tanthDrbar  = tan(thetaWDRbar);
  double cw2DRbar    = sqr(cos(thetaWDRbar));
  double tw2DRbar    = sqr(tanthDrbar);
  double g = displayGaugeCoupling(2), gsq = sqr(g);
  double beta = atan(displayTanb());
  double cb = cos(beta), cos2b = cos(2.0 * beta);
  double sb = sin(beta), sin2b = sin(2.0 * beta);
  double v1 =  displayHvev() * cb, v2 =  displayHvev() * sb; 
  /// LCT: Trilinear CP-even Higgs couplings to s1 in basis HdR HuR SR
  sss1(1, 1) = 0.125 * gsq / cw2DRbar * v1;
  sss1(2, 2) = v1 / 12.0 * (2.0 * lsq - 0.5 * gsq / cw2DRbar);
  sss1(3, 3) = lam / 6.0 * (lam * v1 - kap * v2);
  sss1(1, 2) = v2 / 12.0 * (2.0 * lsq - 0.5 * gsq / cw2DRbar);
  sss1(1, 3) = lsq * s / 6.0;
  sss1(2, 3) = - (al / root2 + lam * kap * s + lam * mupr / root2) / 6.0;
  sss1.symmetrise();
  
  /// LCT: Trilinear CP-odd Higgs couplings to s1 in basis HdI HuI SI
  pps1(1, 1) = 0.125 * gsq / cw2DRbar * v1;
  pps1(2, 2) = 0.25 * v1 * (2.0 * lsq - 0.5 * gsq / cw2DRbar);
  pps1(3, 3) = 0.5 * lam * (lam * v1 + kap * v2);
  pps1(2, 3) = 0.5 * (al / root2 - lam * kap * s - lam * mupr / root2);
  pps1.symmetrise();
  
  /// LCT: Trilinear with charged Higgs. Basis (G+ G- H+ H-)
  hphps1(1, 1) = 0.25 * gsq * (v1 * (1 + tw2DRbar * cos2b) - v2 * sin2b)
     + 0.5 * lsq * sin2b * v2;
  hphps1(2, 2) = 0.25 * gsq * (v1 * (1 - tw2DRbar * cos2b) + v2 * sin2b)
     - 0.5 * lsq * sin2b * v2;
  hphps1(1, 2) = - 0.25 * gsq * (v1 * tw2DRbar * sin2b + v2 * cos2b)
     + 0.5 * v2 * lsq * cos2b;
  hphps1(2, 1) = hphps1(1, 2);
  
}

//PA: Obtains trilnear couplings of s2-higgs-higgs for use in loop functions
void NmssmSoftsusy::getS2HiggsTriCoup(DoubleMatrix & sss2, DoubleMatrix & pps2, DoubleMatrix & hphps2, double thetaWDRbar) const {
  double lam =  displayLambda(), lsq = sqr(lam);
  double s = displaySvev();
  double kap = displayKappa();
  double al = displayTrialambda();
  double mupr = displayMupr();
  double tanthDrbar  = tan(thetaWDRbar);
  double cw2DRbar    = sqr(cos(thetaWDRbar));
  double tw2DRbar    = sqr(tanthDrbar);
  double g = displayGaugeCoupling(2), gsq = sqr(g);
  double beta = atan(displayTanb());
  double cb = cos(beta), cos2b = cos(2.0 * beta);
  double sb = sin(beta), sin2b = sin(2.0 * beta);
  double v1 =  displayHvev() * cb, v2 =  displayHvev() * sb; 
  /// LCT: Trilinear CP-even Higgs couplings to s2 in basis HdR HuR SR
  sss2(1, 1) = v2 * (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 12.0;
  sss2(2, 2) = 0.125 * v2 * gsq / cw2DRbar;
  sss2(3, 3) = lam / 6.0 * (lam * v2 - kap * v1);
  sss2(1, 2) = v1 * (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 12.0;
  sss2(1, 3) = - (al / root2 + lam * kap * s + lam * mupr / root2) / 6.0;
  sss2(2, 3) = lsq * s / 6.0;
  sss2.symmetrise();
  /// LCT: Trilinear CP-odd Higgs couplings to s2 in basis HdI HuI SI
  pps2(1, 1) = 0.25 * v2 * (2.0 * lsq - 0.5 * gsq / cw2DRbar);
  pps2(2, 2) = 0.125 * gsq / cw2DRbar * v2;
  pps2(3, 3) = 0.5 * lam * (lam * v2 + kap * v1);
  pps2(1, 3) = 0.5 * (al / root2 - lam * kap * s - lam * mupr / root2);
  pps2.symmetrise();
  /// LCT: Trilinear coupling with charged Higgs. Basis (G+ G- H+ H-)
  hphps2(1, 1) = 0.25 * gsq * (v2 * (1.0 - tw2DRbar * cos2b) - v1 * sin2b)
     + 0.5 * lsq * v1 * sin2b;
  hphps2(2, 2) = 0.25 * gsq * (v2 * (1.0 + tw2DRbar * cos2b) + v1 * sin2b)
     - 0.5 * lsq * v1 * sin2b;
  hphps2(1, 2) = 0.25 * gsq * (v2 * tw2DRbar * sin2b - v1 * cos2b)
     + 0.5 * lsq * v1 * cos2b;
  hphps2(2, 1) = hphps2(1, 2);

}

//PA: Obtains trilnear couplings of s2-higgs-higgs for use in loop functions
void NmssmSoftsusy::getS3HiggsTriCoup(DoubleMatrix & sss3, DoubleMatrix & pps3, DoubleMatrix & hphps3) const {
  double lam =  displayLambda(), lsq = sqr(lam);
  double s = displaySvev();
  double kap = displayKappa(), ksq = sqr(kap);
  double al = displayTrialambda();
  double ak = displayTriakappa();
  double mupr = displayMupr();
  double beta = atan(displayTanb());
  double cb = cos(beta), cos2b = cos(2.0 * beta);
  double sb = sin(beta), sin2b = sin(2.0 * beta);
  double v1 =  displayHvev() * cb, v2 =  displayHvev() * sb; 
  /// LCT: Trilinear CP-even Higgs couplings to s3 in basis HdR HuR SR
  sss3(1, 1) = lsq / 6.0 * s;
  sss3(2, 2) = lsq / 6.0 * s;
  sss3(3, 3) = (ak / 3.0 + root2 * ksq * s + kap * mupr) / (root2); 
  sss3(1, 2) = -(al / root2 + lam * kap * s + lam * mupr / root2) / (6.0);
  sss3(1, 3) = lam * (lam * v1 - kap * v2) / (6.0);
  sss3(2, 3) = lam * (lam * v2 - kap * v1) / (6.0);
  sss3.symmetrise();

	
  /// LCT: Trilinear CP-odd Higgs couplings to s3 in basis HdI HuI SI
  pps3(1, 1) = 0.5 * lsq * s;
  pps3(2, 2) = pps3(1, 1);
  pps3(3, 3) = (-ak + root2 * ksq * s + kap * mupr) / (root2);
  pps3(1, 2) = 0.5 * (al / root2 + lam * kap * s + lam * mupr / root2);
  pps3(1, 3) = -0.5 * lam * kap * v2;
  pps3(2, 3) = -0.5 * lam * kap * v1;
  pps3.symmetrise();
	
  /// LCT: Trilinear with charged Higgs. Basis (G+ G- H+ H-)
  hphps3(1, 1) = 0.5 * (2.0 * lsq * s
                        - (root2 * al + 2.0 * lam * kap * s) * sin2b);
  hphps3(2, 2) = 0.5 * (2.0 * lsq * s 
                        + (root2 * al + 2.0 * lam * kap * s) * sin2b); 
  hphps3(1, 2) = -0.5 * (root2 * al + 2.0 * lam * kap * s) * cos2b;
  hphps3(2, 1) = hphps3(1, 2);
}





double NmssmSoftsusy::pis1s1Higgs(double p, double q) const {
  double beta = atan(displayTanb()); 
  double thetaWDRbar = asin(calcSinthdrbar());
  double tanthDrbar  = tan(thetaWDRbar);
  double cw2DRbar    = sqr(cos(thetaWDRbar));
  double tw2DRbar    = sqr(tanthDrbar); 
  double lam =  displayLambda(), lsq = sqr(lam);
  double g = displayGaugeCoupling(2), gsq = sqr(g);
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta),
    sinb = sin(beta); 
  double mw = displayMwRun(), mw2 = sqr(mw), mz = displayMzRun(), mz2 = sqr(mz);
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
   /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);
  
/// LCT: Charged Higgs/Goldstone parts unchanged in NMSSM.
  double higgs = 0.0;	
  for (int i=1; i <= 2; i++) {
     higgs = higgs + gsq * 0.5 * sqr(C(i, 1)) * ffn(p, higgsc(i), mw, q);
	}
  
/// LCT: CP-odd states in basis G0 A1 A2
  for (int i = 1; i <= 3; i++) {
     higgs = higgs + gsq * 0.25 / cw2DRbar * sqr(P(i, 1)) 
        * ffn(p, higgsa(i), mz, q);
	} 
	
  /// LCT: Gauge bosons
  higgs = higgs + 1.75 * gsq * cosb2 * (2.0 * mw2 * b0(p, mw, mw, q) 
                                        + mz2 * b0(p, mz, mz, q) / cw2DRbar) 
     + gsq * (2.0 * a0(mw, q) + a0(mz, q) / cw2DRbar);

  //PA: trilinear couplings for s1 to CP even, CP odd and charged Higgs
  DoubleMatrix sss1(3, 3), pps1(3, 3), hphps1(2, 2);
  getS1HiggsTriCoup(sss1, pps1, hphps1, thetaWDRbar);
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
  /// LCT: Quadrilinear CP-even Higgs couplings 
  DoubleMatrix sss1s1(3, 3);
  sss1s1(1, 1) = gsq / (32.0 * cw2DRbar);
  sss1s1(2, 2) = (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 48.0;
  sss1s1(3, 3) = lsq / 24.0;
  sss1s1.symmetrise();
  
  /// LCT: Quadrilinear CP-odd Higgs couplings
  DoubleMatrix pps1s1(3, 3);
  pps1s1(1, 1) = gsq / (16.0 * cw2DRbar);
  pps1s1(2, 2) = (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 8.0;
  pps1s1(3, 3) = 0.25 * lsq;
  pps1s1.symmetrise();
  
  /// LCT: Rotate to mass bases s1 s1 Hi Hi and s1 s1 Ai Ai
  DoubleVector hhs1s1(3), aas1s1(3);
  for (int i = 1; i <=3; i++) {
     for (int a = 1; a <= 3; a++) {
        for (int b = 1; b <=3; b++) {
           hhs1s1(i) = hhs1s1(i) + 6.0 * S(i, a) * S(i, b) * sss1s1(a, b);
           aas1s1(i) = aas1s1(i) + P(i, a) * P(i, b) * pps1s1(a, b);
        }
     }
  }
  
  for (int i=1; i<=3; i++) {
     for (int j=1; j<=3; j++) {
        higgs = higgs + 2.0 * sqr(hhs1(i, j)) * b0(p, higgsm(i), higgsm(j), q);
        higgs = higgs + 2.0 * sqr(aas1(i, j)) * b0(p, higgsa(i), higgsa(j), q);
     }
     higgs = higgs + 2.0 * hhs1s1(i) * a0(higgsm(i), q);
     higgs = higgs + 2.0 * aas1s1(i) * a0(higgsa(i), q);
  }
  
  /// LCT: Quadrilinear (G+ H+)
  DoubleVector hphps1s1(2);
  hphps1s1(1) = gsq * (1.0 + tw2DRbar * cos2b) / 8.0;
  hphps1s1(2) = gsq * (1.0 - tw2DRbar * cos2b) / 8.0;

  for(int i=1; i<=2; i++) {
     for(int j=1; j<=2; j++) {
        higgs = higgs + sqr(hphps1(i, j)) * b0(p, higgsc(i), higgsc(j), q);
     }
     higgs = higgs + 2.0 * hphps1s1(i) * a0(higgsc(i), q);
  }

   return higgs;
}

double NmssmSoftsusy::pis1s2Higgs(double p, double q) const {
  double beta = atan(displayTanb()); 
  double thetaWDRbar = asin(calcSinthdrbar());
  double cw2DRbar    = sqr(cos(thetaWDRbar));
  double lam =  displayLambda(), lsq = sqr(lam);
  double kap = displayKappa();
  double g = displayGaugeCoupling(2), gsq = sqr(g);
  double cosb = cos(beta), sinb = sin(beta), sin2b = sin(2.0 * beta); 
  double mw = displayMwRun(), mw2 = sqr(mw), mz = displayMzRun(), mz2 = sqr(mz);
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
   /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);
  double higgs = 0.0;
  /// LCT: Gauge bosons
  higgs = higgs + 3.5 * gsq * sinb * cosb * (mw2 * b0(p, mw, mw, q) 
        + 0.5 / cw2DRbar * mz2 * b0(p, mz, mz, q));
  /// PA: Charged Higgs-Gauge contribution
  for (int i=1; i <= 2; i++) {
     higgs = higgs - 0.5 * gsq * C(i, 1) * C(i, 2) * ffn(p, higgsc(i), mw, q);
  }
  /// PA: CP-odd Higgs-Gauge contribution
  for (int i=1; i <= 3; i++) {
     higgs = higgs - 0.25 * gsq * P(i, 1) * P(i, 2) 
        * ffn(p, higgsa(i), mz, q) / cw2DRbar;
  }
  //PA: trilinear couplings for s1, s2 to CP even, CP odd and charged Higgs
  DoubleMatrix sss1(3, 3), pps1(3, 3), hphps1(2, 2);
  getS1HiggsTriCoup(sss1, pps1, hphps1, thetaWDRbar);
  DoubleMatrix sss2(3, 3), pps2(3, 3), hphps2(2, 2);
  getS2HiggsTriCoup(sss2, pps2, hphps2, thetaWDRbar);
  /// LCT: Rotate to mass bases s1 Hi Hj, s2 Hi Hj, s1 Ai Aj, s2 Ai Aj
  DoubleMatrix hhs1(3, 3), hhs2(3, 3), aas1(3, 3), aas2(3, 3);
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <= 3; j++) {
        for (int a = 1; a <= 3; a++) {
           for (int b = 1; b <= 3; b++) {
              hhs1(i, j) = hhs1(i, j) + 3.0 * S(i, a) * S(j, b) * sss1(a, b);
              hhs2(i, j) = hhs2(i, j) + 3.0 * S(i, a) * S(j, b) * sss2(a, b);
              aas1(i, j) = aas1(i, j) + P(i, a) * P(j, b) * pps1(a, b);
              aas2(i, j) = aas2(i, j) + P(i, a) * P(j, b) * pps2(a, b);
           }
        }
     }
  }
  
  /// LCT: CP-even quadrilinear Higgs couplings in basis HdR HuR SR
  DoubleMatrix sss1s2(3, 3), pps1s2(3, 3);
  sss1s2(1, 2) = (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 48.0;
  sss1s2(3, 3) = - lam * kap / 24.0;
  sss1s2.symmetrise();
  /// LCT: CP-odd quadrilinear Higgs couplings in basis HdI HuI SI
  pps1s2(3, 3) = 0.25 * lam * kap;
  pps1s2.symmetrise();
  
  /// LCT: Rotate to mass bases hhs1s2 and aas1s2
  DoubleVector hhs1s2(3), aas1s2(3);
  
  for (int i=1; i <= 3; i++) {
     for (int a=1; a <= 3; a++) {
        for (int b=1; b <= 3; b++) {
           hhs1s2(i) = hhs1s2(i) + 6.0 * S(i, a) * S(i, b) * sss1s2(a, b);
           aas1s2(i) = aas1s2(i) + P(i, a) * P(i, b) * pps1s2(a, b);
        }
     }
  }
  /// Trilinear and quadrilinear contributions fron neutral Higgs
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <= 3; j++) {
        higgs = higgs + 2.0 * hhs1(i, j) * hhs2(i, j) 
           * b0(p, higgsm(i), higgsm(j), q) + 2.0 * aas1(i, j) * aas2(i, j) 
           * b0(p, higgsa(i), higgsa(j), q);
     }
     higgs = higgs + 2.0 * (hhs1s2(i) * a0(higgsm(i), q) 
                            + aas1s2(i) * a0(higgsa(i), q));
  }

  /// Quadrilinear couplings for Charged Higgs
  DoubleVector hphps1s2(2);
  hphps1s2(1) = (2.0 * lsq - gsq) * sin2b / 8.0;
  hphps1s2(2) = - hphps1s2(1);
  
  for (int i=1; i <=2; i++) {
     for (int j=1; j <=2; j++) {
        higgs = higgs + hphps1(i, j) * hphps2(j, i) *
           b0(p, higgsc(i), higgsc(j), q);
     }
     higgs = higgs + 2.0 * hphps1s2(i) * a0(higgsc(i), q);
  }
  
  
  return higgs;
}

double NmssmSoftsusy::pis2s2Higgs(double p, double q) const {
  double beta = atan(displayTanb()); 
  double thetaWDRbar = asin(calcSinthdrbar());
  double tanthDrbar  = tan(thetaWDRbar);
  double cw2DRbar    = sqr(cos(thetaWDRbar));
  double tw2DRbar    = sqr(tanthDrbar); 
  double lam =  displayLambda(), lsq = sqr(lam);
  double g = displayGaugeCoupling(2), gsq = sqr(g);
  double cosb = cos(beta), cos2b = cos(2.0 * beta),
     sinb = sin(beta), sinb2 = sqr(sinb); 
  double mw = displayMwRun(), mw2 = sqr(mw), mz = displayMzRun(), mz2 = sqr(mz);
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);
  
  double higgs = 0.0;	
 /// LCT: Charged Higgs/Goldstone parts unchanged in NMSSM.
  for (int i=1; i <= 2; i++) {
     higgs = higgs + gsq * 0.5 * sqr(C(i, 2)) * ffn(p, higgsc(i), mw, q);
	}
	
  /// LCT: CP-odd states in basis G0 A1 A2
  for (int i = 1; i <= 3; i++) {
     higgs = higgs + gsq * 0.25 / cw2DRbar * sqr(P(i, 2)) 
        * ffn(p, higgsa(i), mz, q);
	} 
	
  /// LCT: Gauge bosons
  higgs = higgs + 1.75 * gsq * sinb2 * (2.0 * mw2 * b0(p, mw, mw, q) 
                                        + mz2 * b0(p, mz, mz, q) / cw2DRbar) 
     + gsq * (2.0 * a0(mw, q) + a0(mz, q) / cw2DRbar);
  
  //PA: trilinear couplings for s2 to CP even, CP odd and charged Higgs
  DoubleMatrix sss2(3, 3), pps2(3, 3), hphps2(2, 2);
  getS2HiggsTriCoup(sss2, pps2, hphps2, thetaWDRbar);
  /// LCT: Rotate to mass basis s2 Hi Hj
  DoubleMatrix hhs2(3, 3), aas2(3, 3);
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <= 3; j++) {
        for (int a = 1; a <= 3; a++) {
           for (int b = 1; b <= 3; b++) {
              hhs2(i, j) = hhs2(i, j) + 3.0 * S(i, a) * S(j, b) * sss2(a, b);
              aas2(i, j) = aas2(i, j) + P(i, a) * P(j, b) * pps2(a, b);
           }
        }
     }
  }
  
  /// LCT: Quadrilinear CP-even Higgs couplings 
  DoubleMatrix sss2s2(3, 3);
  sss2s2(1, 1) = (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 48.0;
  sss2s2(2, 2) = gsq / (32.0 * cw2DRbar);
  sss2s2(3, 3) = lsq / 24.0;
  sss2s2.symmetrise();
  
  /// LCT: Quadrilinear CP-odd Higgs couplings
  DoubleMatrix pps2s2(3, 3);
  pps2s2(1, 1) = (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 8.0;
  pps2s2(2, 2) = gsq / (16.0 * cw2DRbar);
  pps2s2(3, 3) = 0.25 * lsq;
  pps2s2.symmetrise();
  
  /// LCT: Rotate to mass bases s1 s1 Hi Hi and s1 s1 Ai Ai
  DoubleVector hhs2s2(3), aas2s2(3);
  for (int i = 1; i <= 3; i++) {
     for (int a = 1; a <= 3; a++) {
        for (int b = 1; b <= 3; b++) {
           hhs2s2(i) = hhs2s2(i) + 6.0 * S(i, a) * S(i, b) * sss2s2(a, b);
           aas2s2(i) = aas2s2(i) + P(i, a) * P(i, b) * pps2s2(a, b);
        }
     }
  }
  
  for (int i=1; i<=3; i++) {
     for (int j=1; j<=3; j++) {
        higgs = higgs + 2.0 * sqr(hhs2(i, j)) * b0(p, higgsm(i), higgsm(j), q);
        higgs = higgs + 2.0 * sqr(aas2(i, j)) * b0(p, higgsa(i), higgsa(j), q);
    }
     higgs = higgs + 2.0 * hhs2s2(i) * a0(higgsm(i), q);
     higgs = higgs + 2.0 * aas2s2(i) * a0(higgsa(i), q);
  }
  
  /// LCT: Quadrilinear (H1+ H1-) (Amended for NMSSM)
  DoubleVector hphps2s2(2);
  hphps2s2(1) = gsq * (1.0 - tw2DRbar * cos2b) / 8.0;
  hphps2s2(2) = gsq * (1.0 + tw2DRbar * cos2b) / 8.0;

  for (int i=1; i <= 2; i++) {
     for (int j=1; j <= 2; j++) {
        higgs = higgs + sqr(hphps2(i, j)) * b0(p, higgsc(i), higgsc(j), q);
		}
     higgs = higgs + 2.0 * hphps2s2(i) * a0(higgsc(i), q);
  }
  
  
  return higgs;
}

double NmssmSoftsusy::pis1s3Higgs(double p, double q) const {
  double lam =  displayLambda(), lsq = sqr(lam);
  // double s = displaySvev();
  double kap = displayKappa();
  double beta = atan(displayTanb()); 
  double cosb = cos(beta), sinb = sin(beta); 
  double thetaWDRbar = asin(calcSinthdrbar());
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);
  
  double higgs = 0.0;	
  //PA: trilinear couplings for s1, s3 to CP even, CP odd and charged Higgs
  DoubleMatrix sss1(3, 3), pps1(3, 3), hphps1(2, 2);
  getS1HiggsTriCoup(sss1, pps1, hphps1, thetaWDRbar);
  DoubleMatrix sss3(3, 3), pps3(3, 3), hphps3(2, 2);
  getS3HiggsTriCoup(sss3, pps3, hphps3);

  /// LCT: Rotate to mass basis s1 Hi Hj, s1 Ai Aj, s3 Hi Hj, s3 Ai Aj
  DoubleMatrix hhs1(3, 3), aas1(3, 3), hhs3(3, 3), aas3(3, 3);
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <=3; j++) {
        for (int a = 1; a <= 3; a++) {
           for (int b = 1; b <= 3; b++) {
              hhs1(i, j) = hhs1(i, j) + 3.0 * S(i, a) * S(j, b) * sss1(a, b);
              aas1(i, j) = aas1(i, j) + P(i, a) * P(j, b) * pps1(a, b);
              hhs3(i, j) = hhs3(i, j) + 3.0 * S(i, a) * S(j, b) * sss3(a, b);
              aas3(i, j) = aas3(i, j) + P(i, a) * P(j, b) * pps3(a, b);
           }
        }
     }
  }
  
  /// LCT: CP-even quadrilinear Higgs couplings to s1s3 in basis HdR HuR SR
  DoubleMatrix sss1s3(3, 3), pps1s3(3, 3);
  sss1s3(1, 3) = lsq / 24.0;
  sss1s3(2, 3) = -lam * kap / 24.0;
  sss1s3.symmetrise();
  /// LCT: CP-odd quadrilinear Higgs couplings to s1s3 in basis HdI HuI SI
  pps1s3(2, 3) = -0.25 * lam * kap;
  pps1s3.symmetrise();
  
  /// LCT: Rotate to mass bases hhs1s3 and aas1s3
  DoubleVector hhs1s3(3), aas1s3(3);
  
  for (int i=1; i <= 3; i++) {
     for (int a=1; a <= 3; a++) {
        for (int b=1; b <= 3; b++) {
           hhs1s3(i) = hhs1s3(i) + 6.0 * S(i, a) * S(i, b) * sss1s3(a, b);
           aas1s3(i) = aas1s3(i) + P(i, a) * P(i, b) * pps1s3(a, b);
        }
     }
  }  
  
  for(int i=1; i <= 3; i++) {
     for(int j=1; j <= 3; j++) {
        higgs = higgs + 2.0 * hhs1(i, j) * hhs3(i, j) 
           * b0(p, higgsm(i), higgsm(j), q)
           + 2.0 * aas1(i, j) * aas3(i, j) * b0(p, higgsa(i), higgsa(j), q);
        
     }
     higgs = higgs + 2.0 * (hhs1s3(i) * a0(higgsm(i), q) 
                            + aas1s3(i) * a0(higgsa(i), q));
  }
  
  /// LCT: No charged quadrilinears 
  //PA: get trilinear charged Higgs contribution
  for (int i=1; i <=2; i++) {
     for (int j=1; j <=2; j++) {
        higgs = higgs + hphps1(i, j) * hphps3(j, i) *
           b0(p, higgsc(i), higgsc(j), q);
     }
  }
return higgs;
}


double NmssmSoftsusy::pis2s3Higgs(double p, double q) const {
  double lam =  displayLambda(), lsq = sqr(lam);
  // double s = displaySvev();
  double kap = displayKappa();
  double beta = atan(displayTanb()); 
  double cosb = cos(beta), sinb = sin(beta); 
  double thetaWDRbar = asin(calcSinthdrbar());
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);
  
  double higgs = 0.0;	
  //PA: trilinear couplings for s2, s3 to CP even, CP odd and charged Higgs
  DoubleMatrix sss2(3, 3), pps2(3, 3), hphps2(2, 2);
  getS2HiggsTriCoup(sss2, pps2, hphps2, thetaWDRbar);
  DoubleMatrix sss3(3, 3), pps3(3, 3), hphps3(2, 2);
  getS3HiggsTriCoup(sss3, pps3, hphps3);

  /// LCT: Rotate to mass basis s2 Hi Hj, s2 Ai Aj, s3 Hi Hj, s3 Ai Aj
  DoubleMatrix hhs2(3, 3), aas2(3, 3), hhs3(3, 3), aas3(3, 3);
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <= 3; j++) {
        for (int a = 1; a <= 3; a++) {
           for (int b = 1; b <= 3; b++) {
              hhs2(i, j) = hhs2(i, j) + 3.0 * S(i, a) * S(j, b) * sss2(a, b);
              aas2(i, j) = aas2(i, j) + P(i, a) * P(j, b) * pps2(a, b);
              hhs3(i, j) = hhs3(i, j) + 3.0 * S(i, a) * S(j, b) * sss3(a, b);
              aas3(i, j) = aas3(i, j) + P(i, a) * P(j, b) * pps3(a, b);
           }
        }
     }
  }
  
  /// LCT: CP-even quadrilinear Higgs couplings to s2s3 in basis HdR HuR SR
  DoubleMatrix sss2s3(3, 3), pps2s3(3, 3);
  sss2s3(1, 3) = -lam * kap / 24.0;
  sss2s3(2, 3) = lsq / 24.0;
  sss2s3.symmetrise();
  /// LCT: CP-odd quadrilinear Higgs couplings to s2s3 in basis HdI HuI SI
  pps2s3(1, 3) = -0.25 * lam * kap;
  pps2s3.symmetrise();
	
  DoubleVector hhs2s3(3), aas2s3(3);
  
  for (int i=1; i <= 3; i++) {
     for (int a=1; a <= 3; a++) {
        for (int b=1; b <= 3; b++) {
           hhs2s3(i) = hhs2s3(i) + 6.0 * S(i, a) * S(i, b) * sss2s3(a, b);
           aas2s3(i) = aas2s3(i) + P(i, a) * P(i, b) * pps2s3(a, b);
        }
     }
  }  
  
  /// Trilinear and quadrilinear contributions
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <= 3; j++) {
        higgs = higgs + 2.0 * hhs2(i, j) * hhs3(i, j) 
           * b0(p, higgsm(i), higgsm(j), q) + 2.0 * aas2(i, j) * aas3(i, j) 
           * b0(p, higgsa(i), higgsa(j), q);
     }
     higgs = higgs + 2.0 * (hhs2s3(i) * a0(higgsm(i), q) 
                            + aas2s3(i) * a0(higgsa(i), q));
  }

  /// LCT: No charged quadrilinears present
  for (int i=1; i <= 2; i++) {
     for (int j=1; j <= 2; j++) {
        higgs = higgs + hphps2(i, j) * hphps3(j, i) *
           b0(p, higgsc(i), higgsc(j), q);
     }
  }

  return higgs;
}


double NmssmSoftsusy::pis3s3Higgs(double p, double q) const {
  double lam =  displayLambda(), lsq = sqr(lam);
  // double s = displaySvev();
  double kap = displayKappa(), ksq = sqr(kap);
  double beta = atan(displayTanb()); 
  double cosb = cos(beta), sinb = sin(beta),  sin2b = sin(2.0 * beta); 
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);
  
  double higgs = 0.0;	
  //PA: trilinear couplings for s3 to CP even, CP odd and charged Higgs
  DoubleMatrix sss3(3, 3), pps3(3, 3), hphps3(2, 2);
  getS3HiggsTriCoup(sss3, pps3, hphps3);

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
  
  /// LCT: Quadrilinear CP-even Higgs couplings to s3s3
  DoubleMatrix sss3s3(3, 3);
  sss3s3(1, 1) = lsq / 24.0;
  sss3s3(2, 2) = lsq / 24.0;
  sss3s3(3, 3) = 0.25 * ksq;
  sss3s3(1, 2) = -lam * kap / 24.0;
  sss3s3.symmetrise();
  
  /// LCT: Quadrilinear CP-odd Higgs couplings to s3s3
  DoubleMatrix pps3s3(3, 3);
  pps3s3(1, 1) = lsq / 4.0;
  pps3s3(2, 2) = lsq / 4.0;
  pps3s3(3, 3) = 0.5 * ksq;
  pps3s3(1, 2) = 0.25 * lam * kap;
  pps3s3.symmetrise();
  
  /// LCT: Rotate to mass bases s3 s3 Hi Hi and s3 s3 Ai Ai
  DoubleVector hhs3s3(3), aas3s3(3);
  for (int i = 1; i <= 3; i++) {
     for (int a = 1; a <= 3; a++) {
        for (int b = 1; b <= 3; b++) {
           hhs3s3(i) = hhs3s3(i) + 6.0 * S(i, a) * S(i, b) * sss3s3(a, b);
           aas3s3(i) = aas3s3(i) + P(i, a) * P(i, b) * pps3s3(a, b);
        }
     }
  }
 
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <= 3; j++) {
        higgs = higgs + 2.0 * sqr(hhs3(i, j)) * b0(p, higgsm(i), higgsm(j), q);
        higgs = higgs + 2.0 * sqr(aas3(i, j)) * b0(p, higgsa(i), higgsa(j), q);
     }
     higgs = higgs + 2.0 * hhs3s3(i) * a0(higgsm(i), q);
     higgs = higgs + 2.0 * aas3s3(i) * a0(higgsa(i), q);
  }

  /// LCT: Quadrilinear (G+ H-) 
  DoubleVector hphps3s3(2);
  hphps3s3(1) = 0.5 * lam * (lam - kap * sin2b);
  hphps3s3(2) = 0.5 * lam * (lam + kap * sin2b);
  
  for (int i=1; i <= 2; i++) {
     for (int j=1; j <= 2; j++) {
        higgs = higgs + sqr(hphps3(i, j)) * b0(p, higgsc(i), higgsc(j), q);
     }
     higgs = higgs + 2.0 * hphps3s3(i) * a0(higgsc(i), q);
  }

  return higgs;
}



void NmssmSoftsusy::getS1NeutralinoCoup(ComplexMatrix & aChi, 
                                        ComplexMatrix & bChi) const {
  double g   = displayGaugeCoupling(2);
  double gp  = displayGaugeCoupling(1) * sqrt(0.6);
  double lam = displayLambda();
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleMatrix aPsi(5, 5);
  aPsi(1, 3) = -gp * 0.5; 
  aPsi(2, 3) = g * 0.5; 
  aPsi(4, 5) = - lam / root2;
  aPsi.symmetrise();
  aChi = n.complexConjugate() * aPsi * n.hermitianConjugate();
  bChi = n * aPsi * n.transpose();
}

void NmssmSoftsusy::getS2NeutralinoCoup(ComplexMatrix & aChi, 
                                        ComplexMatrix & bChi) const {
  double g   = displayGaugeCoupling(2);
  double gp  = displayGaugeCoupling(1) * sqrt(0.6);
  double lam = displayLambda();
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleMatrix aPsi(5, 5);
  aPsi(1, 4) = gp * 0.5; 
  aPsi(2, 4) = -g * 0.5; 
  aPsi(3, 5) = - lam / root2;
  aPsi.symmetrise();
  aChi = n.complexConjugate() * aPsi * n.hermitianConjugate();
  bChi = n * aPsi * n.transpose();
}

void NmssmSoftsusy::getS3NeutralinoCoup(ComplexMatrix & aChi, 
                                        ComplexMatrix & bChi) const {
  double kap = displayKappa();
  double lam = displayLambda();
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleMatrix aPsi(5, 5);
  aPsi(3, 4) = -lam / root2;
  aPsi(5, 5) = kap * root2;
  aPsi.symmetrise();
  aChi = n.complexConjugate() * aPsi * n.hermitianConjugate();
  bChi = n * aPsi * n.transpose();
}


double NmssmSoftsusy::pis1s1Neutralinos(double p, double q) const {
  
  double neutralinos = 0.0;
  ComplexMatrix aChi(5, 5), bChi(5, 5);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  getS1NeutralinoCoup(aChi, bChi); 
  DoubleMatrix fChiChis1s1(5, 5), gChiChis1s1(5, 5);
  for(int i=1; i<=5; i++)
    for (int j=1; j<=5; j++) {
      fChiChis1s1(i, j) = sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod());
      gChiChis1s1(i, j) = (bChi(i, j).conj() * aChi(i, j) + 
	aChi(i, j).conj() * bChi(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
	(fChiChis1s1(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
	 gChiChis1s1(i, j) * mneut(i) * mneut(j) * 
	 b0(p, mneut(i), mneut(j), q));
    }

  return neutralinos;
}

double NmssmSoftsusy::pis1s2Neutralinos(double p, double q) const {
  ComplexMatrix aChi1(5, 5), bChi1(5, 5);
  ComplexMatrix aChi2(5, 5), bChi2(5, 5);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  getS1NeutralinoCoup(aChi1, bChi1);
  getS2NeutralinoCoup(aChi2, bChi2);
  double neutralinos = 0.0;
  DoubleMatrix fChiChis1s2(5, 5), gChiChis1s2(5, 5);
  for(int i=1; i<=5; i++){
     for(int j=1; j<=5; j++) {
      fChiChis1s2(i, j) = (aChi1(i, j).conj() * aChi2(i, j) + 
                           bChi1(i, j).conj() * bChi2(i, j)).real();
      gChiChis1s2(i, j) = (bChi1(i, j).conj() * aChi2(i, j) + 
                           aChi1(i, j).conj() * bChi2(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
         (fChiChis1s2(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
          gChiChis1s2(i, j) * mneut(i) * mneut(j) * 
          b0(p, mneut(i), mneut(j), q));
    }
  }
return neutralinos;
}


double NmssmSoftsusy::pis2s2Neutralinos(double p, double q) const {
  
  double neutralinos = 0.0;
  ComplexMatrix aChi(5, 5), bChi(5, 5);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  getS2NeutralinoCoup(aChi, bChi); 
  DoubleMatrix fChiChis1s1(5, 5), gChiChis1s1(5, 5);
  for(int i=1; i<=5; i++)
    for (int j=1; j<=5; j++) {
      fChiChis1s1(i, j) = sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod());
      gChiChis1s1(i, j) = (bChi(i, j).conj() * aChi(i, j) + 
	aChi(i, j).conj() * bChi(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
	(fChiChis1s1(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
	 gChiChis1s1(i, j) * mneut(i) * mneut(j) * 
	 b0(p, mneut(i), mneut(j), q));
    }

  return neutralinos;
}

double NmssmSoftsusy::pis1s3Neutralinos(double p, double q) const {
  ComplexMatrix aChi1(5, 5), bChi1(5, 5);
  ComplexMatrix aChi3(5, 5), bChi3(5, 5);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  getS1NeutralinoCoup(aChi1, bChi1);
  getS3NeutralinoCoup(aChi3, bChi3);
  double neutralinos = 0.0;
  DoubleMatrix fChiChis1s3(5, 5), gChiChis1s3(5, 5);
  for(int i=1; i<=5; i++)
     for (int j=1; j<=5; j++) {
        fChiChis1s3(i, j) = (aChi1(i, j).conj() * aChi3(i, j) + 
                             bChi1(i, j).conj() * bChi3(i, j)).real();
        gChiChis1s3(i, j) = (bChi1(i, j).conj() * aChi3(i, j) + 
                             aChi1(i, j).conj() * bChi3(i, j)).real();
        neutralinos = neutralinos + 0.5 * 
           (fChiChis1s3(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
            gChiChis1s3(i, j) * mneut(i) * mneut(j) * 
            b0(p, mneut(i), mneut(j), q));
    }

  return neutralinos;
}

double NmssmSoftsusy::pis2s3Neutralinos(double p, double q) const {
  ComplexMatrix aChi2(5, 5), bChi2(5, 5);
  ComplexMatrix aChi3(5, 5), bChi3(5, 5);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  getS2NeutralinoCoup(aChi2, bChi2);
  getS3NeutralinoCoup(aChi3, bChi3);
  double neutralinos = 0.0;
  DoubleMatrix fChiChis2s3(5, 5), gChiChis2s3(5, 5);
  for(int i=1; i<=5; i++)
     for(int j=1; j<=5; j++) {
        fChiChis2s3(i, j) = (aChi2(i, j).conj() * aChi3(i, j) + 
                             bChi2(i, j).conj() * bChi3(i, j)).real();
        gChiChis2s3(i, j) = (bChi2(i, j).conj() * aChi3(i, j) + 
                             aChi2(i, j).conj() * bChi3(i, j)).real();
        neutralinos = neutralinos + 0.5 * 
           (fChiChis2s3(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
            gChiChis2s3(i, j) * mneut(i) * mneut(j) * 
            b0(p, mneut(i), mneut(j), q));
    }

  return neutralinos;
}



double NmssmSoftsusy::pis3s3Neutralinos(double p, double q) const {
  
  double neutralinos = 0.0;
  ComplexMatrix aChi(5, 5), bChi(5, 5);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  getS3NeutralinoCoup(aChi, bChi); 
  DoubleMatrix fChiChis1s1(5, 5), gChiChis1s1(5, 5);
  for(int i=1; i<=5; i++)
    for (int j=1; j<=5; j++) {
      fChiChis1s1(i, j) = sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod());
      gChiChis1s1(i, j) = (bChi(i, j).conj() * aChi(i, j) + 
	aChi(i, j).conj() * bChi(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
	(fChiChis1s1(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
	 gChiChis1s1(i, j) * mneut(i) * mneut(j) * 
	 b0(p, mneut(i), mneut(j), q));
    }

  return neutralinos;
}


double NmssmSoftsusy::pis1s3Charginos(double p, double q) const {
  double chargino = 0.0;
  double g = displayGaugeCoupling(2);
  double lam = displayLambda();
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 
  DoubleMatrix aPsic1(2, 2), aPsic3(2, 2);
  aPsic1(1, 2) = g / root2; 
  ComplexMatrix aChic1(2, 2), bChic1(2, 2);
  ComplexMatrix aChic3(2, 2), bChic3(2, 2);
  aChic1 = v.complexConjugate() * aPsic1 * u.hermitianConjugate();
  bChic1 = u * aPsic1.transpose() * v.transpose();
  aPsic3(2, 2) = lam / root2;
  aChic3 = v.complexConjugate() * aPsic3 * u.hermitianConjugate();
  bChic3 = u * aPsic3.transpose() * v.transpose();
   DoubleMatrix fChiChis1s3(2, 2), gChiChis1s3(2, 2);
  for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      fChiChis1s3(i, j) = (aChic1(i, j).conj() * aChic3(i, j) + 
													 bChic1(i, j).conj() * bChic3(i, j)).real();
      gChiChis1s3(i, j) = (bChic1(i, j).conj() * aChic3(i ,j) + 
													 aChic1(i, j).conj() * bChic3(i, j)).real();
      chargino = chargino + 
			(fChiChis1s3(i, j) * gfn(p, mch(i), mch(j), q) - 2.0 *
			 gChiChis1s3(i, j) * mch(i) * mch(j) * 
			 b0(p, mch(i), mch(j), q));
    }	
			
  return chargino;
}
 
double NmssmSoftsusy::pis2s3Charginos(double p, double q) const {
  double chargino = 0.0;
  double g = displayGaugeCoupling(2);
  double lam = displayLambda();
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 
  DoubleMatrix aPsic2(2, 2), aPsic3(2, 2);
  aPsic2(2, 1) = g / root2; 
  ComplexMatrix aChic2(2, 2), bChic2(2, 2);
  ComplexMatrix aChic3(2, 2), bChic3(2, 2);
  aChic2 = v.complexConjugate() * aPsic2 * u.hermitianConjugate();
  bChic2 = u * aPsic2.transpose() * v.transpose();
  aPsic3(2, 2) = lam / root2;
  aChic3 = v.complexConjugate() * aPsic3 * u.hermitianConjugate();
  bChic3 = u * aPsic3.transpose() * v.transpose();
  DoubleMatrix fChiChis2s3(2, 2), gChiChis2s3(2, 2);
  for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      fChiChis2s3(i, j) = (aChic2(i, j).conj() * aChic3(i, j) + 
													 bChic2(i, j).conj() * bChic3(i, j)).real();
      gChiChis2s3(i, j) = (bChic2(i, j).conj() * aChic3(i ,j) + 
													 aChic2(i, j).conj() * bChic3(i, j)).real();
      chargino = chargino + 
			(fChiChis2s3(i, j) * gfn(p, mch(i), mch(j), q) - 2.0 *
			 gChiChis2s3(i, j) * mch(i) * mch(j) * 
			 b0(p, mch(i), mch(j), q));
    }	
  return chargino;
}


double NmssmSoftsusy::pis3s3Charginos(double p, double q) const {
  double chargino = 0.0;
  double lam = displayLambda();
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 
  DoubleMatrix aPsic(2, 2);
  aPsic(2, 2) = lam / root2; 
  ComplexMatrix aChic(2, 2), bChic(2, 2);
  aChic = v.complexConjugate() * aPsic * u.hermitianConjugate();
  bChic = u * aPsic.transpose() * v.transpose();
   DoubleMatrix fChiChis3s3(2, 2), gChiChis3s3(2, 2);
  for(int i=1; i <= 2; i++)
    for (int j=1; j <= 2; j++) {
      fChiChis3s3(i, j) = sqr(aChic(i, j).mod()) + sqr(bChic(i, j).mod());
      gChiChis3s3(i, j) = (bChic(i, j).conj() * aChic(i, j) + 
													 aChic(i, j).conj() * bChic(i, j)).real();
      chargino = chargino + 
			(fChiChis3s3(i, j) * gfn(p, mch(i), mch(j), q) - 2.0 *
			 gChiChis3s3(i, j) * mch(i) * mch(j) * 
			 b0(p, mch(i), mch(j), q));
    }

 return chargino;
}

double NmssmSoftsusy::pis1s3Sfermions(double p, double q, DoubleMatrix ls1tt, DoubleMatrix ls1bb,  DoubleMatrix ls1tautau, DoubleMatrix ls3tt,  DoubleMatrix ls3bb,  DoubleMatrix ls3tautau) const {
  double lam      = displayLambda();
  double thetat   = displayDrBarPars().thetat ;
  double thetab   = displayDrBarPars().thetab;
  double thetatau = displayDrBarPars().thetatau;
  double st       = sin(thetat) ;
  double ct       = cos(thetat) ;
  double ht       = displayDrBarPars().ht;
  double mstop1   = displayDrBarPars().mu(1, 3);
  double mstop2   = displayDrBarPars().mu(2, 3);
  /// Mix 3rd family up
  ls1tt = rot2d(thetat) * ls1tt * rot2d(-thetat);
  ls1bb = rot2d(thetab) * ls1bb * rot2d(-thetab);
  ls1tautau = rot2d(thetatau) * ls1tautau * rot2d(-thetatau);
  ls3tt = rot2d(thetat) * ls3tt * rot2d(-thetat);
  ls3bb = rot2d(thetab) * ls3bb * rot2d(-thetab);
  ls3tautau = rot2d(thetatau) * ls3tautau * rot2d(-thetatau);
  
  double sfermions = 0.0;
  double stops = -3.0 * ht * lam * ct * st * (a0(mstop1, q) - a0(mstop2, q));
  int i, j; for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      sfermions = sfermions + 3.0 * ls1tt(i, j) * ls3tt(i, j) *
         b0(p, displayDrBarPars().mu(i, 3), displayDrBarPars().mu(j, 3), q);
      sfermions = sfermions + 3.0 * ls1bb(i, j) * ls3bb(i, j) * 
         b0(p, displayDrBarPars().md(i, 3), displayDrBarPars().md(j, 3), q);
      sfermions = sfermions + ls1tautau(i, j) * ls3tautau(i, j) * 
         b0(p, displayDrBarPars().me(i, 3), displayDrBarPars().me(j, 3), q);
      
    }


  sfermions = sfermions + stops;

  return sfermions;
}

double NmssmSoftsusy::pis2s3Sfermions(double p, double q, DoubleMatrix ls2tt, DoubleMatrix ls2bb,  DoubleMatrix ls2tautau, DoubleMatrix ls3tt,  DoubleMatrix ls3bb,  DoubleMatrix ls3tautau) const {
       
  double lam     = displayLambda();
  double msbot1  = displayDrBarPars().md(1, 3);
  double msbot2  = displayDrBarPars().md(2, 3);
  double mstau1  = displayDrBarPars().me(1, 3);
  double mstau2  = displayDrBarPars().me(2, 3);
  double thetat  = displayDrBarPars().thetat ;
  double thetab  = displayDrBarPars().thetab;
  double thetatau= displayDrBarPars().thetatau;
  double cb      = cos(thetab) ;
  double ctau    = cos(thetatau);
  double sb      = sin(thetab) ;
  double stau    = sin(thetatau);
  
  double hb = displayDrBarPars().hb, htau = displayDrBarPars().htau;
  /// Quadrilinear contributions	
  double sbots = -3.0 * hb * lam * cb * sb * (a0(msbot1, q) - a0(msbot2, q));
  double staus =  -1.0 * htau * lam * ctau * stau * (a0(mstau1, q) - a0(mstau2, q));

   /// Mix 3rd family up
  ls2tt = rot2d(thetat) * ls2tt * rot2d(-thetat);
  ls2bb = rot2d(thetab) * ls2bb * rot2d(-thetab);
  ls2tautau = rot2d(thetatau) * ls2tautau * rot2d(-thetatau);
  ls3tt = rot2d(thetat) * ls3tt * rot2d(-thetat);
  ls3bb = rot2d(thetab) * ls3bb * rot2d(-thetab);
  ls3tautau = rot2d(thetatau) * ls3tautau * rot2d(-thetatau);
  
  double sfermions = 0.0;
  for (int i=1; i<=2; i++)
     for (int j=1; j<=2; j++) {
        sfermions = sfermions + 3.0 * ls2tt(i, j) * ls3tt(i, j) 
           * b0(p, displayDrBarPars().mu(i, 3), displayDrBarPars().mu(j, 3), q);
        sfermions = sfermions + 3.0 * ls2bb(i, j) * ls3bb(i, j) 
           * b0(p, displayDrBarPars().md(i, 3), displayDrBarPars().md(j, 3), q);
        sfermions = sfermions + ls2tautau(i, j) * ls3tautau(i, j) 
           * b0(p, displayDrBarPars().me(i, 3), displayDrBarPars().me(j, 3), q);
               }
	
  sfermions =  sfermions + sbots + staus; 


  return sfermions;
}

double NmssmSoftsusy::pis3s3Sfermions(double p, double q, DoubleMatrix ls3tt,  DoubleMatrix ls3bb,  DoubleMatrix ls3tautau) const {
       
 
  double thetat  = displayDrBarPars().thetat ;
  double thetab  = displayDrBarPars().thetab;
  double thetatau= displayDrBarPars().thetatau;

/// Mix 3rd family up
  ls3tt = rot2d(thetat) * ls3tt * rot2d(-thetat);
  ls3bb = rot2d(thetab) * ls3bb * rot2d(-thetab);
  ls3tautau = rot2d(thetatau) * ls3tautau * rot2d(-thetatau);
  
  double sfermions = 0.0;
  int i, j; for (i=1; i <= 2; i++)
    for (j=1; j <= 2; j++) {
      sfermions = sfermions + 3.0 * sqr(ls3tt(i, j)) 
         * b0(p, displayDrBarPars().mu(i, 3), displayDrBarPars().mu(j, 3), q);
      sfermions = sfermions + 3.0 * sqr(ls3bb(i, j)) 
         * b0(p, displayDrBarPars().md(i, 3), displayDrBarPars().md(j, 3), q);
      sfermions = sfermions + sqr(ls3tautau(i, j)) 
         * b0(p, displayDrBarPars().me(i, 3), displayDrBarPars().me(j, 3), q);
    }

 return sfermions;
}
double NmssmSoftsusy::pis1s1(double p, double q) const {
  double    beta    = atan(displayTanb());
  double    mb      = displayDrBarPars().mb;
  double    hb      = displayDrBarPars().hb;
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
    /// minus sign taken into acct here!
  double smu  = -displaySusyMu(); 
  double mz   = displayMzRun();
  double g    = displayGaugeCoupling(2);
  double cosb = cos(beta);
  //PA: get fermion contribution, uses MSSM version
  double fermions = pis1s1Fermions(p, q);
  // sfermion couplings to s1 Higgs state
  DoubleMatrix ls1tt(2, 2), ls1bb(2, 2), ls1tautau(2, 2);
  double gmzOcthW = g * mz / costhDrbar;
  H1SfSfCouplings(ls1tt, ls1bb, ls1tautau, gmzOcthW, smu, cosb, root2*mb/hb);
  //PA: get sfermion contribution
  double sfermions = pis1s1Sfermions(p, q, ls1tt, ls1bb, ls1tautau);
  //PA: get Higgs contribution
  double higgs = pis1s1Higgs(p, q);
  /// Neutralino contribution
  double neutralinos = pis1s1Neutralinos(p, q);
  /// Chargino contribution
  double chargino = pis1s1Charginos(p, q);  

  return 
    (sfermions + 
     fermions + higgs + neutralinos + chargino) / (16.0 * sqr(PI));
}


double NmssmSoftsusy::pis1s2(double p, double q) const {
  double    beta    = atan(displayTanb());
  double    mb   =  displayDrBarPars().mb, hb   =  displayDrBarPars().hb; 
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    smu     = -displaySusyMu(); /// minus sign taken into acct here!
  double    g       = displayGaugeCoupling(2);
  double cosb = cos(beta), sinb = sin(beta);
  double  mz = displayMzRun();
  // sfermion couplings to s1 Higgs state(NMSSM version)
  DoubleMatrix ls1tt(2, 2), ls1bb(2, 2), ls1tautau(2, 2);
  double gmzOcthW = g * mz / costhDrbar;
  H1SfSfCouplings(ls1tt, ls1bb, ls1tautau, gmzOcthW, smu, cosb, root2*mb/hb);
  /// sfermion couplings to s2 Higgs state
  DoubleMatrix ls2tt(2, 2), ls2bb(2, 2), ls2tautau(2, 2);
  H2SfSfCouplings(ls2tt, ls2bb, ls2tautau, gmzOcthW, smu, sinb);
  //PA: get sfermion contribution
  double sfermions = pis1s2Sfermions(p, q, ls1tt, ls1bb, ls1tautau, ls2tt, ls2bb, ls2tautau);
  //PA: get Higgs contribution (NMSSM version)
  double higgs = pis1s2Higgs(p, q);
  /// Neutralino contribution (NMSSM version)
  double neutralinos = pis1s2Neutralinos(p, q); 
  /// Chargino contribution
  double chargino = pis1s2Charginos(p, q);  
 
  return (sfermions + higgs + neutralinos + chargino) 
    / (16.0 * sqr(PI));
}

double NmssmSoftsusy::pis2s2(double p, double q) const {
  double beta = atan(displayTanb());
  double thetaWDRbar = asin(calcSinthdrbar());
  double costhDrbar = cos(thetaWDRbar);
  double smu = -displaySusyMu(); /// minus sign taken into acct here!
  double g = displayGaugeCoupling(2);
  double sinb = sin(beta);
  double mz = displayMzRun();
  double gmzOcthW = g * mz / costhDrbar;
  //PA: get femions (same as MSSM) 
  double fermions = pis2s2Fermions(p, q);
  ///PA: sfermion couplings to s2 Higgs state (NMSSM version)
  DoubleMatrix ls2tt(2, 2), ls2bb(2, 2), ls2tautau(2, 2);
  H2SfSfCouplings(ls2tt, ls2bb, ls2tautau, gmzOcthW, smu, sinb);
  //PA: get sfermions (MSSM routine but with NMSSM couplings input)
  double sfermions = pis2s2Sfermions(p, q, ls2tt, ls2bb, ls2tautau);
  //PA: get higgs (NMSSM version)
  double higgs = pis2s2Higgs(p, q);
  //PA: get neutralinos (NMSSM version)
  double neutralinos = pis2s2Neutralinos(p, q); 
  //PA: get charginos (same as in MSSM version)
  double chargino = pis2s2Charginos(p, q);   

  return (fermions + sfermions + higgs + neutralinos + chargino) 
    / (16.0 * sqr(PI));
}



double NmssmSoftsusy::pis1s3(double p, double q) const {
  double beta = atan(displayTanb());
  double mb = displayDrBarPars().mb;
  double hb = displayDrBarPars().hb;
  double thetaWDRbar = asin(calcSinthdrbar());
  double costhDrbar  = cos(thetaWDRbar);
  double lam = displayLambda();         
    /// minus sign taken into acct here!
  double smu  = -displaySusyMu(); 
  double mz   = displayMzRun();
  double g    = displayGaugeCoupling(2);
  double cosb = cos(beta);
 // sfermion couplings to s1 Higgs state(NMSSM version)
  DoubleMatrix ls1tt(2, 2), ls1bb(2, 2), ls1tautau(2, 2);
  double gmzOcthW = g * mz / costhDrbar;
  H1SfSfCouplings(ls1tt, ls1bb, ls1tautau, gmzOcthW, smu, cosb, root2*mb/hb);
  /// sfermion couplings to s3 Higgs state
  DoubleMatrix ls3tt(2, 2), ls3bb(2, 2), ls3tautau(2, 2);
  SSfSfCouplings(ls3tt, ls3bb, ls3tautau, lam);
  double sfermions = pis1s3Sfermions(p, q, ls1tt, ls1bb, ls1tautau, ls3tt, ls3bb, ls3tautau);
  double higgs =  pis1s3Higgs(p, q);
  double neutralinos = pis1s3Neutralinos(p, q);
  double chargino = pis1s3Charginos(p, q);
return (sfermions + higgs + neutralinos + chargino) 
    / (16.0 * sqr(PI));

}


double NmssmSoftsusy::pis2s3(double p, double q) const {
  double beta = atan(displayTanb());
  double thetaWDRbar = asin(calcSinthdrbar());
  double costhDrbar  = cos(thetaWDRbar);
  double lam = displayLambda();         
    /// minus sign taken into acct here!
  double smu  = -displaySusyMu(); 
  double mz   = displayMzRun();
  double g    = displayGaugeCoupling(2);
  double sinb = sin(beta);
 // sfermion couplings to s1 Higgs state(NMSSM version)
  DoubleMatrix ls2tt(2, 2), ls2bb(2, 2), ls2tautau(2, 2);
  double gmzOcthW = g * mz / costhDrbar;
  H2SfSfCouplings(ls2tt, ls2bb, ls2tautau, gmzOcthW, smu, sinb);
  /// sfermion couplings to s3 Higgs state
  DoubleMatrix ls3tt(2, 2), ls3bb(2, 2), ls3tautau(2, 2);
  SSfSfCouplings(ls3tt, ls3bb, ls3tautau, lam);
  double sfermions = pis2s3Sfermions(p, q, ls2tt, ls2bb, ls2tautau, ls3tt, ls3bb, ls3tautau);
  double higgs =  pis2s3Higgs(p, q);
  double neutralinos = pis2s3Neutralinos(p, q);
  double chargino = pis2s3Charginos(p, q);
  return (sfermions + higgs + neutralinos + chargino) 
    / (16.0 * sqr(PI));
}

double NmssmSoftsusy::pis3s3(double p, double q) const {
  double lam = displayLambda();         
  /// sfermion couplings to s3 Higgs state
  DoubleMatrix ls3tt(2, 2), ls3bb(2, 2), ls3tautau(2, 2);
  SSfSfCouplings(ls3tt, ls3bb, ls3tautau, lam);
  double sfermions = pis3s3Sfermions(p, q, ls3tt, ls3bb, ls3tautau);
  double higgs =  pis3s3Higgs(p, q);
  double neutralinos = pis3s3Neutralinos(p, q);
  double chargino = pis3s3Charginos(p, q);
return (sfermions + higgs + neutralinos + chargino) 
    / (16.0 * sqr(PI));

}



#endif
