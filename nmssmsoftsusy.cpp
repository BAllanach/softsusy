/*file nmssmsoftsusy.cpp
    Project: NMSSMSOFTSUSY
    Author: Ben Allanach, Peter Athron, Lewis Tunstall, Alexander Voigt
    Manual: TBW
    Webpage:  https://github.com/Expander/softsusy.git
*/

#include "nmssmsoftsusy.h"
#ifdef NMSSMSOFTSUSY_H

extern double sw2, gnuL, guL, gdL, geL, guR, gdR, geR, yuL, yuR, ydL,
  ydR, yeL, yeR, ynuL;


double NmssmSoftsusy::testSlavichpiZZT(double g, double gp, double ht, double hb, double htau, double v1, double v2, double p, double Q) const {
   double piZZT = 0.0;
   getpizz_(&g, &gp, &ht, &hb, &htau, &v1, &v2, &p, &Q, &piZZT);
   return piZZT;
 }

double NmssmSoftsusy::testSlavichpiWWT(double g, double gp, double ht, double hb, double htau, double v1, double v2, double p, double Q) const {
   double piWWT = 0.0;
   getpiww_(&g, &gp, &ht, &hb, &htau, &v1, &v2, &p, &Q, &piWWT);
   return piWWT;
 }

double NmssmSoftsusy::getpiSS(double g,double gp, double ll, double kk, double ht, double hb, double htau, double v1, double v2, double xx, double Ak, double Al, double At, double Ab, double Atau, double p, double Q, int i, int j) const {
   double piSS[3][3];
   for(int ii=0; ii<=2; ii++){
      for(int jj=0; jj<=2; jj++){
         piSS[ii][jj] = 0;
      }
   }

getpiss_(&g,&gp, &ll, &kk, &ht, &hb, &htau, &v1, &v2, &xx, &Ak, &Al, &At, &Ab, &Atau, &p, &Q, &piSS); 
  return piSS[i-1][j-1];
 }



double NmssmSoftsusy::getpiPP(double g,double gp, double ll, double kk, double ht, double hb, double htau, double v1, double v2, double xx, double Ak, double Al, double At, double Ab, double Atau, double p, double Q, int i, int j) const {
   double piPP[3][3];
   for(int ii=0; ii<=2; ii++){
      for(int jj=0; jj<=2; jj++){
         piPP[ii][jj] = 0;
      }
   }

   getpipp_(&g,&gp, &ll, &kk, &ht, &hb, &htau, &v1, &v2, &xx, &Ak, &Al, &At, 
	    &Ab, &Atau, &p, &Q, &piPP); 
  return piPP[i-1][j-1];
 }


void NmssmSoftsusy::testSlavichTreeMasses() const {

  double gp = sqrt(0.6) * displayGaugeCoupling(1);
  double g = displayGaugeCoupling(2);
  double lam = displayLambda();
  double kap = displayKappa();
  double Al  = displaySoftAlambda();
  double Ak  = displaySoftAkappa();
  double At  = displaySoftA(UA, 3, 3);
  double Ab  = displaySoftA(DA, 3, 3);
  double Atau  = displaySoftA(EA, 3, 3);
  double mQ3 = sqrt(displaySoftMassSquared(mQl, 3, 3)); 
  double mQ = sqrt(displaySoftMassSquared(mQl, 1, 1));
  if (mQ - sqrt(displaySoftMassSquared(mQl, 2, 2)) > 1e-9){
    cout << "WARNING: called when mQ1 not equal to mQ2" << endl;
  } 
  double mL3 = sqrt(displaySoftMassSquared(mLl, 3, 3)); 
  double mL = sqrt(displaySoftMassSquared(mLl, 1, 1));
  if (mL - sqrt(displaySoftMassSquared(mLl, 2, 2)) > 1e-9){
    cout << "WARNING: called when mL1 not equal to mL2" << endl;
  } 
  double mtr = sqrt(displaySoftMassSquared(mUr, 3, 3)); 
  double mur = sqrt(displaySoftMassSquared(mUr, 1, 1));
  if (mur - sqrt(displaySoftMassSquared(mUr, 2, 2)) > 1e-9){
    cout << "WARNING: called when mU1 not equal to mU2" << endl;
  } 

  double mbr = sqrt(displaySoftMassSquared(mDr, 3, 3));
  double mdr = sqrt(displaySoftMassSquared(mDr, 1, 1));
  if (mdr - sqrt(displaySoftMassSquared(mDr, 2, 2)) > 1e-9){
    cout << "WARNING: called when mD1 not equal to mD2" << endl;
  } 
  double mtaur = sqrt(displaySoftMassSquared(mEr, 3, 3)); 
  double mer = sqrt(displaySoftMassSquared(mEr, 1, 1));
  if (mer - sqrt(displaySoftMassSquared(mEr, 2, 2)) > 1e-9){
    cout << "WARNING: called when mE1 not equal to mE2" << endl;
  } 
  double s   = displaySvev() / root2;
  double ht = displayDrBarPars().ht;
  double hb = displayDrBarPars().hb;
  double htau = displayDrBarPars().htau;
  //divide by root2 for slavich conventions.
  double v1 = displayHvev() * cos(atan(displayTanb())) / root2;
  double v2 = displayHvev() * sin(atan(displayTanb())) / root2;
  bool errmass = false;
  double M1 = displayGaugino(1);
  double M2 = displayGaugino(2);
  double Q = displayMu();

  treemasses_(&g,&gp,&lam,&kap,&ht,&hb,&htau,&v1,&v2,&s,&M1,&M2, &Ak,&Al,&At,&Ab,&Atau,&mQ3,&mtr,&mbr,&mQ,&mur,&mdr,&mL3,&mtaur,&mL,&mer,&Q,&errmass); 

 if(errmass) cout << "Warning tachyon in slavich treemasses."  << endl;
 return;
}


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


//PA: calls routines to calculate all three tadpoles and sets them.
// Currently only works at one loop.  
// Two loop should be added later. 
void NmssmSoftsusy::doTadpoles(double mt, double sinthDRbar) {
   //PA: all the MSSM parts could be replaced by a single call to 
   //Softsusy::doTadpoles and then add MSSM ones
   //But we need to call some 2loop MSSM routines anyway
   // to get 2loop S tadpoles
    calcTadpole1Ms1loop(mt, sinthDRbar);
    calcTadpole2Ms1loop(mt, sinthDRbar);
    calcTadpoleSMs1loop(mt, sinthDRbar);
    //PA: set to one loop first
    double t1OV1 = displayTadpole1Ms1loop();
    double t2OV2 = displayTadpole2Ms1loop();
    double tSOVS = displayTadpoleSMs1loop();
    //  if (numRewsbLoops > 1 && displayProblem().tachyon == none && Z3 == true) {
    if (numRewsbLoops > 1 && Z3 == true) {
       double lam = displayLambda(), s = displaySvev();
       const drBarPars& forLoops = displayDrBarPars();
       /// add the two-loop terms, prepare inputs
       double s1s = 0., s2s = 0., s1t = 0., s2t = 0.,
         gs = displayGaugeCoupling(3),  as = sqr(gs) / (4.0 * PI), 
         rmt = forLoops.mt, rmtsq = sqr(forLoops.mt), 
         scalesq = sqr(displayMu()), vev = displayHvev(),
         vev2 = sqr(vev), tanb = displayTanb(), 
         amu = - lam * s / root2, mg = displayGaugino(3), 
         mAsq = sqr(forLoops.mA0(1)); 
      
      double sxt = sin(forLoops.thetat), cxt = cos(forLoops.thetat);
      double mst1sq = sqr(forLoops.mu(1, 3)), 
	mst2sq = sqr(forLoops.mu(2, 3));
      /// two-loop Higgs corrections: alpha_s alpha_b
      double sxb = sin(forLoops.thetab), 
	cxb = cos(forLoops.thetab);
      double sintau = sin(forLoops.thetatau), 
	costau = cos(forLoops.thetatau);
      double msnusq = sqr(forLoops.msnu(3));
      double msb1sq = sqr(forLoops.md(1, 3)), 
	msb2sq = sqr(forLoops.md(2, 3));
      double mstau1sq = sqr(forLoops.me(1, 3)), 
	mstau2sq = sqr(forLoops.me(2, 3));
      double cotbeta = 1.0 / tanb;
      double rmb = forLoops.mb, rmbsq = sqr(forLoops.mb);
      double rmtausq = sqr(forLoops.mtau);
      double s1b = 0.0, s2b = 0.0, s1tau = 0.0, s2tau = 0.0;
      
      ewsb2loop_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, 
        	 &amu, &tanb, &vev2, &gs, &s1s, &s2s);
      ddstad_(&rmtsq, &rmbsq, &mAsq, &mst1sq, &mst2sq, &msb1sq, &msb2sq, 
              &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &s1t, 
      	      &s2t);
      ewsb2loop_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq,
        	 &amu, &cotbeta, &vev2, &gs, &s2b, &s1b);
      tausqtad_(&rmtausq, &mAsq, &msnusq, &mstau1sq, &mstau2sq, &sintau, 
        	&costau, &scalesq, &amu, &tanb, &vev2, &s1tau, &s2tau);
      
      //rescale T1 to get TS
      double sss = s1s * vev * cos(atan(tanb)) / s;
      double ssb = s1b * vev * sin(atan(tanb)) / s;
   
      if (!testNan(s1s * s1t * s1b * s1tau * s2s * s2t * s2b * s2tau 
                   * sss * ssb)) {
	 t1OV1 += - s1s - s1t - s1b - s1tau;
	 t2OV2 += - s2s - s2t - s2b - s2tau;
         tSOVS += - sss - ssb;   
         /// end of 2-loop bit
      }
      else  {
	flagNoMuConvergence(true);
	if (PRINTOUT > 1) cout << "2-loop tadpoles are nans\n";
      }

    }

    //PA: After one and two loop tadpoles are added they are then set
    setT1OV1Ms(t1OV1); 
    setT2OV2Ms(t2OV2); 
    tSOVSMs = tSOVS;
}

void NmssmSoftsusy::P1SfSfCouplings(DoubleMatrix & lp1tt, DoubleMatrix & lp1bb, DoubleMatrix  & lp1tautau) const {
  double s = displaySvev();
  double lam = displayLambda();
  double ht = displayDrBarPars().ht;
  double ab = displayDrBarPars().ub;
  double atau = displayDrBarPars().utau;
  
  lp1tt(2, 1) = - ht / root2 
    * (displaySusyMu() + lam * s / root2);
  lp1tt(1, 2) = -lp1tt(2, 1);
  
  lp1bb(1, 2) = ab / root2;
  lp1bb(2, 1) = -lp1bb(1, 2);
  
  lp1tautau(1, 2) = atau / root2;
  lp1tautau(2, 1) = -lp1tautau(1, 2); 

}

void NmssmSoftsusy::P2SfSfCouplings(DoubleMatrix & lp2tt, DoubleMatrix & lp2bb, DoubleMatrix  & lp2tautau) const {
  double s = displaySvev();
  double lam = displayLambda();
  double hb = displayDrBarPars().hb;
  double at = displayDrBarPars().ut;
  double htau = displayDrBarPars().htau;

  lp2tt(1, 2) = at / root2;
  lp2tt(2, 1) = -lp2tt(1, 2);
  lp2bb(2, 1) = - hb / root2 * (displaySusyMu() + lam * s / root2);
  lp2bb(1, 2) = -lp2bb(2, 1);
  lp2tautau(2, 1) = - htau / root2 * (displaySusyMu() + lam * s / root2);
  lp2tautau(1, 2) = -lp2tautau(2, 1); 
}

void NmssmSoftsusy::P3SfSfCouplings(DoubleMatrix & lp3tt, DoubleMatrix & lp3bb, DoubleMatrix  & lp3tautau) const {
  double lam = displayLambda();
  double ht = displayDrBarPars().ht;
  double hb = displayDrBarPars().hb;
  double htau = displayDrBarPars().htau;
  double v1 = displayHvev() * cos(atan(displayTanb()));
  double v2 = displayHvev() * sin(atan(displayTanb()));

  lp3tt(2, 1) = - 0.5 * ht * v1 * lam;
  lp3tt(1, 2) = - lp3tt(2, 1);
  lp3bb(2, 1) = - 0.5 * hb * v2 * lam;
  lp3bb(1, 2) = - lp3bb(2, 1);
  lp3tautau(2, 1) = - 0.5 * htau * v2 * lam;
  lp3tautau(1, 2) = - lp3tautau(2, 1);

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
  const drBarPars& fL = displayDrBarPars();
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
  const drBarPars& f = displayDrBarPars();

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
  double sfermions = Softsusy<SoftParsNmssm>::doCalcTad1Sfermions(lTS1Lr, lBS1Lr, lTauS1Lr, costhDRbar);
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
  double sfermions =  Softsusy<SoftParsNmssm>::doCalcTad2Sfermions(lTS2Lr, lBS2Lr, lTauS2Lr, costhDRbar);
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

//PA: functions to set tadpoles.  We could remove these and exploit the virtual nature of these routines, so that when used by an NmssmSoftsusy object, calcTadpole routines would point to the NmssmSoftsusy doCalc routines.  This would work the same as now but the drawback would be the code may be less easy to read.
void NmssmSoftsusy::calcTadpole1Ms1loop(double mt, double sinthDRbar) {
   double t1OV1 = doCalcTadpole1oneLoop(mt, sinthDRbar);
  if (testNan(t1OV1)) {
    flagNoMuConvergence(true);
    t1OV1 = 0.0;
  }
  setT1OV1Ms1loop(t1OV1);
}

void NmssmSoftsusy::calcTadpole2Ms1loop(double mt, double sinthDRbar) {
   double t2OV2 = doCalcTadpole2oneLoop(mt, sinthDRbar);
  if (testNan(t2OV2)) {
    flagNoMuConvergence(true);
    t2OV2 = 0.0;
  }
  setT2OV2Ms1loop(t2OV2);
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
void NmssmSoftsusy::treeCharginos(DoubleMatrix & mass, double beta, double mw) {
  double lam = displayLambda(), svev = displaySvev();

  Softsusy<SoftParsNmssm>::treeCharginos(mass, beta, mw);
  mass(2, 2) = mass(2, 2) + lam * svev / root2;
}

/// LCT: new routine for NMSSM neutralino masses
void NmssmSoftsusy::treeNeutralinos(DoubleMatrix & mass, double beta, double mz, double mw, double sinthDRbar) {
  double lam = displayLambda(), kap = displayKappa();
  double mupr = displayMupr();
  double cosb = cos(beta), sinb = sin(beta);
  double vev = displayHvev(), svev = displaySvev();

  /// Call MSSM 4 x4 neutralino mass matrix
  Softsusy<SoftParsNmssm>::treeNeutralinos(mass, beta, mz, mw, sinthDRbar);

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
  treeCharginos(mCh, beta, mw);
  eg.mch = mCh.asy2by2(eg.thetaL, eg.thetaR);
  eg.mpzCharginos();

  DoubleMatrix mNeut(5, 5);
  treeNeutralinos(mNeut, beta, mz, mw, sinth);
  if (mNeut.diagonaliseSym(eg.mixNeut, eg.mneut) > TOLERANCE *
       1.0e-3) {
      ostringstream ii;
      ii << "accuracy bad in neutralino diagonalisation"<< flush;
      throw ii.str();
   }

  eg.mpzNeutralinos();
}


// PA: fills tree level CP even and CP odd Higgs mass matrices 
// and tree level mHPm all with treel level EWSB conditions used.
// Called in higgs and calcDrBarParsHiggs
// CP even mass matrix in (Hd, Hu, S) basis  
// CP odd Higgs mass matrices mPpr in (Hd, Hu, S) basis 
// and mP2 in rotated basis (A, S) -- goldstone boson removed      
void NmssmSoftsusy::treeHiggs(DoubleMatrix & mS, DoubleMatrix & mPpr, 
                               DoubleMatrix & mP2, double & mHpmsq, 
                               double beta) const {
  double tanb = displayTanb();
  double sinb = sin(beta), cosb = cos(beta);
  double sinb2 = sqr(sinb), cosb2 = sqr(cosb);
  double mz = displayMzRun(), mz2 = sqr(mz);
  double mw = displayMwRun(), mw2 = sqr(mw);
	/// LCT: NMSSM parameters
  double lam = displayLambda(), kap = displayKappa();
  double mupr = displayMupr(), smu = displaySusyMu();
  double al = displayTrialambda(), ak = displayTriakappa();
  double xiF = displayXiF(), xiS = displayXiS();
  double m3sq = displayM3Squared(), mSpsq = displayMspSquared();

  double vev = displayHvev(), v1 = vev * cosb, v2 = vev * sinb;
  double s = displaySvev();

  double mueff = lam * s / root2 + smu;
  double m3hatsq = m3sq + lam * (mupr * s / root2 + xiF);
 
  /// CP-even Higgs in EHT notation and basis (HdR, HuR, SR)
  mS(1, 1) = mz2 * cosb2 + (0.5 * lam * kap * sqr(s)
                            + al * s / root2 + m3hatsq) * tanb;
  mS(1, 2) = - sinb * cosb * mz2 + sqr(lam) * sqr(vev) * sinb * cosb
     - (0.5 * lam * kap * sqr(s) + al * s / root2 + m3hatsq);
  mS(1, 3) = root2 * lam * mueff * v1 
     - (al + lam * kap * s * root2 + lam * mupr) * v2 / root2;
  mS(2, 2) = mz2 * sinb2 + (0.5 * lam * kap * sqr(s)
                            + al * s / root2 + m3hatsq) / tanb;
  mS(2, 3) = root2 * lam * mueff * v2 
     - (al + lam * kap * s * root2 + lam * mupr) * v1 / root2;
  mS(3, 3) = (al + lam * mupr) * v1 * v2 / (root2 * s)
     + s / root2 * (ak + 2.0 * root2 * sqr(kap) * s
                       + 3.0 * kap * mupr) 
     - root2 * (xiS + xiF * mupr + 0.5 * lam * smu * sqr(vev)) / s;

  mS.symmetrise();

  /// CP-odd Higgs in EHT notation and basis (HdI, HuI, SI)
  mPpr(1, 1) = (0.5 * lam * kap * sqr(s) + al * s / root2
                + m3hatsq) * tanb;
  mPpr(1, 2) = 0.5 * lam * kap * sqr(s) + al * s / root2 + m3hatsq;
  mPpr(1, 3) = v2 / root2 * (al - root2 * lam * kap * s - lam * mupr);
  mPpr(2, 2) = (0.5 * lam * kap * sqr(s) + al * s / root2
                + m3hatsq) / tanb;
  mPpr(2, 3) = v1 / root2 * (al - root2 * lam * kap * s - lam * mupr);
  mPpr(3, 3) = (al + lam * kap * s /root2
             + 3.0 * lam * kap * s / root2
             + lam * mupr) * v1 * v2 / (root2 * s) 
     - 3.0 * ak * s / root2 - lam * smu * sqr(vev) / (root2 * s)
     - 2.0 * mSpsq - kap * mupr * s / root2 
     - xiF * (4.0 * kap + root2 * mupr / s) - root2 * xiS / s;

  mPpr.symmetrise();

  DoubleMatrix mP(3, 3);
  // LCT: Rotate CP-odd mass^2 matrix into (G, A, S_I) basis
  mP = rot3d(beta).transpose() * mPpr * rot3d(beta);
  /// LCT: Drop Goldstone from 3 x 3 CP-odd Higgs mass^2 matrix and
  /// construct 2 x 2 matrix in (A, S_I) basis
  mP2(1, 1) = mP(2, 2);
  mP2(1, 2) = mP(2, 3);
  mP2(2, 1) = mP(3, 2);
  mP2(2, 2) = mP(3, 3);

  mHpmsq = mP2(1, 1) + mw2 - 0.5 * sqr(vev) * sqr(lam);
}


// PA: fills tree level CP even and CP odd Higgs mass matrices 
// and tree level mHPm all without treel level EWSB conditions used.
// Same as treeHiggs at tree level *if* tree level EWSB is imposed.
// Useful when s is set to zero and for testing EWSB conditions.
// CP even mass matrix in (Hd, Hu, S) basis  
// CP odd Higgs mass matrices mPpr in (Hd, Hu, S) basis 
// and mP2 in rotated basis (A, S) -- goldstone boson removed      
void NmssmSoftsusy::treeHiggsAlt(DoubleMatrix & mS, DoubleMatrix & mPpr, 
                               DoubleMatrix & mP2, double & mHpmsq, 
                               double beta) const {
  double tanb = displayTanb();
  double sinb = sin(beta), cosb = cos(beta);
  double sinb2 = sqr(sinb), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta);
  double mz = displayMzRun(), mz2 = sqr(mz);
  double mw = displayMwRun(), mw2 = sqr(mw);
  double mHu2 = displayMh2Squared(), mHd2 = displayMh1Squared();
	/// LCT: NMSSM parameters
  double l = displayLambda(), k = displayKappa();
  double mupr = displayMupr(), smu = displaySusyMu();
  double al = displayTrialambda(), ak = displayTriakappa();
  double xiF = displayXiF(), xiS = displayXiS();
  double m3sq = displayM3Squared(), mSpr2 = displayMspSquared();
  double mS2 = displayMsSquared();

  double vev = displayHvev(), v1 = vev * cosb, v2 = vev * sinb;
  double s = displaySvev();

  double mueff = l * s / root2 + smu;
  double m3hatsq = m3sq + l * (mupr * s / root2 + xiF);
 
  /// CP-even Higgs in EHT notation and basis (HdR, HuR, SR)
  mS(1, 1) = mHd2 + sqr(l * s / root2 + smu) + 0.5 * sqr(l * v2)
     + 0.5 * mz2 * cos2b + mz2 * sinb2;
  mS(1, 2) = - sinb * cosb * mz2 + sqr(l) * sqr(vev) * sinb * cosb
     - (0.5 * l * k * sqr(s) + al * s / root2 + m3hatsq);
  mS(1, 3) = root2 * l * mueff * v1 + root2 * l * smu * v1
     - (al + l * k * s * root2 + l * mupr) * v2 / root2;
  mS(2, 2) = mHu2 + sqr(l * s / root2 + smu) + 0.5 * sqr(l * v1) 
     - 0.5 * mz2 * cos2b + mz2 * cosb2;
  mS(2, 3) = root2 * l * mueff * v2 + root2 * l * smu * v2 
     - (al + l * k * s * root2 + l * mupr) * v1 / root2;
  mS(3, 3) = mS2 + mSpr2 + sqr(mupr) + 2.0 * k * xiF + root2 * ak * s
     + 3.0 * sqr(k * s) + 3.0 * root2 * k * s * mupr + 0.5 * sqr(l * vev)
     - l * k * v1 * v2; 

  mS.symmetrise();

  mPpr(1, 1) = mHd2 + sqr(smu + l * s / root2) + sqr(l * v2) 
     + 0.5 * mz2 * cos2b;
  mPpr(1, 2) = 0.5 * l * k * sqr(s) + al * s / root2 + m3hatsq;
  mPpr(1, 3) = v2 / root2 * (al - root2 * l * k * s - l * mupr);
  mPpr(2, 2) = mHu2 + sqr(smu + l * s / root2) + sqr(l * v1) 
     - 0.5 * mz2 * cos2b;
  mPpr(2, 3) = v1 / root2 * (al - root2 * l * k * s - l * mupr);
  mPpr(3, 3) = mS2 - mSpr2 + sqr(mupr) + sqr(k * s) + 0.5 * sqr(l * vev)
     - 2.0 * k * xiF - root2 * ak * s + root2 * k * s * mupr + l * k * v1 * v2;


  mPpr.symmetrise();
  DoubleMatrix mP(3, 3);
  // LCT: Rotate CP-odd mass^2 matrix into (G, A, S_I) basis
  mP = rot3d(beta).transpose() * mPpr * rot3d(beta);
  /// LCT: Drop Goldstone from 3 x 3 CP-odd Higgs mass^2 matrix and
  /// construct 2 x 2 matrix in (A, S_I) basis
  mP2(1, 1) = mP(2, 2);
  mP2(1, 2) = mP(2, 3);
  mP2(2, 1) = mP(3, 2);
  mP2(2, 2) = mP(3, 3);
  mHpmsq = mP2(1, 1) + mw2 - 0.5 * sqr(vev) * sqr(l);
   }


/// LCT: new routine to set tree-level NMSSM CP-even and odd Higgs mass matrices squared
void NmssmSoftsusy::calcDrBarHiggs(double beta, double mz2, double mw2, double sinthDRbar, drBarPars & eg) {
  
  /// LCT: Copied condition on top mass from softsusy.cpp
   if (eg.mt > 200. || eg.mt < 50.) {
    /// Gone badly off-track
    flagProblemThrown(true);
    if (eg.mt > 200.) eg.mt = 200.;
    if (eg.mt < 50.) eg.mt = 50.;
  }
   
   //PA: initialise CP even mass matrix in (Hd, Hu, S) basis
  // CP odd Higgs mass matrices mPpr in (Hd, Hu, S) basis 
   //and mP2 in roatated basis (A, S) -- goldstone boson removed    
   DoubleMatrix mS(3,3), mPpr(3,3), mP2(2,2); 
   double mHpmsq; 
 //PA: fill tree level CP even and CP odd Higgs mass matrices 
  //and tree level mHPm. 
   treeHiggs(mS, mPpr, mP2, mHpmsq, beta);
	/// LCT: Diagonalise 3 x 3 CP-even Higgs mass matrix.
	/// Mass basis (H1, H2, H3) is ordered in increasing mass mH1 < mH2 < mH3
  DoubleVector mhsq(3);
  DoubleMatrix mixh(3,3);
  if (mS.diagonaliseSym(mixh, mhsq) > TOLERANCE *
      1.0e-3) { 
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
 
  /// LCT: Diagonalise
  //PA using thetaH for now since in nmssm this is not used by CP even
  DoubleVector mSq = mP2.sym2by2(eg.thetaH);

  if (mSq(1) < 0. || mSq(2) < 0.) {
     flagTachyon(A0);
     if (PRINTOUT > 1) cout << " mA1/mA2 tachyon";
  }
  DoubleVector temp(mSq.apply(ccbSqrt));
  if (temp(1) > temp(2)) eg.thetaH = eg.thetaH + PI * 0.5;

  int pos;
  eg.mA0(1) = temp.min(pos); eg.mA0(2) = temp.max();
  
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

DoubleMatrix NmssmSoftsusy::addStopHiggs(double p, double mt, DoubleMatrix & higgs) {
  const drBarPars& forLoops = displayDrBarPars(); ///< Contains DRbar tree-level masses
  double    sinthDrbar  = calcSinthdrbar(), sinthDrbar2 = sqr(sinthDrbar);
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double tanthDrbar2 = sqr(tan(asin(sinthDrbar)));
  double    g	    = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    beta    = atan(displayTanb());
  double mw = displayMwRun();
  double mz = displayMzRun();

  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  double svev = displaySvev();
  double lam = displayLambda();
  double mueff = -(displaySusyMu() + lam * svev / root2); ///<< LCT: Note sign!
  double q = displayMu(), 
    ht = forLoops.ht,
    hb = forLoops.hb,
    htsq = sqr(ht), 
    sinb = sin(beta), cosb = cos(beta), 
    hbsq = sqr(hb),
    v1 = displayHvev() * cos(beta),
    v2 = displayHvev() * sin(beta);

  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb;
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Quadrilinear contributions
  ///
  /// LCT: stop-stop-Hp-Hm couplings in (L R) basis.  Order (G+ H+).
  DoubleVector lChHChHstopLstopL(2), lChHChHstopRstopR(2);
  for (int i=1; i<=2; i++) {
    lChHChHstopLstopL(i) = (sqr(g) * (3.0 + tanthDrbar2) 
			    * (sqr(C(i, 1)) - sqr(C(i, 2))) 
			    - 12.0 *  sqr(C(i, 1)) * hbsq) / 12.0;
    lChHChHstopRstopR(i) = - (4.0 * sqr(g) * tanthDrbar2 
			      * (sqr(C(i, 1)) - sqr(C(i, 2))) 
			      + 12.0 *  sqr(C(i, 2)) * htsq) / 12.0;
  }
  /// LCT: Charged contribution
  for (int i=1; i<=2; i++) { 
    higgs(1, 1) += - lChHChHstopLstopL(i) * a0(higgsc(i), q);
    higgs(2, 2) += - lChHChHstopRstopR(i) * a0(higgsc(i), q);
  }

  /// LCT: stop-stop-A0-A0 and stop-stop-h0-h0 couplings in (L R) basis.  
  /// Order (G0 A1 A2)
  double u1ll = - sqr(g) * (3.0 - tanthDrbar2) / 12.0;
  double u1rr = - sqr(g) * tanthDrbar2 / 3.0;
  double u2ll = - htsq;
  double u2rr = u2ll;
  double u3lr = 0.5 * lam * ht;
  DoubleVector lAAstopLstopL(3), lAAstopRstopR(3), lAAstopLstopR(3), lHHstopLstopL(3), lHHstopRstopR(3), lHHstopLstopR(3);
  for (int i=1; i<=3; i++) {
    lAAstopLstopL(i) = u1ll * (sqr(P(i, 1)) - sqr(P(i, 2))) + u2ll * sqr(P(i, 2));
    lAAstopRstopR(i) = u1rr * (sqr(P(i, 1)) - sqr(P(i, 2))) + u2rr * sqr(P(i, 2));
    lAAstopLstopR(i) = - 2.0 * u3lr * P(i, 1) * P(i, 3);

    lHHstopLstopL(i) = u1ll * (sqr(S(i, 1)) - sqr(S(i, 2))) + u2ll * sqr(S(i, 2));
    lHHstopRstopR(i) = u1rr * (sqr(S(i, 1)) - sqr(S(i, 2))) + u2rr * sqr(S(i, 2));
    lHHstopLstopR(i) = 2.0 * u3lr * S(i, 1) * S(i, 3);
  }


  /// LCT: CP-odd and CP-even contribution
 for (int i=1; i<=3; i++) { 
   higgs(1, 1) += - 0.5 * (lAAstopLstopL(i) * a0(higgsa(i), q) 
			   + lHHstopLstopL(i) * a0(higgsm(i), q));
   higgs(2, 2) += - 0.5 * (lAAstopRstopR(i) * a0(higgsa(i), q) 
			   + lHHstopRstopR(i) * a0(higgsm(i), q));
   higgs(1, 2) += - 0.5 * (lAAstopLstopR(i) * a0(higgsa(i), q) 
			   + lHHstopLstopR(i) * a0(higgsm(i), q));
  }

 /// LCT: Trilinear contributions
 ///
 /// LCT: Charged Higgs Feynman rules
  /// (H+ G+, L R) basis
  DoubleMatrix lChHstopLsbotLR(2, 2), lChHstopRsbotLR(2, 2);
  DoubleVector temp(2), temp2(2);
  lChHstopLsbotLR(1, 1) = (g * displayMwRun() * sin(2.0 * beta) 
    - htsq * v2 * cosb - hbsq * v1 * sinb) / root2;
  lChHstopLsbotLR(1, 2) = (mueff * hb * cosb - 
			   forLoops.ub * sinb);
  lChHstopLsbotLR(2, 1) = (-g * displayMwRun() * cos(2.0 * beta) 
    - htsq * v2 * sinb + hbsq * v1 * cosb) / root2;
  lChHstopLsbotLR(2, 2) = hb * mueff * sinb + forLoops.ub * cosb;

  /// (H+ G+, L R) basis
  lChHstopRsbotLR(1, 1) = ht * mueff * sinb - forLoops.ut * cosb;
  lChHstopRsbotLR(1, 2) = ht * hb * (- v1 * cosb - v2 * sinb) / root2;
  lChHstopRsbotLR(2, 1) = -ht * mueff * cosb - forLoops.ut * sinb;

  /// LCT: Flip sign to match Slavich conventions on Hpm
  lChHstopLsbotLR = - lChHstopLsbotLR;
  lChHstopRsbotLR = - lChHstopRsbotLR;

  /// LCT: Rotate sfermions to (1, 2) mass basis
  DoubleMatrix lChHstopLsbot12(2, 2);
  temp(1) = lChHstopLsbotLR(1, 1);
  temp(2) = lChHstopLsbotLR(1, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopLsbot12(1, 1) = temp2(1);
  lChHstopLsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopLsbotLR(2, 1);
  temp(2) = lChHstopLsbotLR(2, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopLsbot12(2, 1) = temp2(1);
  lChHstopLsbot12(2, 2) = temp2(2);

  /// LCT: Rotate sfermions to (1, 2) mass basis
  DoubleMatrix lChHstopRsbot12(2, 2);
  temp(1) = lChHstopRsbotLR(1, 1);
  temp(2) = lChHstopRsbotLR(1, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopRsbot12(1, 1) = temp2(1);
  lChHstopRsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopRsbotLR(2, 1);
  temp(2) = lChHstopRsbotLR(2, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopRsbot12(2, 1) = temp2(1);
  lChHstopRsbot12(2, 2) = temp2(2);

  /// LCT: Rearrange masses to account for change in mixing convention
  DoubleVector Hpm(2);
  Hpm(1) = higgsc(2);
  Hpm(2) = higgsc(1); 

  /// LCT: Charged Higgs contribution
  for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      double b0p = b0(p, Hpm(i), msbot(j), q); 
      higgs(1, 1) += sqr(lChHstopLsbot12(i, j)) * b0p;
      higgs(1, 2) += 
  	lChHstopLsbot12(i, j) * lChHstopRsbot12(i, j) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lChHstopRsbot12(i, j)) * b0p;
    }

  /// LCT: A0-stop-stop couplings in (L R) basis
  double u8lr = forLoops.ut / root2;
  double u9lr = 0.5 * lam * ht;
  DoubleVector lAstopLstopR(3), lAstopRstopL(3);
  for (int i=1; i<=3; i++) {
    lAstopLstopR(i) = u8lr * P(i, 2) + u9lr * (v1 * P(i, 3) + svev * P(i, 1));
    lAstopRstopL(i) = - lAstopLstopR(i);
  }
  /// LCT: Rotate to mixed basis
  DoubleMatrix lAstopLstop12(3, 2), lAstopRstop12(3, 2);
  for (int i=1; i<=3; i++) {
    lAstopLstop12(i, 1) = lAstopLstopR(i) * st;
    lAstopLstop12(i, 2) = lAstopLstopR(i) * ct;
    lAstopRstop12(i, 1) = lAstopRstopL(i) * ct;
    lAstopRstop12(i, 2) = - lAstopRstopL(i) * st;
    }

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=2; j++) {
      higgs(1, 1) += sqr(lAstopLstop12(i, j)) * b0(p, higgsa(i), mstop(j), q);
      higgs(1, 2) += lAstopLstop12(i, j) * lAstopRstop12(i, j) 
	           * b0(p, higgsa(i), mstop(j), q);
      higgs(2, 2) +=  sqr(lAstopRstop12(i, j)) * b0(p, higgsa(i), mstop(j), q);
    }
  }

  /// LCT: (CP-even) Higgs-stop-stop couplings in (L R) basis.  Order (s1 s2 s3)
  //
  /// Taken from Franke & Fraas, arXiv:hep-ph/9512366
  DoubleMatrix lsStopLStopLR(3, 2), lsStopLStop12(3, 2);
  DoubleMatrix lsStopRStopLR(3, 2), lsStopRStop12(3, 2);
  double ut = forLoops.ut;  

  lsStopLStopLR(1, 1) = g * mz / costhDrbar * guL 
  * (S(1, 2) * sinb - S(1, 1) * cosb) - g * sqr(mt) / (mw * sinb) * S(1, 2);
  lsStopLStopLR(1, 2) = - ut / root2 * S(1, 2) 
    + 0.5 * g * mt / (root2 * mw * sinb) * lam 
    * (v1 * S(1, 3) + svev * S(1, 1));
  lsStopLStopLR(2, 1) = g * mz / costhDrbar * guL 
  * (S(2, 2) * sinb - S(2, 1) * cosb) - g * sqr(mt) / (mw * sinb) * S(2, 2);
  lsStopLStopLR(2, 2) = - ut / root2 * S(2, 2) 
    + 0.5 * g * mt / (root2 * mw * sinb) * lam 
    * (v1 * S(2, 3) + svev * S(2, 1));
  lsStopLStopLR(3, 1) = g * mz / costhDrbar * guL 
  * (S(3, 2) * sinb - S(3, 1) * cosb) - g * sqr(mt) / (mw * sinb) * S(3, 2);
  lsStopLStopLR(3, 2) = - ut / root2 * S(3, 2) 
    + 0.5 * g * mt / (root2 * mw * sinb) * lam 
    * (v1 * S(3, 3) + svev * S(3, 1));
  
  lsStopRStopLR(1, 1) = lsStopLStopLR(1, 2);
  lsStopRStopLR(1, 2) = - g * sqr(mt) / (mw * sinb) * S(1, 2) 
    + 2.0 / 3.0 * g * mw * tanthDrbar2 * (S(1, 2) * sinb - S(1, 1) * cosb);
  lsStopRStopLR(2, 1) = lsStopLStopLR(2, 2);
  lsStopRStopLR(2, 2) = - g * sqr(mt) / (mw * sinb) * S(2, 2) 
    + 2.0 / 3.0 * g * mw * tanthDrbar2 * (S(2, 2) * sinb - S(2, 1) * cosb);
  lsStopRStopLR(3, 1) = lsStopLStopLR(3, 2);
  lsStopRStopLR(3, 2) = - g * sqr(mt) / (mw * sinb) * S(3, 2) 
    + 2.0 / 3.0 * g * mw * tanthDrbar2 * (S(3, 2) * sinb - S(3, 1) * cosb);
  
  /// Mix stops up
  for (int i=1; i<=3; i++) {
    temp(1) = lsStopLStopLR(i, 1);
    temp(2) = lsStopLStopLR(i, 2);
    temp2 = rot2d(thetat) * temp;
    lsStopLStop12(i, 1) = temp2(1);
    lsStopLStop12(i, 2) = temp2(2);
    temp(1) = lsStopRStopLR(i, 1);
    temp(2) = lsStopRStopLR(i, 2);
    temp2 = rot2d(thetat) * temp;
    lsStopRStop12(i, 1) = temp2(1);
    lsStopRStop12(i, 2) = temp2(2);
  }  

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=2; j++) {
      higgs(1, 1) += sqr(lsStopLStop12(i, j)) * b0(p, higgsm(i), mstop(j), q);
      higgs(1, 2) += lsStopLStop12(i, j) * lsStopRStop12(i, j) 
      * b0(p, higgsm(i), mstop(j), q);     
      higgs(2, 2) +=  sqr(lsStopRStop12(i, j)) * b0(p, higgsm(i), mstop(j), q);
    }
  }

  higgs(2, 1) = higgs(1, 2);
  return higgs;
}


void NmssmSoftsusy::getNeutralinoFermionStopCoup(DoubleVector & aPsi0TStopr, DoubleVector & bPsi0TStopr, DoubleVector & aPsi0TStopl, DoubleVector & bPsi0TStopl) const {
  double gp = sqrt(0.6) * displayGaugeCoupling(1);
  double g = displayGaugeCoupling(2);
  double lam = displayLambda();

  /// Call MSSM couplings
  Softsusy<SoftParsNmssm>::getNeutralinoFermionStopCoup(aPsi0TStopr, bPsi0TStopr, aPsi0TStopl, bPsi0TStopl);

  /// Fill remaining values
}

DoubleMatrix NmssmSoftsusy::addStopNeutralino(double p, double mt, DoubleMatrix & neutralino) {
  const drBarPars& forLoops = displayDrBarPars();
  double q = displayMu();

  DoubleVector aPsi0TStopr(5), bPsi0TStopr(5), aPsi0TStopl(5), bPsi0TStopl(5);
  getNeutralinoFermionStopCoup(aPsi0TStopr, bPsi0TStopr, aPsi0TStopl, bPsi0TStopl);

  ComplexVector aChi0TStopl(5), bChi0TStopl(5), aChi0TStopr(5), bChi0TStopr(5);
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 

  aChi0TStopl = n.complexConjugate() * aPsi0TStopl;
  bChi0TStopl = n * bPsi0TStopl;
  aChi0TStopr = n.complexConjugate() * aPsi0TStopr;
  bChi0TStopr = n * bPsi0TStopr;

  DoubleVector gChi0TopStopLL(5), fChi0TopStopLL(5);
  DoubleVector gChi0TopStopLR(5), fChi0TopStopLR(5);
  DoubleVector gChi0TopStopRR(5), fChi0TopStopRR(5);

  int i;
  for (i=1; i<=5; i++) {
    fChi0TopStopLL(i) = (aChi0TStopl(i) * aChi0TStopl(i).conj() + 
      bChi0TStopl(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLL(i) = (bChi0TStopl(i).conj() * aChi0TStopl(i) + 
      bChi0TStopl(i) * aChi0TStopl(i).conj()).real();
    fChi0TopStopRR(i) = (aChi0TStopr(i) * aChi0TStopr(i).conj() + 
      bChi0TStopr(i) * bChi0TStopr(i).conj()).real();
    gChi0TopStopRR(i) = (bChi0TStopr(i).conj() * aChi0TStopr(i) + 
      bChi0TStopr(i) * aChi0TStopr(i).conj()).real();
    fChi0TopStopLR(i) = (aChi0TStopr(i) * aChi0TStopl(i).conj() + 
      bChi0TStopr(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLR(i) = (bChi0TStopl(i).conj() * aChi0TStopr(i) + 
      bChi0TStopr(i) * aChi0TStopl(i).conj()).real();
  }

  for (i=1; i<=5; i++) {
    double one = gfn(p, mneut(i), mt, q);
    double two = 2.0 * mneut(i) * mt * b0(p, mneut(i), mt, q);
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0TopStopLL(i) * one - gChi0TopStopLL(i) * two;
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0TopStopRR(i) * one - gChi0TopStopRR(i) * two;
    neutralino(1, 2) = neutralino(1, 2) +
      fChi0TopStopLR(i) * one - gChi0TopStopLR(i) * two;
  }

  return neutralino;
}


void NmssmSoftsusy::addStopCorrection(double p, DoubleMatrix & mass, double mt) {
/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    flagTachyon(stop);
    if (mass(1, 1) < 0.) mass(1, 1) = EPSTOL;
    else mass(2, 2) = EPSTOL;
    return;
  }
  /// one-loop correction matrix
  DoubleMatrix piSq(2, 2); /// Self-energy matrix

  /// Corrections themselves start here
  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2), 
    higgs(2, 2), electroweak(2, 2), chargino(2, 2), neutralino(2, 2);

  /// LCT: Corrections from strong interactions
  addStopQCD(p, mt, strong);
  /// LCT: Corrections from stops
  addStopStop(p, mt, stop);
  /// LCT: Corrections from sbottoms
  addStopSbottom(p, mt, sbottom);
  /// LCT: Corrections from Higgs
  addStopHiggs(p, mt, higgs);
  /// LCT: Electroweak corrections
  addStopEweak(p, electroweak);
  /// LCT: Chargino contributions
  addStopChargino(p, chargino);
  /// LCT: Neutralino contribution
  addStopNeutralino(p, mt, neutralino);

  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (strong + stop + sbottom + electroweak + higgs + chargino + neutralino);

  piSq(2, 1) = piSq(1, 2);
  
  mass = mass - piSq;
}

DoubleMatrix NmssmSoftsusy::addSupHiggs(double p1, double p2, int family, DoubleMatrix & higgs) {
  const drBarPars& forLoops = displayDrBarPars();
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double tanthDrbar2 = sqr(tan(asin(sinthDrbar)));
  double    g	    = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    beta    = atan(displayTanb());
  double mw = displayMwRun();
  double mz = displayMzRun();
  double svev = displaySvev();
  double q = displayMu(), 
    ht = displayDrBarPars().ht,
    hb = displayDrBarPars().hb,
    htsq = sqr(ht), // 
    sinb = sin(beta), cosb = cos(beta), 
    hbsq = sqr(hb),
    v1 = displayHvev() * cos(beta),
    v2 = displayHvev() * sin(beta);
  DoubleVector msd(2);
  msd(1)          = forLoops.md(1, family);
  msd(2)          = forLoops.md(2, family);
  DoubleVector msup(2);
  msup(1)          = forLoops.mu(1, family);
  msup(2)          = forLoops.mu(2, family);

  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb;
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: sup-sup-Hp-Hm quadrilinear couplings in (L R) basis
  DoubleVector lChHChHstopLstopL(2), lChHChHstopRstopR(2);
  for (int i=1; i<=2; i++) {
    lChHChHstopLstopL(i) = (sqr(g) * (3.0 + tanthDrbar2) 
		   * (sqr(C(i, 1)) - sqr(C(i, 2)))) / 12.0;
    lChHChHstopRstopR(i) = - (4.0 * sqr(g) * tanthDrbar2 
		     * (sqr(C(i, 1)) - sqr(C(i, 2)))) / 12.0;
  }
  /// LCT: Charged Higgs contribution
  for (int i=1; i<=2; i++) { 
    double a0p = a0(higgsc(i), q);
    higgs(1, 1) += - lChHChHstopLstopL(i) * a0(higgsc(i), q);
    higgs(2, 2) += - lChHChHstopRstopR(i) * a0(higgsc(i), q);
  }

  /// LCT: stop-stop-A0-A0 and stop-stop-h0-h0 couplings in (L R) basis.  
  /// Order (G0 A1 A2)
  double u1ll = - sqr(g) * (3.0 - tanthDrbar2) / 12.0;
  double u1rr = - sqr(g) * tanthDrbar2 / 3.0;
  double u2ll = 0.0; ///<< LCT: Third family approx
  double u2rr = u2ll;
  DoubleVector lAAstopLstopL(3), lAAstopRstopR(3), lHHstopLstopL(3), lHHstopRstopR(3);
  for (int i=1; i<=3; i++) {
    lAAstopLstopL(i) = u1ll * (sqr(P(i, 1)) - sqr(P(i, 2))) + u2ll * sqr(P(i, 2));
    lAAstopRstopR(i) = u1rr * (sqr(P(i, 1)) - sqr(P(i, 2))) + u2rr * sqr(P(i, 2));
 
    lHHstopLstopL(i) = u1ll * (sqr(S(i, 1)) - sqr(S(i, 2))) + u2ll * sqr(S(i, 2));
    lHHstopRstopR(i) = u1rr * (sqr(S(i, 1)) - sqr(S(i, 2))) + u2rr * sqr(S(i, 2));
  }
  /// LCT: CP-odd and CP-even contribution
 for (int i=1; i<=3; i++) { 
   higgs(1, 1) += - 0.5 * (lAAstopLstopL(i) * a0(higgsa(i), q) 
			   + lHHstopLstopL(i) * a0(higgsm(i), q));
   higgs(2, 2) += - 0.5 * (lAAstopRstopR(i) * a0(higgsa(i), q) 
			   + lHHstopRstopR(i) * a0(higgsm(i), q));
  }

 /// LCT: Trilinear contributions
 ///
 /// LCT: Charged Higgs Feynman rules
 /// Charged Higgs Feynman rules
  DoubleVector temp(2), temp2(2);
  DoubleMatrix lChHstopLsbotLR(2, 2); /// (H+ G+, L R) basis
  lChHstopLsbotLR(1, 1) = g * displayMwRun() * sin(2.0 * beta) / root2;
  lChHstopLsbotLR(2, 1) =-g * displayMwRun() * cos(2.0 * beta) / root2;

  /// LCT: Flip sign to match Slavich conventions on Hpm
  lChHstopLsbotLR = - lChHstopLsbotLR;

  DoubleMatrix lChHstopLsbot12(2, 2);
  temp(1) = lChHstopLsbotLR(1, 1);
  temp(2) = lChHstopLsbotLR(1, 2);
  temp2 = temp;
  lChHstopLsbot12(1, 1) = temp2(1);
  lChHstopLsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopLsbotLR(2, 1);
  temp(2) = lChHstopLsbotLR(2, 2);
  temp2 = temp;
  lChHstopLsbot12(2, 1) = temp2(1);
  lChHstopLsbot12(2, 2) = temp2(2);

  DoubleMatrix lChHstopRsbotLR(2, 2); /// (H+ G+, L R) basis
  /// LCT: Flip sign to match Slavich conventions on Hpm
  lChHstopRsbotLR = - lChHstopRsbotLR;

  DoubleMatrix lChHstopRsbot12(2, 2);
  temp(1) = lChHstopRsbotLR(1, 1);
  temp(2) = lChHstopRsbotLR(1, 2);
  temp2 = temp;
  lChHstopRsbot12(1, 1) = temp2(1);
  lChHstopRsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopRsbotLR(2, 1);
  temp(2) = lChHstopRsbotLR(2, 2);
  temp2 = temp;
  lChHstopRsbot12(2, 1) = temp2(1);
  lChHstopRsbot12(2, 2) = temp2(2);

  /// LCT: Rearrange masses to account for change in mixing convention
  DoubleVector Hpm(2);
  Hpm(1) = higgsc(2);
  Hpm(2) = higgsc(1); 

  /// LCT: Charged Higgs contribution
 for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      double b0p = b0(p1, Hpm(i), msd(j), q); 
      higgs(1, 1) += sqr(lChHstopLsbot12(i, j)) * b0p;
      b0p = b0(p2, Hpm(i), msd(j), q); 
      higgs(2, 2) += sqr(lChHstopRsbot12(i, j)) * b0p;
    }

 /// LCT: A0-sup-sup couplings in (L R) basis
 double u8lr = 0.0; ///<< LCT: Third family approximation
  double u9lr = 0.0;
  DoubleVector lAstopLstopR(3), lAstopRstopL(3);
  for (int i=1; i<=3; i++) {
    lAstopLstopR(i) = u8lr * P(i, 2) + u9lr * (v1 * P(i, 3) + svev * P(i, 1));
    lAstopRstopL(i) = - lAstopLstopR(i);
  }
  /// LCT: Rotate to mixed basis
  DoubleMatrix lAstopLstop12(3, 2), lAstopRstop12(3, 2);
  for (int i=1; i<=3; i++) {
    lAstopLstop12(i, 1) = lAstopLstopR(i);
    lAstopLstop12(i, 2) = lAstopLstopR(i);
    lAstopRstop12(i, 1) = lAstopRstopL(i);
    lAstopRstop12(i, 2) = - lAstopRstopL(i);
    }
 for(int i=1; i<=3; i++)
    for (int j=1; j<=2; j++) {
      double b0p = b0(p1, higgsa(i), msup(j), q); 
      higgs(1, 1) += sqr(lAstopLstop12(i, j)) * b0p;
      b0p = b0(p2, higgsa(i), msup(j), q); 
      higgs(2, 2) += sqr(lAstopRstop12(i, j)) * b0p;
    }
    
  /// LCT: (CP-even) Higgs-stop-stop couplings in (L R) basis.  Order (s1 s2 s3)
  //
  /// Taken from Franke & Fraas, arXiv:hep-ph/9512366
  DoubleMatrix lsStopLStopLR(3, 2), lsStopLStop12(3, 2);
  DoubleMatrix lsStopRStopLR(3, 2), lsStopRStop12(3, 2);  
  
  lsStopLStopLR(1, 1) = g * mz / costhDrbar * guL 
  * (S(1, 2) * sinb - S(1, 1) * cosb);
  lsStopLStopLR(2, 1) = g * mz / costhDrbar * guL 
  * (S(2, 2) * sinb - S(2, 1) * cosb);
  lsStopLStopLR(3, 1) = g * mz / costhDrbar * guL 
  * (S(3, 2) * sinb - S(3, 1) * cosb);
  
  lsStopRStopLR(1, 1) = lsStopLStopLR(1, 2);
  lsStopRStopLR(1, 2) = 2.0 / 3.0 * g * mw * tanthDrbar2 
  * (S(1, 2) * sinb - S(1, 1) * cosb);
  lsStopRStopLR(2, 1) = lsStopLStopLR(2, 2);
  lsStopRStopLR(2, 2) = 2.0 / 3.0 * g * mw * tanthDrbar2 
  * (S(2, 2) * sinb - S(2, 1) * cosb);
  lsStopRStopLR(3, 1) = lsStopLStopLR(3, 2);
  lsStopRStopLR(3, 2) = 2.0 / 3.0 * g * mw * tanthDrbar2 
  * (S(3, 2) * sinb - S(3, 1) * cosb);
  
  /// Mix stops up
  for (int i=1; i<=3; i++) {
    temp(1) = lsStopLStopLR(i, 1);
    temp(2) = lsStopLStopLR(i, 2);
    temp2 = temp;
    lsStopLStop12(i, 1) = temp2(1);
    lsStopLStop12(i, 2) = temp2(2);
    temp(1) = lsStopRStopLR(i, 1);
    temp(2) = lsStopRStopLR(i, 2);
    temp2 = temp;
    lsStopRStop12(i, 1) = temp2(1);
    lsStopRStop12(i, 2) = temp2(2);
  }      
  for(int i=1; i<=3; i++)
    for (int j=1; j<=2; j++) {
      double b0p = b0(p1, higgsm(i), msup(j), q);
      higgs(1, 1) += sqr(lsStopLStop12(i, j)) * b0p;
      b0p = b0(p2, higgsm(i), msup(j), q);
      higgs(2, 2) += sqr(lsStopRStop12(i, j)) * b0p;
    }

  return higgs;
}

DoubleMatrix NmssmSoftsusy::addSupNeutralino(double p1, double p2, int family, DoubleMatrix & neutralino) {
  const drBarPars& forLoops = displayDrBarPars();
  double q = displayMu();

  DoubleVector aPsi0TStopr(5), bPsi0TStopr(5), aPsi0TStopl(5), bPsi0TStopl(5);
  getNeutralinoFermionSupCoup(aPsi0TStopr, bPsi0TStopr, aPsi0TStopl, bPsi0TStopl);

  ComplexVector aChi0TStopl(5), bChi0TStopl(5), aChi0TStopr(5), bChi0TStopr(5);
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 

  aChi0TStopl = n.complexConjugate() * aPsi0TStopl;
  bChi0TStopl = n * bPsi0TStopl;
  aChi0TStopr = n.complexConjugate() * aPsi0TStopr;
  bChi0TStopr = n * bPsi0TStopr;

  DoubleVector gChi0TopStopLL(5), fChi0TopStopLL(5);
  DoubleVector gChi0TopStopLR(5), fChi0TopStopLR(5);
  DoubleVector gChi0TopStopRR(5), fChi0TopStopRR(5);

  int i;
  for (i=1; i<=5; i++) {
    fChi0TopStopLL(i) = (aChi0TStopl(i) * aChi0TStopl(i).conj() + 
      bChi0TStopl(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLL(i) = (bChi0TStopl(i).conj() * aChi0TStopl(i) + 
      bChi0TStopl(i) * aChi0TStopl(i).conj()).real();
    fChi0TopStopRR(i) = (aChi0TStopr(i) * aChi0TStopr(i).conj() + 
      bChi0TStopr(i) * bChi0TStopr(i).conj()).real();
    gChi0TopStopRR(i) = (bChi0TStopr(i).conj() * aChi0TStopr(i) + 
      bChi0TStopr(i) * aChi0TStopr(i).conj()).real();
    fChi0TopStopLR(i) = (aChi0TStopr(i) * aChi0TStopl(i).conj() + 
      bChi0TStopr(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLR(i) = (bChi0TStopl(i).conj() * aChi0TStopr(i) + 
      bChi0TStopr(i) * aChi0TStopl(i).conj()).real();
  }

  for (i=1; i<=5; i++) {
    double one = gfn(p1, mneut(i), 0., q);
    double two = 0.;
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0TopStopLL(i) * one - gChi0TopStopLL(i) * two;
    one = gfn(p2, mneut(i), 0., q);
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0TopStopRR(i) * one - gChi0TopStopRR(i) * two;
  }

  return neutralino;
}

void NmssmSoftsusy::addSupCorrection(DoubleMatrix & mass, int family) {

/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    if (mass(1, 1) < 0.) mass(1, 1) = EPSTOL;
    else mass(2, 2) = EPSTOL;
    if (family == 1) flagTachyon(sup);
    else if (family == 2) flagTachyon(scharm);
    return;
  }
  /// one-loop correction matrix
  DoubleMatrix piSq(2, 2); /// Self-energy matrix
  DoubleVector msup(2);
  msup(1) = displayDrBarPars().mu(1, family);
  msup(2) = displayDrBarPars().mu(2, family);

  DoubleMatrix strong(2, 2),higgs(2, 2), electroweak(2, 2), 
  chargino(2, 2), neutralino(2, 2);
  double p1 = msup(1), p2 = msup(2);

  /// LCT: QCD contribution
  addSupQCD(p1, p2, family, strong);
  /// LCT: Higgs contribution
  addSupHiggs(p1, p2, family, higgs);
  /// LCT: Electroweak contribution
  addSupEweak(p1, p2, family, electroweak);
   /// LCT: Chargino contribution
  addSupChargino(p1, p2, family, chargino);
  /// LCT: Neutralino contribution
  addSupNeutralino(p1, p2, family, neutralino);

  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (strong + higgs + electroweak + chargino + neutralino);

  piSq(2, 1) = piSq(1, 2);
  mass = mass - piSq;
  
}

DoubleMatrix NmssmSoftsusy::addSbotHiggs(double p, double mt, DoubleMatrix & higgs) {
  const drBarPars& forLoops = displayDrBarPars();
  double    mz      = displayMzRun();
  double    mw      = displayMwRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double tanthDrbar2 = sqr(tan(asin(sinthDrbar)));
  double    alpha   = forLoops.thetaH;
  double    g	    = displayGaugeCoupling(2), gsq = sqr(g);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    e       = g * sinthDrbar;
  double    beta    = atan(displayTanb());

  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double    thetatau= forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double    ctau    = cos(thetatau);
  double    stau    = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1)          = forLoops.msnu(1);
  msnu(2)          = forLoops.msnu(2);
  msnu(3)          = forLoops.msnu(3);

  double    mg      = forLoops.mGluino;
  double    smu     = -displaySusyMu();
  double q = displayMu(), g3sq = sqr(displayGaugeCoupling(3)), 
    ht = forLoops.ht,
    hb = forLoops.hb,
    mb = forLoops.mb,
    htau = forLoops.htau,
    htsq = sqr(ht), 
    sinb = sin(beta), cosb = cos(beta), 
    hbsq = sqr(hb);
  double lam = displayLambda();
  double svev = displaySvev();
  double mueff = smu - lam * svev / root2;
  double  v1 = displayHvev() * cosb, v2 = displayHvev() * sinb;

  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb;
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Quadriliner contributions
  /// LCT: sbottom-sbottom-Hp-Hm couplings in (L R) basis.  Order (G+ H+).
  DoubleVector lChHChHsbotLsbotL(2), lChHChHsbotRsbotR(2);
  for (int i=1; i<=2; i++) {
    lChHChHsbotLsbotL(i) = (sqr(C(i, 1)) * sqr(g) * (tanthDrbar2 - 3.0) 
			    - sqr(C(i, 2)) * (sqr(g) * (tanthDrbar2 - 3.0) 
					      + 12.0 * htsq)) / 12.0;
    lChHChHsbotRsbotR(i) = (sqr(C(i, 1)) * (2.0 * sqr(gp) - 12.0 * hbsq) 
			    - sqr(C(i, 2)) * 2.0 * sqr(gp)) / 12.0;
  }
  /// LCT: Charged contribution
  for (int i=1; i<=2; i++) { 
    higgs(1, 1) += - lChHChHsbotLsbotL(i) * a0(higgsc(i), q);
    higgs(2, 2) += - lChHChHsbotRsbotR(i) * a0(higgsc(i), q);
  }

  /// LCT: sbottom-sbottom-A0-A0 couplings in (L R) basis.  Order (G0 A1 A2).
  DoubleVector lAAsbotLsbotL(3), lAAsbotLsbotR(3), lAAsbotRsbotR(3);
  for (int i=1; i<=3; i++) {
    lAAsbotLsbotL(i) = gsq / 12.0 * (3.0 + tanthDrbar2) 
      * (sqr(P(i, 1)) - sqr(P(i, 2))) - hbsq * sqr(P(i, 1));
    lAAsbotRsbotR(i) = gsq * tanthDrbar2 / 6.0 * (sqr(P(i, 1)) - sqr(P(i, 2))) 
      - hbsq * sqr(P(i, 1));
    lAAsbotLsbotR(i) = - lam * hb * P(i, 2) * P(i, 3);
  } 
  /// LCT: sbottom-sbottom-h0-h0 couplings in (L R) basis.  Order (h1 h2 h3)
  DoubleVector lHHsbotLsbotL(3), lHHsbotRsbotR(3), lHHsbotLsbotR(3);
  for (int i=1; i<=3; i++) {
    lHHsbotLsbotL(i) = gsq / 12.0 * (3.0 + tanthDrbar2) 
      * (sqr(S(i, 1)) - sqr(S(i, 2))) - hbsq * sqr(S(i, 1));
    lHHsbotRsbotR(i) = gsq * tanthDrbar2 / 6.0 * (sqr(S(i, 1)) - sqr(S(i, 2))) 
      - hbsq * sqr(S(i, 1));
    lHHsbotLsbotR(i) = lam * hb * S(i, 2) * S(i, 3);
  }

  for (int i=1; i<=3; i++) {
    higgs(1, 1) += - 0.5 * (lAAsbotLsbotL(i) * a0(higgsa(i), q) 
			    + lHHsbotLsbotL(i) * a0(higgsm(i), q));
    higgs(2, 2) += - 0.5 * (lAAsbotRsbotR(i) * a0(higgsa(i), q) 
			    + lHHsbotRsbotR(i) * a0(higgsm(i), q));
    higgs(1, 2) += - 0.5 * (lAAsbotLsbotR(i) * a0(higgsa(i), q) 
			    + lHHsbotLsbotR(i) * a0(higgsm(i), q));
  }

  /// Charged Higgs Feynman rules
  DoubleMatrix lChHsbotLstopLR(2, 2); /// (H+ G+, L R) basis
  lChHsbotLstopLR(2, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2 
    - ht * mt * sinb + hb * mb * cosb;
  lChHsbotLstopLR(2, 2) = -ht * mueff * cosb - forLoops.ut * sinb;
  lChHsbotLstopLR(1, 1) = g * displayMwRun() * sin(2.0 * beta) / root2 
    - ht * mt * cosb - hb * mb * sinb;
  lChHsbotLstopLR(1, 2) = ht * mueff * sinb - forLoops.ut * cosb;

  DoubleVector temp(2), temp2(2);
  DoubleMatrix lChHsbotLstop12(2, 2);
  temp(1) = lChHsbotLstopLR(1, 1);
  temp(2) = lChHsbotLstopLR(1, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotLstop12(1, 1) = temp2(1);
  lChHsbotLstop12(1, 2) = temp2(2);
  temp(1) = lChHsbotLstopLR(2, 1);
  temp(2) = lChHsbotLstopLR(2, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotLstop12(2, 1) = temp2(1);
  lChHsbotLstop12(2, 2) = temp2(2);

  DoubleMatrix lChHsbotRstopLR(2, 2); /// (H+ G+, L R) basis
  lChHsbotRstopLR(1, 1) = hb * mueff * cosb - forLoops.ub * sinb;
  lChHsbotRstopLR(1, 2) = -ht * mb * cosb - hb * mt * sinb;
  lChHsbotRstopLR(2, 1) = hb * mueff * sinb + forLoops.ub * cosb;
  DoubleMatrix lChHsbotRstop12(2, 2);
  temp(1) = lChHsbotRstopLR(1, 1);
  temp(2) = lChHsbotRstopLR(1, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotRstop12(1, 1) = temp2(1);
  lChHsbotRstop12(1, 2) = temp2(2);
  temp(1) = lChHsbotRstopLR(2, 1);
  temp(2) = lChHsbotRstopLR(2, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotRstop12(2, 1) = temp2(1);
  lChHsbotRstop12(2, 2) = temp2(2);

  /// LCT: Flip sign for consistency with Slavich conventions on Hpm mixing
  /// matrix C(i, j)
  lChHsbotLstop12 = - lChHsbotLstop12;
  lChHsbotRstop12 = - lChHsbotRstop12;

 /// LCT: Rearrange masses to account for change in mixing convention
  DoubleVector Hpm(2);
  Hpm(1) = higgsc(2);
  Hpm(2) = higgsc(1); 

  /// LCT: Charged Higgs contribution
  for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      double b0p = b0(p, Hpm(i), mstop(j), q); 
      higgs(1, 1) += sqr(lChHsbotLstop12(i, j)) * b0p;
      higgs(2, 2) += sqr(lChHsbotRstop12(i, j)) * b0p;
      higgs(1, 2) += lChHsbotLstop12(i, j) * lChHsbotRstop12(i, j) * b0p;
    }

  /// LCT: A0-sbottom-sbottom couplings in (L R) basis
  DoubleVector lAsbotLsbotR(3), lAsbotRsbotL(3);
  for (int i=1; i<=3; i++) {
    lAsbotLsbotR(i) = forLoops.ub / root2 * P(i, 1) 
      + 0.5 * lam * hb * (v2 * P(i, 3) + svev * P(i, 2));
    lAsbotRsbotL(i) = - lAsbotLsbotR(i);
  }
  /// LCT: Rotate to mixed basis
  DoubleMatrix lAsbotLsbot12(3, 2), lAsbotRsbot12(3, 2);
  for (int i=1; i<=3; i++) {
    lAsbotLsbot12(i, 1) = lAsbotLsbotR(i) * sb;
    lAsbotLsbot12(i, 2) = lAsbotLsbotR(i) * cb;
    lAsbotRsbot12(i, 1) = lAsbotRsbotL(i) * cb;
    lAsbotRsbot12(i, 2) = - lAsbotRsbotL(i) * sb;
    }
  for (int i=1; i<=3; i++) {
    for (int j=1; j<=2; j++) {
      higgs(1, 1) += sqr(lAsbotLsbot12(i, j)) * b0(p, higgsa(i), msbot(j), q);
      higgs(1, 2) += lAsbotLsbot12(i, j) * lAsbotRsbot12(i, j) 
	* b0(p, higgsa(i), msbot(j), q);
      higgs(2, 2) += sqr(lAsbotRsbot12(i, j)) * b0(p, higgsa(i), msbot(j), q);
    }
  }

  /// LCT: (CP-even) Higgs-sbot-sbot couplings in (L R) basis.  Order (s1 s2 s3)
  //
  /// Taken from Franke & Fraas, arXiv:hep-ph/9512366
  DoubleMatrix lsSbotLSbotLR(3, 2), lsSbotLSbot12(3, 2);
  DoubleMatrix lsSbotRSbotLR(3, 2), lsSbotRSbot12(3, 2);
  double ub = forLoops.ub;  
  
  lsSbotLSbotLR(1, 1) = g * mz / costhDrbar * gdL 
    * (S(1, 2) * sinb - S(1, 1) * cosb) - g * sqr(mb) / (mw * cosb) * S(1, 1);
  lsSbotLSbotLR(1, 2) = - ub / root2 * S(1, 1) 
    + 0.5 * g * mb / (root2 * mw * cosb) * lam * (v2 * S(1, 3) + svev * S(1, 2));
  lsSbotLSbotLR(2, 1) = g * mz / costhDrbar * gdL 
    * (S(2, 2) * sinb - S(2, 1) * cosb) - g * sqr(mb) / (mw * cosb) * S(2, 1);
  lsSbotLSbotLR(2, 2) = - ub / root2 * S(2, 1) 
    + 0.5 * g * mb / (root2 * mw * cosb) * lam * (v2 * S(2, 3) + svev * S(2, 2));
  lsSbotLSbotLR(3, 1) = g * mz / costhDrbar * gdL 
    * (S(3, 2) * sinb - S(3, 1) * cosb) - g * sqr(mb) / (mw * cosb) * S(3, 1);
  lsSbotLSbotLR(3, 2) = - ub / root2 * S(3, 1) 
    + 0.5 * g * mb / (root2 * mw * cosb) * lam * (v2 * S(3, 3) + svev * S(3, 2));
  
  lsSbotRSbotLR(1, 1) = lsSbotLSbotLR(1, 2);
  lsSbotRSbotLR(1, 2) = - g * sqr(mb) / (mw * cosb) * S(1, 1) 
    - 1.0 / 3.0 * g * mw * tanthDrbar2 * (S(1, 2) * sinb - S(1, 1) * cosb);
  lsSbotRSbotLR(2, 1) = lsSbotLSbotLR(2, 2);
  lsSbotRSbotLR(2, 2) = - g * sqr(mb) / (mw * cosb) * S(2, 1) 
    - 1.0 / 3.0 * g * mw * tanthDrbar2 * (S(2, 2) * sinb - S(2, 1) * cosb);
  lsSbotRSbotLR(3, 1) = lsSbotLSbotLR(3, 2);
  lsSbotRSbotLR(3, 2) = - g * sqr(mb) / (mw * cosb) * S(3, 1) 
    - 1.0 / 3.0 * g * mw * tanthDrbar2 * (S(3, 2) * sinb - S(3, 1) * cosb);
  
  /// Mix sbots up
  for (int i=1; i<=3; i++) {
    temp(1) = lsSbotLSbotLR(i, 1);
    temp(2) = lsSbotLSbotLR(i, 2);
    temp2 = rot2d(thetab) * temp;
    lsSbotLSbot12(i, 1) = temp2(1);
    lsSbotLSbot12(i, 2) = temp2(2);
    temp(1) = lsSbotRSbotLR(i, 1);
    temp(2) = lsSbotRSbotLR(i, 2);
    temp2 = rot2d(thetab) * temp;
    lsSbotRSbot12(i, 1) = temp2(1);
    lsSbotRSbot12(i, 2) = temp2(2);
  }  

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=2; j++) {
      higgs(1, 1) += sqr(lsSbotLSbot12(i, j)) * b0(p, higgsm(i), msbot(j), q);
      higgs(2, 2) += sqr(lsSbotRSbot12(i, j)) * b0(p, higgsm(i), msbot(j), q);
      higgs(1, 2) += lsSbotLSbot12(i, j) * lsSbotRSbot12(i, j) 
       * b0(p, higgsm(i), msbot(j), q);
    }
  }

  higgs(2, 1) = higgs(1, 2);

  return higgs;
}

DoubleMatrix NmssmSoftsusy::addSbotNeutralino(double p, double mt, DoubleMatrix & neutralino) {
  const drBarPars& forLoops = displayDrBarPars();
  double q = displayMu(),  mb = forLoops.mb;
  double g = displayGaugeCoupling(2);
  double gp = displayGaugeCoupling(1) * sqrt(0.6);
  double hb = forLoops.hb;
  DoubleVector aPsi0BSbotr(5), bPsi0BSbotr(5), aPsi0BSbotl(5),
    bPsi0BSbotl(5);
 /// Neutralino Feynman rules 
  aPsi0BSbotr(1) = ydR * gp / root2;
  bPsi0BSbotl(1) = gp / (3.0 * root2);
  bPsi0BSbotl(2) = -g / root2;
  aPsi0BSbotl(3) = hb;
  bPsi0BSbotr(3) = hb;

 ComplexVector aChi0BSbotl(5), bChi0BSbotl(5), aChi0BSbotr(5), bChi0BSbotr(5);
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);

  aChi0BSbotl = n.complexConjugate() * aPsi0BSbotl;
  bChi0BSbotl = n * bPsi0BSbotl;
  aChi0BSbotr = n.complexConjugate() * aPsi0BSbotr;
  bChi0BSbotr = n * bPsi0BSbotr;

  DoubleVector gChi0BotSbotLL(5), fChi0BotSbotLL(5);
  DoubleVector gChi0BotSbotLR(5), fChi0BotSbotLR(5);
  DoubleVector gChi0BotSbotRR(5), fChi0BotSbotRR(5);

  int i;
  for (i=1; i<=5; i++) {
    fChi0BotSbotLL(i) = (aChi0BSbotl(i) * aChi0BSbotl(i).conj() + 
      bChi0BSbotl(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLL(i) = (bChi0BSbotl(i).conj() * aChi0BSbotl(i) + 
      bChi0BSbotl(i) * aChi0BSbotl(i).conj()).real();
    fChi0BotSbotRR(i) = (aChi0BSbotr(i) * aChi0BSbotr(i).conj() + 
      bChi0BSbotr(i) * bChi0BSbotr(i).conj()).real();
    gChi0BotSbotRR(i) = (bChi0BSbotr(i).conj() * aChi0BSbotr(i) + 
      bChi0BSbotr(i) * aChi0BSbotr(i).conj()).real();
    fChi0BotSbotLR(i) = (aChi0BSbotr(i) * aChi0BSbotl(i).conj() + 
      bChi0BSbotr(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLR(i) = (bChi0BSbotl(i).conj() * aChi0BSbotr(i) + 
      bChi0BSbotr(i) * aChi0BSbotl(i).conj()).real();
  }

  for (i=1; i<=5; i++) {
    double one = gfn(p, mneut(i), mb, q);
    double two = 2.0 * mneut(i) * mb * b0(p, mneut(i), mb, q);
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0BotSbotLL(i) * one - gChi0BotSbotLL(i) * two;
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0BotSbotRR(i) * one - gChi0BotSbotRR(i) * two;
    neutralino(1, 2) = neutralino(1, 2) +
      fChi0BotSbotLR(i) * one - gChi0BotSbotLR(i) * two;
  }

  return neutralino;
}

void NmssmSoftsusy::addSbotCorrection(double p, DoubleMatrix & mass, double mt) {
/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    flagTachyon(stop);
    if (mass(1, 1) < 0.) mass(1, 1) = EPSTOL;
    else mass(2, 2) = EPSTOL;
    return;
  }

  DoubleMatrix piSq(2, 2); /// Self-energy matrix

  /// Corrections themselves start here
  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2), 
    higgs(2, 2), electroweak(2, 2), chargino(2, 2), neutralino(2, 2);

  /// LCT: Corrections from strong interactions
  addSbotQCD(p, mt, strong);
  /// LCT: Corrections from sfermions
  addSbotSfermion(p, mt, stop, sbottom);
  /// LCT: Corrections from Higgs
  addSbotHiggs(p, mt, higgs);
  /// LCT: Electroweak corrections
  addSbotEweak(p, electroweak);
  // /// LCT: Chargino contributions
  addSbotChargino(p, mt, chargino);
  // /// LCT: Neutralino contribution
  addSbotNeutralino(p, mt, neutralino);

  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (strong + stop + sbottom + electroweak + higgs + chargino + neutralino);

  piSq(2, 1) = piSq(1, 2);

  mass = mass - piSq;
}

DoubleMatrix NmssmSoftsusy::addSdownHiggs(double p1, double p2, int family, DoubleMatrix & higgs) {
  const drBarPars& forLoops = displayDrBarPars();
  double mz      = displayMzRun();
  double mw      = displayMwRun();
  double sinthDrbar = calcSinthdrbar();
  double costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double tanthDrbar2 = sqr(tan(asin(sinthDrbar)));
  double g	    = displayGaugeCoupling(2), gsq = sqr(g);
  double gp      = displayGaugeCoupling(1) * sqrt(0.6), gpsq = sqr(gp);
  double beta    = atan(displayTanb()), sinb = sin(beta), cosb = cos(beta);

  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  DoubleVector mstop(2);
  mstop(1) = forLoops.mu(1, 3);
  mstop(2) = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1) = forLoops.me(1, 3);
  mstau(2) = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1) = forLoops.msnu(1);
  msnu(2) = forLoops.msnu(2);
  msnu(3) = forLoops.msnu(3);
  double smu = -displaySusyMu(); ///<< LCT: Note sign!
  double q = displayMu();
  DoubleVector msup(2);
  msup(1)          = forLoops.mu(1, family);
  msup(2)          = forLoops.mu(2, family);
  DoubleVector msd(2);
  msd(1)          = forLoops.md(1, family);
  msd(2)          = forLoops.md(2, family);

  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb;
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Quadriliner contributions
  /// LCT: sbottom-sbottom-Hp-Hm couplings in (L R) basis.  Order (G+ H+).
  DoubleVector lChHChHsbotLsbotL(2), lChHChHsbotRsbotR(2);
  for (int i=1; i<=2; i++) {
    lChHChHsbotLsbotL(i) = gsq * (tanthDrbar2 - 3.0) / 12.0 
      * (sqr(C(i, 1)) - sqr(C(i, 2)));
    lChHChHsbotRsbotR(i) = gpsq * (sqr(C(i, 1)) - sqr(C(i, 2))) / 6.0;
  }
  /// LCT: Charged contribution
  for (int i=1; i<=2; i++) { 
    higgs(1, 1) += - lChHChHsbotLsbotL(i) * a0(higgsc(i), q);
    higgs(2, 2) += - lChHChHsbotRsbotR(i) * a0(higgsc(i), q);
  }

  /// LCT: sbottom-sbottom-A0-A0 couplings in (L R) basis.  Order (G0 A1 A2).
  DoubleVector lAAsbotLsbotL(3), lAAsbotLsbotR(3), lAAsbotRsbotR(3);
  for (int i=1; i<=3; i++) {
    lAAsbotLsbotL(i) = gsq / 12.0 * (3.0 + tanthDrbar2) 
      * (sqr(P(i, 1)) - sqr(P(i, 2)));
    lAAsbotRsbotR(i) = gsq * tanthDrbar2 / 6.0 * (sqr(P(i, 1)) - sqr(P(i, 2)));
  } 
  /// LCT: sbottom-sbottom-h0-h0 couplings in (L R) basis.  Order (h1 h2 h3)
  DoubleVector lHHsbotLsbotL(3), lHHsbotRsbotR(3), lHHsbotLsbotR(3);
  for (int i=1; i<=3; i++) {
    lHHsbotLsbotL(i) = gsq / 12.0 * (3.0 + tanthDrbar2) 
      * (sqr(S(i, 1)) - sqr(S(i, 2)));
    lHHsbotRsbotR(i) = gsq * tanthDrbar2 / 6.0 * (sqr(S(i, 1)) - sqr(S(i, 2)));
  }

  for (int i=1; i<=3; i++) {
    higgs(1, 1) += - 0.5 * (lAAsbotLsbotL(i) * a0(higgsa(i), q) 
			    + lHHsbotLsbotL(i) * a0(higgsm(i), q));
    higgs(2, 2) += - 0.5 * (lAAsbotRsbotR(i) * a0(higgsa(i), q) 
			    + lHHsbotRsbotR(i) * a0(higgsm(i), q));
  }

 /// Charged Higgs Feynman rules
  DoubleMatrix lChHsbotLstopLR(2, 2); /// (H+ G+, L R) basis
  lChHsbotLstopLR(2, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2;
  lChHsbotLstopLR(1, 1) = g * displayMwRun() * sin(2.0 * beta) / root2;

  DoubleVector temp(2), temp2(2);
  DoubleMatrix lChHsbotLstop12(2, 2);
  temp(1) = lChHsbotLstopLR(1, 1);
  temp(2) = lChHsbotLstopLR(1, 2);
  temp2 = temp;
  lChHsbotLstop12(1, 1) = temp2(1);
  lChHsbotLstop12(1, 2) = temp2(2);
  temp(1) = lChHsbotLstopLR(2, 1);
  temp(2) = lChHsbotLstopLR(2, 2);
  temp2 = temp;
  lChHsbotLstop12(2, 1) = temp2(1);
  lChHsbotLstop12(2, 2) = temp2(2);

  DoubleMatrix lChHsbotRstopLR(2, 2); /// (H+ G+, L R) basis
  DoubleMatrix lChHsbotRstop12(2, 2);
  temp(1) = lChHsbotRstopLR(1, 1);
  temp(2) = lChHsbotRstopLR(1, 2);
  temp2 = temp;
  lChHsbotRstop12(1, 1) = temp2(1);
  lChHsbotRstop12(1, 2) = temp2(2);
  temp(1) = lChHsbotRstopLR(2, 1);
  temp(2) = lChHsbotRstopLR(2, 2);
  temp2 = temp;
  lChHsbotRstop12(2, 1) = temp2(1);
  lChHsbotRstop12(2, 2) = temp2(2);

  /// LCT: Flip sign for consistency with Slavich conventions on Hpm mixing
  /// matrix C(i, j)
  lChHsbotLstop12 = - lChHsbotLstop12;
  lChHsbotRstop12 = - lChHsbotRstop12;

 /// LCT: Rearrange masses to account for change in mixing convention
  DoubleVector Hpm(2);
  Hpm(1) = higgsc(2);
  Hpm(2) = higgsc(1); 

  /// LCT: Charged Higgs contribution
 for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      double b0p = b0(p1, Hpm(i), msup(j), q); 
      higgs(1, 1) += sqr(lChHsbotLstop12(i, j)) * b0p;
      b0p = b0(p2, Hpm(i), msup(j), q); 
      higgs(2, 2) += sqr(lChHsbotRstop12(i, j)) * b0p;
    }

 /// LCT: No sdown-sdown-A0 couplings in third family approximation

 /// LCT: (CP-even) sdown-sdown-h0 couplings in (L R) basis.  Order (s1 s2 s3)
  //
  /// Taken from Franke & Fraas, arXiv:hep-ph/9512366
  DoubleMatrix lsSbotLSbotLR(3, 2), lsSbotLSbot12(3, 2);
  DoubleMatrix lsSbotRSbotLR(3, 2), lsSbotRSbot12(3, 2);
  
  lsSbotLSbotLR(1, 1) = g * mz / costhDrbar * gdL 
    * (S(1, 2) * sinb - S(1, 1) * cosb);
  lsSbotLSbotLR(2, 1) = g * mz / costhDrbar * gdL 
    * (S(2, 2) * sinb - S(2, 1) * cosb);
  lsSbotLSbotLR(3, 1) = g * mz / costhDrbar * gdL 
    * (S(3, 2) * sinb - S(3, 1) * cosb);
  
  lsSbotRSbotLR(1, 1) = lsSbotLSbotLR(1, 2);
  lsSbotRSbotLR(1, 2) = - 1.0 / 3.0 * g * mw * tanthDrbar2 
    * (S(1, 2) * sinb - S(1, 1) * cosb);
  lsSbotRSbotLR(2, 1) = lsSbotLSbotLR(2, 2);
  lsSbotRSbotLR(2, 2) = - 1.0 / 3.0 * g * mw * tanthDrbar2 
    * (S(2, 2) * sinb - S(2, 1) * cosb);
  lsSbotRSbotLR(3, 1) = lsSbotLSbotLR(3, 2);
  lsSbotRSbotLR(3, 2) = - 1.0 / 3.0 * g * mw * tanthDrbar2 
    * (S(3, 2) * sinb - S(3, 1) * cosb);
  
  /// Mix sbots up
  for (int i=1; i<=3; i++) {
    temp(1) = lsSbotLSbotLR(i, 1);
    temp(2) = lsSbotLSbotLR(i, 2);
    temp2 = temp;
    lsSbotLSbot12(i, 1) = temp2(1);
    lsSbotLSbot12(i, 2) = temp2(2);
    temp(1) = lsSbotRSbotLR(i, 1);
    temp(2) = lsSbotRSbotLR(i, 2);
    temp2 = temp;
    lsSbotRSbot12(i, 1) = temp2(1);
    lsSbotRSbot12(i, 2) = temp2(2);
  }  

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=2; j++) {
      double b0p = b0(p1, higgsm(i), msd(j), q);
      higgs(1, 1) += sqr(lsSbotLSbot12(i, j)) * b0p;
      b0p = b0(p2, higgsm(i), msd(j), q);
      higgs(2, 2) += sqr(lsSbotRSbot12(i, j)) * b0p;
    }
  }

  return higgs;
}

DoubleMatrix NmssmSoftsusy::addSdownNeutralino(double p1, double p2, int family, DoubleMatrix & neutralino) {
  const drBarPars& forLoops = displayDrBarPars();
  double g = displayGaugeCoupling(2);
  double gp = displayGaugeCoupling(1) * sqrt(0.6);
  double q = displayMu();

  /// Neutralino Feynman rules
  DoubleVector aPsi0BSbotr(5), bPsi0BSbotr(5), aPsi0BSbotl(5),
    bPsi0BSbotl(5); 
  aPsi0BSbotr(1) = gp * ydR / root2;
  bPsi0BSbotl(1) = gp * ydL / root2;
  bPsi0BSbotl(2) = -0.5 * g * root2;

  ComplexVector aChi0BSbotl(5), bChi0BSbotl(5), aChi0BSbotr(5),
    bChi0BSbotr(5);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);

  aChi0BSbotl = n.complexConjugate() * aPsi0BSbotl;
  bChi0BSbotl = n * bPsi0BSbotl;
  aChi0BSbotr = n.complexConjugate() * aPsi0BSbotr;
  bChi0BSbotr = n * bPsi0BSbotr;

  DoubleVector gChi0BotSbotLL(5), fChi0BotSbotLL(5);
  DoubleVector gChi0BotSbotLR(5), fChi0BotSbotLR(5);
  DoubleVector gChi0BotSbotRR(5), fChi0BotSbotRR(5);
  for (int i=1; i<=5; i++) {
    fChi0BotSbotLL(i) = (aChi0BSbotl(i) * aChi0BSbotl(i).conj() + 
      bChi0BSbotl(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLL(i) = (bChi0BSbotl(i).conj() * aChi0BSbotl(i) + 
      bChi0BSbotl(i) * aChi0BSbotl(i).conj()).real();
    fChi0BotSbotRR(i) = (aChi0BSbotr(i) * aChi0BSbotr(i).conj() + 
      bChi0BSbotr(i) * bChi0BSbotr(i).conj()).real();
    gChi0BotSbotRR(i) = (bChi0BSbotr(i).conj() * aChi0BSbotr(i) + 
      bChi0BSbotr(i) * aChi0BSbotr(i).conj()).real();
    fChi0BotSbotLR(i) = (aChi0BSbotr(i) * aChi0BSbotl(i).conj() + 
      bChi0BSbotr(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLR(i) = (bChi0BSbotl(i).conj() * aChi0BSbotr(i) + 
      bChi0BSbotr(i) * aChi0BSbotl(i).conj()).real();
  }

  for (int i=1; i<=5; i++) {
    double one = gfn(p1, mneut(i), 0., q);
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0BotSbotLL(i) * one;
    one = gfn(p2, mneut(i), 0., q);
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0BotSbotRR(i) * one;
  }

  return neutralino;
}

void NmssmSoftsusy::addSdownCorrection(DoubleMatrix & mass, int family) {

/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    if (family == 1) flagTachyon(sdown); 
    else if (family == 2) flagTachyon(sstrange);
    return;
  }

  const drBarPars& forLoops = displayDrBarPars();
  DoubleMatrix piSq(2, 2); /// Self-energy matrix
  DoubleVector msd(2);
  msd(1) = forLoops.md(1, family);
  msd(2) = forLoops.md(2, family);
  double    p1 = msd(1), p2 = msd(2);
  ///  p1 = p2 = sqrt(msd(1) * msd(2)); 

  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2), higgs(2, 2), 
    chargino(2, 2), neutralino(2, 2), electroweak(2, 2);

  /// LCT: QCD contribution
  addSdownQCD(p1, p2, family, strong);
  /// LCT: Higgs contribution
  addSdownHiggs(p1, p2, family, higgs);
  /// LCT: Electroweak contribution
  addSdownEweak(p1, p2, family, electroweak);
  /// LCT: Chargino contribution
  addSdownChargino(p1, p2, family, chargino);
  /// LCT: Neutralino contribution
  addSdownNeutralino(p1, p2, family, neutralino);
 
  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (strong + higgs + chargino + neutralino + electroweak);

  mass = mass - piSq;	
}

double NmssmSoftsusy::addSnuTauHiggs(double p, double & higgs) {
  const drBarPars& forLoops = displayDrBarPars(); ///< Contains DRbar tree-level masses
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    g	    = displayGaugeCoupling(2), gsq = sqr(g);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    beta    = atan(displayTanb());
  double    htau    = forLoops.htau;
  double mtau = htau * displayHvev() / root2 * cos(beta);
  DoubleVector msnu(3);
  msnu(1)           = forLoops.msnu(1);
  msnu(2)           = forLoops.msnu(2);
  msnu(3)           = forLoops.msnu(3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  double   thetatau = forLoops.thetatau;
  double    smu     = -displaySusyMu();
  double q = displayMu(), 
    htausq = sqr(htau), 
    sinb = sin(beta), cosb = cos(beta), 
    v1 = displayHvev() * cos(beta), v2 = displayHvev() * sin(beta),
    mz = displayMzRun();
  double tanthDrbar2 = sqr(tan(asin(sinthDrbar)));
  double svev = displaySvev();
  double lam = displayLambda();
  double mueff = -(displaySusyMu() + lam * svev / root2); ///<< LCT: Note sign!

  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb;
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Quadrilinear contributions
  ///
  /// LCT: sneutrino-sneutrino-Hp-Hm couplings in (L R) basis.  Order (G+ H+).
  DoubleVector lChHChHsnuLsnuL(2);
  for (int i=1; i<=2; i++) {
    lChHChHsnuLsnuL(i) = (sqr(g) * (1.0 - tanthDrbar2) 
			  * (sqr(C(i, 1)) - sqr(C(i, 2))) - 4.0 *  sqr(C(i, 1)) * htausq) / 4.0;
  }
  /// LCT: Charged contribution
  for (int i=1; i<=2; i++) { 
    higgs += - lChHChHsnuLsnuL(i) * a0(higgsc(i), q);
  }

  /// LCT: sneutrino-sneutrino-A0-A0 and sneutrino-sneutrino-h0-h0 couplings 
  /// in (L R) basis. Order (G0 A1 A2)
  DoubleVector lAAsnuLsnuL(3), lHHsnuLsnuL(3);
  for (int i=1; i<=3; i++) {
    lAAsnuLsnuL(i) = - 0.25 * gsq * (1.0 + tanthDrbar2) * (sqr(P(i, 1)) - sqr(P(i, 2)));
    lHHsnuLsnuL(i) = - 0.25 * gsq * (1.0 + tanthDrbar2) *  (sqr(S(i, 1)) - sqr(S(i, 2)));
  }
  /// LCT: CP-odd and CP-even contribution
 for (int i=1; i<=3; i++) { 
   higgs += - 0.5 * (lAAsnuLsnuL(i) * a0(higgsa(i), q) 
		     + lHHsnuLsnuL(i) * a0(higgsm(i), q));
  }

  /// LCT: Trilinear contributions
  ///
  /// LCT: Charged Higgs Feynman rules
  /// (H+ G+, L R) basis
  DoubleMatrix lChHsnuLstauLR(2, 2); 
  lChHsnuLstauLR(1, 1) = (g * displayMwRun() * sin(2.0 * beta) / root2
    - htau * mtau * sinb);
  lChHsnuLstauLR(1, 2) = (mueff * htau * cosb - forLoops.utau * sinb);
  lChHsnuLstauLR(2, 1) = (-g * displayMwRun() * cos(2.0 * beta) 
    + htausq * v1 * cosb) / root2;
  lChHsnuLstauLR(2, 2) = htau * mueff * sinb + forLoops.utau * cosb;

  DoubleVector temp(2), temp2(2);
  DoubleMatrix lChHsnuLstau12(2, 2);
  temp(1) = lChHsnuLstauLR(1, 1);
  temp(2) = lChHsnuLstauLR(1, 2);
  temp2 = rot2d(thetatau) * temp;
  lChHsnuLstau12(1, 1) = temp2(1);
  lChHsnuLstau12(1, 2) = temp2(2);
  temp(1) = lChHsnuLstauLR(2, 1);
  temp(2) = lChHsnuLstauLR(2, 2);
  temp2 = rot2d(thetatau) * temp;
  lChHsnuLstau12(2, 1) = temp2(1);
  lChHsnuLstau12(2, 2) = temp2(2);

  /// LCT: Flip sign to match Slavich conventions on Hpm
  lChHsnuLstau12 = - lChHsnuLstau12;

  /// LCT: Rearrange masses to account for change in mixing convention
  DoubleVector Hpm(2);
  Hpm(1) = higgsc(2);
  Hpm(2) = higgsc(1); 

  /// LCT: Charged Higgs contribution
  for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      double b0p = b0(p, Hpm(i), mstau(j), q); 
      higgs += sqr(lChHsnuLstau12(i, j)) * b0p;
    }

  /// LCT: (CP-even) Higgs-snuetrino-sneutrino couplings in (L R) basis.  
  /// Order (s1 s2 s3)
  DoubleVector lsSnuLSnuLR(3), lsSnuLSnu12(3);

  for (int i=1; i<=3; i++) {
    lsSnuLSnuLR(i) = - 0.25 * gsq * (1.0 + tanthDrbar2) 
      * (v1 * S(i, 1) - v2 * S(i, 2));
  }
  

  for (int i=1; i<=3; i++) {
      higgs += sqr(lsSnuLSnuLR(i)) * b0(p, higgsm(i), msnu(3), q);
  }

  return higgs;
}

void NmssmSoftsusy::addSnuTauCorrection(double & mass) {

 /// No point adding radiative corrections to tachyonic particles
  if (mass < 0.0) { 
    flagTachyon(snutau);
    mass = EPSTOL;
    return;
  }
  double p = sqrt(mass);
  /// one-loop correction matrix
  double piSq; /// Self-energy matrix
	
  double stop = 0., sbottom = 0., higgs = 0., electroweak = 0.,
    chargino = 0., neutralino = 0.;

  /// LCT: Sfermion contribution
  addSnuTauSfermion(p, stop, sbottom);
  /// LCT: Higgs contribution
  addSnuTauHiggs(p, higgs);
  /// LCT: Electroweak contribution
  addSnuTauEweak(p, electroweak);
  /// LCT: Gaugino contributions    
  addSnuTauGaugino(p, chargino, neutralino);

  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (stop + sbottom + higgs + electroweak + chargino + neutralino);

  mass = mass - piSq;
}

double NmssmSoftsusy::addSnuHiggs(double p, int family, double & higgs) {
  const drBarPars& forLoops = displayDrBarPars(); ///< Contains DRbar tree-level masses
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    alpha   = forLoops.thetaH;
  double    g	    = displayGaugeCoupling(2), gsq = sqr(g);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    beta    = atan(displayTanb());

  DoubleVector msnu(3);
  msnu(1)           = forLoops.msnu(1);
  msnu(2)           = forLoops.msnu(2);
  msnu(3)           = forLoops.msnu(3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double   thetatau = forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double  ctau      = cos(thetatau);
  double  stau      = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  double q = displayMu(), 
    sinb = sin(beta), cosb = cos(beta), 
    mz = displayMzRun(); 
   double tanthDrbar2 = sqr(tan(asin(sinthDrbar)));
   double   v1 = displayHvev() * cos(beta), v2 = displayHvev() * sin(beta);

  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb;
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Quadrilinear contributions
  ///
  /// LCT: sneutrino-sneutrino-Hp-Hm couplings in (L R) basis.  Order (G+ H+).
  DoubleVector lChHChHsnuLsnuL(2);
  for (int i=1; i<=2; i++) {
    lChHChHsnuLsnuL(i) = 0.25 * gsq * (1.0 - tanthDrbar2) 
      * (sqr(C(i, 1)) - sqr(C(i, 2)));
  }
  /// LCT: Charged contribution
  for (int i=1; i<=2; i++) { 
    higgs += - lChHChHsnuLsnuL(i) * a0(higgsc(i), q);
  }

  /// LCT: sneutrino-sneutrino-A0-A0 and sneutrino-sneutrino-h0-h0 couplings 
  /// in (L R) basis. Order (G0 A1 A2)
  DoubleVector lAAsnuLsnuL(3), lHHsnuLsnuL(3);
  for (int i=1; i<=3; i++) {
    lAAsnuLsnuL(i) = - 0.25 * gsq * (1.0 + tanthDrbar2) * (sqr(P(i, 1)) - sqr(P(i, 2)));
    lHHsnuLsnuL(i) = - 0.25 * gsq * (1.0 + tanthDrbar2) *  (sqr(S(i, 1)) - sqr(S(i, 2)));
  }
  /// LCT: CP-odd and CP-even contribution
 for (int i=1; i<=3; i++) { 
   higgs += - 0.5 * (lAAsnuLsnuL(i) * a0(higgsa(i), q) 
		     + lHHsnuLsnuL(i) * a0(higgsm(i), q));
  }

  /// LCT: Trilinear contributions
  ///
  /// LCT: Charged Higgs Feynman rules
  /// (H+ G+, L R) basis
  DoubleMatrix lChHsnuLstauLR(2, 2); 
  lChHsnuLstauLR(2, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2;
  lChHsnuLstauLR(1, 1) = (g * displayMwRun() * sin(2.0 * beta) / root2);

  DoubleVector temp(2), temp2(2);
  DoubleMatrix lChHsnuLstau12(2, 2);
  temp(1) = lChHsnuLstauLR(1, 1);
  temp(2) = lChHsnuLstauLR(1, 2);
  temp2 = temp;
  lChHsnuLstau12(1, 1) = temp2(1);
  lChHsnuLstau12(1, 2) = temp2(2);
  temp(1) = lChHsnuLstauLR(2, 1);
  temp(2) = lChHsnuLstauLR(2, 2);
  temp2 = temp;
  lChHsnuLstau12(2, 1) = temp2(1);
  lChHsnuLstau12(2, 2) = temp2(2);

  /// LCT: Flip sign to match Slavich conventions on Hpm
  lChHsnuLstau12 = - lChHsnuLstau12;

  /// LCT: Rearrange masses to account for change in mixing convention
  DoubleVector Hpm(2);
  Hpm(1) = higgsc(2);
  Hpm(2) = higgsc(1); 

  /// LCT: Charged Higgs contribution
  double meL = forLoops.me(1, family);
  DoubleVector msel(2); 
  msel(1) = forLoops.me(1, family);
  msel(2) = forLoops.me(2, family);
  for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      double b0p = b0(p, Hpm(i), msel(j), q); 
      higgs += sqr(lChHsnuLstau12(i, j)) * b0p;
    }

 /// LCT: (CP-even) Higgs-snuetrino-sneutrino couplings in (L R) basis.  
  /// Order (s1 s2 s3)
  DoubleVector lsSnuLSnuLR(3), lsSnuLSnu12(3);

  for (int i=1; i<=3; i++) {
    lsSnuLSnuLR(i) = - 0.25 * gsq * (1.0 + tanthDrbar2) 
      * (v1 * S(i, 1) - v2 * S(i, 2));
  }
  
  /// LCT: CP-even contribution
  for (int i=1; i<=3; i++) {
      higgs += sqr(lsSnuLSnuLR(i)) * b0(p, higgsm(i), msnu(family), q);
  }

  return higgs;
}

void NmssmSoftsusy::addSnuCorrection(double & mass, int family) {

  /// No point adding radiative corrections to tachyonic particles
  if (mass < 0.0) { 
    if (family == 1) flagTachyon(snue);
    else if (family == 2) flagTachyon(snumu);
    return;
  }

  double p = sqrt(mass);

  /// one-loop correction matrix
  double piSq; /// Self-energy 
	
  double higgs = 0., electroweak = 0., chargino = 0., neutralino = 0.;

  /// LCT: Higg contribution
  addSnuHiggs(p, family, higgs);
  /// LCT: Electroweak contributions
  addSnuEweak(p, family, electroweak);
  /// LCT: Gaugino contributions
  addSnuGaugino(p, family, chargino, neutralino);
 
  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (higgs + electroweak + chargino + neutralino);

  mass = mass - piSq;	  
}

DoubleMatrix NmssmSoftsusy::addStauHiggs(double p, double mtau, DoubleMatrix & higgs) {
  const drBarPars& forLoops = displayDrBarPars();
  double    mz      = displayMzRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    alpha   = forLoops.thetaH;
  double    g	    = displayGaugeCoupling(2), gsq = sqr(g);
  double    beta    = atan(displayTanb());

  double    thetatau= forLoops.thetatau;
  double    ctau    = cos(thetatau);
  double    stau    = sin(thetatau);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1)           = forLoops.msnu(1);
  msnu(2)           = forLoops.msnu(2);
  msnu(3)           = forLoops.msnu(3);

  double    smu     = -displaySusyMu();
  double q = displayMu(), 
    htau = forLoops.htau, 
    sinb = sin(beta), cosb = cos(beta), 
    htausq = sqr(htau);
  double tanthDrbar2 = sqr(tan(asin(sinthDrbar)));
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6), gpsq = sqr(gp);
  double mw = displayMwRun();
  double svev = displaySvev();
  double lam = displayLambda();
  double mueff = -(displaySusyMu() + lam * svev / root2); ///<< LCT: Note sign!
  double v2 = displayHvev() * sin(beta);

  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb;
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Quadrilinear contributions
  ///
  /// LCT: stau-stau-Hp-Hm couplings in (L R) basis.  Order (G+ H+).
  DoubleVector lChHChHstauLstauL(2), lChHChHstauRstauR(2);
  for (int i=1; i<=2; i++) {
    lChHChHstauLstauL(i) = - 0.25 * gsq * (1.0 + tanthDrbar2) 
      * (sqr(C(i, 1)) - sqr(C(i, 2)));
    lChHChHstauRstauR(i) = 0.5 * gpsq * (sqr(C(i, 1)) - sqr(C(i, 2))) 
      - htausq * sqr(C(i, 1));
  }
  /// LCT: Charged contribution
  for (int i=1; i<=2; i++) { 
    higgs(1, 1) += - lChHChHstauLstauL(i) * a0(higgsc(i), q);
    higgs(2, 2) += - lChHChHstauRstauR(i) * a0(higgsc(i), q);
  }

  /// LCT: stau-stau-A0-A0 and stau-stau-h0-h0 couplings in (L R) basis.  
  /// Order (G0 A1 A2) and (h1 h2 h3) respectively
  DoubleVector lAAstauLstauL(3), lAAstauRstauR(3), lAAstauLstauR(3), 
    lHHstauLstauL(3), lHHstauRstauR(3), lHHstauLstauR(3);
  for (int i=1; i<=3; i++) {
    lAAstauLstauL(i) = 0.25 * gsq * (1.0 - tanthDrbar2) 
      * (sqr(P(i, 1)) - sqr(P(i, 2))) - htausq * sqr(P(i, 1));
    lAAstauRstauR(i) = 0.5 * gsq * tanthDrbar2 * (sqr(P(i, 1)) - sqr(P(i, 2))) 
      - htausq * sqr(P(i, 1));
    lAAstauLstauR(i) = - lam * htau * P(i, 2) * P(i, 3);

    lHHstauLstauL(i) = 0.25 * gsq * (1.0 - tanthDrbar2) 
      * (sqr(S(i, 1)) - sqr(S(i, 2))) - htausq * sqr(S(i, 1));
    lHHstauRstauR(i) = 0.5 * gsq * tanthDrbar2 * (sqr(S(i, 1)) - sqr(S(i, 2))) 
      - htausq * sqr(S(i, 1));
    lHHstauLstauR(i) = lam * htau * S(i, 2) * S(i, 3);
  }
  /// LCT: CP-odd and CP-even contribution
 for (int i=1; i<=3; i++) { 
   higgs(1, 1) += - 0.5 * (lAAstauLstauL(i) * a0(higgsa(i), q) 
			   + lHHstauLstauL(i) * a0(higgsm(i), q));
   higgs(2, 2) += - 0.5 * (lAAstauRstauR(i) * a0(higgsa(i), q) 
			   + lHHstauRstauR(i) * a0(higgsm(i), q));
   higgs(1, 2) += - 0.5 * (lAAstauLstauR(i) * a0(higgsa(i), q) 
			   + lHHstauLstauR(i) * a0(higgsm(i), q));
  }

 /// LCT: Trilinear contributions
 /// Charged Higgs
 DoubleMatrix lChHstauLsnu12(2, 2); /// (H+ G+, L R) basis
  lChHstauLsnu12(2, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2 
    + htau * mtau * cosb;
  lChHstauLsnu12(1, 1) = g * displayMwRun() * sin(2.0 * beta) / root2 
    - htau * mtau * sinb;

  DoubleMatrix lChHstauRsnu12(2, 2); /// (H+ G+, L R) basis
  lChHstauRsnu12(2, 1) = htau * mueff * sinb + forLoops.utau * cosb;
  lChHstauRsnu12(1, 1) = htau * mueff * cosb - forLoops.utau * sinb;

  /// Flip sign to match Slavich conventions on H+/G+ mixing matrix C
  lChHstauLsnu12 = - lChHstauLsnu12;
  lChHstauRsnu12 = - lChHstauRsnu12;

  /// LCT: Rearrange masses to account for change in mixing convention
  DoubleVector Hpm(2);
  Hpm(1) = higgsc(2);
  Hpm(2) = higgsc(1); 

  /// LCT: Charged Higgs contribution
  for(int i=1; i<=2; i++) {
      double b0p = b0(p, Hpm(i), msnu(3), q); 
      higgs(1, 1) += sqr(lChHstauLsnu12(i, 1)) * b0p;
      higgs(1, 2) += 
  	lChHstauLsnu12(i, 1) * lChHstauRsnu12(i, 1) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lChHstauRsnu12(i, 1)) * b0p;
  }

 /// LCT: A0-stau-stau couplings in (L R) basis
  DoubleVector lAstauLstauR(3), lAstauRstauL(3);
  for (int i=1; i<=3; i++) {
    lAstauLstauR(i) = forLoops.utau / root2 * P(i, 1) 
      + 0.5 * lam * htau * (v2 * P(i, 3) + svev * P(i, 2));
    lAstauRstauL(i) = - lAstauLstauR(i);
  }
  /// LCT: Rotate to mixed basis
  DoubleMatrix lAstauLstau12(3, 2), lAstauRstau12(3, 2);
  for (int i=1; i<=3; i++) {
    lAstauLstau12(i, 1) = lAstauLstauR(i) * stau;
    lAstauLstau12(i, 2) = lAstauLstauR(i) * ctau;
    lAstauRstau12(i, 1) = lAstauRstauL(i) * ctau;
    lAstauRstau12(i, 2) = - lAstauRstauL(i) * stau;
    }
  for (int i=1; i<=3; i++) {
    for (int j=1; j<=2; j++) {
      higgs(1, 1) += sqr(lAstauLstau12(i, j)) * b0(p, higgsa(i), mstau(j), q);
      higgs(1, 2) += lAstauLstau12(i, j) * lAstauRstau12(i, j) 
	* b0(p, higgsa(i), mstau(j), q);
      higgs(2, 2) += sqr(lAstauRstau12(i, j)) * b0(p, higgsa(i), mstau(j), q);
    }
  }

  /// LCT: (CP-even) Higgs-stau-stau couplings in (L R) basis.  Order (h1 h2 h3)
  //
  /// Taken from Franke & Fraas, arXiv:hep-ph/9512366
  DoubleMatrix lsStauLStauLR(3, 2), lsStauLStau12(3, 2);
  DoubleMatrix lsStauRStauLR(3, 2), lsStauRStau12(3, 2);
  double utau = forLoops.utau;  
  
  lsStauLStauLR(1, 1) = g * mz / costhDrbar * geL 
    * (S(1, 2) * sinb - S(1, 1) * cosb) - g * sqr(mtau) / (mw * cosb) * S(1, 1);
  lsStauLStauLR(1, 2) = - utau / root2 * S(1, 1) 
    + 0.5 * g * mtau / (root2 * mw * cosb) * lam * (v2 * S(1, 3) + svev * S(1, 2));
  lsStauLStauLR(2, 1) = g * mz / costhDrbar * geL 
    * (S(2, 2) * sinb - S(2, 1) * cosb) - g * sqr(mtau) / (mw * cosb) * S(2, 1);
  lsStauLStauLR(2, 2) = - utau / root2 * S(2, 1) 
    + 0.5 * g * mtau / (root2 * mw * cosb) * lam * (v2 * S(2, 3) + svev * S(2, 2));
  lsStauLStauLR(3, 1) = g * mz / costhDrbar * geL 
    * (S(3, 2) * sinb - S(3, 1) * cosb) - g * sqr(mtau) / (mw * cosb) * S(3, 1);
  lsStauLStauLR(3, 2) = - utau / root2 * S(3, 1) 
    + 0.5 * g * mtau / (root2 * mw * cosb) * lam * (v2 * S(3, 3) + svev * S(3, 2));
  
  lsStauRStauLR(1, 1) = lsStauLStauLR(1, 2);
  lsStauRStauLR(1, 2) = - g * sqr(mtau) / (mw * cosb) * S(1, 1) 
    - 1.0 * g * mw * tanthDrbar2 * (S(1, 2) * sinb - S(1, 1) * cosb);
  lsStauRStauLR(2, 1) = lsStauLStauLR(2, 2);
  lsStauRStauLR(2, 2) = - g * sqr(mtau) / (mw * cosb) * S(2, 1) 
    - 1.0 * g * mw * tanthDrbar2 * (S(2, 2) * sinb - S(2, 1) * cosb);
  lsStauRStauLR(3, 1) = lsStauLStauLR(3, 2);
  lsStauRStauLR(3, 2) = - g * sqr(mtau) / (mw * cosb) * S(3, 1) 
    - 1.0 * g * mw * tanthDrbar2 * (S(3, 2) * sinb - S(3, 1) * cosb);
  
  /// Mix staus up
  DoubleVector temp(2), temp2(2);
  for (int i=1; i<=3; i++) {
    temp(1) = lsStauLStauLR(i, 1);
    temp(2) = lsStauLStauLR(i, 2);
    temp2 = rot2d(thetatau) * temp;
    lsStauLStau12(i, 1) = temp2(1);
    lsStauLStau12(i, 2) = temp2(2);
    temp(1) = lsStauRStauLR(i, 1);
    temp(2) = lsStauRStauLR(i, 2);
    temp2 = rot2d(thetatau) * temp;
    lsStauRStau12(i, 1) = temp2(1);
    lsStauRStau12(i, 2) = temp2(2);
  }  

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=2; j++) {
      higgs(1, 1) += sqr(lsStauLStau12(i, j)) * b0(p, higgsm(i), mstau(j), q);
      higgs(2, 2) += sqr(lsStauRStau12(i, j)) * b0(p, higgsm(i), mstau(j), q);
      higgs(1, 2) += lsStauLStau12(i, j) * lsStauRStau12(i, j) 
       * b0(p, higgsm(i), mstau(j), q);
    }
  }

 higgs(2, 1) = higgs(1, 2);

  return higgs;
}

void NmssmSoftsusy::addStauCorrection(double p, DoubleMatrix & mass, 
				     double mtau) {

/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    flagTachyon(stau);
    if (mass(1, 1) < 0.) mass(1, 1) = EPSTOL;
    else mass(2, 2) = EPSTOL;
    return;
  }

  /// one-loop correction matrix
  DoubleMatrix piSq(2, 2); /// Self-energy matrix
	
  DoubleMatrix stop(2, 2), sbottom(2, 2), 
    higgs(2, 2), electroweak(2, 2), chargino(2, 2), neutralino(2, 2);

  /// LCT: Sfermion contribution
  addStauSfermion(p, mtau, stop, sbottom);
  /// LCT: Higgs contributions
  addStauHiggs(p, mtau, higgs);
  /// LCT: Electroweak contribution
  addStauEweak(p, mtau, electroweak);
  /// LCT: Chargino and neutralino contributions
  addStauGaugino(p, mtau, chargino, neutralino);

  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (stop + sbottom + higgs + chargino + neutralino + electroweak);

  piSq(2, 1) = piSq(1, 2);

  
  mass = mass - piSq;	
}

DoubleMatrix NmssmSoftsusy::addSlepHiggs(double p1, double p2, int family, DoubleMatrix & higgs) {
  const drBarPars& forLoops = displayDrBarPars();
  double    sinthDrbar  = calcSinthdrbar();
  double tanthDrbar2 = sqr(tan(asin(sinthDrbar)));
  double    mw      = displayMwRun();
  double    mz      = displayMzRun();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    alpha   = forLoops.thetaH;
  double    g	    = displayGaugeCoupling(2), gsq = sqr(g);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6), gpsq = sqr(gp);
  double    e       = g * sinthDrbar;
  double    beta    = atan(displayTanb());

  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double    thetatau= forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double    ctau    = cos(thetatau);
  double    stau    = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msel(2);
  msel(1)           = forLoops.me(1, family);
  msel(2)           = forLoops.me(2, family);
  DoubleVector msnu(3);
  msnu(1)           = forLoops.msnu(1);
  msnu(2)           = forLoops.msnu(2);
  msnu(3)           = forLoops.msnu(3);
 double q = displayMu(), sinb = sin(beta), cosb = cos(beta);

  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb;
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Quadrilinear contributions
  ///
  /// LCT: stau-stau-Hp-Hm couplings in (L R) basis.  Order (G+ H+).
  DoubleVector lChHChHstauLstauL(2), lChHChHstauRstauR(2);
  for (int i=1; i<=2; i++) {
    lChHChHstauLstauL(i) = - 0.25 * gsq * (1.0 + tanthDrbar2) 
      * (sqr(C(i, 1)) - sqr(C(i, 2)));
    lChHChHstauRstauR(i) = 0.5 * gpsq * (sqr(C(i, 1)) - sqr(C(i, 2)));
  }
  /// LCT: Charged contribution
  for (int i=1; i<=2; i++) { 
    higgs(1, 1) += - lChHChHstauLstauL(i) * a0(higgsc(i), q);
    higgs(2, 2) += - lChHChHstauRstauR(i) * a0(higgsc(i), q);
  }

  /// LCT: stau-stau-A0-A0 and stau-stau-h0-h0 couplings in (L R) basis.  
  /// Order (G0 A1 A2) and (h1 h2 h3) respectively
  DoubleVector lAAstauLstauL(3), lAAstauRstauR(3), lAAstauLstauR(3), 
    lHHstauLstauL(3), lHHstauRstauR(3), lHHstauLstauR(3);
  for (int i=1; i<=3; i++) {
    lAAstauLstauL(i) = 0.25 * gsq * (1.0 - tanthDrbar2) 
      * (sqr(P(i, 1)) - sqr(P(i, 2)));
    lAAstauRstauR(i) = 0.5 * gsq * tanthDrbar2 * (sqr(P(i, 1)) - sqr(P(i, 2)));

    lHHstauLstauL(i) = 0.25 * gsq * (1.0 - tanthDrbar2) 
      * (sqr(S(i, 1)) - sqr(S(i, 2)));
    lHHstauRstauR(i) = 0.5 * gsq * tanthDrbar2 * (sqr(S(i, 1)) - sqr(S(i, 2)));
  }
  /// LCT: CP-odd and CP-even contribution
 for (int i=1; i<=3; i++) { 
   higgs(1, 1) += - 0.5 * (lAAstauLstauL(i) * a0(higgsa(i), q) 
			   + lHHstauLstauL(i) * a0(higgsm(i), q));
   higgs(2, 2) += - 0.5 * (lAAstauRstauR(i) * a0(higgsa(i), q) 
			   + lHHstauRstauR(i) * a0(higgsm(i), q));
  }

/// LCT: Trilinear contributions
 /// Charged Higgs
  /// Charged Higgs Feynman rules
  DoubleMatrix lChHstauLsnu12(2, 2); /// (H+ G+, L R) basis
  lChHstauLsnu12(2, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2;
  lChHstauLsnu12(1, 1) = g * displayMwRun() * sin(2.0 * beta) / root2;

  /// Flip sign to match Slavich conventions on H+/G+ mixing matrix C
  lChHstauLsnu12 = - lChHstauLsnu12;

  /// LCT: Rearrange masses to account for change in mixing convention
  DoubleVector Hpm(2);
  Hpm(1) = higgsc(2);
  Hpm(2) = higgsc(1); 

  /// LCT: Charged Higgs contribution
  for(int i=1; i<=2; i++) {
      double b0p = b0(p1, Hpm(i), msnu(family), q); 
      higgs(1, 1) += sqr(lChHstauLsnu12(i, 1)) * b0p;
  }

  /// LCT: Trilinear contributions
  /// No CP-odd higgs terms

 /// LCT: (CP-even) Higgs-stau-stau couplings in (L R) basis.  Order (h1 h2 h3)
  //
  /// Taken from Franke & Fraas, arXiv:hep-ph/9512366
  DoubleMatrix lsStauLStauLR(3, 2), lsStauLStau12(3, 2);
  DoubleMatrix lsStauRStauLR(3, 2), lsStauRStau12(3, 2);
  double utau = forLoops.utau;  
  
  lsStauLStauLR(1, 1) = g * mz / costhDrbar * geL 
    * (S(1, 2) * sinb - S(1, 1) * cosb);
  lsStauLStauLR(2, 1) = g * mz / costhDrbar * geL 
    * (S(2, 2) * sinb - S(2, 1) * cosb);
  lsStauLStauLR(3, 1) = g * mz / costhDrbar * geL 
    * (S(3, 2) * sinb - S(3, 1) * cosb);
  
  lsStauRStauLR(1, 1) = lsStauLStauLR(1, 2);
  lsStauRStauLR(1, 2) = - 1.0 * g * mw * tanthDrbar2 * (S(1, 2) * sinb - S(1, 1) * cosb);
  lsStauRStauLR(2, 1) = lsStauLStauLR(2, 2);
  lsStauRStauLR(2, 2) = - 1.0 * g * mw * tanthDrbar2 * (S(2, 2) * sinb - S(2, 1) * cosb);
  lsStauRStauLR(3, 1) = lsStauLStauLR(3, 2);
  lsStauRStauLR(3, 2) = - 1.0 * g * mw * tanthDrbar2 * (S(3, 2) * sinb - S(3, 1) * cosb);
  
  /// Mix staus up
  DoubleVector temp(2), temp2(2);
  for (int i=1; i<=3; i++) {
    temp(1) = lsStauLStauLR(i, 1);
    temp(2) = lsStauLStauLR(i, 2);
    temp2 = temp;
    lsStauLStau12(i, 1) = temp2(1);
    lsStauLStau12(i, 2) = temp2(2);
    temp(1) = lsStauRStauLR(i, 1);
    temp(2) = lsStauRStauLR(i, 2);
    temp2 = temp;
    lsStauRStau12(i, 1) = temp2(1);
    lsStauRStau12(i, 2) = temp2(2);
  }  

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=2; j++) {
      double b0p1 = b0(p1, higgsm(i), msel(j), q);
      double b0p2 = b0(p2, higgsm(i), msel(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lsStauLStau12(i, j)) * b0p1;
      higgs(2, 2) = higgs(2, 2) + sqr(lsStauRStau12(i, j)) * b0p2;
    }
  }

 return higgs;
}

void NmssmSoftsusy::addSlepCorrection(DoubleMatrix & mass, int family) {

/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    if (mass(1, 1) < 0.) mass(1, 1) = EPSTOL;
    else mass(2, 2) = EPSTOL;
    if (family == 1) flagTachyon(selectron);
    if (family == 2) flagTachyon(smuon);
    return;
  }

  /// one-loop correction matrix
  DoubleMatrix piSq(2, 2); /// Self-energy matrix

  const drBarPars& forLoops = displayDrBarPars();
  DoubleVector msel(2);
  msel(1)           = forLoops.me(1, family);
  msel(2)           = forLoops.me(2, family);	
  double p1 = msel(1), p2 = msel(2);
  ///  double    p       = sqrt(msel(1) * msel(2));

  DoubleMatrix higgs(2, 2), chargino(2, 2), neutralino(2, 2), 
    electroweak(2, 2); 

  /// LCT: Higgs contribution
  addSlepHiggs(p1, p2, family, higgs);
  /// LCT: Electroweak contribution
  addSlepEweak(p1, p2, family, electroweak);
  /// LCT: Chargino and neutralino contribution
  addSlepGaugino(p1, p2, family, chargino, neutralino);

  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (higgs + chargino + neutralino + electroweak);

  mass = mass - piSq;	
}

void NmssmSoftsusy::getNeutPassarinoVeltman(double p, double q, DoubleMatrix & b0fn, DoubleMatrix & b1fn) {
  const drBarPars& forLoops = displayDrBarPars();
  double    mz = displayMzRun(), mw = displayMwRun();
  DoubleVector mneut(forLoops.mnBpmz);
  DoubleVector mch(forLoops.mchBpmz); 

  const int rank = mneut.displayEnd();
  int k;
     for (k=1; k<=rank; k++) {
       if (k<=2){
	 /// LCT: W+
	 b0fn(k, 1) = b0(p, mch(k), mw, q);
	 b1fn(k, 1) = b1(p, mch(k), mw, q);
	 /// LCT: H+
	 b0fn(k, 2) = b0(p, mch(k), forLoops.mHpm, q); 
	 b1fn(k, 2) = b1(p, mch(k), forLoops.mHpm, q);
       }
 	
       /// LCT: Z0
       b0fn(k, 3) = b0(p, mneut(k), mz, q);
       b1fn(k, 3) = b1(p, mneut(k), mz, q);

       /// LCT: A1 A2
       b0fn(k, 4) = b0(p, mneut(k), forLoops.mA0(1), q);
       b1fn(k, 4) = b1(p, mneut(k), forLoops.mA0(1), q);
       b0fn(k, 5) = b0(p, mneut(k), forLoops.mA0(2), q);
       b1fn(k, 5) = b1(p, mneut(k), forLoops.mA0(2), q);

       /// LCT: h1 h2 h3
       b0fn(k, 6) = b0(p, mneut(k), forLoops.mh0(1), q);
       b1fn(k, 6) = b1(p, mneut(k), forLoops.mh0(1), q);
       b0fn(k, 7) = b0(p, mneut(k), forLoops.mh0(2), q);
       b1fn(k, 7) = b1(p, mneut(k), forLoops.mh0(2), q);
       b0fn(k, 8) = b0(p, mneut(k), forLoops.mh0(3), q);
       b1fn(k, 8) = b1(p, mneut(k), forLoops.mh0(3), q);
     }
}

void NmssmSoftsusy::addNeutLoopHiggs(double p, DoubleMatrix & sigmaL, 
DoubleMatrix & sigmaR, DoubleMatrix & sigmaS) {
  const drBarPars& forLoops = displayDrBarPars();
  double g = displayGaugeCoupling(2), 
    gp = displayGaugeCoupling(1) * sqrt(0.6), 
    ht = displayDrBarPars().ht,
    htau = displayDrBarPars().htau, 
    hb = displayDrBarPars().hb,
    q = displayMu(), 
    tanb = displayTanb(),
    mz = displayMzRun(), mw = displayMwRun();
  double beta = atan(tanb);
  double sinb = sin(beta), cosb = cos(beta);
  double mt = ht * displayHvev() * sin(beta) / root2,
    mb = hb * displayHvev() * cos(beta) / root2,
    mtau = htau * displayHvev() * cos(beta) / root2;
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double    thetatau= forLoops.thetatau;
  double lam = displayLambda();
  double kap = displayKappa();

  DoubleVector msup(2), msdown(2), msel(2), mscharm(2), msstrange(2), 
    msmuon(2), msnumu(2), mstop(2), msbot(2), mstau(2), msnutau(2), msnue(2);
  msup(1) = forLoops.mu(1, 1);      msup(2) = forLoops.mu(2, 1); 
  mscharm(1) = forLoops.mu(1, 2);   mscharm(2) = forLoops.mu(2, 2); 
  mstop(1) = forLoops.mu(1, 3);     mstop(2) = forLoops.mu(2, 3); 
  msdown(1) = forLoops.md(1, 1);    msdown(2) = forLoops.md(2, 1); 
  msstrange(1) = forLoops.md(1, 2); msstrange(2) = forLoops.md(2, 2); 
  msbot(1) = forLoops.md(1, 3);     msbot(2) = forLoops.md(2, 3); 
  msel(1) = forLoops.me(1, 1);      msel(2) = forLoops.me(2, 1); 
  msmuon(1) = forLoops.me(1, 2);    msmuon(2) = forLoops.me(2, 2); 
  mstau(1) = forLoops.me(1, 3);     mstau(2) = forLoops.me(2, 3); 
  msnue(1) = forLoops.msnu(1); 
  msnumu(1) = forLoops.msnu(2);
  msnutau(1) = forLoops.msnu(3);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz);
  const int rank = mneut.displayEnd();

  /// LCT: Neutralino-chargino-Hpm trilinear coupings
  DoubleMatrix aPsiPsiHc1(rank, 2), bPsiPsiHc1(rank, 2);
  DoubleMatrix aPsiPsiHc2(rank, 2), bPsiPsiHc2(rank, 2);
  aPsiPsiHc1(1, 2) = gp / root2;
  bPsiPsiHc2(1, 2) = aPsiPsiHc1(1, 2);
  aPsiPsiHc1(2, 2) = g / root2;
  bPsiPsiHc2(2, 2) = g / root2;
  aPsiPsiHc1(3, 1) = -g;
  bPsiPsiHc2(4, 1) = g;
  aPsiPsiHc2(5, 2) = lam;
  bPsiPsiHc1(5, 2) = -lam;
  ComplexMatrix aPsiChiHc1(rank, 2), bPsiChiHc1(rank, 2);
  ComplexMatrix aPsiChiHc2(rank, 2), bPsiChiHc2(rank, 2);
  aPsiChiHc1 = aPsiPsiHc1 * u.hermitianConjugate();
  aPsiChiHc2 = aPsiPsiHc2 * u.hermitianConjugate();
  bPsiChiHc1 = bPsiPsiHc1 * v.transpose();
  bPsiChiHc2 = bPsiPsiHc2 * v.transpose();
  ComplexMatrix aPsiChiHHp(rank, 2), bPsiChiHHp(rank, 2);
  ComplexMatrix aPsiChiHGp(rank, 2), bPsiChiHGp(rank, 2);
  int i,j,k; for (i=1; i<=rank; i++)
    for (j=1; j<=2; j++) {
      aPsiChiHGp(i, j) = cosb * aPsiChiHc1(i, j) + sinb * aPsiChiHc2(i, j);
      bPsiChiHGp(i, j) = cosb * bPsiChiHc1(i, j) + sinb * bPsiChiHc2(i, j);
      aPsiChiHHp(i, j) =-sinb * aPsiChiHc1(i, j) + cosb * aPsiChiHc2(i, j);
      bPsiChiHHp(i, j) =-sinb * bPsiChiHc1(i, j) + cosb * bPsiChiHc2(i, j);
    }
  
  /// checked this block
  ComplexMatrix aPsiPsis1(rank, rank), aPsiPsis2(rank, rank), 
    aPsiPsis3(rank, rank), aPsiPsip1(rank, rank), aPsiPsip2(rank, rank), 
    aPsiPsip3(rank, rank);
  ComplexMatrix bPsiPsis1(rank, rank), bPsiPsis2(rank, rank), 
    bPsiPsis3(rank, rank), bPsiPsip1(rank, rank), bPsiPsip2(rank, rank), 
    bPsiPsip3(rank, rank);
  ComplexMatrix aPsiChis1(rank, rank), aPsiChis2(rank, rank), 
    aPsiChis3(rank, rank), aPsiChip1(rank, rank), aPsiChip2(rank, rank), 
    aPsiChip3(rank, rank);
  ComplexMatrix bPsiChis1(rank, rank), bPsiChis2(rank, rank), 
    bPsiChis3(rank, rank), bPsiChip1(rank, rank), bPsiChip2(rank, rank), 
    bPsiChip3(rank, rank);
  aPsiPsis1(1, 3) = - gp * 0.5;
  aPsiPsis1(2, 3) = g * 0.5;
  aPsiPsis1(4, 5) = - lam / root2;
  aPsiPsis2(2, 4) = - g * 0.5;
  aPsiPsis2(1, 4) = gp * 0.5;
  aPsiPsis2(3, 5) = - lam / root2;
  aPsiPsis3(3, 4) = - lam / root2;
  aPsiPsis3(5, 5) = root2 * kap;

  aPsiPsip1(1, 3) = - gp * 0.5;
  aPsiPsip1(2, 3) = g * 0.5;
  aPsiPsip1(4, 5) = lam / root2;
  aPsiPsip2(2, 4) = g * 0.5;
  aPsiPsip2(1, 4) = - gp * 0.5;
  aPsiPsip2(3, 5) = - lam / root2; /// LCT: Flipped Staub sign to match BPMZ
  aPsiPsip3(3, 4) = lam / root2;
  aPsiPsip3(5, 5) = - root2 * kap;
 
  aPsiPsis1.symmetrise();
  aPsiPsis2.symmetrise();
  aPsiPsis3.symmetrise();
  aPsiPsip1.symmetrise();
  aPsiPsip2.symmetrise();
  aPsiPsip3.symmetrise();

  bPsiPsis1 = aPsiPsis1;
  bPsiPsis2 = aPsiPsis2;
  bPsiPsis3 = aPsiPsis3;
  bPsiPsip1 = - 1.0 * aPsiPsip1;
  bPsiPsip2 = - 1.0 * aPsiPsip2;
  bPsiPsip3 = - 1.0 * aPsiPsip3;

  aPsiChis1 = aPsiPsis1 * n.hermitianConjugate();
  aPsiChis2 = aPsiPsis2 * n.hermitianConjugate();
  aPsiChis3 = aPsiPsis3 * n.hermitianConjugate();
  aPsiChip1 = aPsiPsip1 * n.hermitianConjugate();
  aPsiChip2 = aPsiPsip2 * n.hermitianConjugate();
  aPsiChip3 = aPsiPsip3 * n.hermitianConjugate();
  bPsiChis1 = bPsiPsis1 * n.transpose();
  bPsiChis2 = bPsiPsis2 * n.transpose();
  bPsiChis3 = bPsiPsis3 * n.transpose();
  bPsiChip1 = bPsiPsip1 * n.transpose();
  bPsiChip2 = bPsiPsip2 * n.transpose();
  bPsiChip3 = bPsiPsip3 * n.transpose();

  /// LCT: Rotate to Higgs mass basis
  ComplexMatrix aPsiChih1(rank, rank), aPsiChih2(rank, rank), 
    aPsiChih3(rank, rank), aPsiChiG(rank, rank), aPsiChiA1(rank, rank), 
    aPsiChiA2(rank, rank);
  ComplexMatrix bPsiChih1(rank, rank), bPsiChih2(rank, rank), 
    bPsiChih3(rank, rank), bPsiChiG(rank, rank), bPsiChiA1(rank, rank), 
    bPsiChiA2(rank, rank);
  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;

  aPsiChih1 = S(1, 1) * aPsiChis1 + S(1, 2) * aPsiChis2 + S(1, 3) * aPsiChis3;
  aPsiChih2 = S(2, 1) * aPsiChis1 + S(2, 2) * aPsiChis2 + S(2, 3) * aPsiChis3;
  aPsiChih3 = S(3, 1) * aPsiChis1 + S(3, 2) * aPsiChis2 + S(3, 3) * aPsiChis3;
  bPsiChih1 = S(1, 1) * bPsiChis1 + S(1, 2) * bPsiChis2 + S(1, 3) * bPsiChis3;
  bPsiChih2 = S(2, 1) * bPsiChis1 + S(2, 2) * bPsiChis2 + S(2, 3) * bPsiChis3;
  bPsiChih3 = S(3, 1) * bPsiChis1 + S(3, 2) * bPsiChis2 + S(3, 3) * bPsiChis3;

  aPsiChiG  = - P(1, 1) * aPsiChip1 + P(1, 2) * aPsiChip2 + P(1, 3) * aPsiChip3;
  aPsiChiA1 =   P(2, 1) * aPsiChip1 - P(2, 2) * aPsiChip2 + P(2, 3) * aPsiChip3;
  aPsiChiA2 =   P(3, 1) * aPsiChip1 - P(3, 2) * aPsiChip2 + P(3, 3) * aPsiChip3;
  bPsiChiG  = - P(1, 1) * bPsiChip1 + P(1, 2) * bPsiChip2 + P(1, 3) * bPsiChip3;
  bPsiChiA1 =   P(2, 1) * bPsiChip1 - P(2, 2) * bPsiChip2 + P(2, 3) * bPsiChip3;
  bPsiChiA2 =   P(3, 1) * bPsiChip1 - P(3, 2) * bPsiChip2 + P(3, 3) * bPsiChip3;

  DoubleMatrix b0fn(rank, 8), b1fn(rank, 8);
  getNeutPassarinoVeltman(p, q, b0fn, b1fn);

  for (i=1; i<=rank; i++) 
    for (j=1; j<=rank; j++) 
      for (k=1; k<=rank; k++) {
	if (k<=2) {
  	/// G+
	sigmaL(i, j) = sigmaL(i, j) + 
	  (aPsiChiHGp(i, k).conj() * aPsiChiHGp(j, k) * b1fn(k, 1)).real();
	sigmaR(i, j) = sigmaR(i, j) + 
	  (bPsiChiHGp(i, k).conj() * bPsiChiHGp(j, k) * b1fn(k, 1)).real();
	sigmaS(i, j) = sigmaS(i, j) + 2.0 * mch(k) * 
	  (bPsiChiHGp(i, k).conj() * aPsiChiHGp(j, k) * b0fn(k, 1)).real();

	/// H+
	sigmaL(i, j) = sigmaL(i, j) + 
	  (aPsiChiHHp(i, k).conj() * aPsiChiHHp(j, k) * b1fn(k, 2)).real();
	sigmaR(i, j) = sigmaR(i, j) + 
	  (bPsiChiHHp(i, k).conj() * bPsiChiHHp(j, k) * b1fn(k, 2)).real();
	sigmaS(i, j) = sigmaS(i, j) + 2.0 * mch(k) * 
	  (bPsiChiHHp(i, k).conj() * aPsiChiHHp(j, k) * b0fn(k, 2)).real();
	}

	  /// h1
  	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
  	  (aPsiChih1(i, k).conj() * aPsiChih1(j, k) * b1fn(k, 6)).real();
  	sigmaR(i, j) = sigmaR(i, j) + 0.5 * 
  	  (bPsiChih1(i, k).conj() * bPsiChih1(j, k) * b1fn(k, 6)).real();
  	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
  	  (bPsiChih1(i, k).conj() * aPsiChih1(j, k) * b0fn(k, 6)).real();

	/// h2
	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
  	  (aPsiChih2(i, k).conj() * aPsiChih2(j, k) * b1fn(k, 7)).real();
  	sigmaR(i, j) = sigmaR(i, j) + 0.5 * 
  	  (bPsiChih2(i, k).conj() * bPsiChih2(j, k) * b1fn(k, 7)).real();
  	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
  	  (bPsiChih2(i, k).conj() * aPsiChih2(j, k) * b0fn(k, 7)).real();

	/// h3
	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
  	  (aPsiChih3(i, k).conj() * aPsiChih3(j, k) * b1fn(k, 8)).real();
  	sigmaR(i, j) = sigmaR(i, j) + 0.5 * 
  	  (bPsiChih3(i, k).conj() * bPsiChih3(j, k) * b1fn(k, 8)).real();
  	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
  	  (bPsiChih3(i, k).conj() * aPsiChih3(j, k) * b0fn(k, 8)).real();

	/// G0
  	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
  	  (aPsiChiG(i, k).conj() * aPsiChiG(j, k) * b1fn(k, 3)).real();
  	sigmaR(i, j) = sigmaR(i, j) + 0.5 *
  	  (bPsiChiG(i, k).conj() * bPsiChiG(j, k) * b1fn(k, 3)).real();
  	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
  	  (bPsiChiG(i, k).conj() * aPsiChiG(j, k) * b0fn(k, 3)).real();

	/// A1
 	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
  	  (aPsiChiA1(i, k).conj() * aPsiChiA1(j, k) * b1fn(k, 4)).real();
  	sigmaR(i, j) = sigmaR(i, j) + 0.5 *
  	  (bPsiChiA1(i, k).conj() * bPsiChiA1(j, k) * b1fn(k, 4)).real();
  	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
  	  (bPsiChiA1(i, k).conj() * aPsiChiA1(j, k) * b0fn(k, 4)).real();

	/// A2
 	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
  	  (aPsiChiA2(i, k).conj() * aPsiChiA2(j, k) * b1fn(k, 5)).real();
  	sigmaR(i, j) = sigmaR(i, j) + 0.5 *
  	  (bPsiChiA2(i, k).conj() * bPsiChiA2(j, k) * b1fn(k, 5)).real();
  	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
  	  (bPsiChiA2(i, k).conj() * aPsiChiA2(j, k) * b0fn(k, 5)).real();
      }
}

void NmssmSoftsusy::addNeutralinoLoop(double p, DoubleMatrix & mass) {
  DoubleMatrix sigmaL(5, 5), sigmaR(5, 5), sigmaS(5, 5);

  /// LCT: Sfermion contribution
  addNeutLoopSfermion(p, sigmaL, sigmaR, sigmaS);
  /// LCT: Gauge contributions
  addNeutLoopGauge(p, sigmaL, sigmaR, sigmaS);
  /// LCT: Higgs contribution
  addNeutLoopHiggs(p, sigmaL, sigmaR, sigmaS);

  DoubleMatrix deltaM(5, 5);
  deltaM = -sigmaR * mass - mass * sigmaL - sigmaS;
  deltaM = (deltaM + deltaM.transpose()) / (32.0 * sqr(PI));

  mass = mass + deltaM;
}

/// Returns mu from rewsb requirement. 
/// returns 1 if there's a problem. Call at MSusy
//PA: To be used in general Z3 violating nmssm 
int NmssmSoftsusy::rewsbMu(int sgnMu, double & mu) const {
  int flag = 0;
   if (abs(sgnMu) != 1) {
    ostringstream ii;     
    ii << "Error: sign mu = " << sgnMu << "\n";
    throw ii.str();
  }
   double mH1sq = displayMh1Squared(), mH2sq = displayMh2Squared(), tanb =
     displayTanb(), tanb2 =  sqr(tanb), lam = displayLambda(), s = displaySvev();
  /// PA: Treelevel relation just like in Mssm,
  //but now \mueff = \mu + lambda s / root2
  double musqeff = (mH1sq - mH2sq * tanb2) 
    / (tanb2 - 1.0) - 0.5 * sqr(displayMz());
  
  double musq = musqeff - 0.5 * sqr(lam) * sqr(s) ;
  //PA real root test
  // bsqM4ac = 2.0 * sqr(lam * s) - 4 * 0.5 * sqr(lam * s) + 4 * musqeff;
  double bsqM4ac =  4.0 * musqeff;
 
  if(bsqM4ac < 0) {
    mu = - lam * s / root2; //PA: take real part
    flag = 1; //PA: mu is not real!
  }
  else{
    mu = - lam * s / root2 + 0.5 * sgnMu * sqrt(bsqM4ac);
  } 
  
  return flag;
  
}

// PA: NMssm rewsb routine which fixes imn much the same way as 
// mu is fixed in the Mssm using mueff = lambda * s / root 
// For use in Z3 constrained version or when other scenarios 
// where mu = 0
int NmssmSoftsusy::rewsbSvev(int sgnMu, double & svev) const {
  int flag = 0;
   if (abs(sgnMu) != 1) {
    ostringstream ii;     
    ii << "Error: sign mu = " << sgnMu << "\n";
    throw ii.str();
  }
   double mH1sq = displayMh1Squared(), mH2sq = displayMh2Squared(), tanb =
     displayTanb(), tanb2 =  sqr(tanb), lam = displayLambda();
  
  /// PA: Tree-level relation just like in Mssm,
  //but now \mu --> \mueff = lambda s / root2
  double musqeff = (mH1sq - mH2sq * tanb2) 
    / (tanb2 - 1.0) - 0.5 * sqr(displayMz());
   if (musqeff < 0.0) flag = 1; /// mu has incorrect sign
     double ssq = 2.0 * musqeff / sqr(lam);
     svev = sgnMu * sqrt(fabs(ssq));

  return flag;
  
}

/// returns 1 if mu < 1.0e-9
//PA:  nmssm version for use in Z3 violating case.  
int NmssmSoftsusy::rewsbM3sq(double mu, double & m3sq) const {
   int flag = 0;

   if (fabs(mu) < 1.0e-9)
     { flag = 1; m3sq = 0.0;} 
   else{
     double lam = displayLambda();
     double s = displaySvev();   
     double vev =  displayHvev();
     double al = displayTrialambda();
     double mupr = displayMupr();
     double kap = displayKappa();
     double mH1sq = displayMh1Squared();
     double mH2sq = displayMh2Squared();
     double t1 = displayTadpole1Ms();
     double t2 = displayTadpole2Ms();
     double s2b = sin(2 * atan(displayTanb()));
     //PA now using:
     // m3sqeff = m3sq  + lam * (mupr * svev / root2 + xiF)
     // + al * svev / root2  +   0.5 * lam * kap * sqr(svev); 
     // and swapping musq on the RHS for musqeff
     double  m3sqeff  =  0.5 * 
        (mH1sq + mH2sq - t2 - t1 + 2.0 * sqr(mu + lam * s / root2)  
         + 0.5 * sqr(lam) * sqr(vev)) * s2b;
     
        m3sq = m3sqeff - 0.5 * (root2 * al * s + lam * kap * s * s)
           - lam * (mupr * s / root2 + displayXiF());         
   }
  
  /// Following means no good rewsb
  if (m3sq < 0.0) flag = 1;
  
  return flag;
}



//PA:: In case of Z3 invariance EWSB outputs kappa instead.
int NmssmSoftsusy::rewsbKap(double & kap) const {
  int flag = 0;
  if(abs(displayLambda()) < 1e-99){
    if(PRINTOUT) cout << "Warning: called with lambd = 0." << endl;
    if(PRINTOUT) cout << "rewsbKap routine rewquires non-zero lambda." << endl;
    flag = 2;
  }
  double lam = displayLambda();
  double s = displaySvev();    
  double vev =  displayHvev();
  double al = displayTrialambda();
  double mupr = displayMupr();
  double xiF = displayXiF();
  double m3sq = displayM3Squared();
  double mu = displaySusyMu();
  double mH1sq = displayMh1Squared();
  double mH2sq = displayMh2Squared();
  double t1 = displayTadpole1Ms();
  double t2 = displayTadpole2Ms();
  double s2b = sin(2 * atan(displayTanb()));
  
    
  double  m3sqeff  =  0.5 * 
    (mH1sq + mH2sq - t2 - t1 + 2.0 * sqr(mu + lam * s / root2)  
     + 0.5 * sqr(lam * vev) ) * s2b;
  
  //PA this doesn't work when lam = 0.  Should add warning.
  kap = (m3sqeff - m3sq - lam * (mupr * s / root2 + xiF)
	 - al * s / root2) * 2.0 / (lam * s * s);
  
  
  return flag;
}

//PA: third EWSB condition (for the singlet Higgs field) 
//new with respect to the MSSM.
int NmssmSoftsusy::rewsbXiS(double & xiS) const {
   double mSsq = displayMsSquared();
   double mSpsq = displayMspSquared();
   double mupr = displayMupr();
   double kap = displayKappa(); 
   double lam = displayLambda();
   double al = displayTrialambda();
   double ak = displayTriakappa();
   double xiF = displayXiF();
   double s = displaySvev();
   double sin2b = sin(2.0 * atan(displayTanb()));
   double vev = displayHvev();
   double smu = displaySusyMu();
   xiS = - s / root2 * (mSsq - displayTadpoleSMs() + mSpsq + sqr(mupr) 
                + 2.0 * kap * xiF + ak * s / root2 + sqr(kap) * sqr(s)
                + 3.0 * kap * s * mupr / root2 + 0.5 * sqr(lam) * sqr(vev)
                          - 0.5 * lam * kap * sqr(vev) * sin2b ) 
      - 0.5 * smu * lam * sqr(vev) - xiF * mupr
      + 0.25 * sqr(vev) * sin2b * (al + lam * mupr);

 return 0;
}

// PA: for Z3 invariant NMSSM where we solve for s, kappa and mS
// Or low energy non-universal Higgs versions
int NmssmSoftsusy::rewsbmSsq(double & mSsq) const {
  double kap   = displayKappa(); 
  double lam   = displayLambda();
  double alam  = displayTrialambda();
  double akap  = displayTriakappa();
  double s     = displaySvev();
  double sin2b = sin(2.0 * atan(displayTanb()));
  double vev   = displayHvev();
  double smu   = displaySusyMu();   
  double mupr  = displayMupr();
  double xiF   = displayXiF();
  double xiS   = displayXiS();
  double mSprsq = displayMspSquared();
  
  //PA: Z3 part first
  mSsq = - akap * s / root2 - sqr(kap * s) - 0.5 * sqr(lam * vev) 
    + 0.5 * lam * kap * sqr(vev) * sin2b 
    + 0.25 * alam * sqr(vev) * sin2b * root2 / s ;
  //PA: now Z3 violating terms added on part added on.
  mSsq = mSsq - mSprsq - sqr(mupr) - 2.0 * kap * xiF  
    - 3.0 * kap * s * mupr / root2 - (xiS + xiF * mupr) * root2 / s 
    + 0.25 * lam * mupr * sqr(vev) * sin2b * root2 / s 
    -  smu * lam * sqr(vev) / (root2 * s) + displayTadpoleSMs();
  
  return 0;
}

//PA: for low energy non-universal Higgs versions 
int NmssmSoftsusy::rewsbmH1sq(double & mH1sq) const {
  double kap = displayKappa(); 
  double lam = displayLambda();
  double alam = displayTrialambda();
  double svev = displaySvev();
  double tb = displayTanb();
  double c2b = cos(2.0 * atan(displayTanb()));
  double vev = displayHvev();
  double vu = vev * sin(atan(displayTanb()));
  double vd = vev * cos(atan(displayTanb()));
  double mupr = displayMupr();
  double mueff = displaySusyMu() + lam * svev / root2;
  double mueff2 = sqr(mueff);
  double mz2 = sqr(displayMzRun());
  double xiF = displayXiF();
  double m3sq = displayM3Squared();
  double m3sqeff = m3sq  + lam * (mupr * svev / root2 + xiF)
     + alam * svev / root2  +   0.5 * lam * kap * sqr(svev); 
  double t1 =  displayTadpole1Ms();

  mH1sq = tb * m3sqeff - 0.5 * mz2 * c2b - 0.5 * sqr(lam * vu) - mueff2 + t1;
 
   return 0;
   

}

//PA: for low energy non-universal Higgs versions 
int NmssmSoftsusy::rewsbmH2sq(double & mH2sq) const {
  double kap = displayKappa(); 
  double lam = displayLambda();
  double alam = displayTrialambda();
  double svev = displaySvev();
  double tb = displayTanb();
  double c2b = cos(2.0 * atan(displayTanb()));
  double vev = displayHvev();
  double vu = vev * sin(atan(displayTanb()));
  double vd = vev * cos(atan(displayTanb()));
  double mupr = displayMupr();
  double mueff = displaySusyMu() + lam * svev / root2;
  double mueff2 = sqr(mueff);
  double mz2 = sqr(displayMzRun());
  double xiF = displayXiF();
  double m3sq = displayM3Squared();
  double m3sqeff = m3sq  + lam * (mupr * svev / root2 + xiF)
    + alam * svev / root2  +   0.5 * lam * kap * sqr(svev);
  double t2 =  displayTadpole2Ms();
  
  mH2sq = m3sqeff / tb + 0.5 * mz2 * c2b - 0.5 * sqr(lam * vd) - mueff2 + t2;
  
  return 0;
}


//PA: Imposes EWSB at the tree level. 
// Curently works for general nmssm mapping
// mu --> mZ, m3sq --> tan beta, s --> XiS  (Z3 = false) 
// ie (mu, m3sq, XiS) --> (mZ, tb, s) 
//and s --> mZ, kappa --> tan beta, mS --> s  (Z3 = true)
//ie (kappa, mS) --> (mZ, tb)   
void NmssmSoftsusy::rewsbTreeLevel(int sgnMu) {
  double mu, m3sq, s, kap;
  double xiS, mSsq;
  //PA: also takes s now for Z3 preserving case with s as output, 
  //but here we set Z3 false anyway
  if(Z3){
    if (rewsbSvev(sgnMu,s)) flagMusqwrongsign(true);
    else flagMusqwrongsign(false);
    setSvev(s);
  }
  else{
    if (rewsbMu(sgnMu, mu)) flagMusqwrongsign(true);
    else flagMusqwrongsign(false);
    setSusyMu(mu);
  }
  if(Z3){ 
    if (rewsbKap(kap)) flagM3sq(true);  
    else flagM3sq(false); 
    setKappa(kap);
  }
  //PA:  again using rewsbM3sq which can work for Z3 violating case  
  else { 
    if (rewsbM3sq(mu, m3sq)) flagM3sq(true);  
    else flagM3sq(false); 
    setM3Squared(m3sq);
  }
  if(Z3 == false){
    rewsbXiS(xiS);
    setXiS(xiS);
  }
  else{
    rewsbmSsq(mSsq);
    setMsSquared(mSsq);
  }
  
  if ((displayMh1Squared() + 2.0 * sqr(displaySusyMu()) +
       displayMh2Squared() - 2.0 * fabs(displayM3Squared())) < 0.0 )
    flagHiggsufb(true);
  else 
    flagHiggsufb(false);
  return;
}


/// Obtains solution of one-loop effective potential minimisation via iteration
/// technique
/// err is 1 if no iteration reached
/// 2 if incorrect rewsb
void NmssmSoftsusy::iterateMu(double & muold, int sgnMu,
			     double mt, int maxTries, double pizzMS,
			     double sinthDRbar, double tol, int & err) {
  static int numTries = 0;
  static double munew = 0.;
  double lam = displayLambda();
  double s = displaySvev();
  if (numTries - 1 > maxTries) { 
    if (PRINTOUT) cout << "iterateMu reached maxtries\n"; 
    numTries = 0; munew = 0.0;
    err = 1; return;
  }
  /// How close to convergence are we?
  double c = 1.0 - minimum(fabs(muold), fabs(munew)) / 
    maximum(fabs(muold), fabs(munew));
  if (PRINTOUT > 2) cout << " diff=" << c;

  if (c < tol) { 
    muold = munew; //err = 0;
    numTries = 0; munew = 0.0;
    if (PRINTOUT > 2) cout << " mu converged\n";
    return; 
  }

  numTries = numTries + 1;
  muold = munew;
 
  double mH1sq = displayMh1Squared(), mH2sq = displayMh2Squared(), 
    tanb = displayTanb(), tanb2 = sqr(tanb);
  
  double treeLevelMusqeff = (mH1sq - mH2sq * tanb2) 
    / (tanb2 - 1.0) - 0.5 * sqr(displayMz());

  try {
    calcDrBarPars();
    double oneLoopMusqeff = treeLevelMusqeff;
    /// calculate the new one-loop tadpoles with old value of mu
    if (numRewsbLoops > 0) {
       doTadpoles(mt, sinthDRbar);

      oneLoopMusqeff = treeLevelMusqeff - 0.5 * pizzMS + 
	(displayTadpole2Ms() * sqr(tanb) - displayTadpole1Ms()) /
	(sqr(tanb) - 1.0); 
    }
     //PA real root test
     // bsqM4ac = 2.0*sqr(lam * s) - 4*0.5*sqr(lam * s) + 4 * oneLoopMusqeff;
     double bsqM4ac =  4 * oneLoopMusqeff;
     if(bsqM4ac < 0) {
        munew = - lam * s / root2; // take real part
        err = 2; 
      if (PRINTOUT > 1) cout << "no real roots for mu";
     }
     else{
        munew = - lam * s / root2 + 0.5 * sgnMu * sqrt(bsqM4ac);
     } 
       
     setSusyMu(munew); 
    double m3sqnew, kapnew;
  
    //PA:  using a rewsbM3sq for the Z3 violating case 
    if (rewsbM3sq(munew, m3sqnew) == 0) {
       flagM3sq(false);
    }
    else{
       flagM3sq(true);
    }
    setM3Squared(m3sqnew);
  
    double xiSnew, mSsqnew;

    rewsbXiS(xiSnew);
    setXiS(xiSnew);

 }

 
  catch(const char *a) {
    numTries = 0;
    throw a;
  }
  catch(const string &a) {
    numTries = 0;
    throw a;
  }
  
  if (PRINTOUT > 2) cout << " mu=" << munew;
  
  iterateMu(muold, sgnMu, mt, maxTries, pizzMS, sinthDRbar, tol, err);
  }


//Routine for iteratively solving for the singlet vev, s = <S>.
// where the EWSB is used to swap (kappa, mS) --> (mZ, tb)   
// and determine s.
void NmssmSoftsusy::iterateSvev(double & sold, int sgnMu,
			     double mt, int maxTries, double pizzMS,
			     double sinthDRbar, double tol, int & err) {
  static int numTries = 0;
  static double snew = 0.0;
  double lam = displayLambda();
  double s = displaySvev();
  if (numTries - 1 > maxTries) { 
    if (PRINTOUT) cout << "iterateMu reached maxtries\n"; 
    numTries = 0; snew = 0.0;
    err = 1; return;
  }
  /// How close to convergence are we?
  double c = 1.0 - minimum(fabs(sold), fabs(snew)) / 
    maximum(fabs(sold), fabs(snew));
  
  if (PRINTOUT > 2) cout << " diff=" << c;

  if (c < tol) { 
    sold = snew; //err = 0;
    numTries = 0; snew = 0.0;
    if (PRINTOUT > 2) cout << " mu converged\n";
    return; 
  }

  numTries = numTries + 1;
 
  sold = snew;
  
  double mH1sq = displayMh1Squared(), mH2sq = displayMh2Squared(), 
    tanb = displayTanb(), tanb2 = sqr(tanb);
  
  double treeLevelMusqeff = (mH1sq - mH2sq * tanb2) 
    / (tanb2 - 1.0) - 0.5 * sqr(displayMz());
  
  try {
   
    calcDrBarPars();

    double oneLoopMusqeff = treeLevelMusqeff;
    /// calculate the new one-loop tadpoles with old value of mu
    if (numRewsbLoops > 0) {
   
       doTadpoles(mt, sinthDRbar);

      oneLoopMusqeff = treeLevelMusqeff - 0.5 * pizzMS + 
	(displayTadpole2Ms() * sqr(tanb) - displayTadpole1Ms()) /
	(sqr(tanb) - 1.0); 
       }
  
    double ssq = 2.0 * oneLoopMusqeff / sqr(lam);
    snew = sgnMu * sqrt(fabs(ssq));
    if (oneLoopMusqeff < 0.0) {
      err = 2; 
      if (PRINTOUT > 1) cout << " mueff^2<0 ";
    }  
     setSvev(snew); 
   
     double m3sqnew, kapnew;
   
    //PA:  using a rewsKap which can work for Z3 symmetric version
    //with kappa as output  
    if (rewsbKap(kapnew) == 0) {
       flagM3sq(false);
    }
    else{
       flagM3sq(true);
    }
    setKappa(kapnew);
  
    double xiSnew, mSsqnew;

     rewsbmSsq(mSsqnew);
     setMsSquared(mSsqnew);
  

  }
  catch(const char *a) {
    numTries = 0;
    throw a;
  }
  catch(const string &a) {
    numTries = 0;
    throw a;
  }
  
  if (PRINTOUT > 2) cout << " singlet VEV =" << snew;
  
  iterateSvev(sold, sgnMu, mt, maxTries, pizzMS, sinthDRbar, tol, err);
  }


/// Organises rewsb: call it at the low scale MS^2 = sqrt(0.5 * (mT1^2 +
/// mT2^2)) is best, or below if it's decoupled from there. 
/// Call with zero, or no mt if you want tree level
void NmssmSoftsusy::rewsb(int sgnMu, double mt, double muOld, double eps) {
  double munew, m3sqnew, kapnew, snew;
  double xiSnew, mSsqnew;
  double sinthDRbarMS = calcSinthdrbar();
  
  calcTadpole1Ms1loop(mt, sinthDRbarMS);  
  calcTadpole2Ms1loop(mt, sinthDRbarMS); 
  calcTadpoleSMs1loop(mt, sinthDRbarMS); 
  
  munew = displaySusyMu();
  snew = displaySvev();
  
  /// Iterate to get a self-consistent mu solution
  int maxTries = 20, err = 0;
  double tol = TOLERANCE * 1.0e-4;
  
  double pizztMS = sqr(displayMzRun()) - sqr(displayMz()); /// resums logs
  //double pizztMS = piZZT(displayMz(), displayMu());
  
  if(Z3){
    iterateSvev(snew, sgnMu, mt, maxTries, pizztMS, sinthDRbarMS,
	    tol, err); 
    if (err == 2) flagMusqwrongsign(true);
    else flagMusqwrongsign(false); 
    if (err == 1) flagNoMuConvergence(true);
    else setSvev(snew);

  }

  else{
    iterateMu(munew, sgnMu, mt, maxTries, pizztMS, sinthDRbarMS,
	      tol, err); 
    
    if (err == 2) flagMusqwrongsign(true);
    else flagMusqwrongsign(false); 
    if (err == 1) flagNoMuConvergence(true);
    else setSusyMu(munew);
    
  
    /// average mu with the input value of muOld, if it isn't the number of the
    /// beast   
    if (muOld > -6.e66) {
      munew = (munew * eps + muOld * (1. - eps));
      setSusyMu(munew);
    }
    
  }
  //PA: using Z3 version of EWSB
  //with kappa as output
  if(Z3){
if (rewsbKap(kapnew) == 0) flagM3sq(false);
  else flagM3sq(true);   
  setKappa(kapnew);
  }
  else{ //PA: use rewsbM3sq which can work for Z3 violating case
  if (rewsbM3sq(munew, m3sqnew) == 0) flagM3sq(false);
  else flagM3sq(true);   
  setM3Squared(m3sqnew);
  }
  if(Z3 == false){
    rewsbXiS(xiSnew);
    setXiS(xiSnew);
  }
 else{
     rewsbmSsq(mSsqnew);
     setMsSquared(mSsqnew);
  }

  if ((displayMh1Squared() + 2.0 * sqr(displaySusyMu()) +
       displayMh2Squared() - 2.0 * fabs(displayM3Squared())) < 0.0 )
    flagHiggsufb(true);
  else 
    flagHiggsufb(false);
}



/// Organises calculation of physical quantities such as sparticle masses etc
/// Call AT MSusy
void NmssmSoftsusy::physical(int accuracy) {
  double sinthDRbarMS, piwwtMS, pizztMS;
 
  calcDrBarPars();
   
  if (accuracy == 0) {
    sinthDRbarMS = 0.0;
    piwwtMS = 0.0;
    pizztMS = 0.0;
  }
  else {
    sinthDRbarMS = calcSinthdrbar();
    piwwtMS = sqr(displayMwRun()) - sqr(displayMw());
    pizztMS = sqr(displayMzRun()) - sqr(displayMz());
  }

   /// Running masses at MSUSY
  double mt = displayDrBarPars().mt;
  double mb = displayDrBarPars().mb;
  double mtau = displayDrBarPars().mtau;

  /// Re-calculate the 1-loop tadpoles for the calculation
  calcTadpole1Ms1loop(mt, sinthDRbarMS);  
  calcTadpole2Ms1loop(mt, sinthDRbarMS);
  calcTadpoleSMs1loop(mt, sinthDRbarMS);

  /// Sfermion masses: all three families in each
  //PA: virtual methods diagonalisation carried out by
  //Softsusy class methods, tree level matrices filled by
  //NmssmSoftsusy methods
  doUpSquarks(mt, pizztMS, sinthDRbarMS, accuracy); 
  doDownSquarks(mb, pizztMS, sinthDRbarMS, accuracy, mt); 
  doChargedSleptons(mtau, pizztMS, sinthDRbarMS, accuracy); 
  doSnu(pizztMS, accuracy);

  //PA:Fill with current values of physpars, 
  //including those set by routines immediately above
  sPhysical phys(displayPhys());
  //PA: set up sPhysical object to be suitable for NMSSM
  phys.mh0.setEnd(3);
  phys.mA0.setEnd(2);
  phys.mixh0.resize(3,3);
  phys.mixA0.resize(2,2);
  phys.mneut.setEnd(5);
  phys.mixNeut.resize(5,5);

  NmssmSoftsusy * ppp;
  ppp = this;
  ppp->higgs(accuracy, piwwtMS, pizztMS, phys); 
  setPhys(phys);

  const int maxHiggsIterations = 20;
  double currentAccuracy = 1.0;
  DoubleVector oldHiggsMasses(6);
  oldHiggsMasses(1) = ppp->displayPhys().mh0(1);   
  oldHiggsMasses(2) = ppp->displayPhys().mh0(2);
  oldHiggsMasses(3) = ppp->displayPhys().mh0(3);
  oldHiggsMasses(4) = ppp->displayPhys().mA0(1);
  oldHiggsMasses(5) = ppp->displayPhys().mA0(2);  
  oldHiggsMasses(6) = ppp->displayPhys().mHpm;

  bool higgsTachyon = false;
  /// Iterate Higgs calculation (unless accuracy=0, in which case we just need
  /// a rough calculation) until the Higgs masses all converge to better than
  /// TOLERANCE fractional accuracy
  int i = 1; while (i < maxHiggsIterations && accuracy > 0 && 
		    currentAccuracy > TOLERANCE) {
     higgsTachyon = ppp->higgs(accuracy, piwwtMS, pizztMS, phys); /// iterate  
     setPhys(phys);
    DoubleVector newHiggsMasses(6);
    newHiggsMasses(1) = ppp->displayPhys().mh0(1);   
    newHiggsMasses(2) = ppp->displayPhys().mh0(2);
    newHiggsMasses(3) = ppp->displayPhys().mh0(3);
    newHiggsMasses(4) = ppp->displayPhys().mA0(1);
    newHiggsMasses(5) = ppp->displayPhys().mA0(2);  
    newHiggsMasses(6) = ppp->displayPhys().mHpm;

    currentAccuracy = oldHiggsMasses.compare(newHiggsMasses);

    oldHiggsMasses = newHiggsMasses;
    
    i++;
  }

  if (higgsTachyon) { flagTachyon(h0); flagTachyon(A0); flagTachyon(hpm); }

  phys.mh0(1) = ppp->displayPhys().mh0(1);
  phys.mh0(2) = ppp->displayPhys().mh0(2);
  phys.mh0(3) = ppp->displayPhys().mh0(3);
  phys.mA0(1) = ppp->displayPhys().mA0(1);
  phys.mA0(2) = ppp->displayPhys().mA0(2);
  phys.mHpm   = ppp->displayPhys().mHpm;

  charginos(accuracy, piwwtMS, phys); 
  neutralinos(accuracy, piwwtMS, pizztMS, phys);
  //PA: now set these values from NMSSM routines
  setPhys(phys);
  gluino(accuracy); 
  return;
}

//PA: Higgs routine for NMSSM
bool NmssmSoftsusy::higgs(int accuracy, double piwwtMS, double pizztMS, 
                          sPhysical & phys) {
  double tanb = displayTanb();
  double cb = cos(atan(tanb)), sb = sin(atan(tanb));     
  double ht = displayYukawaElement(YU, 3, 3);
  double hb = displayYukawaElement(YD, 3, 3);
  double htau = displayYukawaElement(YE, 3, 3);
  double gp = displayGaugeCoupling(1) * sqrt(0.6);
  double gs = displayGaugeCoupling(3);
  double lam = displayLambda();
  double svev = displaySvev();
  double kap = displayKappa();
  double mupr = displayMupr();
  double alam = displayTrialambda();
  double m3sq = displayM3Squared();
  double xiF = displayXiF();
  double vev = displayHvev(), vev2 = sqr(vev);
  double m3sqeff = m3sq  + lam * (mupr * svev / root2 + xiF)
    + alam * svev / root2  +   0.5 * lam * kap * sqr(svev);
  double MAeffsq =  m3sqeff / (sb * cb); ;
   //PA: initialise CP even mass matrix in (Hd, Hu, S) basis
   // CP odd Higgs mass matrices mPpr in (Hd, Hu, S) basis 
   //and mP2 in roatated basis (A, S) -- goldstone boson removed    
   DoubleMatrix mS(3,3), mPpr(3,3), mP2(2,2); 
   double mHpmsq, beta = atan(displayTanb()); 
   //PA: fill tree level CP even and CP odd Higgs mass matrices 
   //and tree level mHPm .
   treeHiggs(mS, mPpr, mP2, mHpmsq, beta);
   
   DoubleMatrix mhAtmH1(mS), mhAtmH2(mS), mhAtmH3(mS);
   DoubleMatrix sigmaMH1(3, 3), sigmaMH2(3, 3), sigmaMH3(3, 3);
   DoubleMatrix maAtmA1(mPpr), maAtmA2(mPpr);
   DoubleMatrix sigmaMA1(3, 3), sigmaMA2(3, 3);
   double q = displayMu(), p; 
   if(accuracy > 0){
     p = phys.mh0(1);    
     sigmaMH1(1, 1) = pis1s1(p, q);
     sigmaMH1(1, 2) = pis1s2(p, q);
     sigmaMH1(1, 3) = pis1s3(p, q);
     sigmaMH1(2, 2) = pis2s2(p, q);
     sigmaMH1(2, 3) = pis2s3(p, q);
     sigmaMH1(3, 3) = pis3s3(p, q);

     p = phys.mh0(2);		
     sigmaMH2(1, 1) = pis1s1(p, q);
     sigmaMH2(1, 2) = pis1s2(p, q);
     sigmaMH2(1, 3) = pis1s3(p, q);
     sigmaMH2(2, 2) = pis2s2(p, q);
     sigmaMH2(2, 3) = pis2s3(p, q);
     sigmaMH2(3, 3) = pis3s3(p, q);
      
     p = phys.mh0(3);    
     sigmaMH3(1, 1) = pis1s1(p, q);
     sigmaMH3(1, 2) = pis1s2(p, q);
     sigmaMH3(1, 3) = pis1s3(p, q);
     sigmaMH3(2, 2) = pis2s2(p, q);
     sigmaMH3(2, 3) = pis2s3(p, q);
     sigmaMH3(3, 3) = pis3s3(p, q);
      
     p = phys.mA0(1);
     sigmaMA1(1, 1) = pip1p1(p, q);
     sigmaMA1(1, 2) = pip1p2(p, q);
     sigmaMA1(1, 3) = pip1p3(p, q);
     sigmaMA1(2, 2) = pip2p2(p, q);
     sigmaMA1(2, 3) = pip2p3(p, q);
     sigmaMA1(3, 3) = pip3p3(p, q);
      
     p = phys.mA0(2);
     sigmaMA2(1, 1) = pip1p1(p, q);
     sigmaMA2(1, 2) = pip1p2(p, q);
     sigmaMA2(1, 3) = pip1p3(p, q);
     sigmaMA2(2, 2) = pip2p2(p, q);
     sigmaMA2(2, 3) = pip2p3(p, q);
     sigmaMA2(3, 3) = pip3p3(p, q);

     double d11 = 0.0;
     double d22 = 0.0;

     if(numHiggsMassLoops > 1) {
       
      
	 double mg = displayGaugino(3),  gs = displayGaugeCoupling(3), 
	   as = sqr(gs) / (4.0 * PI),  rmt = displayDrBarPars().mt, 
	   rmb = displayDrBarPars().mb, 
	   rmtsq = sqr(displayDrBarPars().mt),
	   rmbsq = sqr(displayDrBarPars().mb),
	   mst1sq = sqr(displayDrBarPars().mu(1, 3)), 
	   mst2sq = sqr(displayDrBarPars().mu(2, 3)),
	   msb1sq = sqr(displayDrBarPars().md(1, 3)), 
	   msb2sq = sqr(displayDrBarPars().md(2, 3)); 
	 double sxt = sin(displayDrBarPars().thetat), 
	   cxt = cos(displayDrBarPars().thetat),
	   sxb = sin(displayDrBarPars().thetab), 
	   cxb = cos(displayDrBarPars().thetab);
	 
	 double scalesq = sqr(displayMu()); 
       if(Z3){   
	 double DMS[3][3], DMP[3][3];
	 double DMSB[3][3], DMPB[3][3];
	 for(int i=0; i<=2; i++){
	   for(int j=0; j<=2; j++){
	     DMS[i][j] = 0;
	     DMP[i][j] = 0;
	     DMSB[i][j] = 0;
	     DMPB[i][j] = 0;
	   }
	 }
	 
	 double lamS = lam;
	 double vevS =  vev / root2;
	 double svevS = svev / root2;
	 
	 int loop = 2;
	 //PA: get alpha_s alpha_t pieces
	 effpot_(&loop, &rmt, &mg, &mst1sq, &mst2sq, &sxt, &cxt,
		 &scalesq, &tanb, &vevS, &lamS, &svevS, &as, &DMS, &DMP);
	 //PA: get alpha_s alpha_b pieces
	 double cotb = 1.0 / tanb;
	 effpot_(&loop, &rmb, &mg, &msb1sq, &msb2sq, &sxb, &cxb,
		 &scalesq, &cotb, &vevS, &lamS, &svevS, &as, &DMSB, &DMPB);
	 
	 //PA: Make appropriate substitutions for elements following 0907.4682
	 // bottom of page 9
	 
	 double temp = DMSB[0][0];
	 DMSB[0][0] = DMSB[1][1];
	 DMSB[1][1] = temp;
	 temp = DMSB[0][2];
	 DMSB[0][2] = DMSB[1][2];
	 DMSB[1][2] = temp;


	 for(int i=0; i<=2; i++){
	   for(int j=0; j<=2; j++){
	     DMS[i][j] = DMS[i][j] + DMSB[i][j];
	   }
	 }
	 double amu = - lam * svev / root2;
	 int kkk = 0.0;
	 double s11s = 0.0, s12s = 0.0, s22s = 0.0;
	 double s11b = 0.0, s12b = 0.0, s22b = 0.0;
	 double s11w = 0.0, s12w = 0.0, s22w = 0.0;
	 double s11tau = 0.0, s12tau = 0.0, s22tau = 0.0;
	 double ps2 = 0.0, p2b = 0.0, p2w = 0.0, p2tau = 0.0;
	 double fmasq = fabs(MAeffsq);
	 dszhiggs_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, 
		   &amu, &tanb, &vev2, &gs, &kkk, &s11s, &s22s, &s12s);
	 dszodd_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu,
		 &tanb, &vev2, &gs, &ps2); 
	 dszhiggs_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq,
		   &amu, &cotb, &vev2, &gs, &kkk, &s22b, &s11b, &s12b);
	 dszodd_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq, &amu,
		 &cotb, &vev2, &gs, &p2b);
	 //Corrections as in MSSM, not corrected for NMSSM,
	 //Should be OK for MSSM states when S state is close to decoupled 
	 ddshiggs_(&rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, 
		   &msb2sq, &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, 
		   &vev2, &s11w, &s12w, &s22w);
	 ddsodd_(&rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq, 
		 &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &p2w);
	 
	 double sintau = sin(displayDrBarPars().thetatau),
	   costau = cos(displayDrBarPars().thetatau);
	 double rmtausq = sqr(displayDrBarPars().mtau);
	 int OS = 0;
	 double mstau1sq = sqr(displayDrBarPars().me(1, 3)), 
	   mstau2sq = sqr(displayDrBarPars().me(2, 3));
	 double msnusq = sqr(displayDrBarPars().msnu(3));
	 tausqhiggs_(&rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
		     &costau, &scalesq, &amu, &tanb, &vev2, &OS, &s11tau, 
		     &s22tau, &s12tau);
	 tausqodd_(&rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
		   &costau, &scalesq, &amu, &tanb, &vev2, &p2tau);
       
   	 //PA: Add O(y_t^4 + y_t^2y_b^2 + y_b^4) from mssm two loop parts
	 DMS[0][0] = DMS[0][0] + s11w + sqr(sb) * p2w;
	 DMS[0][1] = DMS[0][1] + s12w - sb * cb * p2w;
	 DMS[1][1] = DMS[1][1] + s22w + sqr(cb) * p2w;
	 //PA: Add alpha_tau^2 and alpha_tau alpha_b from mssm two loop parts  
	 DMS[0][0] = DMS[0][0] + s11tau + sqr(sb) * p2tau;
	 DMS[0][1] = DMS[0][1] + s12tau - sb * cb * p2tau;
	 DMS[1][1] = DMS[1][1] + s22tau + sqr(cb) * p2tau;
	 //PA: check against twoloophiggs.f pieces
	 double dMA = ps2 + p2w + p2b + p2tau;
	 double DMS00 = s11s + s11w + s11b + s11tau + dMA * sqr(sb);
	 double DMS01 = s12s + s12w + s12b + s12tau - dMA * cb * sb;
	 double DMS11 = s22s + s22w + s22b + s22tau + dMA * sqr(cb);
	
	 //PA: Now add two loop parts to the full one loop self energy
	 sigmaMH1(1, 1) = sigmaMH1(1, 1) - DMS[0][0]; 
	 sigmaMH1(1, 2) = sigmaMH1(1, 2) - DMS[0][1]; 
	 sigmaMH1(2, 2) = sigmaMH1(2, 2) - DMS[1][1]; 
	 sigmaMH1(1, 3) = sigmaMH1(1, 3) - DMS[0][2];    
	 sigmaMH1(2, 3) = sigmaMH1(2, 3) - DMS[1][2]; 
	 sigmaMH1(3, 3) = sigmaMH1(3, 3) - DMS[2][2]; 
	 
	 sigmaMH2(1, 1) = sigmaMH2(1, 1) - DMS[0][0]; 
	 sigmaMH2(1, 2) = sigmaMH2(1, 2) - DMS[0][1]; 
	 sigmaMH2(2, 2) = sigmaMH2(2, 2) - DMS[1][1]; 
	 sigmaMH2(1, 3) = sigmaMH2(1, 3) - DMS[0][2];    
	 sigmaMH2(2, 3) = sigmaMH2(2, 3) - DMS[1][2]; 
	 sigmaMH2(3, 3) = sigmaMH2(3, 3) - DMS[2][2]; 
         
	 
	 sigmaMH3(1, 1) = sigmaMH3(1, 1) - DMS[0][0]; 
	 sigmaMH3(1, 2) = sigmaMH3(1, 2) - DMS[0][1]; 
	 sigmaMH3(2, 2) = sigmaMH3(2, 2) - DMS[1][1]; 
	 sigmaMH3(1, 3) = sigmaMH3(1, 3) - DMS[0][2];    
	 sigmaMH3(2, 3) = sigmaMH3(2, 3) - DMS[1][2]; 
	 sigmaMH3(3, 3) = sigmaMH3(3, 3) - DMS[2][2]; 
       }
      
       
     }
     
     sigmaMH1.symmetrise();
     sigmaMH2.symmetrise();
     sigmaMH3.symmetrise();
     sigmaMA1.symmetrise();
     sigmaMA2.symmetrise();
     //PA: adding *one loop* tadpoles
     mhAtmH1(1, 1) =  mhAtmH1(1, 1) + displayTadpole1Ms1loop();
     mhAtmH2(1, 1) =  mhAtmH2(1, 1) + displayTadpole1Ms1loop();
     mhAtmH3(1, 1) =  mhAtmH3(1, 1) + displayTadpole1Ms1loop();

     mhAtmH1(2, 2) = mhAtmH1(2, 2) + displayTadpole2Ms1loop();
     mhAtmH2(2, 2) = mhAtmH2(2, 2) + displayTadpole2Ms1loop();
     mhAtmH3(2, 2) = mhAtmH3(2, 2) + displayTadpole2Ms1loop();
     
     mhAtmH1(3, 3) = mhAtmH1(3, 3) + displayTadpoleSMs1loop();
     mhAtmH2(3, 3) = mhAtmH1(3, 3) + displayTadpoleSMs1loop();
     mhAtmH3(3, 3) = mhAtmH1(3, 3) + displayTadpoleSMs1loop();
     
     mhAtmH1 = mhAtmH1 - sigmaMH1;
     mhAtmH2 = mhAtmH2 - sigmaMH2;
     mhAtmH3 = mhAtmH3 - sigmaMH3;
     
     //PA: adding *one loop* tadpoles
     maAtmA1(1, 1) = mPpr(1, 1) + displayTadpole1Ms1loop();
     maAtmA1(1, 2) = mPpr(1, 2);
     maAtmA1(1, 3) = mPpr(1, 3);
     maAtmA1(2, 1) = maAtmA1(1, 2);
     maAtmA1(2, 2) = mPpr(2, 2) + displayTadpole2Ms1loop();
     maAtmA1(2, 3) = mPpr(2, 3);
     maAtmA1(3, 1) = maAtmA1(1, 3);
     maAtmA1(3, 2) = maAtmA1(2, 3);
     maAtmA1(3, 3) = maAtmA1(3, 3) + displayTadpoleSMs1loop();
     maAtmA1 = maAtmA1 - sigmaMA1;
     
     maAtmA2(1, 1) = mPpr(1, 1) + displayTadpole1Ms1loop();
     maAtmA2(1, 2) = mPpr(1, 2);
     maAtmA2(1, 3) = mPpr(1, 3);
     maAtmA2(2, 1) = maAtmA2(1, 2);
     maAtmA2(2, 2) = mPpr(2, 2) + displayTadpole2Ms1loop();
     maAtmA2(2, 3) = mPpr(2, 3);
     maAtmA2(3, 1) = maAtmA2(1, 3);
     maAtmA2(3, 2) = maAtmA2(2, 3);
     maAtmA2(3, 3) = maAtmA2(3, 3) + displayTadpoleSMs1loop();
     maAtmA2 = maAtmA2 - sigmaMA2;
   }
   
  /// LCT: NMSSM Higgs states.  CP-even first
  DoubleVector temp(3);
  DoubleMatrix mixHiggsLoops(3, 3);
  
  if (mhAtmH1.diagonaliseSym(mixHiggsLoops, temp) > TOLERANCE *
      1.0e-3) { 
    ostringstream ii;
    ii << "accuracy bad in CP-even Higgs diagonalisation"<< flush;
    throw ii.str(); 
	}

  // cout << "In NM higgs accuracy = " << accuracy << endl;
  // cout << "p=mH1 temp: " << temp << endl;
  // cout << "mhAtmH1: " << mhAtmH1 << endl;

 bool h0Htachyon = false;
  if (temp(1) < 0.0 || temp(2) < 0.0 || temp(3) < 0.0) {
    h0Htachyon = true;
    if (PRINTOUT > 2) cout << "H1/H2/H3 tachyon: m^2 = " << temp;
  }
  
  temp = temp.apply(ccbSqrt);
  phys.mixh0 = mixHiggsLoops;
  phys.mh0(1) = temp(1);
  
  /// LCT: Now repeat for p = mH2.  Must reset CP-even mixing matrix
  for (int i=1; i <= 3; i++) {
    for (int j=1; j <= 3; j++) {
      mixHiggsLoops(i, j) = 0.0;
    }
  }
 
  if (mhAtmH2.diagonaliseSym(mixHiggsLoops, temp) > TOLERANCE *
      1.0e-3) {
    ostringstream ii;
    ii << "accuracy bad in CP-even Higgs diagonalisation"<< flush;
    throw ii.str(); 
	}

  // cout << "p=mH2 temp: " << temp << endl;
 
  if (temp(1) < 0.0 || temp(2) < 0.0 || temp(3) < 0.0) {
    h0Htachyon = true;
    if (PRINTOUT > 2) cout << "H1/H2/H3 tachyon: m^2 = " << temp;
  }
  
  temp = temp.apply(ccbSqrt);
  phys.mh0(2) = temp(2);
  
  /// LCT: Now repeat for p = mH3.  Must reset CP-even mixing matrix
  for (int i=1; i <= 3; i++) {
    for (int j=1; j <= 3; j++) {
      mixHiggsLoops(i, j) = 0.0;
    }
  }
  
  if (mhAtmH3.diagonaliseSym(mixHiggsLoops, temp) > TOLERANCE *
      1.0e-3) { 
    ostringstream ii;
    ii << "accuracy bad in CP-even Higgs diagonalisation"<< flush;
    throw ii.str(); 
	}
  
  // cout << "p=mH3 temp: " << temp << endl;

  if (temp(1) < 0.0 || temp(2) < 0.0 || temp(3) < 0.0) {
    h0Htachyon = true;
    if (PRINTOUT > 2) cout << "H1/H2/H3 tachyon: m^2 = " << temp;
  }
  
  temp = temp.apply(ccbSqrt);
  phys.mh0(3) = temp(3);
  
  /// LCT: CP-odd states
  DoubleMatrix mP(3, 3);
  /// LCT: Rotate CP-odd mass^2 matrix into (G, A, S_I) basis
  mP = rot3d(beta).transpose() * maAtmA1 * rot3d(beta);	
  /// LCT: Drop Goldstone from 3 x 3 CP-odd Higgs mass^2 matrix and 
  /// construct 2 x 2 matrix in (A, S_I) basis
  mP2(1, 1) = mP(2, 2);
  mP2(1, 2) = mP(2, 3);
  mP2(2, 1) = mP(3, 2);
  mP2(2, 2) = mP(3, 3);

  /// LCT: Diagonalise	
  DoubleVector Atemp(2);
  double Atheta;
  
  Atemp = mP2.sym2by2(Atheta);

  if (Atemp(1) < 0.0 && Atemp(2) < 0.0) {
     h0Htachyon = true;
     if (PRINTOUT > 2) cout << " A1/A2 tachyon: m^2=" << Atemp;
  }
  Atemp = Atemp.apply(zeroSqrt);
  
  if (Atemp(1) > Atemp(2)) Atheta = Atheta + PI * 0.5; 
  
  phys.thetaH = Atheta; /// Atheta defined for p=mA1  
  int j; double mA1 = Atemp.apply(fabs).min(j);
  
  /// LCT: Rotate CP-odd mass^2 matrix into (G, A, S_I) basis
  mP = rot3d(beta).transpose() * maAtmA2 * rot3d(beta);	
  /// LCT: Drop Goldstone from 3 x 3 CP-odd Higgs mass^2 matrix and 
  /// construct 2 x 2 matrix in (A, S_I) basis
  mP2(1, 1) = mP(2, 2);
  mP2(1, 2) = mP(2, 3);
  mP2(2, 1) = mP(3, 2);
  mP2(2, 2) = mP(3, 3);


  Atemp = mP2.sym2by2(Atheta);

  if (Atemp(1) < 0.0 && Atemp(2) < 0.0) {
    h0Htachyon = true;
    if (PRINTOUT > 2) cout << " A1/A2 tachyon: m^2=" << Atemp;
  }
  Atemp = Atemp.apply(zeroSqrt);
  double mA2 = Atemp.max();
  phys.mA0(1) = mA1;
  phys.mA0(2) = mA2;

  /// LCT: Charged Higgs pole mass accurate to 1-loop (for now)
  double pihphm = piHpHm(phys.mHpm, displayMu());
  
  double mzPole = displayMz(), mzRun2 = sqr(displayMzRun());
  double poleMhcSq = mHpmsq - sqr(displayMwRun())+ sqr(displayMw()) ;
 
 if (accuracy > 0) {
   poleMhcSq = mHpmsq - pihphm - sqr(displayMwRun()) + sqr(displayMw()) 
     + piwwtMS + sqr(sin(beta)) * displayTadpole1Ms1loop() + sqr(cos(beta)) *
     displayTadpole2Ms1loop();
  
 }
 
  phys.mHpm = zeroSqrt(poleMhcSq);
if (poleMhcSq > 0. && !h0Htachyon) return false;
  else {
     if (PRINTOUT) cout << " mHc(phys)^2=" << poleMhcSq 
		        << " but may be first iteration" << endl;
    return true;
  }
}


//PA: adds gauge boson contribitions to the left right and scalar parts 
//of the chargino self energy
void NmssmSoftsusy::addChaLoopHiggs(double p, DoubleMatrix & sigmaL, DoubleMatrix & sigmaR, DoubleMatrix & sigmaS, DoubleMatrix b1pCha, DoubleMatrix b0pCha, DoubleMatrix b1pNeut, DoubleMatrix b0pNeut) const {
  double g = displayGaugeCoupling(2), 
    gp = displayGaugeCoupling(1) * sqrt(0.6), 
    tanb = displayTanb();
  double beta = atan(tanb);
  double sinb = sin(beta), cosb = cos(beta);
  double e = g * calcSinthdrbar(); 
  double lam = displayLambda();

  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 
  const int dimN =  mneut.displayEnd();
  /// checked and corrected
  DoubleMatrix aPsiPsiHc1(dimN, 2), bPsiPsiHc1(dimN, 2);
  DoubleMatrix aPsiPsiHc2(dimN, 2), bPsiPsiHc2(dimN, 2);
  aPsiPsiHc1(1, 2) = gp / root2;
  bPsiPsiHc2(1, 2) = aPsiPsiHc1(1, 2);
  aPsiPsiHc1(2, 2) = g / root2;
  bPsiPsiHc2(2, 2) = g / root2;
  aPsiPsiHc1(3, 1) = -g;
  bPsiPsiHc2(4, 1) = g;
  
  aPsiPsiHc2(5, 2) = lam;
  bPsiPsiHc1(5, 2) = - lam;
  
  ComplexMatrix aChiPsiHc1(dimN, 2), bChiPsiHc1(dimN, 2);
  ComplexMatrix aChiPsiHc2(dimN, 2), bChiPsiHc2(dimN, 2);
  aChiPsiHc1 = n.complexConjugate() * bPsiPsiHc1;
  bChiPsiHc1 = n * aPsiPsiHc1;
  aChiPsiHc2 = n.complexConjugate() * bPsiPsiHc2;
  bChiPsiHc2 = n * aPsiPsiHc2;
  ComplexMatrix aChiPsiHHp(dimN, 2), bChiPsiHHp(dimN, 2);
  ComplexMatrix aChiPsiHGp(dimN, 2), bChiPsiHGp(dimN, 2);
  for (int i=1; i<=dimN; i++)
    for (int j=1; j<=2; j++) {
      aChiPsiHGp(i, j) = cosb * aChiPsiHc1(i, j) + sinb * aChiPsiHc2(i, j);
      bChiPsiHGp(i, j) = cosb * bChiPsiHc1(i, j) + sinb * bChiPsiHc2(i, j);
      aChiPsiHHp(i, j) =-sinb * aChiPsiHc1(i, j) + cosb * aChiPsiHc2(i, j);
      bChiPsiHHp(i, j) =-sinb * bChiPsiHc1(i, j) + cosb * bChiPsiHc2(i, j);
    }

  /// checked this block
  ComplexMatrix aPsiPsis1(2, 2), aPsiPsis2(2, 2), 
    aPsiPsip1(2, 2), aPsiPsip2(2, 2),
    aPsiPsis3(2, 2), aPsiPsip3(2, 2);
  ComplexMatrix bPsiPsis1(2, 2), bPsiPsis2(2, 2), 
    bPsiPsip1(2, 2), bPsiPsip2(2, 2), bPsiPsis3(2, 2), bPsiPsip3(2, 2);
  aPsiPsis1(1, 2) = g / root2;
  aPsiPsis2(2, 1) = g / root2;
  aPsiPsip1(1, 2) = g / root2;
  aPsiPsip2(2, 1) = - g / root2;
  
  aPsiPsis3(2, 2) =  lam / root2; 
  aPsiPsip3(2, 2) = - lam / root2; 
  /// checked and corrected 2/12/08
  bPsiPsis1 = aPsiPsis1.transpose(); 
  bPsiPsis2 = aPsiPsis2.transpose();
  bPsiPsis3 = aPsiPsis3.transpose();
  /// end of correction 2/12/08
  bPsiPsip1 = - 1.0 * aPsiPsip1.transpose();
  bPsiPsip2 = - 1.0 * aPsiPsip2.transpose();
  bPsiPsip3 = - 1.0 * aPsiPsip3.transpose();
  ComplexMatrix aPsiChis1(2, 2), aPsiChis2(2, 2), 
    aPsiChip1(2, 2), aPsiChip2(2, 2),
    aPsiChis3(2, 2),  aPsiChip3(2, 2);
  ComplexMatrix bPsiChis1(2, 2), bPsiChis2(2, 2), 
    bPsiChip1(2, 2), bPsiChip2(2, 2),
    bPsiChis3(2, 2),  bPsiChip3(2, 2);
  aPsiChis1 = aPsiPsis1 * u.hermitianConjugate();
  aPsiChis2 = aPsiPsis2 * u.hermitianConjugate();
  aPsiChis3 = aPsiPsis3 * u.hermitianConjugate();
  
  aPsiChip1 = aPsiPsip1 * u.hermitianConjugate();
  aPsiChip2 = aPsiPsip2 * u.hermitianConjugate();
  aPsiChip3 = aPsiPsip3 * u.hermitianConjugate();
  
  bPsiChis1 = bPsiPsis1 * v.transpose();
  bPsiChis2 = bPsiPsis2 * v.transpose();
  bPsiChis3 = bPsiPsis3 * v.transpose();

  bPsiChip1 = bPsiPsip1 * v.transpose();
  bPsiChip2 = bPsiPsip2 * v.transpose();
  bPsiChip3 = bPsiPsip3 * v.transpose();

  ComplexMatrix aPsiChih1(2, 2), aPsiChih2(2, 2), aPsiChih3(2, 2);
  ComplexMatrix aPsiChiG(2, 2), aPsiChiA1(2, 2), aPsiChiA2(2, 2);
  ComplexMatrix bPsiChih1(2, 2), bPsiChih2(2, 2), bPsiChih3(2, 2); 
  ComplexMatrix bPsiChiG(2, 2), bPsiChiA1(2, 2), bPsiChiA2(2, 2);

  /// LCT: Higgs mixing matrices
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb;
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
  aPsiChih1 = S(1,1) * aPsiChis1 + S(1,2) * aPsiChis2 + S(1,3) * aPsiChis3;
  aPsiChih2 = S(2,1) * aPsiChis1 + S(2,2) * aPsiChis2 + S(2,3) * aPsiChis3;
  aPsiChih3 = S(3,1) * aPsiChis1 + S(3,2) * aPsiChis2 + S(3,3) * aPsiChis3;
  //PA: mixing in slavich conventions, rules should be altered too
  aPsiChiG  = - P(1,1) * aPsiChip1 + P(1,2) * aPsiChip2 + P(1,3) * aPsiChip3;
  aPsiChiA1 =   P(2,1) * aPsiChip1 - P(2,2) * aPsiChip2 + P(2,3) * aPsiChip3;
  aPsiChiA2 =   P(3,1) * aPsiChip1 - P(3,2) * aPsiChip2 + P(3,3) * aPsiChip3;
  
  bPsiChih1 = S(1,1) * bPsiChis1 + S(1,2) * bPsiChis2 + S(1,3) * bPsiChis3;
  bPsiChih2 = S(2,1) * bPsiChis1 + S(2,2) * bPsiChis2 + S(2,3) * bPsiChis3;
  bPsiChih3 = S(3,1) * bPsiChis1 + S(3,2) * bPsiChis2 + S(3,3) * bPsiChis3;
  //PA: mixing in slavich conventions, rules should be altered too
  bPsiChiG  = - P(1,1) * bPsiChip1 + P(1,2) * bPsiChip2 + P(1,3) * bPsiChip3;
  bPsiChiA1 =   P(2,1) * bPsiChip1 - P(2,2) * bPsiChip2 + P(2,3) * bPsiChip3;
  bPsiChiA2 =   P(3,1) * bPsiChip1 - P(3,2) * bPsiChip2 + P(3,3) * bPsiChip3;
   
  /// checked and corrected  
  for (int i=1; i<=2; i++) 
    for (int j=1; j<=2; j++) 
      for (int k=1; k<=dimN; k++) {
	//G+ 
	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	  (aChiPsiHGp(k, i).conj() * aChiPsiHGp(k, j) * b1pNeut(k,1)).real();
	sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	  (bChiPsiHGp(k, i).conj() * bChiPsiHGp(k, j) * b1pNeut(k,1)).real();
	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
	  (bChiPsiHGp(k, i).conj() * aChiPsiHGp(k, j) * b0pNeut(k,1)).real();

	/// H+	
	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	  (aChiPsiHHp(k, i).conj() * aChiPsiHHp(k, j) * b1pNeut(k,2)).real();
	sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	  (bChiPsiHHp(k, i).conj() * bChiPsiHHp(k, j) * b1pNeut(k,2)).real();
	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
	  (bChiPsiHHp(k, i).conj() * aChiPsiHHp(k, j) * b0pNeut(k,2)).real();

	if (k <= 2) {
	 
	  /// h1
	  sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	    (aPsiChih1(i, k).conj() * aPsiChih1(j, k) * b1pCha(k,5)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	    (bPsiChih1(i, k).conj() * bPsiChih1(j, k) * b1pCha(k,5)).real();
	  sigmaS(i, j) = sigmaS(i, j) + mch(k) * 
	    (bPsiChih1(i, k).conj() * aPsiChih1(j, k) * b0pCha(k,5)).real();
  
	  /// h2
	  sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	    (aPsiChih2(i, k).conj() * aPsiChih2(j, k) * b1pCha(k,6)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	    (bPsiChih2(i, k).conj() * bPsiChih2(j, k) * b1pCha(k,6)).real();
	  sigmaS(i, j) = sigmaS(i, j) + mch(k) * 
	    (bPsiChih2(i, k).conj() * aPsiChih2(j, k) * b0pCha(k,6)).real();
	  	  
	   /// h3
	  sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	    (aPsiChih3(i, k).conj() * aPsiChih3(j, k) * b1pCha(k,7)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	    (bPsiChih3(i, k).conj() * bPsiChih3(j, k) * b1pCha(k,7)).real();
	  sigmaS(i, j) = sigmaS(i, j) + mch(k) * 
	   (bPsiChih3(i, k).conj() * aPsiChih3(j, k) * b0pCha(k,7)).real();
    

	  /// G0
	  sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	    (aPsiChiG(i, k).conj() * aPsiChiG(j, k) * b1pCha(k,2)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	    (bPsiChiG(i, k).conj() * bPsiChiG(j, k) * b1pCha(k,2)).real();
	  sigmaS(i, j) = sigmaS(i, j) + mch(k) * 
	    (bPsiChiG(i, k).conj() * aPsiChiG(j, k) * b0pCha(k,2)).real();
	
	  /// A1
	  sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	    (aPsiChiA1(i, k).conj() * aPsiChiA1(j, k) * b1pCha(k,3)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	    (bPsiChiA1(i, k).conj() * bPsiChiA1(j, k) * b1pCha(k,3)).real();
	  sigmaS(i, j) = sigmaS(i, j) + mch(k) * 
	    (bPsiChiA1(i, k).conj() * aPsiChiA1(j, k) * b0pCha(k,3)).real();

	  //A2
	  sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	    (aPsiChiA2(i, k).conj() * aPsiChiA2(j, k) * b1pCha(k,4)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	    (bPsiChiA2(i, k).conj() * bPsiChiA2(j, k) * b1pCha(k,4)).real();
	  sigmaS(i, j) = sigmaS(i, j) + mch(k) * 
	    (bPsiChiA2(i, k).conj() * aPsiChiA2(j, k) * b0pCha(k,4)).real();
	}
      }
  return;
}


void NmssmSoftsusy::addCharginoLoop(double p, DoubleMatrix & mass) {
  double mz = displayMzRun(), mw = displayMwRun(), q = displayMu();
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 
  const int dimN =  mneut.displayEnd();
  DoubleMatrix sigmaL(2, 2), sigmaR(2, 2), sigmaS(2, 2);
  //PA: calculate P-V's for gauge and Higgs at outset 
  //avoids repeated calls for same function
  //1st indice runs over neutralino or charginos
  //2nd runs over higgs and gaige bosons in bases:
  //for Cha - (Z,A1,A1,h1,h2,h3,gamma)
  //for Neut - (W, Hpm) 
  DoubleMatrix b1pCha(2,7), b0pCha(2,7), b1pNeut(dimN,2), b0pNeut(dimN,2); 
  for (int k=1; k<=2; k++) {
    b1pCha(k,1) = b1(p, mch(k), 0., q);	
    b0pCha(k,1) = b0(p, mch(k), 0., q);	
    b1pCha(k,2) = b1(p, mch(k), mz, q);
    b0pCha(k,2) = b0(p, mch(k), mz, q);

    b1pCha(k,3) = b1(p, mch(k), displayDrBarPars().mA0(1), q);
    b0pCha(k,3) = b0(p, mch(k), displayDrBarPars().mA0(1), q);
    b1pCha(k,4) = b1(p, mch(k), displayDrBarPars().mA0(2), q);
    b0pCha(k,4) = b0(p, mch(k), displayDrBarPars().mA0(2), q);
   
    b1pCha(k,5) = b1(p, mch(k), displayDrBarPars().mh0(1), q);
    b0pCha(k,5) = b0(p, mch(k), displayDrBarPars().mh0(1), q);
    b1pCha(k,6) = b1(p, mch(k), displayDrBarPars().mh0(2), q);
    b0pCha(k,6) = b0(p, mch(k), displayDrBarPars().mh0(2), q);
    b1pCha(k,7) = b1(p, mch(k), displayDrBarPars().mh0(3), q);
    b0pCha(k,7) = b0(p, mch(k), displayDrBarPars().mh0(3), q);

}
  
  for (int k=1; k<=dimN; k++) {
    b1pNeut(k,1) = b1(p, mneut(k), mw, q);
    b0pNeut(k,1) = b0(p, mneut(k), mw, q);
    b1pNeut(k,2) = b1(p, mneut(k), displayDrBarPars().mHpm, q);
    b0pNeut(k,2) = b0(p, mneut(k), displayDrBarPars().mHpm, q);
  }

  //PA: sfermion contributions to sigmaL, sigmaR and sigmaS
  //calls MSSM method, no NMSSM alterations
  addChaLoopSfermion(p, sigmaL, sigmaR, sigmaS);
  //Gauge bosons calls MSSM method, no NMSSM alterations!
  addChaLoopGauge(p, sigmaL, sigmaR, sigmaS, b1pCha, b0pCha, b1pNeut, b0pNeut);
  //Higgs bosons, calls NMSSM routine
  addChaLoopHiggs(p, sigmaL, sigmaR, sigmaS, b1pCha, b0pCha, b1pNeut, b0pNeut);

  mass = mass - 1.0 / (16.0 * sqr(PI)) * 
    (sigmaR * mass + mass * sigmaL + sigmaS);
}


void NmssmSoftsusy::charginos(int accuracy, double piwwtMS, sPhysical & phys) {
  DoubleMatrix mCh(2, 2);
  double mw = sqrt( sqr(displayMw()) + piwwtMS) ;   
  double beta = atan(displayTanb());
  treeCharginos(mCh, beta, mw);
  if (accuracy == 0) {
    phys.mch = mCh.asy2by2(phys.thetaL, phys.thetaR);
    return;
  }

  DoubleMatrix mCh2(mCh);

  double p1 = fabs(displayDrBarPars().mch(1)), 
    p2 = fabs(displayDrBarPars().mch(2));
  addCharginoLoop(p1, mCh);
  addCharginoLoop(p2, mCh2);

  double x = 0., y = 0.;
  DoubleVector mch1(mCh.asy2by2(phys.thetaL, phys.thetaR));
  DoubleVector mch2(mCh2.asy2by2(x, y));
  phys.mch(1) = mch1(1);
  /// You should take the sign of the chargino mass to be the same as
  /// got from the chargino_1 determination. Otherwise, if there's a
  /// difference, this will screw things up...
  double sgn_mass = mch1(2) / abs(mch1(2));
  phys.mch(2) = sgn_mass * abs(mch2(2));
}
void NmssmSoftsusy::neutralinos(int accuracy, double piwwtMS, double pizztMS, sPhysical & phys) {
   const drBarPars& forLoops = displayDrBarPars();
   double mw = displayMwRun();
   double mz = displayMzRun();
   double beta = atan(displayTanb());
   double sinth = calcSinthdrbar();
   DoubleMatrix mNeut(5, 5);
   treeNeutralinos(mNeut, beta, mz, mw, sinth);

   if (accuracy == 0) {
      mNeut.diagonaliseSym(phys.mixNeut, phys.mneut);
      return;
   }

   DoubleMatrix mNeut5(mNeut), mNeut4(mNeut), mNeut2(mNeut), mNeut3(mNeut);

  addNeutralinoLoop(fabs(forLoops.mneut(1)), mNeut);
  addNeutralinoLoop(fabs(forLoops.mneut(2)), mNeut2);
  addNeutralinoLoop(fabs(forLoops.mneut(3)), mNeut3);
  addNeutralinoLoop(fabs(forLoops.mneut(4)), mNeut4);
  addNeutralinoLoop(fabs(forLoops.mneut(5)), mNeut5);

  DoubleVector mneut(5), mneut2(5), mneut3(5), mneut4(5), mneut5(5);

  DoubleMatrix dummyMix(5, 5);
  double acceptableTol = TOLERANCE * 1.0e-3; 

  if (mNeut.diagonaliseSym(phys.mixNeut, mneut) > acceptableTol ||
      mNeut2.diagonaliseSym(dummyMix, mneut2) > acceptableTol ||
      mNeut3.diagonaliseSym(dummyMix, mneut3) > acceptableTol ||
       mNeut4.diagonaliseSym(dummyMix, mneut4) > acceptableTol || 
    mNeut5.diagonaliseSym(dummyMix, mneut5) > acceptableTol) { 
    ostringstream ii;
    ii << "accuracy bad in neutralino diagonalisation"<< flush;
    ii << "diagonalising " << phys.mneut << " with "   
       << phys.mixNeut;
    throw ii.str(); 
  }

  /// We should choose sign conventions from the case where the mixing is
  /// defined, in case there is a difference 
  phys.mneut(1) = mneut(1); 
  phys.mneut(2) = mneut(2) / abs(mneut(2)) * abs(mneut2(2));
  phys.mneut(3) = mneut(3) / abs(mneut(3)) * abs(mneut3(3)); 
  phys.mneut(4) = mneut(4) / abs(mneut(4)) * abs(mneut4(4));
  phys.mneut(5) = mneut(5) / abs(mneut(5)) * abs(mneut5(5));
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
  return nmHiggs;
}

double NmssmSoftsusy::piWWTgauginos(double p, double q, double thetaWDRbar) const {
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

  double thetaWDRbar = asin(calcSinthdrbar());
  double g = displayGaugeCoupling(2);

  double ans = 0.0;
  double higgs = piWWTHiggs(p, q, thetaWDRbar);
  double fermions = piWWTfermions(p, q, usePoleMt);   
  double sfermions = piWWTsfermions(p, q);   
  double gauginos = piWWTgauginos(p, q, thetaWDRbar);
  ans = higgs + sfermions + fermions + gauginos;
  double pi = ans * sqr(g) / (16.0 * sqr(PI));

  return pi;
}


//PA: pseudoscalar self energies in basis Im(H_d), Im(H_u), Im(S).
double NmssmSoftsusy::pip1p1(double p, double q) const {
  const drBarPars& tree = displayDrBarPars();
	
  double beta    = atan(displayTanb());
  double mtau    = tree.mtau;
  double mb      = tree.mb;
  double thetaWDRbar = asin(calcSinthdrbar());
  double cw2DRbar    = sqr(cos(thetaWDRbar));
  double tw2DRbar    = sqr(tan(thetaWDRbar));
  double thetat  = tree.thetat ;
  double thetab  = tree.thetab;
  double thetatau = tree.thetatau;
  double msbot1  = tree.md(1, 3);
  double msbot2  = tree.md(2, 3);
  double mstau1  = tree.me(1, 3);
  double mstau2  = tree.me(2, 3);
  double mstop1  = tree.mu(1, 3);
  double mstop2  = tree.mu(2, 3);
  double st      = sin(thetat) ;
  double mz      = displayMzRun();
  double mw	 = displayMwRun(), mw2 = sqr(mw); 
  double sb      = sin(thetab) ;
  double stau    = sin(thetatau);
  double ct      = cos(thetat) ;
  double cb      = cos(thetab) ;
  double ctau    = cos(thetatau);
  double g       = displayGaugeCoupling(2), gsq = sqr(g);
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta);
  double sinb =  sin(beta);
  double hb = tree.hb, htau = tree.htau;
  /// LCT: Extra NMSSM parameters
  double lam = displayLambda(), lsq = sqr(lam);
 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;	
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

   /// Fermions: 3rd family only for now
  double fermions = 3.0 * sqr(hb) *
     (sqr(p) * b0(p, mb, mb, q) - 2.0 * a0(mb, q));
  
  fermions = fermions + sqr(htau) *
     (sqr(p) * b0(p, mtau, mtau, q) - 2.0 * a0(mtau, q));
	
  /// Gauge contributions
  //
  /// LCT: Charged Higgs in basis G+ H+ G- H-
  double higgs = 0.0;
  double GaugeHiggs = 0.0;
  for (int i=1; i <=2; i++) {
     GaugeHiggs += gsq * 0.5 * sqr(C(i, 1)) * ffn(p, higgsc(i), mw, q);
  }
  
  /// LCT: CP-even states in basis H1 H2 H3
  for (int i = 1; i <= 3; i++) {
     GaugeHiggs = GaugeHiggs + gsq * 0.25 / cw2DRbar * 
        sqr(S(i, 1)) * ffn(p, higgsm(i), mz, q);
  } 
  higgs = higgs + GaugeHiggs;
  /// LCT: Gauge bosons
  double GaugeBosons =  gsq * ( 0.5 * cosb2 * mw2 * b0(p, mw, mw, q) 
                         + 2.0 * a0(mw, q) + 1.0 / cw2DRbar * a0(mz, q));
   higgs += GaugeBosons;
  
   /// Upsfermions
  //
  /// Quadrilinears gens 1-2
  /// LH 
  double upsfermions = 0.0;
  upsfermions = upsfermions +	3.0 * (gsq / (cw2DRbar) * 0.5 * guL) 
     * (a0(tree.mu(1, 1), q) + a0(tree.mu(1, 2), q));
  
  /// RH
  upsfermions = upsfermions +	3.0 * (gsq / (cw2DRbar) * 0.5 * guR) 
     * (a0(tree.mu(2, 1), q) + a0(tree.mu(2, 2), q));
  
   /// Downsfermions
  //
  /// Quadrilinears gens 1-2
  /// LH
  double dsfermions = 0.0;
  dsfermions = dsfermions + 
     3.0 * (sqr(g) / (cw2DRbar) * 0.5 * gdL) *
     (a0(tree.md(1, 1), q) + a0(tree.md(1, 2), q));
  
	/// RH
  dsfermions = dsfermions +
	3.0 * (sqr(g) / (cw2DRbar) * 0.5 * gdR) *
	(a0(tree.md(2, 1), q) + a0(tree.md(2, 2), q));
  
  	/// Quadrilinears gens 1-2
	/// LH
	double sleptons = 0.;
  sleptons = (sqr(g) / (cw2DRbar) * 0.5 * geL) *
	(a0(tree.me(1, 1), q) + a0(tree.me(1, 2), q));
	
	/// RH
  sleptons = sleptons +
     (sqr(g) / (cw2DRbar) * 0.5 * geR) *
     (a0(tree.me(2, 1), q) + a0(tree.me(2, 2), q));

  /// Sneutrinos. No trilinear terms
  double sneutrinos = sqr(g) / cw2DRbar * 0.5 * gnuL 
     * (a0(tree.msnu(1), q) + a0(tree.msnu(2), q) 
	+ a0(tree.msnu(3), q));
  
  	/// 3rd Generation
	/// Quadrilinears
  double stops = 0.0;
  stops = stops + 3.0 * (sqr(g) / (cw2DRbar) * 0.5 * guL) 
    * (sqr(ct) * a0(mstop1, q) + sqr(st) * a0(mstop2, q));
	
  stops = stops + 3.0 * sqr(g) / (cw2DRbar) * 0.5 * guR * 
    (sqr(st) * a0(mstop1, q) + sqr(ct) * a0(mstop2, q));
	 
  double sbots = 0.0;
  sbots = sbots + 3.0 * (sqr(hb) + sqr(g) / (cw2DRbar) * 0.5 * gdL) 
    * (sqr(cb) * a0(msbot1, q) + sqr(sb) * a0(msbot2, q));
	
  sbots = sbots + 3.0 *  (sqr(hb) + sqr(g) / (cw2DRbar) * 0.5 * gdR) 
    * (sqr(sb) * a0(msbot1, q) + sqr(cb) * a0(msbot2, q));
	 
  double staus = 0.0;
  staus = staus + 1.0 * (sqr(htau) + sqr(g) / (cw2DRbar) * 0.5 * geL) 
    * (sqr(ctau) * a0(mstau1, q) + sqr(stau) * a0(mstau2, q));
	
  staus = staus + 1.0 * (sqr(htau) +  sqr(g) / (cw2DRbar) * 0.5 * geR) 
    * (sqr(stau) * a0(mstau1, q) + sqr(ctau) * a0(mstau2, q));

  /// Trilinears
  //
  /// stop, sbottom, stau couplings to p1 Higgs state
  DoubleMatrix lp1tt(2, 2), lp1bb(2, 2), lp1tautau(2, 2);
  P1SfSfCouplings(lp1tt, lp1bb, lp1tautau);
  /// Mix 3rd family up
  lp1tt = rot2d(thetat) * lp1tt * rot2d(-thetat);
  lp1bb = rot2d(thetab) * lp1bb * rot2d(-thetab);
  lp1tautau = rot2d(thetatau) * lp1tautau * rot2d(-thetatau);
  
  int i, j; for (i=1; i <= 2; i++) {
    for (j=1; j <= 2; j++) {
      stops = stops + 3.0 * sqr(lp1tt(i, j)) * 
	b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
      sbots = sbots + 3.0 * sqr(lp1bb(i, j)) * 
	b0(p, tree.md(i, 3), tree.md(j, 3), q);
      staus = staus +  sqr(lp1tautau(i, j)) * 
	b0(p, tree.me(i, 3), tree.me(j, 3), q);
    }
  }	

  /// Higgs
  //	
  /// Quadrilinear CP-even Higgs couplings 
  DoubleMatrix ssp1p1(3, 3);
  ssp1p1(1, 1) = gsq / (16.0 * cw2DRbar);
  ssp1p1(2, 2) = (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 8.0;
  ssp1p1(3, 3) = 0.25 * lsq;	
  ssp1p1.symmetrise();
  
  /// Quadrilinear CP-odd Higgs couplings
  DoubleMatrix ppp1p1(3, 3);
  ppp1p1(1, 1) = gsq / (32.0 * cw2DRbar);
  ppp1p1(2, 2) = (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 48.0;
  ppp1p1(3, 3) = lsq / 24.0;
  ppp1p1.symmetrise();
	
  /// LCT: Rotate to mass bases p1 p1 Hi Hi and p1 p1 Ai Ai
  DoubleVector hhp1p1(3), aap1p1(3);
  for (int i = 1; i <= 3; i++) {
     for (int a = 1; a <= 3; a++) {
        for (int b = 1; b <= 3; b++) {
           hhp1p1(i) = hhp1p1(i) + S(i, a) * S(i, b) * ssp1p1(a, b);
           aap1p1(i) = aap1p1(i) + 6.0 * P(i, a) * P(i, b) * ppp1p1(a, b);
        }
     }
  }
  

  /// Trilinear Higgs couplings to p1 in basis (s, p)
  DoubleMatrix spp1(3, 3), hphpp1(2, 2);
  getP1HiggsTriCoup(spp1, hphpp1,cw2DRbar);
  
  /// LCT: Rotate to mass basis p1 Ai Hj
  DoubleMatrix ahp1(3, 3);
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <= 3; j++) {
        for (int a = 1; a <= 3; a++) {
           for (int b = 1; b <= 3; b++) {
              ahp1(i, j) = ahp1(i, j) + 2.0 * S(j, a) * P(i, b) * spp1(a, b);
           }
        }
     }
  }
  double cpeven = 0.0, pseudo = 0.0;
  for (i=1; i <= 3; i++) {
    for (j=1; j <= 3; j++) {
      cpeven = cpeven + sqr(ahp1(i, j)) * b0(p, higgsa(i), higgsm(j), q);
    }
    cpeven = cpeven + 2.0 * hhp1p1(i) * a0(higgsm(i), q);
    pseudo = pseudo + 2.0 * aap1p1(i) * a0(higgsa(i), q);
  }
  
  higgs += cpeven;
  higgs += pseudo;
  /// Quadrilinear charged Higgs.  Basis (G+ G- H+ H-)
  DoubleVector hphpp1p1(2);
  hphpp1p1(1) = gsq * (1.0 + tw2DRbar * cos2b) / 8.0;
  hphpp1p1(2) = gsq * (1.0 - tw2DRbar * cos2b) / 8.0;
   
  double ChargedHiggs = 0.0;
  for (i=1; i <= 2; i++) {
    for (j=1; j <= 2; j++) {
      ChargedHiggs =  ChargedHiggs + sqr(hphpp1(i, j)) * b0(p, higgsc(i), higgsc(j), q);
    }
    ChargedHiggs =  ChargedHiggs + 2.0 * hphpp1p1(i) * a0(higgsc(i), q);
  }
  
  higgs +=  ChargedHiggs;
  
  /// Neutralino contribution
  double neutralinos = 0.0;
	
 
  ComplexMatrix aChi(5, 5), bChi(5, 5);
  DoubleVector mneut(tree.mnBpmz);
  getP1NeutralinoCoup(aChi, bChi);
 	
  DoubleMatrix fChiChip1p1(5, 5), gChiChip1p1(5, 5);
  for(i=1; i <= 5; i++)
    for (j=1; j <= 5; j++) {
      fChiChip1p1(i, j) = sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod());
      gChiChip1p1(i, j) = (bChi(i, j).conj() * aChi(i, j) + 
													 aChi(i, j).conj() * bChi(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
			(fChiChip1p1(i, j) * gfn(p, mneut(i), mneut(j), q) + 2.0 *
			 gChiChip1p1(i, j) * mneut(i) * mneut(j) * 
			 b0(p, mneut(i), mneut(j), q));
    }
	
  /// Chargino contribution
  double chargino = 0.0;
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz); 
  DoubleVector mch(tree.mchBpmz); 
  DoubleMatrix aPsic(2, 2);
  aPsic(1, 2) = -g / root2; 
  ComplexMatrix aChic(2, 2), bChic(2, 2);
  aChic = v.complexConjugate() * aPsic * u.hermitianConjugate();
  bChic = u * aPsic.transpose() * v.transpose();
  for(i=1; i <= 2; i++)
    for (j=1; j <= 2; j++) {
      fChiChip1p1(i, j) = sqr(aChic(i, j).mod()) + sqr(bChic(i, j).mod());
      gChiChip1p1(i, j) = (bChic(i, j).conj() * aChic(i, j) + 
			   aChic(i, j).conj() * bChic(i, j)).real();
      chargino = chargino + 
	(fChiChip1p1(i, j) * gfn(p, mch(i), mch(j), q) + 2.0 *
	 gChiChip1p1(i, j) * mch(i) * mch(j) * 
	 b0(p, mch(i), mch(j), q));
    }	
  
 double sfermions = sneutrinos + sleptons + staus + upsfermions 
     + dsfermions + stops + sbots; 
 
  return 
     (fermions + sfermions + higgs + neutralinos + chargino) / (16.0 * sqr(PI));		
}



double NmssmSoftsusy::pip1p2(double p, double q) const {
  const drBarPars& tree = displayDrBarPars();
	
  double    beta    = atan(displayTanb());
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    mz      = displayMzRun();
  double    mw	    = displayMwRun();
  double    g       = displayGaugeCoupling(2), gsq = sqr(g);
  double cosb = cos(beta), sinb = sin(beta), sin2b = sin(2.0 * beta);
  /// LCT: Extra NMSSM parameters
  double lam = displayLambda(), lsq = sqr(lam);
  double kap = displayKappa();

  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;	
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// No fermion contribution
  
  /// Gauge contributions
  //
  /// LCT: Charged Higgs in basis G+ H+ G- H-
  double higgs = 0.0;
  double GaugeHiggs = 0.0;
  double GaugeBosons = 0.0;
  for (int i=1; i <= 2; i++) {
     GaugeHiggs += gsq * 0.5 * 
        C(i, 1) * C(i, 2) * ffn(p, higgsc(i), displayMwRun(), q);
  }
  
  /// LCT: CP-even states in basis H1 H2 H3
  for (int i = 1; i <= 3; i++) {
     GaugeHiggs = GaugeHiggs - gsq * 0.25 / cw2DRbar * 
        S(i, 1) * S(i, 2) * ffn(p, higgsm(i), mz, q);
  } 
  
  /// LCT: Gauge bosons
  GaugeBosons = GaugeBosons - gsq * 0.5 * cosb * sinb * sqr(mw) 
    * b0(p, mw, mw, q);

  /// Trilinears
  //
  /// stop, sbottom, stau couplings to p1 Higgs state
  DoubleMatrix lp1tt(2, 2), lp1bb(2, 2), lp1tautau(2, 2);
  P1SfSfCouplings(lp1tt, lp1bb, lp1tautau);
   /// stop, sbottom, stau couplings to p2 Higgs state
  DoubleMatrix lp2tt(2, 2), lp2bb(2, 2), lp2tautau(2, 2);
  P2SfSfCouplings(lp2tt, lp2bb, lp2tautau);
  /// Mix 3rd family up
  lp1tt = rot2d(thetat) * lp1tt * rot2d(-thetat);
  lp1bb = rot2d(thetab) * lp1bb * rot2d(-thetab);
  lp1tautau = rot2d(thetatau) * lp1tautau * rot2d(-thetatau);
  
  lp2tt = rot2d(thetat) * lp2tt * rot2d(-thetat);
  lp2bb = rot2d(thetab) * lp2bb * rot2d(-thetab);
  lp2tautau = rot2d(thetatau) * lp2tautau * rot2d(-thetatau);
  

  double stops = 0.0, sbots = 0.0, staus = 0.0;
  int i, j; for (i=1; i <= 2; i++) {
     for (j=1; j <= 2; j++) {
        stops = stops - 3.0 * lp1tt(i, j) * lp2tt(j, i) 
           * b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
        sbots = sbots - 3.0 * lp1bb(i, j) * lp2bb(j, i) 
           * b0(p, tree.md(i, 3), tree.md(j, 3), q);
        staus = staus -  lp1tautau(i, j) * lp2tautau(j, i) 
           * b0(p, tree.me(i, 3), tree.me(j, 3), q);
     }}	

   /// Higgs
  /// Quadrilinear CP-even Higgs couplings 
  DoubleMatrix ssp1p2(3, 3);
  ssp1p2(3, 3) = 0.25 * lam * kap;
  ssp1p2.symmetrise();
  
  /// Quadrilinear CP-odd Higgs couplings
  DoubleMatrix ppp1p2(3, 3);
  ppp1p2(1, 2) = (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 48.0;
  ppp1p2(3, 3) = - lam * kap / 24.0;
  ppp1p2.symmetrise();
  
  /// LCT: Rotate to mass bases p1 p2 Hi Hi and p1 p2 Ai Ai
  DoubleVector hhp1p2(3), aap1p2(3);
  for (int i = 1; i <= 3; i++) {
     for (int a = 1; a <= 3; a++) {
        for (int b = 1; b <= 3; b++) {
           hhp1p2(i) = hhp1p2(i) + S(i, a) * S(i, b) * ssp1p2(a, b);
           aap1p2(i) = aap1p2(i) + 6.0 * P(i, a) * P(i, b) * ppp1p2(a, b);
        }
     }
  }

  DoubleMatrix spp1(3, 3), spp2(3, 3), hphpp2(2, 2),  hphpp1(2, 2);
  getP1HiggsTriCoup(spp1, hphpp1, cw2DRbar);
  getP2HiggsTriCoup(spp2, hphpp2, cw2DRbar);

   /// LCT: Rotate to mass basis p1 Ai Hj and p2 Ai Hj
  DoubleMatrix ahp1(3, 3), ahp2(3, 3);
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <= 3; j++) {
        for (int a = 1; a <= 3; a++) {
           for (int b = 1; b <= 3; b++) {
              ahp1(i, j) = ahp1(i, j) + 2.0 * S(j, a) * P(i, b) * spp1(a, b);
              ahp2(i, j) = ahp2(i, j) + 2.0 * S(j, a) * P(i, b) * spp2(a, b);
           }
        }
     }
  }
  double cpodd = 0.0;
  double cpeven = 0.0;
  for (i=1; i <= 3; i++) {
     for (j=1; j <= 3; j++) {
        cpeven = cpeven + ahp1(i, j) * ahp2(i, j) * b0(p, higgsa(i), higgsm(j), q);
     }
     cpeven = cpeven + 2.0 * hhp1p2(i) * a0(higgsm(i), q);
     cpodd = cpodd + 2.0 * aap1p2(i) * a0(higgsa(i), q);
  }
  
  /// Charged quadrilinear (G+ G- H+ H-)
  DoubleVector hphpp1p2(2);
  hphpp1p2(1) = -(2.0 * lsq - gsq) * sin2b / 8.0;
  hphpp1p2(2) = (2.0 * lsq - gsq) * sin2b / 8.0;

  double ChargedHiggs = 0.0;
  for (i=1; i <= 2; i++) {
     for (j=1; j <= 2; j++) {
        ChargedHiggs = ChargedHiggs - hphpp1(i, j) * hphpp2(j, i) * b0(p, higgsc(i), higgsc(j), q);
     }
     ChargedHiggs = ChargedHiggs  + 2.0 * hphpp1p2(i) * a0(higgsc(i), q);
  }

  higgs = cpodd + cpeven + ChargedHiggs + GaugeHiggs + GaugeBosons;

   /// Neutralino contribution
  double neutralinos = 0.0;
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix aChi1(5, 5), bChi1(5, 5),  aChi2(5, 5), bChi2(5, 5);
  getP1NeutralinoCoup(aChi1, bChi1);
  getP2NeutralinoCoup(aChi2, bChi2);
  DoubleMatrix fChiChip1p2(5, 5), gChiChip1p2(5, 5);
  for(i=1; i <= 5; i++)
    for (j=1; j <= 5; j++) {
      fChiChip1p2(i, j) = (aChi1(i, j).conj() * aChi2(i, j) + 
			   bChi1(i, j).conj() * bChi2(i, j)).real();
      gChiChip1p2(i, j) = (bChi1(i, j).conj() * aChi2(i, j) + 
			   aChi1(i, j).conj() * bChi2(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
	(fChiChip1p2(i, j) * gfn(p, mneut(i), mneut(j), q) + 2.0 *
	 gChiChip1p2(i, j) * mneut(i) * mneut(j) * 
	 b0(p, mneut(i), mneut(j), q));
    }
 /// Chargino contribution
  double chargino = 0.0;
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz); 
  DoubleVector mch(tree.mchBpmz); 
  DoubleMatrix aPsic1(2, 2), aPsic2(2, 2);
  aPsic1(1, 2) = -g / root2; 
  ComplexMatrix aChic1(2, 2), bChic1(2, 2);
  ComplexMatrix aChic2(2, 2), bChic2(2, 2);
  aChic1 = v.complexConjugate() * aPsic1 * u.hermitianConjugate();
  bChic1 = u * aPsic1.transpose() * v.transpose();
  aPsic2(2, 1) = -g / root2;
  aChic2 = v.complexConjugate() * aPsic2 * u.hermitianConjugate();
  bChic2 = u * aPsic2.transpose() * v.transpose();
	
  for(i=1; i <= 2; i++)
     for (j=1; j <= 2; j++) {
        fChiChip1p2(i, j) = (aChic1(i, j).conj() * aChic2(i, j) + 
                             bChic1(i, j).conj() * bChic2(i, j)).real();
        gChiChip1p2(i, j) = (bChic1(i, j).conj() * aChic2(i ,j) + 
                             aChic1(i, j).conj() * bChic2(i, j)).real();
      chargino = chargino + 
         (fChiChip1p2(i, j) * gfn(p, mch(i), mch(j), q) + 2.0 *
          gChiChip1p2(i, j) * mch(i) * mch(j) * 
          b0(p, mch(i), mch(j), q));
     }
  
  double sfermions = staus + stops + sbots; 
  
  return (sfermions + higgs + neutralinos + chargino) / (16.0 * sqr(PI));
}


double NmssmSoftsusy::pip2p2(double p, double q) const {
  const drBarPars& tree = displayDrBarPars();
	
  double beta    = atan(displayTanb());
  double mt      = tree.mt;
  double thetaWDRbar = asin(calcSinthdrbar());
  double cw2DRbar    = sqr(cos(thetaWDRbar));
  double tw2DRbar    = sqr(tan(thetaWDRbar));
  double thetat  = tree.thetat ;
  double thetab  = tree.thetab;
  double thetatau = tree.thetatau;
  double msbot1  = tree.md(1, 3);
  double msbot2  = tree.md(2, 3);
  double mstau1  = tree.me(1, 3);
  double mstau2  = tree.me(2, 3);
  double mstop1  = tree.mu(1, 3);
  double mstop2  = tree.mu(2, 3);
  double st      = sin(thetat) ;
  double mz      = displayMzRun();
  double mw	 = displayMwRun(), mw2 = sqr(mw); 
  double sb      = sin(thetab) ;
  double stau    = sin(thetatau);
  double ct      = cos(thetat) ;
  double cb      = cos(thetab) ;
  double ctau    = cos(thetatau);
  double g       = displayGaugeCoupling(2), gsq = sqr(g);
  double cosb = cos(beta), cos2b = cos(2.0 * beta),
     sinb = sin(beta), sinb2 = sqr(sinb);
  double ht = tree.ht;
 	
	/// LCT: Extra NMSSM parameters
  double lam = displayLambda(), lsq = sqr(lam);
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;	
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// Fermions: 3rd family only for now
  double fermions = 3.0 * sqr(ht) *	(sqr(p) * b0(p, mt, mt, q) - 2.0 * a0(mt, q));
  
  /// Gauge contributions
  /// Charged Higgs in basis (G+ H+ G- H-)
  double higgs = 0.0;	
  double GaugeHiggs = 0.0;
  for (int i=1; i <= 2; i++) {
     GaugeHiggs += gsq * 0.5 * sqr(C(i, 2)) * ffn(p, higgsc(i), mw, q);
  }
  
  /// LCT: CP-even states in basis H1 H2 H3
  for (int i = 1; i <= 3; i++) {
     GaugeHiggs = GaugeHiggs + 0.25 * gsq / cw2DRbar * 
        sqr(S(i, 2)) * ffn(p, higgsm(i), mz, q);
  } 
  
  /// LCT: Gauge bosons
  double GaugeBosons = gsq * ( 0.5 * sinb2 * mw2 * b0(p, mw, mw, q) 
			       + 2.0 * a0(mw, q) + a0(mz, q) / cw2DRbar);

  /// Upsfermions
  /// Quadrilinears gens 1-2
  /// LH 
  double upsfermions = 0.0;
  upsfermions = upsfermions + 3.0 * (-gsq / (cw2DRbar) * 0.5 * guL) 
     * (a0(tree.mu(1, 1), q) + a0(tree.mu(1, 2), q));
  
  /// RH
  upsfermions = upsfermions + 3.0 * (-gsq / (cw2DRbar) * 0.5 * guR) 
              * (a0(tree.mu(2, 1), q) + a0(tree.mu(2, 2), q));

  /// Downsfermions
  /// Quadrilinears gens 1-2
  /// LH
  double dsfermions = 0.0;
  dsfermions = dsfermions + 3.0 * (-gsq / (cw2DRbar) * 0.5 * gdL) 
    * (a0(tree.md(1, 1), q) + a0(tree.md(1, 2), q));
  
  /// RH
  dsfermions = dsfermions + 3.0 * (-sqr(g) / (cw2DRbar) * 0.5 * gdR) 
    * (a0(tree.md(2, 1), q) + a0(tree.md(2, 2), q));

  /// Slepton
  /// Quadrilinears gens 1-2
  /// LH
  double sleptons = 0.;
  sleptons = (-gsq / (cw2DRbar) * 0.5 * geL) 
    * (a0(tree.me(1, 1), q) + a0(tree.me(1, 2), q));
  
  /// RH
  sleptons = sleptons + (-gsq / (cw2DRbar) * 0.5 * geR) 
    * (a0(tree.me(2, 1), q) + a0(tree.me(2, 2), q));
  
  /// Sneutrinos. 
  double sneutrinos = - gsq / cw2DRbar * 0.5 * gnuL * (a0(tree.msnu(1), q) 
                              + a0(tree.msnu(2), q) + a0(tree.msnu(3), q));

   /// 3rd Generation
  /// Quadrilinears
  double stops = 0.0;
  stops = stops + 3.0 * (-gsq / (cw2DRbar) * guL * 0.5 + sqr(ht)) 
     * (sqr(ct) * a0(mstop1, q) + sqr(st) * a0(mstop2, q));
	
  stops = stops + 3.0 * (-gsq / (cw2DRbar) * guR * 0.5 + sqr(ht)) 
        * (sqr(st) * a0(mstop1, q) + sqr(ct) * a0(mstop2, q));
	
  double sbots = 0.0;
  sbots = sbots + 3.0 * (-gsq / (cw2DRbar) * gdL * 0.5) 
     * (sqr(cb) * a0(msbot1, q) + sqr(sb) * a0(msbot2, q));
  
  sbots = sbots + 3.0 * (-gsq / (cw2DRbar) * gdR * 0.5) 
     * (sqr(sb) * a0(msbot1, q) + sqr(cb) * a0(msbot2, q));
  
  double staus = 0.0;
  staus = staus + 1.0 * (-gsq / (cw2DRbar) * geL * 0.5)
     * (sqr(ctau) * a0(mstau1, q) + sqr(stau) * a0(mstau2, q));
  
  staus = staus + 1.0 * (-gsq / (cw2DRbar) * geR * 0.5)
	      * (sqr(stau) * a0(mstau1, q) + sqr(ctau) * a0(mstau2, q));

  /// stop, sbottom, stau couplings to p2 Higgs state
  DoubleMatrix lp2tt(2, 2), lp2bb(2, 2), lp2tautau(2, 2);
  P2SfSfCouplings(lp2tt, lp2bb, lp2tautau);
  /// Mix 3rd family up
  lp2tt = rot2d(thetat) * lp2tt * rot2d(-thetat);
  lp2bb = rot2d(thetab) * lp2bb * rot2d(-thetab);
  lp2tautau = rot2d(thetatau) * lp2tautau * rot2d(-thetatau);
  int i, j; for (i=1; i <= 2; i++) {
     for (j=1; j <= 2; j++) {
        stops = stops + 3.0 * sqr(lp2tt(i, j)) * 
           b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
        sbots = sbots + 3.0 * sqr(lp2bb(i, j)) *
           b0(p, tree.md(i, 3), tree.md(j, 3), q);
        staus = staus + sqr(lp2tautau(i, j)) * 
           b0(p, tree.me(i, 3), tree.me(j, 3), q);
     }
  }	

  /// Higgs
  /// Quadrilinear CP-even Higgs couplings 
  DoubleMatrix ssp2p2(3, 3);
  ssp2p2(1, 1) = (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 8.0;
  ssp2p2(2, 2) = gsq / (16.0 * cw2DRbar);
  ssp2p2(3, 3) = 0.25 * lsq;
  ssp2p2.symmetrise();
  
  /// Quadrilinear CP-odd Higgs couplings
  DoubleMatrix ppp2p2(3, 3);
  ppp2p2(1, 1) = (2.0 * lsq - 0.5 * gsq / cw2DRbar) / 48.0;
  ppp2p2(2, 2) = gsq / (32.0 * cw2DRbar);
  ppp2p2(3, 3) = lsq / 24.0;
  ppp2p2.symmetrise();
  
  /// LCT: Rotate to mass bases p2 p2 Hi Hi and p2 p2 Ai Ai
  DoubleVector hhp2p2(3), aap2p2(3);
  for (int i = 1; i <= 3; i++) {
     for (int a = 1; a <= 3; a++) {
       for (int b = 1; b <= 3; b++) {
	 hhp2p2(i) = hhp2p2(i) + S(i, a) * S(i, b) * ssp2p2(a, b);
	 aap2p2(i) = aap2p2(i) + 6.0 * P(i, a) * P(i, b) * ppp2p2(a, b);
       }
     }
  }
  
  DoubleMatrix spp2(3, 3), hphpp2(2, 2);
  getP2HiggsTriCoup(spp2, hphpp2, cw2DRbar);
  
  /// LCT: Rotate to mass basis p2 Ai Hj
  DoubleMatrix ahp2(3, 3);
  for (int i=1; i <= 3; i++) {
    for (int j=1; j <= 3; j++) {
      for (int a = 1; a <= 3; a++) {
	for (int b = 1; b <= 3; b++) {
	  ahp2(i, j) = ahp2(i, j) + 2.0 * S(j, a) * P(i, b) * spp2(a, b);
	}
      }
    }
  }
  
  double cpodd = 0.0, cpeven = 0.0;
  for (i=1; i <= 3; i++) {
    for (j=1; j <= 3; j++) {
      cpeven = cpeven + sqr(ahp2(i, j)) * b0(p, higgsa(i), higgsm(j), q);
    }
    cpeven = cpeven + 2.0 * hhp2p2(i) * a0(higgsm(i), q);
    cpodd = cpodd + 2.0 * aap2p2(i) * a0(higgsa(i), q);
  }
  
  /// Quadrilinear (G+ G- H+ H-)
  DoubleVector hphpp2p2(2);
  hphpp2p2(1) = gsq * (1.0 - tw2DRbar * cos2b) / 8.0;
  hphpp2p2(2) = gsq * (1.0 + tw2DRbar * cos2b) / 8.0;
  
  double ChargedHiggs = 0.0;
  for (i=1; i <= 2; i++) {
    for (j=1; j <= 2; j++) {
      ChargedHiggs = ChargedHiggs + sqr(hphpp2(i, j)) * b0(p, higgsc(i), higgsc(j), q);
     }
    ChargedHiggs = ChargedHiggs + 2.0 * hphpp2p2(i) * a0(higgsc(i), q);
  }
  
  higgs = GaugeBosons + GaugeHiggs + ChargedHiggs + cpeven + cpodd;

  /// Neutralino contribution
  double neutralinos = 0.0;
  ComplexMatrix aChi(5, 5), bChi(5, 5);
  getP2NeutralinoCoup(aChi, bChi);
  DoubleVector mneut(tree.mnBpmz);
   DoubleMatrix fChiChip2p2(5, 5), gChiChip2p2(5, 5);
  for(i=1; i <= 5; i++)
    for (j=1; j <= 5; j++) {
      fChiChip2p2(i, j) = sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod());
      gChiChip2p2(i, j) = (bChi(i, j).conj() * aChi(i, j) + 
													 aChi(i, j).conj() * bChi(i, j)).real();
      neutralinos = neutralinos + 0.5 * (fChiChip2p2(i, j) 
                  * gfn(p, mneut(i), mneut(j), q) + 2.0 * gChiChip2p2(i, j) 
                  * mneut(i) * mneut(j) * b0(p, mneut(i), mneut(j), q));
    }


  /// Chargino contribution
  double chargino = 0.0;
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz); 
   DoubleVector mch(tree.mchBpmz);
  DoubleMatrix aPsic(2, 2);
  aPsic(2, 1) = -g / root2; 
  ComplexMatrix aChic(2, 2), bChic(2, 2);
  aChic = v.complexConjugate() * aPsic * u.hermitianConjugate();
  bChic = u * aPsic.transpose() * v.transpose();
  for(i=1; i <= 2; i++)
    for (j=1; j <= 2; j++) {
      fChiChip2p2(i, j) = sqr(aChic(i, j).mod()) + sqr(bChic(i, j).mod());
      gChiChip2p2(i, j) = (bChic(i, j).conj() * aChic(i, j) + 
													 aChic(i, j).conj() * bChic(i, j)).real();
      chargino = chargino + (fChiChip2p2(i, j) * gfn(p, mch(i), mch(j), q) 
               + 2.0 * gChiChip2p2(i, j) * mch(i) * mch(j) 
               * b0(p, mch(i), mch(j), q));
    }
	
  double sfermions = sneutrinos + sleptons + staus + upsfermions 
     + dsfermions + stops + sbots; 

  return (fermions + sfermions + higgs + neutralinos 
	  + chargino) / (16.0 * sqr(PI));

}


double NmssmSoftsusy::pip1p3(double p, double q) const {
  const drBarPars& tree = displayDrBarPars();
	
  double    beta    = atan(displayTanb());
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    mstop1  = tree.mu(1, 3);
  double    mstop2  = tree.mu(2, 3);
  double    st      = sin(thetat) ;
  double    ct      = cos(thetat) ;
  double    g       = displayGaugeCoupling(2);
  double cosb = cos(beta), sinb = sin(beta);
  double ht = tree.ht;
 	
  /// LCT: Extra NMSSM parameters
  double lam = displayLambda(), lsq = sqr(lam);
  double kap = displayKappa();

  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;	
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

   /// 3rd Generation sfermions
  /// Quadrilinears
  double stops = 0.0;
  stops = stops + 3.0 * ht * lam * ct * st * (a0(mstop1, q) - a0(mstop2, q));

  /// stop, sbottom, stau couplings to p1 Higgs state
  DoubleMatrix lp1tt(2, 2), lp1bb(2, 2), lp1tautau(2, 2);
  P1SfSfCouplings(lp1tt, lp1bb, lp1tautau);
  /// stop, sbottom, stau couplings to p3 Higgs state
  DoubleMatrix lp3tt(2, 2), lp3bb(2, 2), lp3tautau(2, 2);
  P3SfSfCouplings(lp3tt, lp3bb, lp3tautau);
   /// Mix 3rd family up
  lp1tt = rot2d(thetat) * lp1tt * rot2d(-thetat);
  lp1bb = rot2d(thetab) * lp1bb * rot2d(-thetab);
  lp1tautau = rot2d(thetatau) * lp1tautau * rot2d(-thetatau);
  
  lp3tt = rot2d(thetat) * lp3tt * rot2d(-thetat);
  lp3bb = rot2d(thetab) * lp3bb * rot2d(-thetab);
  lp3tautau = rot2d(thetatau) * lp3tautau * rot2d(-thetatau);
  
  double sbots = 0.0, staus = 0.0;
  int i, j; for (i=1; i <= 2; i++) {
    for (j=1; j<=2; j++) {
      stops = stops - 3.0 * lp1tt(i, j) * lp3tt(j, i) 
	* b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
      sbots = sbots - 3.0 * lp1bb(i, j) * lp3bb(j, i) 
	* b0(p, tree.md(i, 3), tree.md(j, 3), q);
      staus = staus -  lp1tautau(i, j) * lp3tautau(j, i) 
	* b0(p, tree.me(i, 3), tree.me(j, 3), q);
    }
  }	
  
  /// Higgs
  /// Quadrilinear CP-even Higgs couplings 
  DoubleMatrix ssp1p3(3, 3);
  ssp1p3(2, 3) = -0.25 * lam * kap;
  ssp1p3.symmetrise();
  
  /// Quadrilinear CP-odd Higgs couplings
  DoubleMatrix ppp1p3(3, 3);
  ppp1p3(1, 3) = lsq / 24.0;
  ppp1p3(2, 3) = -lam * kap / 24.0;
  ppp1p3.symmetrise();
  /// LCT: Rotate to mass bases p1 p3 Hi Hi and p1 p3 Ai Ai
  DoubleVector hhp1p3(3), aap1p3(3);
  for (int i = 1; i <= 3; i++) {
     for (int a = 1; a <= 3; a++) {
        for (int b = 1; b <= 3; b++) {
           hhp1p3(i) = hhp1p3(i) + S(i, a) * S(i, b) * ssp1p3(a, b);
           aap1p3(i) = aap1p3(i) + 6.0 * P(i, a) * P(i, b) * ppp1p3(a, b);
        }
     }
  }

  DoubleMatrix spp1(3, 3), spp3(3, 3), hphpp1(2, 2),  hphpp3(2, 2);
  getP1HiggsTriCoup(spp1, hphpp1, cw2DRbar);
  getP3HiggsTriCoup(spp3, hphpp3);
  
  /// LCT: Rotate to mass basis p1 Ai Hj and p3 Ai Hj
  DoubleMatrix ahp1(3, 3), ahp3(3, 3);
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <= 3; j++) {
       for (int a = 1; a <= 3; a++) {
	 for (int b = 1; b <= 3; b++) {
	   ahp1(i, j) = ahp1(i, j) + 2.0 * S(j, a) * P(i, b) * spp1(a, b);
	   ahp3(i, j) = ahp3(i, j) + 2.0 * S(j, a) * P(i, b) * spp3(a, b);
           }
       }
     }
  }

  double neutraltri = 0.0, cpeven = 0.0, cpodd = 0.0;
  double higgs = 0.0;
  for (i=1; i <= 3; i++) {
    for (j=1; j <= 3; j++) {
      neutraltri   =  neutraltri + ahp1(i, j) * ahp3(i, j) 
	* b0(p, higgsa(i), higgsm(j), q);
    }
    cpeven = cpeven + 2.0 * hhp1p3(i) * a0(higgsm(i), q);
    cpodd = cpodd + 2.0 * aap1p3(i) * a0(higgsa(i), q);
  }
  
  double ChargedHiggs = 0.0;
  for (i=1; i <= 2; i++) {
    for (j=1; j <= 2; j++) {
       ChargedHiggs = ChargedHiggs - hphpp1(i, j) * hphpp3(j, i) 
	 * b0(p, higgsc(i), higgsc(j), q);
     }
  }
   higgs =  neutraltri + cpodd + cpeven + ChargedHiggs;

   /// Neutralino contribution
  double neutralinos = 0.0;
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix aChi1(5, 5), bChi1(5, 5),  aChi3(5, 5), bChi3(5, 5);
  getP1NeutralinoCoup(aChi1, bChi1);
  getP3NeutralinoCoup(aChi3, bChi3);
  DoubleMatrix fChiChip1p3(5, 5), gChiChip1p3(5, 5);
  for(i=1; i<=5; i++)
    for (j=1; j<=5; j++) {
      fChiChip1p3(i, j) = (aChi1(i, j).conj() * aChi3(i, j) + 
			   bChi1(i, j).conj() * bChi3(i, j)).real();
      gChiChip1p3(i, j) = (bChi1(i, j).conj() * aChi3(i, j) + 
			   aChi1(i, j).conj() * bChi3(i, j)).real();
      neutralinos = neutralinos 
	+ 0.5 * (fChiChip1p3(i, j) * gfn(p, mneut(i), mneut(j), q) 
		 + 2.0 * gChiChip1p3(i, j) 
		 * mneut(i) * mneut(j) * b0(p, mneut(i), mneut(j), q));
     }
   /// Chargino contribution
  double chargino = 0.0;
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz); 
  DoubleVector mch(tree.mchBpmz); 
  DoubleMatrix aPsic1(2, 2), aPsic3(2, 2);
  aPsic1(1, 2) = -g / root2; 
  ComplexMatrix aChic1(2, 2), bChic1(2, 2);
  ComplexMatrix aChic3(2, 2), bChic3(2, 2);
  aChic1 = v.complexConjugate() * aPsic1 * u.hermitianConjugate();
  bChic1 = u * aPsic1.transpose() * v.transpose();
  aPsic3(2, 2) = lam / root2;
  aChic3 = v.complexConjugate() * aPsic3 * u.hermitianConjugate();
  bChic3 = u * aPsic3.transpose() * v.transpose();
  
  for(i=1; i <= 2; i++)
     for (j=1; j <= 2; j++) {
        fChiChip1p3(i, j) = (aChic1(i, j).conj() * aChic3(i, j) + 
                             bChic1(i, j).conj() * bChic3(i, j)).real();
        gChiChip1p3(i, j) = (bChic1(i, j).conj() * aChic3(i ,j) + 
                             aChic1(i, j).conj() * bChic3(i, j)).real();
        chargino = chargino + (fChiChip1p3(i, j) * gfn(p, mch(i), mch(j), q) 
                               + 2.0 * gChiChip1p3(i, j) * mch(i) * mch(j) 
                               * b0(p, mch(i), mch(j), q));
     }
  
  double sfermions = staus + stops + sbots; 

  return (sfermions + higgs + neutralinos + chargino) / (16.0 * sqr(PI));
  
}


double NmssmSoftsusy::pip2p3(double p, double q) const {
  const drBarPars& tree = displayDrBarPars();
  double    beta    = atan(displayTanb());
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    msbot1  = tree.md(1, 3);
  double    msbot2  = tree.md(2, 3);
  double    mstau1  = tree.me(1, 3);
  double    mstau2  = tree.me(2, 3);
  double    sb      = sin(thetab) ;
  double    stau    = sin(thetatau);
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);
  double    g       = displayGaugeCoupling(2);
  double cosb = cos(beta), sinb = sin(beta);
  double hb = tree.hb, htau = tree.htau;
  /// LCT: Extra NMSSM parameters
  double lam = displayLambda();
  double lsq = sqr(lam);
  double kap = displayKappa();

  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;	
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// 3rd Generation sfermions
  /// Quadrilinears
  double sbots = 0.0;
  sbots = sbots + 3.0 * hb * lam * cb * sb  * (a0(msbot1, q) - a0(msbot2, q));
  
  double staus = 0.0;
  staus = staus + htau * lam * ctau * stau * (a0(mstau1, q) - a0(mstau2, q));
  /// Trilinears
  /// stop, sbottom, stau couplings to p1 Higgs state
  DoubleMatrix lp2tt(2, 2), lp2bb(2, 2), lp2tautau(2, 2);
  P2SfSfCouplings(lp2tt, lp2bb, lp2tautau);
  /// stop, sbottom, stau couplings to p3 Higgs state
  DoubleMatrix lp3tt(2, 2), lp3bb(2, 2), lp3tautau(2, 2);
  P3SfSfCouplings(lp3tt, lp3bb, lp3tautau);
  /// Mix 3rd family up
  lp2tt = rot2d(thetat) * lp2tt * rot2d(-thetat);
  lp2bb = rot2d(thetab) * lp2bb * rot2d(-thetab);
  lp2tautau = rot2d(thetatau) * lp2tautau * rot2d(-thetatau);
  
  lp3tt = rot2d(thetat) * lp3tt * rot2d(-thetat);
  lp3bb = rot2d(thetab) * lp3bb * rot2d(-thetab);
  lp3tautau = rot2d(thetatau) * lp3tautau * rot2d(-thetatau);
  
  double stops = 0.0;
  int i, j; for (i=1; i<=2; i++) {
    for (j=1; j<=2; j++) {
       stops = stops - 3.0 * lp2tt(i, j) * lp3tt(j, i) * 
          b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
       sbots = sbots - 3.0 * lp2bb(i, j) * lp3bb(j, i) *
          b0(p, tree.md(i, 3), tree.md(j, 3), q);
       staus = staus -  lp2tautau(i, j) * lp3tautau(j, i) * 
          b0(p, tree.me(i, 3), tree.me(j, 3), q);
    }
  }	
   /// Higgs	
  /// Quadrilinear CP-even Higgs couplings 
  DoubleMatrix ssp2p3(3, 3);
  ssp2p3(1, 3) = - 0.25 * lam * kap;
  ssp2p3.symmetrise();
  
  /// Quadrilinear CP-odd Higgs couplings
  DoubleMatrix ppp2p3(3, 3);
  ppp2p3(1, 3) = - lam * kap / 24.0;
  ppp2p3(2, 3) = lsq / 24.0;
  ppp2p3.symmetrise();
  /// LCT: Rotate to mass bases p2 p3 Hi Hi and p2 p3 Ai Ai
  DoubleVector hhp2p3(3), aap2p3(3);
  for (int i = 1; i <=3; i++) {
     for (int a = 1; a <= 3; a++) {
       for (int b = 1; b <=3; b++) {
	 hhp2p3(i) = hhp2p3(i) + S(i, a) * S(i, b) * ssp2p3(a, b);
	 aap2p3(i) = aap2p3(i) + 6.0 * P(i, a) * P(i, b) * ppp2p3(a, b);
         
        }
     }
  }
  
  DoubleMatrix spp2(3, 3), spp3(3, 3), hphpp2(2, 2),  hphpp3(2, 2);
  getP2HiggsTriCoup(spp2, hphpp2, cw2DRbar);
  getP3HiggsTriCoup(spp3, hphpp3);
 
  /// LCT: Rotate to mass basis p2 Ai Hj and p3 Ai Hj
  DoubleMatrix ahp2(3, 3), ahp3(3, 3);
  for (int i=1; i <= 3; i++) {
    for (int j=1; j <=3; j++) {
      for (int a = 1; a <= 3; a++) {
	for (int b = 1; b <= 3; b++) {
	  ahp2(i, j) =  ahp2(i, j) + 2.0 * S(j, a) * P(i, b) * spp2(a, b);
	  ahp3(i, j) = ahp3(i, j) + 2.0 * S(j, a) * P(i, b) * spp3(a, b);
          

	}
      }
    }
  }
 
  double cpeven = 0.0;   
  double cpodd = 0.0;
  double higgs = 0.0;
  for (i=1; i<=3; i++) {
     for (j=1; j<=3; j++) {
       cpeven = cpeven + ahp2(i, j) * ahp3(i, j) * b0(p, higgsa(i), higgsm(j), q);
     }
     cpeven = cpeven + 2.0 * hhp2p3(i) * a0(higgsm(i), q);
     cpodd = cpodd + 2.0 * aap2p3(i) * a0(higgsa(i), q);
  }

  double ChargedHiggs =0.0;
  for (i=1; i<=2; i++) {
     for (j=1; j<=2; j++) {
        ChargedHiggs = ChargedHiggs - hphpp2(i, j) * hphpp3(j, i) * b0(p, higgsc(i), higgsc(j), q);
     }
  }
  higgs = cpeven + cpodd + ChargedHiggs;

   /// Neutralino contribution
  double neutralinos = 0.0;
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix aChi2(5, 5), bChi2(5, 5),  aChi3(5, 5), bChi3(5, 5);
  getP2NeutralinoCoup(aChi2, bChi2);
  getP3NeutralinoCoup(aChi3, bChi3);
   DoubleMatrix fChiChip2p3(5, 5), gChiChip2p3(5, 5);
  for(i=1; i<=5; i++)
    for (j=1; j<=5; j++) {
      fChiChip2p3(i, j) = (aChi2(i, j).conj() * aChi3(i, j) + 
													 bChi2(i, j).conj() * bChi3(i, j)).real();
      gChiChip2p3(i, j) = (bChi2(i, j).conj() * aChi3(i, j) + 
													 aChi2(i, j).conj() * bChi3(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
			(fChiChip2p3(i, j) * gfn(p, mneut(i), mneut(j), q) + 2.0 *
			 gChiChip2p3(i, j) * mneut(i) * mneut(j) * 
			 b0(p, mneut(i), mneut(j), q));
    }
	

  /// Chargino contribution
  double chargino = 0.0;
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz); 
  DoubleVector mch(tree.mchBpmz); 
  DoubleMatrix aPsic2(2, 2), aPsic3(2, 2);
  aPsic2(2, 1) = -g / root2; 
  ComplexMatrix aChic2(2, 2), bChic2(2, 2);
  ComplexMatrix aChic3(2, 2), bChic3(2, 2);
  aChic2 = v.complexConjugate() * aPsic2 * u.hermitianConjugate();
  bChic2 = u * aPsic2.transpose() * v.transpose();
  aPsic3(2, 2) = lam / root2;
  aChic3 = v.complexConjugate() * aPsic3 * u.hermitianConjugate();
  bChic3 = u * aPsic3.transpose() * v.transpose();
	
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fChiChip2p3(i, j) = (aChic2(i, j).conj() * aChic3(i, j) + 
			   bChic2(i, j).conj() * bChic3(i, j)).real();
      gChiChip2p3(i, j) = (bChic2(i, j).conj() * aChic3(i ,j) + 
			   aChic2(i, j).conj() * bChic3(i, j)).real();
      chargino = chargino + 
	(fChiChip2p3(i, j) * gfn(p, mch(i), mch(j), q) + 2.0 *
	 gChiChip2p3(i, j) * mch(i) * mch(j) * 
	 b0(p, mch(i), mch(j), q));
    }
  
  double sfermions = staus + stops + sbots;  

  return (sfermions + higgs + neutralinos 
	   + chargino) / (16.0 * sqr(PI));
   }

double NmssmSoftsusy::pip3p3(double p, double q) const {
  const drBarPars& tree = displayDrBarPars();
	
  double    beta    = atan(displayTanb());
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double cosb = cos(beta), sinb = sin(beta), sin2b = sin(2.0 * beta);
  /// LCT: Extra NMSSM parameters
  double lam = displayLambda(), lsq =sqr(lam);
  double kap = displayKappa(), ksq = sqr(kap);

  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;	
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);
  
   /// stop, sbottom, stau couplings to p3 Higgs state
  DoubleMatrix lp3tt(2, 2), lp3bb(2, 2), lp3tautau(2, 2);
  P3SfSfCouplings(lp3tt, lp3bb, lp3tautau);
  
  /// Mix 3rd family up
  lp3tt = rot2d(thetat) * lp3tt * rot2d(-thetat);
  lp3bb = rot2d(thetab) * lp3bb * rot2d(-thetab);
  lp3tautau = rot2d(thetatau) * lp3tautau * rot2d(-thetatau);
  
  double stops = 0.0, sbots = 0.0, staus = 0.0;
  int i, j; for (i=1; i<=2; i++) {
     for (j=1; j<=2; j++) {
        stops = stops + 3.0 * sqr(lp3tt(i, j)) * 
           b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
        sbots = sbots + 3.0 * sqr(lp3bb(i, j)) *
           b0(p, tree.md(i, 3), tree.md(j, 3), q);
        staus = staus + sqr(lp3tautau(i, j)) * 
           b0(p, tree.me(i, 3), tree.me(j, 3), q);
     }
  }	

  /// Higgs
  //	
  /// Quadrilinear CP-even Higgs couplings 
  DoubleMatrix ssp3p3(3, 3);
  ssp3p3(1, 1) = 0.25 * lsq;
  ssp3p3(2, 2) = ssp3p3(1, 1);
  ssp3p3(1, 2) = 0.25 * lam * kap;
  ssp3p3(3, 3) = 0.5 * ksq;
  ssp3p3.symmetrise();
  
  /// Quadrilinear CP-odd Higgs couplings
  DoubleMatrix ppp3p3(3, 3);
  ppp3p3(1, 1) = lsq / 24.0;
  ppp3p3(2, 2) = ppp3p3(1, 1);
  ppp3p3(1, 2) = - lam * kap / 24.0;
  ppp3p3(3, 3) = 0.25 * ksq;
  ppp3p3.symmetrise();
  
  /// LCT: Rotate to mass bases p3 p3 Hi Hi and p3 p3 Ai Ai
  DoubleVector hhp3p3(3), aap3p3(3);
  for (int i = 1; i <=3; i++) {
    for (int a = 1; a <= 3; a++) {
      for (int b = 1; b <=3; b++) {
	hhp3p3(i) = hhp3p3(i) + S(i, a) * S(i, b) * ssp3p3(a, b);
	aap3p3(i) = aap3p3(i) + 6.0 * P(i, a) * P(i, b) * ppp3p3(a, b);
      }
    }
  }

  DoubleMatrix spp3(3, 3), hphpp3(2, 2);
  getP3HiggsTriCoup(spp3, hphpp3);
  /// LCT: Rotate to mass basis p3 Ai Hj
  DoubleMatrix ahp3(3, 3);
  for (int i=1; i <= 3; i++) {
     for (int j=1; j <=3; j++) {
        for (int a = 1; a <= 3; a++) {
           for (int b = 1; b <= 3; b++) {
              ahp3(i, j) =  ahp3(i, j) + 2.0 * S(j, a) * P(i, b) * spp3(a, b);
           }
        }
     }
  }
  
   double higgs = 0.0;
  for (i=1; i<=3; i++) {
     for (j=1; j<=3; j++) {
        higgs = higgs + sqr(ahp3(i, j)) * b0(p, higgsa(i), higgsm(j), q);
     }
     higgs = higgs + 2.0 * hhp3p3(i) * a0(higgsm(i), q);
     higgs = higgs + 2.0 * aap3p3(i) * a0(higgsa(i), q);
  }
  
  /// Quadrilinear (H1+ H1- H2+ H2-)
  DoubleVector hphpp3p3(2);
  hphpp3p3(1) = 0.5 * lam * (lam + kap * sin2b);
  hphpp3p3(2) = 0.5 * lam * (lam - kap * sin2b);

  for (i=1; i<=2; i++) {
    for (j=1; j<=2; j++) {
      higgs = higgs + sqr(hphpp3(i, j)) * b0(p, higgsc(i), higgsc(j), q);
    }
    higgs = higgs + 2.0 * hphpp3p3(i) * a0(higgsc(i), q);
  }

  /// Neutralino contribution
  double neutralinos = 0.0;
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix aChi(5, 5), bChi(5, 5);
  getP3NeutralinoCoup(aChi, bChi);
  DoubleMatrix fChiChip3p3(5, 5), gChiChip3p3(5, 5);
  for(i=1; i<=5; i++)
    for (j=1; j<=5; j++) {
      fChiChip3p3(i, j) = sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod());
      gChiChip3p3(i, j) = (bChi(i, j).conj() * aChi(i, j) + 
			   aChi(i, j).conj() * bChi(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
           (fChiChip3p3(i, j) * gfn(p, mneut(i), mneut(j), q) + 2.0 *
            gChiChip3p3(i, j) * mneut(i) * mneut(j) * 
            b0(p, mneut(i), mneut(j), q));
    }
  
  /// Chargino contribution
  double chargino = 0.0;
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz); 
  DoubleVector mch(tree.mchBpmz); 
  DoubleMatrix aPsic(2, 2);
  aPsic(2, 2) = lam / root2; 
  ComplexMatrix aChic(2, 2), bChic(2, 2);
  aChic = v.complexConjugate() * aPsic * u.hermitianConjugate();
  bChic = u * aPsic.transpose() * v.transpose();
  for(i=1; i<=2; i++)
     for (j=1; j<=2; j++) {
        fChiChip3p3(i, j) = sqr(aChic(i, j).mod()) + sqr(bChic(i, j).mod());
        gChiChip3p3(i, j) = (bChic(i, j).conj() * aChic(i, j) + 
                             aChic(i, j).conj() * bChic(i, j)).real();
        chargino = chargino + 
           (fChiChip3p3(i, j) * gfn(p, mch(i), mch(j), q) + 2.0 *
            gChiChip3p3(i, j) * mch(i) * mch(j) * 
            b0(p, mch(i), mch(j), q));
     }
	
  double sfermions = staus + stops + sbots; 
  
  return (sfermions + higgs + neutralinos 
	  + chargino) / (16.0 * sqr(PI));

}

//PA: Obtains trilnear couplings of P1-Pi-Sj and P1-Hpm-Hpm 
//for use in loop functions
void NmssmSoftsusy::getP1HiggsTriCoup(DoubleMatrix & spp1, DoubleMatrix & hphpp1, double cw2DRbar) const {
  double beta = atan(displayTanb());
  double lam =  displayLambda(), lsq = sqr(lam);
  double s = displaySvev();
  double kap = displayKappa();
  double al = displayTrialambda();
  double mupr = displayMupr();
  double g = displayGaugeCoupling(2), gsq = sqr(g);
  double cb = cos(beta), sb = sin(beta);
  double v1 =  displayHvev() * cb, v2 =  displayHvev() * sb; 
  
  spp1(1, 1) = 0.125 * gsq / cw2DRbar * v1;
  spp1(2, 1) = 0.25 * v2 * (2.0 * lsq - 0.5 * gsq / cw2DRbar);
  spp1(3, 1) = 0.5 * lsq * s;
  spp1(3, 3) = - 0.5 * lam * kap * v2;
  spp1(2, 3) = 0.5 * (al / root2 - lam * kap * s - lam * mupr / root2);
  spp1(3, 2) = 0.5 * (al / root2 + lam * kap * s + lam * mupr / root2);


  hphpp1(1, 2) = 0.25 * v2 * (2.0 * lsq - gsq);
  hphpp1(2, 1) = -hphpp1(1, 2);	

}


void NmssmSoftsusy::getP2HiggsTriCoup(DoubleMatrix & spp2, DoubleMatrix & hphpp2, double cw2DRbar) const {
  double beta = atan(displayTanb());
  double lam =  displayLambda(), lsq = sqr(lam);
  double s = displaySvev();
  double kap = displayKappa();
  double al = displayTrialambda();
  double mupr = displayMupr();
  double g = displayGaugeCoupling(2), gsq = sqr(g);
  double cb = cos(beta), sb = sin(beta);
  double v1 =  displayHvev() * cb, v2 =  displayHvev() * sb; 

  spp2(1, 2) = 0.25 * v1 * (2.0 * lsq - 0.5 * gsq / cw2DRbar);
  spp2(1, 3) = 0.5 * (al / root2 - lam * kap * s - lam * mupr / root2);
  spp2(2, 2) = 0.125 * gsq / cw2DRbar * v2;
  spp2(3, 1) = 0.5 * (al / root2 + lam * kap * s + lam * mupr / root2);
  spp2(3, 2) = 0.5 * lsq * s;
  spp2(3, 3) = -0.5 * lam * kap * v1;
  
  hphpp2(1, 2) = 0.25 * v1 * (2.0 * lsq - gsq);
  hphpp2(2, 1) = -hphpp2(1, 2);	


}


void NmssmSoftsusy::getP3HiggsTriCoup(DoubleMatrix & spp3, DoubleMatrix & hphpp3) const {
  double beta = atan(displayTanb());
  double lam =  displayLambda();
  double s = displaySvev();
  double kap = displayKappa(), ksq = sqr(kap);
  double al = displayTrialambda();
  double ak = displayTriakappa();
  double mupr = displayMupr();
  double cb = cos(beta), sb = sin(beta);
  double v1 =  displayHvev() * cb, v2 =  displayHvev() * sb; 

  spp3(1, 2) = 0.5 * (al / root2 - lam * kap * s - lam * mupr / root2);
  spp3(1, 3) = 0.5 * lam * (lam * v1 + kap * v2);
  spp3(2, 1) = spp3(1, 2);
  spp3(2, 3) = 0.5 * lam * (lam * v2 + kap * v1);
  spp3(3, 1) = -0.5 * lam * kap * v2;
  spp3(3, 2) = -0.5 * lam * kap * v1;
  spp3(3, 3) = -ak / root2 + ksq * s + kap * mupr / root2;

  hphpp3(1, 2) =  (al / root2 - lam * kap * s) ;
  hphpp3(2, 1) = - hphpp3(1, 2);

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
//PA: obtains CP odd Higgs-Neutralino couplings
void NmssmSoftsusy::getP1NeutralinoCoup(ComplexMatrix & aChi, 
                                        ComplexMatrix & bChi) const {
  double g   = displayGaugeCoupling(2);
  double gp  = displayGaugeCoupling(1) * sqrt(0.6);
  double lam = displayLambda();
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleMatrix aPsi(5, 5);
  aPsi(1, 3) = gp * 0.5; 
  aPsi(2, 3) = -g * 0.5; 
  aPsi(4, 5) = -lam / root2;
  aPsi.symmetrise();
  aChi = n.complexConjugate() * aPsi * n.hermitianConjugate();
  bChi = n * aPsi * n.transpose();
}
void NmssmSoftsusy::getP2NeutralinoCoup(ComplexMatrix & aChi2, 
                                        ComplexMatrix & bChi2) const {
  double g   = displayGaugeCoupling(2);
  double gp  = displayGaugeCoupling(1) * sqrt(0.6);
  double lam = displayLambda();
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleMatrix aPsi2(5, 5);
  aPsi2(1, 4) = -gp * 0.5; 
  aPsi2(2, 4) = g * 0.5;
  aPsi2(3, 5) = -lam / root2; 
  aPsi2.symmetrise();
  aChi2 = n.complexConjugate() * aPsi2 * n.hermitianConjugate();
  bChi2 = n * aPsi2 * n.transpose();
}


void NmssmSoftsusy::getP3NeutralinoCoup(ComplexMatrix & aChi3, 
                                        ComplexMatrix & bChi3) const {
  double lam = displayLambda();
  double kap = displayKappa();
  
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleMatrix aPsi3(5, 5);
  aPsi3(3, 4) = - lam / root2;
  aPsi3(5, 5) = root2 * kap;
  
  aPsi3.symmetrise();
  aChi3 = n.complexConjugate() * aPsi3 * n.hermitianConjugate();
  bChi3 = n * aPsi3 * n.transpose();

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

void NmssmSoftsusy::getHpWmZCoup(DoubleMatrix & HpWmZ) const {
 double g    = displayGaugeCoupling(2);
 double gp = sqrt(0.6) * displayGaugeCoupling(1);
 double mz   = displayMzRun();
 double sW2  = sqr(calcSinthdrbar());
 double v    = displayHvev();
 double beta = atan(displayTanb()), cosb = cos(beta), sinb = sin(beta);

 HpWmZ(1, 1) = 0.5 * gp * g * v * cosb * sin(calcSinthdrbar()); 
 HpWmZ(2, 2) = - g * mz * sW2 * sinb;
}

void NmssmSoftsusy::getHp1HiggsTriCoup(DoubleMatrix & ahphp1, DoubleMatrix & hhphp1) const {
  double gsq = sqr(displayGaugeCoupling(2));
  double gpsq = 0.6 * sqr(displayGaugeCoupling(1));
  double beta = atan(displayTanb());
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta);
  double sinb = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);
  double lam = displayLambda(), lsq = sqr(lam);
  double kap = displayKappa();
  double mupr = displayMupr();
  double al = displayTrialambda();
  double s = displaySvev();
  double vev = displayHvev(), v1 = vev * cosb, v2 = vev * sinb;
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  for (int i=1; i<=2; i++) {
    for (int j=1; j<=3; j++) {
      ahphp1(i, j) = 0.25 * C(i, 2) * (2.0 * P(j, 3) 
		   * (2.0 * s * kap * lam - root2 * al)
	           + v1 * (gsq - 2.0 * lsq) * P(j, 2) 
	           + v2 * (gsq - 2.0 * lsq) * P(j, 1));
      hhphp1(i, j) = - 0.25 * (C(i, 1) * (S(j, 2) * (gsq - gpsq) * v2 
		   + S(j, 1) * (gsq + gpsq) * v1 + 4.0 * S(j, 3) * s * lsq) 
		   + C(i, 2) * (S(j, 1) * v2 * (gsq - 2.0 * lsq) 
	           + S(j, 2) * v1 * (gsq - 2.0 * lsq) 
		   + 2.0 * S(j, 3) * (2.0 * s * lam * kap + root2 * al))) 
	           + root2 * lam * mupr * S(j, 3) * cosb; // General NMSSM piece
    }
  }
}

void NmssmSoftsusy::getHp2HiggsTriCoup(DoubleMatrix & ahphp2, DoubleMatrix & hhphp2) const {
  double gsq = sqr(displayGaugeCoupling(2));
  double gpsq = 0.6 * sqr(displayGaugeCoupling(1));
  double beta = atan(displayTanb());
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta);
  double sinb = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);
  double lam = displayLambda(), lsq = sqr(lam);
  double kap = displayKappa();
  double mupr = displayMupr();
  double al = displayTrialambda();
  double s = displaySvev();
  double vev = displayHvev(), v1 = vev * cosb, v2 = vev * sinb;
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  for (int i=1; i<=2; i++) {
    for (int j=1; j<=3; j++) {
      ahphp2(i, j) = - 0.25 * C(i, 1) 
	           * (2.0 * (2.0 * s * lam * kap - root2 * al) * P(j, 3)
		   + (gsq - 2.0 * lsq) * (v1 * P(j, 2) + v2 * P(j, 1)));
      hhphp2(i, j) = - 0.25 * (C(i, 1) * ((gsq - 2.0 * lsq) 
		   * (v1 * S(j, 2) + v2 * S(j, 1)) 
                   + 2.0 * S(j, 3) * (2.0 * s * lam * kap + root2 * al))
		   + C(i, 2) * (S(j, 1) * v1 * (gsq - gpsq) 
                   + S(j, 2) * (gsq + gpsq) * v2 + 4.0 * S(j, 3) * s * lsq)) 
	           + root2 * lam * mupr * S(j, 3) * sinb;
    }
  }
}

void NmssmSoftsusy::getNeutralinoCharginoHpmCoup(ComplexMatrix & apph1, ComplexMatrix & apph2, ComplexMatrix & bpph1, ComplexMatrix & bpph2) const {
  double gp = sqrt(0.6) * displayGaugeCoupling(1);
  double g = displayGaugeCoupling(2);
  double lam = displayLambda();

  /// Call MSSM couplings
  Softsusy<SoftParsNmssm>::getNeutralinoCharginoHpmCoup(apph1, apph2, bpph1, bpph2);

  /// Fill remaining values
  apph2(5, 2) = - lam;
  bpph1(5, 2) = lam;
}

void NmssmSoftsusy::getHp1NeutCharCoup(ComplexMatrix & lnchp1, ComplexMatrix & rnchp1) const {
  double g    = displayGaugeCoupling(2), gsq = sqr(g);
  double gp = sqrt(0.6) * displayGaugeCoupling(1);
  double lam = displayLambda();
  ComplexMatrix n(displayDrBarPars().nBpmz);
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz);
  ComplexMatrix nconj = n.complexConjugate(), uconj = u.complexConjugate(), 
                vconj = v.complexConjugate();

  lnchp1(1, 1) = - g * vconj(1, 1) * nconj(1, 3) 
    + vconj(1, 2) * (gp * nconj(1, 1) + g * nconj(1, 2)) / root2;
  lnchp1(1, 2) = - g * vconj(1, 1) * nconj(2, 3) 
    + vconj(1, 2) * (gp * nconj(2, 1) + g * nconj(2, 2)) / root2;
  lnchp1(1, 3) = - g * vconj(1, 1) * nconj(3, 3) 
    + vconj(1, 2) * (gp * nconj(3, 1) + g * nconj(3, 2)) / root2;
  lnchp1(1, 4) = - g * vconj(1, 1) * nconj(4, 3) 
    + vconj(1, 2) * (gp * nconj(4, 1) + g * nconj(4, 2)) / root2;
  lnchp1(1, 5) = - g * vconj(1, 1) * nconj(5, 3) 
    + vconj(1, 2) * (gp * nconj(5, 1) + g * nconj(5, 2)) / root2;

  lnchp1(2, 1) = - g * vconj(2, 1) * nconj(1, 3) 
    + vconj(2, 2) * (gp * nconj(1, 1) + g * nconj(1, 2)) / root2;

}

void NmssmSoftsusy::getGaugeHiggsHpmCoup(DoubleVector & wmhhp1, DoubleVector & wmahp1, DoubleVector & gamhphp1, DoubleVector & zhphp1, DoubleVector & wmhhp2, DoubleVector & wmahp2, DoubleVector & gamhphp2, DoubleVector & zhphp2) const {
  double g    = displayGaugeCoupling(2);
  double gp = sqrt(0.6) * displayGaugeCoupling(1);
  double thetaW = asin(calcSinthdrbar());
  double sW = calcSinthdrbar(), cW = cos(thetaW);
  double beta = atan(displayTanb());
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta);
  double sinb = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
 
  for (int i=1; i<=3; i++) {
    wmhhp1(i) = - 0.5 * g * S(i, 1);
    wmahp1(i) = 0.5 * g * P(i, 1);
    wmhhp2(i) = 0.5 * g * S(i, 2);
    wmahp2(i) = 0.5 * g * P(i, 2);
  }
  for (int i=1; i<=2; i++) {
    gamhphp1(i) = 0.5 * (gp * cW + g * sW) * C(i, 1);
    zhphp1(i) = 0.5 * (-gp * sW + g * cW) * C(i, 1);
    gamhphp2(i) = 0.5 * (gp * cW + g * sW) * C(i, 2);
    zhphp2(i) = 0.5 * (-gp * sW + g * cW) * C(i, 2);
  }
}

double NmssmSoftsusy::piHpm11Higgs(double p, double q) const {
  double gsq = sqr(displayGaugeCoupling(2));
  double gpsq = 0.6 * sqr(displayGaugeCoupling(1));
  double beta    = atan(displayTanb());
  double cosb    = cos(beta);
  double sinb    = sin(beta);
  double lam = displayLambda(), lsq = sqr(lam);
  double al = displayTrialambda();
  double mz = displayMzRun(), mw = displayMwRun();

  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices
  ///
  /// NB: For consistency with the NMSSM conventions of Degrassi and Slavich 
  /// (Nucl. Phys. B825, 119 (2010); arXiv:0907.4682 [hep-ph]), the charged 
  /// Higgs mixing matrix C is taken to be in the mass basis (G+ H+):
  ///
  /// C = [-cos(beta) sin(beta) ]
  ///     [ sin(beta) cos(beta) ]
  ///
  /// Note that compared with BPMZ and MSSM SOFTSUSY, this convention 
  /// introduces a relative minus *sign* in the left-handed 
  /// Neutralino-Chargino-{h1+ h2+} couplings bpph1 and bpph2 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
   /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Quadrilinear charged Higgs couplings (from Staub et al.)
  DoubleVector hp1hp1(2);
  hp1hp1(1) = 0.25 * (- 2.0 * (gsq + gpsq) * sqr(C(1, 1)) 
			+ (gsq + gpsq - 4.0 * lsq) * sqr(C(1, 2)));
  hp1hp1(2) = 0.25 * (- 2.0 * (gsq + gpsq) * sqr(C(2, 1)) 
			+ (gsq + gpsq - 4.0 * lsq) * sqr(C(2, 2)));

  double higgs = 0.0;
  for(int i=1; i<=2; i++) {
      higgs = higgs - hp1hp1(i) * a0(higgsc(i), q);
  }

 /// LCT: Trilinear charged Higgs with CP-odd/-even Higgs hp-A0/h0-Hm
  DoubleMatrix ahphp1(2, 3), hhphp1(2, 3);
  getHp1HiggsTriCoup(ahphp1, hhphp1);

  for (int i=1; i<=2; i++) {
     for (int j=1; j<=3; j++) {
        higgs = higgs + sqr(ahphp1(i, j)) * b0(p, higgsc(i), higgsa(j), q);
	higgs = higgs + sqr(hhphp1(i, j)) * b0(p, higgsc(i), higgsm(j), q);
     }
  }

  /// LCT: Quadrilinear A0A0-hphp and h0h0-hphp couplings
  DoubleVector aahp1hp1(3), hhhp1hp1(3);
  for (int i=1; i<=3; i++) {
    aahp1hp1(i) = - 0.25 * (4.0 * lsq * sqr(P(i, 3)) 
		+ (gsq - gpsq) * sqr(P(i, 2)) + (gsq + gpsq) * sqr(P(i, 1)));
    hhhp1hp1(i) = -0.25 * (4.0 * lsq * sqr(S(i, 3)) 
		+ (gsq - gpsq) * sqr(S(i, 2)) + (gsq + gpsq) * sqr(S(i, 1)));
  }
  for (int i=1; i<=3; i++) {
    higgs = higgs - 0.5 * aahp1hp1(i) * a0(higgsa(i), q);
    higgs = higgs - 0.5 * hhhp1hp1(i) * a0(higgsm(i), q);
  }

  /// LCT: F0 function contributions
  DoubleVector wmhhp1(3), wmahp1(3), gamhphp1(2), zhphp1(2), wmhhp2(3),
    wmahp2(3), gamhphp2(2), zhphp2(2);
  getGaugeHiggsHpmCoup(wmhhp1, wmahp1, gamhphp1, zhphp1, wmhhp2, wmahp2,
  gamhphp2, zhphp2);

  for (int i=1; i<=2; i++) {
    higgs = higgs + sqr(gamhphp1(i)) * ffn(p, higgsc(i), 0, q)
      + sqr(zhphp1(i)) * ffn(p, higgsc(i), mz, q);
  }
  for (int i=1; i<=3; i++) {
    higgs = higgs + sqr(wmahp1(i)) * ffn(p, higgsa(i), mw, q)
      + sqr(wmhhp1(i)) * ffn(p, higgsm(i), mw, q);
  }
  return higgs;
}

double NmssmSoftsusy::piHpm12Higgs(double p, double q) const {
  double gsq = sqr(displayGaugeCoupling(2));
  double gpsq = 0.6 * sqr(displayGaugeCoupling(1));
  double beta    = atan(displayTanb());
  double cosb    = cos(beta);
  double sinb    = sin(beta);
  double lam = displayLambda(), lsq = sqr(lam);
  double kap = displayKappa();
  double al = displayTrialambda();
  double mz = displayMzRun(), mw = displayMwRun();

 /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
   /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Quadrilinear charged Higgs couplings (from Staub et al.)
  DoubleVector hp1hp2(2);
  hp1hp2(1) = 0.25 * (gsq + gpsq - 4.0 * lsq) * C(1, 2) * C(1, 1);
  hp1hp2(2) = 0.25 * (gsq + gpsq - 4.0 * lsq) * C(2, 2) * C(2, 1);

  double higgs = 0.0;
  for(int i=1; i<=2; i++) {
      higgs = higgs - hp1hp2(i) * a0(higgsc(i), q);
  }

 /// LCT: Trilinear charged Higgs with CP-odd/-even Higgs hp-A0/h0-Hm
  DoubleMatrix ahphp1(2, 3), hhphp1(2, 3), ahphp2(2, 3), hhphp2(2, 3);
  getHp1HiggsTriCoup(ahphp1, hhphp1);
  getHp2HiggsTriCoup(ahphp2, hhphp2);

  for (int i=1; i<=2; i++) {
     for (int j=1; j<=3; j++) {
       higgs = higgs + ahphp1(i, j) * ahphp2(i, j) 
             * b0(p, higgsc(i), higgsa(j), q);
       higgs = higgs + hhphp1(i, j) * hhphp2(i, j) 
       	     * b0(p, higgsc(i), higgsm(j), q);
     }
  }

  /// LCT: Quadrilinear A0A0-hphp and h0h0-hphp couplings
  DoubleVector aahp1hp2(3), hhhp1hp2(3);
  for (int i=1; i<=3; i++) {
    aahp1hp2(i) = 0.25 * (4.0 * lam * kap * sqr(P(i, 3)) 
                + 2.0 * (gsq - 2.0 * lsq) * P(i, 1) * P(i, 2));
    hhhp1hp2(i) = - 0.25 * (4.0 * lam * kap * sqr(S(i, 3)) 
		+ 2.0 * (gsq - 2.0 * lsq) * S(i, 1) * S(i, 2));
  }
  for (int i=1; i<=3; i++) {
    higgs = higgs - 0.5 * aahp1hp2(i) * a0(higgsa(i), q);
    higgs = higgs - 0.5 * hhhp1hp2(i) * a0(higgsm(i), q);
  }

  /// LCT: F0 function contributions
  DoubleVector wmhhp1(3), wmahp1(3), gamhphp1(2), zhphp1(2), wmhhp2(3),
    wmahp2(3), gamhphp2(2), zhphp2(2);
  getGaugeHiggsHpmCoup(wmhhp1, wmahp1, gamhphp1, zhphp1, wmhhp2, wmahp2,
  gamhphp2, zhphp2);

  for (int i=1; i<=2; i++) {
    higgs = higgs + gamhphp1(i) * gamhphp2(i) * ffn(p, higgsc(i), 0, q)
          + zhphp1(i) * zhphp2(i) * ffn(p, higgsc(i), mz, q);
  }
  for (int i=1; i<=3; i++) {
    higgs = higgs + wmahp1(i) * wmahp2(i) * ffn(p, higgsa(i), mw, q)
      + wmhhp1(i) * wmhhp2(i) * ffn(p, higgsm(i), mw, q);
  }

  return higgs;
}

double NmssmSoftsusy::piHpm22Higgs(double p, double q) const {
  double gsq = sqr(displayGaugeCoupling(2));
  double gpsq = 0.6 * sqr(displayGaugeCoupling(1));
  double beta    = atan(displayTanb());
  double cosb    = cos(beta);
  double sinb    = sin(beta);
  double lam = displayLambda(), lsq = sqr(lam);
  double al = displayTrialambda();
  double mz = displayMzRun(), mw = displayMwRun();

 /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3), C(2, 2);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;
   /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector higgsm(3), higgsa(3), higgsc(2);
  assignHiggs(higgsm, higgsa, higgsc);

  /// LCT: Quadrilinear charged Higgs couplings (from Staub et al.)
  DoubleVector hp2hp2(2);
  hp2hp2(1) = 0.25 * (- 2.0 * (gsq + gpsq) * sqr(C(1, 2)) 
            + (gsq + gpsq - 4.0 * lsq) * sqr(C(1, 1)));
  hp2hp2(2) = 0.25 * (- 2.0 * (gsq + gpsq) * sqr(C(2, 2)) 
	    + (gsq + gpsq - 4.0 * lsq) * sqr(C(2, 1)));

  double higgs = 0.0;
  for(int i=1; i<=2; i++) {
      higgs = higgs - hp2hp2(i) * a0(higgsc(i), q);
  }

 /// LCT: Trilinear charged Higgs with CP-odd/-even Higgs hp-A0/h0-Hm
  DoubleMatrix ahphp2(2, 3), hhphp2(2, 3);
  getHp2HiggsTriCoup(ahphp2, hhphp2);

  for (int i=1; i<=2; i++) {
     for (int j=1; j<=3; j++) {
        higgs = higgs + sqr(ahphp2(i, j)) * b0(p, higgsc(i), higgsa(j), q);
	higgs = higgs + sqr(hhphp2(i, j)) * b0(p, higgsc(i), higgsm(j), q);
     }
  }

  /// LCT: Quadrilinear A0A0-hphp and h0h0-hphp couplings
  DoubleVector aahp2hp2(3), hhhp2hp2(3);
  for (int i=1; i<=3; i++) {
    aahp2hp2(i) = - 0.25 * (4.0 * lsq * sqr(P(i, 3)) 
		+ (gsq - gpsq) * sqr(P(i, 1)) + (gsq + gpsq) * sqr(P(i, 2)));
    hhhp2hp2(i) = - 0.25 * (4.0 * lsq * sqr(S(i, 3)) 
	        + (gsq - gpsq) * sqr(S(i, 1)) + (gsq + gpsq) * sqr(S(i, 2)));
  }
  for (int i=1; i<=3; i++) {
    higgs = higgs - 0.5 * aahp2hp2(i) * a0(higgsa(i), q);
    higgs = higgs - 0.5 * hhhp2hp2(i) * a0(higgsm(i), q);
  }

  /// LCT: F0 function parts
  DoubleVector wmhhp1(3), wmahp1(3), gamhphp1(2), zhphp1(2), wmhhp2(3),
    wmahp2(3), gamhphp2(2), zhphp2(2);
  getGaugeHiggsHpmCoup(wmhhp1, wmahp1, gamhphp1, zhphp1, wmhhp2, wmahp2,
  gamhphp2, zhphp2);

  for (int i=1; i<=2; i++) {
    higgs = higgs + sqr(gamhphp2(i)) * ffn(p, higgsc(i), 0, q)
      + sqr(zhphp2(i)) * ffn(p, higgsc(i), mz, q);
  }
  for (int i=1; i<=3; i++) {
    higgs = higgs + sqr(wmahp2(i)) * ffn(p, higgsa(i), mw, q)
          + sqr(wmhhp2(i)) * ffn(p, higgsm(i), mw, q);
  }

  return higgs;
}

double NmssmSoftsusy::piHpm11Gaugino(double p, double q) const {
  const drBarPars& tree = displayDrBarPars();
  double gauginos = 0.0;

  /// LCT: Get trilinear couplings
  ComplexMatrix apph1(5, 2), apph2(5, 2), bpph1(5, 2), bpph2(5, 2);
  getNeutralinoCharginoHpmCoup(apph1, apph2, bpph1, bpph2);

  /// LCT: Rotate to physical gaugino e-states
  apph1 = tree.nBpmz.complexConjugate() 
        * apph1 * tree.uBpmz.hermitianConjugate();
  bpph1 = - tree.nBpmz * bpph1 * tree.vBpmz.transpose(); ///<< See piHpm11Higgs 
  /// for origin of minus sign here

  /// LCT: Get neutralino and chargino masses
  DoubleVector mneut(tree.mnBpmz), mch(tree.mchBpmz);

  /// LCT: Contributions start here
  for(int i=1; i<=5; i++) {
    for(int j=1; j<=2; j++) {
      gauginos = gauginos - 2.0 * mch(j) * mneut(i) 
	* 2.0 * (bpph1(i, j).conj() * apph1(i, j)).real()
        * b0(p, mch(j), mneut(i), q) 
        + (sqr(apph1(i, j).mod()) + sqr(bpph1(i, j).mod())) 
	* gfn(p, mch(j), mneut(i), q);
    }
  }
  return gauginos;
}

double NmssmSoftsusy::piHpm12Gaugino(double p, double q) const {
  const drBarPars& tree = displayDrBarPars();
  double gauginos = 0.0;

  /// LCT: Get trilinear couplings
  ComplexMatrix apph1(5, 2), apph2(5, 2), bpph1(5, 2), bpph2(5, 2);
  getNeutralinoCharginoHpmCoup(apph1, apph2, bpph1, bpph2);

  /// LCT: Rotate to physical gaugino e-states
  apph1 = tree.nBpmz.complexConjugate() 
        * apph1 * tree.uBpmz.hermitianConjugate();
  apph2 = tree.nBpmz.complexConjugate() 
        * apph2 * tree.uBpmz.hermitianConjugate();
  bpph1 = - tree.nBpmz * bpph1 * tree.vBpmz.transpose();
  bpph2 = - tree.nBpmz * bpph2 * tree.vBpmz.transpose(); ///<< See piHpm11Higgs 
  /// for origin of minus sign here

  /// LCT: Get neutralino and chargino masses
  DoubleVector mneut(tree.mnBpmz), mch(tree.mchBpmz);

  /// LCT: Contributions start here
  for(int i=1; i<=5; i++) {
    for(int j=1; j<=2; j++) {
      gauginos = gauginos - 2.0 * mch(j) * mneut(i) 
	* (apph2(i, j).conj() * bpph1(i, j) 
        + bpph2(i, j).conj() * apph1(i, j)).real() 
	* b0(p, mch(j), mneut(i), q) 
        + (apph2(i, j).conj() * apph1(i, j) 
	+ bpph2(i, j).conj() * bpph1(i, j)).real() 
	* gfn(p, mch(j), mneut(i), q);
    }
  }
  return gauginos;
}

double NmssmSoftsusy::piHpm22Gaugino(double p, double q) const {
  const drBarPars& tree = displayDrBarPars();
  double gauginos = 0.0;

  /// LCT: Get trilinear couplings
  ComplexMatrix apph1(5, 2), apph2(5, 2), bpph1(5, 2), bpph2(5, 2);
  getNeutralinoCharginoHpmCoup(apph1, apph2, bpph1, bpph2);

  /// LCT: Rotate to physical gaugino e-states
  apph2 = tree.nBpmz.complexConjugate() 
        * apph2 * tree.uBpmz.hermitianConjugate();
  bpph2 = - tree.nBpmz * bpph2 * tree.vBpmz.transpose(); ///<< See piHpm11Higgs 
  /// for origin of minus sign here

  /// LCT: Get neutralino and chargino masses
  DoubleVector mneut(tree.mnBpmz), mch(tree.mchBpmz);

  /// LCT: Contributions start here
  for(int i=1; i<=5; i++) {
    for(int j=1; j<=2; j++) {
      gauginos = gauginos - 2.0 * mch(j) * mneut(i) 
	* 2.0 * (bpph2(i, j).conj() * apph2(i, j)).real()
        * b0(p, mch(j), mneut(i), q)
        + (sqr(apph2(i, j).mod()) + sqr(bpph2(i, j).mod())) 
	* gfn(p, mch(j), mneut(i), q);
    }
  }
  return gauginos;
}

/// LCT: Charged Higgs self-energy in mass basis
double NmssmSoftsusy::piHpHm(double p, double q) const {
  double lam = displayLambda();
  double svev = displaySvev();
  double mueff = -(displaySusyMu() + lam * svev / root2); ///<< Note minus sign
  double beta = atan(displayTanb());
  double cosb = cos(beta),  sinb = sin(beta);
  DoubleMatrix C(2, 2);
  C(1, 1) = - cosb;  C(1, 2) = sinb; 
  C(2, 1) = C(1, 2); C(2, 2) = cosb;

  double gauge = piHpHmGauge(p, q);
  double fermions = piHpHmFermions(p, q);
  double sfermions = piHpHmSfermions(p, q, mueff);
  double higgs = sqr(sinb) * piHpm11Higgs(p, q) + sqr(cosb) * piHpm22Higgs(p, q)               + 2.0 * sinb * cosb * piHpm12Higgs(p, q);
  double gauginos = sqr(sinb) * piHpm11Gaugino(p, q) 
                  + sqr(cosb) * piHpm22Gaugino(p, q)  
                  + 2.0 * sinb * cosb * piHpm12Gaugino(p, q);

  return (gauge + fermions + sfermions + higgs + gauginos) 
    / (16.0 * sqr(PI));
}


double NmssmSoftsusy::calcRunMtHiggs() const {
  double mtpole = displayDataSet().displayPoleMt();
  double p  = mtpole;
  double q  = displayMu();
  double ht = displayDrBarPars().ht;
  double hb = displayDrBarPars().hb;
  double mz = displayMzRun();
  double mw = displayMwRun();
  double mt = displayDrBarPars().mt;
  double mb = displayDrBarPars().mb;
  double g  = displayGaugeCoupling(2);
  double e  = g * calcSinthdrbar();
  double beta = atan(displayTanb());
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta);
  double sinb = sin(beta);
  double thetaWDRbar = asin(calcSinthdrbar());
  double cw2DRbar    = sqr(cos(thetaWDRbar));
  const double  costh   = (displayMw() / displayMz());
  const double    cw2   = sqr(costh) ;
  const double    sw2   = (1.0 - cw2);
  double gtL = 0.5 - 2.0 * sw2 / 3.0, gtR = 2.0 * sw2 / 3.0;
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector mH(3), higgsa(3), higgsc(2);
  assignHiggs(mH, higgsa, higgsc);
  double mHc = higgsc(2);
 /// CP-even Higgs
  double higgs =0.0;
  for (int i = 1; i <= 3; i++) {
     higgs = higgs + 0.5 * sqr(ht) * sqr(S(i,2)) * (b1(p, mt, mH(i), q) 
						    + b0(p, mt, mH(i), q));
  }
  
  //CP-odd Higgs and neutral goldstone
  for (int i = 1; i <= 3; i++) {
    higgs = higgs + 0.5 * sqr(ht) * sqr(P(i, 2)) * (b1(p, mt, higgsa(i), q) 
						    - b0(p, mt, higgsa(i), q));
  }  
  /// Charged Higgs and charged Goldstones contribution
  higgs = higgs + 0.5 * ((sqr(hb) * sqr(sinb) + sqr(ht) * sqr(cosb)) *
			 b1(p, mb, mHc, q) 
                         + (sqr(g) + sqr(hb) * sqr(cosb) + sqr(ht) * sqr(sinb)) 
                         * b1(p, mb, mw, q))
			 + sqr(hb) * sqr(cosb) 
                         * (b0(p, mb, mHc, q) - b0(p, mb, mw, q));
  
  /// Gauge bosons
  higgs = higgs - sqr(e) * 4.0 / 9.0 * (5.0 + 6.0 * log(q / mt)) 
     + sqr(g) / cw2DRbar * ( (sqr(gtL) + sqr(gtR)) * b1(p, mt, mz, q)
                             + 4.0 * gtL * gtR * b0(p, mt, mz, q) );
  return higgs; 
}



double NmssmSoftsusy::calcRunMtNeutralinos() const {
  double q       = displayMu();
  double mtpole  = displayDataSet().displayPoleMt();
  double ht      = displayDrBarPars().ht;
  double g       = displayGaugeCoupling(2);
  double gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double thetat  = displayDrBarPars().thetat ;
  double p = mtpole;
  DoubleMatrix neutralinoContribution(5, 2);
  double neutralinos = 0.0;
  
  /// Neutralino contribution
  DoubleVector aPsi0TStopr(5), bPsi0TStopr(5), aPsi0TStopl(5),
    bPsi0TStopl(5); 
  aPsi0TStopr(1) = - 4 * gp / (3.0 * root2);
  bPsi0TStopl(1) = gp / (3.0 * root2);
  bPsi0TStopl(2) = g / root2;
  aPsi0TStopl(4) = ht;
  bPsi0TStopr(4) = ht;

  ComplexVector aChi0TStopl(5), bChi0TStopl(5), aChi0TStopr(5),
    bChi0TStopr(5);

  /// Neutralinos
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);

  aChi0TStopl = n.complexConjugate() * aPsi0TStopl;
  bChi0TStopl = n * bPsi0TStopl;
  aChi0TStopr = n.complexConjugate() * aPsi0TStopr;
  bChi0TStopr = n * bPsi0TStopr;

  ComplexMatrix aNeutTStop(5, 2), bNeutTStop(5, 2);
  DoubleMatrix fNeutTStop(5, 2), gNeutTStop(5, 2);

  int i, j;
  DoubleMatrix O(2, 2);
  ComplexVector tt(2), t1(2), t2(2);
  O = rot2d(thetat);
  for (i=1; i<=5; i++) {
    tt(1) = aChi0TStopl(i); tt(2) = aChi0TStopr(i);      
    t1 = O * tt;

    tt(1) = bChi0TStopl(i); tt(2) = bChi0TStopr(i);      
    t2 = O * tt;    
    for (j=1; j<=2; j++) {
      aNeutTStop(i, j) = t1(j);
      bNeutTStop(i, j) = t2(j);
      /// functions of couplings needed for loops
      fNeutTStop(i, j) = sqr(aNeutTStop(i, j).mod()) + 
	sqr(bNeutTStop(i, j).mod());

      gNeutTStop(i, j) = 2.0 * 
	(aNeutTStop(i, j) * bNeutTStop(i, j).conj()).real(); 
      
      neutralinoContribution(i, j) = (fNeutTStop(i, j) * 
	 b1(p, mneut(i), displayDrBarPars().mu(j, 3), q) + 
	 gNeutTStop(i, j) * mneut(i) /  mtpole *  
	 b0(p, mneut(i), displayDrBarPars().mu(j, 3), q)) * 0.5;
    
      neutralinos = neutralinos + neutralinoContribution(i, j);
    }
  }


return neutralinos;

}

double NmssmSoftsusy::calcRunningMt() {
  double mtpole  = displayDataSet().displayPoleMt();
  double resigmat = 0.0; 
  double qcd = 0.0, stopGluino = 0.0, higgs = 0.0; 
  //one and two loop qcd
  qcd = Softsusy<SoftParsNmssm>::calcRunMtQCD();
  resigmat = resigmat + qcd;
  /// stop/gluino correction 6% correction
  stopGluino = Softsusy<SoftParsNmssm>::calcRunMtStopGluino();
  resigmat = resigmat + stopGluino;
  /// rest are extra bits from Matchev et al: 2% corrections  
  //Higgs contribution
  higgs = calcRunMtHiggs();
  resigmat = resigmat + higgs;
  /// Neutralino contribution
  double neutralinos = calcRunMtNeutralinos();
  resigmat = resigmat + neutralinos;
  // Chargino contribution
  double charginoContribution = Softsusy<SoftParsNmssm>::calcRunMtCharginos();
  resigmat = resigmat + charginoContribution; 
  
  resigmat = resigmat * mtpole / (16.0 * sqr(PI));  

  return mtpole + resigmat;
}

double NmssmSoftsusy::calcRunMbHiggs() const {
   double mbMZ = displayDataSet().displayMass(mBottom);
  double p = mbMZ;
  double q = displayMu();
  double hb = displayDrBarPars().hb;
  double ht = displayDrBarPars().ht;
  double mb = displayDrBarPars().mb;
  double mt = displayDrBarPars().mt;
  double  cosb    = cos(atan(displayTanb()));
  double  sinb    = sin(atan(displayTanb()));
  double  mz = displayMzRun();
  double  mw = displayMwRun();
  double  thetaWDRbar = asin(calcSinthdrbar());
  double  cw2DRbar    = sqr(cos(thetaWDRbar));
  double g  = displayGaugeCoupling(2);
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector mH(3), higgsa(3), higgsc(2);
  assignHiggs(mH, higgsa, higgsc);
  double mHp = higgsc(2);
  double deltaHiggs = 0.;
  /// CP-even Higgs
  for (int i = 1; i <= 3; i++) {
     deltaHiggs = deltaHiggs 
        + 0.5 * sqr(hb) * (sqr(S(i,1)) * (b1(p, mb, mH(i), q) 
                                          + b0(p, mb, mH(i), q)));
	}
  	
  /// CP-odd Higgs
  for (int i = 1; i <= 3; i++) {
     deltaHiggs = deltaHiggs 
        + 0.5 * sqr(hb) * (sqr(P(i, 1)) * (b1(p, mb, higgsa(i), q) 
                                           - b0(p, mb, higgsa(i), q)));
  }
  /// Charged Higgs/Goldstone + gauge boson contribution
  deltaHiggs =  deltaHiggs 
     + 0.5 * ((sqr(ht) * sqr(cosb) + sqr(hb) * sqr(sinb)) * b1(p, mt, mHp, q) 
              + (sqr(g) + sqr(ht) * sqr(sinb) + sqr(hb) * sqr(cosb)) 
              * b1(p, mt, mw, q)) 
     + sqr(ht) * sqr(sinb) * (b0(p, mt, mHp, q) - b0(p, mt, mw, q)) 
     + sqr(g) / cw2DRbar * ((sqr(gdL) + sqr(gdR)) * b1(p, mb, mz, q) 
                            + 4.0 * gdL * gdR * b0(p, mb, mz, q));
  
   deltaHiggs = - deltaHiggs / (16.0 * sqr(PI));
  return deltaHiggs;
}

double NmssmSoftsusy::calcRunMbNeutralinos() const {
   double mbMZ = displayDataSet().displayMass(mBottom);
  double p = mbMZ;
  double q = displayMu();
  double thetab  = displayDrBarPars().thetab;
  double g       = displayGaugeCoupling(2);
  double g1      = displayGaugeCoupling(1);
  double gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double hb      = displayDrBarPars().hb;
  double mbMSSM  = displayDrBarPars().mb;
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);

  DoubleVector aPsi0Bsbotr(5), bPsi0Bsbotr(5), aPsi0Bsbotl(5),
     bPsi0Bsbotl(5); 
  aPsi0Bsbotr(1) = gp / (root2 * 3.0) * 2.0;
  bPsi0Bsbotl(1) = gp / (root2 * 3.0);
  bPsi0Bsbotl(2) = -root2 * g * 0.5;
  aPsi0Bsbotl(3) = hb;
  bPsi0Bsbotr(3) = hb;

  ComplexVector aChi0Bsbotl(5), bChi0Bsbotl(5), aChi0Bsbotr(5),
    bChi0Bsbotr(5);

  double deltaNeutralino = 0.;
  aChi0Bsbotl = n.complexConjugate() * aPsi0Bsbotl;
  bChi0Bsbotl = n * bPsi0Bsbotl;
  aChi0Bsbotr = n.complexConjugate() * aPsi0Bsbotr;
  bChi0Bsbotr = n * bPsi0Bsbotr;

  ComplexMatrix aNeutBsbot(5, 2), bNeutBsbot(5, 2);
  DoubleMatrix fNeutBsbot(5, 2), gNeutBsbot(5, 2), 
    neutralinoContribution(5, 2);
  int i, j; DoubleMatrix O(2, 2); ComplexVector tt(2), t1(2), t2(2);
  O = rot2d(thetab);
  for (i=1; i<=5; i++) {
    tt(1) = aChi0Bsbotl(i); tt(2) = aChi0Bsbotr(i);      
    t1 = O * tt;

    tt(1) = bChi0Bsbotl(i); tt(2) = bChi0Bsbotr(i);      
    t2 = O * tt;    
    for (j=1; j<=2; j++) {
      aNeutBsbot(i, j) = t1(j);
      bNeutBsbot(i, j) = t2(j);
      /// functions of couplings needed for loops
      fNeutBsbot(i, j) = sqr(aNeutBsbot(i, j).mod()) + 
	sqr(bNeutBsbot(i, j).mod());

      gNeutBsbot(i, j) = 2.0 * 
	(aNeutBsbot(i, j) * bNeutBsbot(i, j).conj()).real(); 
      
      neutralinoContribution(i, j) = (fNeutBsbot(i, j) * 
	 b1(p, mneut(i), displayDrBarPars().md(j, 3), q) + 
	 gNeutBsbot(i, j) * mneut(i) /  mbMSSM *  
	 b0(p, mneut(i), displayDrBarPars().md(j, 3), q)) * 0.5;

      deltaNeutralino = deltaNeutralino + neutralinoContribution(i, j);
    }
  }

  deltaNeutralino = -deltaNeutralino / (16.0 * sqr(PI));
  return deltaNeutralino;
}
double NmssmSoftsusy::calcRunningMb() const {

  if (displayMu() != displayMz()) {
    ostringstream ii;
    ii << "Softsusy<SoftPars>::calcRunningMb called with mu=" <<
      displayMu() << endl; 
    throw ii.str();
  }
  
  double mbMZ = displayDataSet().displayMass(mBottom);
  /// First convert mbMZ into DRbar value from hep-ph/9703293,0207126,9701308
  /// (SM gauge boson contributions)
  mbMZ = mbMZ * Softsusy<SoftParsNmssm>::calcRunMbDrBarConv(); 

  double deltaSquarkGluino = Softsusy<SoftParsNmssm>::calcRunMbSquarkGluino();
  //Chargino-squark loops
  double deltaSquarkChargino = Softsusy<SoftParsNmssm>::calcRunMbChargino();
  /// Higgs
  double deltaHiggs = calcRunMbHiggs();
  /// Neutralinos
  double deltaNeutralino = calcRunMbNeutralinos();

  /// it's NOT clear if this resummation is reliable in the full 1-loop scheme
  /// but it's at least valid to 1 loop. Warning though: if you add higher
  /// loops, you'll have to re-arrange.
  return mbMZ / (1.0 + deltaSquarkGluino + deltaSquarkChargino + deltaHiggs
		 + deltaNeutralino);
}

double NmssmSoftsusy::calcRunMtauHiggs() const {
  double mTauPole = MTAU;
  double p = mTauPole;
  double q = displayMu();
  double mz  = displayMzRun();
  double mw  = displayMwRun();
  double cosb = cos(atan(displayTanb()));
  double sinb = sin(atan(displayTanb()));
  double htau = displayDrBarPars().htau;
  double mtau = displayDrBarPars().mtau;
  double g    = displayGaugeCoupling(2);
  double thetaWDRbar = asin(calcSinthdrbar());
  double cw2DRbar    = sqr(cos(thetaWDRbar));
  double mnu = 0.;
  double sigmaHiggs = 0.0;
  /// LCT: Higgs 3 x 3 CP-even S, CP-odd P, and charged C mixing matrices 
  DoubleMatrix P(3, 3), S(3, 3);
  DegrassiSlavicMix(P);
  S = displayDrBarPars().mixh0;
  /// Define Higgs vector in 't-Hooft Feynman gauge:
  DoubleVector mH(3), higgsa(3), higgsc(2);
  assignHiggs(mH, higgsa, higgsc);
  double mHp = higgsc(2);

  /// CP-even Higgs
  for (int i = 1; i <= 3; i++) {
     sigmaHiggs = sigmaHiggs 
        + 0.5 * sqr(htau) * (sqr(S(i,1)) * (b1(p, mtau, mH(i), q) 
                                            + b0(p, mtau, mH(i), q)));
  }
  for (int i = 1; i <= 3; i++) {
     sigmaHiggs = sigmaHiggs 
        + 0.5 * sqr(htau) * (sqr(P(i, 1)) * (b1(p, mtau, higgsa(i), q) 
                                               -  b0(p, mtau, higgsa(i), q)));
  }
  
  /// Charged Higgs + SM contribution
  sigmaHiggs = sigmaHiggs + 0.5 * (sqr(htau) * sqr(sinb) * b1(p, mnu, mHp, q) + (sqr(g) + sqr(htau) * sqr(cosb)) * b1(p, mnu, mw, q)) +
     sqr(g) / cw2DRbar * ((sqr(geL) + sqr(geR)) * b1(p, mtau, mz, q) 
                          +  4.0 * geL * geR * b0(p, mtau, mz, q));
 
  sigmaHiggs = sigmaHiggs / (16.0 * sqr(PI));

  return  sigmaHiggs;
}

double NmssmSoftsusy::calcRunMtauNeutralinos(double mTauSMMZ) const {
  double thetatau = displayDrBarPars().thetatau;
  double g        = displayGaugeCoupling(2);
  double gp       = displayGaugeCoupling(1) * sqrt(0.6);
  double htau = displayDrBarPars().htau;
  double mtau = displayDrBarPars().mtau;
  double mTauPole = MTAU;
  double p = mTauPole;
  double q = displayMu();
  
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
 /// Neutralinos
  DoubleVector aPsi0TauStaur(5), bPsi0TauStaur(5), aPsi0TauStaul(5),
    bPsi0TauStaul(5); 
  aPsi0TauStaur(1) = gp / root2 * 2.0;
  bPsi0TauStaul(1) = -gp / root2;
  bPsi0TauStaul(2) = -root2 * g * 0.5;
  aPsi0TauStaul(3) = htau;
  bPsi0TauStaur(3) = htau;

  ComplexVector aChi0TauStaul(5), bChi0TauStaul(5), aChi0TauStaur(5),
    bChi0TauStaur(5);

  double sigmaNeutralino = 0.;
  aChi0TauStaul = n.complexConjugate() * aPsi0TauStaul;
  bChi0TauStaul = n * bPsi0TauStaul;
  aChi0TauStaur = n.complexConjugate() * aPsi0TauStaur;
  bChi0TauStaur = n * bPsi0TauStaur;

  ComplexMatrix aNeutTauStau(5, 2), bNeutTauStau(5, 2);
  DoubleMatrix fNeutTauStau(5, 2), gNeutTauStau(5, 2), 
    neutralinoContribution(5, 2);
  
  DoubleMatrix O(2, 2); O = rot2d(thetatau);
  ComplexVector t1(2), t2(2), tt(2);
  for (int i=1; i<=5; i++) {
    tt(1) = aChi0TauStaul(i); tt(2) = aChi0TauStaur(i);      
    t1 = O * tt;

    tt(1) = bChi0TauStaul(i); tt(2) = bChi0TauStaur(i);      
    t2 = O * tt;    
    int j;
    for (int j=1; j<=2; j++) {
      aNeutTauStau(i, j) = t1(j);
      bNeutTauStau(i, j) = t2(j);
      /// functions of couplings needed for loops
      fNeutTauStau(i, j) = sqr(aNeutTauStau(i, j).mod()) + 
	sqr(bNeutTauStau(i, j).mod());

      gNeutTauStau(i, j) = 2.0 * 
	(aNeutTauStau(i, j) * bNeutTauStau(i, j).conj()).real(); 
      
      neutralinoContribution(i, j) = (fNeutTauStau(i, j) * 
	 b1(p, mneut(i), displayDrBarPars().me(j, 3), q) + 
	 gNeutTauStau(i, j) * mneut(i) /  mTauSMMZ *  
	 b0(p, mneut(i), displayDrBarPars().me(j, 3), q)) * 0.5;

      sigmaNeutralino = sigmaNeutralino + neutralinoContribution(i, j);
    }
  }

  sigmaNeutralino = sigmaNeutralino / (16.0 * sqr(PI));
  return sigmaNeutralino;
}

double NmssmSoftsusy::calcRunningMtau() const {
  /// MSbar value
  double mTauSMMZ = displayDataSet().displayMass(mTau);
  double mTauPole = MTAU;
  /// conversion to DRbar
  mTauSMMZ = mTauSMMZ * Softsusy<SoftParsNmssm>::calcRunMtauDrBarConv();
  /// Chargino contribution  
  double sigmaChargino = Softsusy<SoftParsNmssm>::calcRunMtauCharginos(mTauSMMZ);
  /// Higgs
  double sigmaHiggs = calcRunMtauHiggs();
  /// Neutralinos
  double  sigmaNeutralino = calcRunMtauNeutralinos(mTauSMMZ);

  return mTauSMMZ * (1.0 + sigmaNeutralino + sigmaChargino + sigmaHiggs);
}

//PA: returns the mixing of Hu into between h1
double NmssmSoftsusy::h1s2Mix(){
   return displayDrBarPars().mixh0(1,2);
}


double NmssmSoftsusy::predTanb(double MuEff) const  {
  double lam = displayLambda();
  double kap = displayKappa();
  double m3sq = displayM3Squared();
  double mupr = displayMupr();
  double s = displaySvev();
  double xiF = displayXiF();
  double al = displayTrialambda();
  double m3sqeff = m3sq  + lam * (mupr * s / root2 + xiF)
    + al * s / root2  +   0.5 * lam * kap * sqr(s); 
   if (MuEff < -6.e66) MuEff = displaySusyMu() + lam * s / root2;
  double sin2t = 2.0 * m3sqeff / 
    (displayMh1Squared() - displayTadpole1Ms() + 
     displayMh2Squared() - displayTadpole2Ms() + 2.0 *
     sqr(MuEff)); 
  
  /// Note: we want to take inverse sine so that fundamental domain is greater
  /// than pi/4. sin(pi - 2 beta)=sin 2 beta should achieve this.
  /// we also use tan (pi/2 - theta) = 1/tan(theta)
  double theta;
  if (fabs(sin2t) < 1.0) theta = asin(sin2t) * 0.5;
  else return 0.0;
  
  return 1.0 / tan(theta);
}

/// PA: Nmssm version modified to include full mueff = mu + lambda * s / root2
// Apply at scale MSusy:
/// Displays PHYSICAL MZ, ie MZ(q) - piZz^T(q)
double NmssmSoftsusy::predMzsq(double & tanb, double mueffOld, double eps) {
  //PA: new Nmssm parameters needed for mu effective
  double lam = displayLambda();
  double s = displaySvev();
  
  if (fabs(displayTadpole1Ms()) < EPSTOL && 
      fabs(displayTadpole2Ms()) < EPSTOL) {
    double sinthDRbar = calcSinthdrbar();
    calcDrBarPars(); 
    double mt = displayDrBarPars().mt;
    doTadpoles(mt, sinthDRbar);
  }
  //PA: this is now mueff = mu + lam * s /root2
  double susyMu = displaySusyMu() + lam * s /  root2;
  tanb = predTanb(susyMu);
  if (mueffOld > -6.e66) susyMu = susyMu / eps - mueffOld * (1. / eps - 1.);

  double pizztMS = sqr(displayMzRun()) - sqr(displayMz()); ///< resums logs
  double MZsq = 2.0 *
    ((displayMh1Squared() - displayTadpole1Ms() - 
      (displayMh2Squared() - displayTadpole2Ms()) *
      sqr(tanb)) / (sqr(tanb) - 1.0) - sqr(susyMu)) - 
    pizztMS;

  return MZsq;
}


void NmssmSoftsusy::itLowsoft
(int maxTries, double & mx, int sgnMu, double tol, double tanb, 
 void (*boundaryCondition)(NmssmSoftsusy &, const DoubleVector &), 
 const DoubleVector & pars, const DoubleVector nmpars, bool gaugeUnification, bool ewsbBCscale) {

  static NmssmSoftsusy old;
  /// LCT: Resize the drBarPars object to account for NMSSM parameters 
  drBarPars o(old.displayDrBarPars());
  o.mh0.setEnd(3);
  o.mA0.setEnd(2);
  o.mixh0.resize(3,3);
  o.mixA0.resize(2,2);
  o.mneut.setEnd(5);

  o.mixNeut.resize(5,5);
  o.mnBpmz.setEnd(5);
  o.nBpmz.resize(5,5);

  old.setDrBarPars(o);

  static double oldMu = 0.;
  static int numTries = 0;
  double mz = displayMz();

  if (numTries != 0 && sqr(displayMu() / mz - 1.0) > TOLERANCE) {
    cout << "WARNING: itLowsoft called at inappropriate";
    cout << " scale:" << displayMu() << endl; 
    cout << "whereas it should be " << mz << endl; 
  }
  
  if (numTries - 1 > maxTries) {/// Iterating too long: bail out
    flagNoConvergence(true);
    if (PRINTOUT) cout << "itLowsoft reached maxtries\n"; 
    numTries = 0; 
    return;
  }

  if (PRINTOUT > 1) cout << displayProblem(); 

  double mtpole, mtrun;
  
  mtpole = displayDataSet().displayPoleMt();
  /// On first iteration, don't bother with finite corrections
  
  numTries = numTries + 1;
 
  //PA: To check convergence of new parameters.
  double lamold = displayLambda(), kapold = displayKappa(), Svevold = displaySvev(), XiFold = displayXiF(), mu_sold = displayMupr();
 
  //PA: reset new low energy inputs of general nmssm at mz.
  if(Z3){
    setXiF(0.0);                         
    setMupr(0.0);
  }
  else  {
  setKappa(nmpars(2));
  setSvev(nmpars(3));
  setXiF(nmpars(4));                         
  setMupr(nmpars(5));
  }

  try {
    sparticleThresholdCorrections(tanb); 
    
    if (displayProblem().noRhoConvergence) {
      if (PRINTOUT) cout << "No convergence in rhohat\n"; 
    }

    /// precision of running/RGE integration: start off low and increase
    double eps = maximum(exp(double(- numTries) * log(10.0)), tol * 0.01); 
    
    /// first stab at MSUSY: root(mstop1(MZ) mstop2(MZ))
    if (numTries == 1) setMsusy(calcMs()); 

    // cout << "displayMsusy() = " << displayMsusy() << endl; 

    int err = 0;
    err = runto(displayMsusy(), eps);
    double tbIn; double predictedMzSq = 0.;
    predictedMzSq = predMzsq(tbIn);
    setPredMzSq(predictedMzSq); 
    if(!GUTlambda) setLambda(nmpars(1));
                      
    if (!ewsbBCscale) err = runto(mx, eps);

    /// Guard against the top Yukawa fixed point
    if (displayYukawaElement(YU, 3, 3) > 3.0 
	|| displayYukawaElement(YD, 3, 3) > 3.0 
	|| displayYukawaElement(YE, 3, 3) > 3.0) {
      setYukawaElement(YU, 3, 3, minimum(3.0, displayYukawaElement(YU, 3, 3)));
      setYukawaElement(YD, 3, 3, minimum(3.0, displayYukawaElement(YD, 3, 3)));
      setYukawaElement(YE, 3, 3, minimum(3.0, displayYukawaElement(YE, 3, 3)));
      flagIrqfp(true); 
    }
    
    
    if (err) {
      /// problem with running: bail out 
      flagNonperturbative(true);
      if (PRINTOUT) 
	cout << "itLowsoft gone non-perturbative approaching mgut\n"; 
      if (PRINTOUT > 1) printObj();
      numTries = 0; 
      return;
    }
  
    if (gaugeUnification) {
      nmsBrevity (dummy);
      NmssmSusy a(this -> NmssmSusy::beta(dummy));
      
      /// Equal gauge couplings: let them and their derivatives set the boundary
      /// condition scale -- linear approximation
      mx = mx * exp((displayGaugeCoupling(2) - displayGaugeCoupling(1))
		    / (a.displayGaugeCoupling(1) - a.displayGaugeCoupling(2)));

      /// if mx is too high/low, will likely get non-perturbative problems
      if (mx < 1.0e4) {
	mx = 1.0e4;
	if (PRINTOUT > 2) cout << " mx too low ";
	flagMgutOutOfBounds(true);
      }
      if (mx > 5.0e17) {
	if (PRINTOUT > 2) cout << " mx =" << mx <<" too high ";
	mx = 5.0e17;
	flagMgutOutOfBounds(true);
      }
    }
    
    boundaryCondition(*this, pars);
     if(GUTlambda) setLambda(nmpars(1));
    if (!ewsbBCscale) err = runto(displayMsusy(), eps);

    calcDrBarPars();

     if (err) {
      // problem with running: bail out 
      flagNonperturbative(true);
      if (PRINTOUT) cout << "itLowsoft gone non-perturbative on way to MZ\n"; 
      if (PRINTOUT > 1) printObj();
      numTries = 0;
      return;
    }

    setMsusy(calcMs());
    if (ewsbBCscale) mx = displayMsusy();
    if (PRINTOUT > 0) cout << " mgut=" << mx << flush;
    
    mtrun = displayDrBarPars().mt;
    if (numTries < 11) {
      rewsb(sgnMu, mtrun);
    }
    else { ///< After 11 tries, we start averaging old/new mu values
      double epsi = 0.5;
      if (numTries > 20) epsi = 0.2;
      if (numTries > 30) epsi = 0.1;
      rewsb(sgnMu, mtrun, oldMu, epsi);    
      }     

    oldMu = displaySusyMu();
       
    double fD = sumTol(*this, old, numTries);    
    setFracDiff(fD);
    
    if (numTries !=0 && fD < tol) {///< Accuracy achieved: bail out
      numTries = 0; ///< Reset the number of iterations for the next time
      if (PRINTOUT > 0) cout << " sT=" << fD << " " << flush; 
      if (displayProblem().test() && PRINTOUT > 0) 
	cout << " ***problem point***: " << displayProblem() << ".";
      
      return; 
    }
    

    tachyonType nil(none);
    flagAllProblems(false, nil);

    /// Old iteration is 'old': these are the parameters by which convergence is
    /// measured. 
    old.setDrBarPars(displayDrBarPars());
    /// If a print out is desired, print respectively, the difference with the
    /// last iteration (sum tol or sT), the mu parameter and m3^2 from EWSB, and
    /// the predicted MW and MZ boson masses
    if (PRINTOUT > 0) 
      cout << "\n" << numTries << ". sT=" << displayFracDiff()  
	   << " mu=" << displaySusyMu() <<  " m3sq=" << displayM3Squared() 
	   << " MWp=" << displayMw() << " Mzp=" << sqrt(displayPredMzSq())
	   << flush;
   
    if (displayProblem().noMuConvergence) {
      if (PRINTOUT) 
	cout << "itLowsoft doesn't break EWSB\n"; 
      if (PRINTOUT > 1) printObj();
    }
    
    err = runto(mz, eps);
    if (err) {
      /// problem with running: bail out 
      flagNonperturbative(true);
      if (PRINTOUT) cout << "itLowsoft gone non-perturbative on way to MZ\n"; 
      if (PRINTOUT > 1) printObj();
      ///    old = NmssmSoftsusy();
      numTries = 0;
      return;
    }
    
    itLowsoft(maxTries, mx, sgnMu, tol, tanb, boundaryCondition, pars, nmpars,
	      gaugeUnification, ewsbBCscale);
  }
  catch(const char *a) {
    numTries = 0;
    throw a;
  }
  catch(const string &a) {
    numTries = 0;
    throw a;
  }
}

/// Difference between two SOFTSUSY objects in and out: EWSB terms only
/// LCT: Placed here for now.
double sumTol(const NmssmSoftsusy & in, const NmssmSoftsusy & out, int numTries) {

  const drBarPars& inforLoops = in.displayDrBarPars();
  const drBarPars& outforLoops = out.displayDrBarPars();

  const int rank = in.displayDrBarPars().mneut.displayEnd();

  // cout << "rank = " << rank << endl;

  DoubleVector sT(34);
  int k = 1;

  double sTin  = fabs(inforLoops.mh0(1)); double sTout = fabs(outforLoops.mh0(1));

  // cout << "sTin = " << sTin << " .  sTout = " << sTout << endl;

  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mA0(1)); sTout = fabs(outforLoops.mA0(1));
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mh0(2)); sTout = fabs(outforLoops.mh0(2));
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mHpm); sTout = fabs(outforLoops.mHpm);
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;

  // cout << "sumTol af Higgs sT = " << sT << endl;

  int i; for (i=1; i<=3; i++) {
    sTin  = fabs(inforLoops.msnu(i));
    sTout = fabs(outforLoops.msnu(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k++;
  }

  // cout << "sumTol af snu sT = " << sT << endl;

  for (i=1; i<=2; i++) {
    sTin = fabs(inforLoops.mch(i));
    sTout = fabs(outforLoops.mch(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k++;
  }

  // cout << "sumTol af chargino sT = " << sT << endl;

  for (i=1; i<=rank; i++) {
    sTin = fabs(inforLoops.mneut(i));
    sTout = fabs(outforLoops.mneut(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k++;
  }

  // cout << "sumTol af neutralino sT = " << sT << endl;

  sTin = fabs(inforLoops.mGluino);
  sTout = fabs(outforLoops.mGluino);
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
  k++;

  // cout << "sumTol af gluino sT = " << sT << endl;

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

  // cout << "sumTol af sfermions sT = " << sT << endl;

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


/// Provides the first guess at a SUSY object at mt, inputting tanb and oneset
/// (should be at MZ) - it's very crude, doesn't take radiative corrections
/// into account etc. 
NmssmSusy NmssmSoftsusy::guessAtSusyMt(double tanb, DoubleVector nmpars, const QedQcd & oneset) {
   
   //PA: Most of the work is already done by the MSSM
   NmssmSusy t(Softsusy<SoftParsNmssm>::guessAtSusyMt(tanb, oneset));
   //PA: now we just add our new nmssm parameters. 
   //Only lambda directly affects the running of the MSSM Yukawas 
   // at one loop and only kappa additionally at two loop.
   t.setLambda(nmpars(1));
   t.setKappa(nmpars(2));
   t.setSvev(nmpars(3));
   t.setXiF(nmpars(4));       
   t.setMupr(nmpars(5));
   return t;
}



double NmssmSoftsusy::lowOrg
(void (*boundaryCondition)(NmssmSoftsusy &, const DoubleVector &),
 double mxGuess,  const DoubleVector & pars, const DoubleVector nmpars, 
 int sgnMu, double tanb, const QedQcd & oneset, bool gaugeUnification, 
 bool ewsbBCscale) {

  double mx = 0.0;

  try {
     
    double muFirst = displaySusyMu(); /// Remember initial value

    const static NmssmSoftsusy empty;
    bool setTbAtMXflag = displaySetTbAtMX(); 

    setSoftsusy(empty); /// Always starts from an empty object
    /// These are things that are re-written by the new initialisation
    setSetTbAtMX(setTbAtMXflag);
    setData(oneset); 
    setMw(MW); 

    double mz = displayMz();

    if (mxGuess > 0.0) 
      mx = mxGuess; 
    else {
      string ii("Trying to use negative mx in NmssmSoftsusy::lowOrg.\n");
      ii = ii + "Now illegal! Use positive mx for first guess of mx.\n";
      throw ii;
    }
    
    if (oneset.displayMu() != mz) {
      cout << "WARNING: lowOrg in softsusy.cpp called with oneset at scale\n" 
	   << oneset.displayMu() << "\ninstead of " << mz << endl;
    }
    
    /// LCT: Changed maxtries to match softsusy.cpp 8/8/13
    int maxtries = 100;//int(-log(TOLERANCE) / log(10.0) * 10);
    double tol = TOLERANCE;
    
    NmssmSusy t(guessAtSusyMt(tanb, nmpars, oneset));
    t.setLoops(2); /// 2 loops should protect against ht Landau pole 
    t.runto(mx); 
    setSusy(t);

    /// Initial guess: B=0, mu=1st parameter, need better guesses
    boundaryCondition(*this, pars);
    if(GUTlambda) setLambda(nmpars(1));
    if ((sgnMu == 1 || sgnMu == -1) && !ewsbBCscale) {
      /// LCT: Changed sets to match softsusy.cpp 8/8/13
      if(Z3){
	setSusyMu(0.0);
	setM3Squared(0.0);
      }
      else {
	setSusyMu(sgnMu * MZ);
	setM3Squared(1.0e6);
      }
     
    }
    else {
      setSusyMu(muFirst);
      setM3Squared(muFirst); 
    }
   
    run(mx, mz);

    if (sgnMu == 1 || sgnMu == -1) rewsbTreeLevel(sgnMu); 
  
    physical(0);

    setThresholds(3); setLoops(2);
    
    //PA: itLowsoft to be added along with the rest of lowOrg
    itLowsoft(maxtries, mx, sgnMu, tol, tanb, boundaryCondition, pars, 
              nmpars, gaugeUnification, ewsbBCscale);
    
      // if (displayProblem().nonperturbative 
      //   || displayProblem().higgsUfb || displayProblem().tachyon 
      //   || displayProblem().noRhoConvergence)
      // return mx;
     
    
    runto(maximum(displayMsusy(), MZ));
    if (ewsbBCscale) boundaryCondition(*this, pars); 

    physical(3);

    runto(mz);

    if (PRINTOUT > 1) cout << " end of iteration" << endl;
    
  }
  catch(const char *a) {
    ostringstream ii;
    ii << "SOFTSUSY problem: " << a << " pars=" << pars << " tanb=" << tanb 
       << " oneset=" << oneset << endl;
    flagProblemThrown(true);
    throw(ii.str());
  }
  catch(const string & a) {
    ostringstream ii;
    ii << "SOFTSUSY problem: " << a << " pars=" << pars << " tanb=" << tanb 
	 << " oneset=" << oneset << endl;
    flagProblemThrown(true);
    throw ii.str();
  }
  catch(...) {
    ostringstream ii;
    ii << "SOFTSUSY problem: " << endl;
    ii << "pars=" << pars << " tanb=" << tanb
       << " oneset=" << oneset << endl;
    flagProblemThrown(true);
    throw ii.str();
  }

  return mx;
}

/// User supplied routine. Inputs m at the unification scale, and uses
/// inputParameters vector to output m with high energy soft boundary
/// conditions. 
void NmssmMsugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);

  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be a0
  ///  if (m0 < 0.0) m.flagTachyon(true); Deleted on request from A Pukhov
  m.standardSugra(m0, m12, a0);
    
  return;
}
//PA: msugra bcs in the mssm limit of the general nmssm
void MssmMsugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);

  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be a0
  ///  if (m0 < 0.0) m.flagTachyon(true); Deleted on request from A Pukhov
  m.standardsemiSugra(m0, m12, a0, 0.0, 1e-15);
  m.setMspSquared(1e6);  
  return;
}

//PA: semi-msugra bcs for the nmssm
void SemiMsugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);
  //  double mS = inputParameters.display(4);
  double Al = inputParameters.display(4);
  double Ak = inputParameters.display(5);
  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be a0
  ///  if (m0 < 0.0) m.flagTachyon(true); Deleted on request from A Pukhov
  m.standardsemiSugra(m0, m12, a0, Al, Ak);
    
  return;
}

   void generalNmssmBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  NmssmSusy s; SoftParsNmssm r;
  s = m.displaySusy();
  r.set(inputParameters);
  if(Z3==false){
    double m3sq = m.displayM3Squared();
    double XiS = m.displayXiS();
    r.setM3Squared(m3sq);
    r.setXiS(XiS);
  }
  else{
    double mSsq = m.displayMsSquared();
    r.setMsSquared(mSsq);
  }
  
 
  m.setSoftPars(r);
  m.setSusy(s);

  return;
}
#endif
