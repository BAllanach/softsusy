/** \file decays.cpp
   Project:     SOFTSUSY 
   Author:      Tom Cridge, Ben Allanach
   Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   Webpage:     http://softsusy.hepforge.org/

*/

#include "decays.h"

using namespace std;

/// decay global variable declarations
double m1, m2, m3, m4, mq, m5, m6, 
  m7,  m8,  MZboson,  MWboson,  mh,  mH, 
  mA,  mphi,  g1,  g2,  alphamix,  betavac;
int neutralinoj, neutralinoi, AorhorH;
DoubleMatrix NeutMIX(4, 4);

int calculateDecays(ostream & fout, MssmSoftsusy * r,
		    vector<Particle> & decayTable, 
		    const NmssmSoftsusy & nmssm, bool nmssmIsIt) { 
  /// We work on the principle that any problem flags should already be set,
  /// and calculateDecays should not change them
  sProblem rProb = r->displayProblem();
  
  /// Initialise global decay variables
  m1 = 0.; m2 = 0.; m3 = 0.; m4 = 0.; mq = 0.; m5 = 0.; m6 = 0.; 
  m7 = 0.;  m8 = 0.;  MZboson = 0.;  MWboson = 0.;  mh = 0.;  mH = 0.; 
  mA = 0.;  mphi = 0.;  g1 = 0.;  g2 = 0.;  alphamix = 0.;  betavac = 0.,
  neutralinoj = 0;  neutralinoi = 0;  AorhorH = 0;
  for (int i=1; i<=4; i++)
    for (int j=1; j<=4; j++)
      NeutMIX(i, j) = 0.;
  int errorflag = 0;
  
  /// If there is a serious problem with the point, return an error code and
  /// warning
    if (
      (!nmssmIsIt && r->displayProblem().testSeriousProblem()) ||
      (nmssmIsIt && nmssm.displayProblem().testSeriousProblem())
      ) {
    fout << "# Not calculating decays: serious problem with point\n";
    return -1;
  }

    bool flaggluino = 1, flagsupL = 1, flagsupR = 1, flagsdownL = 1, flagsdownR = 1, flagscharmL = 1, flagscharmR = 1, flagsstrangeL = 1, flagsstrangeR = 1, flagstop1 = 1, flagstop2 = 1, flagsbottom1 = 1, flagsbottom2 = 1, flagselectronL = 1, flagselectronR = 1, flagsmuonL = 1, flagsmuonR = 1, flagstau1 = 1, flagstau2 = 1, flagsnueL = 1, flagsnumuL = 1, flagsnutauL = 1, flagneut1 = 1, flagneut2 = 1, flagneut3 = 1, flagneut4 = 1, flagneut5 = 1, flagchar1 = 1, flagchar2 = 1, flagh1 = 1, flagH2 = 1, flagH3 = 1, flagA1 = 1, flagA2 = 1, flagHpm = 1; ///< Flags to turn off decays, default 1 = on, 0 = off

  bool QCDcorr = true; ///Turns on QCD corrections to h->gg and h->qq

  ///Switch on or off 1->3 decays
  bool onetothree = threeBodyDecays; ///Turns on 1->3 decays, reads this in from input file, default is true

  NmssmSoftsusy nmssmrun = nmssm; ///Do this as nmssm is a const to ensure it can't be changed here (e.g. if run in RPV mode where NMSSM is not included don't want this to change NMSSM), copy to nmssmrun to allow me to run it to MSusy, 1000GeV or whatever scale for parameter extraction

  DoubleMatrix S(3,3);
  DoubleMatrix CPEMix(3,3); ///CP even higgses mixing matrix which I will transform to match that in NMSSMTools
   ///initialise CPEMix
   for (int i = 1; i <= 3; i++) {
     for (int j = 1; j <= 3; j++) {
      CPEMix(i,j) = 0;
     }
   }
   double thetaA = 0; double kappa = 0; double lam = 0; double Alambda = 0; double Akappa = 0; double svev = 0; double gs = 0; double gp =0; double g=0; double alphas = 0; double tanthetaW = 0; double greekmu = 0;  double mGluino = 0; DoubleVector mch(2); double mHpm; DoubleVector msnu(3); double thetaL=0; double thetaR=0; DoubleMatrix mu(2,3); DoubleMatrix md(2,3); DoubleMatrix me(2,3);  double mwSoftSusy=0; double runmz=0; double polemw=0; double thetaL2 = 0; double thetaR2 = 0; double mtPole = 0, mbPole = 0;
   DoubleVector mneut(5); DoubleVector mh0(3); DoubleVector mA0(2); DoubleMatrix mixNeut(5,5); DoubleMatrix mixh0(3,3);
   DoubleMatrix Ptemp(3,3);

  DoubleMatrix CPOMix(3,3); ///is pseudoscalar mixing matrix as in P in NMSSMTools (dropped first row as that's goldstones hence 2x3 - just made 3x3 to stop issue with DoubleMatrix, third row is all 0s)
  DoubleMatrix CPOMix2(2,2); ///is P2 in NMSSMTools
   double tanbeta = 0;
   double beta = 0;
   double alpha = 0;
   double runmw = 0;
   double polemz = 0;

   double thetat = 0 , thetab = 0, thetatau = 0;

   double runmt = 0, runmb = 0, runmtau = 0, runmc =0, runms = 0, runmd = 0, runmu = 0, runmel = 0, runmmu = 0, Au = 0, Ad = 0, Ac = 0, As =0, At =0, Ab =0, Atau = 0, Ae =0, Amu =0, mueff=0;

   double alphasAtMA = 0, alphaAtMA = 0, mbAtMA = 0, mtAtMA = 0, mcAtMA = 0, msAtMA = 0, alphasAtMH = 0, alphaAtMH = 0, mbAtMH = 0, mtAtMH = 0, mcAtMH = 0, msAtMH = 0, alphasAtMh = 0, alphaAtMh = 0, mbAtMh = 0, mtAtMh = 0, mcAtMh = 0, msAtMh = 0, alphasAtMA2 = 0, alphaAtMA2 = 0, mbAtMA2 = 0, mtAtMA2 = 0, mcAtMA2 = 0, alphasAtMH3 = 0, alphaAtMH3 = 0, mbAtMH3 = 0, mtAtMH3 = 0, mcAtMH3 = 0; //For running couplings and masses for higgs 1-loop decays
   // msAtMH3 = 0; 
   double g3atmh0 = 0, g3atmH0 = 0, g3atmA0 = 0;
   double mt = 0, mb = 0, mc = 0, ms = 0, mup = 0, mdo = 0, mel = 0, mmu = 0, mtau = 0; ///Quark pole masses for general use

 if(nmssmIsIt == true) {
   
   if(onetothree == true) {
     fout << "# No 1to3 decays included in NMSSM - therefore onetothree set to false" << endl;
   }
   onetothree = false;

   // fout << "MSUSY = " << nmssm.displayMsusy() << endl;

   nmssmrun.runto(nmssm.displayMsusy()); ///Run to scale Msusy for parameter extraction

   nmssmrun.calcDrBarPars(); ///Must redo calcDrBarPars at each scale

   S = nmssmrun.displayDrBarPars().mixh0; /// CP even higgs 3x3 mixing matrix as from softsusy

   CPEMix(1,1) = S(1,2); CPEMix(2,1) = S(2,2); CPEMix(1,2) = S(1,1); CPEMix(2,2) = S(2,1); CPEMix(3,1) = S(3,2); CPEMix(3,2) = S(3,1); CPEMix(1,3) = S(1,3); CPEMix(2,3) = S(2,3); CPEMix(3,3) = S(3,3); ///Transform to match conventions in NMSSMTools


   tanbeta = nmssmrun.displayTanb();
   beta = atan(tanbeta);
   thetaA = nmssmrun.displayDrBarPars().thetaA0; ///CP odd higgses mixing angle

   // fout << "CPEMix = " << CPEMix << std::endl;
   // fout << "beta = " << beta << std::endl;

   Ptemp(1,1) = -cos(beta), Ptemp(1,2) = sin(beta), Ptemp(1,3) = 0;
   Ptemp(2,1) = sin(beta)*cos(thetaA); Ptemp(2,2) = cos(beta)*cos(thetaA); Ptemp(2,3) = sin(thetaA);
   Ptemp(3,1) = sin(beta)*sin(thetaA); Ptemp(3,2) = cos(beta)*sin(thetaA); Ptemp(3,3) = -cos(thetaA);
  
   ///CPOMix is the mixing matrix used in the partial width formulae, not Ptemp or CPOMix2
   CPOMix(1,1) = Ptemp(2,2); CPOMix(1,2) = Ptemp(2,1); CPOMix(1,3) = Ptemp(2,3);
   CPOMix(2,1) = Ptemp(3,2); CPOMix(2,2) = Ptemp(3,1); CPOMix(2,3) = Ptemp(3,3);
   CPOMix(3,1) = 0; CPOMix(3,2) = 0; CPOMix(3,3) = 0; ///0s as never used as CPOMix is actually 2x3 as dropped goldstone row
 
   CPOMix2(1,1) = CPOMix(1,1)/sin(beta);
   CPOMix2(2,1) = CPOMix(2,1)/sin(beta);
   CPOMix2(1,2) = CPOMix(1,3);
   CPOMix2(2,2) = CPOMix(2,3);
   
   // fout << "S = " << S << " CPEMix = " << CPEMix << " CPOMix = " << CPOMix << endl;

   ///Additional NMSSM parameters
   kappa = nmssmrun.displayKappa();
   lam = nmssmrun.displayLambda();
   Alambda = nmssmrun.displaySoftAlambda();
   Akappa = nmssmrun.displaySoftAkappa();
   svev = nmssmrun.displaySvev();

   gs = nmssmrun.displayGaugeCoupling(3);
   gp= nmssmrun.displayGaugeCoupling(1)*pow(0.6,0.5);
   g = nmssmrun.displayGaugeCoupling(2);
   alphas = pow(gs,2)/(4*PI);
   tanthetaW = gp/g;
   //   alphaEm = pow(e,2)/(4*PI);

   ///Set masses from nmssmrun nmssmsoftsusy object passed from softpoint
   mGluino = nmssmrun.displayPhys().mGluino; mneut = nmssmrun.displayPhys().mneut; mch = nmssmrun.displayPhys().mch; mh0 = nmssmrun.displayPhys().mh0; mA0 = nmssmrun.displayPhys().mA0; mHpm = nmssmrun.displayPhys().mHpm; msnu = nmssmrun.displayPhys().msnu; mixNeut = nmssmrun.displayPhys().mixNeut.transpose();
   thetaL = nmssmrun.displayPhys().thetaL; thetaR = nmssmrun.displayPhys().thetaR; mu = nmssmrun.displayPhys().mu; md = nmssmrun.displayPhys().md; me = nmssmrun.displayPhys().me; DoubleMatrix mixh0(nmssmrun.displayPhys().mixh0); mwSoftSusy = nmssmrun.displayMwRun(); runmz = nmssmrun.displayMzRun(); polemw = nmssmrun.displayMw(); polemz = nmssmrun.displayMz();
 
   ///Set mixing angles in my conventions, note shift in thetaL2 and thetaR2, which are the angles used in the PW formulae
   tanbeta = nmssmrun.displayTanb();
   beta = atan(tanbeta);
   alpha =nmssmrun.displayDrBarPars().thetaH;
   thetaL2 = -thetaL + PI/2;
   thetaR2 = -thetaR + PI/2;
   
   runmw = mwSoftSusy; ///"running" W mass, i.e. that used for couplings

  ///Mixing angles taken from softsusy depend on the mass ordering, I have taken the Spheno method whereby the sfermions are always mass ordered, this isn't necessarily true of softsusy so if the mass ordering is reverse I must change the mixing angle and reorder the masses.
  if (mu(1,3) <= mu(2,3)) {
    // thetat = nmssmrun.displayPhys().thetat; 
    thetat = nmssmrun.displayDrBarPars().thetat;
    }
  else if (mu(1,3) > mu(2,3)) {
    // thetat = acos(-sin(nmssmrun.displayPhys().thetat));
    thetat = acos(-sin(nmssmrun.displayDrBarPars().thetat));
    DoubleMatrix mutemp (2,3);
    mutemp(1,1) = mu(1,1); mutemp(1,2) = mu(1,2); mutemp(2,1) = mu(2,1); mutemp(2,2) = mu(2,2); mutemp(1,3) = mu(2,3); mutemp(2,3) = mu(1,3);
    mu = mutemp;
  }
  if (md(1,3) <= md(2,3)) {
    // thetab = nmssmrun.displayPhys().thetab; 
    thetab = nmssmrun.displayDrBarPars().thetab; 
    }
  else if (md(1,3) > md(2,3)) {
    // thetab = acos(-sin(nmssmrun.displayPhys().thetab));
    thetab = acos(-sin(nmssmrun.displayDrBarPars().thetab));
    DoubleMatrix mdtemp (2,3);
    mdtemp(1,1) = md(1,1); mdtemp(1,2) = md(1,2); mdtemp(2,1) = md(2,1); mdtemp(2,2) = md(2,2); mdtemp(1,3) = md(2,3); mdtemp(2,3) = md(1,3);
    md = mdtemp;
  }
  if (me(1,3) <= me(2,3)) {
    // thetatau = nmssmrun.displayPhys().thetatau; 
    thetatau = nmssmrun.displayDrBarPars().thetatau; 
    }
  else if (me(1,3) > me(2,3)) {
    // thetatau = acos(-sin(nmssmrun.displayPhys().thetatau));
    thetatau = acos(-sin(nmssmrun.displayDrBarPars().thetatau));
    DoubleMatrix metemp (2,3);
    metemp(1,1) = me(1,1); metemp(1,2) = me(1,2); metemp(2,1) = me(2,1); metemp(2,2) = me(2,2); metemp(1,3) = me(2,3); metemp(2,3) = me(1,3);
    me = metemp;
  }

  ///Trilinear couplings
  At = nmssmrun.displaySoftA(UA, 3, 3), Ab = nmssmrun.displaySoftA(DA, 3, 3), Atau = nmssmrun.displaySoftA(EA, 3, 3), greekmu = nmssmrun.displaySusyMu();
  Au = nmssmrun.displaySoftA(UA, 1, 1), Ad = nmssmrun.displaySoftA(DA, 1, 1), Ac = nmssmrun.displaySoftA(UA, 2, 2), As = nmssmrun.displaySoftA(DA, 2, 2), Ae = nmssmrun.displaySoftA(EA, 1, 1), Amu = nmssmrun.displaySoftA(EA, 2, 2);

  ///Effective mu parameter of NMSSM - note resolves mu problem of MSSM
  mueff = lam*svev/pow(2,0.5);

  ///Run to mZ for fermion mass extraction
  nmssmrun.runto(nmssmrun.displayMz());

  ///Pole masses for quarks necessary for scheme used for h -> qq QCD corrections
  mtPole = nmssmrun.displayDataSet().displayPoleMt();
  mbPole = nmssmrun.displayDataSet().displayPoleMb();
  
  ///Masses used in general formulae for quarks
  mt = mtPole;
  mb = mbPole;
  mc = nmssmrun.displayDataSet().displayMass(mCharm);
  ms = nmssmrun.displayDataSet().displayMass(mStrange);
  mup = nmssmrun.displayDataSet().displayMass(mUp);
  mdo = nmssmrun.displayDataSet().displayMass(mDown);
  mel = nmssmrun.displayDataSet().displayMass(mElectron);
  mmu = nmssmrun.displayDataSet().displayMass(mMuon);
  mtau = nmssmrun.displayDataSet().displayPoleMtau();
  
  nmssmrun.calcDrBarPars(); ///Must redo calcDrBarPars at new scale

  ///"Running" masses used for couplings, e.g. in setting yukawas
  runmb = nmssmrun.displayDrBarPars().mb;
  runmt = nmssmrun.displayDrBarPars().mt;
  runmtau = nmssmrun.displayDrBarPars().mtau;
  runmc = pow(2,0.5)*runmw*sin(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YU)(2,2)/nmssmrun.displayGaugeCoupling(2);
  runms = pow(2,0.5)*runmw*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YD)(2,2)/nmssmrun.displayGaugeCoupling(2);
  runmd = pow(2,0.5)*runmw*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YD)(1,1)/nmssmrun.displayGaugeCoupling(2);
  runmu = pow(2,0.5)*runmw*sin(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YU)(1,1)/nmssmrun.displayGaugeCoupling(2);
  runmel = pow(2,0.5)*runmw*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YE)(1,1)/nmssmrun.displayGaugeCoupling(2);
  runmmu = pow(2,0.5)*runmw*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YE)(2,2)/nmssmrun.displayGaugeCoupling(2);

  // ///For NMSSM to match up with NMSSMTools I must transform 5x5 neutralino mixing matrix by swapping columns 3 and 4 - DON'T NEED TO FOR SUSY DECAYS, JUST WOULD NEED TO FOR HIGGS DECAYS BUT ACCOUNTED FOR THIS IN THEIR CODES!


  ///Instead of using Softsusy to run gs to different scales (time inefficient) use the one-loop renormalisation group equations for gs directly to run to different scales, this is accurate enough for our purposes here - NOT WHAT IS USED IN CODE, FULL RUNNING OF SOFTSUSY IS USED, just left for potential use:
 double alphasrun (double mu, double mu0, double alphasmu0); ///See later function definition
 double alpharun (double mu, double mu0, double alphamu0); ///See later function definition
 double alphasatmh = alphasrun(mh0(1), polemz, ALPHASMZ);
 double alphasatmH = alphasrun(mh0(2), polemz, ALPHASMZ);
 double alphasatmA = alphasrun(mA0(1), polemz, ALPHASMZ);
 g3atmh0 = pow(4*PI*alphasatmh,0.5);
 g3atmH0 = pow(4*PI*alphasatmH,0.5);
 g3atmA0 = pow(4*PI*alphasatmA,0.5);


 ///Or run fully to the mass of each decaying Higgs: - WHAT IS ACTUALLY USED, these are used for masses of quarks and for gauge couplings for the Higgs loop decays, here the scale at which you choose to set your masses can cause significant differences in the PWs determined:
 nmssmrun.runto(nmssmrun.displayPhys().mA0(2));
 alphasAtMA2 =  pow(nmssmrun.displayGaugeCoupling(3),2)/(4*PI);
 alphaAtMA2 = pow(nmssmrun.displayGaugeCoupling(2)*nmssmrun.calcSinthdrbar(),2)/(4*PI);
 mbAtMA2 = pow(2,0.5)*nmssmrun.displayMwRun()*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YD)(3,3)/nmssmrun.displayGaugeCoupling(2); 
 mtAtMA2 = pow(2,0.5)*nmssmrun.displayMwRun()*sin(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YU)(3,3)/nmssmrun.displayGaugeCoupling(2);
 mcAtMA2 = pow(2,0.5)*nmssmrun.displayMwRun()*sin(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YU)(2,2)/nmssmrun.displayGaugeCoupling(2);

 nmssmrun.runto(nmssmrun.displayPhys().mh0(3));
 alphasAtMH3 =  pow(nmssmrun.displayGaugeCoupling(3),2)/(4*PI);
 alphaAtMH3 = pow(nmssmrun.displayGaugeCoupling(2)*nmssmrun.calcSinthdrbar(),2)/(4*PI);
 mbAtMH3 = pow(2,0.5)*nmssmrun.displayMwRun()*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YD)(3,3)/nmssmrun.displayGaugeCoupling(2);
 mtAtMH3 = pow(2,0.5)*nmssmrun.displayMwRun()*sin(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YU)(3,3)/nmssmrun.displayGaugeCoupling(2);
 mcAtMH3 = pow(2,0.5)*nmssmrun.displayMwRun()*sin(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YU)(2,2)/nmssmrun.displayGaugeCoupling(2);
 // msAtMH3 = pow(2,0.5)*nmssmrun.displayMwRun()*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YD)(2,2)/nmssmrun.displayGaugeCoupling(2);

 nmssmrun.runto(nmssmrun.displayPhys().mA0(1));
 alphasAtMA =  pow(nmssmrun.displayGaugeCoupling(3),2)/(4*PI);
 alphaAtMA = pow(nmssmrun.displayGaugeCoupling(2)*nmssmrun.calcSinthdrbar(),2)/(4*PI);
 mbAtMA = pow(2,0.5)*nmssmrun.displayMwRun()*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YD)(3,3)/nmssmrun.displayGaugeCoupling(2);
 mtAtMA = pow(2,0.5)*nmssmrun.displayMwRun()*sin(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YU)(3,3)/nmssmrun.displayGaugeCoupling(2);
 mcAtMA = pow(2,0.5)*nmssmrun.displayMwRun()*sin(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YU)(2,2)/nmssmrun.displayGaugeCoupling(2);
 msAtMA = pow(2,0.5)*nmssmrun.displayMwRun()*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YD)(2,2)/nmssmrun.displayGaugeCoupling(2);

 nmssmrun.runto(nmssmrun.displayPhys().mh0(2));
 alphasAtMH = pow(nmssmrun.displayGaugeCoupling(3),2)/(4*PI);
 alphaAtMH = pow(nmssmrun.displayGaugeCoupling(2)*nmssmrun.calcSinthdrbar(),2)/(4*PI);
 mbAtMH = pow(2,0.5)*nmssmrun.displayMwRun()*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YD)(3,3)/nmssmrun.displayGaugeCoupling(2);
 mtAtMH = pow(2,0.5)*nmssmrun.displayMwRun()*sin(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YU)(3,3)/nmssmrun.displayGaugeCoupling(2);
 mcAtMH = pow(2,0.5)*nmssmrun.displayMwRun()*sin(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YU)(2,2)/nmssmrun.displayGaugeCoupling(2);
 msAtMH = pow(2,0.5)*nmssmrun.displayMwRun()*cos(atan(nmssmrun.displayTanb()))*nmssmrun.displayYukawaMatrix(YD)(2,2)/nmssmrun.displayGaugeCoupling(2);
 // double gpAtMH = nmssmrun.displayGaugeCoupling(1)*pow(0.6,0.5);
 // double gAtMH = nmssmrun.displayGaugeCoupling(2);

 ///If yukawa matrix elements are 0 (can occur if only include mixing in third generation in softsusy spectrum generator) then use a non-zero mass - avoids issues of inf in foftau for loop decays!
 if (mbAtMH3 < 1.0e-66) {mbAtMH3 = mb;}
 if (mtAtMH3 < 1.0e-66) {mtAtMH3 = mt;}
 if (mcAtMH3 < 1.0e-66) {mcAtMH3 = mc;}
 // if (msAtMH3 < 1.0e-66) {msAtMH3 = ms;}
 if (mbAtMH < 1.0e-66) {mbAtMH = mb;}
 if (mtAtMH < 1.0e-66) {mtAtMH = mt;}
 if (mcAtMH < 1.0e-66) {mcAtMH = mc;}
 if (msAtMH < 1.0e-66) {msAtMH = ms;}
 if (mbAtMA2 < 1.0e-66) {mbAtMA2 = mb;}
 if (mtAtMA2 < 1.0e-66) {mtAtMA2 = mt;}
 if (mcAtMA2 < 1.0e-66) {mcAtMA2 = mc;}
 if (mbAtMA < 1.0e-66) {mbAtMA = mb;}
 if (mtAtMA < 1.0e-66) {mtAtMA = mt;}
 if (mcAtMA < 1.0e-66) {mcAtMA = mc;}
 if (msAtMA < 1.0e-66) {msAtMA = ms;}

///For mh should not use full susy running as susy not around this scale, instead use SM running as in QeDQcD object defined in lowe.h:
 QedQcd a(nmssmrun.displayDataSet());
 
 a.runto(nmssmrun.displayPhys().mh0(1));
 alphasAtMh = a.displayAlpha(ALPHAS);
 alphaAtMh =  a.displayAlpha(ALPHA);
 mbAtMh = a.displayMass()(6);
 mtAtMh = a.displayMass()(3);
 mcAtMh = a.displayMass()(2);
 msAtMh = a.displayMass()(5);

 g3atmh0 = pow(4*PI*alphasAtMh,0.5);
 g3atmH0 = pow(4*PI*alphasAtMH,0.5);
 g3atmA0 = pow(4*PI*alphasAtMA,0.5);

 }
  else if (nmssmIsIt == false) {  ///MSSM

    r->runto(r->displayMsusy()); ///Run to Msusy scale to extract parameters
    r->calcDrBarPars();
 
    gs = r->displayGaugeCoupling(3);
    gp= r->displayGaugeCoupling(1)*pow(0.6,0.5);
    g = r->displayGaugeCoupling(2);

    alphas = pow(gs,2)/(4*PI);
    tanthetaW = gp/g;

    ///Set masses
    mGluino = r->displayPhys().mGluino; mneut = r->displayPhys().mneut; mch = r->displayPhys().mch; mh0 = r->displayPhys().mh0; mA0 = r->displayPhys().mA0; mHpm = r->displayPhys().mHpm; msnu = r->displayPhys().msnu; mixNeut = r->displayPhys().mixNeut.transpose();
    thetaL = r->displayPhys().thetaL; thetaR = r->displayPhys().thetaR; mu = r->displayPhys().mu; md = r->displayPhys().md; me = r->displayPhys().me; mixh0 = r->displayPhys().mixh0; mwSoftSusy = r->displayMwRun(); runmz = r->displayMzRun(); polemw = r->displayMw(); polemz = r->displayMz();
    
    ///Set mixing angles - again note shift in thetaL2 and thetaR2 used in PW formulae
    tanbeta = r->displayTanb();
    beta = atan(tanbeta);
    alpha =r->displayDrBarPars().thetaH;
    thetaL2 = -thetaL + PI/2;
    thetaR2 = -thetaR + PI/2;
    
    runmw = mwSoftSusy; ///"Running" W mass, used for couplings
    
    ///Ensure mass ordering of sfermions, requires shift of sfermion mixing angles if mass order initially wrong
    if (mu(1,3) <= mu(2,3)) {
      // thetat = r->displayPhys().thetat; 
      thetat = r->displayDrBarPars().thetat; 
    }
    else if (mu(1,3) > mu(2,3)) {
    // thetat = acos(-sin(r->displayPhys().thetat)); ///Equivalent of adding pi/2 for the range -pi/4 -> pi/4
      thetat = acos(-sin(r->displayDrBarPars().thetat));
      DoubleMatrix mutemp (2,3);
      mutemp(1,1) = mu(1,1); mutemp(1,2) = mu(1,2); mutemp(2,1) = mu(2,1); mutemp(2,2) = mu(2,2); mutemp(1,3) = mu(2,3); mutemp(2,3) = mu(1,3);
      mu = mutemp;
    }
    
    if (md(1,3) <= md(2,3)) {
      // thetab = r->displayPhys().thetab; 
      thetab = r->displayDrBarPars().thetab; 
    }
    else if (md(1,3) > md(2,3)) {
    // thetab = acos(-sin(r->displayPhys().thetab)); ///Equivalent of adding pi/2 for the range -pi/4 -> pi/4
      thetab = acos(-sin(r->displayDrBarPars().thetab));
      DoubleMatrix mdtemp (2,3);
      mdtemp(1,1) = md(1,1); mdtemp(1,2) = md(1,2); mdtemp(2,1) = md(2,1); mdtemp(2,2) = md(2,2); mdtemp(1,3) = md(2,3); mdtemp(2,3) = md(1,3);
      md = mdtemp;
    }
    if (me(1,3) <= me(2,3)) {
      // thetatau = r->displayPhys().thetatau; 
      thetatau = r->displayDrBarPars().thetatau; 
    }
    else if (me(1,3) > me(2,3)) {
      // thetatau = acos(-sin(r->displayPhys().thetatau)); ///Equivalent of adding pi/2 for the range -pi/4 -> pi/4
      thetatau = acos(-sin(r->displayDrBarPars().thetatau));
      DoubleMatrix metemp (2,3);
      metemp(1,1) = me(1,1); metemp(1,2) = me(1,2); metemp(2,1) = me(2,1); metemp(2,2) = me(2,2); metemp(1,3) = me(2,3); metemp(2,3) = me(1,3);
      me = metemp;
    }
    
    ///Trilinear couplings
    At = r->displaySoftA(UA, 3, 3), Ab = r->displaySoftA(DA, 3, 3), Atau = r->displaySoftA(EA, 3, 3), greekmu = r->displaySusyMu();
    Au = r->displaySoftA(UA, 1, 1), Ad = r->displaySoftA(DA, 1, 1), Ac = r->displaySoftA(UA, 2, 2), As = r->displaySoftA(DA, 2, 2), Ae = r->displaySoftA(EA, 1, 1), Amu = r->displaySoftA(EA, 2, 2);
    
    r->runto(r->displayMz()); ///Run to mZ to extract quark masses
    mtPole = r->displayDataSet().displayPoleMt();
    mbPole = r->displayDataSet().displayPoleMb();
    
    ///Masses used in general formulae for quarks
    mt = mtPole;
    mb = mbPole;
    mc = r->displayDataSet().displayMass(mCharm);
    ms = r->displayDataSet().displayMass(mStrange);
    mup = r->displayDataSet().displayMass(mUp);
    mdo = r->displayDataSet().displayMass(mDown);
    mel = r->displayDataSet().displayMass(mElectron);
    mmu = r->displayDataSet().displayMass(mMuon);
    mtau = r->displayDataSet().displayPoleMtau();
    
    r->calcDrBarPars(); ///Need to recalc DrBarPars at each new scale
    
    ///"Running" masses used for couplings, e.g. yukawa couplings
    runmb = r->displayDrBarPars().mb;
    runmt = r->displayDrBarPars().mt;
    runmtau = r->displayDrBarPars().mtau;
    
    runmc = pow(2,0.5)*polemw*sin(atan(r->displayTanb()))*r->displayYukawaMatrix(YU)(2,2)/r->displayGaugeCoupling(2);
    runms = pow(2,0.5)*polemw*cos(atan(r->displayTanb()))*r->displayYukawaMatrix(YD)(2,2)/r->displayGaugeCoupling(2);
    runmd = pow(2,0.5)*polemw*cos(atan(r->displayTanb()))*r->displayYukawaMatrix(YD)(1,1)/r->displayGaugeCoupling(2);
    runmu = pow(2,0.5)*polemw*sin(atan(r->displayTanb()))*r->displayYukawaMatrix(YU)(1,1)/r->displayGaugeCoupling(2);
    runmel = pow(2,0.5)*runmw*cos(atan(r->displayTanb()))*r->displayYukawaMatrix(YE)(1,1)/r->displayGaugeCoupling(2);
    runmmu = pow(2,0.5)*runmw*cos(atan(r->displayTanb()))*r->displayYukawaMatrix(YE)(2,2)/r->displayGaugeCoupling(2);
    
    
    ///Instead of using Softsusy to run gs to different scales (time inefficient) use the one-loop renormalisation group equations for gs directly to run to different scales, this is accurate enough for our purposes here, NOT WHAT IS DONE - full running is done (see a few lines below) - just left for potential use:
    double alphasrun (double mu, double mu0, double alphasmu0); ///See later function definition
    double alpharun (double mu, double mu0, double alphamu0); ///See later function definition
    double alphasatmh = alphasrun(mh0(1), polemz, ALPHASMZ);
    double alphasatmH = alphasrun(mh0(2), polemz, ALPHASMZ);
    double alphasatmA = alphasrun(mA0(1), polemz, ALPHASMZ);
    g3atmh0 = pow(4*PI*alphasatmh,0.5);
    g3atmH0 = pow(4*PI*alphasatmH,0.5);
    g3atmA0 = pow(4*PI*alphasatmA,0.5);
    
    ///Look into hdecay's alphas values via their running function - NOT USED:
    double alphasrunlambdaQCD (double mu, double LAMBDA, double Nf);
    
    ///Or run fully using softsusy runto: - THIS IS WHAT IS ACTUALLY USED NOT 1-LOOP RENORMALISATION GROUP EQUATIONS APPROXIMATE WAY OF A FEW LINES ABOVE
    ///Determine gauge couplings and quark masses at the mass of the decaying Higgses for use in the Higgs loop decays, here the choice of scale for evaluating the mass of the quarks and couplings can have a significant effect on the PWs
    r->runto(r->displayPhys().mA0(1));
    r->calcDrBarPars();
    alphasAtMA = pow(r->displayGaugeCoupling(3),2)/(4*PI);
    alphaAtMA = pow(r->displayGaugeCoupling(2)*r->calcSinthdrbar(),2)/(4*PI);
    
    mbAtMA = pow(2,0.5)*mwSoftSusy*cos(atan(r->displayTanb()))*r->displayYukawaMatrix(YD)(3,3)/r->displayGaugeCoupling(2);
    mtAtMA = pow(2,0.5)*mwSoftSusy*sin(atan(r->displayTanb()))*r->displayYukawaMatrix(YU)(3,3)/r->displayGaugeCoupling(2);
    mcAtMA = pow(2,0.5)*mwSoftSusy*sin(atan(r->displayTanb()))*r->displayYukawaMatrix(YU)(2,2)/r->displayGaugeCoupling(2);
    msAtMA = pow(2,0.5)*mwSoftSusy*cos(atan(r->displayTanb()))*r->displayYukawaMatrix(YD)(2,2)/r->displayGaugeCoupling(2);
    
    r->runto(r->displayPhys().mh0(2));
    r->calcDrBarPars();
    alphasAtMH = pow(r->displayGaugeCoupling(3),2)/(4*PI);
    alphaAtMH = pow(r->displayGaugeCoupling(2)*r->calcSinthdrbar(),2)/(4*PI);
    mbAtMH = pow(2,0.5)*polemw*cos(beta)*r->displayYukawaMatrix(YD)(3,3)/r->displayGaugeCoupling(2);
    mtAtMH = pow(2,0.5)*polemw*sin(beta)*r->displayYukawaMatrix(YU)(3,3)/r->displayGaugeCoupling(2);
    mcAtMH = pow(2,0.5)*polemw*sin(beta)*r->displayYukawaMatrix(YU)(2,2)/r->displayGaugeCoupling(2);
    msAtMH = pow(2,0.5)*polemw*cos(beta)*r->displayYukawaMatrix(YD)(2,2)/r->displayGaugeCoupling(2);
    
    ///If yukawa matrix elements are 0 (can get if set mixing to 0 in softsusy spectrum generator) then use a non-zero mass - otherwise get inf issue in foftau function for loop decays!
    if (mbAtMH < 1.0e-66) {mbAtMH = mb;}
    if (mtAtMH < 1.0e-66) {mtAtMH = mt;}
    if (mcAtMH < 1.0e-66) {mcAtMH = mc;}
    if (msAtMH < 1.0e-66) {msAtMH = ms;}
    if (mbAtMA < 1.0e-66) {mbAtMA = mb;}
    if (mtAtMA < 1.0e-66) {mtAtMA = mt;}
    if (mcAtMA < 1.0e-66) {mcAtMA = mc;}
    if (msAtMA < 1.0e-66) {msAtMA = ms;}
    
    ///For mh should not use full susy running as susy not around this scale, instead use SM running as in QeDQcD object defined in lowe.h:
    QedQcd a(r->displayDataSet());
    
    a.runto(r->displayPhys().mh0(1));
    alphasAtMh = a.displayAlpha(ALPHAS);
    alphaAtMh =  a.displayAlpha(ALPHA);
    mbAtMh = a.displayMass()(6);
    mtAtMh = a.displayMass()(3);
    mcAtMh = a.displayMass()(2);
    msAtMh = a.displayMass()(5);
    
    g3atmh0 = pow(4*PI*alphasAtMh,0.5);
    g3atmH0 = pow(4*PI*alphasAtMH,0.5);
    g3atmA0 = pow(4*PI*alphasAtMA,0.5);
  }

 ///gravitino stuff - For NLSP decays to LSP gravitino, occur often in GMSB scenarios
 double mgravitino = r->displayGravitino();
 // mgravitino = 5.92500000e-08; ///Doesn't affect functions lsp(m,posi,posj) or nlsp(m,posi,posj) therefore these will read out the lsp and nlsp incorrectly with these, I just used to test gravitino decay formulae.
 int gravonoff = 0; ///a gravitino switch, by default it's off (0), it's turned on (1) a few lines below if the gravitino is the LSP as then decays to it are important, if it's not the LSP decays to it are unimportant as they are Planck suppressed.
 double MPlreduced = 2.4e18; ///Set by hand, could read in a value from softsusy spectrum generator
 int posi = 0, posj = 0;
 double m = mgravitino;
 double LSP = r->lsp(m, posi, posj); ///A number which indicates which type of particle is the LSP
 // LSP = -1; ///Temporary set to gravitino
 if (LSP == -1)     gravonoff = 1; ///< LSP is a gravitino
 else if( LSP > -1 && LSP < 7) gravonoff = 0; ///< LSP is not a gravitino
 else{throw("problem:LSP code not between -1 and 6 - issue reading LSP!");}
 int NLSP = 0; ///default position
 int neutNLSP = 1, upsquNLSP= 1, downsquNLSP = 1, slepNLSP = 1, snuNLSP = 1, gluNLSP = 1; ///Default position is 1 so all SUSY dedcays to LSP gravitino considered unless you only want the NLSP decays to gravitino decays, then if below section uncommented it sets all but NLSP switch to 0.
 
 ///Using the nlsp function from softsusy.cpp, Returns a label which says which particle is NLSP, 0 means NLSP is neutralino, 1=up squark, 2=down squark, 3=sleptons, 4=charginos, 5=sneutrinos, 6=gluino. Uncomment this if you only want the NLSP decays to the LSP gravitino to be considered, note one potential issue with this is particles only slightly heavier than the NLSP may still have only the decay to the LSPgravitino available and this function doesn't take account of this. If this section below is commented all SUSY decays to the LSP gravitino are considered and then they are not output if there BRs are less than BRTol
 NLSP = r->nlsp(m, posi, posj); 
 // fout << "NLSP = " << NLSP << endl;
 // NLSP = 0; /// Temporarily set to neutralino
 // neutNLSP = 0, upsquNLSP= 0, downsquNLSP = 0, slepNLSP = 0, snuNLSP = 0, gluNLSP = 0; /// For scans
 if( NLSP == 0) {
   //    fout << "NLSP is neutralino" << endl;
   upsquNLSP= 0; downsquNLSP = 0; slepNLSP = 0; snuNLSP = 0; gluNLSP = 0;
 }
 else if (NLSP == 1) {
   //    fout << "NLSP is up squark" << endl;
   neutNLSP = 0; downsquNLSP = 0; slepNLSP = 0; snuNLSP = 0; gluNLSP = 0;
 }
 else if (NLSP == 2) {
   //    fout << "NLSP is down squark" << endl;
   neutNLSP = 0; upsquNLSP= 0; slepNLSP = 0; snuNLSP = 0; gluNLSP = 0;
 }
 else if (NLSP == 3) {
   //    fout << "NLSP is slepton" << endl;
   neutNLSP = 0; upsquNLSP= 0; downsquNLSP = 0; snuNLSP = 0; gluNLSP = 0;
 }
 else if (NLSP == 4) {
   //    fout << "NLSP is chargino - WARNING chargino NLSP decays to gravitino LSP not included in program!" << endl;
   neutNLSP = 0; upsquNLSP= 0; downsquNLSP = 0; slepNLSP = 0; snuNLSP = 0; gluNLSP = 0;
 }
 else if (NLSP == 5) {
   //    fout << "NLSP is sneutrino" << endl;
   neutNLSP = 0; upsquNLSP= 0; downsquNLSP = 0; slepNLSP = 0; gluNLSP = 0;
 }
 else if (NLSP == 6) {
   //    fout << "NLSP is gluino" << endl;
   neutNLSP = 0; upsquNLSP= 0; downsquNLSP = 0; slepNLSP = 0; snuNLSP = 0;
 }
 else{throw("problem:NLSP code not between 0 and 6 - issue reading NLSP!");}
 
 ///Obtaining the CKM Matrix from Yukawa Matrices:
 DoubleMatrix YUU(3,3), YDD(3,3), YEE(3,3);
 
 YUU = r->displayYukawaMatrix(YU); YDD = r->displayYukawaMatrix(YD); YEE = r->displayYukawaMatrix(YE);
 
 DoubleMatrix u(3, 3), v(3, 3), vul(3,3), vur(3,3), vdl(3,3), vdr(3,3), vel(3,3), ver(3,3);
 DoubleVector ydDiag(3), yuDiag(3), yeDiag(3);
 r->displayYukawaMatrix(YU).diagonalise(u, v, yuDiag);
 vul = u.transpose(); vur = v.transpose();
 
 r->displayYukawaMatrix(YD).diagonalise(u, v, ydDiag);
 vdl = u.transpose(); vdr = v.transpose();
 
 r->displayYukawaMatrix(YE).diagonalise(u, v, yeDiag);
 vel = u.transpose(); ver = v.transpose();
 
 DoubleMatrix I3(3,3);
 for (int i=1; i<=3; i++) {
   for (int j=1; j<=3; j++) {
     if (i==j) {I3(i,j) = 1;}
     else {I3(i,j) = 0;}
   }
 }

 ///Form hermitian conjugate of vdl - only need to transpose as must be real as taken Yukawa matrices as real:
 DoubleMatrix vdlT = vdl.transpose();
 DoubleMatrix VCKM = vul*vdlT; ///CKM only used so far in H+ decays to q q'bar, otherwise taken as diagonal
 DoubleMatrix Vu(3,3), Uu(3,3), yu(3,3), Vd(3,3), Ud(3,3), yd(3,3), Ve(3,3), Ue(3,3), ye(3,3);
 

 double MCH1=fabs(mch(1)), MCH2=fabs(mch(2));  //BP FIX: not sure why it can be negative??
 
 ///Particle class used to store decay info in one place for ease of output into decay tables
 /// Create object ParticleGluino of class Particle, stores all decay info for gluino decays
 Particle ParticleGluino;
 ParticleGluino.name = "Gluino";
 ParticleGluino.PDG = PDGgluino;
 ParticleGluino.mass = mGluino;
 ParticleGluino.No_1to2_Decays = 24;
 ParticleGluino.No_1to3_Decays = 36;
 ParticleGluino.No_grav_Decays = 1;
 ParticleGluino.No_NMSSM_Decays = 0;
 ParticleGluino.No_of_Decays = ParticleGluino.No_1to2_Decays + ParticleGluino.No_1to3_Decays + ParticleGluino.No_grav_Decays;
 ParticleGluino.Array_Decays.resize(ParticleGluino.No_of_Decays);
 for (int i = 0; i < ParticleGluino.No_of_Decays; i++)
   ParticleGluino.Array_Decays[i].resize(6); 
 ParticleGluino.Array_Comments.resize(ParticleGluino.No_of_Decays);
 ParticleGluino.total_width = 0.0; ///initialise to 0
 ParticleGluino.two_width = 0.0;
 ParticleGluino.three_width = 0.0;
 
 /// Create object ParticleSdownL, stores all decay info for sdownL decays
 Particle ParticleSdownL;
 ParticleSdownL.name = "SdownL";
 ParticleSdownL.PDG = PDGsdownL;
 ParticleSdownL.mass = md(1,1);
 ParticleSdownL.No_1to2_Decays = 7;
 ParticleSdownL.No_1to3_Decays = 0;
 ParticleSdownL.No_grav_Decays = 1;
 ParticleSdownL.No_NMSSM_Decays = 1;
 ParticleSdownL.No_of_Decays = ParticleSdownL.No_1to2_Decays + ParticleSdownL.No_1to3_Decays + ParticleSdownL.No_grav_Decays + ParticleSdownL.No_NMSSM_Decays;
 ParticleSdownL.Array_Decays.resize(ParticleSdownL.No_of_Decays);
 for (int i = 0; i < ParticleSdownL.No_of_Decays; i++)
   ParticleSdownL.Array_Decays[i].resize(6); 
 ParticleSdownL.Array_Comments.resize(ParticleSdownL.No_of_Decays);
 ParticleSdownL.total_width = 0.0; ///initialise to 0
 ParticleSdownL.two_width = 0.0;
 ParticleSdownL.three_width = 0.0;
 
 ///Create object ParticleSdownR, stores all decay info for sdownR decays
 Particle ParticleSdownR;
 ParticleSdownR.name = "SdownR";
 ParticleSdownR.PDG = PDGsdownR;
 ParticleSdownR.mass = md(2,1);
 ParticleSdownR.No_1to2_Decays = 5;
 ParticleSdownR.No_1to3_Decays = 0;
 ParticleSdownR.No_grav_Decays = 1;
 ParticleSdownR.No_NMSSM_Decays = 1;
 ParticleSdownR.No_of_Decays = ParticleSdownR.No_1to2_Decays + ParticleSdownR.No_1to3_Decays + ParticleSdownR.No_grav_Decays + ParticleSdownR.No_NMSSM_Decays;
 ParticleSdownR.Array_Decays.resize(ParticleSdownR.No_of_Decays);
 for (int i = 0; i < ParticleSdownR.No_of_Decays; i++)
   ParticleSdownR.Array_Decays[i].resize(6);
 ParticleSdownR.Array_Comments.resize(ParticleSdownR.No_of_Decays);
 ParticleSdownR.total_width = 0.0; ///initialise to 0
 ParticleSdownR.two_width = 0.0;
 ParticleSdownR.three_width = 0.0;

 ///Create object ParticleSupL, stores all decay info for supL decays
 Particle ParticleSupL;
 ParticleSupL.name = "SupL";
 ParticleSupL.PDG = PDGsupL;
 ParticleSupL.mass = mu(1,1);
 ParticleSupL.No_1to2_Decays = 7;
 ParticleSupL.No_1to3_Decays = 0;
 ParticleSupL.No_grav_Decays = 1;
 ParticleSupL.No_NMSSM_Decays = 1;
 ParticleSupL.No_of_Decays = ParticleSupL.No_1to2_Decays + ParticleSupL.No_1to3_Decays + ParticleSupL.No_grav_Decays + ParticleSupL.No_NMSSM_Decays;
 ParticleSupL.Array_Decays.resize(ParticleSupL.No_of_Decays);
 for (int i = 0; i < ParticleSupL.No_of_Decays; i++)
   ParticleSupL.Array_Decays[i].resize(6);
 ParticleSupL.Array_Comments.resize(ParticleSupL.No_of_Decays);
 ParticleSupL.total_width = 0.0; ///initialise to 0
 ParticleSupL.two_width = 0.0;
 ParticleSupL.three_width = 0.0;
 
 ///Create object ParticleSupR, stores all decay info for supR decays
 Particle ParticleSupR;
 ParticleSupR.name = "SupR";
 ParticleSupR.PDG = PDGsupR;
 ParticleSupR.mass = mu(1,1);
 ParticleSupR.No_1to2_Decays = 5;
 ParticleSupR.No_1to3_Decays = 0;
 ParticleSupR.No_grav_Decays = 1;
 ParticleSupR.No_NMSSM_Decays = 1;
 ParticleSupR.No_of_Decays = ParticleSupR.No_1to2_Decays + ParticleSupR.No_1to3_Decays + ParticleSupR.No_grav_Decays + ParticleSupR.No_NMSSM_Decays;
 ParticleSupR.Array_Decays.resize(ParticleSupR.No_of_Decays);
 for (int i = 0; i < ParticleSupR.No_of_Decays; i++)
   ParticleSupR.Array_Decays[i].resize(6);
 ParticleSupR.Array_Comments.resize(ParticleSupR.No_of_Decays);
 ParticleSupR.total_width = 0.0; ///initialise to 0
 ParticleSupR.two_width = 0.0;
 ParticleSupR.three_width = 0.0;

 ///Create object ParticleSstrangeL, stores all decay info for sstrangeL decays
 Particle ParticleSstrangeL;
 ParticleSstrangeL.name = "SstrangeL";
 ParticleSstrangeL.PDG = PDGsstrangeL;
 ParticleSstrangeL.mass = md(2,1);
 ParticleSstrangeL.No_1to2_Decays = 7;
 ParticleSstrangeL.No_1to3_Decays = 0;
 ParticleSstrangeL.No_grav_Decays = 1;
 ParticleSstrangeL.No_NMSSM_Decays = 1;
 ParticleSstrangeL.No_of_Decays = ParticleSstrangeL.No_1to2_Decays + ParticleSstrangeL.No_1to3_Decays + ParticleSstrangeL.No_grav_Decays+ ParticleSstrangeL.No_NMSSM_Decays;
 ParticleSstrangeL.Array_Decays.resize(ParticleSstrangeL.No_of_Decays);
 for (int i = 0; i < ParticleSstrangeL.No_of_Decays; i++)
   ParticleSstrangeL.Array_Decays[i].resize(6);
 ParticleSstrangeL.Array_Comments.resize(ParticleSstrangeL.No_of_Decays);
 ParticleSstrangeL.total_width = 0.0; //initialise to 0
 ParticleSstrangeL.two_width = 0.0;
 ParticleSstrangeL.three_width = 0.0;

 ///Create object ParticleSstrangeR, stores all decay info for sstrangeR decays
 Particle ParticleSstrangeR;
 ParticleSstrangeR.name = "SstrangeR";
 ParticleSstrangeR.PDG = PDGsstrangeR;
 ParticleSstrangeR.mass = md(2,2);
 ParticleSstrangeR.No_1to2_Decays = 5;
 ParticleSstrangeR.No_1to3_Decays = 0;
 ParticleSstrangeR.No_grav_Decays = 1;
 ParticleSstrangeR.No_NMSSM_Decays = 1;
 ParticleSstrangeR.No_of_Decays = ParticleSstrangeR.No_1to2_Decays + ParticleSstrangeR.No_1to3_Decays + ParticleSstrangeR.No_grav_Decays+ ParticleSstrangeR.No_NMSSM_Decays;
 ParticleSstrangeR.Array_Decays.resize(ParticleSstrangeR.No_of_Decays);
 for (int i = 0; i < ParticleSstrangeR.No_of_Decays; i++)
   ParticleSstrangeR.Array_Decays[i].resize(6); 
 ParticleSstrangeR.Array_Comments.resize(ParticleSstrangeR.No_of_Decays);
 ParticleSstrangeR.total_width = 0.0; //initialise to 0
 ParticleSstrangeR.two_width = 0.0;
 ParticleSstrangeR.three_width = 0.0;
 
 ///Create object ParticleScharmL, stores all decay info for scharmL decays
 Particle ParticleScharmL;
 ParticleScharmL.name = "ScharmL";
 ParticleScharmL.PDG = PDGscharmL;
 ParticleScharmL.mass = mu(1,2);
 ParticleScharmL.No_1to2_Decays = 7;
 ParticleScharmL.No_1to3_Decays = 0;
 ParticleScharmL.No_grav_Decays = 1;
 ParticleScharmL.No_NMSSM_Decays = 1;
 ParticleScharmL.No_of_Decays = ParticleScharmL.No_1to2_Decays + ParticleScharmL.No_1to3_Decays + ParticleScharmL.No_grav_Decays + ParticleScharmL.No_NMSSM_Decays;
 ParticleScharmL.Array_Decays.resize(ParticleScharmL.No_of_Decays);
 for (int i = 0; i < ParticleScharmL.No_of_Decays; i++)
   ParticleScharmL.Array_Decays[i].resize(6);
 ParticleScharmL.Array_Comments.resize(ParticleScharmL.No_of_Decays);
 ParticleScharmL.total_width = 0.0; //initialise to 0
 ParticleScharmL.two_width = 0.0;
 ParticleScharmL.three_width = 0.0;

 ///Create object ParticleScharmR, stores all decay info for scharmR decays
 Particle ParticleScharmR;
 ParticleScharmR.name = "ScharmR";
 ParticleScharmR.PDG = PDGscharmR;
 ParticleScharmR.mass = mu(2,2);
 ParticleScharmR.No_1to2_Decays = 5;
 ParticleScharmR.No_1to3_Decays = 0;
 ParticleScharmR.No_grav_Decays = 1;
 ParticleScharmR.No_NMSSM_Decays = 1;
 ParticleScharmR.No_of_Decays = ParticleScharmR.No_1to2_Decays + ParticleScharmR.No_1to3_Decays + ParticleScharmR.No_grav_Decays + ParticleScharmR.No_NMSSM_Decays;
 ParticleScharmR.Array_Decays.resize(ParticleScharmR.No_of_Decays);
 for (int i = 0; i < ParticleScharmR.No_of_Decays; i++)
   ParticleScharmR.Array_Decays[i].resize(6); 
 ParticleScharmR.Array_Comments.resize(ParticleScharmR.No_of_Decays);
 ParticleScharmR.total_width = 0.0; //initialise to 0
 ParticleScharmR.two_width = 0.0;
 ParticleScharmR.three_width = 0.0;

 ///Create object ParticleSbottom1, stores all decay info for sbottom1 decays
 Particle ParticleSbottom1;
 ParticleSbottom1.name = "Sbottom1";
 ParticleSbottom1.PDG = PDGsbottom1;
 ParticleSbottom1.mass = md(1,3);
 ParticleSbottom1.No_1to2_Decays = 11;
 ParticleSbottom1.No_1to3_Decays = 0;
 ParticleSbottom1.No_grav_Decays = 1;
 ParticleSbottom1.No_NMSSM_Decays = 1;
 ParticleSbottom1.No_of_Decays = ParticleSbottom1.No_1to2_Decays + ParticleSbottom1.No_1to3_Decays + ParticleSbottom1.No_grav_Decays + ParticleSbottom1.No_NMSSM_Decays;
 ParticleSbottom1.Array_Decays.resize(ParticleSbottom1.No_of_Decays);
 for (int i = 0; i < ParticleSbottom1.No_of_Decays; i++)
   ParticleSbottom1.Array_Decays[i].resize(6);
 ParticleSbottom1.Array_Comments.resize(ParticleSbottom1.No_of_Decays);
 ParticleSbottom1.total_width = 0.0; //initialise to 0
 ParticleSbottom1.two_width = 0.0;
 ParticleSbottom1.three_width = 0.0;

 ///Create object ParticleSbottom2, stores all decay info for sbottom2 decays
 Particle ParticleSbottom2;
 ParticleSbottom2.name = "Sbottom2";
 ParticleSbottom2.PDG = PDGsbottom2;
 ParticleSbottom2.mass = md(2,3);
 ParticleSbottom2.No_1to2_Decays = 15;
 ParticleSbottom2.No_1to3_Decays = 0;
 ParticleSbottom2.No_grav_Decays = 1;
 ParticleSbottom2.No_NMSSM_Decays = 3;
 ParticleSbottom2.No_of_Decays = ParticleSbottom2.No_1to2_Decays + ParticleSbottom2.No_1to3_Decays + ParticleSbottom2.No_grav_Decays + ParticleSbottom2.No_NMSSM_Decays;
 ParticleSbottom2.Array_Decays.resize(ParticleSbottom2.No_of_Decays);
 for (int i = 0; i < ParticleSbottom2.No_of_Decays; i++)
   ParticleSbottom2.Array_Decays[i].resize(6);
 ParticleSbottom2.Array_Comments.resize(ParticleSbottom2.No_of_Decays);
 ParticleSbottom2.total_width = 0.0; //initialise to 0
 ParticleSbottom2.two_width = 0.0;
 ParticleSbottom2.three_width = 0.0;

 ///Create object ParticleStop1, stores all decay info for stop1 decays
 Particle ParticleStop1;
 ParticleStop1.name = "Stop1";
 ParticleStop1.PDG = PDGstop1;
 ParticleStop1.mass = mu(1,3);
 ParticleStop1.No_1to2_Decays = 11;
 ParticleStop1.No_1to3_Decays = 0;
 ParticleStop1.No_grav_Decays = 1;
 ParticleStop1.No_NMSSM_Decays = 1;
 ParticleStop1.No_of_Decays = ParticleStop1.No_1to2_Decays + ParticleStop1.No_1to3_Decays + ParticleStop1.No_grav_Decays + ParticleStop1.No_NMSSM_Decays;
 ParticleStop1.Array_Decays.resize(ParticleStop1.No_of_Decays);
 for (int i = 0; i < ParticleStop1.No_of_Decays; i++)
   ParticleStop1.Array_Decays[i].resize(6);
 ParticleStop1.Array_Comments.resize(ParticleStop1.No_of_Decays);
 ParticleStop1.total_width = 0.0; //initialise to 0
 ParticleStop1.two_width = 0.0;
 ParticleStop1.three_width = 0.0;

 ///Create object ParticleStop2, stores all decay info for stop2 decays
 Particle ParticleStop2;
 ParticleStop2.name = "Stop2";
 ParticleStop2.PDG = PDGstop2;
 ParticleStop2.mass = mu(2,3);
 ParticleStop2.No_1to2_Decays = 15;
 ParticleStop2.No_1to3_Decays = 0;
 ParticleStop2.No_grav_Decays = 1;
 ParticleStop2.No_NMSSM_Decays = 3;
 ParticleStop2.No_of_Decays = ParticleStop2.No_1to2_Decays + ParticleStop2.No_1to3_Decays + ParticleStop2.No_grav_Decays + ParticleStop2.No_NMSSM_Decays;
 ParticleStop2.Array_Decays.resize(ParticleStop2.No_of_Decays);
 for (int i = 0; i < ParticleStop2.No_of_Decays; i++)
   ParticleStop2.Array_Decays[i].resize(6);
 ParticleStop2.Array_Comments.resize(ParticleStop2.No_of_Decays);
 ParticleStop2.total_width = 0.0; //initialise to 0
 ParticleStop2.two_width = 0.0;
 ParticleStop2.three_width = 0.0;

 ///Create object ParticleSelectonL, stores all decay info for selectronL decays
 Particle ParticleSelectronL;
 ParticleSelectronL.name = "SelectronL";
 ParticleSelectronL.PDG = PDGselectronL;
 ParticleSelectronL.mass = me(1,1);
 ParticleSelectronL.No_1to2_Decays = 6;
 ParticleSelectronL.No_1to3_Decays = 0;
 ParticleSelectronL.No_grav_Decays = 1;
 ParticleSelectronL.No_NMSSM_Decays = 1;
 ParticleSelectronL.No_of_Decays = ParticleSelectronL.No_1to2_Decays + ParticleSelectronL.No_1to3_Decays + ParticleSelectronL.No_grav_Decays + ParticleSelectronL.No_NMSSM_Decays;
 ParticleSelectronL.Array_Decays.resize(ParticleSelectronL.No_of_Decays);
 for (int i = 0; i < ParticleSelectronL.No_of_Decays; i++)
   ParticleSelectronL.Array_Decays[i].resize(6);
 ParticleSelectronL.Array_Comments.resize(ParticleSelectronL.No_of_Decays);
 ParticleSelectronL.total_width = 0.0;
 ParticleSelectronL.two_width = 0.0;
 ParticleSelectronL.three_width = 0.0;

 ///Create object ParticleSelectronR, stores all decay info for selectronR decays
 Particle ParticleSelectronR;
 ParticleSelectronR.name = "SelectronR";
 ParticleSelectronR.PDG = PDGselectronR;
 ParticleSelectronR.mass = me(2,1);
 ParticleSelectronR.No_1to2_Decays = 4;
 ParticleSelectronR.No_1to3_Decays = 0;
 ParticleSelectronR.No_grav_Decays = 1;
 ParticleSelectronR.No_NMSSM_Decays = 1;
 ParticleSelectronR.No_of_Decays = ParticleSelectronR.No_1to2_Decays + ParticleSelectronR.No_1to3_Decays + ParticleSelectronR.No_grav_Decays + ParticleSelectronR.No_NMSSM_Decays;
 ParticleSelectronR.Array_Decays.resize(ParticleSelectronR.No_of_Decays);
 for (int i = 0; i < ParticleSelectronR.No_of_Decays; i++)
   ParticleSelectronR.Array_Decays[i].resize(6);
 ParticleSelectronR.Array_Comments.resize(ParticleSelectronR.No_of_Decays);
 ParticleSelectronR.total_width = 0.0;
 ParticleSelectronR.two_width = 0.0;
 ParticleSelectronR.three_width = 0.0;

 ///Create object ParticleSmuonL, stores all decay info for smuonL decays
 Particle ParticleSmuonL;
 ParticleSmuonL.name = "SmuonL";
 ParticleSmuonL.PDG = PDGsmuonL;
 ParticleSmuonL.mass = me(1,2);
 ParticleSmuonL.No_1to2_Decays = 6;
 ParticleSmuonL.No_1to3_Decays = 0;
 ParticleSmuonL.No_grav_Decays = 1;
 ParticleSmuonL.No_NMSSM_Decays = 1;
 ParticleSmuonL.No_of_Decays = ParticleSmuonL.No_1to2_Decays + ParticleSmuonL.No_1to3_Decays + ParticleSmuonL.No_grav_Decays + ParticleSmuonL.No_NMSSM_Decays;
 ParticleSmuonL.Array_Decays.resize(ParticleSmuonL.No_of_Decays);
 for (int i = 0; i < ParticleSmuonL.No_of_Decays; i++)
   ParticleSmuonL.Array_Decays[i].resize(6);
 ParticleSmuonL.Array_Comments.resize(ParticleSmuonL.No_of_Decays);
 ParticleSmuonL.total_width = 0.0;
 ParticleSmuonL.two_width = 0.0;
 ParticleSmuonL.three_width = 0.0;

 ///Create object ParticleSmuonR, stores all decay info for smuonR decays
 Particle ParticleSmuonR;
 ParticleSmuonR.name = "SmuonR";
 ParticleSmuonR.PDG = PDGsmuonR;
 ParticleSmuonR.mass = me(2,2);
 ParticleSmuonR.No_1to2_Decays = 4;
 ParticleSmuonR.No_1to3_Decays = 0;
 ParticleSmuonR.No_grav_Decays = 1;
 ParticleSmuonR.No_NMSSM_Decays = 1;
 ParticleSmuonR.No_of_Decays = ParticleSmuonR.No_1to2_Decays + ParticleSmuonR.No_1to3_Decays + ParticleSmuonR.No_grav_Decays + ParticleSmuonR.No_NMSSM_Decays;
 ParticleSmuonR.Array_Decays.resize(ParticleSmuonR.No_of_Decays);
 for (int i = 0; i < ParticleSmuonR.No_of_Decays; i++)
   ParticleSmuonR.Array_Decays[i].resize(6);
 ParticleSmuonR.Array_Comments.resize(ParticleSmuonR.No_of_Decays);
 ParticleSmuonR.total_width = 0.0;
 ParticleSmuonR.two_width = 0.0;
 ParticleSmuonR.three_width = 0.0;

 ///Create object ParticleSnue, stores all decay info for selectron sneutrino decays
 Particle ParticleSnue;
 ParticleSnue.name = "Selectron sneutrino";
 ParticleSnue.PDG = PDGnuselectronL;
 ParticleSnue.mass = msnu(1);
 ParticleSnue.No_1to2_Decays = 6;
 ParticleSnue.No_1to3_Decays = 0;
 ParticleSnue.No_grav_Decays = 1;
 ParticleSnue.No_NMSSM_Decays = 1;
 ParticleSnue.No_of_Decays = ParticleSnue.No_1to2_Decays + ParticleSnue.No_1to3_Decays + ParticleSnue.No_grav_Decays + ParticleSnue.No_NMSSM_Decays;
 ParticleSnue.Array_Decays.resize(ParticleSnue.No_of_Decays);
 for (int i = 0; i < ParticleSnue.No_of_Decays; i++)
    ParticleSnue.Array_Decays[i].resize(6);
 ParticleSnue.Array_Comments.resize(ParticleSnue.No_of_Decays);
 ParticleSnue.total_width = 0.0;
 ParticleSnue.two_width = 0.0;
 ParticleSnue.three_width = 0.0;

 ///Create object ParticleSnumu, stores all decay info for snumu decays
 Particle ParticleSnumu;
 ParticleSnumu.name = "Smuon sneutrino";
 ParticleSnumu.PDG = PDGnusmuonL;
 ParticleSnumu.mass = msnu(2);
 ParticleSnumu.No_1to2_Decays = 6;
 ParticleSnumu.No_1to3_Decays = 0;
 ParticleSnumu.No_grav_Decays = 1;
 ParticleSnumu.No_NMSSM_Decays = 1;
 ParticleSnumu.No_of_Decays = ParticleSnumu.No_1to2_Decays + ParticleSnumu.No_1to3_Decays + ParticleSnumu.No_grav_Decays + ParticleSnumu.No_NMSSM_Decays;
 ParticleSnumu.Array_Decays.resize(ParticleSnumu.No_of_Decays);
 for (int i = 0; i < ParticleSnumu.No_of_Decays; i++)
   ParticleSnumu.Array_Decays[i].resize(6);
 ParticleSnumu.Array_Comments.resize(ParticleSnumu.No_of_Decays);
 ParticleSnumu.total_width = 0.0;
 ParticleSnumu.two_width = 0.0;
 ParticleSnumu.three_width = 0.0;
 
 ///Create object ParticleStau1, stores all decay info for stau1 decays
 Particle ParticleStau1;
 ParticleStau1.name = "Stau1";
 ParticleStau1.PDG = PDGstau1;
 ParticleStau1.mass = me(1,3);
 ParticleStau1.No_1to2_Decays = 8;
 ParticleStau1.No_1to3_Decays = 0;
 ParticleStau1.No_grav_Decays = 1;
 ParticleStau1.No_NMSSM_Decays = 1;
 ParticleStau1.No_of_Decays = ParticleStau1.No_1to2_Decays + ParticleStau1.No_1to3_Decays + ParticleStau1.No_grav_Decays + ParticleStau1.No_NMSSM_Decays;
 ParticleStau1.Array_Decays.resize(ParticleStau1.No_of_Decays);
 for (int i = 0; i < ParticleStau1.No_of_Decays; i++)
   ParticleStau1.Array_Decays[i].resize(6);
 ParticleStau1.Array_Comments.resize(ParticleStau1.No_of_Decays);
 ParticleStau1.total_width = 0.0;
 ParticleStau1.two_width = 0.0;
 ParticleStau1.three_width = 0.0;

 ///Create object ParticleStau2, stores all decay info for stau2 decays
 Particle ParticleStau2;
 ParticleStau2.name = "Stau2";
 ParticleStau2.PDG = PDGstau2;
 ParticleStau2.mass = me(2,3);
 ParticleStau2.No_1to2_Decays = 12;
 ParticleStau2.No_1to3_Decays = 0;
 ParticleStau2.No_grav_Decays = 1;
 ParticleStau2.No_NMSSM_Decays = 3;
 ParticleStau2.No_of_Decays = ParticleStau2.No_1to2_Decays + ParticleStau2.No_1to3_Decays + ParticleStau2.No_grav_Decays + ParticleStau2.No_NMSSM_Decays;
 ParticleStau2.Array_Decays.resize(ParticleStau2.No_of_Decays);
 for (int i = 0; i < ParticleStau2.No_of_Decays; i++)
   ParticleStau2.Array_Decays[i].resize(6);
 ParticleStau2.Array_Comments.resize(ParticleStau2.No_of_Decays);
 ParticleStau2.total_width = 0.0;
 ParticleStau2.two_width = 0.0;
 ParticleStau2.three_width = 0.0;
  
 ///Create object ParticleSnutau, stores all decay info for snutau decays
 Particle ParticleSnutau;
 ParticleSnutau.name = "Stau sneutrino";
 ParticleSnutau.PDG = PDGnustauL;
 ParticleSnutau.mass = msnu(3);
 ParticleSnutau.No_1to2_Decays = 10;
 ParticleSnutau.No_1to3_Decays = 0;
 ParticleSnutau.No_grav_Decays = 1;
 ParticleSnutau.No_NMSSM_Decays = 1;
 ParticleSnutau.No_of_Decays = ParticleSnutau.No_1to2_Decays + ParticleSnutau.No_1to3_Decays + ParticleSnutau.No_grav_Decays + ParticleSnutau.No_NMSSM_Decays;
 ParticleSnutau.Array_Decays.resize(ParticleSnutau.No_of_Decays);
 for (int i = 0; i < ParticleSnutau.No_of_Decays; i++)
   ParticleSnutau.Array_Decays[i].resize(6);
 ParticleSnutau.Array_Comments.resize(ParticleSnutau.No_of_Decays);
 ParticleSnutau.total_width = 0.0;
 ParticleSnutau.two_width = 0.0;
 ParticleSnutau.three_width = 0.0;
  
 ///Create object ParticleChargino1, stores all decay info for chargino1 decays
 Particle ParticleChargino1;
 ParticleChargino1.name = "Chargino 1+ (lightest)";
 ParticleChargino1.PDG = PDGchargino1;
 ParticleChargino1.mass = MCH1;
 /// Ben: added in the charged pion/neutralino decay here
 ParticleChargino1.No_1to2_Decays = 24; ///We consider decays of the W1+, the W1- decays then just follow with the same amplitudes but often particles swapped for their anitparticles
 /// Added pi^+ pi^0 one here
 ParticleChargino1.No_1to3_Decays = 21;
 ParticleChargino1.No_grav_Decays = 0;
 ParticleChargino1.No_NMSSM_Decays = 2;
 ParticleChargino1.No_1to4_Decays = 1;
 ParticleChargino1.No_of_Decays = ParticleChargino1.No_1to2_Decays + ParticleChargino1.No_1to3_Decays + ParticleChargino1.No_grav_Decays + ParticleChargino1.No_NMSSM_Decays +
   ParticleChargino1.No_1to4_Decays; 
 ParticleChargino1.Array_Decays.resize(ParticleChargino1.No_of_Decays);
 for (int i = 0; i < ParticleChargino1.No_of_Decays; i++)
   ParticleChargino1.Array_Decays[i].resize(6);
 ParticleChargino1.Array_Comments.resize(ParticleChargino1.No_of_Decays);
 ParticleChargino1.total_width = 0.0;
 ParticleChargino1.two_width = 0.0;
 ParticleChargino1.three_width = 0.0;
 ParticleChargino1.four_width = 0.0;

 ///Create object ParticleChargino2, stores all decay info for chargino2 decays
 Particle ParticleChargino2;
 ParticleChargino2.name = "Chargino 2+ (heaviest)";
 ParticleChargino2.PDG = PDGchargino2;
 ParticleChargino2.mass = MCH2;
 ParticleChargino2.No_1to2_Decays = 27; /// has 4 additional decays cf chargino 1 as chargino2 -> chargino1 + Z/h/H/A as well
 ParticleChargino2.No_1to3_Decays = 20;
 ParticleChargino2.No_grav_Decays = 0;
 ParticleChargino2.No_NMSSM_Decays = 4;
 ParticleChargino2.No_of_Decays = ParticleChargino2.No_1to2_Decays + ParticleChargino2.No_1to3_Decays + ParticleChargino2.No_grav_Decays + ParticleChargino2.No_NMSSM_Decays;
 ParticleChargino2.Array_Decays.resize(ParticleChargino2.No_of_Decays);
 for (int i = 0; i < ParticleChargino2.No_of_Decays; i++)
   ParticleChargino2.Array_Decays[i].resize(6);
 ParticleChargino2.Array_Comments.resize(ParticleChargino2.No_of_Decays);
 ParticleChargino2.total_width = 0.0;
 ParticleChargino2.two_width = 0.0;
 ParticleChargino2.three_width = 0.0;
 
 ///Create object ParticleNeutralino1, stores all decay info for neutralino1 decays (remember neutralinos are mass-ordered)
 Particle ParticleNeutralino1;
 ParticleNeutralino1.name = "Neutralino1";
 ParticleNeutralino1.PDG = PDGneutralino1;
 ParticleNeutralino1.mass = mneut(1);
 ParticleNeutralino1.No_1to2_Decays = 62;
 ParticleNeutralino1.No_1to3_Decays = 10;
 ParticleNeutralino1.No_grav_Decays = 5;
 ParticleNeutralino1.No_NMSSM_Decays = 0;
 ParticleNeutralino1.No_of_Decays = ParticleNeutralino1.No_1to2_Decays + ParticleNeutralino1.No_1to3_Decays + ParticleNeutralino1.No_grav_Decays + ParticleNeutralino1.No_NMSSM_Decays;
 ParticleNeutralino1.Array_Decays.resize(ParticleNeutralino1.No_of_Decays);
 for (int i = 0; i < ParticleNeutralino1.No_of_Decays; i++)
   ParticleNeutralino1.Array_Decays[i].resize(6);
 ParticleNeutralino1.Array_Comments.resize(ParticleNeutralino1.No_of_Decays);
 ParticleNeutralino1.total_width = 0.0;
 ParticleNeutralino1.two_width = 0.0;
 ParticleNeutralino1.three_width = 0.0;

  ///Create object ParticleNeutralino2, stores all decay info for neutralino2 decays (remember neutralinos are mass-ordered)
 Particle ParticleNeutralino2;
 ParticleNeutralino2.name = "Neutralino2";
 ParticleNeutralino2.PDG = PDGneutralino2;
 ParticleNeutralino2.mass = mneut(2);
 ParticleNeutralino2.No_1to2_Decays = 62;
 ParticleNeutralino2.No_1to3_Decays = 12+10; 
 ParticleNeutralino2.No_grav_Decays = 5;
 ParticleNeutralino2.No_NMSSM_Decays = 2;
 ParticleNeutralino2.No_of_Decays = ParticleNeutralino2.No_1to2_Decays + ParticleNeutralino2.No_1to3_Decays + ParticleNeutralino2.No_grav_Decays + ParticleNeutralino2.No_NMSSM_Decays;
 ParticleNeutralino2.Array_Decays.resize(ParticleNeutralino2.No_of_Decays);
 for (int i = 0; i < ParticleNeutralino2.No_of_Decays; i++)
   ParticleNeutralino2.Array_Decays[i].resize(6); 
 ParticleNeutralino2.Array_Comments.resize(ParticleNeutralino2.No_of_Decays);
 ParticleNeutralino2.total_width = 0.0;
 ParticleNeutralino2.two_width = 0.0;
 ParticleNeutralino2.three_width = 0.0;
 
 ///Create object ParticleNeutralino3, stores all decay info for neutralino3 decays (remember neutralinos are mass-ordered)
 Particle ParticleNeutralino3;
 ParticleNeutralino3.name = "Neutralino3";
 ParticleNeutralino3.PDG = PDGneutralino3;
 ParticleNeutralino3.mass = mneut(3);
 ParticleNeutralino3.No_1to2_Decays = 62;
 ParticleNeutralino3.No_1to3_Decays = 24+10;
 ParticleNeutralino3.No_grav_Decays = 5;
 ParticleNeutralino3.No_NMSSM_Decays = 4;
 ParticleNeutralino3.No_of_Decays = ParticleNeutralino3.No_1to2_Decays + ParticleNeutralino3.No_1to3_Decays + ParticleNeutralino3.No_grav_Decays + ParticleNeutralino3.No_NMSSM_Decays;
 ParticleNeutralino3.Array_Decays.resize(ParticleNeutralino3.No_of_Decays);
 for (int i = 0; i < ParticleNeutralino3.No_of_Decays; i++)
   ParticleNeutralino3.Array_Decays[i].resize(6);
 ParticleNeutralino3.Array_Comments.resize(ParticleNeutralino3.No_of_Decays);
 ParticleNeutralino3.total_width = 0.0;
 ParticleNeutralino3.two_width = 0.0;
 ParticleNeutralino3.three_width = 0.0;

 ///Create object ParticleNeutralin4, stores all decay info for neutralino4 decays (remember neutralinos are mass-ordered)
 Particle ParticleNeutralino4;
 ParticleNeutralino4.name = "Neutralino4";
 ParticleNeutralino4.PDG = PDGneutralino4;
 ParticleNeutralino4.mass = mneut(4);
 ParticleNeutralino4.No_1to2_Decays = 62;
 ParticleNeutralino4.No_1to3_Decays = 46;
 ParticleNeutralino4.No_grav_Decays = 5;
 ParticleNeutralino4.No_NMSSM_Decays = 6;
 ParticleNeutralino4.No_of_Decays = ParticleNeutralino4.No_1to2_Decays + ParticleNeutralino4.No_1to3_Decays + ParticleNeutralino4.No_grav_Decays + ParticleNeutralino4.No_NMSSM_Decays;
 ParticleNeutralino4.Array_Decays.resize(ParticleNeutralino4.No_of_Decays);
 for (int i = 0; i < ParticleNeutralino4.No_of_Decays; i++)
   ParticleNeutralino4.Array_Decays[i].resize(6);
 ParticleNeutralino4.Array_Comments.resize(ParticleNeutralino4.No_of_Decays);
 ParticleNeutralino4.total_width = 0.0;
 ParticleNeutralino4.two_width = 0.0;
 ParticleNeutralino4.three_width = 0.0;
 
 ///Create object Particlehiggsl, stores all decay info for higgsl (lightest CP even neutral Higgs) decays
 Particle Particlehiggsl;
 Particlehiggsl.name = "light higgs";
 Particlehiggsl.PDG = PDGh0;
 Particlehiggsl.mass = mh0(1);
 Particlehiggsl.No_1to2_Decays = 68;
 Particlehiggsl.No_1to3_Decays = 2;
 Particlehiggsl.No_grav_Decays = 0;
 Particlehiggsl.No_NMSSM_Decays = 9;
 Particlehiggsl.No_of_Decays = Particlehiggsl.No_1to2_Decays + Particlehiggsl.No_1to3_Decays + Particlehiggsl.No_grav_Decays + Particlehiggsl.No_NMSSM_Decays;
 Particlehiggsl.Array_Decays.resize(Particlehiggsl.No_of_Decays);
 for (int i = 0; i < Particlehiggsl.No_of_Decays; i++)
   Particlehiggsl.Array_Decays[i].resize(6);
 Particlehiggsl.Array_Comments.resize(Particlehiggsl.No_of_Decays);
 Particlehiggsl.total_width = 0.0;
 Particlehiggsl.two_width = 0.0;
 Particlehiggsl.three_width = 0.0;

 ///Create object ParticleHiggsH, stores all decay info for HiggsH (heaviest CP even neutral Higgs of MSSM or second heaviest in NMSSM)
 Particle ParticleHiggsH;
 ParticleHiggsH.name = "heavy higgs";
 ParticleHiggsH.PDG = PDGH0;
 ParticleHiggsH.mass = mh0(2);
 ParticleHiggsH.No_1to2_Decays = 70;
 ParticleHiggsH.No_1to3_Decays = 2;
 ParticleHiggsH.No_grav_Decays = 0;
 ParticleHiggsH.No_NMSSM_Decays = 9;
 ParticleHiggsH.No_of_Decays = ParticleHiggsH.No_1to2_Decays + ParticleHiggsH.No_1to3_Decays + ParticleHiggsH.No_grav_Decays + ParticleHiggsH.No_NMSSM_Decays;
 ParticleHiggsH.Array_Decays.resize(ParticleHiggsH.No_of_Decays);
 for (int i = 0; i < ParticleHiggsH.No_of_Decays; i++)
   ParticleHiggsH.Array_Decays[i].resize(6);
 ParticleHiggsH.Array_Comments.resize(ParticleHiggsH.No_of_Decays);
 ParticleHiggsH.total_width = 0.0;
 ParticleHiggsH.two_width = 0.0;
 ParticleHiggsH.three_width = 0.0;

 ///Create object ParticleHiggsA, stores all decay info for HiggsA (CP odd neutral Higgs, or lightest CP odd neutral Higgs in NMSSM)
 Particle ParticleHiggsA;
 ParticleHiggsA.name = "pseudoscalar higgs";
 ParticleHiggsA.PDG = PDGA0;
 ParticleHiggsA.mass = mA0(1);
 ParticleHiggsA.No_1to2_Decays = 47;///Note A cannot decay into alike sfermion antisfermion paris because of CP conservation
 ParticleHiggsA.No_1to3_Decays = 0;
 ParticleHiggsA.No_grav_Decays = 0;
 ParticleHiggsA.No_NMSSM_Decays = 6;
 ParticleHiggsA.No_of_Decays = ParticleHiggsA.No_1to2_Decays + ParticleHiggsA.No_1to3_Decays + ParticleHiggsA.No_grav_Decays + ParticleHiggsA.No_NMSSM_Decays;
 ParticleHiggsA.Array_Decays.resize(ParticleHiggsA.No_of_Decays);
 for (int i = 0; i < ParticleHiggsA.No_of_Decays; i++)
   ParticleHiggsA.Array_Decays[i].resize(6);
 ParticleHiggsA.Array_Comments.resize(ParticleHiggsA.No_of_Decays);
 ParticleHiggsA.total_width = 0.0;
 ParticleHiggsA.two_width = 0.0;
 ParticleHiggsA.three_width = 0.0;

 ///Create object ParticleHiggsplus, stores all decay info for charged Higgs H+
 Particle ParticleHiggsplus;
 ParticleHiggsplus.name = "Charged higgs+";
 ParticleHiggsplus.PDG = PDGHplus;
 ParticleHiggsplus.mass = mHpm;
 ParticleHiggsplus.No_1to2_Decays = 39;
 ParticleHiggsplus.No_1to3_Decays = 0;
 ParticleHiggsplus.No_grav_Decays = 0;
 ParticleHiggsplus.No_NMSSM_Decays = 6;
 ParticleHiggsplus.No_of_Decays = ParticleHiggsplus.No_1to2_Decays + ParticleHiggsplus.No_1to3_Decays + ParticleHiggsplus.No_grav_Decays+ParticleHiggsplus.No_NMSSM_Decays;
 ParticleHiggsplus.Array_Decays.resize(ParticleHiggsplus.No_of_Decays);
 for (int i = 0; i < ParticleHiggsplus.No_of_Decays; i++)
   ParticleHiggsplus.Array_Decays[i].resize(6);
 ParticleHiggsplus.Array_Comments.resize(ParticleHiggsplus.No_of_Decays);
 ParticleHiggsplus.total_width = 0.0;
 ParticleHiggsplus.two_width = 0.0;
 ParticleHiggsplus.three_width = 0.0;
 
 ///Create also objects for NMSSM particles:
 ///Create object ParticleHiggsA2, stores all decay info for HiggsA2 (heaviest CP odd neutral Higgs of NMSSM)
 Particle ParticleHiggsA2;
 ParticleHiggsA2.name = "NMSSM second pseudoscalar higgs";
 ParticleHiggsA2.PDG = PDGA2;
 if (nmssmIsIt) ParticleHiggsA2.mass = mA0(2);
 ParticleHiggsA2.No_1to2_Decays = 0;
 ParticleHiggsA2.No_1to3_Decays = 0;
 ParticleHiggsA2.No_grav_Decays = 0;
 if (nmssmIsIt == false) { ParticleHiggsA2.No_NMSSM_Decays = 0;}
 else if (nmssmIsIt == true) { ParticleHiggsA2.No_NMSSM_Decays = 56;} ///Number of extra 1->2 decays when in NMSSM (47 like A1 then 5 similarly extra decays to neutralinos + A2 -> h3 Z + 3 A2 -> A1 h decays?????)
 ParticleHiggsA2.No_of_Decays = ParticleHiggsA2.No_1to2_Decays + ParticleHiggsA2.No_1to3_Decays + ParticleHiggsA2.No_grav_Decays + ParticleHiggsA2.No_NMSSM_Decays;
 ParticleHiggsA2.Array_Decays.resize(ParticleHiggsA2.No_of_Decays);
 for (int i = 0; i < ParticleHiggsA2.No_of_Decays; i++)
   ParticleHiggsA2.Array_Decays[i].resize(6);
 ParticleHiggsA2.Array_Comments.resize(ParticleHiggsA2.No_of_Decays);
 ParticleHiggsA2.total_width = 0.0;
 ParticleHiggsA2.two_width = 0.0;
 ParticleHiggsA2.three_width = 0.0;
 
 ///Create object ParticleHiggsH3, stores all decay info for HiggsH3 (Heaviest of the 3 CP even neutral Higgses of the NMSSM)
 Particle ParticleHiggsH3;
 ParticleHiggsH3.name = "NMSSM third scalar higgs";
 ParticleHiggsH3.PDG = PDGH3;
 if (nmssmIsIt) ParticleHiggsH3.mass = mh0(3);
 ParticleHiggsH3.No_1to2_Decays = 0;
 ParticleHiggsH3.No_1to3_Decays = 2;
 ParticleHiggsH3.No_grav_Decays = 0;
 if (nmssmIsIt == false) { ParticleHiggsH3.No_NMSSM_Decays = 0;}
 else if (nmssmIsIt == true) { ParticleHiggsH3.No_NMSSM_Decays = 81;}
 ParticleHiggsH3.No_of_Decays = ParticleHiggsH3.No_1to2_Decays + ParticleHiggsH3.No_1to3_Decays + ParticleHiggsH3.No_grav_Decays  + ParticleHiggsH3.No_NMSSM_Decays;
 ParticleHiggsH3.Array_Decays.resize(ParticleHiggsH3.No_of_Decays);
 for (int i = 0; i < ParticleHiggsH3.No_of_Decays; i++)
   ParticleHiggsH3.Array_Decays[i].resize(6);
 ParticleHiggsH3.Array_Comments.resize(ParticleHiggsH3.No_of_Decays);
 ParticleHiggsH3.total_width = 0.0;
 ParticleHiggsH3.two_width = 0.0;
 ParticleHiggsH3.three_width = 0.0;

 ///Create object ParticleNeutralino5, stores all decay info for neutralino5 decays (remember neutralinos are mass-ordered, so this is amix of singlino of NMSSM with the usual 4 neutralinos of MSSM)
 Particle ParticleNeutralino5;
 ParticleNeutralino5.name = "Neutralino5";
 ParticleNeutralino5.PDG = PDGneutralino5;
 if (nmssmIsIt) ParticleNeutralino5.mass = mneut(5);
 ParticleNeutralino5.No_1to2_Decays = 0;
 ParticleNeutralino5.No_1to3_Decays = 0; ///1->3 decays not included in NMSSM, only in MSSM
 ParticleNeutralino5.No_grav_Decays = 0; ///Decays to gravitinos not included in NMSSM, only in MSSM
 if (nmssmIsIt == false) { ParticleNeutralino5.No_NMSSM_Decays = 0;}
 else if (nmssmIsIt == true) { ParticleNeutralino5.No_NMSSM_Decays = 74;} 
 ParticleNeutralino5.No_of_Decays = ParticleNeutralino5.No_1to2_Decays + ParticleNeutralino5.No_1to3_Decays + ParticleNeutralino5.No_grav_Decays + ParticleNeutralino5.No_NMSSM_Decays;
 ParticleNeutralino5.Array_Decays.resize(ParticleNeutralino5.No_of_Decays);
 for (int i = 0; i < ParticleNeutralino5.No_of_Decays; i++)
   ParticleNeutralino5.Array_Decays[i].resize(6);
 ParticleNeutralino5.Array_Comments.resize(ParticleNeutralino5.No_of_Decays);
 ParticleNeutralino5.total_width = 0.0;
 ParticleNeutralino5.two_width = 0.0;
 ParticleNeutralino5.three_width = 0.0;
 
  g1 = g; g2 = gp; alphamix = alpha; betavac = beta;
  ///Gluino Decays
  if (flaggluino == 1) {

 /// Now need to calculate the partial decays of the gluino
  double gluinoamplitudeupantisupL = 0, gluinoamplitudeupantisupR = 0, gluinoamplitudeantiupsupL = 0, gluinoamplitudeantiupsupR = 0;
  gluinoamplitudeupantisupL = gluinoamplitudedecay(mGluino, mup, mu(1,1), alphas);
  gluinoamplitudeupantisupR = gluinoamplitudedecay (mGluino, mup, mu(2,1), alphas);
  gluinoamplitudeantiupsupL = gluinoamplitudedecay (mGluino, mup, mu(1,1), alphas);
  gluinoamplitudeantiupsupR = gluinoamplitudedecay (mGluino, mup, mu(2,1), alphas);
  
 double gluinoamplitudedownantisdownL = 0, gluinoamplitudedownantisdownR = 0, gluinoamplitudeantidownsdownL = 0, gluinoamplitudeantidownsdownR = 0;
 gluinoamplitudedownantisdownL = gluinoamplitudedecay(mGluino, mdo, md(1,1), alphas);
 gluinoamplitudedownantisdownR = gluinoamplitudedecay(mGluino, mdo, md(2,1), alphas);
 gluinoamplitudeantidownsdownL = gluinoamplitudedecay(mGluino, mdo, md(1,1), alphas);
 gluinoamplitudeantidownsdownR = gluinoamplitudedecay(mGluino, mdo, md(2,1), alphas);
 
 double gluinoamplitudecharmantischarmL = 0, gluinoamplitudecharmantischarmR = 0, gluinoamplitudeanticharmscharmL = 0, 	gluinoamplitudeanticharmscharmR = 0;
 gluinoamplitudecharmantischarmL = gluinoamplitudedecay(mGluino, mc, mu(1,2), alphas);
 gluinoamplitudecharmantischarmR = gluinoamplitudedecay(mGluino, mc, mu(2,2), alphas);
 gluinoamplitudeanticharmscharmL = gluinoamplitudedecay(mGluino, mc, mu(1,2), alphas);
 gluinoamplitudeanticharmscharmR = gluinoamplitudedecay(mGluino, mc, mu(2,2), alphas);
 
 double gluinoamplitudestrangeantisstrangeL = 0, gluinoamplitudestrangeantisstrangeR = 0, gluinoamplitudeantistrangesstrangeL = 0, 	gluinoamplitudeantistrangesstrangeR = 0;
 gluinoamplitudestrangeantisstrangeL = gluinoamplitudedecay(mGluino, ms, md(1,2), alphas);
 gluinoamplitudestrangeantisstrangeR = gluinoamplitudedecay(mGluino, ms, md(2,2), alphas);
 gluinoamplitudeantistrangesstrangeL = gluinoamplitudedecay(mGluino, ms, md(1,2), alphas);
 gluinoamplitudeantistrangesstrangeR = gluinoamplitudedecay(mGluino, ms, md(2,2), alphas); 

 double gluinoamplitudetopantistop1 = 0, gluinoamplitudetopantistop2 = 0, gluinoamplitudeantitopstop1 = 0, 	gluinoamplitudeantitopstop2 = 0;
 gluinoamplitudetopantistop1 = gluinoamplitudedecaymix(mGluino, mt, mu(1,3), alphas, 1, thetat);
 gluinoamplitudetopantistop2 = gluinoamplitudedecaymix(mGluino, mt, mu(2,3), alphas, 2, thetat);
 gluinoamplitudeantitopstop1 = gluinoamplitudedecaymix(mGluino, mt, mu(1,3), alphas, 1, thetat);
 gluinoamplitudeantitopstop2 = gluinoamplitudedecaymix(mGluino, mt, mu(2,3), alphas, 2, thetat);

 double gluinoamplitudebottomantisbottom1 = 0, gluinoamplitudebottomantisbottom2 = 0, gluinoamplitudeantibottomsbottom1 = 0, 	gluinoamplitudeantibottomsbottom2 = 0;
 gluinoamplitudebottomantisbottom1 = gluinoamplitudedecaymix(mGluino, mb, md(1,3), alphas, 1, thetab);
 gluinoamplitudebottomantisbottom2 = gluinoamplitudedecaymix(mGluino, mb, md(2,3), alphas, 2, thetab);
 gluinoamplitudeantibottomsbottom1 = gluinoamplitudedecaymix(mGluino, mb, md(1,3), alphas, 1, thetab);
 gluinoamplitudeantibottomsbottom2 = gluinoamplitudedecaymix(mGluino, mb, md(2,3), alphas, 2, thetab); 

 double gluinoamplitudeneut1uubar = 0, gluinoamplitudeneut2uubar = 0, gluinoamplitudeneut3uubar = 0, gluinoamplitudeneut4uubar = 0, gluinoamplitudeneut1ddbar = 0, gluinoamplitudeneut2ddbar = 0, gluinoamplitudeneut3ddbar = 0, gluinoamplitudeneut4ddbar = 0, gluinoamplitudeneut1ccbar = 0, gluinoamplitudeneut2ccbar = 0, gluinoamplitudeneut3ccbar = 0, gluinoamplitudeneut4ccbar = 0, gluinoamplitudeneut1ssbar = 0, gluinoamplitudeneut2ssbar = 0, gluinoamplitudeneut3ssbar = 0, gluinoamplitudeneut4ssbar = 0, gluinoamplitudeneut1ttbar = 0, gluinoamplitudeneut2ttbar = 0, gluinoamplitudeneut3ttbar = 0, gluinoamplitudeneut4ttbar = 0, gluinoamplitudeneut1bbbar = 0, gluinoamplitudeneut2bbbar = 0, gluinoamplitudeneut3bbbar = 0, gluinoamplitudeneut4bbbar = 0, gluinoamplitudegravitinogluon = 0;
 
 ///1 to 3 decays of gluinos to neutralinos and first two gen quarks via dgauss method:
  // gluinoamplitudeneut1uubar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(1), mu(1,1), mu(2,1), mup, g, gp, mixNeut, alphas, 'u', 1, onetothree);
 gluinoamplitudeneut1uubar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogenlimit (mGluino, mneut(1), mu(1,1), mu(2,1), mup, g, gp, mixNeut, alphas, 'u', 1, onetothree);
  gluinoamplitudeneut1ddbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogenlimit (mGluino, mneut(1), md(1,1), md(2,1), mdo, g, gp, mixNeut, alphas, 'd', 1, onetothree);
     // May get a PW < 0 in cases of very compressed spectra due to numerical precision, in this case take limit of compressed spectra so can avoid numerical issues associated with fine cancellations
 

 gluinoamplitudeneut2uubar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(2), mu(1,1), mu(2,1), mup, g, gp, mixNeut, alphas, 'u', 2, onetothree);
 gluinoamplitudeneut3uubar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(3), mu(1,1), mu(2,1), mup, g, gp, mixNeut, alphas, 'u', 3, onetothree);
 gluinoamplitudeneut4uubar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(4), mu(1,1), mu(2,1), mup, g, gp, mixNeut, alphas, 'u', 4, onetothree);
 // gluinoamplitudeneut1ddbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(1), md(1,1), md(2,1), mdo, g, gp, mixNeut, alphas, 'd', 1, onetothree);
  if (gluinoamplitudeneut1ddbar < 0) {
    gluinoamplitudeneut1ddbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(1), md(1,1), md(2,1), 0.0, g, gp, mixNeut, alphas, 'd', 1, onetothree); //Reason as for the gluinoamplitudeneut1uubar case
 }

 gluinoamplitudeneut2ddbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(2), md(1,1), md(2,1), mdo, g, gp, mixNeut, alphas, 'd', 2, onetothree);
 gluinoamplitudeneut3ddbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(3), md(1,1), md(2,1), mdo, g, gp, mixNeut, alphas, 'd', 3, onetothree);
 gluinoamplitudeneut4ddbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(4), md(1,1), md(2,1), mdo, g, gp, mixNeut, alphas, 'd', 4, onetothree);
 gluinoamplitudeneut1ccbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(1), mu(1,2), mu(2,2), mc, g, gp, mixNeut, alphas, 'u', 1, onetothree);
 gluinoamplitudeneut2ccbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(2), mu(1,2), mu(2,2), mc, g, gp, mixNeut, alphas, 'u', 2, onetothree);
 gluinoamplitudeneut3ccbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(3), mu(1,2), mu(2,2), mc, g, gp, mixNeut, alphas, 'u', 3, onetothree);
 gluinoamplitudeneut4ccbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(4), mu(1,2), mu(2,2), mc, g, gp, mixNeut, alphas, 'u', 4, onetothree);
 gluinoamplitudeneut1ssbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(1), md(1,2), md(2,2), ms, g, gp, mixNeut, alphas, 'd', 1, onetothree);
 gluinoamplitudeneut2ssbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(2), md(1,2), md(2,2), ms, g, gp, mixNeut, alphas, 'd', 2, onetothree);
 gluinoamplitudeneut3ssbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(3), md(1,2), md(2,2), ms, g, gp, mixNeut, alphas, 'd', 3, onetothree);
 gluinoamplitudeneut4ssbar = gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (mGluino, mneut(4), md(1,2), md(2,2), ms, g, gp, mixNeut, alphas, 'd', 4, onetothree);
 
 gluinoamplitudeneut1ttbar = gluinoamplitudedecaydgaussneutralinottbar (mGluino, mu(1,3), mu(2,3), mneut(1), mt, runmw, g, gp, thetat, beta, alphas, mixNeut, runmt, 1, onetothree, 't');
 gluinoamplitudeneut2ttbar = gluinoamplitudedecaydgaussneutralinottbar (mGluino, mu(1,3), mu(2,3), mneut(2), mt, runmw, g, gp, thetat, beta, alphas, mixNeut, runmt, 2, onetothree, 't'); 
 gluinoamplitudeneut3ttbar = gluinoamplitudedecaydgaussneutralinottbar (mGluino, mu(1,3), mu(2,3), mneut(3), mt, runmw, g, gp, thetat, beta, alphas, mixNeut, runmt, 3, onetothree, 't');
 gluinoamplitudeneut4ttbar = gluinoamplitudedecaydgaussneutralinottbar (mGluino, mu(1,3), mu(2,3), mneut(4), mt, runmw, g, gp, thetat, beta, alphas, mixNeut, runmt, 4, onetothree, 't');
 gluinoamplitudeneut1bbbar = gluinoamplitudedecaydgaussneutralinottbar (mGluino, md(1,3), md(2,3), mneut(1), mb, runmw, g, gp, thetab, beta, alphas, mixNeut, runmb, 1, onetothree, 'b');
 gluinoamplitudeneut2bbbar = gluinoamplitudedecaydgaussneutralinottbar (mGluino, md(1,3), md(2,3), mneut(2), mb, runmw, g, gp, thetab, beta, alphas, mixNeut, runmb, 2, onetothree, 'b');
 gluinoamplitudeneut3bbbar = gluinoamplitudedecaydgaussneutralinottbar (mGluino, md(1,3), md(2,3), mneut(3), mb, runmw, g, gp, thetab, beta, alphas, mixNeut, runmb, 3, onetothree, 'b'); 
 gluinoamplitudeneut4bbbar = gluinoamplitudedecaydgaussneutralinottbar (mGluino, md(1,3), md(2,3), mneut(4), mb, runmw, g, gp, thetab, beta, alphas, mixNeut, runmb, 4, onetothree, 'b');

 ///Note no 1to3 decays calculated in NMSSM, only in MSSM

 gluinoamplitudegravitinogluon = gluinoamplitudedecaygravitino (mGluino, mgravitino, MPlreduced, gravonoff, gluNLSP);
  
 double gluinoamplitudechar1udbar = 0, gluinoamplitudechar2udbar = 0, gluinoamplitudechar1csbar = 0, gluinoamplitudechar2csbar = 0, gluinoamplitudechar1tbbar = 0, gluinoamplitudechar2tbbar = 0;

 ///1 to 3 decays of gluinos to charginos and q qpbar using dgauss method
 gluinoamplitudechar1udbar = gluinoamplitudedecaydgausscharginoqqpbarfirsttwogen (mGluino, MCH1, mup, mdo, mu(1,1), md(1,1), g, thetaL2, thetaR2, alphas, 1, onetothree);
 if (gluinoamplitudechar1udbar < 0) {
   gluinoamplitudechar1udbar = gluinoamplitudedecaydgausscharginoqqpbarfirsttwogen (mGluino, MCH1, 0, 0, mu(1,1), md(1,1), g, thetaL2, thetaR2, alphas, 1, onetothree); //For very compressed spectra these partial widths contain very fine cancellations, which for small quark masses of the first generation can mean numerical precision issues give negative PWs, to avoid this we re-evaluate the PWs in this limit in the masssless quark limit for the first generation quarks, this has no real effect on the PW as for very compressed regions the regions of phase space where the finite quark mass would have an effect are suppressed
 }
   
 gluinoamplitudechar2udbar = gluinoamplitudedecaydgausscharginoqqpbarfirsttwogen (mGluino, MCH2, mup, mdo, mu(1,1), md(1,1), g, thetaL2, thetaR2, alphas, 2, onetothree);
 gluinoamplitudechar1csbar = gluinoamplitudedecaydgausscharginoqqpbarfirsttwogen (mGluino, MCH1, mc, ms, mu(1,2), md(1,2), g, thetaL2, thetaR2, alphas, 1, onetothree);
 gluinoamplitudechar2csbar = gluinoamplitudedecaydgausscharginoqqpbarfirsttwogen (mGluino, MCH2, mc, ms, mu(1,2), md(1,2), g, thetaL2, thetaR2, alphas, 2, onetothree);

  gluinoamplitudechar1tbbar = gluinoamplitudedecaydgausschartbbar (mGluino, mu(1,3), mu(2,3), md(1,3), md(2,3), mt, mb, MCH1, alphas, thetat, thetab, runmw, g, gp, thetaL2, thetaR2, beta, runmt, runmb, 1, onetothree);
  gluinoamplitudechar2tbbar = gluinoamplitudedecaydgausschartbbar (mGluino, mu(1,3), mu(2,3), md(1,3), md(2,3), mt, mb, MCH2, alphas, thetat, thetab, runmw, g, gp, thetaL2, thetaR2, beta, runmt, runmb, 2, onetothree);

 ///Now fill up Array_Decays for the Gluino
 ParticleGluino.Array_Decays[4][0] = PDGup; ParticleGluino.Array_Decays[4][1] = -PDGsupL; ParticleGluino.Array_Decays[4][2] =  gluinoamplitudeupantisupL; ParticleGluino.Array_Decays[4][3] = 2; ParticleGluino.Array_Comments[4] = "# ~g -> u ~u_L*"; 
 ParticleGluino.Array_Decays[5][0] = -PDGup; ParticleGluino.Array_Decays[5][1] = PDGsupL; ParticleGluino.Array_Decays[5][2] = gluinoamplitudeantiupsupL;ParticleGluino.Array_Decays[5][3] = 2; ParticleGluino.Array_Comments[5] = "# ~g -> ub ~u_L";
 ParticleGluino.Array_Decays[6][0] = PDGup;  ParticleGluino.Array_Decays[6][1] = -PDGsupR;  ParticleGluino.Array_Decays[6][2] = gluinoamplitudeupantisupR; ParticleGluino.Array_Decays[6][3] = 2; ParticleGluino.Array_Comments[6] = "# ~g -> u ~u_R*";
 ParticleGluino.Array_Decays[7][0] = -PDGup;  ParticleGluino.Array_Decays[7][1] = PDGsupR;  ParticleGluino.Array_Decays[7][2] = gluinoamplitudeantiupsupR; ParticleGluino.Array_Decays[7][3] = 2; ParticleGluino.Array_Comments[7] = "# ~g -> ub ~u_R";
 ParticleGluino.Array_Decays[0][0] = PDGdown;  ParticleGluino.Array_Decays[0][1] = -PDGsdownL;  ParticleGluino.Array_Decays[0][2] = gluinoamplitudedownantisdownL; ParticleGluino.Array_Decays[0][3] = 2; ParticleGluino.Array_Comments[0] = "# ~g -> d ~d_L*";
 ParticleGluino.Array_Decays[1][0] = -PDGdown;  ParticleGluino.Array_Decays[1][1] = PDGsdownL;  ParticleGluino.Array_Decays[1][2] = gluinoamplitudeantidownsdownL; ParticleGluino.Array_Decays[1][3] = 2; ParticleGluino.Array_Comments[1] = "# ~g -> db ~d_L";
 ParticleGluino.Array_Decays[2][0] = PDGdown;  ParticleGluino.Array_Decays[2][1] = -PDGsdownR;  ParticleGluino.Array_Decays[2][2] = gluinoamplitudedownantisdownR; ParticleGluino.Array_Decays[2][3] = 2; ParticleGluino.Array_Comments[2] = "# ~g -> d ~d_R*";
 ParticleGluino.Array_Decays[3][0] = -PDGdown;  ParticleGluino.Array_Decays[3][1] = PDGsdownR;  ParticleGluino.Array_Decays[3][2] = gluinoamplitudeantidownsdownR; ParticleGluino.Array_Decays[3][3] = 2; ParticleGluino.Array_Comments[3] = "# ~g -> db ~d_R";
 ParticleGluino.Array_Decays[12][0] = PDGcharm;  ParticleGluino.Array_Decays[12][1] = -PDGscharmL;  ParticleGluino.Array_Decays[12][2] = gluinoamplitudecharmantischarmL; ParticleGluino.Array_Decays[12][3] = 2; ParticleGluino.Array_Comments[12] = "# ~g -> c ~c_L*";
 ParticleGluino.Array_Decays[13][0] = -PDGcharm;  ParticleGluino.Array_Decays[13][1] = PDGscharmL;  ParticleGluino.Array_Decays[13][2] = gluinoamplitudeanticharmscharmL; ParticleGluino.Array_Decays[13][3] = 2; ParticleGluino.Array_Comments[13] = "# ~g -> cb ~c_L";
 ParticleGluino.Array_Decays[14][0] = PDGcharm;  ParticleGluino.Array_Decays[14][1] = -PDGscharmR;  ParticleGluino.Array_Decays[14][2] = gluinoamplitudecharmantischarmR; ParticleGluino.Array_Decays[14][3] = 2; ParticleGluino.Array_Comments[14] = "# ~g -> c ~c_R*";
 ParticleGluino.Array_Decays[15][0] = -PDGcharm;  ParticleGluino.Array_Decays[15][1] = PDGscharmR;  ParticleGluino.Array_Decays[15][2] = gluinoamplitudeanticharmscharmR; ParticleGluino.Array_Decays[15][3] = 2; ParticleGluino.Array_Comments[15] = "# ~g -> cb ~c_R";
 ParticleGluino.Array_Decays[8][0] = PDGstrange;  ParticleGluino.Array_Decays[8][1] = -PDGsstrangeL;  ParticleGluino.Array_Decays[8][2]=gluinoamplitudestrangeantisstrangeL; ParticleGluino.Array_Decays[8][3] = 2; ParticleGluino.Array_Comments[8] = "# ~g -> s ~s_L*";
 ParticleGluino.Array_Decays[9][0] = -PDGstrange;  ParticleGluino.Array_Decays[9][1] = PDGsstrangeL;  ParticleGluino.Array_Decays[9][2]=gluinoamplitudeantistrangesstrangeL; ParticleGluino.Array_Decays[9][3] = 2; ParticleGluino.Array_Comments[9] = "# ~g -> sb ~s_L";
 ParticleGluino.Array_Decays[10][0] = PDGstrange;  ParticleGluino.Array_Decays[10][1] = -PDGsstrangeR;  ParticleGluino.Array_Decays[10][2] = gluinoamplitudestrangeantisstrangeR; ParticleGluino.Array_Decays[10][3] = 2; ParticleGluino.Array_Comments[10] = "# ~g -> s ~s_R*";
 ParticleGluino.Array_Decays[11][0] = -PDGstrange;  ParticleGluino.Array_Decays[11][1] = PDGsstrangeR;  ParticleGluino.Array_Decays[11][2] = gluinoamplitudeantistrangesstrangeR; ParticleGluino.Array_Decays[11][3] = 2; ParticleGluino.Array_Comments[11] = "# ~g -> sb ~s_R";
 ParticleGluino.Array_Decays[20][0] = PDGtop;  ParticleGluino.Array_Decays[20][1] = -PDGstop1;  ParticleGluino.Array_Decays[20][2] = gluinoamplitudetopantistop1; ParticleGluino.Array_Decays[20][3] = 2; ParticleGluino.Array_Comments[20] = "# ~g -> t ~t_1*";
 ParticleGluino.Array_Decays[21][0] = -PDGtop;  ParticleGluino.Array_Decays[21][1] = PDGstop1;  ParticleGluino.Array_Decays[21][2] = gluinoamplitudeantitopstop1; ParticleGluino.Array_Decays[21][3] = 2; ParticleGluino.Array_Comments[21] = "# ~g -> tb ~t_1";
 ParticleGluino.Array_Decays[22][0] = PDGtop;  ParticleGluino.Array_Decays[22][1] = -PDGstop2;  ParticleGluino.Array_Decays[22][2] = gluinoamplitudetopantistop2; ParticleGluino.Array_Decays[22][3] = 2; ParticleGluino.Array_Comments[22] = "# ~g -> t ~t_2*";
 ParticleGluino.Array_Decays[23][0] = -PDGtop;  ParticleGluino.Array_Decays[23][1] = PDGstop2;  ParticleGluino.Array_Decays[23][2] = gluinoamplitudeantitopstop2; ParticleGluino.Array_Decays[23][3] = 2; ParticleGluino.Array_Comments[23] = "# ~g -> tb ~t_2";
 ParticleGluino.Array_Decays[16][0] = PDGbottom;  ParticleGluino.Array_Decays[16][1] = -PDGsbottom1;  ParticleGluino.Array_Decays[16][2] = gluinoamplitudebottomantisbottom1; ParticleGluino.Array_Decays[16][3] = 2; ParticleGluino.Array_Comments[16] = "# ~g -> b ~b_1*";
 ParticleGluino.Array_Decays[17][0] = -PDGbottom;  ParticleGluino.Array_Decays[17][1] = PDGsbottom1;  ParticleGluino.Array_Decays[17][2] = gluinoamplitudeantibottomsbottom1; ParticleGluino.Array_Decays[17][3] = 2; ParticleGluino.Array_Comments[17] = "# ~g -> bb ~b_1";
 ParticleGluino.Array_Decays[18][0] = PDGbottom;  ParticleGluino.Array_Decays[18][1] = -PDGsbottom2;  ParticleGluino.Array_Decays[18][2] = gluinoamplitudebottomantisbottom2; ParticleGluino.Array_Decays[18][3] = 2; ParticleGluino.Array_Comments[18] = "# ~g -> b ~b_2*";
 ParticleGluino.Array_Decays[19][0] = -PDGbottom;  ParticleGluino.Array_Decays[19][1] = PDGsbottom2;  ParticleGluino.Array_Decays[19][2] = gluinoamplitudeantibottomsbottom2; ParticleGluino.Array_Decays[19][3] = 2; ParticleGluino.Array_Comments[19] = "# ~g -> bb ~b_2";
 ParticleGluino.Array_Decays[24][0] = PDGgluon; ParticleGluino.Array_Decays[24][1] = PDGgravitino; ParticleGluino.Array_Decays[24][4] = 0; ParticleGluino.Array_Decays[24][2] =  gluinoamplitudegravitinogluon; ParticleGluino.Array_Decays[24][3] = 2; ParticleGluino.Array_Comments[24] = "# ~g -> g ~G";
 ParticleGluino.Array_Decays[25][0] = PDGneutralino1;ParticleGluino.Array_Decays[25][1] = PDGup; ParticleGluino.Array_Decays[25][4] = -PDGup; ParticleGluino.Array_Decays[25][2] = gluinoamplitudeneut1uubar; ParticleGluino.Array_Decays[25][3] = 3; ParticleGluino.Array_Comments[25] = "# ~g -> ~Z1 u ub";
 ParticleGluino.Array_Decays[26][0] = PDGneutralino2;ParticleGluino.Array_Decays[26][1] = PDGup; ParticleGluino.Array_Decays[26][4] = -PDGup; ParticleGluino.Array_Decays[26][2] = gluinoamplitudeneut2uubar; ParticleGluino.Array_Decays[26][3] = 3; ParticleGluino.Array_Comments[26] = "# ~g -> ~Z2 u ub";
 ParticleGluino.Array_Decays[27][0] = PDGneutralino3;ParticleGluino.Array_Decays[27][1] = PDGup; ParticleGluino.Array_Decays[27][4] = -PDGup; ParticleGluino.Array_Decays[27][2] = gluinoamplitudeneut3uubar; ParticleGluino.Array_Decays[27][3] = 3; ParticleGluino.Array_Comments[27] = "# ~g -> ~Z3 u ub";
 ParticleGluino.Array_Decays[28][0] = PDGneutralino4;ParticleGluino.Array_Decays[28][1] = PDGup; ParticleGluino.Array_Decays[28][4] = -PDGup; ParticleGluino.Array_Decays[28][2] = gluinoamplitudeneut4uubar; ParticleGluino.Array_Decays[28][3] = 3; ParticleGluino.Array_Comments[28] = "# ~g -> ~Z4 u ub";
 ParticleGluino.Array_Decays[29][0] = PDGneutralino1;ParticleGluino.Array_Decays[29][1] = PDGdown; ParticleGluino.Array_Decays[29][4] = -PDGdown; ParticleGluino.Array_Decays[29][2] = gluinoamplitudeneut1ddbar; ParticleGluino.Array_Decays[29][3] = 3; ParticleGluino.Array_Comments[29] = "# ~g -> ~Z1 d db";
 ParticleGluino.Array_Decays[30][0] = PDGneutralino2;ParticleGluino.Array_Decays[30][1] = PDGdown; ParticleGluino.Array_Decays[30][4] = -PDGdown; ParticleGluino.Array_Decays[30][2] = gluinoamplitudeneut2ddbar; ParticleGluino.Array_Decays[30][3] = 3; ParticleGluino.Array_Comments[30] = "# ~g -> ~Z2 d db";
 ParticleGluino.Array_Decays[31][0] = PDGneutralino3;ParticleGluino.Array_Decays[31][1] = PDGdown; ParticleGluino.Array_Decays[31][4] = -PDGdown; ParticleGluino.Array_Decays[31][2] = gluinoamplitudeneut3ddbar; ParticleGluino.Array_Decays[31][3] = 3; ParticleGluino.Array_Comments[31] = "# ~g -> ~Z3 d db"; 
 ParticleGluino.Array_Decays[32][0] = PDGneutralino4;ParticleGluino.Array_Decays[32][1] = PDGdown; ParticleGluino.Array_Decays[32][4] = -PDGdown; ParticleGluino.Array_Decays[32][2] = gluinoamplitudeneut4ddbar; ParticleGluino.Array_Decays[32][3] = 3; ParticleGluino.Array_Comments[32] = "# ~g -> ~Z4 d db";
 ParticleGluino.Array_Decays[33][0] = PDGneutralino1;ParticleGluino.Array_Decays[33][1] = PDGcharm; ParticleGluino.Array_Decays[33][4] = -PDGcharm; ParticleGluino.Array_Decays[33][2] = gluinoamplitudeneut1ccbar; ParticleGluino.Array_Decays[33][3] = 3; ParticleGluino.Array_Comments[33] = "# ~g -> ~Z1 c cb";
 ParticleGluino.Array_Decays[34][0] = PDGneutralino2;ParticleGluino.Array_Decays[34][1] = PDGcharm; ParticleGluino.Array_Decays[34][4] = -PDGcharm; ParticleGluino.Array_Decays[34][2] = gluinoamplitudeneut2ccbar; ParticleGluino.Array_Decays[34][3] = 3; ParticleGluino.Array_Comments[34] = "# ~g -> ~Z2 c cb";
 ParticleGluino.Array_Decays[35][0] = PDGneutralino3;ParticleGluino.Array_Decays[35][1] = PDGcharm; ParticleGluino.Array_Decays[35][4] = -PDGcharm; ParticleGluino.Array_Decays[35][2] = gluinoamplitudeneut3ccbar; ParticleGluino.Array_Decays[35][3] = 3; ParticleGluino.Array_Comments[35] = "# ~g -> ~Z3 c cb";
 ParticleGluino.Array_Decays[36][0] = PDGneutralino4;ParticleGluino.Array_Decays[36][1] = PDGcharm; ParticleGluino.Array_Decays[36][4] = -PDGcharm; ParticleGluino.Array_Decays[36][2] = gluinoamplitudeneut4ccbar; ParticleGluino.Array_Decays[36][3] = 3; ParticleGluino.Array_Comments[36] = "# ~g -> ~Z4 c cb";
 ParticleGluino.Array_Decays[37][0] = PDGneutralino1;ParticleGluino.Array_Decays[37][1] = PDGstrange; ParticleGluino.Array_Decays[37][4] = -PDGstrange; ParticleGluino.Array_Decays[37][2] = gluinoamplitudeneut1ssbar; ParticleGluino.Array_Decays[37][3] = 3; ParticleGluino.Array_Comments[37] = "# ~g -> ~Z1 s sb";
 ParticleGluino.Array_Decays[38][0] = PDGneutralino2;ParticleGluino.Array_Decays[38][1] = PDGstrange; ParticleGluino.Array_Decays[38][4] = -PDGstrange; ParticleGluino.Array_Decays[38][2] = gluinoamplitudeneut2ssbar; ParticleGluino.Array_Decays[38][3] = 3; ParticleGluino.Array_Comments[38] = "# ~g -> ~Z2 s sb";
 ParticleGluino.Array_Decays[39][0] = PDGneutralino3;ParticleGluino.Array_Decays[39][1] = PDGstrange; ParticleGluino.Array_Decays[39][4] = -PDGstrange; ParticleGluino.Array_Decays[39][2] = gluinoamplitudeneut3ssbar; ParticleGluino.Array_Decays[39][3] = 3; ParticleGluino.Array_Comments[39] = "# ~g -> ~Z3 s sb";
 ParticleGluino.Array_Decays[40][0] = PDGneutralino4;ParticleGluino.Array_Decays[40][1] = PDGstrange; ParticleGluino.Array_Decays[40][4] = -PDGstrange; ParticleGluino.Array_Decays[40][2] = gluinoamplitudeneut4ssbar; ParticleGluino.Array_Decays[40][3] = 3; ParticleGluino.Array_Comments[40] = "# ~g -> ~Z4 s sb";
 ParticleGluino.Array_Decays[41][0] = -PDGchargino1;ParticleGluino.Array_Decays[41][1] = PDGup; ParticleGluino.Array_Decays[41][4] = -PDGdown; ParticleGluino.Array_Decays[41][2] = gluinoamplitudechar1udbar; ParticleGluino.Array_Decays[41][3] = 3; ParticleGluino.Array_Comments[41] = "# ~g -> ~W1- u db";
 ParticleGluino.Array_Decays[42][0] = PDGchargino1;ParticleGluino.Array_Decays[42][1] = PDGdown; ParticleGluino.Array_Decays[42][4] = -PDGup; ParticleGluino.Array_Decays[42][2] = gluinoamplitudechar1udbar; ParticleGluino.Array_Decays[42][3] = 3; ParticleGluino.Array_Comments[42] = "# ~g -> ~W1+ ub d";
 ParticleGluino.Array_Decays[43][0] = -PDGchargino1;ParticleGluino.Array_Decays[43][1] = PDGcharm; ParticleGluino.Array_Decays[43][4] = -PDGstrange; ParticleGluino.Array_Decays[43][2] = gluinoamplitudechar1csbar; ParticleGluino.Array_Decays[43][3] = 3; ParticleGluino.Array_Comments[43] = "# ~g -> ~W1- c sb";
 ParticleGluino.Array_Decays[44][0] = PDGchargino1;ParticleGluino.Array_Decays[44][1] = PDGstrange; ParticleGluino.Array_Decays[44][4] = -PDGcharm; ParticleGluino.Array_Decays[44][2] = gluinoamplitudechar1csbar; ParticleGluino.Array_Decays[44][3] = 3; ParticleGluino.Array_Comments[44] = "# ~g -> ~W1+ cb s";
 ParticleGluino.Array_Decays[45][0] = -PDGchargino2;ParticleGluino.Array_Decays[45][1] = PDGup; ParticleGluino.Array_Decays[45][4] = -PDGdown; ParticleGluino.Array_Decays[45][2] = gluinoamplitudechar2udbar; ParticleGluino.Array_Decays[45][3] = 3; ParticleGluino.Array_Comments[45] = "# ~g -> ~W2- u db";
 ParticleGluino.Array_Decays[46][0] = PDGchargino2;ParticleGluino.Array_Decays[46][1] = PDGdown; ParticleGluino.Array_Decays[46][4] = -PDGup; ParticleGluino.Array_Decays[46][2] = gluinoamplitudechar2udbar; ParticleGluino.Array_Decays[46][3] = 3; ParticleGluino.Array_Comments[46] = "# ~g -> ~W2+ ub d";
 ParticleGluino.Array_Decays[47][0] = -PDGchargino2;ParticleGluino.Array_Decays[47][1] = PDGcharm; ParticleGluino.Array_Decays[47][4] = -PDGstrange; ParticleGluino.Array_Decays[47][2] = gluinoamplitudechar2csbar; ParticleGluino.Array_Decays[47][3] = 3; ParticleGluino.Array_Comments[47] = "# ~g -> ~W2- c sb";
 ParticleGluino.Array_Decays[48][0] = PDGchargino2;ParticleGluino.Array_Decays[48][1] = PDGstrange; ParticleGluino.Array_Decays[48][4] = -PDGcharm; ParticleGluino.Array_Decays[48][2] = gluinoamplitudechar2csbar; ParticleGluino.Array_Decays[48][3] = 3; ParticleGluino.Array_Comments[48] = "# ~g -> ~W2+ cb s";  
 ParticleGluino.Array_Decays[49][0] = PDGneutralino1; ParticleGluino.Array_Decays[49][1] = PDGtop; ParticleGluino.Array_Decays[49][4] = -PDGtop;  ParticleGluino.Array_Decays[49][2] = gluinoamplitudeneut1ttbar;  ParticleGluino.Array_Decays[49][3] = 3; ParticleGluino.Array_Comments[49] = "# ~g -> ~Z1 t tb";
 ParticleGluino.Array_Decays[50][0] = PDGneutralino2; ParticleGluino.Array_Decays[50][1] = PDGtop; ParticleGluino.Array_Decays[50][4] = -PDGtop;  ParticleGluino.Array_Decays[50][2] = gluinoamplitudeneut2ttbar;  ParticleGluino.Array_Decays[50][3] = 3; ParticleGluino.Array_Comments[50] = "# ~g -> ~Z2 t tb";
 ParticleGluino.Array_Decays[51][0] = PDGneutralino3; ParticleGluino.Array_Decays[51][1] = PDGtop; ParticleGluino.Array_Decays[51][4] = -PDGtop;  ParticleGluino.Array_Decays[51][2] = gluinoamplitudeneut3ttbar;  ParticleGluino.Array_Decays[51][3] = 3; ParticleGluino.Array_Comments[51] = "# ~g -> ~Z3 t tb";
 ParticleGluino.Array_Decays[52][0] = PDGneutralino4; ParticleGluino.Array_Decays[52][1] = PDGtop; ParticleGluino.Array_Decays[52][4] = -PDGtop;  ParticleGluino.Array_Decays[52][2] = gluinoamplitudeneut4ttbar;  ParticleGluino.Array_Decays[52][3] = 3; ParticleGluino.Array_Comments[52] = "# ~g -> ~Z4 t tb";
ParticleGluino.Array_Decays[53][0] = PDGneutralino1; ParticleGluino.Array_Decays[53][1] = PDGbottom; ParticleGluino.Array_Decays[53][4] = -PDGbottom;  ParticleGluino.Array_Decays[53][2] = gluinoamplitudeneut1bbbar;  ParticleGluino.Array_Decays[53][3] = 3; ParticleGluino.Array_Comments[53] = "# ~g -> ~Z1 b bb";
 ParticleGluino.Array_Decays[54][0] = PDGneutralino2; ParticleGluino.Array_Decays[54][1] = PDGbottom; ParticleGluino.Array_Decays[54][4] = -PDGbottom;  ParticleGluino.Array_Decays[54][2] = gluinoamplitudeneut2bbbar;  ParticleGluino.Array_Decays[54][3] = 3; ParticleGluino.Array_Comments[54] = "# ~g -> ~Z2 b bb";
 ParticleGluino.Array_Decays[55][0] = PDGneutralino3; ParticleGluino.Array_Decays[55][1] = PDGbottom; ParticleGluino.Array_Decays[55][4] = -PDGbottom;  ParticleGluino.Array_Decays[55][2] = gluinoamplitudeneut3bbbar;  ParticleGluino.Array_Decays[55][3] = 3; ParticleGluino.Array_Comments[55] = "# ~g -> ~Z3 b bb";
 ParticleGluino.Array_Decays[56][0] = PDGneutralino4; ParticleGluino.Array_Decays[56][1] = PDGbottom; ParticleGluino.Array_Decays[56][4] = -PDGbottom;  ParticleGluino.Array_Decays[56][2] = gluinoamplitudeneut4bbbar;  ParticleGluino.Array_Decays[56][3] = 3; ParticleGluino.Array_Comments[56] = "# ~g -> ~Z4 b bb";
 ParticleGluino.Array_Decays[57][0] = PDGchargino1; ParticleGluino.Array_Decays[57][1] = -PDGtop; ParticleGluino.Array_Decays[57][4] = PDGbottom; ParticleGluino.Array_Decays[57][2] =  gluinoamplitudechar1tbbar; ParticleGluino.Array_Decays[57][3] = 3; ParticleGluino.Array_Comments[57] = "# ~g -> ~W1+ tb b";
 ParticleGluino.Array_Decays[58][0] = -PDGchargino1; ParticleGluino.Array_Decays[58][1] = PDGtop; ParticleGluino.Array_Decays[58][4] = -PDGbottom; ParticleGluino.Array_Decays[58][2] =  gluinoamplitudechar1tbbar; ParticleGluino.Array_Decays[58][3] = 3; ParticleGluino.Array_Comments[58] = "# ~g -> ~W1- t bb";
 ParticleGluino.Array_Decays[59][0] = PDGchargino2; ParticleGluino.Array_Decays[59][1] = -PDGtop; ParticleGluino.Array_Decays[59][4] = PDGbottom; ParticleGluino.Array_Decays[59][2] =  gluinoamplitudechar2tbbar; ParticleGluino.Array_Decays[59][3] = 3; ParticleGluino.Array_Comments[59] = "# ~g -> ~W2+ tb b";
 ParticleGluino.Array_Decays[60][0] = -PDGchargino2; ParticleGluino.Array_Decays[60][1] = PDGtop; ParticleGluino.Array_Decays[60][4] = -PDGbottom; ParticleGluino.Array_Decays[60][2] =  gluinoamplitudechar2tbbar; ParticleGluino.Array_Decays[60][3] = 3; ParticleGluino.Array_Comments[60] = "# ~g -> ~W2- t bb";

 // std::cout << "Checking for Negative Partial Widths... " << std::endl;
 for(int i = 0; i<ParticleGluino.No_of_Decays; i++) {
   if (ParticleGluino.Array_Decays[i][2] < 0) {
     fout << "#warning! Partial Width for " << ParticleGluino.Array_Comments[i] << " is negative = " << ParticleGluino.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
     ParticleGluino.Array_Decays[i][2] = 0;
     errorflag = -1;
   }
 }

 double Gluino_No_1to2_Decays = 0;
 ParticleGluino.two_width = 0;
 ParticleGluino.three_width = 0;

  Gluino_No_1to2_Decays = ParticleGluino.No_1to2_Decays + ParticleGluino.No_grav_Decays;
   
 for (int j = 0; j<Gluino_No_1to2_Decays; j++) {
   ParticleGluino.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
 } 

 for (int j=0; j<Gluino_No_1to2_Decays; j++) {
    ParticleGluino.two_width = ParticleGluino.two_width + ParticleGluino.Array_Decays[j][2];
  }
 for (int j=Gluino_No_1to2_Decays; j<ParticleGluino.No_of_Decays; j++) {
   ParticleGluino.three_width = ParticleGluino.three_width + ParticleGluino.Array_Decays[j][2];
   // std::cout << "ParticleGluino.Array_Decays[" << j << "][2] = "<< ParticleGluino.Array_Decays[j][2] << std::endl;
 }
  
  if ( ParticleGluino.three_width != ParticleGluino.three_width) /// Tests for a nan as only nans aren't equal to themselves
    {
      fout << "# Three body decays give nan for gluino - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
      errorflag = -1;
      ParticleGluino.No_of_Decays = Gluino_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
      ParticleGluino.total_width = ParticleGluino.two_width;
    }
  else {
    ParticleGluino.total_width = ParticleGluino.two_width + ParticleGluino.three_width;
  }

  if ( ParticleGluino.total_width != ParticleGluino.total_width) /// Tests for a nan as only nans aren't equal to themselves
    {
      errorflag = -1;
      // for (int i = 0; i<ParticleGluino.No_of_Decays; i++) {
      //   fout << i << " " << ParticleGluino.Array_Decays[i][2] << endl;
      // }	  
      throw( "nan in gluino total width \n");
    }
  }

 
 ///Squark Decays
 ///SdownL Decays

 double sdownLamplitudegluinodown=0, sdownLamplitudedownneutralinoZ1=0, sdownLamplitudedownneutralinoZ2=0, sdownLamplitudedownneutralinoZ3=0, sdownLamplitudedownneutralinoZ4=0, sdownLamplitudecharginoW1up=0, sdownLamplitudecharginoW2up=0, sdownLamplitudedowngravitino=0;

 double sdownLamplitudedownneutralinoZ5 = 0;

 if (flagsdownL == 1) {
   sdownLamplitudegluinodown =  squarkamplitudedecaygluino (md(1,1), mdo, mGluino, alphas);
   sdownLamplitudecharginoW1up = squarkamplitudedecaycharginoW1 (md(1,1), mup, MCH1, g, thetaL2);
   sdownLamplitudecharginoW2up = squarkamplitudedecaycharginoW2 (md(1,1), mup, MCH2, g, thetaL2);
   if (nmssmIsIt == false) {
     sdownLamplitudedownneutralinoZ1 = squarkLamplitudedecayneutralino(md(1,1), mdo, mneut(1), g, gp, mixNeut, 1, -1);
     sdownLamplitudedownneutralinoZ2 = squarkLamplitudedecayneutralino(md(1,1), mdo, mneut(2), g, gp, mixNeut, 2, -1);
     sdownLamplitudedownneutralinoZ3 = squarkLamplitudedecayneutralino(md(1,1), mdo, mneut(3), g, gp, mixNeut, 3, -1);
     sdownLamplitudedownneutralinoZ4 = squarkLamplitudedecayneutralino(md(1,1), mdo, mneut(4), g, gp, mixNeut, 4, -1);
     sdownLamplitudedowngravitino = squarkamplitudedecaygravitino(md(1,1), mgravitino, mdo, MPlreduced, gravonoff, downsquNLSP); ///Gravitino decays only in MSSM
   }
   else if (nmssmIsIt == true) { ///Extended neutralino sector of NMSSM means decays to neutralinos are different in NMSSM
     sdownLamplitudedownneutralinoZ1 = squarkamplitudedecayquarkneutralinoNMSSM (md(1,1), mdo, mneut(1), g, gp, mixNeut, 'd', 'L', 1); 
     sdownLamplitudedownneutralinoZ2 = squarkamplitudedecayquarkneutralinoNMSSM (md(1,1), mdo, mneut(2), g, gp, mixNeut, 'd', 'L', 2);
     sdownLamplitudedownneutralinoZ3 = squarkamplitudedecayquarkneutralinoNMSSM (md(1,1), mdo, mneut(3), g, gp, mixNeut, 'd', 'L', 3);
     sdownLamplitudedownneutralinoZ4 = squarkamplitudedecayquarkneutralinoNMSSM (md(1,1), mdo, mneut(4), g, gp, mixNeut, 'd', 'L', 4);
     sdownLamplitudedownneutralinoZ5 = squarkamplitudedecayquarkneutralinoNMSSM (md(1,1), mdo, mneut(5), g, gp, mixNeut, 'd', 'L', 5);
   }


   ParticleSdownL.Array_Decays[0][0] = PDGdown; ParticleSdownL.Array_Decays[0][1] = PDGgluino; ParticleSdownL.Array_Decays[0][2] = sdownLamplitudegluinodown; ParticleSdownL.Array_Decays[0][3] = 2; ParticleSdownL.Array_Comments[0] = "# ~d_L -> d ~g";
   ParticleSdownL.Array_Decays[1][0] = PDGup; ParticleSdownL.Array_Decays[1][1] = -PDGchargino1; ParticleSdownL.Array_Decays[1][2] = sdownLamplitudecharginoW1up; ParticleSdownL.Array_Decays[1][3] = 2; ParticleSdownL.Array_Comments[1] = "# ~d_L -> u ~chi_1-" ;
   ParticleSdownL.Array_Decays[2][0] = PDGup; ParticleSdownL.Array_Decays[2][1] = -PDGchargino2; ParticleSdownL.Array_Decays[2][2] = sdownLamplitudecharginoW2up; ParticleSdownL.Array_Decays[2][3] = 2; ParticleSdownL.Array_Comments[2] = "# ~d_L -> u ~chi_2-";
   ParticleSdownL.Array_Decays[3][0] = PDGdown; ParticleSdownL.Array_Decays[3][1] = PDGneutralino1; ParticleSdownL.Array_Decays[3][2] = sdownLamplitudedownneutralinoZ1; ParticleSdownL.Array_Decays[3][3] = 2; ParticleSdownL.Array_Comments[3] = "# ~d_L -> d ~chi_10";
   ParticleSdownL.Array_Decays[4][0] = PDGdown; ParticleSdownL.Array_Decays[4][1] = PDGneutralino2; ParticleSdownL.Array_Decays[4][2] = sdownLamplitudedownneutralinoZ2; ParticleSdownL.Array_Decays[4][3] = 2; ParticleSdownL.Array_Comments[4] = "# ~d_L -> d ~chi_20";
   ParticleSdownL.Array_Decays[5][0] = PDGdown; ParticleSdownL.Array_Decays[5][1] = PDGneutralino3; ParticleSdownL.Array_Decays[5][2] = sdownLamplitudedownneutralinoZ3; ParticleSdownL.Array_Decays[5][3] = 2; ParticleSdownL.Array_Comments[5] = "# ~d_L -> d ~chi_30";
   ParticleSdownL.Array_Decays[6][0] = PDGdown; ParticleSdownL.Array_Decays[6][1] = PDGneutralino4; ParticleSdownL.Array_Decays[6][2] = sdownLamplitudedownneutralinoZ4; ParticleSdownL.Array_Decays[6][3] = 2; ParticleSdownL.Array_Comments[6] = "# ~d_L -> d ~chi_40";
   ParticleSdownL.Array_Decays[7][0] = PDGdown; ParticleSdownL.Array_Decays[7][1] = PDGneutralino5; ParticleSdownL.Array_Decays[7][2] = sdownLamplitudedownneutralinoZ5; ParticleSdownL.Array_Decays[7][3] = 2; ParticleSdownL.Array_Comments[7] = "# ~d_L -> d ~chi_50";
   ParticleSdownL.Array_Decays[8][0] = PDGdown; ParticleSdownL.Array_Decays[8][1] = PDGgravitino; ParticleSdownL.Array_Decays[8][2] = sdownLamplitudedowngravitino; ParticleSdownL.Array_Decays[8][3] = 2; ParticleSdownL.Array_Comments[8] = "# ~d_L -> d ~G";

   for(int i = 0; i<ParticleSdownL.No_of_Decays; i++) {
     if (ParticleSdownL.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSdownL.Array_Comments[i] << " is negative = " << ParticleSdownL.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSdownL.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }

   
 double SdownL_No_1to2_Decays = 0;
 
 SdownL_No_1to2_Decays = ParticleSdownL.No_1to2_Decays + ParticleSdownL.No_grav_Decays + ParticleSdownL.No_NMSSM_Decays;
 
 for (int j = 0; j<SdownL_No_1to2_Decays; j++) {
   ParticleSdownL.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
 }
 
 for (int j=0; j<SdownL_No_1to2_Decays; j++) {
   ParticleSdownL.two_width = ParticleSdownL.two_width + ParticleSdownL.Array_Decays[j][2];
 }
 for (int j=SdownL_No_1to2_Decays; j<ParticleSdownL.No_of_Decays; j++) {
   ParticleSdownL.three_width = ParticleSdownL.three_width + ParticleSdownL.Array_Decays[j][2];
 }
 
 ///Note currently no squark three body decays included - may change in future versions
 if ( ParticleSdownL.three_width != ParticleSdownL.three_width) /// Tests for a nan as only nans aren't equal to themselves
   {
     fout << "# Three body decays give nan for sdownL - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
     errorflag = -1;
     ParticleSdownL.No_of_Decays = SdownL_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
     ParticleSdownL.total_width = ParticleSdownL.two_width;
   }
 else {
   ParticleSdownL.total_width = ParticleSdownL.two_width + ParticleSdownL.three_width;
 }
 if ( ParticleSdownL.total_width != ParticleSdownL.total_width) /// Tests for a nan as only nans aren't equal to themselves
   {
     errorflag = -1;
     // for (int i = 0; i<ParticleSdownL.No_of_Decays; i++) {
     //   fout << i << " " << ParticleSdownL.Array_Decays[i][2] << endl;
     // }	  
     throw( "nan in sdownL total width \n");
   }
 
 }

  
 ///SdownR Decays

 double sdownRamplitudegluinodown=0, sdownRamplitudedownneutralinoZ1=0, sdownRamplitudedownneutralinoZ2=0, sdownRamplitudedownneutralinoZ3=0, sdownRamplitudedownneutralinoZ4=0, sdownRamplitudedowngravitino=0;
 
 double sdownRamplitudedownneutralinoZ5 = 0;
 
 if (flagsdownR == 1) {
   
   sdownRamplitudegluinodown = squarkamplitudedecaygluino (md(2,1), mdo, mGluino, alphas);
   if (nmssmIsIt == false) {
     sdownRamplitudedownneutralinoZ1 = squarkRamplitudedecayneutralino (md(2,1), mdo, mneut(1), g, gp, mixNeut, 1, -1);
     sdownRamplitudedownneutralinoZ2 = squarkRamplitudedecayneutralino (md(2,1), mdo, mneut(2), g, gp, mixNeut, 2, -1);
     sdownRamplitudedownneutralinoZ3 = squarkRamplitudedecayneutralino (md(2,1), mdo, mneut(3), g, gp, mixNeut, 3, -1);
     sdownRamplitudedownneutralinoZ4 = squarkRamplitudedecayneutralino (md(2,1), mdo, mneut(4), g, gp, mixNeut, 4, -1);
     sdownRamplitudedowngravitino = squarkamplitudedecaygravitino(md(2,1), mgravitino, mdo, MPlreduced, gravonoff, downsquNLSP);
   }
   else if (nmssmIsIt == true) {
     sdownRamplitudedownneutralinoZ1 = squarkamplitudedecayquarkneutralinoNMSSM (md(2,1), mdo, mneut(1), g, gp, mixNeut, 'd', 'R', 1); 
     sdownRamplitudedownneutralinoZ2 = squarkamplitudedecayquarkneutralinoNMSSM (md(2,1), mdo, mneut(2), g, gp, mixNeut, 'd', 'R', 2);
     sdownRamplitudedownneutralinoZ3 = squarkamplitudedecayquarkneutralinoNMSSM (md(2,1), mdo, mneut(3), g, gp, mixNeut, 'd', 'R', 3);
     sdownRamplitudedownneutralinoZ4 = squarkamplitudedecayquarkneutralinoNMSSM (md(2,1), mdo, mneut(4), g, gp, mixNeut, 'd', 'R', 4);
     sdownRamplitudedownneutralinoZ5 = squarkamplitudedecayquarkneutralinoNMSSM (md(2,1), mdo, mneut(5), g, gp, mixNeut, 'd', 'R', 5);
   }
   
   ParticleSdownR.Array_Decays[0][0] = PDGdown; ParticleSdownR.Array_Decays[0][1] = PDGgluino; ParticleSdownR.Array_Decays[0][2] = sdownRamplitudegluinodown; ParticleSdownR.Array_Decays[0][3] = 2; ParticleSdownR.Array_Comments[0] = "# ~d_R -> d ~g";
   ParticleSdownR.Array_Decays[1][0] = PDGdown; ParticleSdownR.Array_Decays[1][1] = PDGneutralino1; ParticleSdownR.Array_Decays[1][2] = sdownRamplitudedownneutralinoZ1; ParticleSdownR.Array_Decays[1][3] = 2; ParticleSdownR.Array_Comments[1] = "# ~d_R -> d ~chi_10";
   ParticleSdownR.Array_Decays[2][0] = PDGdown; ParticleSdownR.Array_Decays[2][1] = PDGneutralino2; ParticleSdownR.Array_Decays[2][2] = sdownRamplitudedownneutralinoZ2; ParticleSdownR.Array_Decays[2][3] = 2; ParticleSdownR.Array_Comments[2] = "# ~d_R -> d ~chi_20";
   ParticleSdownR.Array_Decays[3][0] = PDGdown; ParticleSdownR.Array_Decays[3][1] = PDGneutralino3; ParticleSdownR.Array_Decays[3][2] = sdownRamplitudedownneutralinoZ3; ParticleSdownR.Array_Decays[3][3] = 2; ParticleSdownR.Array_Comments[3] = "# ~d_R -> d ~chi_30";
   ParticleSdownR.Array_Decays[4][0] = PDGdown; ParticleSdownR.Array_Decays[4][1] = PDGneutralino4; ParticleSdownR.Array_Decays[4][2] = sdownRamplitudedownneutralinoZ4; ParticleSdownR.Array_Decays[4][3] = 2; ParticleSdownR.Array_Comments[4] = "# ~d_R -> d ~chi_40";
   ParticleSdownR.Array_Decays[5][0] = PDGdown; ParticleSdownR.Array_Decays[5][1] = PDGneutralino4; ParticleSdownR.Array_Decays[5][2] = sdownRamplitudedownneutralinoZ5; ParticleSdownR.Array_Decays[5][3] = 2; ParticleSdownR.Array_Comments[5] = "# ~d_R -> d ~chi_50";
   
   ParticleSdownR.Array_Decays[6][0] = PDGdown; ParticleSdownR.Array_Decays[6][1] = PDGgravitino; ParticleSdownR.Array_Decays[6][2] = sdownRamplitudedowngravitino; ParticleSdownR.Array_Decays[6][3] = 2; ParticleSdownR.Array_Comments[6] = "# ~d_R -> d ~G";

   for(int i = 0; i<ParticleSdownR.No_of_Decays; i++) {
     if (ParticleSdownR.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSdownR.Array_Comments[i] << " is negative = " << ParticleSdownR.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSdownR.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   

   double SdownR_No_1to2_Decays = 0;
   
   SdownR_No_1to2_Decays = ParticleSdownR.No_1to2_Decays + ParticleSdownR.No_grav_Decays + ParticleSdownR.No_NMSSM_Decays;
   
   for (int j = 0; j<SdownR_No_1to2_Decays; j++) {
     ParticleSdownR.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<SdownR_No_1to2_Decays; j++) {
     ParticleSdownR.two_width = ParticleSdownR.two_width + ParticleSdownR.Array_Decays[j][2];
   }
   for (int j=SdownR_No_1to2_Decays; j<ParticleSdownR.No_of_Decays; j++) {
     ParticleSdownR.three_width = ParticleSdownR.three_width + ParticleSdownR.Array_Decays[j][2];
   }
   
   ///Note currently no squark three body decays included - may change in future versions
   if ( ParticleSdownR.three_width != ParticleSdownR.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for sdownR - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSdownR.No_of_Decays = SdownR_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSdownR.total_width = ParticleSdownR.two_width;
     }
   else {
     ParticleSdownR.total_width = ParticleSdownR.two_width + ParticleSdownR.three_width;
   }
   
   if ( ParticleSdownR.total_width != ParticleSdownR.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSdownR.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSdownR.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in sdownR total width \n");
     }
   
 }
 
 ///SupL Decays
 
 double supLamplitudegluinoup=0, supLamplitudeupneutralinoZ1=0, supLamplitudeupneutralinoZ2=0, supLamplitudeupneutralinoZ3=0, supLamplitudeupneutralinoZ4=0, supLamplitudecharginoW1down=0, supLamplitudecharginoW2down=0, supLamplitudeupgravitino=0; 
 double supLamplitudeupneutralinoZ5 = 0;
 
 if (flagsupL == 1) {
   supLamplitudegluinoup = squarkamplitudedecaygluino (mu(1,1), mup, mGluino, alphas);
   supLamplitudecharginoW1down = squarkamplitudedecaycharginoW1 (mu(1,1), mdo, MCH1, g, thetaR2);
   supLamplitudecharginoW2down = squarkamplitudedecaycharginoW2 (mu(1,1), mdo, MCH2, g, thetaR2);
   if (nmssmIsIt == false) {
     supLamplitudeupneutralinoZ1 = squarkLamplitudedecayneutralino (mu(1,1), mup, mneut(1), g, gp, mixNeut, 1, 1);
     supLamplitudeupneutralinoZ2 = squarkLamplitudedecayneutralino (mu(1,1), mup, mneut(2), g, gp, mixNeut, 2, 1);
     supLamplitudeupneutralinoZ3 = squarkLamplitudedecayneutralino (mu(1,1), mup, mneut(3), g, gp, mixNeut, 3, 1);
     supLamplitudeupneutralinoZ4 = squarkLamplitudedecayneutralino (mu(1,1), mup, mneut(4), g, gp, mixNeut, 4, 1);
     supLamplitudeupgravitino = squarkamplitudedecaygravitino (mu(1,1), mgravitino, mup, MPlreduced, gravonoff, upsquNLSP);
   }
   else if (nmssmIsIt == true) {
     supLamplitudeupneutralinoZ1 = squarkamplitudedecayquarkneutralinoNMSSM (mu(1,1), mup, mneut(1), g, gp, mixNeut, 'u', 'L', 1); 
     supLamplitudeupneutralinoZ2 = squarkamplitudedecayquarkneutralinoNMSSM (mu(1,1), mup, mneut(2), g, gp, mixNeut, 'u', 'L', 2);
     supLamplitudeupneutralinoZ3 = squarkamplitudedecayquarkneutralinoNMSSM (mu(1,1), mup, mneut(3), g, gp, mixNeut, 'u', 'L', 3);
     supLamplitudeupneutralinoZ4 = squarkamplitudedecayquarkneutralinoNMSSM (mu(1,1), mup, mneut(4), g, gp, mixNeut, 'u', 'L', 4);
     supLamplitudeupneutralinoZ5 = squarkamplitudedecayquarkneutralinoNMSSM (mu(1,1), mup, mneut(5), g, gp, mixNeut, 'u', 'L', 5);
   }
   
   ParticleSupL.Array_Decays[0][0] = PDGup; ParticleSupL.Array_Decays[0][1] = PDGgluino; ParticleSupL.Array_Decays[0][2] = supLamplitudegluinoup; ParticleSupL.Array_Decays[0][3] = 2; ParticleSupL.Array_Comments[0] = "# ~u_L -> u ~g";
   ParticleSupL.Array_Decays[1][0] = PDGdown; ParticleSupL.Array_Decays[1][1] = PDGchargino1; ParticleSupL.Array_Decays[1][2] = supLamplitudecharginoW1down; ParticleSupL.Array_Decays[1][3] = 2; ParticleSupL.Array_Comments[1] = "# ~u_L -> d ~chi_1+"; 
   ParticleSupL.Array_Decays[2][0] = PDGdown; ParticleSupL.Array_Decays[2][1] = PDGchargino2; ParticleSupL.Array_Decays[2][2] = supLamplitudecharginoW2down; ParticleSupL.Array_Decays[2][3] = 2; ParticleSupL.Array_Comments[2] = "# ~u_L -> d ~chi_2+"; 
   ParticleSupL.Array_Decays[3][0] = PDGup; ParticleSupL.Array_Decays[3][1] = PDGneutralino1; ParticleSupL.Array_Decays[3][2] = supLamplitudeupneutralinoZ1; ParticleSupL.Array_Decays[3][3] = 2; ParticleSupL.Array_Comments[3] = "# ~u_L -> u ~chi_10";
   ParticleSupL.Array_Decays[4][0] = PDGup; ParticleSupL.Array_Decays[4][1] = PDGneutralino2; ParticleSupL.Array_Decays[4][2] = supLamplitudeupneutralinoZ2; ParticleSupL.Array_Decays[4][3] = 2; ParticleSupL.Array_Comments[4] = "# ~u_L -> u ~chi_20";
   ParticleSupL.Array_Decays[5][0] = PDGup; ParticleSupL.Array_Decays[5][1] = PDGneutralino3; ParticleSupL.Array_Decays[5][2] = supLamplitudeupneutralinoZ3; ParticleSupL.Array_Decays[5][3] = 2; ParticleSupL.Array_Comments[5] = "# ~u_L -> u ~chi_30";
   ParticleSupL.Array_Decays[6][0] = PDGup; ParticleSupL.Array_Decays[6][1] = PDGneutralino4; ParticleSupL.Array_Decays[6][2] = supLamplitudeupneutralinoZ4; ParticleSupL.Array_Decays[6][3] = 2; ParticleSupL.Array_Comments[6] = "# ~u_L -> u ~chi_40";
   ParticleSupL.Array_Decays[7][0] = PDGup; ParticleSupL.Array_Decays[7][1] = PDGneutralino5; ParticleSupL.Array_Decays[7][2] = supLamplitudeupneutralinoZ5; ParticleSupL.Array_Decays[7][3] = 2; ParticleSupL.Array_Comments[7] = "# ~u_L -> u ~chi_50";
   
   ParticleSupL.Array_Decays[8][0] = PDGup; ParticleSupL.Array_Decays[8][1] = PDGgravitino; ParticleSupL.Array_Decays[8][2] = supLamplitudeupgravitino; ParticleSupL.Array_Decays[8][3] = 2; ParticleSupL.Array_Comments[8] = "# ~u_L -> u ~G";

   for(int i = 0; i<ParticleSupL.No_of_Decays; i++) {
     if (ParticleSupL.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSupL.Array_Comments[i] << " is negative = " << ParticleSupL.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSupL.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }
   
   
   double SupL_No_1to2_Decays = 0;
   
   SupL_No_1to2_Decays = ParticleSupL.No_1to2_Decays + ParticleSupL.No_grav_Decays + ParticleSupL.No_NMSSM_Decays;
   
   for (int j = 0; j<SupL_No_1to2_Decays; j++) {
     ParticleSupL.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<SupL_No_1to2_Decays; j++) {
     ParticleSupL.two_width = ParticleSupL.two_width + ParticleSupL.Array_Decays[j][2];
   }
   for (int j=SupL_No_1to2_Decays; j<ParticleSupL.No_of_Decays; j++) {
     ParticleSupL.three_width = ParticleSupL.three_width + ParticleSupL.Array_Decays[j][2];
   }
   
   ///Note currently no squark three body decays included - may change in future versions
   if ( ParticleSupL.three_width != ParticleSupL.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for supL - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSupL.No_of_Decays = SupL_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSupL.total_width = ParticleSupL.two_width;
     }
   else {
     ParticleSupL.total_width = ParticleSupL.two_width + ParticleSupL.three_width;
   }
   
   if ( ParticleSupL.total_width != ParticleSupL.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSupL.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSupL.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in supL total width \n");
     }
 
 }

 ///SupR Decays
 double supRamplitudegluinoup=0, supRamplitudeupneutralinoZ1=0, supRamplitudeupneutralinoZ2=0, supRamplitudeupneutralinoZ3=0, supRamplitudeupneutralinoZ4=0, supRamplitudeupgravitino=0;
 double supRamplitudeupneutralinoZ5=0;

 if(flagsupR == 1) {
   
   supRamplitudegluinoup = squarkamplitudedecaygluino (mu(2,1), mup, mGluino, alphas);
   if (nmssmIsIt == false) {
     supRamplitudeupneutralinoZ1 = squarkRamplitudedecayneutralino (mu(2,1), mup, mneut(1), g, gp, mixNeut, 1, 1);
     supRamplitudeupneutralinoZ2 = squarkRamplitudedecayneutralino (mu(2,1), mup, mneut(2), g, gp, mixNeut, 2, 1);
     supRamplitudeupneutralinoZ3 = squarkRamplitudedecayneutralino (mu(2,1), mup, mneut(3), g, gp, mixNeut, 3, 1);
     supRamplitudeupneutralinoZ4 = squarkRamplitudedecayneutralino (mu(2,1), mup, mneut(4), g, gp, mixNeut, 4, 1);
     supRamplitudeupgravitino = squarkamplitudedecaygravitino(mu(2,1), mgravitino, mup, MPlreduced, gravonoff, upsquNLSP);
   }
   else if (nmssmIsIt == true) {
     supRamplitudeupneutralinoZ1 = squarkamplitudedecayquarkneutralinoNMSSM (mu(2,1), mup, mneut(1), g, gp, mixNeut, 'u', 'R', 1); 
     supRamplitudeupneutralinoZ2 = squarkamplitudedecayquarkneutralinoNMSSM (mu(2,1), mup, mneut(2), g, gp, mixNeut, 'u', 'R', 2);
     supRamplitudeupneutralinoZ3 = squarkamplitudedecayquarkneutralinoNMSSM (mu(2,1), mup, mneut(3), g, gp, mixNeut, 'u', 'R', 3);
     supRamplitudeupneutralinoZ4 = squarkamplitudedecayquarkneutralinoNMSSM (mu(2,1), mup, mneut(4), g, gp, mixNeut, 'u', 'R', 4);
     supRamplitudeupneutralinoZ5 = squarkamplitudedecayquarkneutralinoNMSSM (mu(2,1), mup, mneut(5), g, gp, mixNeut, 'u', 'R', 5);
   }
   
   ParticleSupR.Array_Decays[0][0] = PDGup; ParticleSupR.Array_Decays[0][1] = PDGgluino; ParticleSupR.Array_Decays[0][2] = supRamplitudegluinoup; ParticleSupR.Array_Decays[0][3] = 2; ParticleSupR.Array_Comments[0] = "# ~u_R -> u ~g";
   ParticleSupR.Array_Decays[1][0] = PDGup; ParticleSupR.Array_Decays[1][1] = PDGneutralino1; ParticleSupR.Array_Decays[1][2] = supRamplitudeupneutralinoZ1; ParticleSupR.Array_Decays[1][3] = 2; ParticleSupR.Array_Comments[1] = "# ~u_R -> u ~chi_10";
   ParticleSupR.Array_Decays[2][0] = PDGup; ParticleSupR.Array_Decays[2][1] = PDGneutralino2; ParticleSupR.Array_Decays[2][2] = supRamplitudeupneutralinoZ2; ParticleSupR.Array_Decays[2][3] = 2; ParticleSupR.Array_Comments[2] = "# ~u_R -> u ~chi_20";
   ParticleSupR.Array_Decays[3][0] = PDGup; ParticleSupR.Array_Decays[3][1] = PDGneutralino3; ParticleSupR.Array_Decays[3][2] = supRamplitudeupneutralinoZ3; ParticleSupR.Array_Decays[3][3] = 2; ParticleSupR.Array_Comments[3] = "# ~u_R -> u ~chi_30";
   ParticleSupR.Array_Decays[4][0] = PDGup; ParticleSupR.Array_Decays[4][1] = PDGneutralino4; ParticleSupR.Array_Decays[4][2] = supRamplitudeupneutralinoZ4; ParticleSupR.Array_Decays[4][3] = 2; ParticleSupR.Array_Comments[4] = "# ~u_R -> u ~chi_40";
   ParticleSupR.Array_Decays[5][0] = PDGup; ParticleSupR.Array_Decays[5][1] = PDGneutralino5; ParticleSupR.Array_Decays[5][2] = supRamplitudeupneutralinoZ5; ParticleSupR.Array_Decays[5][3] = 2; ParticleSupR.Array_Comments[5] = "# ~u_R -> u ~chi_50";
   
   ParticleSupR.Array_Decays[6][0] = PDGup; ParticleSupR.Array_Decays[6][1] = PDGgravitino; ParticleSupR.Array_Decays[6][2] = supRamplitudeupgravitino; ParticleSupR.Array_Decays[6][3] = 2; ParticleSupR.Array_Comments[6] = "# ~u_R -> u ~G";

   for(int i = 0; i<ParticleSupR.No_of_Decays; i++) {
     if (ParticleSupR.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSupR.Array_Comments[i] << " is negative = " << ParticleSupR.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSupR.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   

   double SupR_No_1to2_Decays = 0;
   
   SupR_No_1to2_Decays = ParticleSupR.No_1to2_Decays + ParticleSupR.No_grav_Decays + ParticleSupR.No_NMSSM_Decays;
 
   for (int j = 0; j<SupR_No_1to2_Decays; j++) {
     ParticleSupR.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<SupR_No_1to2_Decays; j++) {
     ParticleSupR.two_width = ParticleSupR.two_width + ParticleSupR.Array_Decays[j][2];
   }
   for (int j=SupR_No_1to2_Decays; j<ParticleSupR.No_of_Decays; j++) {
     ParticleSupR.three_width = ParticleSupR.three_width + ParticleSupR.Array_Decays[j][2];
   }
   ///Note currently no squark three body decays included - may change in future versions
   if ( ParticleSupR.three_width != ParticleSupR.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for supR - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSupR.No_of_Decays = SupR_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSupR.total_width = ParticleSupR.two_width;
     }
   else {
     ParticleSupR.total_width = ParticleSupR.two_width + ParticleSupR.three_width;
   }
   
   if ( ParticleSupR.total_width != ParticleSupR.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSupR.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSupR.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in supR total width \n");
     }
  
 }

 ///SstrangeL Decays

 double sstrangeLamplitudegluinostrange=0, sstrangeLamplitudestrangeneutralinoZ1=0, sstrangeLamplitudestrangeneutralinoZ2=0, sstrangeLamplitudestrangeneutralinoZ3=0, sstrangeLamplitudestrangeneutralinoZ4=0, sstrangeLamplitudecharginoW1charm=0, sstrangeLamplitudecharginoW2charm=0, sstrangeLamplitudestrangegravitino=0;
 double sstrangeLamplitudestrangeneutralinoZ5 = 0;
 
 if(flagsstrangeL == 1) {
   sstrangeLamplitudegluinostrange = squarkamplitudedecaygluino (md(1,2), ms, mGluino, alphas);
   sstrangeLamplitudecharginoW1charm = squarkamplitudedecaycharginoW1 (md(1,2), mc, MCH1, g, thetaL2);
   sstrangeLamplitudecharginoW2charm = squarkamplitudedecaycharginoW2 (md(1,2), mc, MCH2, g, thetaL2);
   if (nmssmIsIt == false) {
     sstrangeLamplitudestrangeneutralinoZ1 = squarkLamplitudedecayneutralino (md(1,2), ms, mneut(1), g, gp, mixNeut, 1, -1);
     sstrangeLamplitudestrangeneutralinoZ2 = squarkLamplitudedecayneutralino (md(1,2), ms, mneut(2), g, gp, mixNeut, 2, -1);
     sstrangeLamplitudestrangeneutralinoZ3 = squarkLamplitudedecayneutralino (md(1,2), ms, mneut(3), g, gp, mixNeut, 3, -1);
     sstrangeLamplitudestrangeneutralinoZ4 = squarkLamplitudedecayneutralino (md(1,2), ms, mneut(4), g, gp, mixNeut, 4, -1);
     
     sstrangeLamplitudestrangegravitino = squarkamplitudedecaygravitino (md(1,2), mgravitino, ms, MPlreduced, gravonoff, downsquNLSP);
   }
   else if (nmssmIsIt == true) {
     sstrangeLamplitudestrangeneutralinoZ1 = squarkamplitudedecayquarkneutralinoNMSSM (md(1,2), ms, mneut(1), g, gp, mixNeut, 'd', 'L', 1); 
     sstrangeLamplitudestrangeneutralinoZ2 = squarkamplitudedecayquarkneutralinoNMSSM (md(1,2), ms, mneut(2), g, gp, mixNeut, 'd', 'L', 2);
     sstrangeLamplitudestrangeneutralinoZ3 = squarkamplitudedecayquarkneutralinoNMSSM (md(1,2), ms, mneut(3), g, gp, mixNeut, 'd', 'L', 3);
     sstrangeLamplitudestrangeneutralinoZ4 = squarkamplitudedecayquarkneutralinoNMSSM (md(1,2), ms, mneut(4), g, gp, mixNeut, 'd', 'L', 4);
     sstrangeLamplitudestrangeneutralinoZ5 = squarkamplitudedecayquarkneutralinoNMSSM (md(1,2), ms, mneut(5), g, gp, mixNeut, 'd', 'L', 5);
   }
   
   ParticleSstrangeL.Array_Decays[0][0] = PDGstrange; ParticleSstrangeL.Array_Decays[0][1] = PDGgluino; ParticleSstrangeL.Array_Decays[0][2] = sstrangeLamplitudegluinostrange; ParticleSstrangeL.Array_Decays[0][3] = 2; ParticleSstrangeL.Array_Comments[0] = "# ~s_L -> s ~g";
   ParticleSstrangeL.Array_Decays[1][0] = PDGcharm; ParticleSstrangeL.Array_Decays[1][1] = -PDGchargino1; ParticleSstrangeL.Array_Decays[1][2] = sstrangeLamplitudecharginoW1charm; ParticleSstrangeL.Array_Decays[1][3] = 2; ParticleSstrangeL.Array_Comments[1] = "# ~s_L -> c ~chi_1-" ;
   ParticleSstrangeL.Array_Decays[2][0] = PDGcharm; ParticleSstrangeL.Array_Decays[2][1] = -PDGchargino2; ParticleSstrangeL.Array_Decays[2][2] = sstrangeLamplitudecharginoW2charm; ParticleSstrangeL.Array_Decays[2][3] = 2; ParticleSstrangeL.Array_Comments[2] = "# ~s_L -> c ~chi_2-";
   ParticleSstrangeL.Array_Decays[3][0] = PDGstrange; ParticleSstrangeL.Array_Decays[3][1] = PDGneutralino1; ParticleSstrangeL.Array_Decays[3][2] = sstrangeLamplitudestrangeneutralinoZ1; ParticleSstrangeL.Array_Decays[3][3] = 2; ParticleSstrangeL.Array_Comments[3] = "# ~s_L -> s ~chi_10";
   ParticleSstrangeL.Array_Decays[4][0] = PDGstrange; ParticleSstrangeL.Array_Decays[4][1] = PDGneutralino2; ParticleSstrangeL.Array_Decays[4][2] = sstrangeLamplitudestrangeneutralinoZ2; ParticleSstrangeL.Array_Decays[4][3] = 2; ParticleSstrangeL.Array_Comments[4] = "# ~s_L -> s ~chi_20";
   ParticleSstrangeL.Array_Decays[5][0] = PDGstrange; ParticleSstrangeL.Array_Decays[5][1] = PDGneutralino3; ParticleSstrangeL.Array_Decays[5][2] = sstrangeLamplitudestrangeneutralinoZ3; ParticleSstrangeL.Array_Decays[5][3] = 2; ParticleSstrangeL.Array_Comments[5] = "# ~s_L -> s ~chi_30";
   ParticleSstrangeL.Array_Decays[6][0] = PDGstrange; ParticleSstrangeL.Array_Decays[6][1] = PDGneutralino4; ParticleSstrangeL.Array_Decays[6][2] = sstrangeLamplitudestrangeneutralinoZ4; ParticleSstrangeL.Array_Decays[6][3] = 2; ParticleSstrangeL.Array_Comments[6] = "# ~s_L -> s ~chi_40";
   ParticleSstrangeL.Array_Decays[7][0] = PDGstrange; ParticleSstrangeL.Array_Decays[7][1] = PDGneutralino5; ParticleSstrangeL.Array_Decays[7][2] = sstrangeLamplitudestrangeneutralinoZ5; ParticleSstrangeL.Array_Decays[7][3] = 2; ParticleSstrangeL.Array_Comments[7] = "# ~s_L -> s ~chi_50";

   ParticleSstrangeL.Array_Decays[8][0] = PDGstrange; ParticleSstrangeL.Array_Decays[8][1] = PDGgravitino; ParticleSstrangeL.Array_Decays[8][2] = sstrangeLamplitudestrangegravitino; ParticleSstrangeL.Array_Decays[8][3] = 2; ParticleSstrangeL.Array_Comments[8] = "# ~s_L -> s ~G";

   for(int i = 0; i<ParticleSstrangeL.No_of_Decays; i++) {
     if (ParticleSstrangeL.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSstrangeL.Array_Comments[i] << " is negative = " << ParticleSstrangeL.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSstrangeL.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }
      
   double SstrangeL_No_1to2_Decays = 0;

   SstrangeL_No_1to2_Decays = ParticleSstrangeL.No_1to2_Decays + ParticleSstrangeL.No_grav_Decays + ParticleSstrangeL.No_NMSSM_Decays;
 
   for (int j = 0; j<SstrangeL_No_1to2_Decays; j++) {
     ParticleSstrangeL.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<SstrangeL_No_1to2_Decays; j++) {
     ParticleSstrangeL.two_width = ParticleSstrangeL.two_width + ParticleSstrangeL.Array_Decays[j][2];
   }
   for (int j=SstrangeL_No_1to2_Decays; j<ParticleSstrangeL.No_of_Decays; j++) {
     ParticleSstrangeL.three_width = ParticleSstrangeL.three_width + ParticleSstrangeL.Array_Decays[j][2];
   }
   ///Note currently no squark three body decays included - may change in future versions 
   if ( ParticleSstrangeL.three_width != ParticleSstrangeL.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for sstrangeL - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSstrangeL.No_of_Decays = SstrangeL_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSstrangeL.total_width = ParticleSstrangeL.two_width;
     }
   else {
     ParticleSstrangeL.total_width = ParticleSstrangeL.two_width + ParticleSstrangeL.three_width;
   }
   
   if ( ParticleSstrangeL.total_width != ParticleSstrangeL.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSstrangeL.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSstrangeL.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in sstrangeL total width \n");
     }
   
 }


 ///SstrangeR decays

 double sstrangeRamplitudegluinostrange=0, sstrangeRamplitudestrangeneutralinoZ1=0, sstrangeRamplitudestrangeneutralinoZ2=0, sstrangeRamplitudestrangeneutralinoZ3=0, sstrangeRamplitudestrangeneutralinoZ4=0, sstrangeRamplitudestrangegravitino=0;

 double sstrangeRamplitudestrangeneutralinoZ5 = 0;

 if (flagsstrangeR == 1) {
   sstrangeRamplitudegluinostrange = squarkamplitudedecaygluino (md(2,2), ms, mGluino, alphas);
   if (nmssmIsIt == false) {
     sstrangeRamplitudestrangeneutralinoZ1 = squarkRamplitudedecayneutralino (md(2,2), ms, mneut(1), g, gp, mixNeut, 1, -1);
     sstrangeRamplitudestrangeneutralinoZ2 = squarkRamplitudedecayneutralino (md(2,2), ms, mneut(2), g, gp, mixNeut, 2, -1);
     sstrangeRamplitudestrangeneutralinoZ3 = squarkRamplitudedecayneutralino (md(2,2), ms, mneut(3), g, gp, mixNeut, 3, -1);
     sstrangeRamplitudestrangeneutralinoZ4 = squarkRamplitudedecayneutralino (md(2,2), ms, mneut(4), g, gp, mixNeut, 4, -1);
     
     sstrangeRamplitudestrangegravitino = squarkamplitudedecaygravitino(md(2,2), mgravitino, ms, MPlreduced, gravonoff, downsquNLSP);
   }
   else if (nmssmIsIt == true) {
     sstrangeRamplitudestrangeneutralinoZ1 = squarkamplitudedecayquarkneutralinoNMSSM (md(2,2), ms, mneut(1), g, gp, mixNeut, 'd', 'R', 1); 
     sstrangeRamplitudestrangeneutralinoZ2 = squarkamplitudedecayquarkneutralinoNMSSM (md(2,2), ms, mneut(2), g, gp, mixNeut, 'd', 'R', 2);
     sstrangeRamplitudestrangeneutralinoZ3 = squarkamplitudedecayquarkneutralinoNMSSM (md(2,2), ms, mneut(3), g, gp, mixNeut, 'd', 'R', 3);
     sstrangeRamplitudestrangeneutralinoZ4 = squarkamplitudedecayquarkneutralinoNMSSM (md(2,2), ms, mneut(4), g, gp, mixNeut, 'd', 'R', 4);
     sstrangeRamplitudestrangeneutralinoZ5 = squarkamplitudedecayquarkneutralinoNMSSM (md(2,2), ms, mneut(5), g, gp, mixNeut, 'd', 'R', 5);
   }
   
   ParticleSstrangeR.Array_Decays[0][0] = PDGstrange; ParticleSstrangeR.Array_Decays[0][1] = PDGgluino; ParticleSstrangeR.Array_Decays[0][2] = sstrangeRamplitudegluinostrange; ParticleSstrangeR.Array_Decays[0][3] = 2; ParticleSstrangeR.Array_Comments[0] = "# ~s_R -> s ~g";
   ParticleSstrangeR.Array_Decays[1][0] = PDGstrange; ParticleSstrangeR.Array_Decays[1][1] = PDGneutralino1; ParticleSstrangeR.Array_Decays[1][2] = sstrangeRamplitudestrangeneutralinoZ1; ParticleSstrangeR.Array_Decays[1][3] = 2; ParticleSstrangeR.Array_Comments[1] = "# ~s_R -> s ~chi_10";
   ParticleSstrangeR.Array_Decays[2][0] = PDGstrange; ParticleSstrangeR.Array_Decays[2][1] = PDGneutralino2; ParticleSstrangeR.Array_Decays[2][2] = sstrangeRamplitudestrangeneutralinoZ2; ParticleSstrangeR.Array_Decays[2][3] = 2; ParticleSstrangeR.Array_Comments[2] = "# ~s_R -> s ~chi_20";
   ParticleSstrangeR.Array_Decays[3][0] = PDGstrange; ParticleSstrangeR.Array_Decays[3][1] = PDGneutralino3; ParticleSstrangeR.Array_Decays[3][2] = sstrangeRamplitudestrangeneutralinoZ3; ParticleSstrangeR.Array_Decays[3][3] = 2; ParticleSstrangeR.Array_Comments[3] = "# ~s_R -> s ~chi_30";
   ParticleSstrangeR.Array_Decays[4][0] = PDGstrange; ParticleSstrangeR.Array_Decays[4][1] = PDGneutralino4; ParticleSstrangeR.Array_Decays[4][2] = sstrangeRamplitudestrangeneutralinoZ4; ParticleSstrangeR.Array_Decays[4][3] = 2; ParticleSstrangeR.Array_Comments[4] = "# ~s_R -> s ~chi_40";
   ParticleSstrangeR.Array_Decays[5][0] = PDGstrange; ParticleSstrangeR.Array_Decays[5][1] = PDGneutralino5; ParticleSstrangeR.Array_Decays[5][2] = sstrangeRamplitudestrangeneutralinoZ5; ParticleSstrangeR.Array_Decays[5][3] = 2; ParticleSstrangeR.Array_Comments[5] = "# ~s_R -> s ~chi_50";

   ParticleSstrangeR.Array_Decays[6][0] = PDGstrange; ParticleSstrangeR.Array_Decays[6][1] = PDGgravitino; ParticleSstrangeR.Array_Decays[6][2] = sstrangeRamplitudestrangegravitino; ParticleSstrangeR.Array_Decays[6][3] = 2; ParticleSstrangeR.Array_Comments[6] = "# ~s_R -> s ~G";

   for(int i = 0; i<ParticleSstrangeR.No_of_Decays; i++) {
     if (ParticleSstrangeR.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSstrangeR.Array_Comments[i] << " is negative = " << ParticleSstrangeR.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSstrangeR.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }
   
   double SstrangeR_No_1to2_Decays = 0;

   SstrangeR_No_1to2_Decays = ParticleSstrangeR.No_1to2_Decays + ParticleSstrangeR.No_grav_Decays + ParticleSstrangeR.No_NMSSM_Decays;
 
   for (int j = 0; j<SstrangeR_No_1to2_Decays; j++) {
     ParticleSstrangeR.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<SstrangeR_No_1to2_Decays; j++) {
     ParticleSstrangeR.two_width = ParticleSstrangeR.two_width + ParticleSstrangeR.Array_Decays[j][2];
   }
   for (int j=SstrangeR_No_1to2_Decays; j<ParticleSstrangeR.No_of_Decays; j++) {
     ParticleSstrangeR.three_width = ParticleSstrangeR.three_width + ParticleSstrangeR.Array_Decays[j][2];
   }
   ///Note currently no squark three body decays included - may change in future versions
   if ( ParticleSstrangeR.three_width != ParticleSstrangeR.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for sstrangeR - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSstrangeR.No_of_Decays = SstrangeR_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSstrangeR.total_width = ParticleSstrangeR.two_width;
     }
   else {
     ParticleSstrangeR.total_width = ParticleSstrangeR.two_width + ParticleSstrangeR.three_width;
   }

   if ( ParticleSstrangeR.total_width != ParticleSstrangeR.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSstrangeR.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSstrangeR.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in sstrangeR total width \n");
     }
    
 }

 ///ScharmL Decays
 
 double scharmLamplitudegluinocharm=0, scharmLamplitudecharmneutralinoZ1=0, scharmLamplitudecharmneutralinoZ2=0, scharmLamplitudecharmneutralinoZ3=0, scharmLamplitudecharmneutralinoZ4=0, scharmLamplitudecharginoW1strange=0, scharmLamplitudecharginoW2strange=0, scharmLamplitudecharmgravitino=0;
 double scharmLamplitudecharmneutralinoZ5 = 0;

 if (flagscharmL == 1) { 
   scharmLamplitudegluinocharm = squarkamplitudedecaygluino (mu(1,2), mc, mGluino, alphas);
   scharmLamplitudecharginoW1strange = squarkamplitudedecaycharginoW1 (mu(1,2), ms, MCH1, g, thetaR2);
   scharmLamplitudecharginoW2strange = squarkamplitudedecaycharginoW2 (mu(1,2), ms, MCH2, g, thetaR2);
   if (nmssmIsIt == false) {
     scharmLamplitudecharmneutralinoZ1 = squarkLamplitudedecayneutralino (mu(1,2), mc, mneut(1), g, gp, mixNeut, 1, 1);
     scharmLamplitudecharmneutralinoZ2 = squarkLamplitudedecayneutralino (mu(1,2), mc, mneut(2), g, gp, mixNeut, 2, 1);
     scharmLamplitudecharmneutralinoZ3 = squarkLamplitudedecayneutralino (mu(1,2), mc, mneut(3), g, gp, mixNeut, 3, 1);
     scharmLamplitudecharmneutralinoZ4 = squarkLamplitudedecayneutralino (mu(1,2), mc, mneut(4), g, gp, mixNeut, 4, 1);
     
     scharmLamplitudecharmgravitino = squarkamplitudedecaygravitino(mu(1,2), mgravitino, mc, MPlreduced, gravonoff, upsquNLSP);
   }
   else if (nmssmIsIt == true) {
     scharmLamplitudecharmneutralinoZ1 = squarkamplitudedecayquarkneutralinoNMSSM (mu(1,2), mc, mneut(1), g, gp, mixNeut, 'u', 'L', 1); 
     scharmLamplitudecharmneutralinoZ2 = squarkamplitudedecayquarkneutralinoNMSSM (mu(1,2), mc, mneut(2), g, gp, mixNeut, 'u', 'L', 2);
     scharmLamplitudecharmneutralinoZ3 = squarkamplitudedecayquarkneutralinoNMSSM (mu(1,2), mc, mneut(3), g, gp, mixNeut, 'u', 'L', 3);
     scharmLamplitudecharmneutralinoZ4 = squarkamplitudedecayquarkneutralinoNMSSM (mu(1,2), mc, mneut(4), g, gp, mixNeut, 'u', 'L', 4);
     scharmLamplitudecharmneutralinoZ5 = squarkamplitudedecayquarkneutralinoNMSSM (mu(1,2), mc, mneut(5), g, gp, mixNeut, 'u', 'L', 5);
   }
 
   ParticleScharmL.Array_Decays[0][0] = PDGcharm; ParticleScharmL.Array_Decays[0][1] = PDGgluino; ParticleScharmL.Array_Decays[0][2] = scharmLamplitudegluinocharm; ParticleScharmL.Array_Decays[0][3] = 2; ParticleScharmL.Array_Comments[0] = "# ~c_L -> c ~g";
   ParticleScharmL.Array_Decays[1][0] = PDGstrange; ParticleScharmL.Array_Decays[1][1] = PDGchargino1; ParticleScharmL.Array_Decays[1][2] = scharmLamplitudecharginoW1strange; ParticleScharmL.Array_Decays[1][3] = 2; ParticleScharmL.Array_Comments[1] = "# ~c_L -> s ~chi_1+";
   ParticleScharmL.Array_Decays[2][0] = PDGstrange; ParticleScharmL.Array_Decays[2][1] = PDGchargino2; ParticleScharmL.Array_Decays[2][2] = scharmLamplitudecharginoW2strange; ParticleScharmL.Array_Decays[2][3] = 2; ParticleScharmL.Array_Comments[2] = "# ~c_L -> s ~chi_2+";
   ParticleScharmL.Array_Decays[3][0] = PDGcharm; ParticleScharmL.Array_Decays[3][1] = PDGneutralino1; ParticleScharmL.Array_Decays[3][2] = scharmLamplitudecharmneutralinoZ1; ParticleScharmL.Array_Decays[3][3] = 2; ParticleScharmL.Array_Comments[3] = "# ~c_L -> c ~chi_10";
   ParticleScharmL.Array_Decays[4][0] = PDGcharm; ParticleScharmL.Array_Decays[4][1] = PDGneutralino2; ParticleScharmL.Array_Decays[4][2] = scharmLamplitudecharmneutralinoZ2; ParticleScharmL.Array_Decays[4][3] = 2; ParticleScharmL.Array_Comments[4] = "# ~c_L -> c ~chi_20";
   ParticleScharmL.Array_Decays[5][0] = PDGcharm; ParticleScharmL.Array_Decays[5][1] = PDGneutralino3; ParticleScharmL.Array_Decays[5][2] = scharmLamplitudecharmneutralinoZ3; ParticleScharmL.Array_Decays[5][3] = 2; ParticleScharmL.Array_Comments[5] = "# ~c_L -> c ~chi_30";
   ParticleScharmL.Array_Decays[6][0] = PDGcharm; ParticleScharmL.Array_Decays[6][1] = PDGneutralino4; ParticleScharmL.Array_Decays[6][2] = scharmLamplitudecharmneutralinoZ4; ParticleScharmL.Array_Decays[6][3] = 2; ParticleScharmL.Array_Comments[6] = "# ~c_L -> c ~chi_40";
   ParticleScharmL.Array_Decays[7][0] = PDGcharm; ParticleScharmL.Array_Decays[7][1] = PDGneutralino5; ParticleScharmL.Array_Decays[7][2] = scharmLamplitudecharmneutralinoZ5; ParticleScharmL.Array_Decays[7][3] = 2; ParticleScharmL.Array_Comments[7] = "# ~c_L -> c ~chi_50";
   
   ParticleScharmL.Array_Decays[8][0] = PDGcharm; ParticleScharmL.Array_Decays[8][1] = PDGgravitino; ParticleScharmL.Array_Decays[8][2] = scharmLamplitudecharmgravitino; ParticleScharmL.Array_Decays[8][3] = 2; ParticleScharmL.Array_Comments[8] = "# ~c_L -> c ~G";

   for(int i = 0; i<ParticleScharmL.No_of_Decays; i++) {
     if (ParticleScharmL.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleScharmL.Array_Comments[i] << " is negative = " << ParticleScharmL.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleScharmL.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   

   double ScharmL_No_1to2_Decays = 0;

   ScharmL_No_1to2_Decays = ParticleScharmL.No_1to2_Decays + ParticleScharmL.No_grav_Decays + ParticleScharmL.No_NMSSM_Decays;
 
   for (int j = 0; j<ScharmL_No_1to2_Decays; j++) {
     ParticleScharmL.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<ScharmL_No_1to2_Decays; j++) {
     ParticleScharmL.two_width = ParticleScharmL.two_width + ParticleScharmL.Array_Decays[j][2];
   }
   for (int j=ScharmL_No_1to2_Decays; j<ParticleScharmL.No_of_Decays; j++) {
     ParticleScharmL.three_width = ParticleScharmL.three_width + ParticleScharmL.Array_Decays[j][2];
   }
   ///Note currently no squark three body decays included - may change in future versions
   if ( ParticleScharmL.three_width != ParticleScharmL.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for scharmL - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleScharmL.No_of_Decays = ScharmL_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleScharmL.total_width = ParticleScharmL.two_width;
     }
   else {
     ParticleScharmL.total_width = ParticleScharmL.two_width + ParticleScharmL.three_width;
   }

   if ( ParticleScharmL.total_width != ParticleScharmL.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleScharmL.No_of_Decays; i++) {
       //   fout << i << " " << ParticleScharmL.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in scharmL total width \n");
     }
 
 }
 
 
 ///ScharmR decays
 
 double scharmRamplitudegluinocharm=0, scharmRamplitudecharmneutralinoZ1=0, scharmRamplitudecharmneutralinoZ2=0, scharmRamplitudecharmneutralinoZ3=0, scharmRamplitudecharmneutralinoZ4=0, scharmRamplitudecharmgravitino=0;
 double scharmRamplitudecharmneutralinoZ5 = 0;

 if (flagscharmR == 1) {
   scharmRamplitudegluinocharm = squarkamplitudedecaygluino (mu(2,2), mc, mGluino, alphas);
   if (nmssmIsIt == false){
     scharmRamplitudecharmneutralinoZ1 = squarkRamplitudedecayneutralino (mu(2,2), mc, mneut(1), g, gp, mixNeut, 1, 1);
     scharmRamplitudecharmneutralinoZ2 = squarkRamplitudedecayneutralino (mu(2,2), mc, mneut(2), g, gp, mixNeut, 2, 1);
     scharmRamplitudecharmneutralinoZ3 = squarkRamplitudedecayneutralino (mu(2,2), mc, mneut(3), g, gp, mixNeut, 3, 1);
     scharmRamplitudecharmneutralinoZ4 = squarkRamplitudedecayneutralino (mu(2,2), mc, mneut(4), g, gp, mixNeut, 4, 1);
     
     scharmRamplitudecharmgravitino = squarkamplitudedecaygravitino(mu(2,2), mgravitino, mc, MPlreduced,gravonoff, upsquNLSP);
   }
   else if (nmssmIsIt == true) {
     scharmRamplitudecharmneutralinoZ1 = squarkamplitudedecayquarkneutralinoNMSSM (mu(2,2), mc, mneut(1), g, gp, mixNeut, 'u', 'R', 1); 
     scharmRamplitudecharmneutralinoZ2 = squarkamplitudedecayquarkneutralinoNMSSM (mu(2,2), mc, mneut(2), g, gp, mixNeut, 'u', 'R', 2);
     scharmRamplitudecharmneutralinoZ3 = squarkamplitudedecayquarkneutralinoNMSSM (mu(2,2), mc, mneut(3), g, gp, mixNeut, 'u', 'R', 3);
     scharmRamplitudecharmneutralinoZ4 = squarkamplitudedecayquarkneutralinoNMSSM (mu(2,2), mc, mneut(4), g, gp, mixNeut, 'u', 'R', 4);
     scharmRamplitudecharmneutralinoZ5 = squarkamplitudedecayquarkneutralinoNMSSM (mu(2,2), mc, mneut(5), g, gp, mixNeut, 'u', 'R', 5);
   }
   
   ParticleScharmR.Array_Decays[0][0] = PDGcharm; ParticleScharmR.Array_Decays[0][1] = PDGgluino; ParticleScharmR.Array_Decays[0][2] = scharmRamplitudegluinocharm; ParticleScharmR.Array_Decays[0][3] = 2; ParticleScharmR.Array_Comments[0] = "# ~c_R -> c ~g";
   ParticleScharmR.Array_Decays[1][0] = PDGcharm; ParticleScharmR.Array_Decays[1][1] = PDGneutralino1; ParticleScharmR.Array_Decays[1][2] = scharmRamplitudecharmneutralinoZ1; ParticleScharmR.Array_Decays[1][3] = 2; ParticleScharmR.Array_Comments[1] = "# ~c_R -> c ~chi_10";
   ParticleScharmR.Array_Decays[2][0] = PDGcharm; ParticleScharmR.Array_Decays[2][1] = PDGneutralino2; ParticleScharmR.Array_Decays[2][2] = scharmRamplitudecharmneutralinoZ2; ParticleScharmR.Array_Decays[2][3] = 2; ParticleScharmR.Array_Comments[2] = "# ~c_R -> c ~chi_20";
   ParticleScharmR.Array_Decays[3][0] = PDGcharm; ParticleScharmR.Array_Decays[3][1] = PDGneutralino3; ParticleScharmR.Array_Decays[3][2] = scharmRamplitudecharmneutralinoZ3; ParticleScharmR.Array_Decays[3][3] = 2; ParticleScharmR.Array_Comments[3] = "# ~c_R -> c ~chi_30";
   ParticleScharmR.Array_Decays[4][0] = PDGcharm; ParticleScharmR.Array_Decays[4][1] = PDGneutralino4; ParticleScharmR.Array_Decays[4][2] = scharmRamplitudecharmneutralinoZ4; ParticleScharmR.Array_Decays[4][3] = 2; ParticleScharmR.Array_Comments[4] = "# ~c_R -> c ~chi_40";
   ParticleScharmR.Array_Decays[5][0] = PDGcharm; ParticleScharmR.Array_Decays[5][1] = PDGneutralino5; ParticleScharmR.Array_Decays[5][2] = scharmRamplitudecharmneutralinoZ5; ParticleScharmR.Array_Decays[5][3] = 2; ParticleScharmR.Array_Comments[5] = "# ~c_R -> c ~chi_50";
   
   ParticleScharmR.Array_Decays[6][0] = PDGcharm; ParticleScharmR.Array_Decays[6][1] = PDGgravitino; ParticleScharmR.Array_Decays[6][2] = scharmRamplitudecharmgravitino; ParticleScharmR.Array_Decays[6][3] = 2; ParticleScharmR.Array_Comments[6] = "# ~c_R -> c ~G";

   for(int i = 0; i<ParticleScharmR.No_of_Decays; i++) {
     if (ParticleScharmR.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleScharmR.Array_Comments[i] << " is negative = " << ParticleScharmR.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleScharmR.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }
   

   double ScharmR_No_1to2_Decays = 0;

   ScharmR_No_1to2_Decays = ParticleScharmR.No_1to2_Decays + ParticleScharmR.No_grav_Decays + ParticleScharmR.No_NMSSM_Decays;
 
   for (int j = 0; j<ScharmR_No_1to2_Decays; j++) {
     ParticleScharmR.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<ScharmR_No_1to2_Decays; j++) {
     ParticleScharmR.two_width = ParticleScharmR.two_width + ParticleScharmR.Array_Decays[j][2];
   }
   for (int j=ScharmR_No_1to2_Decays; j<ParticleScharmR.No_of_Decays; j++) {
     ParticleScharmR.three_width = ParticleScharmR.three_width + ParticleScharmR.Array_Decays[j][2];
   }
   ///Note currently no squark three body decays included - may change in future versions 
   if ( ParticleScharmR.three_width != ParticleScharmR.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for scharmR - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleScharmR.No_of_Decays = ScharmR_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleScharmR.total_width = ParticleScharmR.two_width;
     }
   else {
     ParticleScharmR.total_width = ParticleScharmR.two_width + ParticleScharmL.three_width;
   }
   
   if ( ParticleScharmR.total_width != ParticleScharmR.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleScharmR.No_of_Decays; i++) {
       //   fout << i << " " << ParticleScharmR.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in scharmR total width \n");
     }
 }
 
 ///Sbottom1 Decays
 double sbottom1amplitudegluinobottom=0, sbottom1amplitudebottomneutralinoZ1=0, sbottom1amplitudebottomneutralinoZ2=0, sbottom1amplitudebottomneutralinoZ3=0, sbottom1amplitudebottomneutralinoZ4=0, sbottom1amplitudetopcharginoW1=0, sbottom1amplitudetopcharginoW2=0, sbottom1amplitudeWbosonstop1=0, sbottom1amplitudeWbosonstop2=0, sbottom1amplitudeHminusstop1=0, sbottom1amplitudeHminusstop2=0, sbottom1amplitudebottomgravitino=0;
 double sbottom1amplitudebottomneutralinoZ5 = 0;

 if (flagsbottom1 == 1) {
   sbottom1amplitudegluinobottom = squarkamplitudedecaygluinomix (md(1,3), mb, mGluino, alphas, 1, thetab);
   sbottom1amplitudetopcharginoW1 = squark1amplitudedecaycharginoW1mix (md(1,3), mt, MCH1, g, thetaL2, thetaR2, thetab, beta, runmw, runmt, runmb, 2);
   sbottom1amplitudetopcharginoW2 = squark1amplitudedecaycharginoW2mix (md(1,3), mt, MCH2, g, thetaL2, thetaR2, thetab, beta, runmw, runmt, runmb, 2);
   if (nmssmIsIt == false) {
     sbottom1amplitudebottomneutralinoZ1 = squark3amplitudedecayneutralino (md(1,3), mb, mneut(1), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 1);
     sbottom1amplitudebottomneutralinoZ2 = squark3amplitudedecayneutralino (md(1,3), mb, mneut(2), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 2);
     sbottom1amplitudebottomneutralinoZ3 = squark3amplitudedecayneutralino (md(1,3), mb, mneut(3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 3);
     sbottom1amplitudebottomneutralinoZ4 = squark3amplitudedecayneutralino (md(1,3), mb, mneut(4), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 4);
     sbottom1amplitudebottomgravitino = squarkamplitudedecaygravitino(md(1,3), mgravitino, mb, MPlreduced, gravonoff, downsquNLSP);
   }
   else if (nmssmIsIt == true) {
     sbottom1amplitudebottomneutralinoZ1 = sbottomamplitudedecaybottomneutralinoNMSSM (md(1,3), mb, mneut(1), g, gp, thetab, mixNeut, runmb, runmw, beta, 1, 1);
     sbottom1amplitudebottomneutralinoZ2 = sbottomamplitudedecaybottomneutralinoNMSSM (md(1,3), mb, mneut(2), g, gp, thetab, mixNeut, runmb, runmw, beta, 1, 2);
     sbottom1amplitudebottomneutralinoZ3 = sbottomamplitudedecaybottomneutralinoNMSSM (md(1,3), mb, mneut(3), g, gp, thetab, mixNeut, runmb, runmw, beta, 1, 3);
     sbottom1amplitudebottomneutralinoZ4 = sbottomamplitudedecaybottomneutralinoNMSSM (md(1,3), mb, mneut(4), g, gp, thetab, mixNeut, runmb, runmw, beta, 1, 4);
     sbottom1amplitudebottomneutralinoZ5 = sbottomamplitudedecaybottomneutralinoNMSSM (md(1,3), mb, mneut(5), g, gp, thetab, mixNeut, runmb, runmw, beta, 1, 5);
   }

   sbottom1amplitudeWbosonstop1 = squark3amplitudedecaysquark3Wboson (md(1,3), polemw, mu(1,3), g, thetat, thetab, 2, 1, 1, 1);
   sbottom1amplitudeWbosonstop2 = squark3amplitudedecaysquark3Wboson (md(1,3), polemw, mu(2,3), g, thetat, thetab, 2, 1, 1, 2);
   sbottom1amplitudeHminusstop1 = squark3amplitudedecaychargedHiggssquark3 (md(1,3), mHpm, mu(1,3), g, runmw, beta, thetat, thetab, greekmu, At, Ab, runmt, runmb, 1, 1);
   sbottom1amplitudeHminusstop2 = squark3amplitudedecaychargedHiggssquark3 (md(1,3), mHpm, mu(2,3), g, runmw, beta, thetat, thetab, greekmu, At, Ab, runmt, runmb, 2, 1);

   ParticleSbottom1.Array_Decays[0][0] = PDGbottom; ParticleSbottom1.Array_Decays[0][1] = PDGgluino; ParticleSbottom1.Array_Decays[0][2] = sbottom1amplitudegluinobottom; ParticleSbottom1.Array_Decays[0][3] = 2; ParticleSbottom1.Array_Comments[0] = "# ~b_1 -> b ~g";
   ParticleSbottom1.Array_Decays[1][0] = PDGtop; ParticleSbottom1.Array_Decays[1][1] = -PDGchargino1; ParticleSbottom1.Array_Decays[1][2] = sbottom1amplitudetopcharginoW1; ParticleSbottom1.Array_Decays[1][3] = 2; ParticleSbottom1.Array_Comments[1] = "# ~b_1 -> t ~chi_1-"; 
   ParticleSbottom1.Array_Decays[2][0] = PDGtop; ParticleSbottom1.Array_Decays[2][1] = -PDGchargino2; ParticleSbottom1.Array_Decays[2][2] = sbottom1amplitudetopcharginoW2; ParticleSbottom1.Array_Decays[2][3] = 2; ParticleSbottom1.Array_Comments[2] = "# ~b_1 -> t ~chi_2-";

   ParticleSbottom1.Array_Decays[3][0] = -PDGWplus; ParticleSbottom1.Array_Decays[3][1] = PDGstop1 ; ParticleSbottom1.Array_Decays[3][2] = sbottom1amplitudeWbosonstop1; ParticleSbottom1.Array_Decays[3][3] = 2; ParticleSbottom1.Array_Comments[3] = "# ~b_1 -> W- ~t_1";
   ParticleSbottom1.Array_Decays[4][0] = -PDGWplus; ParticleSbottom1.Array_Decays[4][1] = PDGstop2 ; ParticleSbottom1.Array_Decays[4][2] = sbottom1amplitudeWbosonstop2; ParticleSbottom1.Array_Decays[4][3] = 2; ParticleSbottom1.Array_Comments[4] = "# ~b_1 -> W- ~t_2";
   ParticleSbottom1.Array_Decays[5][0] = -PDGHplus; ParticleSbottom1.Array_Decays[5][1] = PDGstop1; ParticleSbottom1.Array_Decays[5][2] = sbottom1amplitudeHminusstop1; ParticleSbottom1.Array_Decays[5][3] = 2; ParticleSbottom1.Array_Comments[5] = "# ~b_1 -> H- ~t_1";
   ParticleSbottom1.Array_Decays[6][0] = -PDGHplus; ParticleSbottom1.Array_Decays[6][1] = PDGstop2; ParticleSbottom1.Array_Decays[6][2] = sbottom1amplitudeHminusstop2; ParticleSbottom1.Array_Decays[6][3] = 2; ParticleSbottom1.Array_Comments[6] = "# ~b_1 -> H- ~t_2";

   ParticleSbottom1.Array_Decays[7][0] = PDGbottom; ParticleSbottom1.Array_Decays[7][1] = PDGneutralino1; ParticleSbottom1.Array_Decays[7][2] = sbottom1amplitudebottomneutralinoZ1; ParticleSbottom1.Array_Decays[7][3] = 2; ParticleSbottom1.Array_Comments[7] = "# ~b_1 -> b ~chi_10";
   ParticleSbottom1.Array_Decays[8][0] = PDGbottom; ParticleSbottom1.Array_Decays[8][1] = PDGneutralino2; ParticleSbottom1.Array_Decays[8][2] = sbottom1amplitudebottomneutralinoZ2; ParticleSbottom1.Array_Decays[8][3] = 2; ParticleSbottom1.Array_Comments[8] = "# ~b_1 -> b ~chi_20";
   ParticleSbottom1.Array_Decays[9][0] = PDGbottom; ParticleSbottom1.Array_Decays[9][1] = PDGneutralino3; ParticleSbottom1.Array_Decays[9][2] = sbottom1amplitudebottomneutralinoZ3; ParticleSbottom1.Array_Decays[9][3] = 2; ParticleSbottom1.Array_Comments[9] = "# ~b_1 -> b ~chi_30";
   ParticleSbottom1.Array_Decays[10][0] = PDGbottom; ParticleSbottom1.Array_Decays[10][1] = PDGneutralino4; ParticleSbottom1.Array_Decays[10][2] = sbottom1amplitudebottomneutralinoZ4; ParticleSbottom1.Array_Decays[10][3] = 2; ParticleSbottom1.Array_Comments[10] = "# ~b_1 -> b ~chi_40";
   ParticleSbottom1.Array_Decays[11][0] = PDGbottom; ParticleSbottom1.Array_Decays[11][1] = PDGneutralino5; ParticleSbottom1.Array_Decays[11][2] = sbottom1amplitudebottomneutralinoZ5; ParticleSbottom1.Array_Decays[11][3] = 2; ParticleSbottom1.Array_Comments[11] = "# ~b_1 -> b ~chi_50";

   ParticleSbottom1.Array_Decays[12][0] = PDGbottom; ParticleSbottom1.Array_Decays[12][1] = PDGgravitino; ParticleSbottom1.Array_Decays[12][2] = sbottom1amplitudebottomgravitino; ParticleSbottom1.Array_Decays[12][3] = 2; ParticleSbottom1.Array_Comments[12] = "# ~b_1 -> b ~G";

   for(int i = 0; i<ParticleSbottom1.No_of_Decays; i++) {
     if (ParticleSbottom1.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSbottom1.Array_Comments[i] << " is negative = " << ParticleSbottom1.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSbottom1.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }

   double Sbottom1_No_1to2_Decays = 0;

   Sbottom1_No_1to2_Decays = ParticleSbottom1.No_1to2_Decays + ParticleSbottom1.No_grav_Decays + ParticleSbottom1.No_NMSSM_Decays;
 
   for (int j = 0; j<Sbottom1_No_1to2_Decays; j++) {
     ParticleSbottom1.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<Sbottom1_No_1to2_Decays; j++) {
     ParticleSbottom1.two_width = ParticleSbottom1.two_width + ParticleSbottom1.Array_Decays[j][2];
   }
   for (int j=Sbottom1_No_1to2_Decays; j<ParticleSbottom1.No_of_Decays; j++) {
     ParticleSbottom1.three_width = ParticleSbottom1.three_width + ParticleSbottom1.Array_Decays[j][2];
   }
   ///Note currently no squark three body decays included - may change in future versions
   if ( ParticleSbottom1.three_width != ParticleSbottom1.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for sbottom1 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       ParticleSbottom1.No_of_Decays = Sbottom1_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSbottom1.total_width = ParticleSbottom1.two_width;
     }
   else {
     ParticleSbottom1.total_width = ParticleSbottom1.two_width + ParticleSbottom1.three_width;
   }

   if ( ParticleSbottom1.total_width != ParticleSbottom1.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSbottom1.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSbottom1.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in sbottom1 total width \n");
     }
 }

 ///Sbottom2 Decays

 double sbottom2amplitudegluinobottom=0, sbottom2amplitudebottomneutralinoZ1=0, sbottom2amplitudebottomneutralinoZ2=0, sbottom2amplitudebottomneutralinoZ3=0, sbottom2amplitudebottomneutralinoZ4=0, sbottom2amplitudetopcharginoW1=0, sbottom2amplitudetopcharginoW2=0, sbottom2amplitudeWbosonstop1=0, sbottom2amplitudeWbosonstop2=0, sbottom2amplitudeHminusstop1=0, sbottom2amplitudeHminusstop2=0, sbottom2amplitudeZbosonsbottom1=0, sbottom2amplitudehsbottom1=0, sbottom2amplitudeAsbottom1=0, sbottom2amplitudeHsbottom1=0, sbottom2amplitudebottomgravitino=0;

 double sbottom2amplitudeH3sbottom1 = 0, sbottom2amplitudeA2sbottom1 = 0, sbottom2amplitudebottomneutralinoZ5 = 0;
 
 if (flagsbottom2 == 1) {
   sbottom2amplitudegluinobottom =  squarkamplitudedecaygluinomix (md(2,3), mb, mGluino, alphas, 2, thetab);
   sbottom2amplitudetopcharginoW1 = squark2amplitudedecaycharginoW1mix (md(2,3), mt, MCH1, g, thetaL2, thetaR2, thetab, beta, runmw, runmt, runmb,2);
   sbottom2amplitudetopcharginoW2 = squark2amplitudedecaycharginoW2mix (md(2,3), mt, MCH2, g, thetaL2, thetaR2, thetab, beta, runmw, runmt, runmb,2); 
   sbottom2amplitudebottomneutralinoZ1 = squark3amplitudedecayneutralino (md(2,3), mb, mneut(1), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 1);
   sbottom2amplitudebottomneutralinoZ2 = squark3amplitudedecayneutralino (md(2,3), mb, mneut(2), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 2);
   sbottom2amplitudebottomneutralinoZ3 = squark3amplitudedecayneutralino (md(2,3), mb, mneut(3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 3);
   sbottom2amplitudebottomneutralinoZ4 = squark3amplitudedecayneutralino (md(2,3), mb, mneut(4), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 4);
   sbottom2amplitudeWbosonstop1 = squark3amplitudedecaysquark3Wboson (md(2,3), polemw, mu(1,3), g, thetat, thetab, 2, 2, 1, 1);
   sbottom2amplitudeWbosonstop2 = squark3amplitudedecaysquark3Wboson (md(2,3), polemw, mu(2,3), g, thetat, thetab, 2, 2, 1, 2);
   sbottom2amplitudeHminusstop1 = squark3amplitudedecaychargedHiggssquark3 (md(2,3), mHpm, mu(1,3), g, runmw, beta, thetat, thetab, greekmu, At, Ab, runmt, runmb, 1, 2);
   sbottom2amplitudeHminusstop2 = squark3amplitudedecaychargedHiggssquark3 (md(2,3), mHpm, mu(2,3), g, runmw, beta, thetat, thetab, greekmu, At, Ab, runmt, runmb, 2, 2);
   sbottom2amplitudeZbosonsbottom1 = squark32amplitudedecaysquark3Zboson (md(2,3), md(1,3), polemz, g, gp, thetab);
   if (nmssmIsIt == false) {
     sbottom2amplitudehsbottom1 = squark32amplitudedecayneutralHiggssquark3 (md(2,3), mh0(1), md(1,3), g, gp, runmw, beta, alpha, thetat, thetab, greekmu, At, Ab, runmt, runmb, 2, 'h');
     sbottom2amplitudeHsbottom1 = squark32amplitudedecayneutralHiggssquark3 (md(2,3), mh0(2), md(1,3), g, gp, runmw, beta, alpha, thetat, thetab, greekmu, At, Ab, runmt, runmb, 2, 'H');
     sbottom2amplitudeAsbottom1 = squark32amplitudedecayneutralHiggssquark3 (md(2,3), mA0(1), md(1,3), g, gp, runmw, beta, alpha, thetat, thetab, greekmu, At, Ab, runmt, runmb, 2, 'A');
     sbottom2amplitudebottomgravitino = squarkamplitudedecaygravitino(md(2,3), mgravitino, mb, MPlreduced, gravonoff, downsquNLSP);
   }
   else if (nmssmIsIt == true) {
     sbottom2amplitudehsbottom1 = sbottom2amplitudedecaysbottom1CPevenhiggsNMSSM (md(2,3), md(1,3), mh0(1), runmb, thetab, CPEMix, beta, runmw, g, gp, Ab, mueff, lam, 1);
     sbottom2amplitudeHsbottom1 = sbottom2amplitudedecaysbottom1CPevenhiggsNMSSM (md(2,3), md(1,3), mh0(2), runmb, thetab, CPEMix, beta, runmw, g, gp, Ab, mueff, lam, 2);
     sbottom2amplitudeH3sbottom1 = sbottom2amplitudedecaysbottom1CPevenhiggsNMSSM (md(2,3), md(1,3), mh0(3), runmb, thetab, CPEMix, beta, runmw, g, gp, Ab, mueff, lam, 3);
     sbottom2amplitudeAsbottom1 = sbottom2amplitudedecaysbottom1CPoddhiggsNMSSM (md(2,3), md(1,3), mA0(1), runmb , thetab, CPOMix, beta, runmw, g, Ab, mueff, lam, 1); 
     sbottom2amplitudeA2sbottom1 = sbottom2amplitudedecaysbottom1CPoddhiggsNMSSM (md(2,3), md(1,3), mA0(2), runmb , thetab, CPOMix, beta, runmw, g, Ab, mueff, lam, 2);
     sbottom2amplitudebottomneutralinoZ1 = sbottomamplitudedecaybottomneutralinoNMSSM (md(2,3), mb, mneut(1), g, gp, thetab, mixNeut, runmb, runmw, beta, 2, 1);
     sbottom2amplitudebottomneutralinoZ2 = sbottomamplitudedecaybottomneutralinoNMSSM (md(2,3), mb, mneut(2), g, gp, thetab, mixNeut, runmb, runmw, beta, 2, 2);
     sbottom2amplitudebottomneutralinoZ3 = sbottomamplitudedecaybottomneutralinoNMSSM (md(2,3), mb, mneut(3), g, gp, thetab, mixNeut, runmb, runmw, beta, 2, 3);
     sbottom2amplitudebottomneutralinoZ4 = sbottomamplitudedecaybottomneutralinoNMSSM (md(2,3), mb, mneut(4), g, gp, thetab, mixNeut, runmb, runmw, beta, 2, 4);
     sbottom2amplitudebottomneutralinoZ5 = sbottomamplitudedecaybottomneutralinoNMSSM (md(2,3), mb, mneut(5), g, gp, thetab, mixNeut, runmb, runmw, beta, 2, 5);   
   }
   
   ParticleSbottom2.Array_Decays[0][0] = PDGbottom; ParticleSbottom2.Array_Decays[0][1] = PDGgluino; ParticleSbottom2.Array_Decays[0][2] = sbottom2amplitudegluinobottom; ParticleSbottom2.Array_Decays[0][3] = 2; ParticleSbottom2.Array_Comments[0] = "# ~b_2 -> b ~g";
   ParticleSbottom2.Array_Decays[1][0] = PDGtop; ParticleSbottom2.Array_Decays[1][1] = -PDGchargino1; ParticleSbottom2.Array_Decays[1][2] = sbottom2amplitudetopcharginoW1; ParticleSbottom2.Array_Decays[1][3] = 2; ParticleSbottom2.Array_Comments[1] = "# ~b_2 -> t ~chi_1-";
   ParticleSbottom2.Array_Decays[2][0] = PDGtop; ParticleSbottom2.Array_Decays[2][1] = -PDGchargino2; ParticleSbottom2.Array_Decays[2][2] = sbottom2amplitudetopcharginoW2; ParticleSbottom2.Array_Decays[2][3] = 2; ParticleSbottom2.Array_Comments[2] = "# ~b_2 -> t ~chi_2-";

   ParticleSbottom2.Array_Decays[3][0] = -PDGWplus; ParticleSbottom2.Array_Decays[3][1] = PDGstop1 ; ParticleSbottom2.Array_Decays[3][2] = sbottom2amplitudeWbosonstop1; ParticleSbottom2.Array_Decays[3][3] = 2; ParticleSbottom2.Array_Comments[3] = "# ~b_2 -> W- ~t_1";
   ParticleSbottom2.Array_Decays[4][0] = -PDGWplus; ParticleSbottom2.Array_Decays[4][1] = PDGstop2 ; ParticleSbottom2.Array_Decays[4][2] = sbottom2amplitudeWbosonstop2; ParticleSbottom2.Array_Decays[4][3] = 2; ParticleSbottom2.Array_Comments[4] = "# ~b_2 -> W- ~t_2";
   ParticleSbottom2.Array_Decays[5][0] = -PDGHplus; ParticleSbottom2.Array_Decays[5][1] = PDGstop1; ParticleSbottom2.Array_Decays[5][2] = sbottom2amplitudeHminusstop1; ParticleSbottom2.Array_Decays[5][3] = 2; ParticleSbottom2.Array_Comments[5] = "# ~b_2 -> H- ~t_1";
   ParticleSbottom2.Array_Decays[6][0] = -PDGHplus; ParticleSbottom2.Array_Decays[6][1] = PDGstop2; ParticleSbottom2.Array_Decays[6][2] = sbottom2amplitudeHminusstop2; ParticleSbottom2.Array_Decays[6][3] = 2; ParticleSbottom2.Array_Comments[6] = "# ~b_2 -> H- ~t_2";
   ParticleSbottom2.Array_Decays[7][0] = PDGh0; ParticleSbottom2.Array_Decays[7][1] = PDGsbottom1; ParticleSbottom2.Array_Decays[7][2] = sbottom2amplitudehsbottom1; ParticleSbottom2.Array_Decays[7][3] = 2; ParticleSbottom2.Array_Comments[7] = "# ~b_2 -> h ~b_1";
   ParticleSbottom2.Array_Decays[8][0] = PDGH0; ParticleSbottom2.Array_Decays[8][1] = PDGsbottom1; ParticleSbottom2.Array_Decays[8][2] = sbottom2amplitudeHsbottom1; ParticleSbottom2.Array_Decays[8][3] = 2; ParticleSbottom2.Array_Comments[8] = "# ~b_2 -> H ~b_1";
   ParticleSbottom2.Array_Decays[9][0] = PDGA0; ParticleSbottom2.Array_Decays[9][1] = PDGsbottom1; ParticleSbottom2.Array_Decays[9][2] = sbottom2amplitudeAsbottom1; ParticleSbottom2.Array_Decays[9][3] = 2; ParticleSbottom2.Array_Comments[9] = "# ~b_2 -> A ~b_1";
   ParticleSbottom2.Array_Decays[10][0] = PDGsbottom1; ParticleSbottom2.Array_Decays[10][1] = PDGZboson; ParticleSbottom2.Array_Decays[10][2] = sbottom2amplitudeZbosonsbottom1; ParticleSbottom2.Array_Decays[10][3] = 2; ParticleSbottom2.Array_Comments[10] = "# ~b_2 -> Z ~b_1";

   ParticleSbottom2.Array_Decays[11][0] = PDGbottom; ParticleSbottom2.Array_Decays[11][1] = PDGneutralino1; ParticleSbottom2.Array_Decays[11][2] = sbottom2amplitudebottomneutralinoZ1; ParticleSbottom2.Array_Decays[11][3] = 2; ParticleSbottom2.Array_Comments[11] = "# ~b_2 -> b ~chi_10";
   ParticleSbottom2.Array_Decays[12][0] = PDGbottom; ParticleSbottom2.Array_Decays[12][1] = PDGneutralino2; ParticleSbottom2.Array_Decays[12][2] = sbottom2amplitudebottomneutralinoZ2; ParticleSbottom2.Array_Decays[12][3] = 2; ParticleSbottom2.Array_Comments[12] = "# ~b_2 -> b ~chi_20";
   ParticleSbottom2.Array_Decays[13][0] = PDGbottom; ParticleSbottom2.Array_Decays[13][1] = PDGneutralino3; ParticleSbottom2.Array_Decays[13][2] = sbottom2amplitudebottomneutralinoZ3; ParticleSbottom2.Array_Decays[13][3] = 2; ParticleSbottom2.Array_Comments[13] = "# ~b_2 -> b ~chi_30";
   ParticleSbottom2.Array_Decays[14][0] = PDGbottom; ParticleSbottom2.Array_Decays[14][1] = PDGneutralino4; ParticleSbottom2.Array_Decays[14][2] = sbottom2amplitudebottomneutralinoZ4; ParticleSbottom2.Array_Decays[14][3] = 2; ParticleSbottom2.Array_Comments[14] = "# ~b_2 -> b ~chi_40";
   ParticleSbottom2.Array_Decays[15][0] = PDGbottom; ParticleSbottom2.Array_Decays[15][1] = PDGneutralino5; ParticleSbottom2.Array_Decays[15][2] = sbottom2amplitudebottomneutralinoZ5; ParticleSbottom2.Array_Decays[15][3] = 2; ParticleSbottom2.Array_Comments[15] = "# ~b_2 -> b ~chi_50";

   ParticleSbottom2.Array_Decays[16][0] = PDGbottom; ParticleSbottom2.Array_Decays[16][1] = PDGgravitino; ParticleSbottom2.Array_Decays[16][2] = sbottom2amplitudebottomgravitino; ParticleSbottom2.Array_Decays[16][3] = 2; ParticleSbottom2.Array_Comments[16] = "# ~b_2 -> b ~G";
   ParticleSbottom2.Array_Decays[17][0] = PDGH3; ParticleSbottom2.Array_Decays[17][1] = PDGsbottom1; ParticleSbottom2.Array_Decays[17][2] = sbottom2amplitudeH3sbottom1; ParticleSbottom2.Array_Decays[17][3] = 2; ParticleSbottom2.Array_Comments[17] = "# ~b_2 -> H3 ~b_1";
   ParticleSbottom2.Array_Decays[18][0] = PDGA2; ParticleSbottom2.Array_Decays[18][1] = PDGsbottom1; ParticleSbottom2.Array_Decays[18][2] = sbottom2amplitudeA2sbottom1; ParticleSbottom2.Array_Decays[18][3] = 2; ParticleSbottom2.Array_Comments[18] = "# ~b_2 -> A2 ~b_1";

   for(int i = 0; i<ParticleSbottom2.No_of_Decays; i++) {
     if (ParticleSbottom2.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSbottom2.Array_Comments[i] << " is negative = " << ParticleSbottom2.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSbottom2.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }
   
   double Sbottom2_No_1to2_Decays = 0;

   Sbottom2_No_1to2_Decays = ParticleSbottom2.No_1to2_Decays + ParticleSbottom2.No_grav_Decays + ParticleSbottom2.No_NMSSM_Decays;

   for (int j = 0; j<Sbottom2_No_1to2_Decays; j++) {
     ParticleSbottom2.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<Sbottom2_No_1to2_Decays; j++) {
     ParticleSbottom2.two_width = ParticleSbottom2.two_width + ParticleSbottom2.Array_Decays[j][2];
   }
   for (int j=Sbottom2_No_1to2_Decays; j<ParticleSbottom2.No_of_Decays; j++) {
     ParticleSbottom2.three_width = ParticleSbottom2.three_width + ParticleSbottom2.Array_Decays[j][2];
   }
   ///Note currently no squark three body decays included - may change in future versions
   if ( ParticleSbottom2.three_width != ParticleSbottom2.three_width) /// Tests for a nan as only nans aren't equal to themselves
   {
     fout << "# Three body decays give nan for sbottom2 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
     errorflag = -1;
     ParticleSbottom2.No_of_Decays = Sbottom2_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
     ParticleSbottom2.total_width = ParticleSbottom2.two_width;
   }
   else {
     ParticleSbottom2.total_width = ParticleSbottom2.two_width + ParticleSbottom2.three_width;
   }
   
   if ( ParticleSbottom2.total_width != ParticleSbottom2.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSbottom2.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSbottom2.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in sbottom2 total width \n");
     }
   
 }

 ///Stop1 Decays
 double stop1amplitudegluinotop=0, stop1amplitudetopneutralinoZ1=0, stop1amplitudetopneutralinoZ2=0, stop1amplitudetopneutralinoZ3=0, stop1amplitudetopneutralinoZ4=0, stop1amplitudebottomcharginoW1=0, stop1amplitudebottomcharginoW2=0, stop1amplitudeWbosonsbottom1=0, stop1amplitudeWbosonsbottom2=0, stop1amplitudeHplussbottom1=0, stop1amplitudeHplussbottom2=0, stop1amplitudetopgravitino=0;
 
 double stop1amplitudetopneutralinoZ5 = 0;

 if (flagstop1 == 1) {
   stop1amplitudegluinotop = squarkamplitudedecaygluinomix (mu(1,3), mt, mGluino, alphas, 1, thetat);
   stop1amplitudebottomcharginoW1 = squark1amplitudedecaycharginoW1mix (mu(1,3), mb, MCH1, g, thetaL2, thetaR2, thetat, beta, runmw, runmt, runmb, 1);
   stop1amplitudebottomcharginoW2 = squark1amplitudedecaycharginoW2mix (mu(1,3), mb, MCH2, g, thetaL2, thetaR2, thetat, beta, runmw, runmt, runmb, 1);
   if (nmssmIsIt == false) {
     stop1amplitudetopneutralinoZ1 = squark3amplitudedecayneutralino (mu(1,3), mt, mneut(1), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 1);
     stop1amplitudetopneutralinoZ2 = squark3amplitudedecayneutralino (mu(1,3), mt, mneut(2), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 2);
     stop1amplitudetopneutralinoZ3 = squark3amplitudedecayneutralino (mu(1,3), mt, mneut(3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 3);
     stop1amplitudetopneutralinoZ4 = squark3amplitudedecayneutralino (mu(1,3), mt, mneut(4), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 4);
     stop1amplitudetopgravitino = squarkamplitudedecaygravitino(mu(1,3), mgravitino, mt, MPlreduced, gravonoff, upsquNLSP);
   }
   else if (nmssmIsIt == true) {
     stop1amplitudetopneutralinoZ1 = stopamplitudedecaytopneutralinoNMSSM (mu(1,3), mt, mneut(1), g, gp, thetat, mixNeut, runmt, runmw, beta, 1, 1);
     stop1amplitudetopneutralinoZ2 = stopamplitudedecaytopneutralinoNMSSM (mu(1,3), mt, mneut(2), g, gp, thetat, mixNeut, runmt, runmw, beta, 1, 2);
     stop1amplitudetopneutralinoZ3 = stopamplitudedecaytopneutralinoNMSSM (mu(1,3), mt, mneut(3), g, gp, thetat, mixNeut, runmt, runmw, beta, 1, 3);
     stop1amplitudetopneutralinoZ4 = stopamplitudedecaytopneutralinoNMSSM (mu(1,3), mt, mneut(4), g, gp, thetat, mixNeut, runmt, runmw, beta, 1, 4);
     stop1amplitudetopneutralinoZ5 = stopamplitudedecaytopneutralinoNMSSM (mu(1,3), mt, mneut(5), g, gp, thetat, mixNeut, runmt, runmw, beta, 1, 5);
   }
   
   stop1amplitudeWbosonsbottom1 = squark3amplitudedecaysquark3Wboson (mu(1,3), polemw, md(1,3), g, thetat, thetab, 1, 1, 2, 1);
   stop1amplitudeWbosonsbottom2 = squark3amplitudedecaysquark3Wboson (mu(1,3), polemw, md(2,3), g, thetat, thetab, 1, 1, 2, 2);
   stop1amplitudeHplussbottom1 = squark3amplitudedecaychargedHiggssquark3 (mu(1,3), mHpm, md(1,3), g, runmw, beta, thetat, thetab, greekmu, At, Ab, runmt, runmb, 1, 1);
   stop1amplitudeHplussbottom2 = squark3amplitudedecaychargedHiggssquark3 (mu(1,3), mHpm, md(2,3), g, runmw, beta, thetat, thetab, greekmu, At, Ab, runmt, runmb, 1, 2);
 
   ParticleStop1.Array_Decays[0][0] = PDGtop; ParticleStop1.Array_Decays[0][1] = PDGgluino; ParticleStop1.Array_Decays[0][2] = stop1amplitudegluinotop; ParticleStop1.Array_Decays[0][3] = 2; ParticleStop1.Array_Comments[0] = "# ~t_1 -> t ~g";
   ParticleStop1.Array_Decays[1][0] = PDGbottom; ParticleStop1.Array_Decays[1][1] = PDGchargino1; ParticleStop1.Array_Decays[1][2] = stop1amplitudebottomcharginoW1; ParticleStop1.Array_Decays[1][3] = 2; ParticleStop1.Array_Comments[1] = "# ~t_1 -> b ~chi_1+";
   ParticleStop1.Array_Decays[2][0] = PDGbottom; ParticleStop1.Array_Decays[2][1] = PDGchargino2; ParticleStop1.Array_Decays[2][2] = stop1amplitudebottomcharginoW2; ParticleStop1.Array_Decays[2][3] = 2; ParticleStop1.Array_Comments[2] = "# ~t_1 -> b ~chi_2+";
   
   ParticleStop1.Array_Decays[3][0] = PDGWplus; ParticleStop1.Array_Decays[3][1] = PDGsbottom1 ; ParticleStop1.Array_Decays[3][2] = stop1amplitudeWbosonsbottom1; ParticleStop1.Array_Decays[3][3] = 2; ParticleStop1.Array_Comments[3] = "# ~t_1 -> W+ ~b_1";
   ParticleStop1.Array_Decays[4][0] = PDGWplus; ParticleStop1.Array_Decays[4][1] = PDGsbottom2 ; ParticleStop1.Array_Decays[4][2] = stop1amplitudeWbosonsbottom2; ParticleStop1.Array_Decays[4][3] = 2; ParticleStop1.Array_Comments[4] = "# ~t_1 -> W+ ~b_2";
   ParticleStop1.Array_Decays[5][0] = PDGHplus; ParticleStop1.Array_Decays[5][1] = PDGsbottom1; ParticleStop1.Array_Decays[5][2] = stop1amplitudeHplussbottom1; ParticleStop1.Array_Decays[5][3] = 2; ParticleStop1.Array_Comments[5] = "# ~t_1 -> H+ ~b_1";
   ParticleStop1.Array_Decays[6][0] = PDGHplus; ParticleStop1.Array_Decays[6][1] = PDGsbottom2; ParticleStop1.Array_Decays[6][2] = stop1amplitudeHplussbottom2; ParticleStop1.Array_Decays[6][3] = 2; ParticleStop1.Array_Comments[6] = "# ~t_1 -> H+ ~b_2"; 

   ParticleStop1.Array_Decays[7][0] = PDGtop; ParticleStop1.Array_Decays[7][1] = PDGneutralino1; ParticleStop1.Array_Decays[7][2] = stop1amplitudetopneutralinoZ1; ParticleStop1.Array_Decays[7][3] = 2; ParticleStop1.Array_Comments[7] = "# ~t_1 -> t ~chi_10";
   ParticleStop1.Array_Decays[8][0] = PDGtop; ParticleStop1.Array_Decays[8][1] = PDGneutralino2; ParticleStop1.Array_Decays[8][2] = stop1amplitudetopneutralinoZ2; ParticleStop1.Array_Decays[8][3] = 2; ParticleStop1.Array_Comments[8] = "# ~t_1 -> t ~chi_20";
   ParticleStop1.Array_Decays[9][0] = PDGtop; ParticleStop1.Array_Decays[9][1] = PDGneutralino3; ParticleStop1.Array_Decays[9][2] = stop1amplitudetopneutralinoZ3; ParticleStop1.Array_Decays[9][3] = 2; ParticleStop1.Array_Comments[9] = "# ~t_1 -> t ~chi_30";
   ParticleStop1.Array_Decays[10][0] = PDGtop; ParticleStop1.Array_Decays[10][1] = PDGneutralino4; ParticleStop1.Array_Decays[10][2] = stop1amplitudetopneutralinoZ4; ParticleStop1.Array_Decays[10][3] = 2; ParticleStop1.Array_Comments[10] = "# ~t_1 -> t ~chi_40";
   ParticleStop1.Array_Decays[11][0] = PDGtop; ParticleStop1.Array_Decays[11][1] = PDGneutralino5; ParticleStop1.Array_Decays[11][2] = stop1amplitudetopneutralinoZ5; ParticleStop1.Array_Decays[11][3] = 2; ParticleStop1.Array_Comments[11] = "# ~t_1 -> t ~chi_50";

   ParticleStop1.Array_Decays[12][0] = PDGtop; ParticleStop1.Array_Decays[12][1] = PDGgravitino; ParticleStop1.Array_Decays[12][2] = stop1amplitudetopgravitino; ParticleStop1.Array_Decays[12][3] = 2; ParticleStop1.Array_Comments[12] = "# ~t_1 -> t ~G";

   for(int i = 0; i<ParticleStop1.No_of_Decays; i++) {
     if (ParticleStop1.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleStop1.Array_Comments[i] << " is negative = " << ParticleStop1.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleStop1.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   

   double Stop1_No_1to2_Decays = 0;

   Stop1_No_1to2_Decays = ParticleStop1.No_1to2_Decays + ParticleStop1.No_grav_Decays + ParticleStop1.No_NMSSM_Decays;
 
   for (int j = 0; j<Stop1_No_1to2_Decays; j++) {
     ParticleStop1.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<Stop1_No_1to2_Decays; j++) {
     ParticleStop1.two_width = ParticleStop1.two_width + ParticleStop1.Array_Decays[j][2];
   }
   for (int j=Stop1_No_1to2_Decays; j<ParticleStop1.No_of_Decays; j++) {
     ParticleStop1.three_width = ParticleStop1.three_width + ParticleStop1.Array_Decays[j][2];
   }
   ///Note currently no squark three body decays included - may change in future versions
   if ( ParticleStop1.three_width != ParticleStop1.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for stop1 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleStop1.No_of_Decays = Stop1_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleStop1.total_width = ParticleStop1.two_width;
     }
   else {
     ParticleStop1.total_width = ParticleStop1.two_width + ParticleStop1.three_width;
   }

   if ( ParticleStop1.total_width != ParticleStop1.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleStop1.No_of_Decays; i++) {
       //   fout << i << " " << ParticleStop1.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in stop1 total width \n");
     }
  }

 ///Stop2 Decays
 
 double stop2amplitudegluinotop=0, stop2amplitudetopneutralinoZ1=0, stop2amplitudetopneutralinoZ2=0, stop2amplitudetopneutralinoZ3=0, stop2amplitudetopneutralinoZ4=0, stop2amplitudebottomcharginoW1=0, stop2amplitudebottomcharginoW2=0, stop2amplitudeWbosonsbottom1=0, stop2amplitudeWbosonsbottom2=0, stop2amplitudeHplussbottom1=0, stop2amplitudeHplussbottom2=0, stop2amplitudeZbosonstop1=0, stop2amplitudehstop1=0, stop2amplitudeHstop1=0, stop2amplitudeAstop1=0, stop2amplitudetopgravitino=0;
 double stop2amplitudeH3stop1 = 0, stop2amplitudeA2stop1 = 0, stop2amplitudetopneutralinoZ5 = 0;

 if (flagstop2 == 1) {
   stop2amplitudegluinotop = squarkamplitudedecaygluinomix (mu(2,3), mt, mGluino, alphas, 2, thetat);
   stop2amplitudebottomcharginoW1 = squark2amplitudedecaycharginoW1mix (mu(2,3), mb, MCH1, g, thetaL2, thetaR2, thetat, beta, runmw, runmt, runmb,1);
   stop2amplitudebottomcharginoW2 = squark2amplitudedecaycharginoW2mix (mu(2,3), mb, MCH2, g, thetaL2, thetaR2, thetat, beta, runmw, runmt, runmb,1);
   
   stop2amplitudeWbosonsbottom1 = squark3amplitudedecaysquark3Wboson (mu(2,3), polemw, md(1,3), g, thetat, thetab, 1, 2, 2, 1);
   stop2amplitudeWbosonsbottom2 = squark3amplitudedecaysquark3Wboson (mu(2,3), polemw, md(2,3), g, thetat, thetab, 1, 2, 2, 2);
   stop2amplitudeHplussbottom1 = squark3amplitudedecaychargedHiggssquark3 (mu(2,3), mHpm, md(1,3), g, runmw, beta, thetat, thetab, greekmu, At, Ab, runmt, runmb, 2, 1);
   stop2amplitudeHplussbottom2 = squark3amplitudedecaychargedHiggssquark3 (mu(2,3), mHpm, md(2,3), g, runmw, beta, thetat, thetab, greekmu, At, Ab, runmt, runmb, 2, 2);

   stop2amplitudeZbosonstop1 = squark32amplitudedecaysquark3Zboson (mu(2,3), mu(1,3), polemz, g, gp, thetat); 
   
   if (nmssmIsIt == false) {
     stop2amplitudehstop1 = squark32amplitudedecayneutralHiggssquark3 (mu(2,3), mh0(1), mu(1,3), g, gp, runmw, beta, alpha, thetat, thetab, greekmu, At, Ab, runmt, runmb, 1, 'h');
     stop2amplitudeHstop1 = squark32amplitudedecayneutralHiggssquark3 (mu(2,3), mh0(2), mu(1,3), g, gp, runmw, beta, alpha, thetat, thetab, greekmu, At, Ab, runmt, runmb, 1, 'H'); 
     stop2amplitudeAstop1 = squark32amplitudedecayneutralHiggssquark3 (mu(2,3), mA0(1), mu(1,3), g, gp, runmw, beta, alpha, thetat, thetab, greekmu, At, Ab, runmt, runmb, 1, 'A');
     stop2amplitudetopneutralinoZ1 = squark3amplitudedecayneutralino (mu(2,3), mt, mneut(1), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 1);
     stop2amplitudetopneutralinoZ2 = squark3amplitudedecayneutralino (mu(2,3), mt, mneut(2), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 2);
     stop2amplitudetopneutralinoZ3 = squark3amplitudedecayneutralino (mu(2,3), mt, mneut(3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 3);
     stop2amplitudetopneutralinoZ4 = squark3amplitudedecayneutralino (mu(2,3), mt, mneut(4), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 4);
     stop2amplitudetopgravitino =squarkamplitudedecaygravitino(mu(2,3), mgravitino, mt, MPlreduced, gravonoff, upsquNLSP);
   }
   else if (nmssmIsIt == true) {
     stop2amplitudehstop1 = stop2amplitudedecaystop1CPevenhiggsNMSSM (mu(2,3), mu(1,3), mh0(1), runmt, thetat, CPEMix, beta, runmw, g, gp, At, mueff, lam, 1);
     stop2amplitudeHstop1 = stop2amplitudedecaystop1CPevenhiggsNMSSM (mu(2,3), mu(1,3), mh0(2), runmt, thetat, CPEMix, beta, runmw, g, gp, At, mueff, lam, 2);
     stop2amplitudeH3stop1 = stop2amplitudedecaystop1CPevenhiggsNMSSM (mu(2,3), mu(1,3), mh0(3), runmt, thetat, CPEMix, beta, runmw, g, gp, At, mueff, lam, 3);
     stop2amplitudeAstop1 = stop2amplitudedecaystop1CPoddhiggsNMSSM (mu(2,3), mu(1,3), mA0(1), runmt , thetat, CPOMix, beta, runmw, g, At, mueff, lam, 1);
     stop2amplitudeA2stop1 = stop2amplitudedecaystop1CPoddhiggsNMSSM (mu(2,3), mu(1,3), mA0(2), runmt , thetat, CPOMix, beta, runmw, g, At, mueff, lam, 2);
     stop2amplitudetopneutralinoZ1 = stopamplitudedecaytopneutralinoNMSSM (mu(2,3), mt, mneut(1), g, gp, thetat, mixNeut, runmt, runmw, beta, 2, 1);
     stop2amplitudetopneutralinoZ2 = stopamplitudedecaytopneutralinoNMSSM (mu(2,3), mt, mneut(2), g, gp, thetat, mixNeut, runmt, runmw, beta, 2, 2);
     stop2amplitudetopneutralinoZ3 = stopamplitudedecaytopneutralinoNMSSM (mu(2,3), mt, mneut(3), g, gp, thetat, mixNeut, runmt, runmw, beta, 2, 3);
     stop2amplitudetopneutralinoZ4 = stopamplitudedecaytopneutralinoNMSSM (mu(2,3), mt, mneut(4), g, gp, thetat, mixNeut, runmt, runmw, beta, 2, 4);
     stop2amplitudetopneutralinoZ5 = stopamplitudedecaytopneutralinoNMSSM (mu(2,3), mt, mneut(5), g, gp, thetat, mixNeut, runmt, runmw, beta, 2, 5);
   }
 
   ParticleStop2.Array_Decays[0][0] = PDGtop; ParticleStop2.Array_Decays[0][1] = PDGgluino; ParticleStop2.Array_Decays[0][2] = stop2amplitudegluinotop; ParticleStop2.Array_Decays[0][3] = 2; ParticleStop2.Array_Comments[0] = "# ~t_2 -> t ~g";
   ParticleStop2.Array_Decays[1][0] = PDGbottom; ParticleStop2.Array_Decays[1][1] = PDGchargino1; ParticleStop2.Array_Decays[1][2] = stop2amplitudebottomcharginoW1; ParticleStop2.Array_Decays[1][3] = 2; ParticleStop2.Array_Comments[1] = "# ~t_2 -> b ~chi_1+";
   ParticleStop2.Array_Decays[2][0] = PDGbottom; ParticleStop2.Array_Decays[2][1]= PDGchargino2; ParticleStop2.Array_Decays[2][2] = stop2amplitudebottomcharginoW2; ParticleStop2.Array_Decays[2][3] = 2; ParticleStop2.Array_Comments[2] = "# ~t_2 -> b ~chi_2+";
   
   ParticleStop2.Array_Decays[3][0] = PDGWplus; ParticleStop2.Array_Decays[3][1] = PDGsbottom1 ; ParticleStop2.Array_Decays[3][2] = stop2amplitudeWbosonsbottom1; ParticleStop2.Array_Decays[3][3] = 2; ParticleStop2.Array_Comments[3] = "# ~t_2 -> W+ ~b_1";
   ParticleStop2.Array_Decays[4][0] = PDGWplus; ParticleStop2.Array_Decays[4][1] = PDGsbottom2 ; ParticleStop2.Array_Decays[4][2] = stop2amplitudeWbosonsbottom2; ParticleStop2.Array_Decays[4][3] = 2; ParticleStop2.Array_Comments[4] = "# ~t_2 -> W+ ~b_2";
   ParticleStop2.Array_Decays[5][0] = PDGHplus; ParticleStop2.Array_Decays[5][1] = PDGsbottom1; ParticleStop2.Array_Decays[5][2] = stop2amplitudeHplussbottom1; ParticleStop2.Array_Decays[5][3] = 2; ParticleStop2.Array_Comments[5] = "# ~t_2 -> H+ ~b_1";
   ParticleStop2.Array_Decays[6][0] = PDGHplus; ParticleStop2.Array_Decays[6][1] = PDGsbottom2; ParticleStop2.Array_Decays[6][2] = stop2amplitudeHplussbottom2; ParticleStop2.Array_Decays[6][3] = 2; ParticleStop2.Array_Comments[6] = "# ~t_2 -> H+ ~b_2";
   ParticleStop2.Array_Decays[7][0] = PDGh0; ParticleStop2.Array_Decays[7][1] = PDGstop1; ParticleStop2.Array_Decays[7][2] = stop2amplitudehstop1; ParticleStop2.Array_Decays[7][3] = 2; ParticleStop2.Array_Comments[7] = "# ~t_2 -> h ~t_1";
   ParticleStop2.Array_Decays[8][0] = PDGH0; ParticleStop2.Array_Decays[8][1] = PDGstop1; ParticleStop2.Array_Decays[8][2] = stop2amplitudeHstop1; ParticleStop2.Array_Decays[8][3] = 2; ParticleStop2.Array_Comments[8] = "# ~t_2 -> H ~t_1";
   ParticleStop2.Array_Decays[9][0] = PDGA0; ParticleStop2.Array_Decays[9][1] = PDGstop1; ParticleStop2.Array_Decays[9][2] = stop2amplitudeAstop1; ParticleStop2.Array_Decays[9][3] = 2; ParticleStop2.Array_Comments[9] = "# ~t_2 -> A ~t_1";
   ParticleStop2.Array_Decays[10][0] = PDGstop1; ParticleStop2.Array_Decays[10][1] = PDGZboson; ParticleStop2.Array_Decays[10][2] = stop2amplitudeZbosonstop1; ParticleStop2.Array_Decays[10][3] = 2; ParticleStop2.Array_Comments[10] = "# ~t_2 -> Z ~t_1";

   ParticleStop2.Array_Decays[11][0] = PDGH3; ParticleStop2.Array_Decays[11][1] = PDGstop1; ParticleStop2.Array_Decays[11][2] = stop2amplitudeH3stop1; ParticleStop2.Array_Decays[11][3] = 2; ParticleStop2.Array_Comments[11] = "# ~t_2 -> H3 ~t_1";
   ParticleStop2.Array_Decays[12][0] = PDGA2; ParticleStop2.Array_Decays[12][1] = PDGstop1; ParticleStop2.Array_Decays[12][2] = stop2amplitudeA2stop1; ParticleStop2.Array_Decays[12][3] = 2; ParticleStop2.Array_Comments[12] = "# ~t_2 -> A2 ~t_1";
   ParticleStop2.Array_Decays[13][0] = PDGtop; ParticleStop2.Array_Decays[13][1] = PDGneutralino1; ParticleStop2.Array_Decays[13][2] = stop2amplitudetopneutralinoZ1; ParticleStop2.Array_Decays[13][3] = 2; ParticleStop2.Array_Comments[13] = "# ~t_2 -> t ~chi_10";
   ParticleStop2.Array_Decays[14][0] = PDGtop; ParticleStop2.Array_Decays[14][1] = PDGneutralino2; ParticleStop2.Array_Decays[14][2] = stop2amplitudetopneutralinoZ2; ParticleStop2.Array_Decays[14][3] = 2; ParticleStop2.Array_Comments[14] = "# ~t_2 -> t ~chi_20";
   ParticleStop2.Array_Decays[15][0] = PDGtop; ParticleStop2.Array_Decays[15][1] = PDGneutralino3; ParticleStop2.Array_Decays[15][2] = stop2amplitudetopneutralinoZ3; ParticleStop2.Array_Decays[15][3] = 2; ParticleStop2.Array_Comments[15] = "# ~t_2 -> t ~chi_30";
   ParticleStop2.Array_Decays[16][0] = PDGtop; ParticleStop2.Array_Decays[16][1] = PDGneutralino4; ParticleStop2.Array_Decays[16][2] = stop2amplitudetopneutralinoZ4; ParticleStop2.Array_Decays[16][3] = 2; ParticleStop2.Array_Comments[16] = "# ~t_2 -> t ~chi_40";
   ParticleStop2.Array_Decays[17][0] = PDGtop; ParticleStop2.Array_Decays[17][1] = PDGneutralino5; ParticleStop2.Array_Decays[17][2] = stop2amplitudetopneutralinoZ5; ParticleStop2.Array_Decays[17][3] = 2; ParticleStop2.Array_Comments[17] = "# ~t_2 -> t ~chi_50";
   
   ParticleStop2.Array_Decays[18][0] = PDGtop; ParticleStop2.Array_Decays[18][1] = PDGgravitino; ParticleStop2.Array_Decays[18][2] = stop2amplitudetopgravitino; ParticleStop2.Array_Decays[18][3] = 2; ParticleStop2.Array_Comments[18] = "# ~t_2 -> t ~G";

   for(int i = 0; i<ParticleStop2.No_of_Decays; i++) {
     if (ParticleStop2.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleStop2.Array_Comments[i] << " is negative = " << ParticleStop2.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleStop2.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }
    
   double Stop2_No_1to2_Decays = 0;

   Stop2_No_1to2_Decays = ParticleStop2.No_1to2_Decays + ParticleStop2.No_grav_Decays + ParticleStop2.No_NMSSM_Decays;
   
   for (int j = 0; j<Stop2_No_1to2_Decays; j++) {
     ParticleStop2.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<Stop2_No_1to2_Decays; j++) {
     ParticleStop2.two_width = ParticleStop2.two_width + ParticleStop2.Array_Decays[j][2];
   }
   for (int j=Stop2_No_1to2_Decays; j<ParticleStop2.No_of_Decays; j++) {
     ParticleStop2.three_width = ParticleStop2.three_width + ParticleStop2.Array_Decays[j][2];
   }
   ///Note currently no squark three body decays included - may change in future versions 
   if ( ParticleStop2.three_width != ParticleStop2.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for stop2 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleStop2.No_of_Decays = Stop2_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleStop2.total_width = ParticleStop2.two_width;
     }
   else {
     ParticleStop2.total_width = ParticleStop2.two_width + ParticleStop2.three_width;
   }
   
   if ( ParticleStop2.total_width != ParticleStop2.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleStop2.No_of_Decays; i++) {
       //   fout << i << " " << ParticleStop2.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in stop2 total width \n");
     }
   
 }

 /// Slepton Decays

 ///SelectonL decays
 
 double selectronLamplitudeeneutralinoZ1=0, selectronLamplitudeeneutralinoZ2=0, selectronLamplitudeeneutralinoZ3=0, selectronLamplitudeeneutralinoZ4=0, selectronLamplitudenuecharginoW1=0, selectronLamplitudenuecharginoW2=0, selectronLamplitudeelectrongravitino=0;
 double selectronLamplitudeeneutralinoZ5 = 0;

 if (flagselectronL == 1) {
   
   if (nmssmIsIt == false) {
     selectronLamplitudeeneutralinoZ1 = sleptonamplitudedecayleptonneutralino(me(1,1), mel, mneut(1), g, gp, mixNeut, 'L', 1);
     selectronLamplitudeeneutralinoZ2 = sleptonamplitudedecayleptonneutralino(me(1,1), mel, mneut(2), g, gp, mixNeut, 'L', 2);
     selectronLamplitudeeneutralinoZ3 = sleptonamplitudedecayleptonneutralino(me(1,1), mel, mneut(3), g, gp, mixNeut, 'L', 3);
     selectronLamplitudeeneutralinoZ4 = sleptonamplitudedecayleptonneutralino(me(1,1), mel, mneut(4), g, gp, mixNeut, 'L', 4);
     selectronLamplitudeelectrongravitino = squarkamplitudedecaygravitino(me(1,1), mgravitino, mel, MPlreduced, gravonoff, slepNLSP);
   }
   else if (nmssmIsIt == true) {
     selectronLamplitudeeneutralinoZ1 = sleptonamplitudedecayleptonneutralinoNMSSM (me(1,1), mel, mneut(1), g, gp, mixNeut, 'd', 'L', 1); 
     selectronLamplitudeeneutralinoZ2 = sleptonamplitudedecayleptonneutralinoNMSSM (me(1,1), mel, mneut(2), g, gp, mixNeut, 'd', 'L', 2);
     selectronLamplitudeeneutralinoZ3 = sleptonamplitudedecayleptonneutralinoNMSSM (me(1,1), mel, mneut(3), g, gp, mixNeut, 'd', 'L', 3);
     selectronLamplitudeeneutralinoZ4 = sleptonamplitudedecayleptonneutralinoNMSSM (me(1,1), mel, mneut(4), g, gp, mixNeut, 'd', 'L', 4);
     selectronLamplitudeeneutralinoZ5 = sleptonamplitudedecayleptonneutralinoNMSSM (me(1,1), mel, mneut(5), g, gp, mixNeut, 'd', 'L', 5);
   }
   selectronLamplitudenuecharginoW1 = sleptonamplitudedecaychargino (me(1,1), 0, MCH1, g, thetaL2, 1);
   selectronLamplitudenuecharginoW2 = sleptonamplitudedecaychargino (me(1,1), 0, MCH2, g, thetaL2, 2);
   
   ParticleSelectronL.Array_Decays[0][0] = PDGnuelectron; ParticleSelectronL.Array_Decays[0][1] = -PDGchargino1; ParticleSelectronL.Array_Decays[0][2] = selectronLamplitudenuecharginoW1; ParticleSelectronL.Array_Decays[0][3] = 2; ParticleSelectronL.Array_Comments[0] = "# ~e_L -> nu_e ~chi_1-";
   ParticleSelectronL.Array_Decays[1][0] = PDGnuelectron; ParticleSelectronL.Array_Decays[1][1] = -PDGchargino2; ParticleSelectronL.Array_Decays[1][2] = selectronLamplitudenuecharginoW2; ParticleSelectronL.Array_Decays[1][3] = 2; ParticleSelectronL.Array_Comments[1] = "# ~e_L -> nu_e ~chi_2-";
   ParticleSelectronL.Array_Decays[2][0] = PDGelectron; ParticleSelectronL.Array_Decays[2][1] = PDGneutralino1; ParticleSelectronL.Array_Decays[2][2] = selectronLamplitudeeneutralinoZ1; ParticleSelectronL.Array_Decays[2][3] = 2; ParticleSelectronL.Array_Comments[2] = "# ~e_L -> e- ~chi_10";
   ParticleSelectronL.Array_Decays[3][0] = PDGelectron; ParticleSelectronL.Array_Decays[3][1] = PDGneutralino2; ParticleSelectronL.Array_Decays[3][2] = selectronLamplitudeeneutralinoZ2; ParticleSelectronL.Array_Decays[3][3] = 2; ParticleSelectronL.Array_Comments[3] = "# ~e_L -> e- ~chi_20";
   ParticleSelectronL.Array_Decays[4][0] = PDGelectron; ParticleSelectronL.Array_Decays[4][1] = PDGneutralino3; ParticleSelectronL.Array_Decays[4][2] = selectronLamplitudeeneutralinoZ3; ParticleSelectronL.Array_Decays[4][3] = 2; ParticleSelectronL.Array_Comments[4] = "# ~e_L -> e- ~chi_30";
   ParticleSelectronL.Array_Decays[5][0] = PDGelectron; ParticleSelectronL.Array_Decays[5][1] = PDGneutralino4; ParticleSelectronL.Array_Decays[5][2] = selectronLamplitudeeneutralinoZ4; ParticleSelectronL.Array_Decays[5][3] = 2; ParticleSelectronL.Array_Comments[5] = "# ~e_L -> e- ~chi_40";
   ParticleSelectronL.Array_Decays[6][0] = PDGelectron; ParticleSelectronL.Array_Decays[6][1] = PDGneutralino5; ParticleSelectronL.Array_Decays[6][2] = selectronLamplitudeeneutralinoZ5; ParticleSelectronL.Array_Decays[6][3] = 2; ParticleSelectronL.Array_Comments[6] = "# ~e_L -> e- ~chi_50";

   ParticleSelectronL.Array_Decays[7][0] = PDGelectron; ParticleSelectronL.Array_Decays[7][1] = PDGgravitino; ParticleSelectronL.Array_Decays[7][2] = selectronLamplitudeelectrongravitino; ParticleSelectronL.Array_Decays[7][3] = 2; ParticleSelectronL.Array_Comments[7] = "# ~e_L -> e- ~G";

   for(int i = 0; i<ParticleSelectronL.No_of_Decays; i++) {
     if (ParticleSelectronL.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSelectronL.Array_Comments[i] << " is negative = " << ParticleSelectronL.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSelectronL.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   

   double SelectronL_No_1to2_Decays = 0;
   
   SelectronL_No_1to2_Decays = ParticleSelectronL.No_1to2_Decays + ParticleSelectronL.No_grav_Decays + ParticleSelectronL.No_NMSSM_Decays;
 
   for (int j = 0; j<SelectronL_No_1to2_Decays; j++) {
     ParticleSelectronL.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<SelectronL_No_1to2_Decays; j++) {
     ParticleSelectronL.two_width = ParticleSelectronL.two_width + ParticleSelectronL.Array_Decays[j][2];
   }
   for (int j=SelectronL_No_1to2_Decays; j<ParticleSelectronL.No_of_Decays; j++) {
     ParticleSelectronL.three_width = ParticleSelectronL.three_width + ParticleSelectronL.Array_Decays[j][2];
   }
   ///Note currently no slepton three body decays included - may change in future versions
   if ( ParticleSelectronL.three_width != ParticleSelectronL.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for selectronL - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSelectronL.No_of_Decays = SelectronL_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSelectronL.total_width = ParticleSelectronL.two_width;
     }
   else {
     ParticleSelectronL.total_width = ParticleSelectronL.two_width + ParticleSelectronL.three_width;
   }

   if ( ParticleSelectronL.total_width != ParticleSelectronL.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSelectronL.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSelectronL.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in selectronL total width \n");
     }
   
 }


 ///SelectonR decays
 
 double selectronRamplitudeeneutralinoZ1=0, selectronRamplitudeeneutralinoZ2=0, selectronRamplitudeeneutralinoZ3=0, selectronRamplitudeeneutralinoZ4=0, selectronRamplitudeelectrongravitino=0;

 double selectronRamplitudeeneutralinoZ5 = 0;

 if (flagselectronR == 1) {

   if (nmssmIsIt == false) {
     selectronRamplitudeeneutralinoZ1 = sleptonamplitudedecayleptonneutralino(me(2,1), mel, mneut(1), g, gp, mixNeut, 'R', 1);
     selectronRamplitudeeneutralinoZ2 = sleptonamplitudedecayleptonneutralino(me(2,1), mel, mneut(2), g, gp, mixNeut, 'R', 2);
     selectronRamplitudeeneutralinoZ3 = sleptonamplitudedecayleptonneutralino(me(2,1), mel, mneut(3), g, gp, mixNeut, 'R', 3);
     selectronRamplitudeeneutralinoZ4 = sleptonamplitudedecayleptonneutralino(me(2,1), mel, mneut(4), g, gp, mixNeut, 'R', 4);
     selectronRamplitudeelectrongravitino = squarkamplitudedecaygravitino(me(2,1), mgravitino, mel, MPlreduced, gravonoff, slepNLSP);
   }
   else if (nmssmIsIt == true) {
     selectronRamplitudeeneutralinoZ1 = sleptonamplitudedecayleptonneutralinoNMSSM (me(2,1), mel, mneut(1), g, gp, mixNeut, 'd', 'R', 1);     selectronRamplitudeeneutralinoZ2 = sleptonamplitudedecayleptonneutralinoNMSSM (me(2,1), mel, mneut(2), g, gp, mixNeut, 'd', 'R', 2);
     selectronRamplitudeeneutralinoZ3 = sleptonamplitudedecayleptonneutralinoNMSSM (me(2,1), mel, mneut(3), g, gp, mixNeut, 'd', 'R', 3);
     selectronRamplitudeeneutralinoZ4 = sleptonamplitudedecayleptonneutralinoNMSSM (me(2,1), mel, mneut(4), g, gp, mixNeut, 'd', 'R', 4);
     selectronRamplitudeeneutralinoZ5 = sleptonamplitudedecayleptonneutralinoNMSSM (me(2,1), mel, mneut(5), g, gp, mixNeut, 'd', 'R', 5);
   }

   ParticleSelectronR.Array_Decays[0][0] = PDGelectron; ParticleSelectronR.Array_Decays[0][1] = PDGneutralino1; ParticleSelectronR.Array_Decays[0][2] = selectronRamplitudeeneutralinoZ1; ParticleSelectronR.Array_Decays[0][3] = 2; ParticleSelectronR.Array_Comments[0] = "# ~e_R -> e- ~chi_10";
   ParticleSelectronR.Array_Decays[1][0] = PDGelectron; ParticleSelectronR.Array_Decays[1][1] = PDGneutralino2; ParticleSelectronR.Array_Decays[1][2] = selectronRamplitudeeneutralinoZ2; ParticleSelectronR.Array_Decays[1][3] = 2; ParticleSelectronR.Array_Comments[1] = "# ~e_R -> e- ~chi_20";
   ParticleSelectronR.Array_Decays[2][0] = PDGelectron; ParticleSelectronR.Array_Decays[2][1] = PDGneutralino3; ParticleSelectronR.Array_Decays[2][2] = selectronRamplitudeeneutralinoZ3; ParticleSelectronR.Array_Decays[2][3] = 2; ParticleSelectronR.Array_Comments[2] = "# ~e_R -> e- ~chi_30";
   ParticleSelectronR.Array_Decays[3][0] = PDGelectron; ParticleSelectronR.Array_Decays[3][1] = PDGneutralino4; ParticleSelectronR.Array_Decays[3][2] = selectronRamplitudeeneutralinoZ4; ParticleSelectronR.Array_Decays[3][3] = 2; ParticleSelectronR.Array_Comments[3] = "# ~e_R -> e- ~chi_40";
   ParticleSelectronR.Array_Decays[4][0] = PDGelectron; ParticleSelectronR.Array_Decays[4][1] = PDGneutralino5; ParticleSelectronR.Array_Decays[4][2] = selectronRamplitudeeneutralinoZ5; ParticleSelectronR.Array_Decays[4][3] = 2; ParticleSelectronR.Array_Comments[4] = "# ~e_R -> e- ~chi_50";
 
   ParticleSelectronR.Array_Decays[5][0] = PDGelectron; ParticleSelectronR.Array_Decays[5][1] = PDGgravitino; ParticleSelectronR.Array_Decays[5][2] = selectronRamplitudeelectrongravitino; ParticleSelectronR.Array_Decays[5][3] = 2; ParticleSelectronR.Array_Comments[5] = "# ~e_R -> e- ~G";

   for(int i = 0; i<ParticleSelectronR.No_of_Decays; i++) {
     if (ParticleSelectronR.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSelectronR.Array_Comments[i] << " is negative = " << ParticleSelectronR.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSelectronR.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   

   double SelectronR_No_1to2_Decays = 0;

   SelectronR_No_1to2_Decays = ParticleSelectronR.No_1to2_Decays + ParticleSelectronR.No_grav_Decays + ParticleSelectronR.No_NMSSM_Decays;
 
   for (int j = 0; j<SelectronR_No_1to2_Decays; j++) {
     ParticleSelectronR.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<SelectronR_No_1to2_Decays; j++) {
     ParticleSelectronR.two_width = ParticleSelectronR.two_width + ParticleSelectronR.Array_Decays[j][2];
   }
   for (int j=SelectronR_No_1to2_Decays; j<ParticleSelectronR.No_of_Decays; j++) {
     ParticleSelectronR.three_width = ParticleSelectronR.three_width + ParticleSelectronR.Array_Decays[j][2];
   }
   ///Note currently no slepton three body decays included - may change in future versions
   if ( ParticleSelectronR.three_width != ParticleSelectronR.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for selectronR - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSelectronR.No_of_Decays = SelectronR_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSelectronR.total_width = ParticleSelectronR.two_width;
     }
   else {
     ParticleSelectronR.total_width = ParticleSelectronR.two_width + ParticleSelectronR.three_width;
   }

   if ( ParticleSelectronR.total_width != ParticleSelectronR.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSelectronR.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSelectronR.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in selectronR total width \n");
     }
  
 }


 ///SmuonL decays

 double smuonLamplitudemuneutralinoZ1=0, smuonLamplitudemuneutralinoZ2=0, smuonLamplitudemuneutralinoZ3=0, smuonLamplitudemuneutralinoZ4=0, smuonLamplitudenumucharginoW1=0, smuonLamplitudenumucharginoW2=0, smuonLamplitudemuongravitino=0;
 double smuonLamplitudemuneutralinoZ5 = 0;

 if(flagsmuonL == 1) {

   if (nmssmIsIt == false) {
     smuonLamplitudemuneutralinoZ1 = sleptonamplitudedecayleptonneutralino(me(1,2), mmu, mneut(1), g, gp, mixNeut, 'L', 1);
     smuonLamplitudemuneutralinoZ2 = sleptonamplitudedecayleptonneutralino(me(1,2), mmu, mneut(2), g, gp, mixNeut, 'L', 2);
     smuonLamplitudemuneutralinoZ3 = sleptonamplitudedecayleptonneutralino(me(1,2), mmu, mneut(3), g, gp, mixNeut, 'L', 3);
     smuonLamplitudemuneutralinoZ4 = sleptonamplitudedecayleptonneutralino(me(1,2), mmu, mneut(4), g, gp, mixNeut, 'L', 4);
     smuonLamplitudemuongravitino = squarkamplitudedecaygravitino(me(1,2), mgravitino, mmu, MPlreduced, gravonoff, slepNLSP);
   }
   else if (nmssmIsIt == true) {
     smuonLamplitudemuneutralinoZ1 = sleptonamplitudedecayleptonneutralinoNMSSM (me(1,2), mmu, mneut(1), g, gp, mixNeut, 'd', 'L', 1); 
     smuonLamplitudemuneutralinoZ2 = sleptonamplitudedecayleptonneutralinoNMSSM (me(1,2), mmu, mneut(2), g, gp, mixNeut, 'd', 'L', 2);
     smuonLamplitudemuneutralinoZ3 = sleptonamplitudedecayleptonneutralinoNMSSM (me(1,2), mmu, mneut(3), g, gp, mixNeut, 'd', 'L', 3);
     smuonLamplitudemuneutralinoZ4 = sleptonamplitudedecayleptonneutralinoNMSSM (me(1,2), mmu, mneut(4), g, gp, mixNeut, 'd', 'L', 4);
     smuonLamplitudemuneutralinoZ5 = sleptonamplitudedecayleptonneutralinoNMSSM (me(1,2), mmu, mneut(5), g, gp, mixNeut, 'd', 'L', 5);
   }
   smuonLamplitudenumucharginoW1 = sleptonamplitudedecaychargino (me(1,2), 0, MCH1, g, thetaL2, 1);
   smuonLamplitudenumucharginoW2 = sleptonamplitudedecaychargino (me(1,2), 0, MCH2, g, thetaL2, 2);

   ParticleSmuonL.Array_Decays[0][0] = PDGnumuon; ParticleSmuonL.Array_Decays[0][1] = -PDGchargino1; ParticleSmuonL.Array_Decays[0][2] = smuonLamplitudenumucharginoW1; ParticleSmuonL.Array_Decays[0][3] = 2; ParticleSmuonL.Array_Comments[0] = "# ~mu_L -> nu_mu ~chi_1-";
   ParticleSmuonL.Array_Decays[1][0] = PDGnumuon; ParticleSmuonL.Array_Decays[1][1] = -PDGchargino2; ParticleSmuonL.Array_Decays[1][2] = smuonLamplitudenumucharginoW2; ParticleSmuonL.Array_Decays[1][3] = 2; ParticleSmuonL.Array_Comments[1] = "# ~mu_L -> nu_mu ~chi_2-";
   ParticleSmuonL.Array_Decays[2][0] = PDGmuon; ParticleSmuonL.Array_Decays[2][1] = PDGneutralino1; ParticleSmuonL.Array_Decays[2][2] = smuonLamplitudemuneutralinoZ1; ParticleSmuonL.Array_Decays[2][3] = 2; ParticleSmuonL.Array_Comments[2] = "# ~mu_L -> mu- ~chi_10";
   ParticleSmuonL.Array_Decays[3][0] = PDGmuon; ParticleSmuonL.Array_Decays[3][1] = PDGneutralino2; ParticleSmuonL.Array_Decays[3][2] = smuonLamplitudemuneutralinoZ2; ParticleSmuonL.Array_Decays[3][3] = 2; ParticleSmuonL.Array_Comments[3] = "# ~mu_L -> mu- ~chi_20";
   ParticleSmuonL.Array_Decays[4][0] = PDGmuon; ParticleSmuonL.Array_Decays[4][1] = PDGneutralino3; ParticleSmuonL.Array_Decays[4][2] = smuonLamplitudemuneutralinoZ3; ParticleSmuonL.Array_Decays[4][3] = 2; ParticleSmuonL.Array_Comments[4] = "# ~mu_L -> mu- ~chi_30";
   ParticleSmuonL.Array_Decays[5][0] = PDGmuon; ParticleSmuonL.Array_Decays[5][1] = PDGneutralino4; ParticleSmuonL.Array_Decays[5][2] = smuonLamplitudemuneutralinoZ4; ParticleSmuonL.Array_Decays[5][3] = 2; ParticleSmuonL.Array_Comments[5] = "# ~mu_L -> mu- ~chi_40";
   ParticleSmuonL.Array_Decays[6][0] = PDGmuon; ParticleSmuonL.Array_Decays[6][1] = PDGneutralino4; ParticleSmuonL.Array_Decays[6][2] = smuonLamplitudemuneutralinoZ5; ParticleSmuonL.Array_Decays[6][3] = 2; ParticleSmuonL.Array_Comments[6] = "# ~mu_L -> mu- ~chi_50";

   ParticleSmuonL.Array_Decays[7][0] = PDGmuon; ParticleSmuonL.Array_Decays[7][1] = PDGgravitino; ParticleSmuonL.Array_Decays[7][2] = smuonLamplitudemuongravitino; ParticleSmuonL.Array_Decays[7][3] = 2; ParticleSmuonL.Array_Comments[7] = "# ~mu_L -> mu- ~G";

   for(int i = 0; i<ParticleSmuonL.No_of_Decays; i++) {
     if (ParticleSmuonL.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSmuonL.Array_Comments[i] << " is negative = " << ParticleSmuonL.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSmuonL.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   
   
   double SmuonL_No_1to2_Decays = 0;

   SmuonL_No_1to2_Decays = ParticleSmuonL.No_1to2_Decays + ParticleSmuonL.No_grav_Decays + ParticleSmuonL.No_NMSSM_Decays;
 
   for (int j = 0; j<SmuonL_No_1to2_Decays; j++) {
     ParticleSmuonL.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<SmuonL_No_1to2_Decays; j++) {
     ParticleSmuonL.two_width = ParticleSmuonL.two_width + ParticleSmuonL.Array_Decays[j][2];
   }
   for (int j=SmuonL_No_1to2_Decays; j<ParticleSmuonL.No_of_Decays; j++) {
     ParticleSmuonL.three_width = ParticleSmuonL.three_width + ParticleSmuonL.Array_Decays[j][2];
   }
   ///Note currently no slepton three body decays included - may change in future versions
   if ( ParticleSmuonL.three_width != ParticleSmuonL.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for smuonL - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSmuonL.No_of_Decays = SmuonL_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSmuonL.total_width = ParticleSmuonL.two_width;
     }
   else {
     ParticleSmuonL.total_width = ParticleSmuonL.two_width + ParticleSmuonL.three_width;
   }

   if ( ParticleSmuonL.total_width != ParticleSmuonL.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSmuonL.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSmuonL.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in smuonL total width \n");
     }
  
 }

 ///SmuonR decays
 
 double smuonRamplitudemuneutralinoZ1=0, smuonRamplitudemuneutralinoZ2=0, smuonRamplitudemuneutralinoZ3=0, smuonRamplitudemuneutralinoZ4=0, smuonRamplitudemuongravitino=0;
 double smuonRamplitudemuneutralinoZ5 = 0;
 
 if (flagsmuonR == 1) {
   if (nmssmIsIt == false) {
     smuonRamplitudemuneutralinoZ1 = sleptonamplitudedecayleptonneutralino(me(2,2), mmu, mneut(1), g, gp, mixNeut, 'R', 1);
     smuonRamplitudemuneutralinoZ2 = sleptonamplitudedecayleptonneutralino(me(2,2), mmu, mneut(2), g, gp, mixNeut, 'R', 2);
     smuonRamplitudemuneutralinoZ3 = sleptonamplitudedecayleptonneutralino(me(2,2), mmu, mneut(3), g, gp, mixNeut, 'R', 3);
     smuonRamplitudemuneutralinoZ4 = sleptonamplitudedecayleptonneutralino(me(2,2), mmu, mneut(4), g, gp, mixNeut, 'R', 4);
     smuonRamplitudemuongravitino = squarkamplitudedecaygravitino (me(2,2), mgravitino, mmu, MPlreduced, gravonoff, slepNLSP);
   }
   else if (nmssmIsIt == true) {
     smuonRamplitudemuneutralinoZ1 = sleptonamplitudedecayleptonneutralinoNMSSM (me(2,2), mmu, mneut(1), g, gp, mixNeut, 'd', 'R', 1); 
     smuonRamplitudemuneutralinoZ2 = sleptonamplitudedecayleptonneutralinoNMSSM (me(2,2), mmu, mneut(2), g, gp, mixNeut, 'd', 'R', 2);
     smuonRamplitudemuneutralinoZ3 = sleptonamplitudedecayleptonneutralinoNMSSM (me(2,2), mmu, mneut(3), g, gp, mixNeut, 'd', 'R', 3);
     smuonRamplitudemuneutralinoZ4 = sleptonamplitudedecayleptonneutralinoNMSSM (me(2,2), mmu, mneut(4), g, gp, mixNeut, 'd', 'R', 4);
     smuonRamplitudemuneutralinoZ5 = sleptonamplitudedecayleptonneutralinoNMSSM (me(2,2), mmu, mneut(5), g, gp, mixNeut, 'd', 'R', 5);
   }

   ParticleSmuonR.Array_Decays[0][0] = PDGmuon; ParticleSmuonR.Array_Decays[0][1] = PDGneutralino1; ParticleSmuonR.Array_Decays[0][2] = smuonRamplitudemuneutralinoZ1; ParticleSmuonR.Array_Decays[0][3] = 2; ParticleSmuonR.Array_Comments[0] = "# ~mu_R -> mu- ~chi_10";
   ParticleSmuonR.Array_Decays[1][0] = PDGmuon; ParticleSmuonR.Array_Decays[1][1] = PDGneutralino2; ParticleSmuonR.Array_Decays[1][2] = smuonRamplitudemuneutralinoZ2; ParticleSmuonR.Array_Decays[1][3] = 2; ParticleSmuonR.Array_Comments[1] = "# ~mu_R -> mu- ~chi_20";
   ParticleSmuonR.Array_Decays[2][0] = PDGmuon; ParticleSmuonR.Array_Decays[2][1] = PDGneutralino3; ParticleSmuonR.Array_Decays[2][2] = smuonRamplitudemuneutralinoZ3; ParticleSmuonR.Array_Decays[2][3] = 2; ParticleSmuonR.Array_Comments[2] = "# ~mu_R -> mu- ~chi_30";
   ParticleSmuonR.Array_Decays[3][0] = PDGmuon; ParticleSmuonR.Array_Decays[3][1] = PDGneutralino4; ParticleSmuonR.Array_Decays[3][2] = smuonRamplitudemuneutralinoZ4; ParticleSmuonR.Array_Decays[3][3] = 2; ParticleSmuonR.Array_Comments[3] = "# ~mu_R -> mu- ~chi_40";
   ParticleSmuonR.Array_Decays[4][0] = PDGmuon; ParticleSmuonR.Array_Decays[4][1] = PDGneutralino5; ParticleSmuonR.Array_Decays[4][2] = smuonRamplitudemuneutralinoZ5; ParticleSmuonR.Array_Decays[4][3] = 2; ParticleSmuonR.Array_Comments[4] = "# ~mu_R -> mu- ~chi_50";
 
   ParticleSmuonR.Array_Decays[5][0] = PDGmuon; ParticleSmuonR.Array_Decays[5][1] = PDGgravitino; ParticleSmuonR.Array_Decays[5][2] = smuonRamplitudemuongravitino; ParticleSmuonR.Array_Decays[5][3] = 2; ParticleSmuonR.Array_Comments[5] = "# ~mu_R -> mu- ~G";

   for(int i = 0; i<ParticleSmuonR.No_of_Decays; i++) {
     if (ParticleSmuonR.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSmuonR.Array_Comments[i] << " is negative = " << ParticleSmuonR.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSmuonR.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }     
 
   double SmuonR_No_1to2_Decays = 0;

   SmuonR_No_1to2_Decays = ParticleSmuonR.No_1to2_Decays + ParticleSmuonR.No_grav_Decays + ParticleSmuonR.No_NMSSM_Decays;
 
   for (int j = 0; j<SmuonR_No_1to2_Decays; j++) {
     ParticleSmuonR.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<SmuonR_No_1to2_Decays; j++) {
     ParticleSmuonR.two_width = ParticleSmuonR.two_width + ParticleSmuonR.Array_Decays[j][2];
   }
   for (int j=SmuonR_No_1to2_Decays; j<ParticleSmuonR.No_of_Decays; j++) {
     ParticleSmuonR.three_width = ParticleSmuonR.three_width + ParticleSmuonL.Array_Decays[j][2];
   }
   ///Note currently no slepton three body decays included - may change in future versions
   if ( ParticleSmuonR.three_width != ParticleSmuonR.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for smuonR - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSmuonR.No_of_Decays = SmuonR_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSmuonR.total_width = ParticleSmuonR.two_width;
     }
   else {
     ParticleSmuonR.total_width = ParticleSmuonR.two_width + ParticleSmuonR.three_width;
   }

   if ( ParticleSmuonR.total_width != ParticleSmuonR.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSmuonR.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSmuonR.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in smuonR total width \n");
     }
 }

 ///Selectron sneutrino decays

 double snueamplitudenueneutralinoZ1=0, snueamplitudenueneutralinoZ2=0, snueamplitudenueneutralinoZ3=0, snueamplitudenueneutralinoZ4=0, snueamplitudeecharginoW1=0, snueamplitudeecharginoW2=0, snueamplitudenuegravitino=0;
 double snueamplitudenueneutralinoZ5 = 0;

 if (flagsnueL == 1) {
   if (nmssmIsIt == false) {
     snueamplitudenueneutralinoZ1 = sneutrinoamplitudedecayneutrinoneutralino (msnu(1), 0, mneut(1), g, gp, mixNeut, 1);
     snueamplitudenueneutralinoZ2 = sneutrinoamplitudedecayneutrinoneutralino (msnu(1), 0, mneut(2), g, gp, mixNeut, 2);
     snueamplitudenueneutralinoZ3 = sneutrinoamplitudedecayneutrinoneutralino (msnu(1), 0, mneut(3), g, gp, mixNeut, 3);
     snueamplitudenueneutralinoZ4 = sneutrinoamplitudedecayneutrinoneutralino (msnu(1), 0, mneut(4), g, gp, mixNeut, 4);
     snueamplitudenuegravitino = squarkamplitudedecaygravitino(msnu(1),mgravitino, 0, MPlreduced, gravonoff, snuNLSP);
   }
   else if (nmssmIsIt == true) {
     snueamplitudenueneutralinoZ1 = sleptonamplitudedecayleptonneutralinoNMSSM (msnu(1), 0, mneut(1), g, gp, mixNeut, 'u', 'L', 1); 
     snueamplitudenueneutralinoZ2 = sleptonamplitudedecayleptonneutralinoNMSSM (msnu(1), 0, mneut(2), g, gp, mixNeut, 'u', 'L', 2);
     snueamplitudenueneutralinoZ3 = sleptonamplitudedecayleptonneutralinoNMSSM (msnu(1), 0, mneut(3), g, gp, mixNeut, 'u', 'L', 3);
     snueamplitudenueneutralinoZ4 = sleptonamplitudedecayleptonneutralinoNMSSM (msnu(1), 0, mneut(4), g, gp, mixNeut, 'u', 'L', 4);
     snueamplitudenueneutralinoZ5 = sleptonamplitudedecayleptonneutralinoNMSSM (msnu(1), 0, mneut(5), g, gp, mixNeut, 'u', 'L', 5);
   }

   snueamplitudeecharginoW1 = sleptonamplitudedecaychargino (msnu(1), mel, MCH1, g, thetaR2, 1);
   snueamplitudeecharginoW2 = sleptonamplitudedecaychargino (msnu(1), mel, MCH2, g, thetaR2, 2);

   ParticleSnue.Array_Decays[0][0] = PDGelectron; ParticleSnue.Array_Decays[0][1] = PDGchargino1; ParticleSnue.Array_Decays[0][2] = snueamplitudeecharginoW1; ParticleSnue.Array_Decays[0][3] = 2; ParticleSnue.Array_Comments[0] = "# ~nu_eL -> e- ~chi_1+";
   ParticleSnue.Array_Decays[1][0] = PDGelectron; ParticleSnue.Array_Decays[1][1] = PDGchargino2; ParticleSnue.Array_Decays[1][2] = snueamplitudeecharginoW2; ParticleSnue.Array_Decays[1][3] = 2; ParticleSnue.Array_Comments[1] = "# ~nu_eL -> e- ~chi_2+";
   ParticleSnue.Array_Decays[2][0] = PDGnuelectron; ParticleSnue.Array_Decays[2][1] = PDGneutralino1; ParticleSnue.Array_Decays[2][2] = snueamplitudenueneutralinoZ1; ParticleSnue.Array_Decays[2][3] = 2; ParticleSnue.Array_Comments[2] = "# ~nu_eL -> nu_e ~chi_10" ;
   ParticleSnue.Array_Decays[3][0] = PDGnuelectron; ParticleSnue.Array_Decays[3][1] = PDGneutralino2; ParticleSnue.Array_Decays[3][2] = snueamplitudenueneutralinoZ2; ParticleSnue.Array_Decays[3][3] = 2; ParticleSnue.Array_Comments[3] = "# ~nu_eL -> nu_e ~chi_20";
   ParticleSnue.Array_Decays[4][0] = PDGnuelectron; ParticleSnue.Array_Decays[4][1] = PDGneutralino3; ParticleSnue.Array_Decays[4][2] = snueamplitudenueneutralinoZ3; ParticleSnue.Array_Decays[4][3] = 2; ParticleSnue.Array_Comments[4] = "# ~nu_eL -> nu_e ~chi_30";
   ParticleSnue.Array_Decays[5][0] = PDGnuelectron; ParticleSnue.Array_Decays[5][1] = PDGneutralino4; ParticleSnue.Array_Decays[5][2] = snueamplitudenueneutralinoZ4; ParticleSnue.Array_Decays[5][3] = 2; ParticleSnue.Array_Comments[5] = "# ~nu_eL -> nu_e ~chi_40";
   ParticleSnue.Array_Decays[6][0] = PDGnuelectron; ParticleSnue.Array_Decays[6][1] = PDGneutralino5; ParticleSnue.Array_Decays[6][2] = snueamplitudenueneutralinoZ5; ParticleSnue.Array_Decays[6][3] = 2; ParticleSnue.Array_Comments[6] = "# ~nu_eL -> nu_e ~chi_50";

   ParticleSnue.Array_Decays[7][0] = PDGnuelectron; ParticleSnue.Array_Decays[7][1] = PDGgravitino; ParticleSnue.Array_Decays[7][2] = snueamplitudenuegravitino; ParticleSnue.Array_Decays[7][3] = 2; ParticleSnue.Array_Comments[7] = "# ~nu_eL -> nu_e ~G";

   for(int i = 0; i<ParticleSnue.No_of_Decays; i++) {
     if (ParticleSnue.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSnue.Array_Comments[i] << " is negative = " << ParticleSnue.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSnue.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }     

   double Snue_No_1to2_Decays = 0;

   Snue_No_1to2_Decays = ParticleSnue.No_1to2_Decays + ParticleSnue.No_grav_Decays + ParticleSnue.No_NMSSM_Decays;
  
   for (int j = 0; j<Snue_No_1to2_Decays; j++) {
     ParticleSnue.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<Snue_No_1to2_Decays; j++) {
     ParticleSnue.two_width = ParticleSnue.two_width + ParticleSnue.Array_Decays[j][2];
   }
   for (int j=Snue_No_1to2_Decays; j<ParticleSnue.No_of_Decays; j++) {
     ParticleSnue.three_width = ParticleSnue.three_width + ParticleSnue.Array_Decays[j][2];
   }
   ///Note currently no slepton three body decays included - may change in future versions
   if ( ParticleSnue.three_width != ParticleSnue.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for snue - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSnue.No_of_Decays = Snue_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSnue.total_width = ParticleSnue.two_width;
     }
   else {
     ParticleSnue.total_width = ParticleSnue.two_width + ParticleSnue.three_width;
   }

   if ( ParticleSnue.total_width != ParticleSnue.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSnue.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSnue.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in snue total width \n");
     }  
 }
 

 ///Smuon sneutrino decays

 double snumuamplitudenumuneutralinoZ1=0, snumuamplitudenumuneutralinoZ2=0, snumuamplitudenumuneutralinoZ3=0, snumuamplitudenumuneutralinoZ4=0, snumuamplitudemucharginoW1=0, snumuamplitudemucharginoW2=0, snumuamplitudenumugravitino=0;
 double snumuamplitudenumuneutralinoZ5 = 0;

 if (flagsnumuL == 1) {
   if (nmssmIsIt == false) { 
     snumuamplitudenumuneutralinoZ1 = sneutrinoamplitudedecayneutrinoneutralino (msnu(2), 0, mneut(1), g, gp, mixNeut, 1);
     snumuamplitudenumuneutralinoZ2 = sneutrinoamplitudedecayneutrinoneutralino (msnu(2), 0, mneut(2), g, gp, mixNeut, 2);
     snumuamplitudenumuneutralinoZ3 = sneutrinoamplitudedecayneutrinoneutralino (msnu(2), 0, mneut(3), g, gp, mixNeut, 3);
     snumuamplitudenumuneutralinoZ4 = sneutrinoamplitudedecayneutrinoneutralino (msnu(2), 0, mneut(4), g, gp, mixNeut, 4);
     snumuamplitudenumugravitino = squarkamplitudedecaygravitino (msnu(2), mgravitino, 0, MPlreduced, gravonoff, snuNLSP);
   }
   else if (nmssmIsIt == true) {
     snumuamplitudenumuneutralinoZ1 = sleptonamplitudedecayleptonneutralinoNMSSM (msnu(2), 0, mneut(1), g, gp, mixNeut, 'u', 'L', 1); 
     snumuamplitudenumuneutralinoZ2 = sleptonamplitudedecayleptonneutralinoNMSSM (msnu(2), 0, mneut(2), g, gp, mixNeut, 'u', 'L', 2);
     snumuamplitudenumuneutralinoZ3 = sleptonamplitudedecayleptonneutralinoNMSSM (msnu(2), 0, mneut(3), g, gp, mixNeut, 'u', 'L', 3);
     snumuamplitudenumuneutralinoZ4 = sleptonamplitudedecayleptonneutralinoNMSSM (msnu(2), 0, mneut(4), g, gp, mixNeut, 'u', 'L', 4);
     snumuamplitudenumuneutralinoZ5 = sleptonamplitudedecayleptonneutralinoNMSSM (msnu(2), 0, mneut(5), g, gp, mixNeut, 'u', 'L', 5);
   }
 
   snumuamplitudemucharginoW1 = sleptonamplitudedecaychargino (msnu(2), mmu, MCH1, g, thetaR2, 1);
   snumuamplitudemucharginoW2 = sleptonamplitudedecaychargino (msnu(2), mmu, MCH2, g, thetaR2, 2);

   ParticleSnumu.Array_Decays[0][0] = PDGmuon; ParticleSnumu.Array_Decays[0][1] = PDGchargino1; ParticleSnumu.Array_Decays[0][2] = snumuamplitudemucharginoW1; ParticleSnumu.Array_Decays[0][3] = 2; ParticleSnumu.Array_Comments[0] = "# ~nu_muL -> mu- ~chi_1+";
   ParticleSnumu.Array_Decays[1][0] = PDGmuon; ParticleSnumu.Array_Decays[1][1] = PDGchargino2; ParticleSnumu.Array_Decays[1][2] = snumuamplitudemucharginoW2; ParticleSnumu.Array_Decays[1][3] = 2; ParticleSnumu.Array_Comments[1] = "# ~nu_muL -> mu- ~chi_2+";
   ParticleSnumu.Array_Decays[2][0] = PDGnumuon; ParticleSnumu.Array_Decays[2][1] = PDGneutralino1; ParticleSnumu.Array_Decays[2][2] = snumuamplitudenumuneutralinoZ1; ParticleSnumu.Array_Decays[2][3] = 2; ParticleSnumu.Array_Comments[2] = "# ~nu_muL -> nu_mu ~chi_10";
   ParticleSnumu.Array_Decays[3][0] = PDGnumuon; ParticleSnumu.Array_Decays[3][1] = PDGneutralino2; ParticleSnumu.Array_Decays[3][2] = snumuamplitudenumuneutralinoZ2; ParticleSnumu.Array_Decays[3][3] = 2; ParticleSnumu.Array_Comments[3] = "# ~nu_muL -> nu_mu ~chi_20";
   ParticleSnumu.Array_Decays[4][0] = PDGnumuon; ParticleSnumu.Array_Decays[4][1] = PDGneutralino3; ParticleSnumu.Array_Decays[4][2] = snumuamplitudenumuneutralinoZ3; ParticleSnumu.Array_Decays[4][3] = 2; ParticleSnumu.Array_Comments[4] = "# ~nu_muL -> nu_mu ~chi_30";
   ParticleSnumu.Array_Decays[5][0] = PDGnumuon; ParticleSnumu.Array_Decays[5][1] = PDGneutralino4; ParticleSnumu.Array_Decays[5][2] = snumuamplitudenumuneutralinoZ4; ParticleSnumu.Array_Decays[5][3] = 2; ParticleSnumu.Array_Comments[5] = "# ~nu_muL -> nu_mu ~chi_40";
   ParticleSnumu.Array_Decays[6][0] = PDGnumuon; ParticleSnumu.Array_Decays[6][1] = PDGneutralino5; ParticleSnumu.Array_Decays[6][2] = snumuamplitudenumuneutralinoZ5; ParticleSnumu.Array_Decays[6][3] = 2; ParticleSnumu.Array_Comments[6] = "# ~nu_muL -> nu_mu ~chi_50";

   ParticleSnumu.Array_Decays[7][0] = PDGnumuon; ParticleSnumu.Array_Decays[7][1] = PDGgravitino; ParticleSnumu.Array_Decays[7][2] = snumuamplitudenumugravitino; ParticleSnumu.Array_Decays[7][3] = 2; ParticleSnumu.Array_Comments[7] = "# ~nu_muL -> nu_mu ~G";

   for(int i = 0; i<ParticleSnumu.No_of_Decays; i++) {
     if (ParticleSnumu.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSnumu.Array_Comments[i] << " is negative = " << ParticleSnumu.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSnumu.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }  
   
   double Snumu_No_1to2_Decays = 0;

   Snumu_No_1to2_Decays = ParticleSnumu.No_1to2_Decays + ParticleSnumu.No_grav_Decays + ParticleSnumu.No_NMSSM_Decays;
 
   for (int j = 0; j<Snumu_No_1to2_Decays; j++) {
     ParticleSnumu.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<Snumu_No_1to2_Decays; j++) {
     ParticleSnumu.two_width = ParticleSnumu.two_width + ParticleSnumu.Array_Decays[j][2];
   }
   for (int j=Snumu_No_1to2_Decays; j<ParticleSnumu.No_of_Decays; j++) {
     ParticleSnumu.three_width = ParticleSnumu.three_width + ParticleSnumu.Array_Decays[j][2];
   }
   ///Note currently no slepton three body decays included - may change in future versions
   if ( ParticleSnumu.three_width != ParticleSnumu.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for snumu - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSnumu.No_of_Decays = Snumu_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSnumu.total_width = ParticleSnumu.two_width;
     }
   else {
     ParticleSnumu.total_width = ParticleSnumu.two_width + ParticleSnumu.three_width;
   }
   
   if ( ParticleSnumu.total_width != ParticleSnumu.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSnumu.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSnumu.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in snumu total width \n");
     }
   
 }


 
 ///Stau1 decays
 double stau1amplitudetauneutralinoZ1=0, stau1amplitudetauneutralinoZ2=0, stau1amplitudetauneutralinoZ3=0, stau1amplitudetauneutralinoZ4=0, stau1amplitudetauneutrinocharginoW1=0, stau1amplitudetauneutrinocharginoW2=0, stau1amplitudetausneutrinoHminus=0, stau1amplitudesnustauWboson=0, stau1amplitudetaugravitino=0;
 double stau1amplitudetauneutralinoZ5 = 0;
 
 if (flagstau1 == 1) {
   if (nmssmIsIt == false) {
     stau1amplitudetauneutralinoZ1 = stauamplitudedecaytauneutralino (me(1,3), mtau, mneut(1), g, gp, runmw, mixNeut, thetatau, beta, 1, 1);
     stau1amplitudetauneutralinoZ2 = stauamplitudedecaytauneutralino (me(1,3), mtau, mneut(2), g, gp, runmw, mixNeut, thetatau, beta, 1, 2);
     stau1amplitudetauneutralinoZ3 = stauamplitudedecaytauneutralino (me(1,3), mtau, mneut(3), g, gp, runmw, mixNeut, thetatau, beta, 1, 3);
     stau1amplitudetauneutralinoZ4 = stauamplitudedecaytauneutralino (me(1,3), mtau, mneut(4), g, gp, runmw, mixNeut, thetatau, beta, 1, 4);
     stau1amplitudetaugravitino = squarkamplitudedecaygravitino(me(1,3), mgravitino, mtau, MPlreduced, gravonoff, slepNLSP);  
 }
   else if (nmssmIsIt == true) {
     stau1amplitudetauneutralinoZ1 = stauamplitudedecaytauneutralinoNMSSM(me(1,3), mtau, mneut(1), g, gp, thetatau, mixNeut, runmtau, runmw, beta, 1, 1);
     stau1amplitudetauneutralinoZ2 = stauamplitudedecaytauneutralinoNMSSM(me(1,3), mtau, mneut(2), g, gp, thetatau, mixNeut, runmtau, runmw, beta, 1, 2);
     stau1amplitudetauneutralinoZ3 = stauamplitudedecaytauneutralinoNMSSM(me(1,3), mtau, mneut(3), g, gp, thetatau, mixNeut, runmtau, runmw, beta, 1, 3);
     stau1amplitudetauneutralinoZ4 = stauamplitudedecaytauneutralinoNMSSM(me(1,3), mtau, mneut(4), g, gp, thetatau, mixNeut, runmtau, runmw, beta, 1, 4);
     stau1amplitudetauneutralinoZ5 = stauamplitudedecaytauneutralinoNMSSM(me(1,3), mtau, mneut(5), g, gp, thetatau, mixNeut, runmtau, runmw, beta, 1, 5);
   }

   stau1amplitudetauneutrinocharginoW1 = stauamplitudedecaynutauchargino (me(1,3), 0, MCH1, g, runmw, thetatau, thetaL2, beta, runmtau, 1, 1);
   stau1amplitudetauneutrinocharginoW2 = stauamplitudedecaynutauchargino (me(1,3), 0, MCH2, g, runmw, thetatau, thetaL2, beta, runmtau, 1, 2);
   stau1amplitudetausneutrinoHminus = stauamplitudedecaysnustauHminus (me(1,3), msnu(3), mHpm, g, runmw, beta, thetatau, runmtau, greekmu, Atau, 1);
   stau1amplitudesnustauWboson = stauamplitudedecaysnustauWboson (me(1,3), msnu(3), polemw, g, thetatau, 1);
   
   ParticleStau1.Array_Decays[0][0] = PDGnutau; ParticleStau1.Array_Decays[0][1] = -PDGchargino1; ParticleStau1.Array_Decays[0][2] = stau1amplitudetauneutrinocharginoW1; ParticleStau1.Array_Decays[0][3] = 2; ParticleStau1.Array_Comments[0] = "# ~tau_1- -> nu_tau ~chi_1-";
   ParticleStau1.Array_Decays[1][0] = PDGnutau; ParticleStau1.Array_Decays[1][1] = -PDGchargino2; ParticleStau1.Array_Decays[1][2] = stau1amplitudetauneutrinocharginoW2; ParticleStau1.Array_Decays[1][3] = 2; ParticleStau1.Array_Comments[1] = "# ~tau_1- -> nu_tau ~chi_2-";
   ParticleStau1.Array_Decays[2][0] = PDGnustauL; ParticleStau1.Array_Decays[2][1] = -PDGHplus; ParticleStau1.Array_Decays[2][2] = stau1amplitudetausneutrinoHminus; ParticleStau1.Array_Decays[2][3] = 2; ParticleStau1.Array_Comments[2] = "# ~tau_1- -> H- ~nu_tauL";
   ParticleStau1.Array_Decays[3][0] = PDGnustauL; ParticleStau1.Array_Decays[3][1] = -PDGWplus; ParticleStau1.Array_Decays[3][2] = stau1amplitudesnustauWboson; ParticleStau1.Array_Decays[3][3] = 2; ParticleStau1.Array_Comments[3] = "# ~tau_1- -> W- ~nu_tauL";
   ParticleStau1.Array_Decays[4][0] = PDGtau; ParticleStau1.Array_Decays[4][1] = PDGneutralino1; ParticleStau1.Array_Decays[4][2] = stau1amplitudetauneutralinoZ1; ParticleStau1.Array_Decays[4][3] = 2; ParticleStau1.Array_Comments[4] = "# ~tau_1- -> tau- ~chi_10";
   ParticleStau1.Array_Decays[5][0] = PDGtau; ParticleStau1.Array_Decays[5][1] = PDGneutralino2; ParticleStau1.Array_Decays[5][2] = stau1amplitudetauneutralinoZ2; ParticleStau1.Array_Decays[5][3] = 2; ParticleStau1.Array_Comments[5] = "# ~tau_1- -> tau- ~chi_20";
   ParticleStau1.Array_Decays[6][0] = PDGtau; ParticleStau1.Array_Decays[6][1] = PDGneutralino3; ParticleStau1.Array_Decays[6][2] = stau1amplitudetauneutralinoZ3; ParticleStau1.Array_Decays[6][3] = 2; ParticleStau1.Array_Comments[6] = "# ~tau_1- -> tau- ~chi_30";
   ParticleStau1.Array_Decays[7][0] = PDGtau; ParticleStau1.Array_Decays[7][1] = PDGneutralino4; ParticleStau1.Array_Decays[7][2] = stau1amplitudetauneutralinoZ4; ParticleStau1.Array_Decays[7][3] = 2; ParticleStau1.Array_Comments[7] = "# ~tau_1- -> tau- ~chi_40";
   ParticleStau1.Array_Decays[8][0] = PDGtau; ParticleStau1.Array_Decays[8][1] = PDGneutralino5; ParticleStau1.Array_Decays[8][2] = stau1amplitudetauneutralinoZ5; ParticleStau1.Array_Decays[8][3] = 2; ParticleStau1.Array_Comments[8] = "# ~tau_1- -> tau- ~chi_50";
 
   ParticleStau1.Array_Decays[9][0] = PDGtau; ParticleStau1.Array_Decays[9][1] = PDGgravitino; ParticleStau1.Array_Decays[9][2] = stau1amplitudetaugravitino; ParticleStau1.Array_Decays[9][3] = 2; ParticleStau1.Array_Comments[9] = "# ~tau_1- -> tau- ~G";

   for(int i = 0; i<ParticleStau1.No_of_Decays; i++) {
     if (ParticleStau1.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleStau1.Array_Comments[i] << " is negative = " << ParticleStau1.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleStau1.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }     
   
   double Stau1_No_1to2_Decays = 0;
   
   Stau1_No_1to2_Decays = ParticleStau1.No_1to2_Decays + ParticleStau1.No_grav_Decays + ParticleStau1.No_NMSSM_Decays;
 
   for (int j = 0; j<Stau1_No_1to2_Decays; j++) {
     ParticleStau1.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<Stau1_No_1to2_Decays; j++) {
     ParticleStau1.two_width = ParticleStau1.two_width + ParticleStau1.Array_Decays[j][2];
   }
   for (int j=Stau1_No_1to2_Decays; j<ParticleStau1.No_of_Decays; j++) {
     ParticleStau1.three_width = ParticleStau1.three_width + ParticleStau1.Array_Decays[j][2];
   }
   ///Note currently no slepton three body decays included - may change in future versions
   if ( ParticleStau1.three_width != ParticleStau1.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for stau1 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleStau1.No_of_Decays = Stau1_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleStau1.total_width = ParticleStau1.two_width;
     }
   else {
     ParticleStau1.total_width = ParticleStau1.two_width + ParticleStau1.three_width;
   }
   if ( ParticleStau1.total_width != ParticleStau1.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleStau1.No_of_Decays; i++) {
       //   fout << i << " " << ParticleStau1.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in stau1 total width \n");
     }
  
 }


///Stau2 decays

 double stau2amplitudetauneutralinoZ1=0, stau2amplitudetauneutralinoZ2=0, stau2amplitudetauneutralinoZ3=0, stau2amplitudetauneutralinoZ4=0, stau2amplitudetauneutrinocharginoW1=0, stau2amplitudetauneutrinocharginoW2=0, stau2amplitudestausneutrinoHminus=0, stau2amplitudestausneutrinoWboson=0, stau2amplitudestau1Zboson=0, stau2amplitudestau1h=0, stau2amplitudestau1H=0, stau2amplitudestau1A=0, stau2amplitudetaugravitino=0;
 double stau2amplitudestau1H3 = 0, stau2amplitudestau1A2 = 0, stau2amplitudetauneutralinoZ5 = 0;
 
 if (flagstau2 == 1) {
   stau2amplitudetauneutrinocharginoW1 = stauamplitudedecaynutauchargino (me(2,3), 0, MCH1, g, runmw, thetatau, thetaL2, beta, runmtau, 2, 1);
   stau2amplitudetauneutrinocharginoW2 = stauamplitudedecaynutauchargino (me(2,3), 0, MCH2, g, runmw, thetatau, thetaL2, beta, runmtau, 2, 2);
   stau2amplitudestausneutrinoHminus = stauamplitudedecaysnustauHminus (me(2,3), msnu(3), mHpm, g, runmw, beta, thetatau, runmtau, greekmu, Atau, 2);
   stau2amplitudestausneutrinoWboson = stauamplitudedecaysnustauWboson (me(2,3), msnu(3), polemw, g, thetatau, 2);
   stau2amplitudestau1Zboson = stau2amplitudedecaystau1Zboson (me(2,3), me(1,3), polemz, g, gp, thetatau);
   if (nmssmIsIt == false) {
     stau2amplitudestau1h = stau2amplitudedecaystau1phi (me(2,3), me(1,3), mh0(1), g, gp, thetatau, beta, alpha, runmw, runmtau, greekmu, Atau, 'h');
     stau2amplitudestau1H = stau2amplitudedecaystau1phi (me(2,3), me(1,3), mh0(2), g, gp, thetatau, beta, alpha, runmw, runmtau, greekmu, Atau, 'H');
     stau2amplitudestau1A = stau2amplitudedecaystau1phi (me(2,3), me(1,3), mA0(1), g, gp, thetatau, beta, alpha, runmw, runmtau, greekmu, Atau, 'A');
     stau2amplitudetauneutralinoZ1 = stauamplitudedecaytauneutralino (me(2,3), mtau, mneut(1), g, gp, runmw, mixNeut, thetatau, beta, 2, 1);
     stau2amplitudetauneutralinoZ2 = stauamplitudedecaytauneutralino (me(2,3), mtau, mneut(2), g, gp, runmw, mixNeut, thetatau, beta, 2, 2);
     stau2amplitudetauneutralinoZ3 = stauamplitudedecaytauneutralino (me(2,3), mtau, mneut(3), g, gp, runmw, mixNeut, thetatau, beta, 2, 3);
     stau2amplitudetauneutralinoZ4 = stauamplitudedecaytauneutralino (me(2,3), mtau, mneut(4), g, gp, runmw, mixNeut, thetatau, beta, 2, 4);
     stau2amplitudetaugravitino = squarkamplitudedecaygravitino(me(2,3), mgravitino, mtau, MPlreduced, gravonoff, slepNLSP);
   }
   else if (nmssmIsIt == true) {
     stau2amplitudestau1h = stau2amplitudedecaystau1CPevenhiggsNMSSM (me(2,3), me(1,3), mh0(1), runmtau, thetatau, CPEMix, beta, runmw, g, gp, Atau, mueff, lam, 1); 
     stau2amplitudestau1H = stau2amplitudedecaystau1CPevenhiggsNMSSM (me(2,3), me(1,3), mh0(2), runmtau, thetatau, CPEMix, beta, runmw, g, gp, Atau, mueff, lam, 2); 
     stau2amplitudestau1H3 = stau2amplitudedecaystau1CPevenhiggsNMSSM (me(2,3), me(1,3), mh0(3), runmtau, thetatau, CPEMix, beta, runmw, g, gp, Atau, mueff, lam, 3); 
     stau2amplitudestau1A = stau2amplitudedecaystau1CPoddhiggsNMSSM (me(2,3), me(1,3), mA0(1), runmtau, thetatau, CPOMix, beta, runmw, g , Atau, mueff, lam, 1); 
     stau2amplitudestau1A2 = stau2amplitudedecaystau1CPoddhiggsNMSSM (me(2,3), me(1,3), mA0(2), runmtau, thetatau, CPOMix, beta, runmw, g , Atau, mueff, lam, 2);
     stau2amplitudetauneutralinoZ1 = stauamplitudedecaytauneutralinoNMSSM(me(2,3), mtau, mneut(1), g, gp, thetatau, mixNeut, runmtau, runmw, beta, 2, 1);
     stau2amplitudetauneutralinoZ2 = stauamplitudedecaytauneutralinoNMSSM(me(2,3), mtau, mneut(2), g, gp, thetatau, mixNeut, runmtau, runmw, beta, 2, 2);
     stau2amplitudetauneutralinoZ3 = stauamplitudedecaytauneutralinoNMSSM(me(2,3), mtau, mneut(3), g, gp, thetatau, mixNeut, runmtau, runmw, beta, 2, 3);
     stau2amplitudetauneutralinoZ4 = stauamplitudedecaytauneutralinoNMSSM(me(2,3), mtau, mneut(4), g, gp, thetatau, mixNeut, runmtau, runmw, beta, 2, 4);
     stau2amplitudetauneutralinoZ5 = stauamplitudedecaytauneutralinoNMSSM(me(2,3), mtau, mneut(5), g, gp, thetatau, mixNeut, runmtau, runmw, beta, 2, 5);
   }
   
   ParticleStau2.Array_Decays[0][0] = PDGnutau; ParticleStau2.Array_Decays[0][1] = -PDGchargino1; ParticleStau2.Array_Decays[0][2] = stau2amplitudetauneutrinocharginoW1; ParticleStau2.Array_Decays[0][3] = 2; ParticleStau2.Array_Comments[0] = "# ~tau_2- -> nu_tau ~chi_1-";
   ParticleStau2.Array_Decays[1][0] = PDGnutau; ParticleStau2.Array_Decays[1][1] = -PDGchargino2; ParticleStau2.Array_Decays[1][2] = stau2amplitudetauneutrinocharginoW2; ParticleStau2.Array_Decays[1][3] = 2; ParticleStau2.Array_Comments[1] = "# ~tau_2- -> nu_tau ~chi_2-";
   ParticleStau2.Array_Decays[2][0] = PDGnustauL; ParticleStau2.Array_Decays[2][1] = -PDGHplus; ParticleStau2.Array_Decays[2][2] = stau2amplitudestausneutrinoHminus; ParticleStau2.Array_Decays[2][3] = 2; ParticleStau2.Array_Comments[2] = "# ~tau_2- -> H- ~nu_tauL";
   ParticleStau2.Array_Decays[3][0] = PDGnustauL; ParticleStau2.Array_Decays[3][1] = -PDGWplus; ParticleStau2.Array_Decays[3][2] = stau2amplitudestausneutrinoWboson; ParticleStau2.Array_Decays[3][3] = 2; ParticleStau2.Array_Comments[3] = "# ~tau_2- -> W- ~nu_tauL";
   ParticleStau2.Array_Decays[4][0] = PDGstau1; ParticleStau2.Array_Decays[4][1] = PDGZboson; ParticleStau2.Array_Decays[4][2] = stau2amplitudestau1Zboson; ParticleStau2.Array_Decays[4][3] = 2; ParticleStau2.Array_Comments[4] = "# ~tau_2- -> Z ~tau_1-";
   ParticleStau2.Array_Decays[5][0] = PDGstau1; ParticleStau2.Array_Decays[5][1] = PDGh0; ParticleStau2.Array_Decays[5][2] = stau2amplitudestau1h; ParticleStau2.Array_Decays[5][3] = 2; ParticleStau2.Array_Comments[5] = "# ~tau_2- -> h ~tau_1-";
   ParticleStau2.Array_Decays[6][0] = PDGstau1; ParticleStau2.Array_Decays[6][1] = PDGH0; ParticleStau2.Array_Decays[6][2] = stau2amplitudestau1H; ParticleStau2.Array_Decays[6][3] = 2; ParticleStau2.Array_Comments[6] = "# ~tau_2- -> H ~tau_1-";
   ParticleStau2.Array_Decays[7][0] = PDGstau1; ParticleStau2.Array_Decays[7][1] = PDGA0; ParticleStau2.Array_Decays[7][2] = stau2amplitudestau1A; ParticleStau2.Array_Decays[7][3] = 2; ParticleStau2.Array_Comments[7] = "# ~tau_2- -> A ~tau_1-";

   ParticleStau2.Array_Decays[8][0] = PDGH3; ParticleStau2.Array_Decays[8][1] = PDGstau1; ParticleStau2.Array_Decays[8][2] = stau2amplitudestau1H3; ParticleStau2.Array_Decays[8][3] = 2; ParticleStau2.Array_Comments[8] = "# ~tau_2 -> H3 ~tau_1-";
   ParticleStau2.Array_Decays[9][0] = PDGA2; ParticleStau2.Array_Decays[9][1] = PDGstau1; ParticleStau2.Array_Decays[9][2] = stau2amplitudestau1A2; ParticleStau2.Array_Decays[9][3] = 2; ParticleStau2.Array_Comments[9] = "# ~tau_2 -> A2 ~tau_1-"; 
   ParticleStau2.Array_Decays[10][0] = PDGtau; ParticleStau2.Array_Decays[10][1] = PDGneutralino1; ParticleStau2.Array_Decays[10][2] = stau2amplitudetauneutralinoZ1; ParticleStau2.Array_Decays[10][3] = 2; ParticleStau2.Array_Comments[10] = "# ~tau_2- -> tau- ~chi_10";
   ParticleStau2.Array_Decays[11][0] = PDGtau; ParticleStau2.Array_Decays[11][1] = PDGneutralino2; ParticleStau2.Array_Decays[11][2] = stau2amplitudetauneutralinoZ2; ParticleStau2.Array_Decays[11][3] = 2; ParticleStau2.Array_Comments[11] = "# ~tau_2- -> tau- ~chi_20";
   ParticleStau2.Array_Decays[12][0] = PDGtau; ParticleStau2.Array_Decays[12][1] = PDGneutralino3; ParticleStau2.Array_Decays[12][2] = stau2amplitudetauneutralinoZ3; ParticleStau2.Array_Decays[12][3] = 2; ParticleStau2.Array_Comments[12] = "# ~tau_2- -> tau- ~chi_30";
   ParticleStau2.Array_Decays[13][0] = PDGtau; ParticleStau2.Array_Decays[13][1] = PDGneutralino4; ParticleStau2.Array_Decays[13][2] = stau2amplitudetauneutralinoZ4; ParticleStau2.Array_Decays[13][3] = 2; ParticleStau2.Array_Comments[13] = "# ~tau_2- -> tau- ~chi_40";
   ParticleStau2.Array_Decays[14][0] = PDGtau; ParticleStau2.Array_Decays[14][1] = PDGneutralino5; ParticleStau2.Array_Decays[14][2] = stau2amplitudetauneutralinoZ5; ParticleStau2.Array_Decays[14][3] = 2; ParticleStau2.Array_Comments[14] = "# ~tau_2- -> tau- ~chi_50";
   
   ParticleStau2.Array_Decays[15][0] = PDGtau; ParticleStau2.Array_Decays[15][1] = PDGgravitino; ParticleStau2.Array_Decays[15][2] = stau2amplitudetaugravitino; ParticleStau2.Array_Decays[15][3] = 2; ParticleStau2.Array_Comments[15] = "# ~tau_2- -> tau- ~G";

   for(int i = 0; i<ParticleStau2.No_of_Decays; i++) {
     if (ParticleStau2.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleStau2.Array_Comments[i] << " is negative = " << ParticleStau2.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleStau2.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   

   double Stau2_No_1to2_Decays = 0;

   Stau2_No_1to2_Decays = ParticleStau2.No_1to2_Decays + ParticleStau2.No_grav_Decays + ParticleStau2.No_NMSSM_Decays;
   
   for (int j = 0; j<Stau2_No_1to2_Decays; j++) {
     ParticleStau2.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<Stau2_No_1to2_Decays; j++) {
     ParticleStau2.two_width = ParticleStau2.two_width + ParticleStau2.Array_Decays[j][2];
   }
   for (int j=Stau2_No_1to2_Decays; j<ParticleStau2.No_of_Decays; j++) {
     ParticleStau2.three_width = ParticleStau2.three_width + ParticleStau2.Array_Decays[j][2];
   }
   ///Note currently no slepton three body decays included - may change in future versions
   if ( ParticleStau2.three_width != ParticleStau2.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for stau2 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleStau2.No_of_Decays = Stau2_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleStau2.total_width = ParticleStau2.two_width;
     }
   else {
     ParticleStau2.total_width = ParticleStau2.two_width + ParticleStau2.three_width;
   }

   if ( ParticleStau2.total_width != ParticleStau2.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleStau2.No_of_Decays; i++) {
       //   fout << i << " " << ParticleStau2.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in stau2 total width \n");
     }
 }
 
 ///Stau sneutrino decays

 double snutauamplitudenutauneutralinoZ1=0, snutauamplitudenutauneutralinoZ2=0, snutauamplitudenutauneutralinoZ3=0, snutauamplitudenutauneutralinoZ4=0, snutauamplitudetaucharginoW1=0, snutauamplitudetaucharginoW2=0, snutauamplitudeWbosonstau1=0, snutauamplitudeWbosonstau2=0, snutauamplitudeHplusstau1=0, snutauamplitudeHplusstau2=0, snutauamplitudenutaugravitino=0;
 double snutauamplitudenutauneutralinoZ5 = 0;
 
 if (flagsnutauL == 1) {
   snutauamplitudetaucharginoW1 = stausneutrinoamplitudedecaytauchargino (msnu(3), mtau, MCH1, g, runmw, beta, thetaL2, thetaR2, 1);
   snutauamplitudetaucharginoW2 = stausneutrinoamplitudedecaytauchargino (msnu(3), mtau, MCH2, g, runmw, beta, thetaL2, thetaR2, 2);
   snutauamplitudeHplusstau1 = stauamplitudedecaysnustauHminus (msnu(3), me(1,3), mHpm, g, runmw, beta, thetatau, runmtau, greekmu, Atau, 1);
   snutauamplitudeHplusstau2 = stauamplitudedecaysnustauHminus (msnu(3), me(2,3), mHpm, g, runmw, beta, thetatau, runmtau, greekmu, Atau, 2);
   snutauamplitudeWbosonstau1 = stauamplitudedecaysnustauWboson (msnu(3), me(1,3), polemw, g, thetatau, 1);
   snutauamplitudeWbosonstau2 = stauamplitudedecaysnustauWboson (msnu(3), me(2,3), polemw, g, thetatau, 2);
   if (nmssmIsIt == false) {
     snutauamplitudenutauneutralinoZ1 = stausneutrinoamplitudedecaytauneutrinoneutralino (msnu(3), 0, mneut(1), g, gp, mixNeut, 1);
     snutauamplitudenutauneutralinoZ2 = stausneutrinoamplitudedecaytauneutrinoneutralino (msnu(3), 0, mneut(2), g, gp, mixNeut, 2);
     snutauamplitudenutauneutralinoZ3 = stausneutrinoamplitudedecaytauneutrinoneutralino (msnu(3), 0, mneut(3), g, gp, mixNeut, 3);
     snutauamplitudenutauneutralinoZ4 = stausneutrinoamplitudedecaytauneutrinoneutralino (msnu(3), 0, mneut(4), g, gp, mixNeut, 4);
     snutauamplitudenutaugravitino = squarkamplitudedecaygravitino(msnu(3), mgravitino, 0, MPlreduced, gravonoff, snuNLSP);
   }
   else if (nmssmIsIt == true) {
     snutauamplitudenutauneutralinoZ1 = snutauamplitudedecaynutauneutralinoNMSSM (msnu(3), mneut(1), g, gp, mixNeut, 1);
     snutauamplitudenutauneutralinoZ2 = snutauamplitudedecaynutauneutralinoNMSSM (msnu(3), mneut(2), g, gp, mixNeut, 2);
     snutauamplitudenutauneutralinoZ3 = snutauamplitudedecaynutauneutralinoNMSSM (msnu(3), mneut(3), g, gp, mixNeut, 3);
     snutauamplitudenutauneutralinoZ4 = snutauamplitudedecaynutauneutralinoNMSSM (msnu(3), mneut(4), g, gp, mixNeut, 4);
     snutauamplitudenutauneutralinoZ5 = snutauamplitudedecaynutauneutralinoNMSSM (msnu(3), mneut(5), g, gp, mixNeut, 5);
   }
   
   ParticleSnutau.Array_Decays[0][0] = PDGnutau; ParticleSnutau.Array_Decays[0][1] = PDGneutralino1; ParticleSnutau.Array_Decays[0][2] = snutauamplitudenutauneutralinoZ1; ParticleSnutau.Array_Decays[0][3] = 2; ParticleSnutau.Array_Comments[0] = "# ~nu_tauL -> nu_tau ~chi_10";
   ParticleSnutau.Array_Decays[1][0] = PDGnutau; ParticleSnutau.Array_Decays[1][1] = PDGneutralino2; ParticleSnutau.Array_Decays[1][2] = snutauamplitudenutauneutralinoZ2; ParticleSnutau.Array_Decays[1][3] = 2; ParticleSnutau.Array_Comments[1] = "# ~nu_tauL -> nu_tau ~chi_20";
   ParticleSnutau.Array_Decays[2][0] = PDGnutau; ParticleSnutau.Array_Decays[2][1] = PDGneutralino3; ParticleSnutau.Array_Decays[2][2] = snutauamplitudenutauneutralinoZ3; ParticleSnutau.Array_Decays[2][3] = 2; ParticleSnutau.Array_Comments[2] = "# ~nu_tauL -> nu_tau ~chi_30";
   ParticleSnutau.Array_Decays[3][0] = PDGnutau; ParticleSnutau.Array_Decays[3][1] = PDGneutralino4; ParticleSnutau.Array_Decays[3][2] = snutauamplitudenutauneutralinoZ4; ParticleSnutau.Array_Decays[3][3] = 2; ParticleSnutau.Array_Comments[3] = "# ~nu_tauL -> nu_tau ~chi_40";
   ParticleSnutau.Array_Decays[4][0] = PDGtau; ParticleSnutau.Array_Decays[4][1] = PDGchargino1; ParticleSnutau.Array_Decays[4][2] = snutauamplitudetaucharginoW1; ParticleSnutau.Array_Decays[4][3] = 2; ParticleSnutau.Array_Comments[4] = "# ~nu_tauL -> tau- ~chi_1+";
   ParticleSnutau.Array_Decays[5][0] = PDGtau; ParticleSnutau.Array_Decays[5][1] = PDGchargino2; ParticleSnutau.Array_Decays[5][2] = snutauamplitudetaucharginoW2; ParticleSnutau.Array_Decays[5][3] = 2; ParticleSnutau.Array_Comments[5] = "# ~nu_tauL -> tau- ~chi_2+";
   ParticleSnutau.Array_Decays[6][0] = PDGstau1; ParticleSnutau.Array_Decays[6][1] = PDGHplus; ParticleSnutau.Array_Decays[6][2] = snutauamplitudeHplusstau1; ParticleSnutau.Array_Decays[6][3] = 2; ParticleSnutau.Array_Comments[6] = "# ~nu_tauL -> H+ ~tau_1-";
   ParticleSnutau.Array_Decays[7][0] = PDGstau2; ParticleSnutau.Array_Decays[7][1] = PDGHplus; ParticleSnutau.Array_Decays[7][2] = snutauamplitudeHplusstau2; ParticleSnutau.Array_Decays[7][3] = 2; ParticleSnutau.Array_Comments[7] = "# ~nu_tauL -> H+ ~tau_2-";
   ParticleSnutau.Array_Decays[8][0] = PDGWplus; ParticleSnutau.Array_Decays[8][1] = PDGstau1; ParticleSnutau.Array_Decays[8][2] = snutauamplitudeWbosonstau1; ParticleSnutau.Array_Decays[8][3] = 2; ParticleSnutau.Array_Comments[8] = "# ~nu_tauL -> W+ ~tau_1-";
   ParticleSnutau.Array_Decays[9][0] = PDGWplus; ParticleSnutau.Array_Decays[9][1] = PDGstau2; ParticleSnutau.Array_Decays[9][2] = snutauamplitudeWbosonstau2; ParticleSnutau.Array_Decays[9][3] = 2; ParticleSnutau.Array_Comments[9] = "# ~nu_tauL -> W+ ~tau_2-";

   ParticleSnutau.Array_Decays[10][0] = PDGnutau; ParticleSnutau.Array_Decays[10][1] = PDGneutralino5; ParticleSnutau.Array_Decays[10][2] = snutauamplitudenutauneutralinoZ5; ParticleSnutau.Array_Decays[10][3] = 2; ParticleSnutau.Array_Comments[10] = "# ~nu_tauL -> nu_tau ~chi_50";

   ParticleSnutau.Array_Decays[11][0] = PDGnutau; ParticleSnutau.Array_Decays[11][1] = PDGgravitino; ParticleSnutau.Array_Decays[11][2] = snutauamplitudenutaugravitino; ParticleSnutau.Array_Decays[11][3] = 2; ParticleSnutau.Array_Comments[11] = "# ~nu_tauL -> nu_tau ~G";

   for(int i = 0; i<ParticleSnutau.No_of_Decays; i++) {
     if (ParticleSnutau.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleSnutau.Array_Comments[i] << " is negative = " << ParticleSnutau.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleSnutau.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   

   double Snutau_No_1to2_Decays = 0;

   Snutau_No_1to2_Decays = ParticleSnutau.No_1to2_Decays + ParticleSnutau.No_grav_Decays + ParticleSnutau.No_NMSSM_Decays;
 
   for (int j = 0; j<Snutau_No_1to2_Decays; j++) {
     ParticleSnutau.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<Snutau_No_1to2_Decays; j++) {
     ParticleSnutau.two_width = ParticleSnutau.two_width + ParticleSnutau.Array_Decays[j][2];
   }
   for (int j=Snutau_No_1to2_Decays; j<ParticleSnutau.No_of_Decays; j++) {
     ParticleSnutau.three_width = ParticleSnutau.three_width + ParticleSnutau.Array_Decays[j][2];
   }
   ///Note currently no slepton three body decays included - may change in future versions
   if ( ParticleSnutau.three_width != ParticleSnutau.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for snutau - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleSnutau.No_of_Decays = Snutau_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleSnutau.total_width = ParticleSnutau.two_width;
     }
   else {
     ParticleSnutau.total_width = ParticleSnutau.two_width + ParticleSnutau.three_width;
   }
   
   if ( ParticleSnutau.total_width != ParticleSnutau.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleSnutau.No_of_Decays; i++) {
       //   fout << i << " " << ParticleSnutau.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in snutau total width \n");
     }
   
 }


 ///Chargino Decays
 ///Chargino1 Decays
 
 double chargino1amplitudeWbosonneutralinoZ1=0, chargino1amplitudeWbosonneutralinoZ2=0, chargino1amplitudeWbosonneutralinoZ3=0, chargino1amplitudeWbosonneutralinoZ4=0, chargino1amplitudeHminusneutralinoZ1=0, chargino1amplitudeHminusneutralinoZ2=0, chargino1amplitudeHminusneutralinoZ3=0, chargino1amplitudeHminusneutralinoZ4=0, chargino1amplitudesupLd=0, chargino1amplitudesdownLu=0, chargino1amplitudescharmLs=0, chargino1amplitudesstrangeLc=0, chargino1amplitudestop1b=0, chargino1amplitudestop2b=0, chargino1amplitudesbottom1t=0, chargino1amplitudesbottom2t=0, chargino1amplitudesnuee=0, chargino1amplitudeselectronLnue=0, chargino1amplitudesnumumu=0, chargino1amplitudesmuonLnumu=0, chargino1amplitudesnutautau=0, chargino1amplitudenutaustau1=0, chargino1amplitudenutaustau2=0, chargino1amplitudeneut1udbar=0, chargino1amplitudeneut1csbar=0, chargino1amplitudeneut1nueebar=0, chargino1amplitudeneut1numumubar=0, chargino1amplitudeneut1nutautaubar=0, chargino1amplitudeneut2udbar=0, chargino1amplitudeneut2csbar=0, chargino1amplitudeneut2nueebar=0, chargino1amplitudeneut2numumubar=0, chargino1amplitudeneut2nutautaubar=0, chargino1amplitudeneut3udbar=0, chargino1amplitudeneut3csbar=0, chargino1amplitudeneut3nueebar=0, chargino1amplitudeneut3numumubar=0, chargino1amplitudeneut3nutautaubar=0, chargino1amplitudeneut4udbar=0, chargino1amplitudeneut4csbar=0, chargino1amplitudeneut4nueebar=0, chargino1amplitudeneut4numumubar=0, chargino1amplitudeneut4nutautaubar=0;
 double chargino1amplitudeHminusneutralinoZ5 = 0, chargino1amplitudeWbosonneutralinoZ5 = 0;

 if (flagchar1 == 1) {
   chargino1amplitudesupLd = charginoamplitudedecayquarksquarkL (MCH1, mdo, mu(1,1), g, thetaR2, 1);
   chargino1amplitudesdownLu = charginoamplitudedecayquarksquarkL (MCH1, mup, md(1,1), g, thetaL2, 1);
   chargino1amplitudescharmLs = charginoamplitudedecayquarksquarkL (MCH1, ms, mu(1,2), g , thetaR2, 1);
   chargino1amplitudesstrangeLc = charginoamplitudedecayquarksquarkL (MCH1, mc, md(1,2), g, thetaL2, 1);
   chargino1amplitudestop1b = charginoamplitudedecayquarksquarkmix (MCH1, mb, mu(1,3), g, thetat, thetaL2, thetaR2, beta, runmt, runmb, runmw, 1, 1, 1);
   chargino1amplitudestop2b = charginoamplitudedecayquarksquarkmix (MCH1, mb, mu(2,3), g, thetat, thetaL2, thetaR2, beta, runmt, runmb,  runmw, 1, 1, 2); 
   chargino1amplitudesbottom1t = charginoamplitudedecayquarksquarkmix (MCH1, mt, md(1,3), g, thetab, thetaL2, thetaR2, beta, runmt, runmb, runmw, 1, 2, 1);
   chargino1amplitudesbottom2t = charginoamplitudedecayquarksquarkmix (MCH1, mt, md(2,3), g, thetab, thetaL2, thetaR2, beta, runmt, runmb, runmw, 1, 2, 2);
   chargino1amplitudesnuee = charginoamplitudedecayleptonsleptonL (MCH1, mel, msnu(1), g, thetaR2, 1);
   chargino1amplitudeselectronLnue = charginoamplitudedecayleptonsleptonL (MCH1, 0, me(1,1), g, thetaL2, 1);
   chargino1amplitudesnumumu = charginoamplitudedecayleptonsleptonL (MCH1, mmu, msnu(2), g, thetaR2, 1);
   chargino1amplitudesmuonLnumu = charginoamplitudedecayleptonsleptonL (MCH1, 0, me(1,2), g, thetaL2, 1); 
   chargino1amplitudesnutautau = charginoamplitudedecaysnutautau (MCH1, mtau, msnu(3), g, thetaL2, thetaR2, beta, runmw, 1);
   chargino1amplitudenutaustau1 = charginoamplitudedecaystaunutau (MCH1, 0, me(1,3), g, thetaL2, thetaR2, thetatau, beta, runmw, runmtau, 1, 1);
   chargino1amplitudenutaustau2 = charginoamplitudedecaystaunutau (MCH1, 0, me(2,3), g, thetaL2, thetaR2, thetatau, beta, runmw, runmtau, 2, 1);
   if (nmssmIsIt == false) {
     chargino1amplitudeWbosonneutralinoZ1 = charginoamplitudedecayWbosonneutralino (MCH1, polemw, mneut(1), g, thetaL2, thetaR2, mixNeut, 1, 1);
     chargino1amplitudeWbosonneutralinoZ2 = charginoamplitudedecayWbosonneutralino (MCH1, polemw, mneut(2), g, thetaL2, thetaR2, mixNeut, 1, 2);
     chargino1amplitudeWbosonneutralinoZ3 = charginoamplitudedecayWbosonneutralino (MCH1, polemw, mneut(3), g, thetaL2, thetaR2, mixNeut, 1, 3);
     chargino1amplitudeWbosonneutralinoZ4 = charginoamplitudedecayWbosonneutralino (MCH1, polemw, mneut(4), g, thetaL2, thetaR2, mixNeut, 1, 4);
     chargino1amplitudeHminusneutralinoZ1 = charginoamplitudedecayHminusneutralino (MCH1, mHpm, mneut(1), g, gp, thetaL2, thetaR2, beta, mixNeut, 1, 1);
     chargino1amplitudeHminusneutralinoZ2 = charginoamplitudedecayHminusneutralino (MCH1, mHpm, mneut(2), g, gp, thetaL2, thetaR2, beta, mixNeut, 1, 2);
     chargino1amplitudeHminusneutralinoZ3 = charginoamplitudedecayHminusneutralino (MCH1, mHpm, mneut(3), g, gp, thetaL2, thetaR2, beta, mixNeut, 1, 3);
     chargino1amplitudeHminusneutralinoZ4 = charginoamplitudedecayHminusneutralino (MCH1, mHpm, mneut(4), g, gp, thetaL2, thetaR2, beta, mixNeut, 1, 4); 
   }
   else if (nmssmIsIt == true) {
     chargino1amplitudeHminusneutralinoZ1 = charginoiamplitudedecayneutralinojHpmNMSSM (MCH1, mneut(1), mHpm, g, gp, thetaL2, thetaR2, beta, mixNeut, lam, 1, 1);
     chargino1amplitudeHminusneutralinoZ2 = charginoiamplitudedecayneutralinojHpmNMSSM (MCH1, mneut(2), mHpm, g, gp, thetaL2, thetaR2, beta, mixNeut, lam, 1, 2);
     chargino1amplitudeHminusneutralinoZ3 = charginoiamplitudedecayneutralinojHpmNMSSM (MCH1, mneut(3), mHpm, g, gp, thetaL2, thetaR2, beta, mixNeut, lam, 1, 3);
     chargino1amplitudeHminusneutralinoZ4 = charginoiamplitudedecayneutralinojHpmNMSSM (MCH1, mneut(4), mHpm, g, gp, thetaL2, thetaR2, beta, mixNeut, lam, 1, 4);
     chargino1amplitudeHminusneutralinoZ5 = charginoiamplitudedecayneutralinojHpmNMSSM (MCH1, mneut(5), mHpm, g, gp, thetaL2, thetaR2, beta, mixNeut, lam, 1, 5);
     chargino1amplitudeWbosonneutralinoZ1 = charginoiamplitudedecayneutralinojWNMSSM (MCH1, mneut(1), polemw, g, gp, thetaL2, thetaR2, mixNeut, 1, 1);
     chargino1amplitudeWbosonneutralinoZ2 = charginoiamplitudedecayneutralinojWNMSSM (MCH1, mneut(2), polemw, g, gp, thetaL2, thetaR2, mixNeut, 1, 2);
     chargino1amplitudeWbosonneutralinoZ3 = charginoiamplitudedecayneutralinojWNMSSM (MCH1, mneut(3), polemw, g, gp, thetaL2, thetaR2, mixNeut, 1, 3);
     chargino1amplitudeWbosonneutralinoZ4 = charginoiamplitudedecayneutralinojWNMSSM (MCH1, mneut(4), polemw, g, gp, thetaL2, thetaR2, mixNeut, 1, 4);
     chargino1amplitudeWbosonneutralinoZ5 = charginoiamplitudedecayneutralinojWNMSSM (MCH1, mneut(5), polemw, g, gp, thetaL2, thetaR2, mixNeut, 1, 5);
   }

   chargino1amplitudeneut1udbar = neutralinoamplitudedecaycharginoffprimebar (MCH1, mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mneut(1), mup, mdo, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 1, 1, onetothree, 'q', 'c');
   chargino1amplitudeneut1csbar = neutralinoamplitudedecaycharginoffprimebar (MCH1, mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mneut(1), mc, ms, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 1, 1, onetothree, 'q', 'c');
   chargino1amplitudeneut1nueebar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(1), 10000000000, me(1,1), me(2,1), polemw, mHpm, mneut(1), 0, mel, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 1, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut1numumubar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(2), 10000000000, me(1,2), me(2,2), polemw, mHpm, mneut(1), 0, mmu, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 1, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut1nutautaubar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(3), 10000000000, me(1,3), me(2,3), polemw, mHpm, mneut(1), 0, mtau, 0, thetatau-PI/2, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 1, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut2udbar = neutralinoamplitudedecaycharginoffprimebar (MCH1, mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mneut(2), mup, mdo, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 2, 1, onetothree, 'q', 'c');
   chargino1amplitudeneut2csbar = neutralinoamplitudedecaycharginoffprimebar (MCH1, mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mneut(2), mc, ms, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 2, 1, onetothree, 'q', 'c');
   chargino1amplitudeneut2nueebar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(1), 10000000000, me(1,1), me(2,1), polemw, mHpm, mneut(2), 0, mel, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 2, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut2numumubar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(2), 10000000000, me(1,2), me(2,2), polemw, mHpm, mneut(2), 0, mmu, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 2, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut2nutautaubar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(3), 10000000000, me(1,3), me(2,3), polemw, mHpm, mneut(2), 0, mtau, 0, thetatau-PI/2, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 2, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut3udbar = neutralinoamplitudedecaycharginoffprimebar (MCH1, mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mneut(3), mup, mdo, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 3, 1, onetothree, 'q', 'c');
   chargino1amplitudeneut3csbar = neutralinoamplitudedecaycharginoffprimebar (MCH1, mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mneut(3), mc, ms, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 3, 1, onetothree, 'q', 'c');
   chargino1amplitudeneut3nueebar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(1), 10000000000, me(1,1), me(2,1), polemw, mHpm, mneut(3), 0, mel, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 3, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut3numumubar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(2), 10000000000, me(1,2), me(2,2), polemw, mHpm, mneut(3), 0, mmu, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 3, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut3nutautaubar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(3), 10000000000, me(1,3), me(2,3), polemw, mHpm, mneut(3), 0, mtau, 0, thetatau-PI/2, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 3, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut4udbar = neutralinoamplitudedecaycharginoffprimebar (MCH1, mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mneut(4), mup, mdo, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 4, 1, onetothree, 'q', 'c');
   chargino1amplitudeneut4csbar = neutralinoamplitudedecaycharginoffprimebar (MCH1, mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mneut(4), mc, ms, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 4, 1, onetothree, 'q', 'c');
   chargino1amplitudeneut4nueebar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(1), 10000000000, me(1,1), me(2,1), polemw, mHpm, mneut(4), 0, mel, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 4, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut4numumubar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(2), 10000000000, me(1,2), me(2,2), polemw, mHpm, mneut(4), 0, mmu, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 4, 1, onetothree, 'l', 'c');
   chargino1amplitudeneut4nutautaubar = neutralinoamplitudedecaycharginoffprimebar (MCH1, msnu(3), 10000000000, me(1,3), me(2,3), polemw, mHpm, mneut(4), 0, mtau, 0, thetatau-PI/2, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 4, 1, onetothree, 'l', 'c');
 
   ParticleChargino1.Array_Decays[0][0] = -PDGdown; ParticleChargino1.Array_Decays[0][1] = PDGsupL; ParticleChargino1.Array_Decays[0][2] = chargino1amplitudesupLd; ParticleChargino1.Array_Decays[0][3] = 2; ParticleChargino1.Array_Comments[0] = "# ~chi_1+ -> db ~u_L";
   ParticleChargino1.Array_Decays[1][0] = PDGup; ParticleChargino1.Array_Decays[1][1] = -PDGsdownL; ParticleChargino1.Array_Decays[1][2] = chargino1amplitudesdownLu; ParticleChargino1.Array_Decays[1][3] = 2; ParticleChargino1.Array_Comments[1] = "# ~chi_1+ -> u ~d_L*"; 
   ParticleChargino1.Array_Decays[2][0] = -PDGstrange; ParticleChargino1.Array_Decays[2][1] = PDGscharmL; ParticleChargino1.Array_Decays[2][2] = chargino1amplitudescharmLs; ParticleChargino1.Array_Decays[2][3] = 2; ParticleChargino1.Array_Comments[2] = "# ~chi_1+ -> sb ~s_L";
   ParticleChargino1.Array_Decays[3][0] = PDGcharm; ParticleChargino1.Array_Decays[3][1] = -PDGsstrangeL; ParticleChargino1.Array_Decays[3][2] = chargino1amplitudesstrangeLc; ParticleChargino1.Array_Decays[3][3] = 2; ParticleChargino1.Array_Comments[3] = "# ~chi_1+ -> c ~s_L*";
   ParticleChargino1.Array_Decays[4][0] = -PDGbottom; ParticleChargino1.Array_Decays[4][1] = PDGstop1; ParticleChargino1.Array_Decays[4][2] = chargino1amplitudestop1b; ParticleChargino1.Array_Decays[4][3] = 2; ParticleChargino1.Array_Comments[4] = "# ~chi_1+ -> bb ~t_1";
   ParticleChargino1.Array_Decays[5][0] = -PDGbottom; ParticleChargino1.Array_Decays[5][1] = PDGstop2; ParticleChargino1.Array_Decays[5][2] = chargino1amplitudestop2b; ParticleChargino1.Array_Decays[5][3] = 2; ParticleChargino1.Array_Comments[5] = "# ~chi_1+ -> bb ~t2";
   ParticleChargino1.Array_Decays[6][0] = PDGtop; ParticleChargino1.Array_Decays[6][1] = PDGsbottom1; ParticleChargino1.Array_Decays[6][2] = chargino1amplitudesbottom1t; ParticleChargino1.Array_Decays[6][3] = 2; ParticleChargino1.Array_Comments[6] = "# ~chi_1+ -> t ~b_1";
   ParticleChargino1.Array_Decays[7][0] = PDGtop; ParticleChargino1.Array_Decays[7][1] = PDGsbottom2; ParticleChargino1.Array_Decays[7][2] = chargino1amplitudesbottom2t; ParticleChargino1.Array_Decays[7][3] = 2; ParticleChargino1.Array_Comments[7] = "# ~chi_1+ -> t ~b_2";
   ParticleChargino1.Array_Decays[8][0] = PDGnuselectronL; ParticleChargino1.Array_Decays[8][1] = -PDGelectron; ParticleChargino1.Array_Decays[8][2] = chargino1amplitudesnuee; ParticleChargino1.Array_Decays[8][3] = 2; ParticleChargino1.Array_Comments[8] = "# ~chi_1+ -> e+ ~nu_eL";
   ParticleChargino1.Array_Decays[9][0] = -PDGselectronL; ParticleChargino1.Array_Decays[9][1] = PDGnuelectron; ParticleChargino1.Array_Decays[9][2] = chargino1amplitudeselectronLnue; ParticleChargino1.Array_Decays[9][3] = 2; ParticleChargino1.Array_Comments[9] = "# ~chi_1+ -> nu_e ~e_L+";
   ParticleChargino1.Array_Decays[10][0] = PDGnusmuonL; ParticleChargino1.Array_Decays[10][1] = -PDGmuon; ParticleChargino1.Array_Decays[10][2] = chargino1amplitudesnumumu; ParticleChargino1.Array_Decays[10][3] = 2; ParticleChargino1.Array_Comments[10] = "# ~chi_1+ -> mu+ ~nu_muL";
   ParticleChargino1.Array_Decays[11][0] = -PDGsmuonL; ParticleChargino1.Array_Decays[11][1] = PDGnumuon; ParticleChargino1.Array_Decays[11][2] = chargino1amplitudesmuonLnumu; ParticleChargino1.Array_Decays[11][3] = 2; ParticleChargino1.Array_Comments[11] = "# ~chi_1+ -> nu_mu ~muL+";
   ParticleChargino1.Array_Decays[12][0] = -PDGtau; ParticleChargino1.Array_Decays[12][1] = PDGnustauL; ParticleChargino1.Array_Decays[12][2] = chargino1amplitudesnutautau; ParticleChargino1.Array_Decays[12][3] = 2; ParticleChargino1.Array_Comments[12] = "# ~chi_1+ -> tau+ ~nu_tauL";
   ParticleChargino1.Array_Decays[13][0] = PDGnutau; ParticleChargino1.Array_Decays[13][1] = -PDGstau1; ParticleChargino1.Array_Decays[13][2] = chargino1amplitudenutaustau1; ParticleChargino1.Array_Decays[13][3] = 2; ParticleChargino1.Array_Comments[13] = "# ~chi_1+ -> nu_tau ~tau_1+";
   ParticleChargino1.Array_Decays[14][0] = PDGnutau; ParticleChargino1.Array_Decays[14][1] = -PDGstau2; ParticleChargino1.Array_Decays[14][2] = chargino1amplitudenutaustau2; ParticleChargino1.Array_Decays[14][3] = 2; ParticleChargino1.Array_Comments[14] = "# ~chi_1+ -> nu_tau ~tau_2+";
   ParticleChargino1.Array_Decays[15][0] = PDGWplus; ParticleChargino1.Array_Decays[15][1] = PDGneutralino1; ParticleChargino1.Array_Decays[15][2] = chargino1amplitudeWbosonneutralinoZ1; ParticleChargino1.Array_Decays[15][3] = 2; ParticleChargino1.Array_Comments[15] = "# ~chi_1+ -> W+ ~chi_10";
   ParticleChargino1.Array_Decays[16][0] = PDGWplus; ParticleChargino1.Array_Decays[16][1] = PDGneutralino2; ParticleChargino1.Array_Decays[16][2] = chargino1amplitudeWbosonneutralinoZ2; ParticleChargino1.Array_Decays[16][3] = 2; ParticleChargino1.Array_Comments[16] = "# ~chi_1+ -> W+ ~chi_20";
   ParticleChargino1.Array_Decays[17][0] = PDGWplus; ParticleChargino1.Array_Decays[17][1] = PDGneutralino3; ParticleChargino1.Array_Decays[17][2] = chargino1amplitudeWbosonneutralinoZ3; ParticleChargino1.Array_Decays[17][3] = 2; ParticleChargino1.Array_Comments[17] = "# ~chi_1+ -> W+ ~chi_30";
   ParticleChargino1.Array_Decays[18][0] = PDGWplus; ParticleChargino1.Array_Decays[18][1] = PDGneutralino4; ParticleChargino1.Array_Decays[18][2] = chargino1amplitudeWbosonneutralinoZ4; ParticleChargino1.Array_Decays[18][3] = 2; ParticleChargino1.Array_Comments[18] = "# ~chi_1+ -> W+ ~chi_40";
   ParticleChargino1.Array_Decays[19][0] = PDGHplus; ParticleChargino1.Array_Decays[19][1] = PDGneutralino1; ParticleChargino1.Array_Decays[19][2] = chargino1amplitudeHminusneutralinoZ1; ParticleChargino1.Array_Decays[19][3] = 2; ParticleChargino1.Array_Comments[19] = "# ~chi_1+ -> H+ ~chi_10";
   ParticleChargino1.Array_Decays[20][0] = PDGHplus; ParticleChargino1.Array_Decays[20][1] = PDGneutralino2; ParticleChargino1.Array_Decays[20][2] = chargino1amplitudeHminusneutralinoZ2; ParticleChargino1.Array_Decays[20][3] = 2; ParticleChargino1.Array_Comments[20] = "# ~chi_1+ -> H+ ~chi_20";
   ParticleChargino1.Array_Decays[21][0] = PDGHplus; ParticleChargino1.Array_Decays[21][1] = PDGneutralino3; ParticleChargino1.Array_Decays[21][2] = chargino1amplitudeHminusneutralinoZ3; ParticleChargino1.Array_Decays[21][3] = 2; ParticleChargino1.Array_Comments[21] = "# ~chi_1+ -> H+ ~chi_30";
   ParticleChargino1.Array_Decays[22][0] = PDGHplus; ParticleChargino1.Array_Decays[22][1] = PDGneutralino4; ParticleChargino1.Array_Decays[22][2] = chargino1amplitudeHminusneutralinoZ4; ParticleChargino1.Array_Decays[22][3] = 2; ParticleChargino1.Array_Comments[22] = "# ~chi_1+ -> H+ ~chi_40";

   /// Somehow, you need to sneak in the pion decay here and re-number the
   /// other decays
   double gammaChi1ToPiNeut1 = charginoToNeutralino1pion(r);
   ParticleChargino1.Array_Decays[23][0] = PDGpiPlus; ParticleChargino1.Array_Decays[23][1] = PDGneutralino1; ParticleChargino1.Array_Decays[23][2] = gammaChi1ToPiNeut1; ParticleChargino1.Array_Decays[23][3] = 2; ParticleChargino1.Array_Comments[23] = "# ~chi_1+ -> pi+ ~chi_10";   
   
   ParticleChargino1.Array_Decays[24][0] = PDGHplus; ParticleChargino1.Array_Decays[24][1] = PDGneutralino5; ParticleChargino1.Array_Decays[24][2] = chargino1amplitudeHminusneutralinoZ5; ParticleChargino1.Array_Decays[24][3] = 2; ParticleChargino1.Array_Comments[24] = "# ~chi_1+ -> H+ ~chi_50";
   ParticleChargino1.Array_Decays[25][0] = PDGWplus; ParticleChargino1.Array_Decays[25][1] = PDGneutralino5; ParticleChargino1.Array_Decays[25][2] = chargino1amplitudeWbosonneutralinoZ5; ParticleChargino1.Array_Decays[25][3] = 2; ParticleChargino1.Array_Comments[25] = "# ~chi_1+ -> W+ ~chi_50";
   
   
   ParticleChargino1.Array_Decays[45][0] = PDGneutralino1; ParticleChargino1.Array_Decays[45][1] = PDGup; ParticleChargino1.Array_Decays[45][4] = -PDGdown; ParticleChargino1.Array_Decays[45][2] = chargino1amplitudeneut1udbar; ParticleChargino1.Array_Decays[45][3] = 3; ParticleChargino1.Array_Comments[45] = "# ~chi_1+ -> chi_10 u dbar";
   ParticleChargino1.Array_Decays[26][0] = PDGneutralino1; ParticleChargino1.Array_Decays[26][1] = PDGcharm; ParticleChargino1.Array_Decays[26][4] = -PDGstrange; ParticleChargino1.Array_Decays[26][2] = chargino1amplitudeneut1csbar; ParticleChargino1.Array_Decays[26][3] = 3; ParticleChargino1.Array_Comments[26] = "# ~chi_1+ -> chi_10 c sbar";
   ParticleChargino1.Array_Decays[27][0] = PDGneutralino1; ParticleChargino1.Array_Decays[27][1] = PDGnuelectron; ParticleChargino1.Array_Decays[27][4] = -PDGelectron; ParticleChargino1.Array_Decays[27][2] = chargino1amplitudeneut1nueebar; ParticleChargino1.Array_Decays[27][3] = 3; ParticleChargino1.Array_Comments[27] = "# ~chi_1+ -> chi_10 nu_e e+";
   ParticleChargino1.Array_Decays[28][0] = PDGneutralino1; ParticleChargino1.Array_Decays[28][1] = PDGnumuon; ParticleChargino1.Array_Decays[28][4] = -PDGmuon; ParticleChargino1.Array_Decays[28][2] = chargino1amplitudeneut1numumubar; ParticleChargino1.Array_Decays[28][3] = 3; ParticleChargino1.Array_Comments[28] = "# ~chi_1+ -> chi_10 nu_mu mu+";
   ParticleChargino1.Array_Decays[29][0] = PDGneutralino1; ParticleChargino1.Array_Decays[29][1] = PDGnutau; ParticleChargino1.Array_Decays[29][4] = -PDGtau; ParticleChargino1.Array_Decays[29][2] = chargino1amplitudeneut1nutautaubar; ParticleChargino1.Array_Decays[29][3] = 3; ParticleChargino1.Array_Comments[29] = "# ~chi_1+ -> chi_10 nu_tau tau+";
   ParticleChargino1.Array_Decays[30][0] = PDGneutralino2; ParticleChargino1.Array_Decays[30][1] = PDGup; ParticleChargino1.Array_Decays[30][4] = -PDGdown; ParticleChargino1.Array_Decays[30][2] = chargino1amplitudeneut2udbar; ParticleChargino1.Array_Decays[30][3] = 3; ParticleChargino1.Array_Comments[30] = "# ~chi_1+ -> chi_20 u dbar";
   ParticleChargino1.Array_Decays[31][0] = PDGneutralino2; ParticleChargino1.Array_Decays[31][1] = PDGcharm; ParticleChargino1.Array_Decays[31][4] = -PDGstrange; ParticleChargino1.Array_Decays[31][2] = chargino1amplitudeneut2csbar; ParticleChargino1.Array_Decays[31][3] = 3; ParticleChargino1.Array_Comments[31] = "# ~chi_1+ -> chi_20 c sbar";
   ParticleChargino1.Array_Decays[32][0] = PDGneutralino2; ParticleChargino1.Array_Decays[32][1] = PDGnuelectron; ParticleChargino1.Array_Decays[32][4] = -PDGelectron; ParticleChargino1.Array_Decays[32][2] = chargino1amplitudeneut2nueebar; ParticleChargino1.Array_Decays[32][3] = 3; ParticleChargino1.Array_Comments[32] = "# ~chi_1+ -> chi_20 nu_e e+";
   ParticleChargino1.Array_Decays[33][0] = PDGneutralino2; ParticleChargino1.Array_Decays[33][1] = PDGnumuon; ParticleChargino1.Array_Decays[33][4] = -PDGmuon; ParticleChargino1.Array_Decays[33][2] = chargino1amplitudeneut2numumubar; ParticleChargino1.Array_Decays[33][3] = 3; ParticleChargino1.Array_Comments[33] = "# ~chi_1+ -> chi_20 nu_mu mu+";
   ParticleChargino1.Array_Decays[34][0] = PDGneutralino2; ParticleChargino1.Array_Decays[34][1] = PDGnutau; ParticleChargino1.Array_Decays[34][4] = -PDGtau; ParticleChargino1.Array_Decays[34][2] = chargino1amplitudeneut2nutautaubar; ParticleChargino1.Array_Decays[34][3] = 3; ParticleChargino1.Array_Comments[34] = "# ~chi_1+ -> chi_20 nu_tau tau+";
   ParticleChargino1.Array_Decays[35][0] = PDGneutralino3; ParticleChargino1.Array_Decays[35][1] = PDGup; ParticleChargino1.Array_Decays[35][4] = -PDGdown; ParticleChargino1.Array_Decays[35][2] = chargino1amplitudeneut3udbar; ParticleChargino1.Array_Decays[35][3] = 3; ParticleChargino1.Array_Comments[35] = "# ~chi_1+ -> chi_30 u dbar";
   ParticleChargino1.Array_Decays[36][0] = PDGneutralino3; ParticleChargino1.Array_Decays[36][1] = PDGcharm; ParticleChargino1.Array_Decays[36][4] = -PDGstrange; ParticleChargino1.Array_Decays[36][2] = chargino1amplitudeneut3csbar; ParticleChargino1.Array_Decays[36][3] = 3; ParticleChargino1.Array_Comments[36] = "# ~chi_1+ -> chi_30 c sbar";
   ParticleChargino1.Array_Decays[37][0] = PDGneutralino3; ParticleChargino1.Array_Decays[37][1] = PDGnuelectron; ParticleChargino1.Array_Decays[37][4] = -PDGelectron; ParticleChargino1.Array_Decays[37][2] = chargino1amplitudeneut3nueebar; ParticleChargino1.Array_Decays[37][3] = 3; ParticleChargino1.Array_Comments[37] = "# ~chi_1+ -> chi_30 nu_e e+";
   ParticleChargino1.Array_Decays[38][0] = PDGneutralino3; ParticleChargino1.Array_Decays[38][1] = PDGnumuon; ParticleChargino1.Array_Decays[38][4] = -PDGmuon; ParticleChargino1.Array_Decays[38][2] = chargino1amplitudeneut3numumubar; ParticleChargino1.Array_Decays[38][3] = 3; ParticleChargino1.Array_Comments[38] = "# ~chi_1+ -> chi_30 nu_mu mu+";
   ParticleChargino1.Array_Decays[39][0] = PDGneutralino3; ParticleChargino1.Array_Decays[39][1] = PDGnutau; ParticleChargino1.Array_Decays[39][4] = -PDGtau; ParticleChargino1.Array_Decays[39][2] = chargino1amplitudeneut3nutautaubar; ParticleChargino1.Array_Decays[39][3] = 3; ParticleChargino1.Array_Comments[39] = "# ~chi_1+ -> chi_30 nu_tau tau+";
   ParticleChargino1.Array_Decays[40][0] = PDGneutralino4; ParticleChargino1.Array_Decays[40][1] = PDGup; ParticleChargino1.Array_Decays[40][4] = -PDGdown; ParticleChargino1.Array_Decays[40][2] = chargino1amplitudeneut4udbar; ParticleChargino1.Array_Decays[40][3] = 3; ParticleChargino1.Array_Comments[40] = "# ~chi_1+ -> chi_40 u dbar";
   ParticleChargino1.Array_Decays[41][0] = PDGneutralino4; ParticleChargino1.Array_Decays[41][1] = PDGcharm; ParticleChargino1.Array_Decays[41][4] = -PDGstrange; ParticleChargino1.Array_Decays[41][2] = chargino1amplitudeneut4csbar; ParticleChargino1.Array_Decays[41][3] = 3; ParticleChargino1.Array_Comments[41] = "# ~chi_1+ -> chi_40 c sbar";
   ParticleChargino1.Array_Decays[42][0] = PDGneutralino4; ParticleChargino1.Array_Decays[42][1] = PDGnuelectron; ParticleChargino1.Array_Decays[42][4] = -PDGelectron; ParticleChargino1.Array_Decays[42][2] = chargino1amplitudeneut4nueebar; ParticleChargino1.Array_Decays[42][3] = 3; ParticleChargino1.Array_Comments[42] = "# ~chi_1+ -> chi_40 nu_e e+";
   ParticleChargino1.Array_Decays[43][0] = PDGneutralino4; ParticleChargino1.Array_Decays[43][1] = PDGnumuon; ParticleChargino1.Array_Decays[43][4] = -PDGmuon; ParticleChargino1.Array_Decays[43][2] = chargino1amplitudeneut4numumubar; ParticleChargino1.Array_Decays[43][3] = 3; ParticleChargino1.Array_Comments[43] = "# ~chi_1+ -> chi_40 nu_mu mu+";
   ParticleChargino1.Array_Decays[44][0] = PDGneutralino4; ParticleChargino1.Array_Decays[44][1] = PDGnutau; ParticleChargino1.Array_Decays[44][4] = -PDGtau; ParticleChargino1.Array_Decays[44][2] = chargino1amplitudeneut4nutautaubar; ParticleChargino1.Array_Decays[44][3] = 3; ParticleChargino1.Array_Comments[44] = "# ~chi_1+ -> chi_40 nu_tau tau+";
   /// Ben: Added this one
   ParticleChargino1.Array_Decays[46][0] = PDGneutralino1; ParticleChargino1.Array_Decays[46][1] = PDGpiPlus; ParticleChargino1.Array_Decays[46][4] = PDGpi0; ParticleChargino1.Array_Decays[46][2] = charginoToNeutralino2pion(r); ParticleChargino1.Array_Decays[46][3] = 3; ParticleChargino1.Array_Comments[46] = "# ~chi_1+ -> chi_10 pi+ pi0";
  /// Ben: Added this one
   /*   ParticleChargino1.Array_Decays[47][0] = PDGneutralino2; ParticleChargino1.Array_Decays[47][1] = PDGpiPlus; ParticleChargino1.Array_Decays[47][4] = PDGpi0; ParticleChargino1.Array_Decays[47][2] = charginoToNeutralino21pion(r); ParticleChargino1.Array_Decays[47][3] = 3; ParticleChargino1.Array_Comments[47] = "# ~chi_1+ -> chi_20 pi+ pi0";*/

   /// If at too low mass difference, use pions rather than quarks
   if (fabs(MCH1) - fabs(mneut(1)) < hadronicScale) {
     ParticleChargino1.Array_Decays[26][2] = 0.;
     ParticleChargino1.Array_Decays[45][2] = 0.;
   } else {
     ParticleChargino1.Array_Decays[23][2] = 0.;
     ParticleChargino1.Array_Decays[46][2] = 0.;
   }

   for(int i = 0; i<ParticleChargino1.No_of_Decays; i++) {
     if (ParticleChargino1.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleChargino1.Array_Comments[i] << " is negative = " << ParticleChargino1.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleChargino1.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }

   int Chargino1_No_1to2_Decays = 0;
 
   Chargino1_No_1to2_Decays = ParticleChargino1.No_1to2_Decays + ParticleChargino1.No_grav_Decays + ParticleChargino1.No_NMSSM_Decays;
 
   for (int j = 0; j<Chargino1_No_1to2_Decays; j++) {
     ParticleChargino1.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<Chargino1_No_1to2_Decays; j++) {
     ParticleChargino1.two_width = ParticleChargino1.two_width + ParticleChargino1.Array_Decays[j][2];
   }
   for (int j=Chargino1_No_1to2_Decays; j<ParticleChargino1.No_of_Decays; j++) {
     ParticleChargino1.three_width = ParticleChargino1.three_width + ParticleChargino1.Array_Decays[j][2];
   }
   if ( ParticleChargino1.three_width != ParticleChargino1.three_width) /// Tests for a nan as only nans aren't equal to themselves
  {
    fout << "# Three body decays give nan for chargino1 - warning! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
     errorflag = -1;
     ParticleChargino1.No_of_Decays = Chargino1_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
     ParticleChargino1.total_width = ParticleChargino1.two_width;
   }
   else {
     ParticleChargino1.total_width = ParticleChargino1.two_width + ParticleChargino1.three_width;
   }

   if ( ParticleChargino1.total_width != ParticleChargino1.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleChargino1.No_of_Decays; i++) {
       //   fout << i << " " << ParticleChargino1.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in chargino1 total width \n");
     }
 }

 ///Chargino2 Decays
 double chargino2amplitudeWbosonneutralinoZ1=0, chargino2amplitudeWbosonneutralinoZ2=0, chargino2amplitudeWbosonneutralinoZ3=0, chargino2amplitudeWbosonneutralinoZ4=0, chargino2amplitudeHminusneutralinoZ1=0, chargino2amplitudeHminusneutralinoZ2=0, chargino2amplitudeHminusneutralinoZ3=0, chargino2amplitudeHminusneutralinoZ4=0, chargino2amplitudesupLd=0, chargino2amplitudesdownLu=0, chargino2amplitudescharmLs=0, chargino2amplitudesstrangeLc=0, chargino2amplitudestop1b=0, chargino2amplitudestop2b=0, chargino2amplitudesbottom1t=0, chargino2amplitudesbottom2t=0, chargino2amplitudesnuee=0, chargino2amplitudeselectronLnue=0, chargino2amplitudesnumumu=0, chargino2amplitudesmuonLnumu=0, chargino2amplitudesnutautau=0, chargino2amplitudenutaustau1=0, chargino2amplitudenutaustau2=0, chargino2amplitudechargino1Zboson=0, chargino2amplitudechargino1h=0, chargino2amplitudechargino1H=0, chargino2amplitudechargino1A=0, chargino2amplitudeneut1udbar=0, chargino2amplitudeneut1csbar=0, chargino2amplitudeneut1nueebar=0, chargino2amplitudeneut1numumubar=0, chargino2amplitudeneut1nutautaubar=0, chargino2amplitudeneut2udbar=0, chargino2amplitudeneut2csbar=0, chargino2amplitudeneut2nueebar=0, chargino2amplitudeneut2numumubar=0, chargino2amplitudeneut2nutautaubar=0, chargino2amplitudeneut3udbar=0, chargino2amplitudeneut3csbar=0, chargino2amplitudeneut3nueebar=0, chargino2amplitudeneut3numumubar=0, chargino2amplitudeneut3nutautaubar=0, chargino2amplitudeneut4udbar=0, chargino2amplitudeneut4csbar=0, chargino2amplitudeneut4nueebar=0, chargino2amplitudeneut4numumubar=0, chargino2amplitudeneut4nutautaubar=0;

 double chargino2amplitudechargino1H3 = 0, chargino2amplitudechargino1A2 = 0, chargino2amplitudeHminusneutralinoZ5 = 0, chargino2amplitudeWbosonneutralinoZ5 = 0;

 if (flagchar2 == 1) {
   chargino2amplitudesupLd = charginoamplitudedecayquarksquarkL (MCH2, mdo, mu(1,1), g, thetaR2, 2);
   chargino2amplitudesdownLu = charginoamplitudedecayquarksquarkL (MCH2, mup, md(1,1), g, thetaL2, 2);
   chargino2amplitudescharmLs = charginoamplitudedecayquarksquarkL (MCH2, ms, mu(1,2), g , thetaR2, 2);
   chargino2amplitudesstrangeLc = charginoamplitudedecayquarksquarkL (MCH2, mc, md(1,2), g, thetaL2, 2);
   chargino2amplitudestop1b = charginoamplitudedecayquarksquarkmix (MCH2, mb, mu(1,3), g, thetat, thetaL2, thetaR2, beta, runmt, runmb, runmw, 2, 1, 1);
   chargino2amplitudestop2b = charginoamplitudedecayquarksquarkmix (MCH2, mb, mu(2,3), g, thetat, thetaL2, thetaR2, beta, runmt, runmb, runmw, 2, 1, 2);
   chargino2amplitudesbottom1t = charginoamplitudedecayquarksquarkmix (MCH2, mt, md(1,3), g, thetab, thetaL2, thetaR2, beta, runmt, runmb, runmw, 2, 2, 1);
   chargino2amplitudesbottom2t = charginoamplitudedecayquarksquarkmix (MCH2, mt, md(2,3), g, thetab, thetaL2, thetaR2, beta, runmt, runmb, runmw, 2, 2, 2);
   chargino2amplitudesnuee = charginoamplitudedecayleptonsleptonL (MCH2, mel, msnu(1), g, thetaR2, 2);
   chargino2amplitudeselectronLnue = charginoamplitudedecayleptonsleptonL (MCH2, 0, me(1,1), g, thetaL2, 2);
   chargino2amplitudesnumumu = charginoamplitudedecayleptonsleptonL (MCH2, mmu, msnu(2), g, thetaR2, 2);
   chargino2amplitudesmuonLnumu = charginoamplitudedecayleptonsleptonL (MCH2, 0, me(1,2), g, thetaL2, 2); 
   chargino2amplitudesnutautau = charginoamplitudedecaysnutautau (MCH2, mtau, msnu(3), g, thetaL2, thetaR2, beta, runmw, 2);
   chargino2amplitudenutaustau1 = charginoamplitudedecaystaunutau (MCH2, 0, me(1,3), g, thetaL2, thetaR2, thetatau, beta, runmw, runmtau, 1, 2);
   chargino2amplitudenutaustau2 = charginoamplitudedecaystaunutau (MCH2, 0, me(2,3), g, thetaL2, thetaR2, thetatau, beta, runmw, runmtau, 2, 2);
   chargino2amplitudechargino1Zboson = chargino2amplitudedecaychargino1Zboson (MCH2, polemz, MCH1, g, gp, thetaL2, thetaR2);
   if (nmssmIsIt == false) {
     chargino2amplitudechargino1h = chargino2amplitudedecaychargino1neutHiggs (MCH2, mh0(1), MCH1, g, gp, thetaL2, thetaR2, beta, alpha, 'h');
     chargino2amplitudechargino1H = chargino2amplitudedecaychargino1neutHiggs (MCH2, mh0(2), MCH1, g, gp, thetaL2, thetaR2, beta, alpha, 'H');
     chargino2amplitudechargino1A = chargino2amplitudedecaychargino1neutHiggs (MCH2, mA0(1), MCH1, g, gp, thetaL2, thetaR2, beta, alpha, 'A');
     chargino2amplitudeWbosonneutralinoZ1 = charginoamplitudedecayWbosonneutralino (MCH2, polemw, mneut(1), g, thetaL2, thetaR2, mixNeut, 2, 1);
     chargino2amplitudeWbosonneutralinoZ2 = charginoamplitudedecayWbosonneutralino (MCH2, polemw, mneut(2), g, thetaL2, thetaR2, mixNeut, 2, 2);
     chargino2amplitudeWbosonneutralinoZ3 = charginoamplitudedecayWbosonneutralino (MCH2, polemw, mneut(3), g, thetaL2, thetaR2, mixNeut, 2, 3);
     chargino2amplitudeWbosonneutralinoZ4 = charginoamplitudedecayWbosonneutralino (MCH2, polemw, mneut(4), g, thetaL2, thetaR2, mixNeut, 2, 4);
     chargino2amplitudeHminusneutralinoZ1 = charginoamplitudedecayHminusneutralino (MCH2, mHpm, mneut(1), g, gp, thetaL2, thetaR2, beta, mixNeut, 2, 1);
     chargino2amplitudeHminusneutralinoZ2 = charginoamplitudedecayHminusneutralino (MCH2, mHpm, mneut(2), g, gp, thetaL2, thetaR2, beta, mixNeut, 2, 2);
     chargino2amplitudeHminusneutralinoZ3 = charginoamplitudedecayHminusneutralino (MCH2, mHpm, mneut(3), g, gp, thetaL2, thetaR2, beta, mixNeut, 2, 3);
     chargino2amplitudeHminusneutralinoZ4 = charginoamplitudedecayHminusneutralino (MCH2, mHpm, mneut(4), g, gp, thetaL2, thetaR2, beta, mixNeut, 2, 4);
   }
   else if (nmssmIsIt == true) {
     chargino2amplitudechargino1h = chargino2amplitudedecaychargino1CPevenhiggsNMSSM (MCH2, MCH1, mh0(1), g, lam, thetaL2, thetaR2, CPEMix, 1);
     chargino2amplitudechargino1H = chargino2amplitudedecaychargino1CPevenhiggsNMSSM (MCH2, MCH1, mh0(2), g, lam, thetaL2, thetaR2, CPEMix, 2);
     chargino2amplitudechargino1H3 = chargino2amplitudedecaychargino1CPevenhiggsNMSSM (MCH2, MCH1, mh0(3), g, lam, thetaL2, thetaR2, CPEMix, 3);
     chargino2amplitudechargino1A = chargino2amplitudedecaychargino1CPoddhiggsNMSSM (MCH2, MCH1, mA0(1), g, lam, thetaL2, thetaR2, CPOMix, 1);
     chargino2amplitudechargino1A2 = chargino2amplitudedecaychargino1CPoddhiggsNMSSM (MCH2, MCH1, mA0(2), g, lam, thetaL2, thetaR2, CPOMix, 2); 
     chargino2amplitudeHminusneutralinoZ1 = charginoiamplitudedecayneutralinojHpmNMSSM (MCH2, mneut(1), mHpm, g, gp, thetaL2, thetaR2, beta, mixNeut, lam, 2, 1);
     chargino2amplitudeHminusneutralinoZ2 = charginoiamplitudedecayneutralinojHpmNMSSM (MCH2, mneut(2), mHpm, g, gp, thetaL2, thetaR2, beta, mixNeut, lam, 2, 2);
     chargino2amplitudeHminusneutralinoZ3 = charginoiamplitudedecayneutralinojHpmNMSSM (MCH2, mneut(3), mHpm, g, gp, thetaL2, thetaR2, beta, mixNeut, lam, 2, 3);
     chargino2amplitudeHminusneutralinoZ4 = charginoiamplitudedecayneutralinojHpmNMSSM (MCH2, mneut(4), mHpm, g, gp, thetaL2, thetaR2, beta, mixNeut, lam, 2, 4);
     chargino2amplitudeHminusneutralinoZ5 = charginoiamplitudedecayneutralinojHpmNMSSM (MCH2, mneut(5), mHpm, g, gp, thetaL2, thetaR2, beta, mixNeut, lam, 2, 5);
     chargino2amplitudeWbosonneutralinoZ1 = charginoiamplitudedecayneutralinojWNMSSM (MCH2, mneut(1), polemw, g, gp, thetaL2, thetaR2, mixNeut, 2, 1);
     chargino2amplitudeWbosonneutralinoZ2 = charginoiamplitudedecayneutralinojWNMSSM (MCH2, mneut(2), polemw, g, gp, thetaL2, thetaR2, mixNeut, 2, 2);
     chargino2amplitudeWbosonneutralinoZ3 = charginoiamplitudedecayneutralinojWNMSSM (MCH2, mneut(3), polemw, g, gp, thetaL2, thetaR2, mixNeut, 2, 3);
     chargino2amplitudeWbosonneutralinoZ4 = charginoiamplitudedecayneutralinojWNMSSM (MCH2, mneut(4), polemw, g, gp, thetaL2, thetaR2, mixNeut, 2, 4);
     chargino2amplitudeWbosonneutralinoZ5 = charginoiamplitudedecayneutralinojWNMSSM (MCH2, mneut(5), polemw, g, gp, thetaL2, thetaR2, mixNeut, 2, 5);
   }

   chargino2amplitudeneut1udbar = neutralinoamplitudedecaycharginoffprimebar (MCH2, mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mneut(1), mup, mdo, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 1, 2, onetothree, 'q', 'c');
   chargino2amplitudeneut1csbar = neutralinoamplitudedecaycharginoffprimebar (MCH2, mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mneut(1), mc, ms, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 1, 2, onetothree, 'q', 'c');
   chargino2amplitudeneut1nueebar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(1), 0, me(1,1), me(2,1), polemw, mHpm, mneut(1), 0, mel, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 1, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut1numumubar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(2), 0, me(1,2), me(2,2), polemw, mHpm, mneut(1), 0, mmu, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 1, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut1nutautaubar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(3), 0, me(1,3), me(2,3), polemw, mHpm, mneut(1), 0, mtau, 0, thetatau-PI/2, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 1, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut2udbar = neutralinoamplitudedecaycharginoffprimebar (MCH2, mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mneut(2), mup, mdo, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 2, 2, onetothree, 'q', 'c');
   chargino2amplitudeneut2csbar = neutralinoamplitudedecaycharginoffprimebar (MCH2, mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mneut(2), mc, ms, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 2, 2, onetothree, 'q', 'c');
   chargino2amplitudeneut2nueebar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(1), 0, me(1,1), me(2,1), polemw, mHpm, mneut(2), 0, mel, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 2, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut2numumubar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(2), 0, me(1,2), me(2,2), polemw, mHpm, mneut(2), 0, mmu, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 2, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut2nutautaubar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(3), 0, me(1,3), me(2,3), polemw, mHpm, mneut(2), 0, mtau, 0, thetatau-PI/2, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 2, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut3udbar = neutralinoamplitudedecaycharginoffprimebar (MCH2, mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mneut(3), mup, mdo, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 3, 2, onetothree, 'q', 'c');
   chargino2amplitudeneut3csbar = neutralinoamplitudedecaycharginoffprimebar (MCH2, mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mneut(3), mc, ms, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 3, 2, onetothree, 'q', 'c');
   chargino2amplitudeneut3nueebar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(1), 0, me(1,1), me(2,1), polemw, mHpm, mneut(3), 0, mel, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 3, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut3numumubar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(2), 0, me(1,2), me(2,2), polemw, mHpm, mneut(3), 0, mmu, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 3, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut3nutautaubar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(3), 0, me(1,3), me(2,3), polemw, mHpm, mneut(3), 0, mtau, 0, thetatau-PI/2, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 3, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut4udbar = neutralinoamplitudedecaycharginoffprimebar (MCH2, mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mneut(4), mup, mdo, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 4, 2, onetothree, 'q', 'c');
   chargino2amplitudeneut4csbar = neutralinoamplitudedecaycharginoffprimebar (MCH2, mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mneut(4), mc, ms, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 4, 2, onetothree, 'q', 'c');
   chargino2amplitudeneut4nueebar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(1), 0, me(1,1), me(2,1), polemw, mHpm, mneut(4), 0, mel, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 4, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut4numumubar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(2), 0, me(1,2), me(2,2), polemw, mHpm, mneut(4), 0, mmu, 0, 0, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 4, 2, onetothree, 'l', 'c');
   chargino2amplitudeneut4nutautaubar = neutralinoamplitudedecaycharginoffprimebar (MCH2, msnu(3), 0, me(1,3), me(2,3), polemw, mHpm, mneut(4), 0, mtau, 0, thetatau-PI/2, g, gp, alpha, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 4, 2, onetothree, 'l', 'c');
 
   ParticleChargino2.Array_Decays[0][0] = -PDGdown; ParticleChargino2.Array_Decays[0][1] = PDGsupL; ParticleChargino2.Array_Decays[0][2] = chargino2amplitudesupLd; ParticleChargino2.Array_Decays[0][3] = 2; ParticleChargino2.Array_Comments[0] = "# ~chi_2+ -> db ~u_L";
   ParticleChargino2.Array_Decays[1][0] = PDGup; ParticleChargino2.Array_Decays[1][1] = -PDGsdownL; ParticleChargino2.Array_Decays[1][2] = chargino2amplitudesdownLu; ParticleChargino2.Array_Decays[1][3] = 2; ParticleChargino2.Array_Comments[1] = "# ~chi_2+ -> u ~d_L*";
   ParticleChargino2.Array_Decays[2][0] = -PDGstrange; ParticleChargino2.Array_Decays[2][1] = PDGscharmL; ParticleChargino2.Array_Decays[2][2] = chargino2amplitudescharmLs; ParticleChargino2.Array_Decays[2][3] = 2; ParticleChargino2.Array_Comments[2] = "# ~chi_2+ -> sb ~c_L";
   ParticleChargino2.Array_Decays[3][0] = PDGcharm; ParticleChargino2.Array_Decays[3][1] = -PDGsstrangeL; ParticleChargino2.Array_Decays[3][2] = chargino2amplitudesstrangeLc; ParticleChargino2.Array_Decays[3][3] = 2; ParticleChargino2.Array_Comments[3] = "# ~chi_2+ -> c ~s_L*";
   ParticleChargino2.Array_Decays[4][0] = -PDGbottom; ParticleChargino2.Array_Decays[4][1] = PDGstop1; ParticleChargino2.Array_Decays[4][2] = chargino2amplitudestop1b; ParticleChargino2.Array_Decays[4][3] = 2; ParticleChargino2.Array_Comments[4] = "# ~chi_2+ -> bb ~t_1";
   ParticleChargino2.Array_Decays[5][0] = -PDGbottom; ParticleChargino2.Array_Decays[5][1] = PDGstop2; ParticleChargino2.Array_Decays[5][2] = chargino2amplitudestop2b; ParticleChargino2.Array_Decays[5][3] = 2; ParticleChargino2.Array_Comments[5] = "# ~chi_2+ -> bb ~t_2";
   ParticleChargino2.Array_Decays[6][0] = PDGtop; ParticleChargino2.Array_Decays[6][1] = -PDGsbottom1; ParticleChargino2.Array_Decays[6][2] = chargino2amplitudesbottom1t; ParticleChargino2.Array_Decays[6][3] = 2; ParticleChargino2.Array_Comments[6] = "# ~chi_2+ -> t ~b_1*";
   ParticleChargino2.Array_Decays[7][0] = PDGtop; ParticleChargino2.Array_Decays[7][1] = -PDGsbottom2; ParticleChargino2.Array_Decays[7][2] = chargino2amplitudesbottom2t; ParticleChargino2.Array_Decays[7][3] = 2; ParticleChargino2.Array_Comments[7] = "# ~chi_2+ -> t ~b_2*";
   ParticleChargino2.Array_Decays[8][0] = PDGnuselectronL; ParticleChargino2.Array_Decays[8][1] = -PDGelectron; ParticleChargino2.Array_Decays[8][2] = chargino2amplitudesnuee; ParticleChargino2.Array_Decays[8][3] = 2; ParticleChargino2.Array_Comments[8] = "# ~chi_2+ -> e+ ~nu_eL";
   ParticleChargino2.Array_Decays[9][0] = -PDGselectronL; ParticleChargino2.Array_Decays[9][1] = PDGnuelectron; ParticleChargino2.Array_Decays[9][2] = chargino2amplitudeselectronLnue; ParticleChargino2.Array_Decays[9][3] = 2; ParticleChargino2.Array_Comments[9] = "# ~chi_2+ -> nu_e ~e_L+";
   ParticleChargino2.Array_Decays[10][0] = PDGnusmuonL; ParticleChargino2.Array_Decays[10][1] = -PDGmuon; ParticleChargino2.Array_Decays[10][2] = chargino2amplitudesnumumu; ParticleChargino2.Array_Decays[10][3] = 2; ParticleChargino2.Array_Comments[10] = "# ~chi_2+ -> mu+ ~nu_muL";
   ParticleChargino2.Array_Decays[11][0] = -PDGsmuonL; ParticleChargino2.Array_Decays[11][1] = PDGnumuon; ParticleChargino2.Array_Decays[11][2] = chargino2amplitudesmuonLnumu; ParticleChargino2.Array_Decays[11][3] = 2; ParticleChargino2.Array_Comments[11] = "# ~chi_2+ -> nu_mu ~muL+";
   ParticleChargino2.Array_Decays[12][0] = -PDGtau; ParticleChargino2.Array_Decays[12][1] = PDGnustauL; ParticleChargino2.Array_Decays[12][2] = chargino2amplitudesnutautau; ParticleChargino2.Array_Decays[12][3] = 2; ParticleChargino2.Array_Comments[12] = "# ~chi_2+ -> tau+ ~nu_tauL";
   ParticleChargino2.Array_Decays[13][0] = PDGnutau; ParticleChargino2.Array_Decays[13][1] = -PDGstau1; ParticleChargino2.Array_Decays[13][2] = chargino2amplitudenutaustau1; ParticleChargino2.Array_Decays[13][3] = 2; ParticleChargino2.Array_Comments[13] = "# ~chi_2+ -> nu_tau ~tau_1+";
   ParticleChargino2.Array_Decays[14][0] = PDGnutau; ParticleChargino2.Array_Decays[14][1] = -PDGstau2; ParticleChargino2.Array_Decays[14][2] = chargino2amplitudenutaustau2; ParticleChargino2.Array_Decays[14][3] = 2; ParticleChargino2.Array_Comments[14] = "# ~chi_2+ -> nu_tau ~tau_2+";
   ParticleChargino2.Array_Decays[15][0] = PDGWplus; ParticleChargino2.Array_Decays[15][1] = PDGneutralino1; ParticleChargino2.Array_Decays[15][2] = chargino2amplitudeWbosonneutralinoZ1; ParticleChargino2.Array_Decays[15][3] = 2; ParticleChargino2.Array_Comments[15] = "# ~chi_2+ -> W+ ~chi_10";
   ParticleChargino2.Array_Decays[16][0] = PDGWplus; ParticleChargino2.Array_Decays[16][1] = PDGneutralino2; ParticleChargino2.Array_Decays[16][2] = chargino2amplitudeWbosonneutralinoZ2; ParticleChargino2.Array_Decays[16][3] = 2; ParticleChargino2.Array_Comments[16] = "# ~chi_2+ -> W+ ~chi_20";
   ParticleChargino2.Array_Decays[17][0] = PDGWplus; ParticleChargino2.Array_Decays[17][1] = PDGneutralino3; ParticleChargino2.Array_Decays[17][2] = chargino2amplitudeWbosonneutralinoZ3; ParticleChargino2.Array_Decays[17][3] = 2; ParticleChargino2.Array_Comments[17] = "# ~chi_2+ -> W+ ~chi_30";
   ParticleChargino2.Array_Decays[18][0] = PDGWplus; ParticleChargino2.Array_Decays[18][1] = PDGneutralino4; ParticleChargino2.Array_Decays[18][2] = chargino2amplitudeWbosonneutralinoZ4; ParticleChargino2.Array_Decays[18][3] = 2; ParticleChargino2.Array_Comments[18] = "# ~chi_2+ -> W+ ~chi_40";
   ParticleChargino2.Array_Decays[19][0] = PDGHplus; ParticleChargino2.Array_Decays[19][1] = PDGneutralino1; ParticleChargino2.Array_Decays[19][2] = chargino2amplitudeHminusneutralinoZ1; ParticleChargino2.Array_Decays[19][3] = 2; ParticleChargino2.Array_Comments[19] = "# ~chi_2+ -> H+ ~chi_10";
   ParticleChargino2.Array_Decays[20][0] = PDGHplus; ParticleChargino2.Array_Decays[20][1] = PDGneutralino2; ParticleChargino2.Array_Decays[20][2] = chargino2amplitudeHminusneutralinoZ2; ParticleChargino2.Array_Decays[20][3] = 2; ParticleChargino2.Array_Comments[20] = "# ~chi_2+ -> H+ ~chi_20";
   ParticleChargino2.Array_Decays[21][0] = PDGHplus; ParticleChargino2.Array_Decays[21][1] = PDGneutralino3; ParticleChargino2.Array_Decays[21][2] = chargino2amplitudeHminusneutralinoZ3; ParticleChargino2.Array_Decays[21][3] = 2; ParticleChargino2.Array_Comments[21] = "# ~chi_2+ -> H+ ~chi_30";
   ParticleChargino2.Array_Decays[22][0] = PDGHplus; ParticleChargino2.Array_Decays[22][1] = PDGneutralino4; ParticleChargino2.Array_Decays[22][2] = chargino2amplitudeHminusneutralinoZ4; ParticleChargino2.Array_Decays[22][3] = 2; ParticleChargino2.Array_Comments[22] = "# ~chi_2+ -> H+ ~chi_40";
   ParticleChargino2.Array_Decays[23][0] = PDGchargino1; ParticleChargino2.Array_Decays[23][1] = PDGZboson; ParticleChargino2.Array_Decays[23][2] = chargino2amplitudechargino1Zboson; ParticleChargino2.Array_Decays[23][3] = 2; ParticleChargino2.Array_Comments[23] = "# ~chi_2+ -> Z ~chi_1+";
   ParticleChargino2.Array_Decays[24][0] = PDGchargino1; ParticleChargino2.Array_Decays[24][1] = PDGh0; ParticleChargino2.Array_Decays[24][2] = chargino2amplitudechargino1h; ParticleChargino2.Array_Decays[24][3] = 2; ParticleChargino2.Array_Comments[24] = "# ~chi_2+ -> h ~chi_1+";
   ParticleChargino2.Array_Decays[25][0] = PDGchargino1; ParticleChargino2.Array_Decays[25][1] = PDGH0; ParticleChargino2.Array_Decays[25][2] = chargino2amplitudechargino1H; ParticleChargino2.Array_Decays[25][3] = 2; ParticleChargino2.Array_Comments[25] = "# ~chi_2+ -> H ~chi_1+";
   ParticleChargino2.Array_Decays[26][0] = PDGchargino1; ParticleChargino2.Array_Decays[26][1] = PDGA0; ParticleChargino2.Array_Decays[26][2] = chargino2amplitudechargino1A; ParticleChargino2.Array_Decays[26][3] = 2; ParticleChargino2.Array_Comments[26] = "# ~chi_2+ -> A ~chi_1+";

   ParticleChargino2.Array_Decays[27][0] = PDGHplus; ParticleChargino2.Array_Decays[27][1] = PDGneutralino5; ParticleChargino2.Array_Decays[27][2] = chargino2amplitudeHminusneutralinoZ5; ParticleChargino2.Array_Decays[27][3] = 2; ParticleChargino2.Array_Comments[27] = "# ~chi_2+ -> H+ ~chi_50";
   ParticleChargino2.Array_Decays[28][0] = PDGWplus; ParticleChargino2.Array_Decays[28][1] = PDGneutralino5; ParticleChargino2.Array_Decays[28][2] = chargino2amplitudeWbosonneutralinoZ5; ParticleChargino2.Array_Decays[28][3] = 2; ParticleChargino2.Array_Comments[28] = "# ~chi_2+ -> W+ ~chi_50";
 
   ParticleChargino2.Array_Decays[29][0] = PDGchargino1; ParticleChargino2.Array_Decays[29][1] = PDGH3; ParticleChargino2.Array_Decays[29][4] = 0; ParticleChargino2.Array_Decays[29][2] = chargino2amplitudechargino1H3; ParticleChargino2.Array_Decays[29][3] = 2; ParticleChargino2.Array_Comments[29] = "# ~chi_2+ -> ~chi_1+ H3";
   ParticleChargino2.Array_Decays[30][0] = PDGchargino1; ParticleChargino2.Array_Decays[30][1] = PDGA2; ParticleChargino2.Array_Decays[30][4] = 0; ParticleChargino2.Array_Decays[30][2] = chargino2amplitudechargino1A2; ParticleChargino2.Array_Decays[30][3] = 2; ParticleChargino2.Array_Comments[30] = "# ~chi_2+ -> ~chi_1+ A2";
   
   ParticleChargino2.Array_Decays[31][0] = PDGneutralino1; ParticleChargino2.Array_Decays[31][1] = PDGup; ParticleChargino2.Array_Decays[31][4] = -PDGdown; ParticleChargino2.Array_Decays[31][2] = chargino2amplitudeneut1udbar; ParticleChargino2.Array_Decays[31][3] = 3; ParticleChargino2.Array_Comments[31] = "# ~chi_2+ -> ~chi_10 u dbar";
   ParticleChargino2.Array_Decays[32][0] = PDGneutralino1; ParticleChargino2.Array_Decays[32][1] = PDGcharm; ParticleChargino2.Array_Decays[32][4] = -PDGstrange; ParticleChargino2.Array_Decays[32][2] = chargino2amplitudeneut1csbar; ParticleChargino2.Array_Decays[32][3] = 3; ParticleChargino2.Array_Comments[32] = "# ~chi_2+ -> ~chi_10 c sbar";
   ParticleChargino2.Array_Decays[33][0] = PDGneutralino1; ParticleChargino2.Array_Decays[33][1] = PDGnuelectron; ParticleChargino2.Array_Decays[33][4] = -PDGelectron; ParticleChargino2.Array_Decays[33][2] = chargino2amplitudeneut1nueebar; ParticleChargino2.Array_Decays[33][3] = 3; ParticleChargino2.Array_Comments[33] = "# ~chi_2+ -> ~chi_10 nu_e e+";
   ParticleChargino2.Array_Decays[34][0] = PDGneutralino1; ParticleChargino2.Array_Decays[34][1] = PDGnumuon; ParticleChargino2.Array_Decays[34][4] = -PDGmuon; ParticleChargino2.Array_Decays[34][2] = chargino2amplitudeneut1numumubar; ParticleChargino2.Array_Decays[34][3] = 3; ParticleChargino2.Array_Comments[34] = "# ~chi_2+ -> ~chi_10 nu_mu mu+";
   ParticleChargino2.Array_Decays[35][0] = PDGneutralino1; ParticleChargino2.Array_Decays[35][1] = PDGnutau; ParticleChargino2.Array_Decays[35][4] = -PDGtau; ParticleChargino2.Array_Decays[35][2] = chargino2amplitudeneut1nutautaubar; ParticleChargino2.Array_Decays[35][3] = 3; ParticleChargino2.Array_Comments[35] = "# ~chi_2+ -> ~chi_10 nu_tau tau+";
   ParticleChargino2.Array_Decays[36][0] = PDGneutralino2; ParticleChargino2.Array_Decays[36][1] = PDGup; ParticleChargino2.Array_Decays[36][4] = -PDGdown; ParticleChargino2.Array_Decays[36][2] = chargino2amplitudeneut2udbar; ParticleChargino2.Array_Decays[36][3] = 3; ParticleChargino2.Array_Comments[36] = "# ~chi_2+ -> ~chi_20 u dbar";
   ParticleChargino2.Array_Decays[37][0] = PDGneutralino2; ParticleChargino2.Array_Decays[37][1] = PDGcharm; ParticleChargino2.Array_Decays[37][4] = -PDGstrange; ParticleChargino2.Array_Decays[37][2] = chargino2amplitudeneut2csbar; ParticleChargino2.Array_Decays[37][3] = 3; ParticleChargino2.Array_Comments[37] = "# ~chi_2+ -> ~chi_20 c sbar";
   ParticleChargino2.Array_Decays[38][0] = PDGneutralino2; ParticleChargino2.Array_Decays[38][1] = PDGnuelectron; ParticleChargino2.Array_Decays[38][4] = -PDGelectron; ParticleChargino2.Array_Decays[38][2] = chargino2amplitudeneut2nueebar; ParticleChargino2.Array_Decays[38][3] = 3; ParticleChargino2.Array_Comments[38] = "# ~chi_2+ -> ~chi_20 nu_e e+";
   ParticleChargino2.Array_Decays[39][0] = PDGneutralino2; ParticleChargino2.Array_Decays[39][1] = PDGnumuon; ParticleChargino2.Array_Decays[39][4] = -PDGmuon; ParticleChargino2.Array_Decays[39][2] = chargino2amplitudeneut2numumubar; ParticleChargino2.Array_Decays[39][3] = 3; ParticleChargino2.Array_Comments[39] = "# ~chi_2+ -> ~chi_20 nu_mu mu+";
   ParticleChargino2.Array_Decays[40][0] = PDGneutralino2; ParticleChargino2.Array_Decays[40][1] = PDGnutau; ParticleChargino2.Array_Decays[40][4] = -PDGtau; ParticleChargino2.Array_Decays[40][2] = chargino2amplitudeneut2nutautaubar; ParticleChargino2.Array_Decays[40][3] = 3; ParticleChargino2.Array_Comments[40] = "# ~chi_2+ -> ~chi_20 nu_tau tau+";
   ParticleChargino2.Array_Decays[41][0] = PDGneutralino3; ParticleChargino2.Array_Decays[41][1] = PDGup; ParticleChargino2.Array_Decays[41][4] = -PDGdown; ParticleChargino2.Array_Decays[41][2] = chargino2amplitudeneut3udbar; ParticleChargino2.Array_Decays[41][3] = 3; ParticleChargino2.Array_Comments[41] = "# ~chi_2+ -> ~chi_30 u dbar";
   ParticleChargino2.Array_Decays[42][0] = PDGneutralino3; ParticleChargino2.Array_Decays[42][1] = PDGcharm; ParticleChargino2.Array_Decays[42][4] = -PDGstrange; ParticleChargino2.Array_Decays[42][2] = chargino2amplitudeneut3csbar; ParticleChargino2.Array_Decays[42][3] = 3; ParticleChargino2.Array_Comments[42] = "# ~chi_2+ -> ~chi_30 c sbar";
   ParticleChargino2.Array_Decays[43][0] = PDGneutralino3; ParticleChargino2.Array_Decays[43][1] = PDGnuelectron; ParticleChargino2.Array_Decays[43][4] = -PDGelectron; ParticleChargino2.Array_Decays[43][2] = chargino2amplitudeneut3nueebar; ParticleChargino2.Array_Decays[43][3] = 3; ParticleChargino2.Array_Comments[43] = "# ~chi_2+ -> ~chi_30 nu_e e+";
   ParticleChargino2.Array_Decays[44][0] = PDGneutralino3; ParticleChargino2.Array_Decays[44][1] = PDGnumuon; ParticleChargino2.Array_Decays[44][4] = -PDGmuon; ParticleChargino2.Array_Decays[44][2] = chargino2amplitudeneut3numumubar; ParticleChargino2.Array_Decays[44][3] = 3; ParticleChargino2.Array_Comments[44] = "# ~chi_2+ -> ~chi_30 nu_mu mu+";
   ParticleChargino2.Array_Decays[45][0] = PDGneutralino3; ParticleChargino2.Array_Decays[45][1] = PDGnutau; ParticleChargino2.Array_Decays[45][4] = -PDGtau; ParticleChargino2.Array_Decays[45][2] = chargino2amplitudeneut3nutautaubar; ParticleChargino2.Array_Decays[45][3] = 3; ParticleChargino2.Array_Comments[45] = "# ~chi_2+ -> ~chi_30 nu_tau tau+";
   ParticleChargino2.Array_Decays[46][0] = PDGneutralino4; ParticleChargino2.Array_Decays[46][1] = PDGup; ParticleChargino2.Array_Decays[46][4] = -PDGdown; ParticleChargino2.Array_Decays[46][2] = chargino2amplitudeneut4udbar; ParticleChargino2.Array_Decays[46][3] = 3; ParticleChargino2.Array_Comments[46] = "# ~chi_2+ -> ~chi_40 u dbar";
   ParticleChargino2.Array_Decays[47][0] = PDGneutralino4; ParticleChargino2.Array_Decays[47][1] = PDGcharm; ParticleChargino2.Array_Decays[47][4] = -PDGstrange; ParticleChargino2.Array_Decays[47][2] = chargino2amplitudeneut4csbar; ParticleChargino2.Array_Decays[47][3] = 3; ParticleChargino2.Array_Comments[47] = "# ~chi_2+ -> ~chi_40 c sbar";
   ParticleChargino2.Array_Decays[48][0] = PDGneutralino4; ParticleChargino2.Array_Decays[48][1] = PDGnuelectron; ParticleChargino2.Array_Decays[48][4] = -PDGelectron; ParticleChargino2.Array_Decays[48][2] = chargino2amplitudeneut4nueebar; ParticleChargino2.Array_Decays[48][3] = 3; ParticleChargino2.Array_Comments[48] = "# ~chi_2+ -> ~chi_40 nu_e e+";
   ParticleChargino2.Array_Decays[49][0] = PDGneutralino4; ParticleChargino2.Array_Decays[49][1] = PDGnumuon; ParticleChargino2.Array_Decays[49][4] = -PDGmuon; ParticleChargino2.Array_Decays[49][2] = chargino2amplitudeneut4numumubar; ParticleChargino2.Array_Decays[49][3] = 3; ParticleChargino2.Array_Comments[49] = "# ~chi_2+ -> ~chi_40 nu_mu mu+";
   ParticleChargino2.Array_Decays[50][0] = PDGneutralino4; ParticleChargino2.Array_Decays[50][1] = PDGnutau; ParticleChargino2.Array_Decays[50][4] = -PDGtau; ParticleChargino2.Array_Decays[50][2] = chargino2amplitudeneut4nutautaubar; ParticleChargino2.Array_Decays[50][3] = 3; ParticleChargino2.Array_Comments[50] = "# ~chi_2+ -> ~chi_40 nu_tau tau+";

   for(int i = 0; i<ParticleChargino2.No_of_Decays; i++) {
     if (ParticleChargino2.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleChargino2.Array_Comments[i] << " is negative = " << ParticleChargino2.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleChargino2.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }

   double Chargino2_No_1to2_Decays = 0;
   Chargino2_No_1to2_Decays = ParticleChargino2.No_1to2_Decays + ParticleChargino2.No_grav_Decays + ParticleChargino2.No_NMSSM_Decays;
 
   for (int j = 0; j<Chargino2_No_1to2_Decays; j++) {
     ParticleChargino2.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
 
   for (int j=0; j<Chargino2_No_1to2_Decays; j++) {
     ParticleChargino2.two_width = ParticleChargino2.two_width + ParticleChargino2.Array_Decays[j][2];
   }
   for (int j=Chargino2_No_1to2_Decays; j<ParticleChargino2.No_of_Decays; j++) {
     ParticleChargino2.three_width = ParticleChargino2.three_width + ParticleChargino2.Array_Decays[j][2];
   }
 
   if ( ParticleChargino2.three_width != ParticleChargino2.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for chargino2 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleChargino2.No_of_Decays = Chargino2_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleChargino2.total_width = ParticleChargino2.two_width;
     }
   else {
     ParticleChargino2.total_width = ParticleChargino2.two_width + ParticleChargino2.three_width;
   }
   
   if ( ParticleChargino2.total_width != ParticleChargino2.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleChargino2.No_of_Decays; i++) {
       //   fout << i << " " << ParticleChargino2.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in chargino2 total width \n");
     }
 
 }

 ///Neutralino Decays
 ///Neutralino1 Decays
  
 double neutralino1amplitudeuLubar=0, neutralino1amplitudeuRubar=0, neutralino1amplitudeuLbaru=0, neutralino1amplitudeuRbaru=0, neutralino1amplitudedLdbar=0, neutralino1amplitudedRdbar=0, neutralino1amplitudedLbard=0, neutralino1amplitudedRbard=0, neutralino1amplitudecLcbar=0, neutralino1amplitudecRcbar=0, neutralino1amplitudecLbarc=0, neutralino1amplitudecRbarc=0, neutralino1amplitudesLsbar=0, neutralino1amplitudesRsbar=0, neutralino1amplitudesLbars=0, neutralino1amplitudesRbars=0, neutralino1amplitudeeLebar=0, neutralino1amplitudeeRebar=0, neutralino1amplitudeeLbare=0, neutralino1amplitudeeRbare=0, neutralino1amplitudemuLmubar=0, neutralino1amplitudemuRmubar=0, neutralino1amplitudemuLbarmu=0, neutralino1amplitudemuRbarmu=0, neutralino1amplitudesnuenuebar=0, neutralino1amplitudesnuebarnue=0, neutralino1amplitudesnumunumubar=0, neutralino1amplitudesnumubarnumu=0, neutralino1amplitudetopstop1bar=0, neutralino1amplitudetopstop2bar=0, neutralino1amplitudetopbarstop1=0, neutralino1amplitudetopbarstop2=0, neutralino1amplitudebottomsbottom1bar=0, neutralino1amplitudebottomsbottom2bar=0, neutralino1amplitudebottombarsbottom1=0, neutralino1amplitudebottombarsbottom2=0, neutralino1amplitudetaustau1bar=0, neutralino1amplitudetaustau2bar=0, neutralino1amplitudetaubarstau1=0, neutralino1amplitudetaubarstau2=0, neutralino1amplitudenutausnutaubar=0, neutralino1amplitudenutaubarsnutau=0, neutralino1amplitudeWbosonpluscharginoW1=0, neutralino1amplitudeWbosonpluscharginoW2=0, neutralino1amplitudeWbosonminuscharginoW1=0, neutralino1amplitudeWbosonminuscharginoW2=0, neutralino1amplitudeHpluscharginoW1=0, neutralino1amplitudeHpluscharginoW2=0, neutralino1amplitudeHminuscharginoW1=0, neutralino1amplitudeHminuscharginoW2=0, neutralino1amplitudeZbosonneutralino2=0, neutralino1amplitudeZbosonneutralino3=0, neutralino1amplitudeZbosonneutralino4=0, neutralino1amplitudehneutralino2=0, neutralino1amplitudehneutralino3=0, neutralino1amplitudehneutralino4=0, neutralino1amplitudeHneutralino2=0, neutralino1amplitudeHneutralino3=0, neutralino1amplitudeHneutralino4=0, neutralino1amplitudeAneutralino2=0, neutralino1amplitudeAneutralino3=0, neutralino1amplitudeAneutralino4=0, neutralino1amplitudephotongravitino=0, neutralino1amplitudeZgravitino=0, neutralino1amplitudehgravitino=0, neutralino1amplitudeHgravitino=0, neutralino1amplitudeAgravitino=0;

 double neutralino1amplitudechargino1udbar=0, neutralino1amplitudechargino1csbar=0, neutralino1amplitudechargino1enuebar=0, neutralino1amplitudechargino1munumubar=0, neutralino1amplitudechargino1taunutaubar=0, neutralino1amplitudechargino2udbar=0, neutralino1amplitudechargino2csbar=0, neutralino1amplitudechargino2enuebar=0, neutralino1amplitudechargino2munumubar=0, neutralino1amplitudechargino2taunutaubar=0;

 if (flagneut1 == 1) {
   if (nmssmIsIt == false) {
     neutralino1amplitudeuLubar = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mup, mu(1,1), g, gp, mixNeut, 1, 'L', 1);
     neutralino1amplitudeuRubar = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mup, mu(2,1), g, gp, mixNeut, 1, 'R', 1);
     neutralino1amplitudeuLbaru = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mup, mu(1,1), g, gp, mixNeut, 1, 'L', 1); 
     neutralino1amplitudeuRbaru = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mup, mu(2,1), g, gp, mixNeut, 1, 'R', 1);
     neutralino1amplitudedLdbar = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mdo, md(1,1), g, gp, mixNeut, 2, 'L', 1);
     neutralino1amplitudedRdbar = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mdo, md(2,1), g, gp, mixNeut, 2, 'R', 1);
     neutralino1amplitudedLbard = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mdo, md(1,1), g, gp, mixNeut, 2, 'L', 1);
     neutralino1amplitudedRbard = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mdo, md(2,1), g, gp, mixNeut, 2, 'R', 1);
     neutralino1amplitudecLcbar = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mc, mu(1,2), g, gp, mixNeut, 1, 'L', 1);
     neutralino1amplitudecRcbar = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mc, mu(2,2), g, gp, mixNeut, 1, 'R', 1);
     neutralino1amplitudecLbarc = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mc, mu(1,2), g, gp, mixNeut, 1, 'L', 1);
     neutralino1amplitudecRbarc = neutralinoamplitudedecayquarksquarkLorR (mneut(1), mc, mu(2,2), g, gp, mixNeut, 1, 'R', 1);
     neutralino1amplitudesLsbar = neutralinoamplitudedecayquarksquarkLorR (mneut(1), ms, md(1,2), g, gp, mixNeut, 2, 'L', 1);
     neutralino1amplitudesRsbar = neutralinoamplitudedecayquarksquarkLorR (mneut(1), ms, md(2,2), g, gp, mixNeut, 2, 'R', 1);
     neutralino1amplitudesLbars = neutralinoamplitudedecayquarksquarkLorR (mneut(1), ms, md(1,2), g, gp, mixNeut, 2, 'L', 1);
     neutralino1amplitudesRbars = neutralinoamplitudedecayquarksquarkLorR (mneut(1), ms, md(2,2), g, gp, mixNeut, 2, 'R', 1);
     neutralino1amplitudeeLebar = neutralinoamplitudedecayleptonsleptonLorR (mneut(1), mel, me(1,1), g, gp, mixNeut, 'L', 1);
     neutralino1amplitudeeRebar = neutralinoamplitudedecayleptonsleptonLorR (mneut(1), mel, me(2,1), g, gp, mixNeut, 'R', 1);
     neutralino1amplitudeeLbare = neutralinoamplitudedecayleptonsleptonLorR (mneut(1), mel, me(1,1), g, gp, mixNeut, 'L', 1);
     neutralino1amplitudeeRbare = neutralinoamplitudedecayleptonsleptonLorR (mneut(1), mel, me(2,1), g, gp, mixNeut, 'R', 1);
     neutralino1amplitudemuLmubar = neutralinoamplitudedecayleptonsleptonLorR (mneut(1), mmu, me(1,2), g, gp, mixNeut, 'L', 1);
     neutralino1amplitudemuRmubar = neutralinoamplitudedecayleptonsleptonLorR (mneut(1), mmu, me(2,2), g, gp, mixNeut, 'R', 1);
     neutralino1amplitudemuLbarmu = neutralinoamplitudedecayleptonsleptonLorR (mneut(1), mmu, me(1,2), g, gp, mixNeut, 'L', 1);
     neutralino1amplitudemuRbarmu = neutralinoamplitudedecayleptonsleptonLorR (mneut(1), mmu, me(2,2), g, gp, mixNeut, 'R', 1);
     neutralino1amplitudesnuenuebar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(1), 0, msnu(1), g, gp, mixNeut, 1);
     neutralino1amplitudesnuebarnue = neutralinoamplitudedecayneutrinosneutrinoL (mneut(1), 0, msnu(1), g, gp, mixNeut, 1);
     neutralino1amplitudesnumunumubar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(1), 0, msnu(2), g, gp, mixNeut, 1);
     neutralino1amplitudesnumubarnumu = neutralinoamplitudedecayneutrinosneutrinoL (mneut(1), 0, msnu(2), g, gp, mixNeut, 1);
     neutralino1amplitudetopstop1bar = neutralinoamplitudedecaysquark3quarkmix (mneut(1), mt, mu(1,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 1);
     neutralino1amplitudetopstop2bar = neutralinoamplitudedecaysquark3quarkmix (mneut(1), mt, mu(2,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 1);
     neutralino1amplitudetopbarstop1 = neutralinoamplitudedecaysquark3quarkmix (mneut(1), mt, mu(1,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 1);
     neutralino1amplitudetopbarstop2 = neutralinoamplitudedecaysquark3quarkmix (mneut(1), mt, mu(2,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 1);
     neutralino1amplitudebottomsbottom1bar = neutralinoamplitudedecaysquark3quarkmix (mneut(1), mb, md(1,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 1);
     neutralino1amplitudebottomsbottom2bar = neutralinoamplitudedecaysquark3quarkmix (mneut(1), mb, md(2,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 1);
     neutralino1amplitudebottombarsbottom1 = neutralinoamplitudedecaysquark3quarkmix (mneut(1), mb, md(1,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 1);
     neutralino1amplitudebottombarsbottom2 = neutralinoamplitudedecaysquark3quarkmix (mneut(1), mb, md(2,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 1);
     neutralino1amplitudetaustau1bar = neutralinoamplitudedecaystautau (mneut(1), mtau, me(1,3), runmw, thetatau, beta, mixNeut, g, gp, 1, 1);
     neutralino1amplitudetaustau2bar = neutralinoamplitudedecaystautau (mneut(1), mtau, me(2,3), runmw, thetatau, beta, mixNeut, g, gp, 2, 1);
     neutralino1amplitudetaubarstau1 = neutralinoamplitudedecaystautau (mneut(1), mtau, me(1,3), runmw, thetatau, beta, mixNeut, g, gp, 1, 1);
     neutralino1amplitudetaubarstau2 = neutralinoamplitudedecaystautau (mneut(1), mtau, me(2,3), runmw, thetatau, beta, mixNeut, g, gp, 2, 1);
     neutralino1amplitudenutausnutaubar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(1), 0, msnu(3), g, gp, mixNeut, 1);
     neutralino1amplitudenutaubarsnutau = neutralinoamplitudedecayneutrinosneutrinoL (mneut(1), 0, msnu(3), g, gp, mixNeut, 1);
     neutralino1amplitudeWbosonpluscharginoW1 = neutralinoamplitudedecaycharginoWboson (mneut(1), polemw, MCH1, g, thetaL2, thetaR2, mixNeut, 1, 1);
     neutralino1amplitudeWbosonpluscharginoW2 = neutralinoamplitudedecaycharginoWboson (mneut(1), polemw, MCH2, g, thetaL2, thetaR2, mixNeut, 1, 2);
     neutralino1amplitudeWbosonminuscharginoW1 = neutralino1amplitudeWbosonpluscharginoW1;
     neutralino1amplitudeWbosonminuscharginoW2 = neutralino1amplitudeWbosonpluscharginoW2;
     neutralino1amplitudeHpluscharginoW1 = neutralinoamplitudedecaycharginoHplus (mneut(1), mHpm, MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 1, 1);
     neutralino1amplitudeHpluscharginoW2 = neutralinoamplitudedecaycharginoHplus (mneut(1), mHpm, MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 1, 2);
     neutralino1amplitudeHminuscharginoW1 = neutralino1amplitudeHpluscharginoW1;
     neutralino1amplitudeHminuscharginoW2 = neutralino1amplitudeHpluscharginoW2;
     
     neutralino1amplitudechargino1udbar = neutralinoamplitudedecaycharginoffprimebar (mneut(1), mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mch(1), mup, mdo, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 1, 1, onetothree, 'q', 'n');
     neutralino1amplitudechargino1csbar = neutralinoamplitudedecaycharginoffprimebar (mneut(1), mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mch(1), mc, ms, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 1, 1, onetothree, 'q', 'n');
     neutralino1amplitudechargino1enuebar = neutralinoamplitudedecaycharginoffprimebar (mneut(1), me(1,1), me(2,1), msnu(1), 100000000000, polemw, mHpm, mch(1), 0, mel, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 1, 1, onetothree, 'l', 'n');
     neutralino1amplitudechargino1munumubar = neutralinoamplitudedecaycharginoffprimebar (mneut(1), me(1,2), me(2,2), msnu(2), 100000000000, polemw, mHpm, mch(1), 0, mmu, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 1, 1, onetothree, 'l', 'n');
     neutralino1amplitudechargino1taunutaubar = neutralinoamplitudedecaycharginoffprimebar (mneut(1), me(1,3), me(2,3), msnu(3), 100000000000, polemw, mHpm, mch(1), 0, mtau, 0, thetatau-PI/2, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 1, 1, onetothree, 'l', 'n');
     neutralino1amplitudechargino2udbar = neutralinoamplitudedecaycharginoffprimebar (mneut(1), mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mch(2), mup, mdo, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 1, 2, onetothree, 'q', 'n');
     neutralino1amplitudechargino2csbar = neutralinoamplitudedecaycharginoffprimebar (mneut(1), mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mch(2), mc, ms, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 1, 2, onetothree, 'q', 'n');
     neutralino1amplitudechargino2enuebar = neutralinoamplitudedecaycharginoffprimebar (mneut(1), me(1,1), me(2,1), msnu(1), 100000000000, polemw, mHpm, mch(2), 0, mel, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 1, 2, onetothree, 'l', 'n');
     neutralino1amplitudechargino2munumubar = neutralinoamplitudedecaycharginoffprimebar (mneut(1), me(1,2), me(2,2), msnu(2), 100000000000, polemw, mHpm, mch(2), 0, mmu, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 1, 2, onetothree, 'l', 'n');
     neutralino1amplitudechargino2taunutaubar = neutralinoamplitudedecaycharginoffprimebar (mneut(1), me(1,3), me(2,3), msnu(3), 100000000000, polemw, mHpm, mch(2), 0, mtau, 0, thetatau-PI/2, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 1, 2, onetothree, 'l', 'n');
     
     neutralino1amplitudephotongravitino = neutralinoamplitudedecayphotongravitino(mneut(1), mgravitino, MPlreduced, mixNeut, g, gp, 1, gravonoff, neutNLSP);
     neutralino1amplitudeZgravitino = neutralinoamplitudedecayZgravitino(mneut(1), polemz, mgravitino, MPlreduced, mixNeut, g, gp, beta, 1, gravonoff, neutNLSP);
     neutralino1amplitudehgravitino = neutralinoamplitudedecayphigravitino(mneut(1), mh0(1), mgravitino, MPlreduced, mixNeut, alpha, beta, 1, gravonoff, 'h', neutNLSP);
     neutralino1amplitudeHgravitino = neutralinoamplitudedecayphigravitino(mneut(1), mh0(2), mgravitino, MPlreduced, mixNeut, alpha, beta, 1, gravonoff, 'H', neutNLSP);
     neutralino1amplitudeAgravitino = neutralinoamplitudedecayphigravitino(mneut(1), mA0(1), mgravitino, MPlreduced, mixNeut, alpha, beta, 1, gravonoff, 'A', neutNLSP);
   }
   else if (nmssmIsIt == true){
       neutralino1amplitudeWbosonpluscharginoW1 = neutralinoamplitudedecaycharginoWNMSSM (mneut(1), MCH1, polemw, g, thetaL2, thetaR2, mixNeut, 1, 1);
       neutralino1amplitudeWbosonpluscharginoW2 = neutralinoamplitudedecaycharginoWNMSSM (mneut(1), MCH2, polemw, g, thetaL2, thetaR2, mixNeut, 1, 2);
       neutralino1amplitudeWbosonminuscharginoW1 = neutralino1amplitudeWbosonpluscharginoW1;
       neutralino1amplitudeWbosonminuscharginoW2 = neutralino1amplitudeWbosonpluscharginoW2;
       
       neutralino1amplitudeHpluscharginoW1 = neutralinoamplitudecharginoHpmNMSSM (mneut(1), MCH1, mHpm, g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 1, 1);
       neutralino1amplitudeHpluscharginoW2 = neutralinoamplitudecharginoHpmNMSSM (mneut(1), MCH2, mHpm, g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 1, 2);
       neutralino1amplitudeHminuscharginoW1 = neutralino1amplitudeHpluscharginoW1;
       neutralino1amplitudeHminuscharginoW2 = neutralino1amplitudeHpluscharginoW2;
       
       neutralino1amplitudeuLubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), mu(1,1), mup, g, gp, mixNeut, 1, 'u', 'L');
       neutralino1amplitudeuRubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), mu(2,1), mup, g, gp, mixNeut, 1, 'u', 'R');
       neutralino1amplitudeuLbaru = neutralino1amplitudeuLubar;
       neutralino1amplitudeuRbaru = neutralino1amplitudeuRubar;
       neutralino1amplitudedLdbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), md(1,1), mdo, g, gp, mixNeut, 1, 'd', 'L');
       neutralino1amplitudedRdbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), md(2,1), mdo, g, gp, mixNeut, 1, 'd', 'R');
       neutralino1amplitudedLbard = neutralino1amplitudedLdbar;
       neutralino1amplitudedRbard = neutralino1amplitudedRdbar;
       neutralino1amplitudecLcbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), mu(1,2), mc, g, gp, mixNeut, 1, 'u', 'L');
       neutralino1amplitudecRcbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), mu(2,2), mc, g, gp, mixNeut, 1, 'u', 'R');
       neutralino1amplitudecLbarc = neutralino1amplitudecLcbar;
       neutralino1amplitudecRbarc = neutralino1amplitudecRcbar;
       neutralino1amplitudesLsbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), md(1,2), ms, g, gp, mixNeut, 1, 'd', 'L');
       neutralino1amplitudesRsbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), md(2,2), ms, g, gp, mixNeut, 1, 'd', 'R');
       neutralino1amplitudesLbars = neutralino1amplitudesLsbar;
       neutralino1amplitudesRbars = neutralino1amplitudesRsbar;
       neutralino1amplitudeeLebar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), me(1,1), mel, g, gp, mixNeut, 1, 'l', 'L');
       neutralino1amplitudeeRebar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), me(2,1), mel, g, gp, mixNeut, 1, 'l', 'R');
       neutralino1amplitudeeLbare = neutralino1amplitudeeLebar;
       neutralino1amplitudeeRbare = neutralino1amplitudeeRebar;
       neutralino1amplitudemuLmubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), me(1,2), mmu, g, gp, mixNeut, 1, 'l', 'L');
       neutralino1amplitudemuRmubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(1), me(2,2), mmu, g, gp, mixNeut, 1, 'l', 'R');
       neutralino1amplitudemuLbarmu = neutralino1amplitudemuLmubar;
       neutralino1amplitudemuRbarmu = neutralino1amplitudemuRmubar;
       
       neutralino1amplitudetopstop1bar = neutralinoamplitudestoptopNMSSM (mneut(1), mu(1,3), mt, g, gp, thetat, beta, runmw, mixNeut, runmt, 1, 1);
       neutralino1amplitudetopstop2bar = neutralinoamplitudestoptopNMSSM (mneut(1), mu(2,3), mt, g, gp, thetat, beta, runmw, mixNeut, runmt, 1, 2);
       neutralino1amplitudetopbarstop1 = neutralino1amplitudetopstop1bar;
       neutralino1amplitudetopbarstop2 = neutralino1amplitudetopstop2bar;
       neutralino1amplitudebottomsbottom1bar = neutralinoamplitudesbottombottomNMSSM (mneut(1), md(1,3), mb, g, gp, thetab, beta, runmw, mixNeut, runmb, 1, 1);
       neutralino1amplitudebottomsbottom2bar = neutralinoamplitudesbottombottomNMSSM (mneut(1), md(2,3), mb, g, gp, thetab, beta, runmw, mixNeut, runmb, 1, 2);
       neutralino1amplitudebottombarsbottom1 = neutralino1amplitudebottomsbottom1bar;
       neutralino1amplitudebottombarsbottom2 = neutralino1amplitudebottomsbottom2bar;
       
       neutralino1amplitudetaustau1bar = neutralinoamplitudestautauNMSSM (mneut(1), me(1,3), mtau, g, gp, thetatau, beta, runmw, mixNeut, runmtau, 1, 1);
       neutralino1amplitudetaustau2bar = neutralinoamplitudestautauNMSSM (mneut(1), me(2,3), mtau, g, gp, thetatau, beta, runmw, mixNeut, runmtau, 1, 2);
       neutralino1amplitudetaubarstau1 = neutralino1amplitudetaustau1bar;
       neutralino1amplitudetaubarstau2 = neutralino1amplitudetaustau2bar;
       
       neutralino1amplitudesnuenuebar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(1), msnu(1), 0, g, gp, mixNeut, 1);
       neutralino1amplitudesnuebarnue = neutralino1amplitudesnuenuebar;
       neutralino1amplitudesnumunumubar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(1), msnu(2), 0, g, gp, mixNeut, 1);
       neutralino1amplitudesnumubarnumu = neutralino1amplitudesnumunumubar;
       neutralino1amplitudenutausnutaubar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(1), msnu(3), 0, g, gp, mixNeut, 1);
       neutralino1amplitudenutaubarsnutau = neutralino1amplitudenutausnutaubar;
   }
 
   ParticleNeutralino1.Array_Decays[0][0] = -PDGup; ParticleNeutralino1.Array_Decays[0][1] = PDGsupL; ParticleNeutralino1.Array_Decays[0][2] = neutralino1amplitudeuLubar; ParticleNeutralino1.Array_Decays[0][3] = 2; ParticleNeutralino1.Array_Comments[0] = "# ~chi_10 -> ub ~u_L";
   ParticleNeutralino1.Array_Decays[1][0] = -PDGup; ParticleNeutralino1.Array_Decays[1][1] = PDGsupR; ParticleNeutralino1.Array_Decays[1][2] = neutralino1amplitudeuRubar; ParticleNeutralino1.Array_Decays[1][3] = 2; ParticleNeutralino1.Array_Comments[1] = "# ~chi_10 -> ub ~u_R"; 
   ParticleNeutralino1.Array_Decays[2][0] = PDGup; ParticleNeutralino1.Array_Decays[2][1] = -PDGsupL; ParticleNeutralino1.Array_Decays[2][2] = neutralino1amplitudeuLbaru; ParticleNeutralino1.Array_Decays[2][3] = 2; ParticleNeutralino1.Array_Comments[2] = "# ~chi_10 -> u ~u_L*";
   ParticleNeutralino1.Array_Decays[3][0] = PDGup; ParticleNeutralino1.Array_Decays[3][1] = -PDGsupR; ParticleNeutralino1.Array_Decays[3][2] = neutralino1amplitudeuRbaru; ParticleNeutralino1.Array_Decays[3][3] = 2; ParticleNeutralino1.Array_Comments[3] = "# ~chi_10 -> u ~u_R*";
   ParticleNeutralino1.Array_Decays[4][0] = -PDGdown; ParticleNeutralino1.Array_Decays[4][1] = PDGsdownL; ParticleNeutralino1.Array_Decays[4][2] = neutralino1amplitudedLdbar; ParticleNeutralino1.Array_Decays[4][3] = 2; ParticleNeutralino1.Array_Comments[4] = "# ~chi_10 -> db ~d_L";
   ParticleNeutralino1.Array_Decays[5][0] = -PDGdown; ParticleNeutralino1.Array_Decays[5][1] = PDGsdownR; ParticleNeutralino1.Array_Decays[5][2] = neutralino1amplitudedRdbar; ParticleNeutralino1.Array_Decays[5][3] = 2; ParticleNeutralino1.Array_Comments[5] = "# ~chi_10 -> db ~d_R";
   ParticleNeutralino1.Array_Decays[6][0] = PDGdown; ParticleNeutralino1.Array_Decays[6][1] = -PDGsdownL; ParticleNeutralino1.Array_Decays[6][2] = neutralino1amplitudedLbard; ParticleNeutralino1.Array_Decays[6][3] = 2; ParticleNeutralino1.Array_Comments[6] = "# ~chi_10 -> d ~d_L*";
   ParticleNeutralino1.Array_Decays[7][0] = PDGdown; ParticleNeutralino1.Array_Decays[7][1] = -PDGsdownR; ParticleNeutralino1.Array_Decays[7][2] = neutralino1amplitudedRbard; ParticleNeutralino1.Array_Decays[7][3] = 2; ParticleNeutralino1.Array_Comments[7] = "# ~chi_10 -> d ~d_R*";
   ParticleNeutralino1.Array_Decays[8][0] = -PDGcharm; ParticleNeutralino1.Array_Decays[8][1] = PDGscharmL; ParticleNeutralino1.Array_Decays[8][2] = neutralino1amplitudecLcbar; ParticleNeutralino1.Array_Decays[8][3] = 2; ParticleNeutralino1.Array_Comments[8] = "# ~chi_10 -> cb ~c_L";
   ParticleNeutralino1.Array_Decays[9][0] = -PDGcharm; ParticleNeutralino1.Array_Decays[9][1] = PDGscharmR; ParticleNeutralino1.Array_Decays[9][2] = neutralino1amplitudecRcbar; ParticleNeutralino1.Array_Decays[9][3] = 2; ParticleNeutralino1.Array_Comments[9] = "# ~chi_10 -> cb ~c_R";
   ParticleNeutralino1.Array_Decays[10][0] = PDGcharm; ParticleNeutralino1.Array_Decays[10][1] = -PDGscharmL; ParticleNeutralino1.Array_Decays[10][2] = neutralino1amplitudecLbarc; ParticleNeutralino1.Array_Decays[10][3] = 2; ParticleNeutralino1.Array_Comments[10] = "# ~chi_10 -> c ~c_L*";
   ParticleNeutralino1.Array_Decays[11][0] = PDGcharm; ParticleNeutralino1.Array_Decays[11][1] = -PDGscharmR; ParticleNeutralino1.Array_Decays[11][2] = neutralino1amplitudecRbarc; ParticleNeutralino1.Array_Decays[11][3] = 2; ParticleNeutralino1.Array_Comments[11] = "# ~chi_10 -> c ~c_R*";
   ParticleNeutralino1.Array_Decays[12][0] = -PDGstrange; ParticleNeutralino1.Array_Decays[12][1] = PDGsstrangeL; ParticleNeutralino1.Array_Decays[12][2] = neutralino1amplitudesLsbar; ParticleNeutralino1.Array_Decays[12][3] = 2; ParticleNeutralino1.Array_Comments[12] = "# ~chi_10 -> sb ~s_L";
   ParticleNeutralino1.Array_Decays[13][0] = -PDGstrange; ParticleNeutralino1.Array_Decays[13][1] = PDGsstrangeR; ParticleNeutralino1.Array_Decays[13][2] = neutralino1amplitudesRsbar; ParticleNeutralino1.Array_Decays[13][3] = 2; ParticleNeutralino1.Array_Comments[13] = "# ~chi_10 -> sb ~s_R";
   ParticleNeutralino1.Array_Decays[14][0] = PDGstrange; ParticleNeutralino1.Array_Decays[14][1] = -PDGsstrangeL; ParticleNeutralino1.Array_Decays[14][2] = neutralino1amplitudesLbars; ParticleNeutralino1.Array_Decays[14][3] = 2; ParticleNeutralino1.Array_Comments[14] = "# ~chi_10 -> s ~s_L*";
   ParticleNeutralino1.Array_Decays[15][0] = PDGstrange; ParticleNeutralino1.Array_Decays[15][1] = -PDGsstrangeR; ParticleNeutralino1.Array_Decays[15][2] = neutralino1amplitudesRbars; ParticleNeutralino1.Array_Decays[15][3] = 2; ParticleNeutralino1.Array_Comments[15] = "# ~chi_10 -> s ~s_R*";
   ParticleNeutralino1.Array_Decays[16][0] = -PDGelectron; ParticleNeutralino1.Array_Decays[16][1] = PDGselectronL; ParticleNeutralino1.Array_Decays[16][2] = neutralino1amplitudeeLebar; ParticleNeutralino1.Array_Decays[16][3] = 2; ParticleNeutralino1.Array_Comments[16] = "# ~chi_10 -> e+ ~e_L-";
   ParticleNeutralino1.Array_Decays[17][0] = -PDGelectron; ParticleNeutralino1.Array_Decays[17][1] = PDGselectronR; ParticleNeutralino1.Array_Decays[17][2] = neutralino1amplitudeeRebar; ParticleNeutralino1.Array_Decays[17][3] = 2; ParticleNeutralino1.Array_Comments[17] = "# ~chi_10 -> e+ ~e_R-";
   ParticleNeutralino1.Array_Decays[18][0] = PDGelectron; ParticleNeutralino1.Array_Decays[18][1] = -PDGselectronL; ParticleNeutralino1.Array_Decays[18][2] = neutralino1amplitudeeLbare; ParticleNeutralino1.Array_Decays[18][3] = 2; ParticleNeutralino1.Array_Comments[18] = "# ~chi_10 -> e- ~e_L+";
   ParticleNeutralino1.Array_Decays[19][0] = PDGelectron; ParticleNeutralino1.Array_Decays[19][1] = -PDGselectronR; ParticleNeutralino1.Array_Decays[19][2] = neutralino1amplitudeeRbare; ParticleNeutralino1.Array_Decays[19][3] = 2; ParticleNeutralino1.Array_Comments[19] = "# ~chi_10 -> e- ~e_R+";
   ParticleNeutralino1.Array_Decays[20][0] = -PDGmuon; ParticleNeutralino1.Array_Decays[20][1] = PDGsmuonL; ParticleNeutralino1.Array_Decays[20][2] = neutralino1amplitudemuLmubar; ParticleNeutralino1.Array_Decays[20][3] = 2; ParticleNeutralino1.Array_Comments[20] = "# ~chi_10 -> mu+ ~mu_L-";
   ParticleNeutralino1.Array_Decays[21][0] = -PDGmuon; ParticleNeutralino1.Array_Decays[21][1] = PDGsmuonR; ParticleNeutralino1.Array_Decays[21][2] = neutralino1amplitudemuRmubar; ParticleNeutralino1.Array_Decays[21][3] = 2; ParticleNeutralino1.Array_Comments[21] = "# ~chi_10 -> mu+ ~mu_R-";
   ParticleNeutralino1.Array_Decays[22][0] = PDGmuon; ParticleNeutralino1.Array_Decays[22][1] = -PDGsmuonL; ParticleNeutralino1.Array_Decays[22][2] = neutralino1amplitudemuLbarmu; ParticleNeutralino1.Array_Decays[22][3] = 2; ParticleNeutralino1.Array_Comments[22] = "# ~chi_10 -> mu- ~mu_L+";
   ParticleNeutralino1.Array_Decays[23][0] = PDGmuon; ParticleNeutralino1.Array_Decays[23][1] = -PDGsmuonR; ParticleNeutralino1.Array_Decays[23][2] = neutralino1amplitudemuRbarmu; ParticleNeutralino1.Array_Decays[23][3] = 2; ParticleNeutralino1.Array_Comments[23] = "# ~chi_10 -> mu- ~mu_R+";
   ParticleNeutralino1.Array_Decays[24][0] = PDGnuelectron; ParticleNeutralino1.Array_Decays[24][1] = -PDGnuselectronL; ParticleNeutralino1.Array_Decays[24][2] = neutralino1amplitudesnuebarnue; ParticleNeutralino1.Array_Decays[24][3] = 2; ParticleNeutralino1.Array_Comments[24] = "# ~chi_10 -> nu_e ~nu_eL*";
   ParticleNeutralino1.Array_Decays[25][0] = -PDGnuelectron; ParticleNeutralino1.Array_Decays[25][1] = PDGnuselectronL; ParticleNeutralino1.Array_Decays[25][2] = neutralino1amplitudesnuenuebar; ParticleNeutralino1.Array_Decays[25][3] = 2; ParticleNeutralino1.Array_Comments[25] = "# ~chi_10 -> nu_eb ~nu_eL";
   ParticleNeutralino1.Array_Decays[26][0] = PDGnumuon; ParticleNeutralino1.Array_Decays[26][1] = -PDGnusmuonL; ParticleNeutralino1.Array_Decays[26][2] = neutralino1amplitudesnumubarnumu; ParticleNeutralino1.Array_Decays[26][3] = 2; ParticleNeutralino1.Array_Comments[26] = "# ~chi_10 -> nu_mu ~nu_muL*";
   ParticleNeutralino1.Array_Decays[27][0] = -PDGnumuon; ParticleNeutralino1.Array_Decays[27][1] = PDGnusmuonL; ParticleNeutralino1.Array_Decays[27][2] = neutralino1amplitudesnumunumubar; ParticleNeutralino1.Array_Decays[27][3] = 2; ParticleNeutralino1.Array_Comments[27] = "# ~chi_10 -> nu_mub ~nu_muL";
   ParticleNeutralino1.Array_Decays[28][0] = PDGtop; ParticleNeutralino1.Array_Decays[28][1] = -PDGstop1; ParticleNeutralino1.Array_Decays[28][2] = neutralino1amplitudetopstop1bar; ParticleNeutralino1.Array_Decays[28][3] = 2; ParticleNeutralino1.Array_Comments[28] = "# ~chi_10 -> t ~t_1*";
   ParticleNeutralino1.Array_Decays[29][0] = PDGtop; ParticleNeutralino1.Array_Decays[29][1] = -PDGstop2; ParticleNeutralino1.Array_Decays[29][2] = neutralino1amplitudetopstop2bar; ParticleNeutralino1.Array_Decays[29][3] = 2; ParticleNeutralino1.Array_Comments[29] = "# ~chi_10 -> t ~t_2*";
   ParticleNeutralino1.Array_Decays[30][0] = -PDGtop; ParticleNeutralino1.Array_Decays[30][1] = PDGstop1; ParticleNeutralino1.Array_Decays[30][2] = neutralino1amplitudetopbarstop1; ParticleNeutralino1.Array_Decays[30][3] = 2; ParticleNeutralino1.Array_Comments[30] = "# ~chi_10 -> tb ~t_1";
   ParticleNeutralino1.Array_Decays[31][0] = -PDGtop; ParticleNeutralino1.Array_Decays[31][1] = PDGstop2; ParticleNeutralino1.Array_Decays[31][2] = neutralino1amplitudetopbarstop2; ParticleNeutralino1.Array_Decays[31][3] = 2; ParticleNeutralino1.Array_Comments[31] = "# ~chi_10 -> tb ~t_2";
   ParticleNeutralino1.Array_Decays[32][0] = PDGbottom; ParticleNeutralino1.Array_Decays[32][1] = -PDGsbottom1; ParticleNeutralino1.Array_Decays[32][2] = neutralino1amplitudebottomsbottom1bar; ParticleNeutralino1.Array_Decays[32][3] = 2; ParticleNeutralino1.Array_Comments[32] = "# ~chi_10 -> b ~b_1*";
   ParticleNeutralino1.Array_Decays[33][0] = PDGbottom; ParticleNeutralino1.Array_Decays[33][1] = -PDGsbottom2; ParticleNeutralino1.Array_Decays[33][2] = neutralino1amplitudebottomsbottom2bar; ParticleNeutralino1.Array_Decays[33][3] = 2; ParticleNeutralino1.Array_Comments[33] = "# ~chi_10 -> b ~b_2*";
   ParticleNeutralino1.Array_Decays[34][0] = -PDGbottom; ParticleNeutralino1.Array_Decays[34][1] = PDGsbottom1; ParticleNeutralino1.Array_Decays[34][2] = neutralino1amplitudebottombarsbottom1; ParticleNeutralino1.Array_Decays[34][3] = 2; ParticleNeutralino1.Array_Comments[34] = "# ~chi_10 -> bb ~b_1";
   ParticleNeutralino1.Array_Decays[35][0] = -PDGbottom; ParticleNeutralino1.Array_Decays[35][1] = PDGsbottom2; ParticleNeutralino1.Array_Decays[35][2] = neutralino1amplitudebottombarsbottom2; ParticleNeutralino1.Array_Decays[35][3] = 2; ParticleNeutralino1.Array_Comments[35] = "# ~chi_10 -> bb ~b_2";
   ParticleNeutralino1.Array_Decays[36][0] = -PDGstau1; ParticleNeutralino1.Array_Decays[36][1] = PDGtau; ParticleNeutralino1.Array_Decays[36][2] = neutralino1amplitudetaustau1bar; ParticleNeutralino1.Array_Decays[36][3] = 2; ParticleNeutralino1.Array_Comments[36] = "# ~chi_10 -> tau- ~tau_1+";
   ParticleNeutralino1.Array_Decays[37][0] = -PDGstau2; ParticleNeutralino1.Array_Decays[37][1] = PDGtau; ParticleNeutralino1.Array_Decays[37][2] = neutralino1amplitudetaustau2bar; ParticleNeutralino1.Array_Decays[37][3] = 2; ParticleNeutralino1.Array_Comments[37] = "# ~chi_10 -> tau- ~tau_2+";
   ParticleNeutralino1.Array_Decays[38][0] = PDGstau1; ParticleNeutralino1.Array_Decays[38][1] = -PDGtau; ParticleNeutralino1.Array_Decays[38][2] = neutralino1amplitudetaubarstau1; ParticleNeutralino1.Array_Decays[38][3] = 2; ParticleNeutralino1.Array_Comments[38] = "# ~chi_10 -> tau+ ~tau_1-";
   ParticleNeutralino1.Array_Decays[39][0] = PDGstau2; ParticleNeutralino1.Array_Decays[39][1] = -PDGtau; ParticleNeutralino1.Array_Decays[39][2] = neutralino1amplitudetaubarstau2; ParticleNeutralino1.Array_Decays[39][3] = 2; ParticleNeutralino1.Array_Comments[39] = "# ~chi_10 -> tau+ ~tau_2-";
   ParticleNeutralino1.Array_Decays[40][0] = PDGnutau; ParticleNeutralino1.Array_Decays[40][1] = -PDGnustauL; ParticleNeutralino1.Array_Decays[40][2] = neutralino1amplitudenutausnutaubar; ParticleNeutralino1.Array_Decays[40][3] = 2; ParticleNeutralino1.Array_Comments[40] = "# ~chi_10 -> nu_tau ~nu_tauL*";
   ParticleNeutralino1.Array_Decays[41][0] = -PDGnutau; ParticleNeutralino1.Array_Decays[41][1] = PDGnustauL; ParticleNeutralino1.Array_Decays[41][2] = neutralino1amplitudenutaubarsnutau; ParticleNeutralino1.Array_Decays[41][3] = 2; ParticleNeutralino1.Array_Comments[41] = "# ~chi_10 -> nu_taub ~nu_tauL";
   ParticleNeutralino1.Array_Decays[42][0] = PDGWplus; ParticleNeutralino1.Array_Decays[42][1] = -PDGchargino1; ParticleNeutralino1.Array_Decays[42][2] = neutralino1amplitudeWbosonpluscharginoW1; ParticleNeutralino1.Array_Decays[42][3] = 2; ParticleNeutralino1.Array_Comments[42] = "# ~chi_10 -> W+ ~chi_1-";
   ParticleNeutralino1.Array_Decays[43][0] = PDGWplus; ParticleNeutralino1.Array_Decays[43][1] = -PDGchargino2; ParticleNeutralino1.Array_Decays[43][2] = neutralino1amplitudeWbosonpluscharginoW2; ParticleNeutralino1.Array_Decays[43][3] = 2; ParticleNeutralino1.Array_Comments[43] = "# ~chi_10 -> W+ ~chi_2-";
   ParticleNeutralino1.Array_Decays[44][0] = -PDGWplus; ParticleNeutralino1.Array_Decays[44][1] = PDGchargino1; ParticleNeutralino1.Array_Decays[44][2] = neutralino1amplitudeWbosonminuscharginoW1; ParticleNeutralino1.Array_Decays[44][3] = 2; ParticleNeutralino1.Array_Comments[44] = "# ~chi_10 -> W- ~chi_1+";
   ParticleNeutralino1.Array_Decays[45][0] = -PDGWplus; ParticleNeutralino1.Array_Decays[45][1] = PDGchargino2; ParticleNeutralino1.Array_Decays[45][2] = neutralino1amplitudeWbosonminuscharginoW2; ParticleNeutralino1.Array_Decays[45][3] = 2; ParticleNeutralino1.Array_Comments[45] = "# ~chi_10 -> W- ~chi_2+";
   ParticleNeutralino1.Array_Decays[46][0] = PDGHplus; ParticleNeutralino1.Array_Decays[46][1] = -PDGchargino1; ParticleNeutralino1.Array_Decays[46][2] = neutralino1amplitudeHpluscharginoW1; ParticleNeutralino1.Array_Decays[46][3] = 2; ParticleNeutralino1.Array_Comments[46] = "# ~chi_10 -> H+ ~chi_1-";
   ParticleNeutralino1.Array_Decays[47][0] = PDGHplus; ParticleNeutralino1.Array_Decays[47][1] = -PDGchargino2; ParticleNeutralino1.Array_Decays[47][2] = neutralino1amplitudeHpluscharginoW2; ParticleNeutralino1.Array_Decays[47][3] = 2; ParticleNeutralino1.Array_Comments[47] = "# ~chi_10 -> H+ ~chi_2-";
   ParticleNeutralino1.Array_Decays[48][0] = -PDGHplus; ParticleNeutralino1.Array_Decays[48][1] = PDGchargino1; ParticleNeutralino1.Array_Decays[48][2] = neutralino1amplitudeHminuscharginoW1; ParticleNeutralino1.Array_Decays[48][3] = 2; ParticleNeutralino1.Array_Comments[48] = "# ~chi_10 -> H- ~chi_1+";
   ParticleNeutralino1.Array_Decays[49][0] = -PDGHplus; ParticleNeutralino1.Array_Decays[49][1] = PDGchargino2; ParticleNeutralino1.Array_Decays[49][2] = neutralino1amplitudeHminuscharginoW2; ParticleNeutralino1.Array_Decays[49][3] = 2; ParticleNeutralino1.Array_Comments[49] = "# ~chi_10 -> H- ~chi_2+";
   ParticleNeutralino1.Array_Decays[50][0] = PDGZboson; ParticleNeutralino1.Array_Decays[50][1] = PDGneutralino2; ParticleNeutralino1.Array_Decays[50][2] = neutralino1amplitudeZbosonneutralino2; ParticleNeutralino1.Array_Decays[50][3] = 2; ParticleNeutralino1.Array_Comments[50] = "# ~chi_10 -> Z ~chi_20";
   ParticleNeutralino1.Array_Decays[51][0] = PDGZboson; ParticleNeutralino1.Array_Decays[51][1] = PDGneutralino3; ParticleNeutralino1.Array_Decays[51][2] = neutralino1amplitudeZbosonneutralino3; ParticleNeutralino1.Array_Decays[51][3] = 2; ParticleNeutralino1.Array_Comments[51] = "# ~chi_10 -> Z ~chi_30";
   ParticleNeutralino1.Array_Decays[52][0] = PDGZboson; ParticleNeutralino1.Array_Decays[52][1] = PDGneutralino4; ParticleNeutralino1.Array_Decays[52][2] = neutralino1amplitudeZbosonneutralino4; ParticleNeutralino1.Array_Decays[52][3] = 2; ParticleNeutralino1.Array_Comments[52] = "# ~chi_10 -> Z ~chi_40";
   ParticleNeutralino1.Array_Decays[53][0] = PDGh0; ParticleNeutralino1.Array_Decays[53][1] = PDGneutralino2; ParticleNeutralino1.Array_Decays[53][2] = neutralino1amplitudehneutralino2; ParticleNeutralino1.Array_Decays[53][3] = 2; ParticleNeutralino1.Array_Comments[53] = "# ~chi_10 -> h ~chi_20";
   ParticleNeutralino1.Array_Decays[54][0] = PDGh0; ParticleNeutralino1.Array_Decays[54][1] = PDGneutralino3; ParticleNeutralino1.Array_Decays[54][2] = neutralino1amplitudehneutralino3;  ParticleNeutralino1.Array_Decays[54][3] = 2; ParticleNeutralino1.Array_Comments[54] = "# ~chi_10 -> h ~chi_30";
   ParticleNeutralino1.Array_Decays[55][0] = PDGh0; ParticleNeutralino1.Array_Decays[55][1] = PDGneutralino4; ParticleNeutralino1.Array_Decays[55][2] = neutralino1amplitudehneutralino4; ParticleNeutralino1.Array_Decays[55][3] = 2; ParticleNeutralino1.Array_Comments[55] = "# ~chi_10 -> h ~chi_40";
   ParticleNeutralino1.Array_Decays[56][0] = PDGH0; ParticleNeutralino1.Array_Decays[56][1] = PDGneutralino2; ParticleNeutralino1.Array_Decays[56][2] = neutralino1amplitudeHneutralino2; ParticleNeutralino1.Array_Decays[56][3] = 2; ParticleNeutralino1.Array_Comments[56] = "# ~chi_10 -> H ~chi_20";
   ParticleNeutralino1.Array_Decays[57][0] = PDGH0; ParticleNeutralino1.Array_Decays[57][1] = PDGneutralino3; ParticleNeutralino1.Array_Decays[57][2] = neutralino1amplitudeHneutralino3; ParticleNeutralino1.Array_Decays[57][3] = 2; ParticleNeutralino1.Array_Comments[57] = "# ~chi_10 -> H ~chi_30";
   ParticleNeutralino1.Array_Decays[58][0] = PDGH0; ParticleNeutralino1.Array_Decays[58][1] = PDGneutralino4; ParticleNeutralino1.Array_Decays[58][2] = neutralino1amplitudeHneutralino4; ParticleNeutralino1.Array_Decays[58][3] = 2; ParticleNeutralino1.Array_Comments[58] = "# ~chi_10 -> H ~chi_40";
   ParticleNeutralino1.Array_Decays[59][0] = PDGA0; ParticleNeutralino1.Array_Decays[59][1] = PDGneutralino2; ParticleNeutralino1.Array_Decays[59][2] = neutralino1amplitudeAneutralino2; ParticleNeutralino1.Array_Decays[59][3] = 2; ParticleNeutralino1.Array_Comments[59] = "# ~chi_10 -> A ~chi_20";
   ParticleNeutralino1.Array_Decays[60][0] = PDGA0; ParticleNeutralino1.Array_Decays[60][1] = PDGneutralino3; ParticleNeutralino1.Array_Decays[60][2] = neutralino1amplitudeAneutralino3; ParticleNeutralino1.Array_Decays[60][3] = 2; ParticleNeutralino1.Array_Comments[60] = "# ~chi_10 -> A ~chi_30";
   ParticleNeutralino1.Array_Decays[61][0] = PDGA0; ParticleNeutralino1.Array_Decays[61][1] = PDGneutralino4; ParticleNeutralino1.Array_Decays[61][2] = neutralino1amplitudeAneutralino4; ParticleNeutralino1.Array_Decays[61][3] = 2; ParticleNeutralino1.Array_Comments[61] = "# ~chi_10 -> A ~chi_40";

   ParticleNeutralino1.Array_Decays[62][0] = PDGphoton; ParticleNeutralino1.Array_Decays[62][1] = PDGgravitino; ParticleNeutralino1.Array_Decays[62][2] = neutralino1amplitudephotongravitino; ParticleNeutralino1.Array_Decays[62][3] = 2; ParticleNeutralino1.Array_Comments[62] = "# ~chi_10 -> gamma ~G";
   ParticleNeutralino1.Array_Decays[63][0] = PDGZboson; ParticleNeutralino1.Array_Decays[63][1] = PDGgravitino; ParticleNeutralino1.Array_Decays[63][2] = neutralino1amplitudeZgravitino; ParticleNeutralino1.Array_Decays[63][3] = 2; ParticleNeutralino1.Array_Comments[63] = "# ~chi_10 -> Z ~G";
   ParticleNeutralino1.Array_Decays[64][0] = PDGh0; ParticleNeutralino1.Array_Decays[64][1] = PDGgravitino; ParticleNeutralino1.Array_Decays[64][2] = neutralino1amplitudehgravitino; ParticleNeutralino1.Array_Decays[64][3] = 2; ParticleNeutralino1.Array_Comments[64] = "# ~chi_10 -> h ~G";
   ParticleNeutralino1.Array_Decays[65][0] = PDGH0; ParticleNeutralino1.Array_Decays[65][1] = PDGgravitino; ParticleNeutralino1.Array_Decays[65][2] = neutralino1amplitudeHgravitino; ParticleNeutralino1.Array_Decays[65][3] = 2; ParticleNeutralino1.Array_Comments[65] = "# ~chi_10 -> H ~G";
   ParticleNeutralino1.Array_Decays[66][0] = PDGA0; ParticleNeutralino1.Array_Decays[66][1] = PDGgravitino; ParticleNeutralino1.Array_Decays[66][2] = neutralino1amplitudeAgravitino; ParticleNeutralino1.Array_Decays[66][3] = 2; ParticleNeutralino1.Array_Comments[66] = "# ~chi_10 -> A ~G";
   
   ParticleNeutralino1.Array_Decays[67][0] = PDGchargino1; ParticleNeutralino1.Array_Decays[67][1] = PDGup; ParticleNeutralino1.Array_Decays[67][4] = -PDGdown; ParticleNeutralino1.Array_Decays[67][2] = neutralino1amplitudechargino1udbar; ParticleNeutralino1.Array_Decays[67][3] = 3; ParticleNeutralino1.Array_Comments[67] = "# ~chi_10 -> chi_1- u db";
   ParticleNeutralino1.Array_Decays[68][0] = PDGchargino1; ParticleNeutralino1.Array_Decays[68][1] = PDGcharm; ParticleNeutralino1.Array_Decays[68][4] = -PDGstrange; ParticleNeutralino1.Array_Decays[68][2] = neutralino1amplitudechargino1csbar; ParticleNeutralino1.Array_Decays[68][3] = 3; ParticleNeutralino1.Array_Comments[68] = "# ~chi_10 -> chi_1- c sb";
   ParticleNeutralino1.Array_Decays[69][0] = PDGchargino1; ParticleNeutralino1.Array_Decays[69][1] = PDGnuelectron; ParticleNeutralino1.Array_Decays[69][4] = -PDGelectron; ParticleNeutralino1.Array_Decays[69][2] = neutralino1amplitudechargino1enuebar; ParticleNeutralino1.Array_Decays[69][3] = 3; ParticleNeutralino1.Array_Comments[69] = "# ~chi_10 -> chi_1- nu_e eb";
   ParticleNeutralino1.Array_Decays[70][0] = PDGchargino1; ParticleNeutralino1.Array_Decays[70][1] = PDGnumuon; ParticleNeutralino1.Array_Decays[70][4] = -PDGmuon; ParticleNeutralino1.Array_Decays[70][2] = neutralino1amplitudechargino1munumubar; ParticleNeutralino1.Array_Decays[70][3] = 3; ParticleNeutralino1.Array_Comments[70] = "# ~chi_10 -> chi_1- nu_mu mub";
   ParticleNeutralino1.Array_Decays[71][0] = PDGchargino1; ParticleNeutralino1.Array_Decays[71][1] = PDGnutau; ParticleNeutralino1.Array_Decays[71][4] = -PDGtau; ParticleNeutralino1.Array_Decays[71][2] = neutralino1amplitudechargino1taunutaubar; ParticleNeutralino1.Array_Decays[71][3] = 3; ParticleNeutralino1.Array_Comments[71] = "# ~chi_10 -> chi_1- nu_tau taub";
   ParticleNeutralino1.Array_Decays[72][0] = PDGchargino2; ParticleNeutralino1.Array_Decays[72][1] = PDGup; ParticleNeutralino1.Array_Decays[72][4] = -PDGdown; ParticleNeutralino1.Array_Decays[72][2] = neutralino1amplitudechargino2udbar; ParticleNeutralino1.Array_Decays[72][3] = 3; ParticleNeutralino1.Array_Comments[72] = "# ~chi_10 -> chi_2- u dbar";
   ParticleNeutralino1.Array_Decays[73][0] = PDGchargino2; ParticleNeutralino1.Array_Decays[73][1] = PDGcharm; ParticleNeutralino1.Array_Decays[73][4] = -PDGstrange; ParticleNeutralino1.Array_Decays[73][2] = neutralino1amplitudechargino2csbar; ParticleNeutralino1.Array_Decays[73][3] = 3; ParticleNeutralino1.Array_Comments[73] = "# ~chi_10 -> chi_2- c sbar";
   ParticleNeutralino1.Array_Decays[74][0] = PDGchargino2; ParticleNeutralino1.Array_Decays[74][1] = PDGnuelectron; ParticleNeutralino1.Array_Decays[74][4] = -PDGelectron; ParticleNeutralino1.Array_Decays[74][2] = neutralino1amplitudechargino2enuebar; ParticleNeutralino1.Array_Decays[74][3] = 3; ParticleNeutralino1.Array_Comments[74] = "# ~chi_10 -> chi_2- nu_e eb";
   ParticleNeutralino1.Array_Decays[75][0] = PDGchargino2; ParticleNeutralino1.Array_Decays[75][1] = PDGnumuon; ParticleNeutralino1.Array_Decays[75][4] = -PDGmuon; ParticleNeutralino1.Array_Decays[75][2] = neutralino1amplitudechargino2munumubar; ParticleNeutralino1.Array_Decays[75][3] = 3; ParticleNeutralino1.Array_Comments[75] = "# ~chi_10 -> chi_2- nu_mu mub";
   ParticleNeutralino1.Array_Decays[76][0] = PDGchargino2; ParticleNeutralino1.Array_Decays[76][1] = PDGnutau; ParticleNeutralino1.Array_Decays[67][4] = -PDGtau; ParticleNeutralino1.Array_Decays[76][2] = neutralino1amplitudechargino2taunutaubar; ParticleNeutralino1.Array_Decays[76][3] = 3; ParticleNeutralino1.Array_Comments[76] = "# ~chi_10 -> chi_2- nu_tau taubar";

   for(int i = 0; i<ParticleNeutralino1.No_of_Decays; i++) {
     if (ParticleNeutralino1.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleNeutralino1.Array_Comments[i] << " is negative = " << ParticleNeutralino1.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleNeutralino1.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }
   
   double Neut1_No_1to2_Decays = 0;
   
   Neut1_No_1to2_Decays = ParticleNeutralino1.No_1to2_Decays + ParticleNeutralino1.No_grav_Decays;
   
   for (int j = 0; j<Neut1_No_1to2_Decays; j++) {
     ParticleNeutralino1.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<Neut1_No_1to2_Decays; j++) {
     ParticleNeutralino1.two_width = ParticleNeutralino1.two_width + ParticleNeutralino1.Array_Decays[j][2];
   }
   for (int j=Neut1_No_1to2_Decays; j<ParticleNeutralino1.No_of_Decays; j++) {
     ParticleNeutralino1.three_width = ParticleNeutralino1.three_width + ParticleNeutralino1.Array_Decays[j][2];
   }
   
   if ( ParticleNeutralino1.three_width != ParticleNeutralino1.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for neutralino 1 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleNeutralino1.No_of_Decays = Neut1_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleNeutralino1.total_width = ParticleNeutralino1.two_width;
     }
   else {
     ParticleNeutralino1.total_width = ParticleNeutralino1.two_width + ParticleNeutralino1.three_width;
   }

   if ( ParticleNeutralino1.total_width != ParticleNeutralino1.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleNeutralino1.No_of_Decays; i++) {
       //   fout << i << " " << ParticleNeutralino1.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in Neutralino1 total width \n");
     }
   
 }
///Neutralino2 Decays
 double neutralino2amplitudeuLubar=0, neutralino2amplitudeuRubar=0, neutralino2amplitudeuLbaru=0, neutralino2amplitudeuRbaru=0, neutralino2amplitudedLdbar=0, neutralino2amplitudedRdbar=0, neutralino2amplitudedLbard=0, neutralino2amplitudedRbard=0, neutralino2amplitudecLcbar=0, neutralino2amplitudecRcbar=0, neutralino2amplitudecLbarc=0, neutralino2amplitudecRbarc=0, neutralino2amplitudesLsbar=0, neutralino2amplitudesRsbar=0, neutralino2amplitudesLbars=0, neutralino2amplitudesRbars=0, neutralino2amplitudeeLebar=0, neutralino2amplitudeeRebar=0, neutralino2amplitudeeLbare=0, neutralino2amplitudeeRbare=0, neutralino2amplitudemuLmubar=0, neutralino2amplitudemuRmubar=0, neutralino2amplitudemuLbarmu=0, neutralino2amplitudemuRbarmu=0, neutralino2amplitudesnuenuebar=0, neutralino2amplitudesnuebarnue=0, neutralino2amplitudesnumunumubar=0, neutralino2amplitudesnumubarnumu=0, neutralino2amplitudetopstop1bar=0, neutralino2amplitudetopstop2bar=0, neutralino2amplitudetopbarstop1=0, neutralino2amplitudetopbarstop2=0, neutralino2amplitudebottomsbottom1bar=0, neutralino2amplitudebottomsbottom2bar=0, neutralino2amplitudebottombarsbottom1=0, neutralino2amplitudebottombarsbottom2=0, neutralino2amplitudetaustau1bar=0, neutralino2amplitudetaustau2bar=0, neutralino2amplitudetaubarstau1=0, neutralino2amplitudetaubarstau2=0, neutralino2amplitudenutausnutaubar=0, neutralino2amplitudenutaubarsnutau=0, neutralino2amplitudeWbosonpluscharginoW1=0, neutralino2amplitudeWbosonpluscharginoW2=0, neutralino2amplitudeWbosonminuscharginoW1=0, neutralino2amplitudeWbosonminuscharginoW2=0, neutralino2amplitudeHpluscharginoW1=0, neutralino2amplitudeHpluscharginoW2=0, neutralino2amplitudeHminuscharginoW1=0, neutralino2amplitudeHminuscharginoW2=0, neutralino2amplitudeZbosonneutralino1=0, neutralino2amplitudeZbosonneutralino3=0, neutralino2amplitudeZbosonneutralino4=0, neutralino2amplitudehneutralino1=0, neutralino2amplitudehneutralino3=0, neutralino2amplitudehneutralino4=0, neutralino2amplitudeHneutralino1=0, neutralino2amplitudeHneutralino3=0, neutralino2amplitudeHneutralino4=0, neutralino2amplitudeAneutralino1=0, neutralino2amplitudeAneutralino3=0, neutralino2amplitudeAneutralino4=0, neutralino2amplitudephotongravitino=0, neutralino2amplitudeZgravitino=0, neutralino2amplitudehgravitino=0, neutralino2amplitudeHgravitino=0, neutralino2amplitudeAgravitino=0;

 double neutralino2amplitudeneut1uubar=0, neutralino2amplitudeneut1ddbar=0, neutralino2amplitudeneut1ccbar=0, neutralino2amplitudeneut1ssbar=0, neutralino2amplitudeneut1ttbar=0, neutralino2amplitudeneut1bbbar=0, neutralino2amplitudeneut1eebar=0, neutralino2amplitudeneut1mumubar=0, neutralino2amplitudeneut1tautaubar=0, neutralino2amplitudeneut1nuenuebar=0, neutralino2amplitudeneut1numunumubar=0, neutralino2amplitudeneut1nutaunutaubar=0, neutralino2amplitudechargino1udbar=0, neutralino2amplitudechargino1csbar=0, neutralino2amplitudechargino1enuebar=0, neutralino2amplitudechargino1munumubar=0, neutralino2amplitudechargino1taunutaubar=0, neutralino2amplitudechargino2udbar=0, neutralino2amplitudechargino2csbar=0, neutralino2amplitudechargino2enuebar=0, neutralino2amplitudechargino2munumubar=0, neutralino2amplitudechargino2taunutaubar=0;

 double neutralino2amplitudeH3neutralino1 = 0, neutralino2amplitudeA2neutralino1 = 0;

 if (flagneut2 == 1) {
   if (nmssmIsIt == false) {
     neutralino2amplitudeuLubar = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mup, mu(1,1), g, gp, mixNeut, 1, 'L', 2);
     neutralino2amplitudeuRubar = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mup, mu(2,1), g, gp, mixNeut, 1, 'R', 2);
     neutralino2amplitudeuLbaru = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mup, mu(1,1), g, gp, mixNeut, 1, 'L', 2); 
     neutralino2amplitudeuRbaru = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mup, mu(2,1), g, gp, mixNeut, 1, 'R', 2);
     neutralino2amplitudedLdbar = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mdo, md(1,1), g, gp, mixNeut, 2, 'L', 2);
     neutralino2amplitudedRdbar = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mdo, md(2,1), g, gp, mixNeut, 2, 'R', 2);
     neutralino2amplitudedLbard = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mdo, md(1,1), g, gp, mixNeut, 2, 'L', 2);
     neutralino2amplitudedRbard = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mdo, md(2,1), g, gp, mixNeut, 2, 'R', 2);
     neutralino2amplitudecLcbar = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mc, mu(1,2), g, gp, mixNeut, 1, 'L', 2);
     neutralino2amplitudecRcbar = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mc, mu(2,2), g, gp, mixNeut, 1, 'R', 2);
     neutralino2amplitudecLbarc = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mc, mu(1,2), g, gp, mixNeut, 1, 'L', 2);
     neutralino2amplitudecRbarc = neutralinoamplitudedecayquarksquarkLorR (mneut(2), mc, mu(2,2), g, gp, mixNeut, 1, 'R', 2);
     neutralino2amplitudesLsbar = neutralinoamplitudedecayquarksquarkLorR (mneut(2), ms, md(1,2), g, gp, mixNeut, 2, 'L', 2);
     neutralino2amplitudesRsbar = neutralinoamplitudedecayquarksquarkLorR (mneut(2), ms, md(2,2), g, gp, mixNeut, 2, 'R', 2);
     neutralino2amplitudesLbars = neutralinoamplitudedecayquarksquarkLorR (mneut(2), ms, md(1,2), g, gp, mixNeut, 2, 'L', 2);
     neutralino2amplitudesRbars = neutralinoamplitudedecayquarksquarkLorR (mneut(2), ms, md(2,2), g, gp, mixNeut, 2, 'R', 2);
     neutralino2amplitudeeLebar = neutralinoamplitudedecayleptonsleptonLorR (mneut(2), mel, me(1,1), g, gp, mixNeut, 'L', 2);
     neutralino2amplitudeeRebar = neutralinoamplitudedecayleptonsleptonLorR (mneut(2), mel, me(2,1), g, gp, mixNeut, 'R', 2);
     neutralino2amplitudeeLbare = neutralinoamplitudedecayleptonsleptonLorR (mneut(2), mel, me(1,1), g, gp, mixNeut, 'L', 2);
     neutralino2amplitudeeRbare = neutralinoamplitudedecayleptonsleptonLorR (mneut(2), mel, me(2,1), g, gp, mixNeut, 'R', 2);
     neutralino2amplitudemuLmubar = neutralinoamplitudedecayleptonsleptonLorR (mneut(2), mmu, me(1,2), g, gp, mixNeut, 'L', 2);
     neutralino2amplitudemuRmubar = neutralinoamplitudedecayleptonsleptonLorR (mneut(2), mmu, me(2,2), g, gp, mixNeut, 'R', 2);
     neutralino2amplitudemuLbarmu = neutralinoamplitudedecayleptonsleptonLorR (mneut(2), mmu, me(1,2), g, gp, mixNeut, 'L', 2);
     neutralino2amplitudemuRbarmu = neutralinoamplitudedecayleptonsleptonLorR (mneut(2), mmu, me(2,2), g, gp, mixNeut, 'R', 2);
     neutralino2amplitudesnuenuebar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(2), 0, msnu(1), g, gp, mixNeut, 2);
     neutralino2amplitudesnuebarnue = neutralinoamplitudedecayneutrinosneutrinoL (mneut(2), 0, msnu(1), g, gp, mixNeut, 2);
     neutralino2amplitudesnumunumubar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(2), 0, msnu(2), g, gp, mixNeut, 2);
     neutralino2amplitudesnumubarnumu = neutralinoamplitudedecayneutrinosneutrinoL (mneut(2), 0, msnu(2), g, gp, mixNeut, 2);
     neutralino2amplitudetopstop1bar = neutralinoamplitudedecaysquark3quarkmix (mneut(2), mt, mu(1,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 2);
     neutralino2amplitudetopstop2bar = neutralinoamplitudedecaysquark3quarkmix (mneut(2), mt, mu(2,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 2);
     neutralino2amplitudetopbarstop1 = neutralinoamplitudedecaysquark3quarkmix (mneut(2), mt, mu(1,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 2);
     neutralino2amplitudetopbarstop2 = neutralinoamplitudedecaysquark3quarkmix (mneut(2), mt, mu(2,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 2);
     neutralino2amplitudebottomsbottom1bar = neutralinoamplitudedecaysquark3quarkmix (mneut(2), mb, md(1,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 2);
     neutralino2amplitudebottomsbottom2bar = neutralinoamplitudedecaysquark3quarkmix (mneut(2), mb, md(2,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 2);
     neutralino2amplitudebottombarsbottom1 = neutralinoamplitudedecaysquark3quarkmix (mneut(2), mb, md(1,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 2);
     neutralino2amplitudebottombarsbottom2 = neutralinoamplitudedecaysquark3quarkmix (mneut(2), mb, md(2,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 2);
     neutralino2amplitudetaustau1bar = neutralinoamplitudedecaystautau (mneut(2), mtau, me(1,3), runmw, thetatau, beta, mixNeut, g, gp, 1, 2);
     neutralino2amplitudetaustau2bar = neutralinoamplitudedecaystautau (mneut(2), mtau, me(2,3), runmw, thetatau, beta, mixNeut, g, gp, 2, 2);
     neutralino2amplitudetaubarstau1 = neutralinoamplitudedecaystautau (mneut(2), mtau, me(1,3), runmw, thetatau, beta, mixNeut, g, gp, 1, 2);
     neutralino2amplitudetaubarstau2 = neutralinoamplitudedecaystautau (mneut(2), mtau, me(2,3), runmw, thetatau, beta, mixNeut, g, gp, 2, 2);
     neutralino2amplitudenutausnutaubar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(2), 0, msnu(3), g, gp, mixNeut, 2);
     neutralino2amplitudenutaubarsnutau = neutralinoamplitudedecayneutrinosneutrinoL (mneut(2), 0, msnu(3), g, gp, mixNeut, 2);
     neutralino2amplitudeWbosonpluscharginoW1 = neutralinoamplitudedecaycharginoWboson (mneut(2), polemw, MCH1, g, thetaL2, thetaR2, mixNeut, 2, 1);
     neutralino2amplitudeWbosonpluscharginoW2 = neutralinoamplitudedecaycharginoWboson (mneut(2), polemw, MCH2, g, thetaL2, thetaR2, mixNeut, 2, 2);
     neutralino2amplitudeWbosonminuscharginoW1 = neutralinoamplitudedecaycharginoWboson (mneut(2), polemw, MCH1, g, thetaL2, thetaR2, mixNeut, 2, 1);
     neutralino2amplitudeWbosonminuscharginoW2 = neutralinoamplitudedecaycharginoWboson (mneut(2), polemw, MCH2, g, thetaL2, thetaR2, mixNeut, 2, 2);
     neutralino2amplitudeHpluscharginoW1 = neutralinoamplitudedecaycharginoHplus (mneut(2), mHpm, MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 2, 1);
     neutralino2amplitudeHpluscharginoW2 = neutralinoamplitudedecaycharginoHplus (mneut(2), mHpm, MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 2, 2);
     neutralino2amplitudeHminuscharginoW1 = neutralinoamplitudedecaycharginoHplus (mneut(2), mHpm, MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 2, 1);
     neutralino2amplitudeHminuscharginoW2 = neutralinoamplitudedecaycharginoHplus (mneut(2), mHpm, MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 2, 2);
     neutralino2amplitudeZbosonneutralino1 = neutralinoamplitudedecayneutralinoZboson (mneut(2), polemz, mneut(1), g, gp, mixNeut, 2, 1);
     
     neutralino2amplitudehneutralino1 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(2), mh0(1), mneut(1), g, gp, mixNeut, alpha, 2, 1, 'h');
     
     neutralino2amplitudeHneutralino1 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(2), mh0(2), mneut(1), g, gp, mixNeut, alpha, 2, 1, 'H');
     neutralino2amplitudeAneutralino1 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(2), mA0(1), mneut(1), g, gp, mixNeut, beta, 2, 1, 'A');
     
     neutralino2amplitudephotongravitino = neutralinoamplitudedecayphotongravitino(mneut(2), mgravitino, MPlreduced, mixNeut, g, gp, 2, gravonoff, neutNLSP);
     neutralino2amplitudeZgravitino = neutralinoamplitudedecayZgravitino(mneut(2), polemz, mgravitino, MPlreduced, mixNeut, g, gp, beta, 2, gravonoff, neutNLSP);
     neutralino2amplitudehgravitino = neutralinoamplitudedecayphigravitino(mneut(2), mh0(1), mgravitino, MPlreduced, mixNeut, alpha, beta, 2, gravonoff, 'h', neutNLSP);
     neutralino2amplitudeHgravitino = neutralinoamplitudedecayphigravitino(mneut(2), mh0(2), mgravitino, MPlreduced, mixNeut, alpha, beta, 2, gravonoff, 'H', neutNLSP);
     neutralino2amplitudeAgravitino = neutralinoamplitudedecayphigravitino(mneut(2), mA0(1), mgravitino, MPlreduced, mixNeut, alpha, beta, 2, gravonoff, 'A', neutNLSP);
     
     neutralino2amplitudeneut1uubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), mu(1,1), mu(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mup, alphas, 0, runmw, g, gp, alpha, beta, runmu, mixNeut, 2, 1, onetothree, 'u');
     neutralino2amplitudeneut1ddbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), md(1,1), md(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mdo, alphas, 0, runmw, g, gp, alpha, beta, runmd, mixNeut, 2, 1, onetothree, 'd');
     neutralino2amplitudeneut1ccbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), mu(1,2), mu(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mc, alphas, 0, runmw, g, gp, alpha, beta, runmc, mixNeut, 2, 1, onetothree, 'u');
     neutralino2amplitudeneut1ssbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), md(1,2), md(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(1), ms, alphas, 0, runmw, g, gp, alpha, beta, runms, mixNeut, 2, 1, onetothree, 'd');
     neutralino2amplitudeneut1ttbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), mu(1,3), mu(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mt, alphas, thetat, runmw, g, gp, alpha, beta, runmt, mixNeut, 2, 1, onetothree, 'u');
     neutralino2amplitudeneut1bbbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), md(1,3), md(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mb, alphas, thetab, runmw, g, gp, alpha, beta, runmb, mixNeut, 2, 1, onetothree, 'd');
     neutralino2amplitudeneut1eebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), me(1,1), me(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mel, alphas, 0, runmw, g, gp, alpha, beta, runmel, mixNeut, 2, 1, onetothree, 'l');
     neutralino2amplitudeneut1mumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), me(1,2), me(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mmu, alphas, 0, runmw, g, gp, alpha, beta, runmmu, mixNeut, 2, 1, onetothree, 'l');

										     neutralino2amplitudeneut1nuenuebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), msnu(1), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(1), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 2, 1, onetothree, 'n');
     neutralino2amplitudeneut1numunumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), msnu(2), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(1), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 2, 1, onetothree, 'n');
     neutralino2amplitudeneut1nutaunutaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(2), msnu(3), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(1), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 2, 1, onetothree, 'n'); ///Note set msf(2) very large as there is no msnuR so need this intermediate to decouple and not be present
     neutralino2amplitudechargino1udbar = neutralinoamplitudedecaycharginoffprimebar (mneut(2), mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mch(1), mup, mdo, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 2, 1, onetothree, 'q', 'n');
     neutralino2amplitudechargino1csbar = neutralinoamplitudedecaycharginoffprimebar (mneut(2), mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mch(1), mc, ms, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 2, 1, onetothree, 'q', 'n');
     neutralino2amplitudechargino1enuebar = neutralinoamplitudedecaycharginoffprimebar (mneut(2), me(1,1), me(2,1), msnu(1), 100000000000, polemw, mHpm, mch(1), 0, mel, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 2, 1, onetothree, 'l', 'n');
     neutralino2amplitudechargino1munumubar = neutralinoamplitudedecaycharginoffprimebar (mneut(2), me(1,2), me(2,2), msnu(2), 100000000000, polemw, mHpm, mch(1), 0, mmu, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 2, 1, onetothree, 'l', 'n');
     neutralino2amplitudechargino1taunutaubar = neutralinoamplitudedecaycharginoffprimebar (mneut(2), me(1,3), me(2,3), msnu(3), 100000000000, polemw, mHpm, mch(1), 0, mtau, 0, thetatau-PI/2, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 2, 1, onetothree, 'l', 'n');
     neutralino2amplitudechargino2udbar = neutralinoamplitudedecaycharginoffprimebar (mneut(2), mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mch(2), mup, mdo, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 2, 2, onetothree, 'q', 'n');
     neutralino2amplitudechargino2csbar = neutralinoamplitudedecaycharginoffprimebar (mneut(2), mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mch(2), mc, ms, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 2, 2, onetothree, 'q', 'n');
     neutralino2amplitudechargino2enuebar = neutralinoamplitudedecaycharginoffprimebar (mneut(2), me(1,1), me(2,1), msnu(1), 100000000000, polemw, mHpm, mch(2), 0, mel, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 2, 2, onetothree, 'l', 'n');
     neutralino2amplitudechargino2munumubar = neutralinoamplitudedecaycharginoffprimebar (mneut(2), me(1,2), me(2,2), msnu(2), 100000000000, polemw, mHpm, mch(2), 0, mmu, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 2, 2, onetothree, 'l', 'n');
     neutralino2amplitudechargino2taunutaubar = neutralinoamplitudedecaycharginoffprimebar (mneut(2), me(1,3), me(2,3), msnu(3), 100000000000, polemw, mHpm, mch(2), 0, mtau, 0, thetatau-PI/2, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 2, 2, onetothree, 'l', 'n');
   }
   
   else if (nmssmIsIt == true){
     neutralino2amplitudeWbosonpluscharginoW1 = neutralinoamplitudedecaycharginoWNMSSM (mneut(2), MCH1, polemw, g, thetaL2, thetaR2, mixNeut, 2, 1);
     neutralino2amplitudeWbosonpluscharginoW2 = neutralinoamplitudedecaycharginoWNMSSM (mneut(2), MCH2, polemw, g, thetaL2, thetaR2, mixNeut, 2, 2);
     neutralino2amplitudeWbosonminuscharginoW1 = neutralino2amplitudeWbosonpluscharginoW1;
     neutralino2amplitudeWbosonminuscharginoW2 = neutralino2amplitudeWbosonpluscharginoW2;
     neutralino2amplitudeZbosonneutralino1 = neutralinoamplitudedecayneutralinoZNMSSM (mneut(2), mneut(1), polemz, g, gp, mixNeut, 2, 1);
     
     neutralino2amplitudeHpluscharginoW1 = neutralinoamplitudecharginoHpmNMSSM (mneut(2), MCH1, mHpm, g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 2, 1);
     neutralino2amplitudeHpluscharginoW2 = neutralinoamplitudecharginoHpmNMSSM (mneut(2), MCH2, mHpm, g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 2, 2);
     neutralino2amplitudeHminuscharginoW1 = neutralino2amplitudeHpluscharginoW1;
     neutralino2amplitudeHminuscharginoW2 = neutralino2amplitudeHpluscharginoW2;
     
     neutralino2amplitudehneutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(2), mneut(1), mh0(1), g, gp, lam, kappa, mixNeut, CPEMix, 2, 1, 1);
     neutralino2amplitudeHneutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(2), mneut(1), mh0(2), g, gp, lam, kappa, mixNeut, CPEMix, 2, 1, 2);
     neutralino2amplitudeH3neutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(2), mneut(1), mh0(3), g, gp, lam, kappa, mixNeut, CPEMix, 2, 1, 3);
     
     neutralino2amplitudeAneutralino1 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(2), mneut(1), mA0(1), g, gp, lam, kappa, mixNeut, CPOMix, 2, 1, 1);
     neutralino2amplitudeA2neutralino1 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(2), mneut(1), mA0(2), g, gp, lam, kappa, mixNeut, CPOMix, 2, 1, 2);
     
     neutralino2amplitudeuLubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), mu(1,1), mup, g, gp, mixNeut, 2, 'u', 'L');
     neutralino2amplitudeuRubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), mu(2,1), mup, g, gp, mixNeut, 2, 'u', 'R');
     neutralino2amplitudeuLbaru = neutralino2amplitudeuLubar;
     neutralino2amplitudeuRbaru = neutralino2amplitudeuRubar;
     neutralino2amplitudedLdbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), md(1,1), mdo, g, gp, mixNeut, 2, 'd', 'L');
     neutralino2amplitudedRdbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), md(2,1), mdo, g, gp, mixNeut, 2, 'd', 'R');
     neutralino2amplitudedLbard = neutralino2amplitudedLdbar;
     neutralino2amplitudedRbard = neutralino2amplitudedRdbar;
     neutralino2amplitudecLcbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), mu(1,2), mc, g, gp, mixNeut, 2, 'u', 'L');
     neutralino2amplitudecRcbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), mu(2,2), mc, g, gp, mixNeut, 2, 'u', 'R');
     neutralino2amplitudecLbarc = neutralino2amplitudecLcbar;
     neutralino2amplitudecRbarc = neutralino2amplitudecRcbar;
     neutralino2amplitudesLsbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), md(1,2), ms, g, gp, mixNeut, 2, 'd', 'L');
     neutralino2amplitudesRsbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), md(2,2), ms, g, gp, mixNeut, 2, 'd', 'R');
     neutralino2amplitudesLbars = neutralino2amplitudesLsbar;
     neutralino2amplitudesRbars = neutralino2amplitudesRsbar;
     neutralino2amplitudeeLebar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), me(1,1), mel, g, gp, mixNeut, 2, 'l', 'L');
     neutralino2amplitudeeRebar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), me(2,1), mel, g, gp, mixNeut, 2, 'l', 'R');
     neutralino2amplitudeeLbare = neutralino2amplitudeeLebar;
     neutralino2amplitudeeRbare = neutralino2amplitudeeRebar;
     neutralino2amplitudemuLmubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), me(1,2), mmu, g, gp, mixNeut, 2, 'l', 'L');
     neutralino2amplitudemuRmubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(2), me(2,2), mmu, g, gp, mixNeut, 2, 'l', 'R');
     neutralino2amplitudemuLbarmu = neutralino2amplitudemuLmubar;
     neutralino2amplitudemuRbarmu = neutralino2amplitudemuRmubar;
     
     neutralino2amplitudetopstop1bar = neutralinoamplitudestoptopNMSSM (mneut(2), mu(1,3), mt, g, gp, thetat, beta, runmw, mixNeut, runmt, 2, 1);
     neutralino2amplitudetopstop2bar = neutralinoamplitudestoptopNMSSM (mneut(2), mu(2,3), mt, g, gp, thetat, beta, runmw, mixNeut, runmt, 2, 2);
     neutralino2amplitudetopbarstop1 = neutralino2amplitudetopstop1bar;
     neutralino2amplitudetopbarstop2 = neutralino2amplitudetopstop2bar;
     neutralino2amplitudebottomsbottom1bar = neutralinoamplitudesbottombottomNMSSM (mneut(2), md(1,3), mb, g, gp, thetab, beta, runmw, mixNeut, runmb, 2, 1);
     neutralino2amplitudebottomsbottom2bar = neutralinoamplitudesbottombottomNMSSM (mneut(2), md(2,3), mb, g, gp, thetab, beta, runmw, mixNeut, runmb, 2, 2);
     neutralino2amplitudebottombarsbottom1 = neutralino2amplitudebottomsbottom1bar;
     neutralino2amplitudebottombarsbottom2 = neutralino2amplitudebottomsbottom2bar;
     neutralino2amplitudetaustau1bar = neutralinoamplitudestautauNMSSM (mneut(2), me(1,3), mtau, g, gp, thetatau, beta, runmw, mixNeut, runmtau, 2, 1);
     neutralino2amplitudetaustau2bar = neutralinoamplitudestautauNMSSM (mneut(2), me(2,3), mtau, g, gp, thetatau, beta, runmw, mixNeut, runmtau, 2, 2);
     neutralino2amplitudetaubarstau1 = neutralino2amplitudetaustau1bar;
     neutralino2amplitudetaubarstau2 = neutralino2amplitudetaustau2bar;
     neutralino2amplitudesnuenuebar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(2), msnu(1), 0, g, gp, mixNeut, 2);
     neutralino2amplitudesnuebarnue = neutralino2amplitudesnuenuebar;
     neutralino2amplitudesnumunumubar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(2), msnu(2), 0, g, gp, mixNeut, 2);
     neutralino2amplitudesnumubarnumu = neutralino2amplitudesnumunumubar;
     neutralino2amplitudenutausnutaubar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(2), msnu(3), 0, g, gp, mixNeut, 2);
     neutralino2amplitudenutaubarsnutau = neutralino2amplitudenutausnutaubar;
   }
    
   ParticleNeutralino2.Array_Decays[0][0] = -PDGup; ParticleNeutralino2.Array_Decays[0][1] = PDGsupL; ParticleNeutralino2.Array_Decays[0][2] = neutralino2amplitudeuLubar; ParticleNeutralino2.Array_Decays[0][3] = 2; ParticleNeutralino2.Array_Comments[0] = "# ~chi_20 -> ub ~u_L";
   ParticleNeutralino2.Array_Decays[1][0] = -PDGup; ParticleNeutralino2.Array_Decays[1][1] = PDGsupR; ParticleNeutralino2.Array_Decays[1][2] = neutralino2amplitudeuRubar; ParticleNeutralino2.Array_Decays[1][3] = 2; ParticleNeutralino2.Array_Comments[1] = "# ~chi_20 -> ub ~u_R";
   ParticleNeutralino2.Array_Decays[2][0] = PDGup; ParticleNeutralino2.Array_Decays[2][1] = -PDGsupL; ParticleNeutralino2.Array_Decays[2][2] = neutralino2amplitudeuLbaru; ParticleNeutralino2.Array_Decays[2][3] = 2; ParticleNeutralino2.Array_Comments[2] = "# ~chi_20 -> u ~u_L*";
   ParticleNeutralino2.Array_Decays[3][0] = PDGup; ParticleNeutralino2.Array_Decays[3][1] = -PDGsupR; ParticleNeutralino2.Array_Decays[3][2] = neutralino2amplitudeuRbaru; ParticleNeutralino2.Array_Decays[3][3] = 2; ParticleNeutralino2.Array_Comments[3] = "# ~chi_20 -> u ~u_R*";
   ParticleNeutralino2.Array_Decays[4][0] = -PDGdown; ParticleNeutralino2.Array_Decays[4][1] = PDGsdownL; ParticleNeutralino2.Array_Decays[4][2] = neutralino2amplitudedLdbar; ParticleNeutralino2.Array_Decays[4][3] = 2; ParticleNeutralino2.Array_Comments[4] = "# ~chi_20 -> db ~d_L";
   ParticleNeutralino2.Array_Decays[5][0] = -PDGdown; ParticleNeutralino2.Array_Decays[5][1] = PDGsdownR; ParticleNeutralino2.Array_Decays[5][2] = neutralino2amplitudedRdbar; ParticleNeutralino2.Array_Decays[5][3] = 2; ParticleNeutralino2.Array_Comments[5] = "# ~chi_20 -> db ~d_R";
   ParticleNeutralino2.Array_Decays[6][0] = PDGdown; ParticleNeutralino2.Array_Decays[6][1] = -PDGsdownL; ParticleNeutralino2.Array_Decays[6][2] = neutralino2amplitudedLbard; ParticleNeutralino2.Array_Decays[6][3] = 2; ParticleNeutralino2.Array_Comments[6] = "# ~chi_20 -> d ~d_L*";
   ParticleNeutralino2.Array_Decays[7][0] = PDGdown; ParticleNeutralino2.Array_Decays[7][1] = -PDGsdownR; ParticleNeutralino2.Array_Decays[7][2] = neutralino2amplitudedRbard; ParticleNeutralino2.Array_Decays[7][3] = 2; ParticleNeutralino2.Array_Comments[7] = "# ~chi_20 -> d ~d_R*";
   ParticleNeutralino2.Array_Decays[8][0] = -PDGcharm; ParticleNeutralino2.Array_Decays[8][1] = PDGscharmL; ParticleNeutralino2.Array_Decays[8][2] = neutralino2amplitudecLcbar; ParticleNeutralino2.Array_Decays[8][3] = 2; ParticleNeutralino2.Array_Comments[8] = "# ~chi_20 -> cb ~c_L";
   ParticleNeutralino2.Array_Decays[9][0] = -PDGcharm; ParticleNeutralino2.Array_Decays[9][1] = PDGscharmR; ParticleNeutralino2.Array_Decays[9][2] = neutralino2amplitudecRcbar; ParticleNeutralino2.Array_Decays[9][3] = 2; ParticleNeutralino2.Array_Comments[9] = "# ~chi_20 -> cb ~c_R";
   ParticleNeutralino2.Array_Decays[10][0] = PDGcharm; ParticleNeutralino2.Array_Decays[10][1] = -PDGscharmL; ParticleNeutralino2.Array_Decays[10][2] = neutralino2amplitudecLbarc; ParticleNeutralino2.Array_Decays[10][3] = 2; ParticleNeutralino2.Array_Comments[10] = "# ~chi_20 -> c ~c_L*";
   ParticleNeutralino2.Array_Decays[11][0] = PDGcharm; ParticleNeutralino2.Array_Decays[11][1] = -PDGscharmR; ParticleNeutralino2.Array_Decays[11][2] = neutralino2amplitudecRbarc; ParticleNeutralino2.Array_Decays[11][3] = 2; ParticleNeutralino2.Array_Comments[11] = "# ~chi_20 -> c ~c_R*";
   ParticleNeutralino2.Array_Decays[12][0] = -PDGstrange; ParticleNeutralino2.Array_Decays[12][1] = PDGsstrangeL; ParticleNeutralino2.Array_Decays[12][2] = neutralino2amplitudesLsbar; ParticleNeutralino2.Array_Decays[12][3] = 2; ParticleNeutralino2.Array_Comments[12] = "# ~chi_20 -> sb ~s_L";
   ParticleNeutralino2.Array_Decays[13][0] = -PDGstrange; ParticleNeutralino2.Array_Decays[13][1] = PDGsstrangeR; ParticleNeutralino2.Array_Decays[13][2] = neutralino2amplitudesRsbar; ParticleNeutralino2.Array_Decays[13][3] = 2; ParticleNeutralino2.Array_Comments[13] = "# ~chi_20 -> sb ~s_R";
   ParticleNeutralino2.Array_Decays[14][0] = PDGstrange; ParticleNeutralino2.Array_Decays[14][1] = -PDGsstrangeL; ParticleNeutralino2.Array_Decays[14][2] = neutralino2amplitudesLbars; ParticleNeutralino2.Array_Decays[14][3] = 2; ParticleNeutralino2.Array_Comments[14] = "# ~chi_20 -> s ~s_L*";
   ParticleNeutralino2.Array_Decays[15][0] = PDGstrange; ParticleNeutralino2.Array_Decays[15][1] = -PDGsstrangeR; ParticleNeutralino2.Array_Decays[15][2] = neutralino2amplitudesRbars; ParticleNeutralino2.Array_Decays[15][3] = 2; ParticleNeutralino2.Array_Comments[15] = "# ~chi_20 -> s ~s_R*";
   ParticleNeutralino2.Array_Decays[16][0] = -PDGelectron; ParticleNeutralino2.Array_Decays[16][1] = PDGselectronL; ParticleNeutralino2.Array_Decays[16][2] = neutralino2amplitudeeLebar; ParticleNeutralino2.Array_Decays[16][3] = 2; ParticleNeutralino2.Array_Comments[16] = "# ~chi_20 -> e+ ~e_L-";
   ParticleNeutralino2.Array_Decays[17][0] = -PDGelectron; ParticleNeutralino2.Array_Decays[17][1] = PDGselectronR; ParticleNeutralino2.Array_Decays[17][2] = neutralino2amplitudeeRebar; ParticleNeutralino2.Array_Decays[17][3] = 2; ParticleNeutralino2.Array_Comments[17] = "# ~chi_20 -> e+ ~e_R-";
   ParticleNeutralino2.Array_Decays[18][0] = PDGelectron; ParticleNeutralino2.Array_Decays[18][1] = -PDGselectronL; ParticleNeutralino2.Array_Decays[18][2] = neutralino2amplitudeeLbare; ParticleNeutralino2.Array_Decays[18][3] = 2; ParticleNeutralino2.Array_Comments[18] = "# ~chi_20 -> e- ~e_L+";
   ParticleNeutralino2.Array_Decays[19][0] = PDGelectron; ParticleNeutralino2.Array_Decays[19][1] = -PDGselectronR; ParticleNeutralino2.Array_Decays[19][2] = neutralino2amplitudeeRbare; ParticleNeutralino2.Array_Decays[19][3] = 2; ParticleNeutralino2.Array_Comments[19] = "# ~chi_20 -> e- ~e_R+";   
   ParticleNeutralino2.Array_Decays[20][0] = -PDGmuon; ParticleNeutralino2.Array_Decays[20][1] = PDGsmuonL; ParticleNeutralino2.Array_Decays[20][2] = neutralino2amplitudemuLmubar; ParticleNeutralino2.Array_Decays[20][3] = 2; ParticleNeutralino2.Array_Comments[20] = "# ~chi_20 -> mu+ ~mu_L-";
   ParticleNeutralino2.Array_Decays[21][0] = -PDGmuon; ParticleNeutralino2.Array_Decays[21][1] = PDGsmuonR; ParticleNeutralino2.Array_Decays[21][2] = neutralino2amplitudemuRmubar; ParticleNeutralino2.Array_Decays[21][3] = 2; ParticleNeutralino2.Array_Comments[21] = "# ~chi_20 -> mu+ ~mu_R-";
   ParticleNeutralino2.Array_Decays[22][0] = PDGmuon; ParticleNeutralino2.Array_Decays[22][1] = -PDGsmuonL; ParticleNeutralino2.Array_Decays[22][2] = neutralino2amplitudemuLbarmu; ParticleNeutralino2.Array_Decays[22][3] = 2; ParticleNeutralino2.Array_Comments[22] = "# ~chi_20 -> mu- ~mu_L+";
   ParticleNeutralino2.Array_Decays[23][0] = PDGmuon; ParticleNeutralino2.Array_Decays[23][1] = -PDGsmuonR; ParticleNeutralino2.Array_Decays[23][2] = neutralino2amplitudemuRbarmu; ParticleNeutralino2.Array_Decays[23][3] = 2; ParticleNeutralino2.Array_Comments[23] = "# ~chi_20 -> mu- ~mu_R+";
   ParticleNeutralino2.Array_Decays[24][0] = PDGnuelectron; ParticleNeutralino2.Array_Decays[24][1] = -PDGnuselectronL; ParticleNeutralino2.Array_Decays[24][2] = neutralino2amplitudesnuebarnue; ParticleNeutralino2.Array_Decays[24][3] = 2; ParticleNeutralino2.Array_Comments[24] = "# ~chi_20 -> nu_e ~nu_eL*";
   ParticleNeutralino2.Array_Decays[25][0] = -PDGnuelectron; ParticleNeutralino2.Array_Decays[25][1] = PDGnuselectronL; ParticleNeutralino2.Array_Decays[25][2] = neutralino2amplitudesnuenuebar; ParticleNeutralino2.Array_Decays[25][3] = 2; ParticleNeutralino2.Array_Comments[25] = "# ~chi_20 -> nu_eb ~nu_eL";
   ParticleNeutralino2.Array_Decays[26][0] = PDGnumuon; ParticleNeutralino2.Array_Decays[26][1] = -PDGnusmuonL; ParticleNeutralino2.Array_Decays[26][2] = neutralino2amplitudesnumubarnumu; ParticleNeutralino2.Array_Decays[26][3] = 2; ParticleNeutralino2.Array_Comments[26] = "# ~chi_20 -> nu_mu ~nu_muL*";
   ParticleNeutralino2.Array_Decays[27][0] = -PDGnumuon; ParticleNeutralino2.Array_Decays[27][1] = PDGnusmuonL; ParticleNeutralino2.Array_Decays[27][2] = neutralino2amplitudesnumunumubar; ParticleNeutralino2.Array_Decays[27][3] = 2; ParticleNeutralino2.Array_Comments[27] = "# ~chi_20 -> nu_mub ~nu_muL";
   ParticleNeutralino2.Array_Decays[28][0] = PDGtop; ParticleNeutralino2.Array_Decays[28][1] = -PDGstop1; ParticleNeutralino2.Array_Decays[28][2] = neutralino2amplitudetopstop1bar; ParticleNeutralino2.Array_Decays[28][3] = 2; ParticleNeutralino2.Array_Comments[28] = "# ~chi_20 -> t ~t_1*";
   ParticleNeutralino2.Array_Decays[29][0] = PDGtop; ParticleNeutralino2.Array_Decays[29][1] = -PDGstop2; ParticleNeutralino2.Array_Decays[29][2] = neutralino2amplitudetopstop2bar; ParticleNeutralino2.Array_Decays[29][3] = 2; ParticleNeutralino2.Array_Comments[29] = "# ~chi_20 -> t ~t_2*";
   ParticleNeutralino2.Array_Decays[30][0] = -PDGtop; ParticleNeutralino2.Array_Decays[30][1] = PDGstop1; ParticleNeutralino2.Array_Decays[30][2] = neutralino2amplitudetopbarstop1; ParticleNeutralino2.Array_Decays[30][3] = 2; ParticleNeutralino2.Array_Comments[30] = "# ~chi_20 -> tb ~t_1";
   ParticleNeutralino2.Array_Decays[31][0] = -PDGtop; ParticleNeutralino2.Array_Decays[31][1] = PDGstop2; ParticleNeutralino2.Array_Decays[31][2] = neutralino2amplitudetopbarstop2; ParticleNeutralino2.Array_Decays[31][3] = 2; ParticleNeutralino2.Array_Comments[31] = "# ~chi_20 -> tb ~t_2";
   ParticleNeutralino2.Array_Decays[32][0] = PDGbottom; ParticleNeutralino2.Array_Decays[32][1] = -PDGsbottom1; ParticleNeutralino2.Array_Decays[32][2] = neutralino2amplitudebottomsbottom1bar; ParticleNeutralino2.Array_Decays[32][3] = 2; ParticleNeutralino2.Array_Comments[32] = "# ~chi_20 -> b ~b_1*";
   ParticleNeutralino2.Array_Decays[33][0] = PDGbottom; ParticleNeutralino2.Array_Decays[33][1] = -PDGsbottom2; ParticleNeutralino2.Array_Decays[33][2] = neutralino2amplitudebottomsbottom2bar; ParticleNeutralino2.Array_Decays[33][3] = 2; ParticleNeutralino2.Array_Comments[33] = "# ~chi_20 -> b ~b_2*";
   ParticleNeutralino2.Array_Decays[34][0] = -PDGbottom; ParticleNeutralino2.Array_Decays[34][1] = PDGsbottom1; ParticleNeutralino2.Array_Decays[34][2] = neutralino2amplitudebottombarsbottom1; ParticleNeutralino2.Array_Decays[34][3] = 2; ParticleNeutralino2.Array_Comments[34] = "# ~chi_20 -> bb ~b_1";
   ParticleNeutralino2.Array_Decays[35][0] = -PDGbottom; ParticleNeutralino2.Array_Decays[35][1] = PDGsbottom2; ParticleNeutralino2.Array_Decays[35][2] = neutralino2amplitudebottombarsbottom2; ParticleNeutralino2.Array_Decays[35][3] = 2; ParticleNeutralino2.Array_Comments[35] = "# ~chi_20 -> bb ~b_2";
   ParticleNeutralino2.Array_Decays[36][0] = -PDGstau1; ParticleNeutralino2.Array_Decays[36][1] = PDGtau; ParticleNeutralino2.Array_Decays[36][2] = neutralino2amplitudetaustau1bar; ParticleNeutralino2.Array_Decays[36][3] = 2; ParticleNeutralino2.Array_Comments[36] = "# ~chi_20 -> tau- ~tau_1+";
   ParticleNeutralino2.Array_Decays[37][0] = -PDGstau2; ParticleNeutralino2.Array_Decays[37][1] = PDGtau; ParticleNeutralino2.Array_Decays[37][2] = neutralino2amplitudetaustau2bar; ParticleNeutralino2.Array_Decays[37][3] = 2; ParticleNeutralino2.Array_Comments[37] = "# ~chi_20 -> tau- ~tau_2+";
   ParticleNeutralino2.Array_Decays[38][0] = PDGstau1; ParticleNeutralino2.Array_Decays[38][1] = -PDGtau; ParticleNeutralino2.Array_Decays[38][2] = neutralino2amplitudetaubarstau1; ParticleNeutralino2.Array_Decays[38][3] = 2; ParticleNeutralino2.Array_Comments[38] = "# ~chi_20 -> tau+ ~tau_1-";
   ParticleNeutralino2.Array_Decays[39][0] = PDGstau2; ParticleNeutralino2.Array_Decays[39][1] = -PDGtau; ParticleNeutralino2.Array_Decays[39][2] = neutralino2amplitudetaubarstau2; ParticleNeutralino2.Array_Decays[39][3] = 2; ParticleNeutralino2.Array_Comments[39] = "# ~chi_20 -> tau+ ~tau_2-";
   ParticleNeutralino2.Array_Decays[40][0] = PDGnutau; ParticleNeutralino2.Array_Decays[40][1] = -PDGnustauL; ParticleNeutralino2.Array_Decays[40][2] = neutralino2amplitudenutausnutaubar; ParticleNeutralino2.Array_Decays[40][3] = 2; ParticleNeutralino2.Array_Comments[40] = "# ~chi_20 -> nu_tau ~nu_tauL*";
   ParticleNeutralino2.Array_Decays[41][0] = -PDGnutau; ParticleNeutralino2.Array_Decays[41][1] = PDGnustauL; ParticleNeutralino2.Array_Decays[41][2] = neutralino2amplitudenutaubarsnutau; ParticleNeutralino2.Array_Decays[41][3] = 2; ParticleNeutralino2.Array_Comments[41] = "# ~chi_20 -> nu_taub ~nu_tauL";
   ParticleNeutralino2.Array_Decays[42][0] = PDGWplus; ParticleNeutralino2.Array_Decays[42][1] = -PDGchargino1; ParticleNeutralino2.Array_Decays[42][2] = neutralino2amplitudeWbosonpluscharginoW1; ParticleNeutralino2.Array_Decays[42][3] = 2; ParticleNeutralino2.Array_Comments[42] = "# ~chi_20 -> W+ ~chi_1-";
   ParticleNeutralino2.Array_Decays[43][0] = PDGWplus; ParticleNeutralino2.Array_Decays[43][1] = -PDGchargino2; ParticleNeutralino2.Array_Decays[43][2] = neutralino2amplitudeWbosonpluscharginoW2; ParticleNeutralino2.Array_Decays[43][3] = 2; ParticleNeutralino2.Array_Comments[43] = "# ~chi_20 -> W+ ~chi_2-";
   ParticleNeutralino2.Array_Decays[44][0] = -PDGWplus; ParticleNeutralino2.Array_Decays[44][1] = PDGchargino1; ParticleNeutralino2.Array_Decays[44][2] = neutralino2amplitudeWbosonminuscharginoW1; ParticleNeutralino2.Array_Decays[44][3] = 2; ParticleNeutralino2.Array_Comments[44] = "# ~chi_20 -> W- ~chi_1+";
   ParticleNeutralino2.Array_Decays[45][0] = -PDGWplus; ParticleNeutralino2.Array_Decays[45][1] = PDGchargino2; ParticleNeutralino2.Array_Decays[45][2] = neutralino2amplitudeWbosonminuscharginoW2; ParticleNeutralino2.Array_Decays[45][3] = 2; ParticleNeutralino2.Array_Comments[45] = "# ~chi_20 -> W- ~chi_2+";
   ParticleNeutralino2.Array_Decays[46][0] = PDGHplus; ParticleNeutralino2.Array_Decays[46][1] = -PDGchargino1; ParticleNeutralino2.Array_Decays[46][2] = neutralino2amplitudeHpluscharginoW1; ParticleNeutralino2.Array_Decays[46][3] = 2; ParticleNeutralino2.Array_Comments[46] = "# ~chi_20 -> H+ ~chi_1-";
   ParticleNeutralino2.Array_Decays[47][0] = PDGHplus; ParticleNeutralino2.Array_Decays[47][1] = -PDGchargino2; ParticleNeutralino2.Array_Decays[47][2] = neutralino2amplitudeHpluscharginoW2; ParticleNeutralino2.Array_Decays[47][3] = 2; ParticleNeutralino2.Array_Comments[47] = "# ~chi_20 -> H+ ~chi_2-";
   ParticleNeutralino2.Array_Decays[48][0] = -PDGHplus; ParticleNeutralino2.Array_Decays[48][1] = PDGchargino1; ParticleNeutralino2.Array_Decays[48][2] = neutralino2amplitudeHminuscharginoW1; ParticleNeutralino2.Array_Decays[48][3] = 2; ParticleNeutralino2.Array_Comments[48] = "# ~chi_20 -> H- ~chi1+";
   ParticleNeutralino2.Array_Decays[49][0] = -PDGHplus; ParticleNeutralino2.Array_Decays[49][1] = PDGchargino2; ParticleNeutralino2.Array_Decays[49][2] = neutralino2amplitudeHminuscharginoW2; ParticleNeutralino2.Array_Decays[49][3] = 2; ParticleNeutralino2.Array_Comments[49] = "# ~chi_20 -> H- ~chi_2+";
   ParticleNeutralino2.Array_Decays[50][0] = PDGZboson; ParticleNeutralino2.Array_Decays[50][1] = PDGneutralino1; ParticleNeutralino2.Array_Decays[50][2] = neutralino2amplitudeZbosonneutralino1; ParticleNeutralino2.Array_Decays[50][3] = 2; ParticleNeutralino2.Array_Comments[50] = "# ~chi_20 -> Z ~chi_10";
   ParticleNeutralino2.Array_Decays[51][0] = PDGZboson; ParticleNeutralino2.Array_Decays[51][1] = PDGneutralino3; ParticleNeutralino2.Array_Decays[51][2] = neutralino2amplitudeZbosonneutralino3; ParticleNeutralino2.Array_Decays[51][3] = 2; ParticleNeutralino2.Array_Comments[51] = "# ~chi_20 -> Z ~chi_30";
   ParticleNeutralino2.Array_Decays[52][0] = PDGZboson; ParticleNeutralino2.Array_Decays[52][1] = PDGneutralino4; ParticleNeutralino2.Array_Decays[52][2] = neutralino2amplitudeZbosonneutralino4; ParticleNeutralino2.Array_Decays[52][3] = 2; ParticleNeutralino2.Array_Comments[52] = "# ~chi_20 -> Z ~chi_40";
   ParticleNeutralino2.Array_Decays[53][0] = PDGh0; ParticleNeutralino2.Array_Decays[53][1] = PDGneutralino1; ParticleNeutralino2.Array_Decays[53][2] = neutralino2amplitudehneutralino1; ParticleNeutralino2.Array_Decays[53][3] = 2; ParticleNeutralino2.Array_Comments[53] = "# ~chi_20 -> h ~chi_10";
   ParticleNeutralino2.Array_Decays[54][0] = PDGh0; ParticleNeutralino2.Array_Decays[54][1] = PDGneutralino3; ParticleNeutralino2.Array_Decays[54][2] = neutralino2amplitudehneutralino3; ParticleNeutralino2.Array_Decays[54][3] = 2; ParticleNeutralino2.Array_Comments[54] = "# ~chi_20 -> h ~chi_30";
   ParticleNeutralino2.Array_Decays[55][0] = PDGh0; ParticleNeutralino2.Array_Decays[55][1] = PDGneutralino4; ParticleNeutralino2.Array_Decays[55][2] = neutralino2amplitudehneutralino4; ParticleNeutralino2.Array_Decays[55][3] = 2; ParticleNeutralino2.Array_Comments[55] = "# ~chi_20 -> h ~chi_40";
   ParticleNeutralino2.Array_Decays[56][0] = PDGH0; ParticleNeutralino2.Array_Decays[56][1] = PDGneutralino1; ParticleNeutralino2.Array_Decays[56][2] = neutralino2amplitudeHneutralino1; ParticleNeutralino2.Array_Decays[56][3] = 2; ParticleNeutralino2.Array_Comments[56] = "# ~chi_20 -> H ~chi_10";
   ParticleNeutralino2.Array_Decays[57][0] = PDGH0; ParticleNeutralino2.Array_Decays[57][1] = PDGneutralino3; ParticleNeutralino2.Array_Decays[57][2] = neutralino2amplitudeHneutralino3; ParticleNeutralino2.Array_Decays[57][3] = 2; ParticleNeutralino2.Array_Comments[57] = "# ~chi_20 -> H ~chi_30";
   ParticleNeutralino2.Array_Decays[58][0] = PDGH0; ParticleNeutralino2.Array_Decays[58][1] = PDGneutralino4; ParticleNeutralino2.Array_Decays[58][2] = neutralino2amplitudeHneutralino4; ParticleNeutralino2.Array_Decays[58][3] = 2; ParticleNeutralino2.Array_Comments[58] = "# ~chi_20 -> H ~chi_40";
   ParticleNeutralino2.Array_Decays[59][0] = PDGA0; ParticleNeutralino2.Array_Decays[59][1] = PDGneutralino1; ParticleNeutralino2.Array_Decays[59][2] = neutralino2amplitudeAneutralino1; ParticleNeutralino2.Array_Decays[59][3] = 2; ParticleNeutralino2.Array_Comments[59] = "# ~chi_20 -> A ~chi_10";
   ParticleNeutralino2.Array_Decays[60][0] = PDGA0; ParticleNeutralino2.Array_Decays[60][1] = PDGneutralino3; ParticleNeutralino2.Array_Decays[60][2] = neutralino2amplitudeAneutralino3; ParticleNeutralino2.Array_Decays[60][3] = 2; ParticleNeutralino2.Array_Comments[60] = "# ~chi_20 -> A ~chi_30";
   ParticleNeutralino2.Array_Decays[61][0] = PDGA0; ParticleNeutralino2.Array_Decays[61][1] = PDGneutralino4; ParticleNeutralino2.Array_Decays[61][2] = neutralino2amplitudeAneutralino4; ParticleNeutralino2.Array_Decays[61][3] = 2; ParticleNeutralino2.Array_Comments[61] = "# ~chi_20 -> A ~chi_40";
   
   ParticleNeutralino2.Array_Decays[62][0] = PDGH3; ParticleNeutralino2.Array_Decays[62][1] = PDGneutralino1; ParticleNeutralino2.Array_Decays[62][2] = neutralino2amplitudeH3neutralino1; ParticleNeutralino2.Array_Decays[62][3] = 2; ParticleNeutralino2.Array_Comments[62] = "# ~chi_20 -> H3 ~chi_10";
   ParticleNeutralino2.Array_Decays[63][0] = PDGA2; ParticleNeutralino2.Array_Decays[63][1] = PDGneutralino1; ParticleNeutralino2.Array_Decays[63][2] = neutralino2amplitudeA2neutralino1; ParticleNeutralino2.Array_Decays[63][3] = 2; ParticleNeutralino2.Array_Comments[63] = "# ~chi_20 -> A2 ~chi_10";
   
   ParticleNeutralino2.Array_Decays[64][0] = PDGphoton; ParticleNeutralino2.Array_Decays[64][1] = PDGgravitino; ParticleNeutralino2.Array_Decays[64][2] = neutralino2amplitudephotongravitino; ParticleNeutralino2.Array_Decays[64][3] = 2; ParticleNeutralino2.Array_Decays[64][4]=0; ParticleNeutralino2.Array_Comments[64] = "# ~chi_20 -> gamma ~G";
   ParticleNeutralino2.Array_Decays[65][0] = PDGZboson; ParticleNeutralino2.Array_Decays[65][1] = PDGgravitino; ParticleNeutralino2.Array_Decays[65][2] = neutralino2amplitudeZgravitino; ParticleNeutralino2.Array_Decays[65][3] = 2; ParticleNeutralino2.Array_Decays[65][4]=0; ParticleNeutralino2.Array_Comments[65] = "# ~chi_20 -> Z ~G";
   ParticleNeutralino2.Array_Decays[66][0] = PDGh0; ParticleNeutralino2.Array_Decays[66][1] = PDGgravitino; ParticleNeutralino2.Array_Decays[66][2] = neutralino2amplitudehgravitino; ParticleNeutralino2.Array_Decays[66][3] = 2; ParticleNeutralino2.Array_Decays[66][4] = 0; ParticleNeutralino2.Array_Comments[66] = "# ~chi_20 -> h ~G";
   ParticleNeutralino2.Array_Decays[67][0] = PDGH0; ParticleNeutralino2.Array_Decays[67][1] = PDGgravitino; ParticleNeutralino2.Array_Decays[67][2] = neutralino2amplitudeHgravitino; ParticleNeutralino2.Array_Decays[67][3] = 2; ParticleNeutralino2.Array_Decays[67][4] = 0; ParticleNeutralino2.Array_Comments[67] = "# ~chi_20 -> H ~G";
   ParticleNeutralino2.Array_Decays[68][0] = PDGA0; ParticleNeutralino2.Array_Decays[68][1] = PDGgravitino; ParticleNeutralino2.Array_Decays[68][2] = neutralino2amplitudeAgravitino; ParticleNeutralino2.Array_Decays[68][3] = 2; ParticleNeutralino2.Array_Decays[68][4] = 0; ParticleNeutralino2.Array_Comments[68] = "# ~chi_20 -> A ~G";
   
   ParticleNeutralino2.Array_Decays[69][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[69][1] = PDGup; ParticleNeutralino2.Array_Decays[69][4] = -PDGup; ParticleNeutralino2.Array_Decays[69][2] = neutralino2amplitudeneut1uubar; ParticleNeutralino2.Array_Decays[69][3] = 3; ParticleNeutralino2.Array_Comments[69] = "# ~chi_20 -> ~chi_10 u ubar";
   ParticleNeutralino2.Array_Decays[70][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[70][1] = PDGdown; ParticleNeutralino2.Array_Decays[70][4] = -PDGdown; ParticleNeutralino2.Array_Decays[70][2] = neutralino2amplitudeneut1ddbar; ParticleNeutralino2.Array_Decays[70][3] = 3; ParticleNeutralino2.Array_Comments[70] = "# ~chi_20 -> ~chi_10 d dbar";
   ParticleNeutralino2.Array_Decays[71][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[71][1] = PDGcharm; ParticleNeutralino2.Array_Decays[71][4] = -PDGcharm; ParticleNeutralino2.Array_Decays[71][2] = neutralino2amplitudeneut1ccbar; ParticleNeutralino2.Array_Decays[71][3] = 3; ParticleNeutralino2.Array_Comments[71] = "# ~chi_20 -> ~chi_10 c cbar";
   ParticleNeutralino2.Array_Decays[72][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[72][1] = PDGstrange; ParticleNeutralino2.Array_Decays[72][4] = -PDGstrange; ParticleNeutralino2.Array_Decays[72][2] = neutralino2amplitudeneut1ssbar; ParticleNeutralino2.Array_Decays[72][3] = 3; ParticleNeutralino2.Array_Comments[72] = "# ~chi_20 -> ~chi_10 s sbar";
   ParticleNeutralino2.Array_Decays[73][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[73][1] = PDGtop; ParticleNeutralino2.Array_Decays[73][4] = -PDGtop; ParticleNeutralino2.Array_Decays[73][2] = neutralino2amplitudeneut1ttbar; ParticleNeutralino2.Array_Decays[73][3] = 3; ParticleNeutralino2.Array_Comments[73] = "# ~chi_20 -> ~chi_10 t tbar";
   ParticleNeutralino2.Array_Decays[74][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[74][1] = PDGbottom; ParticleNeutralino2.Array_Decays[74][4] = -PDGbottom; ParticleNeutralino2.Array_Decays[74][2] = neutralino2amplitudeneut1bbbar; ParticleNeutralino2.Array_Decays[74][3] = 3; ParticleNeutralino2.Array_Comments[74] = "# ~chi_20 -> ~chi_10 b bbar";
   ParticleNeutralino2.Array_Decays[75][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[75][1] = PDGelectron; ParticleNeutralino2.Array_Decays[75][4] = -PDGelectron; ParticleNeutralino2.Array_Decays[75][2] = neutralino2amplitudeneut1eebar; ParticleNeutralino2.Array_Decays[75][3] = 3; ParticleNeutralino2.Array_Comments[75] = "# ~chi_20 -> ~chi_10 e- e+";
   ParticleNeutralino2.Array_Decays[76][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[76][1] = PDGmuon; ParticleNeutralino2.Array_Decays[76][4] = -PDGmuon; ParticleNeutralino2.Array_Decays[76][2] = neutralino2amplitudeneut1mumubar; ParticleNeutralino2.Array_Decays[76][3] = 3; ParticleNeutralino2.Array_Comments[76] = "# ~chi_20 -> ~chi_10 mu- mu+";
   ParticleNeutralino2.Array_Decays[77][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[77][1] = PDGtau; ParticleNeutralino2.Array_Decays[77][4] = -PDGtau; ParticleNeutralino2.Array_Decays[77][2] = neutralino2amplitudeneut1tautaubar; ParticleNeutralino2.Array_Decays[77][3] = 3; ParticleNeutralino2.Array_Comments[77] = "# ~chi_20 -> ~chi_10 tau- tau+";
   ParticleNeutralino2.Array_Decays[78][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[78][1] = PDGnuelectron; ParticleNeutralino2.Array_Decays[78][4] = -PDGnuelectron; ParticleNeutralino2.Array_Decays[78][2] = neutralino2amplitudeneut1nuenuebar; ParticleNeutralino2.Array_Decays[78][3] = 3; ParticleNeutralino2.Array_Comments[78] = "# ~chi_20 -> ~chi_10 nue nuebar";
   ParticleNeutralino2.Array_Decays[79][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[79][1] = PDGnumuon; ParticleNeutralino2.Array_Decays[79][4] = -PDGnumuon; ParticleNeutralino2.Array_Decays[79][2] = neutralino2amplitudeneut1numunumubar; ParticleNeutralino2.Array_Decays[79][3] = 3; ParticleNeutralino2.Array_Comments[79] = "# ~chi_20 -> ~chi_10 numu numubar";
   ParticleNeutralino2.Array_Decays[80][0] = PDGneutralino1; ParticleNeutralino2.Array_Decays[80][1] = PDGnutau; ParticleNeutralino2.Array_Decays[80][4] = -PDGnutau; ParticleNeutralino2.Array_Decays[80][2] = neutralino2amplitudeneut1nutaunutaubar; ParticleNeutralino2.Array_Decays[80][3] = 3; ParticleNeutralino2.Array_Comments[80] = "# ~chi_20 -> ~chi_10 nutau nutaubar";
   
   ParticleNeutralino2.Array_Decays[81][0] = PDGchargino1; ParticleNeutralino2.Array_Decays[81][1] = PDGup; ParticleNeutralino2.Array_Decays[81][4] = -PDGdown; ParticleNeutralino2.Array_Decays[81][2] = neutralino2amplitudechargino1udbar; ParticleNeutralino2.Array_Decays[81][3] = 3; ParticleNeutralino2.Array_Comments[81] = "# ~chi_20 -> chi_1- u db";
   ParticleNeutralino2.Array_Decays[82][0] = PDGchargino1; ParticleNeutralino2.Array_Decays[82][1] = PDGcharm; ParticleNeutralino2.Array_Decays[82][4] = -PDGstrange; ParticleNeutralino2.Array_Decays[82][2] = neutralino2amplitudechargino1csbar; ParticleNeutralino2.Array_Decays[82][3] = 3; ParticleNeutralino2.Array_Comments[82] = "# ~chi_20 -> chi_1- c sb";
   ParticleNeutralino2.Array_Decays[83][0] = PDGchargino1; ParticleNeutralino2.Array_Decays[83][1] = PDGnuelectron; ParticleNeutralino2.Array_Decays[83][4] = -PDGelectron; ParticleNeutralino2.Array_Decays[83][2] = neutralino2amplitudechargino1enuebar; ParticleNeutralino2.Array_Decays[83][3] = 3; ParticleNeutralino2.Array_Comments[83] = "# ~chi_20 -> chi_1- nu_e eb";
   ParticleNeutralino2.Array_Decays[84][0] = PDGchargino1; ParticleNeutralino2.Array_Decays[84][1] = PDGnumuon; ParticleNeutralino2.Array_Decays[84][4] = -PDGmuon; ParticleNeutralino2.Array_Decays[84][2] = neutralino2amplitudechargino1munumubar; ParticleNeutralino2.Array_Decays[84][3] = 3; ParticleNeutralino2.Array_Comments[84] = "# ~chi_20 -> chi_1- nu_mu mub";
   ParticleNeutralino2.Array_Decays[85][0] = PDGchargino1; ParticleNeutralino2.Array_Decays[85][1] = PDGnutau; ParticleNeutralino2.Array_Decays[85][4] = -PDGtau; ParticleNeutralino2.Array_Decays[85][2] = neutralino2amplitudechargino1taunutaubar; ParticleNeutralino2.Array_Decays[85][3] = 3; ParticleNeutralino2.Array_Comments[85] = "# ~chi_20 -> chi_1- nu_tau taub";
   ParticleNeutralino2.Array_Decays[86][0] = PDGchargino2; ParticleNeutralino2.Array_Decays[86][1] = PDGup; ParticleNeutralino2.Array_Decays[86][4] = -PDGdown; ParticleNeutralino2.Array_Decays[86][2] = neutralino2amplitudechargino2udbar; ParticleNeutralino2.Array_Decays[86][3] = 3; ParticleNeutralino2.Array_Comments[86] = "# ~chi_20 -> chi_2- u dbar";
   ParticleNeutralino2.Array_Decays[87][0] = PDGchargino2; ParticleNeutralino2.Array_Decays[87][1] = PDGcharm; ParticleNeutralino2.Array_Decays[87][4] = -PDGstrange; ParticleNeutralino2.Array_Decays[87][2] = neutralino2amplitudechargino2csbar; ParticleNeutralino2.Array_Decays[87][3] = 3; ParticleNeutralino2.Array_Comments[87] = "# ~chi_20 -> chi_2- c sbar";
   ParticleNeutralino2.Array_Decays[88][0] = PDGchargino2; ParticleNeutralino2.Array_Decays[88][1] = PDGnuelectron; ParticleNeutralino2.Array_Decays[88][4] = -PDGelectron; ParticleNeutralino2.Array_Decays[88][2] = neutralino2amplitudechargino2enuebar; ParticleNeutralino2.Array_Decays[88][3] = 3; ParticleNeutralino2.Array_Comments[88] = "# ~chi_20 -> chi_2- nu_e eb";
   ParticleNeutralino2.Array_Decays[89][0] = PDGchargino2; ParticleNeutralino2.Array_Decays[89][1] = PDGnumuon; ParticleNeutralino2.Array_Decays[89][4] = -PDGmuon; ParticleNeutralino2.Array_Decays[89][2] = neutralino2amplitudechargino2munumubar; ParticleNeutralino2.Array_Decays[89][3] = 3; ParticleNeutralino2.Array_Comments[89] = "# ~chi_20 -> chi_2- nu_mu mub";
   ParticleNeutralino2.Array_Decays[90][0] = PDGchargino2; ParticleNeutralino2.Array_Decays[90][1] = PDGnutau; ParticleNeutralino2.Array_Decays[90][4] = -PDGtau; ParticleNeutralino2.Array_Decays[90][2] = neutralino2amplitudechargino2taunutaubar; ParticleNeutralino2.Array_Decays[90][3] = 3; ParticleNeutralino2.Array_Comments[90] = "# ~chi_20 -> chi_2- nu_tau taubar";

   for(int i = 0; i<ParticleNeutralino2.No_of_Decays; i++) {
     if (ParticleNeutralino2.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleNeutralino2.Array_Comments[i] << " is negative = " << ParticleNeutralino2.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleNeutralino2.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }
   
   double Neut2_No_1to2_Decays = 0;
   
   Neut2_No_1to2_Decays = ParticleNeutralino2.No_1to2_Decays + ParticleNeutralino2.No_grav_Decays + ParticleNeutralino2.No_NMSSM_Decays;
   
   for (int j = 0; j<Neut2_No_1to2_Decays; j++) {
     ParticleNeutralino2.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<Neut2_No_1to2_Decays; j++) {
     ParticleNeutralino2.two_width = ParticleNeutralino2.two_width + ParticleNeutralino2.Array_Decays[j][2];
   }
   for (int j=Neut2_No_1to2_Decays; j<ParticleNeutralino2.No_of_Decays; j++) {
     ParticleNeutralino2.three_width = ParticleNeutralino2.three_width + ParticleNeutralino2.Array_Decays[j][2];
   }
   
   if ( ParticleNeutralino2.three_width != ParticleNeutralino2.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for neutralino 2 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleNeutralino2.No_of_Decays = Neut2_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleNeutralino2.total_width = ParticleNeutralino2.two_width;
     }
   else {
     ParticleNeutralino2.total_width = ParticleNeutralino2.two_width + ParticleNeutralino2.three_width;
   }
   
   if ( ParticleNeutralino2.total_width != ParticleNeutralino2.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleNeutralino2.No_of_Decays; i++) {
       //   fout << i << " " << ParticleNeutralino2.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in Neutralino2 total width \n");
     }
   
 }

///Neutralino3 Decays
 
 double neutralino3amplitudeuLubar=0, neutralino3amplitudeuRubar=0, neutralino3amplitudeuLbaru=0, neutralino3amplitudeuRbaru=0, neutralino3amplitudedLdbar=0, neutralino3amplitudedRdbar=0, neutralino3amplitudedLbard=0, neutralino3amplitudedRbard=0, neutralino3amplitudecLcbar=0, neutralino3amplitudecRcbar=0, neutralino3amplitudecLbarc=0, neutralino3amplitudecRbarc=0, neutralino3amplitudesLsbar=0, neutralino3amplitudesRsbar=0, neutralino3amplitudesLbars=0, neutralino3amplitudesRbars=0, neutralino3amplitudeeLebar=0, neutralino3amplitudeeRebar=0, neutralino3amplitudeeLbare=0, neutralino3amplitudeeRbare=0, neutralino3amplitudemuLmubar=0, neutralino3amplitudemuRmubar=0, neutralino3amplitudemuLbarmu=0, neutralino3amplitudemuRbarmu=0, neutralino3amplitudesnuenuebar=0, neutralino3amplitudesnuebarnue=0, neutralino3amplitudesnumunumubar=0, neutralino3amplitudesnumubarnumu=0, neutralino3amplitudetopstop1bar=0, neutralino3amplitudetopstop2bar=0, neutralino3amplitudetopbarstop1=0, neutralino3amplitudetopbarstop2=0, neutralino3amplitudebottomsbottom1bar=0, neutralino3amplitudebottomsbottom2bar=0, neutralino3amplitudebottombarsbottom1=0, neutralino3amplitudebottombarsbottom2=0, neutralino3amplitudetaustau1bar=0, neutralino3amplitudetaustau2bar=0, neutralino3amplitudetaubarstau1=0, neutralino3amplitudetaubarstau2=0, neutralino3amplitudenutausnutaubar=0, neutralino3amplitudenutaubarsnutau=0, neutralino3amplitudeWbosonpluscharginoW1=0, neutralino3amplitudeWbosonpluscharginoW2=0, neutralino3amplitudeWbosonminuscharginoW1=0, neutralino3amplitudeWbosonminuscharginoW2=0, neutralino3amplitudeHpluscharginoW1=0, neutralino3amplitudeHpluscharginoW2=0, neutralino3amplitudeHminuscharginoW1=0, neutralino3amplitudeHminuscharginoW2=0, neutralino3amplitudeZbosonneutralino1=0, neutralino3amplitudeZbosonneutralino2=0, neutralino3amplitudeZbosonneutralino4=0, neutralino3amplitudehneutralino1=0, neutralino3amplitudehneutralino2=0, neutralino3amplitudehneutralino4=0, neutralino3amplitudeHneutralino1=0, neutralino3amplitudeHneutralino2=0, neutralino3amplitudeHneutralino4=0, neutralino3amplitudeAneutralino1=0, neutralino3amplitudeAneutralino2=0, neutralino3amplitudeAneutralino4=0, neutralino3amplitudephotongravitino=0, neutralino3amplitudeZgravitino=0, neutralino3amplitudehgravitino=0, neutralino3amplitudeHgravitino=0, neutralino3amplitudeAgravitino=0;

 double neutralino3amplitudeneut1uubar=0, neutralino3amplitudeneut1ddbar=0, neutralino3amplitudeneut1ccbar=0, neutralino3amplitudeneut1ssbar=0, neutralino3amplitudeneut1ttbar=0, neutralino3amplitudeneut1bbbar=0, neutralino3amplitudeneut1eebar=0, neutralino3amplitudeneut1mumubar=0, neutralino3amplitudeneut1tautaubar=0, neutralino3amplitudeneut1nuenuebar=0, neutralino3amplitudeneut1numunumubar=0, neutralino3amplitudeneut1nutaunutaubar=0, neutralino3amplitudeneut2uubar=0, neutralino3amplitudeneut2ddbar=0, neutralino3amplitudeneut2ccbar=0, neutralino3amplitudeneut2ssbar=0, neutralino3amplitudeneut2ttbar=0, neutralino3amplitudeneut2bbbar=0, neutralino3amplitudeneut2eebar=0, neutralino3amplitudeneut2mumubar=0, neutralino3amplitudeneut2tautaubar=0, neutralino3amplitudeneut2nuenuebar=0, neutralino3amplitudeneut2numunumubar=0, neutralino3amplitudeneut2nutaunutaubar=0, neutralino3amplitudechargino1udbar=0, neutralino3amplitudechargino1csbar=0, neutralino3amplitudechargino1enuebar=0, neutralino3amplitudechargino1munumubar=0, neutralino3amplitudechargino1taunutaubar=0, neutralino3amplitudechargino2udbar=0, neutralino3amplitudechargino2csbar=0, neutralino3amplitudechargino2enuebar=0, neutralino3amplitudechargino2munumubar=0, neutralino3amplitudechargino2taunutaubar=0;

 double neutralino3amplitudeH3neutralino1 = 0, neutralino3amplitudeH3neutralino2 = 0, neutralino3amplitudeA2neutralino1 = 0, neutralino3amplitudeA2neutralino2 = 0;

 if (flagneut3 == 1) {
   if (nmssmIsIt == false) {
     neutralino3amplitudeuLubar = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mup, mu(1,1), g, gp, mixNeut, 1, 'L', 3);
     neutralino3amplitudeuRubar = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mup, mu(2,1), g, gp, mixNeut, 1, 'R', 3);
     neutralino3amplitudeuLbaru = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mup, mu(1,1), g, gp, mixNeut, 1, 'L', 3); 
     neutralino3amplitudeuRbaru = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mup, mu(2,1), g, gp, mixNeut, 1, 'R', 3);
     neutralino3amplitudedLdbar = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mdo, md(1,1), g, gp, mixNeut, 2, 'L', 3);
     neutralino3amplitudedRdbar = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mdo, md(2,1), g, gp, mixNeut, 2, 'R', 3);
     neutralino3amplitudedLbard = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mdo, md(1,1), g, gp, mixNeut, 2, 'L', 3);
     neutralino3amplitudedRbard = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mdo, md(2,1), g, gp, mixNeut, 2, 'R', 3);
     neutralino3amplitudecLcbar = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mc, mu(1,2), g, gp, mixNeut, 1, 'L', 3);
     neutralino3amplitudecRcbar = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mc, mu(2,2), g, gp, mixNeut, 1, 'R', 3);
     neutralino3amplitudecLbarc = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mc, mu(1,2), g, gp, mixNeut, 1, 'L', 3);
     neutralino3amplitudecRbarc = neutralinoamplitudedecayquarksquarkLorR (mneut(3), mc, mu(2,2), g, gp, mixNeut, 1, 'R', 3);
     neutralino3amplitudesLsbar = neutralinoamplitudedecayquarksquarkLorR (mneut(3), ms, md(1,2), g, gp, mixNeut, 2, 'L', 3);
     neutralino3amplitudesRsbar = neutralinoamplitudedecayquarksquarkLorR (mneut(3), ms, md(2,2), g, gp, mixNeut, 2, 'R', 3);
     neutralino3amplitudesLbars = neutralinoamplitudedecayquarksquarkLorR (mneut(3), ms, md(1,2), g, gp, mixNeut, 2, 'L', 3);
     neutralino3amplitudesRbars = neutralinoamplitudedecayquarksquarkLorR (mneut(3), ms, md(2,2), g, gp, mixNeut, 2, 'R', 3);
     neutralino3amplitudeeLebar = neutralinoamplitudedecayleptonsleptonLorR (mneut(3), mel, me(1,1), g, gp, mixNeut, 'L', 3);
     neutralino3amplitudeeRebar = neutralinoamplitudedecayleptonsleptonLorR (mneut(3), mel, me(2,1), g, gp, mixNeut, 'R', 3);
     neutralino3amplitudeeLbare = neutralinoamplitudedecayleptonsleptonLorR (mneut(3), mel, me(1,1), g, gp, mixNeut, 'L', 3);
     neutralino3amplitudeeRbare = neutralinoamplitudedecayleptonsleptonLorR (mneut(3), mel, me(2,1), g, gp, mixNeut, 'R', 3);
     neutralino3amplitudemuLmubar = neutralinoamplitudedecayleptonsleptonLorR (mneut(3), mmu, me(1,2), g, gp, mixNeut, 'L', 3);
     neutralino3amplitudemuRmubar = neutralinoamplitudedecayleptonsleptonLorR (mneut(3), mmu, me(2,2), g, gp, mixNeut, 'R', 3);
     neutralino3amplitudemuLbarmu = neutralinoamplitudedecayleptonsleptonLorR (mneut(3), mmu, me(1,2), g, gp, mixNeut, 'L', 3);
     neutralino3amplitudemuRbarmu = neutralinoamplitudedecayleptonsleptonLorR (mneut(3), mmu, me(2,2), g, gp, mixNeut, 'R', 3);
     neutralino3amplitudesnuenuebar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(3), 0, msnu(1), g, gp, mixNeut, 3);
     neutralino3amplitudesnuebarnue = neutralinoamplitudedecayneutrinosneutrinoL (mneut(3), 0, msnu(1), g, gp, mixNeut, 3);
     neutralino3amplitudesnumunumubar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(3), 0, msnu(2), g, gp, mixNeut, 3);
     neutralino3amplitudesnumubarnumu = neutralinoamplitudedecayneutrinosneutrinoL (mneut(3), 0, msnu(2), g, gp, mixNeut, 3);
     neutralino3amplitudetopstop1bar = neutralinoamplitudedecaysquark3quarkmix (mneut(3), mt, mu(1,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 3);
     neutralino3amplitudetopstop2bar = neutralinoamplitudedecaysquark3quarkmix (mneut(3), mt, mu(2,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 3);
     neutralino3amplitudetopbarstop1 = neutralinoamplitudedecaysquark3quarkmix (mneut(3), mt, mu(1,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 3);
     neutralino3amplitudetopbarstop2 = neutralinoamplitudedecaysquark3quarkmix (mneut(3), mt, mu(2,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 3);
     neutralino3amplitudebottomsbottom1bar = neutralinoamplitudedecaysquark3quarkmix (mneut(3), mb, md(1,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 3);
     neutralino3amplitudebottomsbottom2bar = neutralinoamplitudedecaysquark3quarkmix (mneut(3), mb, md(2,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 3);
     neutralino3amplitudebottombarsbottom1 = neutralinoamplitudedecaysquark3quarkmix (mneut(3), mb, md(1,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 3);
     neutralino3amplitudebottombarsbottom2 = neutralinoamplitudedecaysquark3quarkmix (mneut(3), mb, md(2,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 3);
     neutralino3amplitudetaustau1bar = neutralinoamplitudedecaystautau (mneut(3), mtau, me(1,3), runmw, thetatau, beta, mixNeut, g, gp, 1, 3);
     neutralino3amplitudetaustau2bar = neutralinoamplitudedecaystautau (mneut(3), mtau, me(2,3), runmw, thetatau, beta, mixNeut, g, gp, 2, 3);
     neutralino3amplitudetaubarstau1 = neutralinoamplitudedecaystautau (mneut(3), mtau, me(1,3), runmw, thetatau, beta, mixNeut, g, gp, 1, 3);
     neutralino3amplitudetaubarstau2 = neutralinoamplitudedecaystautau (mneut(3), mtau, me(2,3), runmw, thetatau, beta, mixNeut, g, gp, 2, 3);
     neutralino3amplitudenutausnutaubar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(3), 0, msnu(3), g, gp, mixNeut, 3);
     neutralino3amplitudenutaubarsnutau = neutralinoamplitudedecayneutrinosneutrinoL (mneut(3), 0, msnu(3), g, gp, mixNeut, 3);
     neutralino3amplitudeWbosonpluscharginoW1 = neutralinoamplitudedecaycharginoWboson (mneut(3), polemw, MCH1, g, thetaL2, thetaR2, mixNeut, 3, 1);
     neutralino3amplitudeWbosonpluscharginoW2 = neutralinoamplitudedecaycharginoWboson (mneut(3), polemw, MCH2, g, thetaL2, thetaR2, mixNeut, 3, 2);
     neutralino3amplitudeWbosonminuscharginoW1 = neutralinoamplitudedecaycharginoWboson (mneut(3), polemw, MCH1, g, thetaL2, thetaR2, mixNeut, 3, 1);
     neutralino3amplitudeWbosonminuscharginoW2 = neutralinoamplitudedecaycharginoWboson (mneut(3), polemw, MCH2, g, thetaL2, thetaR2, mixNeut, 3, 2);
     neutralino3amplitudeHpluscharginoW1 = neutralinoamplitudedecaycharginoHplus (mneut(3), mHpm, MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 3, 1);
     neutralino3amplitudeHpluscharginoW2 = neutralinoamplitudedecaycharginoHplus (mneut(3), mHpm, MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 3, 2);
     neutralino3amplitudeHminuscharginoW1 = neutralinoamplitudedecaycharginoHplus (mneut(3), mHpm, MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 3, 1);
     neutralino3amplitudeHminuscharginoW2 = neutralinoamplitudedecaycharginoHplus (mneut(3), mHpm, MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 3, 2);

     neutralino3amplitudeZbosonneutralino1 = neutralinoamplitudedecayneutralinoZboson (mneut(3), polemz, mneut(1), g, gp, mixNeut, 3, 1);
     neutralino3amplitudeZbosonneutralino2 = neutralinoamplitudedecayneutralinoZboson (mneut(3), polemz, mneut(2), g, gp, mixNeut, 3, 2);
     neutralino3amplitudehneutralino1 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(3), mh0(1), mneut(1), g, gp, mixNeut, alpha, 3, 1, 'h');
     neutralino3amplitudehneutralino2 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(3), mh0(1), mneut(2), g, gp, mixNeut, alpha, 3, 2, 'h');
     neutralino3amplitudeHneutralino1 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(3), mh0(2), mneut(1), g, gp, mixNeut, alpha, 3, 1, 'H');
     neutralino3amplitudeHneutralino2 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(3), mh0(2), mneut(2), g, gp, mixNeut, alpha, 3, 2, 'H');
     neutralino3amplitudeAneutralino1 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(3), mA0(1), mneut(1), g, gp, mixNeut, beta, 3, 1, 'A');
     neutralino3amplitudeAneutralino2 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(3), mA0(1), mneut(2), g, gp, mixNeut, beta, 3, 2, 'A');
     
     neutralino3amplitudephotongravitino = neutralinoamplitudedecayphotongravitino(mneut(3), mgravitino, MPlreduced, mixNeut, g, gp, 3, gravonoff, neutNLSP);
     neutralino3amplitudeZgravitino = neutralinoamplitudedecayZgravitino(mneut(3), polemz, mgravitino, MPlreduced, mixNeut, g, gp, beta, 3, gravonoff, neutNLSP);
     neutralino3amplitudehgravitino = neutralinoamplitudedecayphigravitino(mneut(3), mh0(1), mgravitino, MPlreduced, mixNeut, alpha, beta, 3, gravonoff, 'h', neutNLSP);
     neutralino3amplitudeHgravitino = neutralinoamplitudedecayphigravitino(mneut(3), mh0(2), mgravitino, MPlreduced, mixNeut, alpha, beta, 3, gravonoff, 'H', neutNLSP);
     neutralino3amplitudeAgravitino = neutralinoamplitudedecayphigravitino(mneut(3), mA0(1), mgravitino, MPlreduced, mixNeut, alpha, beta, 3, gravonoff, 'A', neutNLSP);
     
     neutralino3amplitudeneut1uubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), mu(1,1), mu(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mup, alphas, 0, runmw, g, gp, alpha, beta, runmu, mixNeut, 3, 1, onetothree, 'u');
     neutralino3amplitudeneut1ddbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), md(1,1), md(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mdo, alphas, 0, runmw, g, gp, alpha, beta, runmd, mixNeut, 3, 1, onetothree, 'd');
     neutralino3amplitudeneut1ccbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), mu(1,2), mu(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mc, alphas, 0, runmw, g, gp, alpha, beta, runmc, mixNeut, 3, 1, onetothree, 'u');
     neutralino3amplitudeneut1ssbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), md(1,2), md(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(1), ms, alphas, 0, runmw, g, gp, alpha, beta, runms, mixNeut, 3, 1, onetothree, 'd');
     neutralino3amplitudeneut1ttbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), mu(1,3), mu(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mt, alphas, thetat, runmw, g, gp, alpha, beta, runmt, mixNeut, 3, 1, onetothree, 'u');
     neutralino3amplitudeneut1bbbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), md(1,3), md(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mb, alphas, thetab, runmw, g, gp, alpha, beta, runmb, mixNeut, 3, 1, onetothree, 'd');
     neutralino3amplitudeneut1eebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), me(1,1), me(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mel, alphas, 0, runmw, g, gp, alpha, beta, runmel, mixNeut, 3, 1, onetothree, 'l');
     neutralino3amplitudeneut1mumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), me(1,2), me(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mmu, alphas, 0, runmw, g, gp, alpha, beta, runmmu, mixNeut, 3, 1, onetothree, 'l');
     neutralino3amplitudeneut1tautaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), me(1,3), me(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mtau, alphas, thetatau-PI/2, runmw, g, gp, alpha, beta, runmtau, mixNeut, 3, 1, onetothree, 'l');
     neutralino3amplitudeneut1nuenuebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), msnu(1), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(1), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 3, 1, onetothree, 'n');
     neutralino3amplitudeneut1numunumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), msnu(2), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(1), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 3, 1, onetothree, 'n');
     neutralino3amplitudeneut1nutaunutaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), msnu(3), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(1), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 3, 1, onetothree, 'n');
     neutralino3amplitudeneut2uubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), mu(1,1), mu(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mup, alphas, 0, runmw, g, gp, alpha, beta, runmu, mixNeut, 3, 2, onetothree, 'u');
     neutralino3amplitudeneut2ddbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), md(1,1), md(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mdo, alphas, 0, runmw, g, gp, alpha, beta, runmd, mixNeut, 3, 2, onetothree, 'd');
     neutralino3amplitudeneut2ccbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), mu(1,2), mu(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mc, alphas, 0, runmw, g, gp, alpha, beta, runmc, mixNeut, 3, 2, onetothree, 'u');
     neutralino3amplitudeneut2ssbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), md(1,2), md(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(2), ms, alphas, 0, runmw, g, gp, alpha, beta, runms, mixNeut, 3, 2, onetothree, 'd');
     neutralino3amplitudeneut2ttbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), mu(1,3), mu(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mt, alphas, thetat, runmw, g, gp, alpha, beta, runmt, mixNeut, 3, 2, onetothree, 'u');
     neutralino3amplitudeneut2bbbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), md(1,3), md(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mb, alphas, thetab, runmw, g, gp, alpha, beta, runmb, mixNeut, 3, 2, onetothree, 'd');
     neutralino3amplitudeneut2eebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), me(1,1), me(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mel, alphas, 0, runmw, g, gp, alpha, beta, runmel, mixNeut, 3, 2, onetothree, 'l');
     neutralino3amplitudeneut2mumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), me(1,2), me(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mmu, alphas, 0, runmw, g, gp, alpha, beta, runmmu, mixNeut, 3, 2, onetothree, 'l');
     neutralino3amplitudeneut2tautaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), me(1,3), me(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mtau, alphas, thetatau-PI/2, runmw, g, gp, alpha, beta, runmtau, mixNeut, 3, 2, onetothree, 'l');
     neutralino3amplitudeneut2nuenuebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), msnu(1), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(2), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 3, 2, onetothree, 'n');
     neutralino3amplitudeneut2numunumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), msnu(2), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(2), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 3, 2, onetothree, 'n');
     neutralino3amplitudeneut2nutaunutaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(3), msnu(3), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(2), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 3, 2, onetothree, 'n');
     
     neutralino3amplitudechargino1udbar = neutralinoamplitudedecaycharginoffprimebar (mneut(3), mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mch(1), mup, mdo, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 3, 1, onetothree, 'q', 'n');
     neutralino3amplitudechargino1csbar = neutralinoamplitudedecaycharginoffprimebar (mneut(3), mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mch(1), mc, ms, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 3, 1, onetothree, 'q', 'n');
     neutralino3amplitudechargino1enuebar = neutralinoamplitudedecaycharginoffprimebar (mneut(3), me(1,1), me(2,1), msnu(1), 100000000000, polemw, mHpm, mch(1), 0, mel, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 3, 1, onetothree, 'l', 'n');
     neutralino3amplitudechargino1munumubar = neutralinoamplitudedecaycharginoffprimebar (mneut(3), me(1,2), me(2,2), msnu(2), 100000000000, polemw, mHpm, mch(1), 0, mmu, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 3, 1, onetothree, 'l', 'n');
     neutralino3amplitudechargino1taunutaubar = neutralinoamplitudedecaycharginoffprimebar (mneut(3), me(1,3), me(2,3), msnu(3), 100000000000, polemw, mHpm, mch(1), 0, mtau, 0, thetatau-PI/2, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 3, 1, onetothree, 'l', 'n');
     neutralino3amplitudechargino2udbar = neutralinoamplitudedecaycharginoffprimebar (mneut(3), mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mch(2), mup, mdo, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 3, 2, onetothree, 'q', 'n');
     neutralino3amplitudechargino2csbar = neutralinoamplitudedecaycharginoffprimebar (mneut(3), mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mch(2), mc, ms, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 3, 2, onetothree, 'q', 'n');
     neutralino3amplitudechargino2enuebar = neutralinoamplitudedecaycharginoffprimebar (mneut(3), me(1,1), me(2,1), msnu(1), 100000000000, polemw, mHpm, mch(2), 0, mel, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 3, 2, onetothree, 'l', 'n');
     neutralino3amplitudechargino2munumubar = neutralinoamplitudedecaycharginoffprimebar (mneut(3), me(1,2), me(2,2), msnu(2), 100000000000, polemw, mHpm, mch(2), 0, mmu, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 3, 2, onetothree, 'l', 'n');
     neutralino3amplitudechargino2taunutaubar = neutralinoamplitudedecaycharginoffprimebar (mneut(3), me(1,3), me(2,3), msnu(3), 100000000000, polemw, mHpm, mch(2), 0, mtau, 0, thetatau-PI/2, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 3, 2, onetothree, 'l', 'n');
     
   }
   else if (nmssmIsIt == true){
     neutralino3amplitudeWbosonpluscharginoW1 = neutralinoamplitudedecaycharginoWNMSSM (mneut(3), MCH1, polemw, g, thetaL2, thetaR2, mixNeut, 3, 1);
     neutralino3amplitudeWbosonpluscharginoW2 = neutralinoamplitudedecaycharginoWNMSSM (mneut(3), MCH2, polemw, g, thetaL2, thetaR2, mixNeut, 3, 2);
     neutralino3amplitudeWbosonminuscharginoW1 = neutralino3amplitudeWbosonpluscharginoW1;
     neutralino3amplitudeWbosonminuscharginoW2 = neutralino3amplitudeWbosonpluscharginoW2;
     neutralino3amplitudeZbosonneutralino1 = neutralinoamplitudedecayneutralinoZNMSSM (mneut(3), mneut(1), polemz, g, gp, mixNeut, 3, 1);
     neutralino3amplitudeZbosonneutralino2 = neutralinoamplitudedecayneutralinoZNMSSM (mneut(3), mneut(2), polemz, g, gp, mixNeut, 3, 2);
     
     neutralino3amplitudeHpluscharginoW1 = neutralinoamplitudecharginoHpmNMSSM (mneut(3), MCH1, mHpm, g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 3, 1);
     neutralino3amplitudeHpluscharginoW2 = neutralinoamplitudecharginoHpmNMSSM (mneut(3), MCH2, mHpm, g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 3, 2);
     neutralino3amplitudeHminuscharginoW1 = neutralino3amplitudeHpluscharginoW1;
     neutralino3amplitudeHminuscharginoW2 = neutralino3amplitudeHpluscharginoW2;
     
     neutralino3amplitudehneutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(3), mneut(1), mh0(1), g, gp, lam, kappa, mixNeut, CPEMix, 3, 1, 1);
     neutralino3amplitudeHneutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(3), mneut(1), mh0(2), g, gp, lam, kappa, mixNeut, CPEMix, 3, 1, 2);
     neutralino3amplitudeH3neutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(3), mneut(1), mh0(3), g, gp, lam, kappa, mixNeut, CPEMix, 3, 1, 3);
     neutralino3amplitudehneutralino2 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(3), mneut(2), mh0(1), g, gp, lam, kappa, mixNeut, CPEMix, 3, 2, 1);
     neutralino3amplitudeHneutralino2 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(3), mneut(2), mh0(2), g, gp, lam, kappa, mixNeut, CPEMix, 3, 2, 2);
     neutralino3amplitudeH3neutralino2 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(3), mneut(2), mh0(3), g, gp, lam, kappa, mixNeut, CPEMix, 3, 2, 3);
     
     neutralino3amplitudeAneutralino1 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(3), mneut(1), mA0(1), g, gp, lam, kappa, mixNeut, CPOMix, 3, 1, 1);
     neutralino3amplitudeA2neutralino1 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(3), mneut(1), mA0(2), g, gp, lam, kappa, mixNeut, CPOMix, 3, 1, 2);
     neutralino3amplitudeAneutralino2 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(3), mneut(2), mA0(1), g, gp, lam, kappa, mixNeut, CPOMix, 3, 2, 1);
     neutralino3amplitudeA2neutralino2 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(3), mneut(2), mA0(2), g, gp, lam, kappa, mixNeut, CPOMix, 3, 2, 2);
     
     neutralino3amplitudeuLubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), mu(1,1), mup, g, gp, mixNeut, 3, 'u', 'L');
     neutralino3amplitudeuRubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), mu(2,1), mup, g, gp, mixNeut, 3, 'u', 'R');
     neutralino3amplitudeuLbaru = neutralino3amplitudeuLubar;
     neutralino3amplitudeuRbaru = neutralino3amplitudeuRubar;
     neutralino3amplitudedLdbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), md(1,1), mdo, g, gp, mixNeut, 3, 'd', 'L');
     neutralino3amplitudedRdbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), md(2,1), mdo, g, gp, mixNeut, 3, 'd', 'R');
     neutralino3amplitudedLbard = neutralino3amplitudedLdbar;
     neutralino3amplitudedRbard = neutralino3amplitudedRdbar;
     neutralino3amplitudecLcbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), mu(1,2), mc, g, gp, mixNeut, 3, 'u', 'L');
     neutralino3amplitudecRcbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), mu(2,2), mc, g, gp, mixNeut, 3, 'u', 'R');
     neutralino3amplitudecLbarc = neutralino3amplitudecLcbar;
     neutralino3amplitudecRbarc = neutralino3amplitudecRcbar;
     neutralino3amplitudesLsbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), md(1,2), ms, g, gp, mixNeut, 3, 'd', 'L');
     neutralino3amplitudesRsbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), md(2,2), ms, g, gp, mixNeut, 3, 'd', 'R');
     neutralino3amplitudesLbars = neutralino3amplitudesLsbar;
     neutralino3amplitudesRbars = neutralino3amplitudesRsbar;
     neutralino3amplitudeeLebar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), me(1,1), mel, g, gp, mixNeut, 3, 'l', 'L');
     neutralino3amplitudeeRebar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), me(2,1), mel, g, gp, mixNeut, 3, 'l', 'R');
     neutralino3amplitudeeLbare = neutralino3amplitudeeLebar;
     neutralino3amplitudeeRbare = neutralino3amplitudeeRebar;
     neutralino3amplitudemuLmubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), me(1,2), mmu, g, gp, mixNeut, 3, 'l', 'L');
     neutralino3amplitudemuRmubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(3), me(2,2), mmu, g, gp, mixNeut, 3, 'l', 'R');
     neutralino3amplitudemuLbarmu = neutralino3amplitudemuLmubar;
     neutralino3amplitudemuRbarmu = neutralino3amplitudemuRmubar;
     
     neutralino3amplitudetopstop1bar = neutralinoamplitudestoptopNMSSM (mneut(3), mu(1,3), mt, g, gp, thetat, beta, runmw, mixNeut, runmt, 3, 1);
     neutralino3amplitudetopstop2bar = neutralinoamplitudestoptopNMSSM (mneut(3), mu(2,3), mt, g, gp, thetat, beta, runmw, mixNeut, runmt, 3, 2);
     neutralino3amplitudetopbarstop1 = neutralino3amplitudetopstop1bar;
     neutralino3amplitudetopbarstop2 = neutralino3amplitudetopstop2bar;
     neutralino3amplitudebottomsbottom1bar = neutralinoamplitudesbottombottomNMSSM (mneut(3), md(1,3), mb, g, gp, thetab, beta, runmw, mixNeut, runmb, 3, 1);
     neutralino3amplitudebottomsbottom2bar = neutralinoamplitudesbottombottomNMSSM (mneut(3), md(2,3), mb, g, gp, thetab, beta, runmw, mixNeut, runmb, 3, 2);
     neutralino3amplitudebottombarsbottom1 = neutralino3amplitudebottomsbottom1bar;
     neutralino3amplitudebottombarsbottom2 = neutralino3amplitudebottomsbottom2bar;
     
     neutralino3amplitudetaustau1bar = neutralinoamplitudestautauNMSSM (mneut(3), me(1,3), mtau, g, gp, thetatau, beta, runmw, mixNeut, runmtau, 3, 1);
     neutralino3amplitudetaustau2bar = neutralinoamplitudestautauNMSSM (mneut(3), me(2,3), mtau, g, gp, thetatau, beta, runmw, mixNeut, runmtau, 3, 2);
     neutralino3amplitudetaubarstau1 = neutralino3amplitudetaustau1bar;
     neutralino3amplitudetaubarstau2 = neutralino3amplitudetaustau2bar;
     
     neutralino3amplitudesnuenuebar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(3), msnu(1), 0, g, gp, mixNeut, 3);
     neutralino3amplitudesnuebarnue = neutralino3amplitudesnuenuebar;
     neutralino3amplitudesnumunumubar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(3), msnu(2), 0, g, gp, mixNeut, 3);
     neutralino3amplitudesnumubarnumu = neutralino3amplitudesnumunumubar;
     neutralino3amplitudenutausnutaubar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(3), msnu(3), 0, g, gp, mixNeut, 3);
     neutralino3amplitudenutaubarsnutau = neutralino3amplitudenutausnutaubar;
   }
   
   ParticleNeutralino3.Array_Decays[0][0] = -PDGup; ParticleNeutralino3.Array_Decays[0][1] = PDGsupL; ParticleNeutralino3.Array_Decays[0][2] = neutralino3amplitudeuLubar; ParticleNeutralino3.Array_Decays[0][3] = 2; ParticleNeutralino3.Array_Comments[0] = "# ~chi_30 -> ub ~u_L";
   ParticleNeutralino3.Array_Decays[1][0] = -PDGup; ParticleNeutralino3.Array_Decays[1][1] = PDGsupR; ParticleNeutralino3.Array_Decays[1][2] = neutralino3amplitudeuRubar; ParticleNeutralino3.Array_Decays[1][3] = 2; ParticleNeutralino3.Array_Comments[1] = "# ~chi_30 -> ub ~u_R";
   ParticleNeutralino3.Array_Decays[2][0] = PDGup; ParticleNeutralino3.Array_Decays[2][1] = -PDGsupL; ParticleNeutralino3.Array_Decays[2][2] = neutralino3amplitudeuLbaru; ParticleNeutralino3.Array_Decays[2][3] = 2; ParticleNeutralino3.Array_Comments[2] = "# ~chi_30 -> u ~u_L*";
   ParticleNeutralino3.Array_Decays[3][0] = PDGup; ParticleNeutralino3.Array_Decays[3][1] = -PDGsupR; ParticleNeutralino3.Array_Decays[3][2] = neutralino3amplitudeuRbaru; ParticleNeutralino3.Array_Decays[3][3] = 2; ParticleNeutralino3.Array_Comments[3] = "# ~chi_30 -> u ~u_R*";
   ParticleNeutralino3.Array_Decays[4][0] = -PDGdown; ParticleNeutralino3.Array_Decays[4][1] = PDGsdownL; ParticleNeutralino3.Array_Decays[4][2] = neutralino3amplitudedLdbar; ParticleNeutralino3.Array_Decays[4][3] = 2; ParticleNeutralino3.Array_Comments[4] = "# ~chi_30 -> db ~d_L";
   ParticleNeutralino3.Array_Decays[5][0] = -PDGdown; ParticleNeutralino3.Array_Decays[5][1] = PDGsdownR; ParticleNeutralino3.Array_Decays[5][2] = neutralino3amplitudedRdbar; ParticleNeutralino3.Array_Decays[5][3] = 2; ParticleNeutralino3.Array_Comments[5] = "# ~chi_30 -> db ~d_R";
   ParticleNeutralino3.Array_Decays[6][0] = PDGdown; ParticleNeutralino3.Array_Decays[6][1] = -PDGsdownL; ParticleNeutralino3.Array_Decays[6][2] = neutralino3amplitudedLbard; ParticleNeutralino3.Array_Decays[6][3] = 2; ParticleNeutralino3.Array_Comments[6] = "# ~chi_30 -> d ~d_L*";
   ParticleNeutralino3.Array_Decays[7][0] = PDGdown; ParticleNeutralino3.Array_Decays[7][1] = -PDGsdownR; ParticleNeutralino3.Array_Decays[7][2] = neutralino3amplitudedRbard; ParticleNeutralino3.Array_Decays[7][3] = 2; ParticleNeutralino3.Array_Comments[7] = "# ~chi_30 -> d ~d_R*";
   ParticleNeutralino3.Array_Decays[8][0] = -PDGcharm; ParticleNeutralino3.Array_Decays[8][1] = PDGscharmL; ParticleNeutralino3.Array_Decays[8][2] = neutralino3amplitudecLcbar; ParticleNeutralino3.Array_Decays[8][3] = 2; ParticleNeutralino3.Array_Comments[8] = "# ~chi_30 -> cb ~c_L";
   ParticleNeutralino3.Array_Decays[9][0] = -PDGcharm; ParticleNeutralino3.Array_Decays[9][1] = PDGscharmR; ParticleNeutralino3.Array_Decays[9][2] = neutralino3amplitudecRcbar; ParticleNeutralino3.Array_Decays[9][3] = 2; ParticleNeutralino3.Array_Comments[9] = "# ~chi_30 -> cb ~c_R";
   ParticleNeutralino3.Array_Decays[10][0] = PDGcharm; ParticleNeutralino3.Array_Decays[10][1] = -PDGscharmL; ParticleNeutralino3.Array_Decays[10][2] = neutralino3amplitudecLbarc; ParticleNeutralino3.Array_Decays[10][3] = 2; ParticleNeutralino3.Array_Comments[10] = "# ~chi_30 -> c ~c_L*"; 
   ParticleNeutralino3.Array_Decays[11][0] = PDGcharm; ParticleNeutralino3.Array_Decays[11][1] = -PDGscharmR; ParticleNeutralino3.Array_Decays[11][2] = neutralino3amplitudecRbarc; ParticleNeutralino3.Array_Decays[11][3] = 2; ParticleNeutralino3.Array_Comments[11] = "# ~chi_30 -> c ~c_R*";
   ParticleNeutralino3.Array_Decays[12][0] = -PDGstrange; ParticleNeutralino3.Array_Decays[12][1] = PDGsstrangeL; ParticleNeutralino3.Array_Decays[12][2] = neutralino3amplitudesLsbar; ParticleNeutralino3.Array_Decays[12][3] = 2; ParticleNeutralino3.Array_Comments[12] = "# ~chi_30 -> sb ~s_L";
   ParticleNeutralino3.Array_Decays[13][0] = -PDGstrange; ParticleNeutralino3.Array_Decays[13][1] = PDGsstrangeR; ParticleNeutralino3.Array_Decays[13][2] = neutralino3amplitudesRsbar; ParticleNeutralino3.Array_Decays[13][3] = 2; ParticleNeutralino3.Array_Comments[13] = "# ~chi_30 -> sb ~s_R";
   ParticleNeutralino3.Array_Decays[14][0] = PDGstrange; ParticleNeutralino3.Array_Decays[14][1] = -PDGsstrangeL; ParticleNeutralino3.Array_Decays[14][2] = neutralino3amplitudesLbars; ParticleNeutralino3.Array_Decays[14][3] = 2; ParticleNeutralino3.Array_Comments[14] = "# ~chi_30 -> s ~s_L*";
   ParticleNeutralino3.Array_Decays[15][0] = PDGstrange; ParticleNeutralino3.Array_Decays[15][1] = -PDGsstrangeR; ParticleNeutralino3.Array_Decays[15][2] = neutralino3amplitudesRbars; ParticleNeutralino3.Array_Decays[15][3] = 2; ParticleNeutralino3.Array_Comments[15] = "# ~chi_30 -> s ~s_R*";
   ParticleNeutralino3.Array_Decays[16][0] = -PDGelectron; ParticleNeutralino3.Array_Decays[16][1] = PDGselectronL; ParticleNeutralino3.Array_Decays[16][2] = neutralino3amplitudeeLebar; ParticleNeutralino3.Array_Decays[16][3] = 2; ParticleNeutralino3.Array_Comments[16] = "# ~chi_30 -> e+ ~e_L-";
   ParticleNeutralino3.Array_Decays[17][0] = -PDGelectron; ParticleNeutralino3.Array_Decays[17][1] = PDGselectronR; ParticleNeutralino3.Array_Decays[17][2] = neutralino3amplitudeeRebar; ParticleNeutralino3.Array_Decays[17][3] = 2; ParticleNeutralino3.Array_Comments[17] = "# ~chi_30 -> e+ ~e_R-";
   ParticleNeutralino3.Array_Decays[18][0] = PDGelectron; ParticleNeutralino3.Array_Decays[18][1] = -PDGselectronL; ParticleNeutralino3.Array_Decays[18][2] = neutralino3amplitudeeLbare; ParticleNeutralino3.Array_Decays[18][3] = 2; ParticleNeutralino3.Array_Comments[18] = "# ~chi_30 -> e- ~e_L+";
   ParticleNeutralino3.Array_Decays[19][0] = PDGelectron; ParticleNeutralino3.Array_Decays[19][1] = -PDGselectronR; ParticleNeutralino3.Array_Decays[19][2] = neutralino3amplitudeeRbare; ParticleNeutralino3.Array_Decays[19][3] = 2; ParticleNeutralino3.Array_Comments[19] = "# ~chi_30 -> e- ~e_R+";   
   ParticleNeutralino3.Array_Decays[20][0] = -PDGmuon; ParticleNeutralino3.Array_Decays[20][1] = PDGsmuonL; ParticleNeutralino3.Array_Decays[20][2] = neutralino3amplitudemuLmubar; ParticleNeutralino3.Array_Decays[20][3] = 2; ParticleNeutralino3.Array_Comments[20] = "# ~chi_30 -> mu+ ~mu_L-";
   ParticleNeutralino3.Array_Decays[21][0] = -PDGmuon; ParticleNeutralino3.Array_Decays[21][1] = PDGsmuonR; ParticleNeutralino3.Array_Decays[21][2] = neutralino3amplitudemuRmubar; ParticleNeutralino3.Array_Decays[21][3] = 2; ParticleNeutralino3.Array_Comments[21] = "# ~chi_30 -> mu+ ~mu_R-";
   ParticleNeutralino3.Array_Decays[22][0] = PDGmuon; ParticleNeutralino3.Array_Decays[22][1] = -PDGsmuonL; ParticleNeutralino3.Array_Decays[22][2] = neutralino3amplitudemuLbarmu; ParticleNeutralino3.Array_Decays[22][3] = 2; ParticleNeutralino3.Array_Comments[22] = "# ~chi_30 -> mu- ~mu_L+";
   ParticleNeutralino3.Array_Decays[23][0] = PDGmuon; ParticleNeutralino3.Array_Decays[23][1] = -PDGsmuonR; ParticleNeutralino3.Array_Decays[23][2] = neutralino3amplitudemuRbarmu; ParticleNeutralino3.Array_Decays[23][3] = 2; ParticleNeutralino3.Array_Comments[23] = "# ~chi_30 -> mu- ~mu_R+";
   ParticleNeutralino3.Array_Decays[24][0] = PDGnuelectron; ParticleNeutralino3.Array_Decays[24][1] = -PDGnuselectronL; ParticleNeutralino3.Array_Decays[24][2] = neutralino3amplitudesnuebarnue; ParticleNeutralino3.Array_Decays[24][3] = 2; ParticleNeutralino3.Array_Comments[24] = "# ~chi_30 -> nu_e ~nu_eL*";
   ParticleNeutralino3.Array_Decays[25][0] = -PDGnuelectron; ParticleNeutralino3.Array_Decays[25][1] = PDGnuselectronL; ParticleNeutralino3.Array_Decays[25][2] = neutralino3amplitudesnuenuebar; ParticleNeutralino3.Array_Decays[25][3] = 2; ParticleNeutralino3.Array_Comments[25] = "# ~chi_30 -> nu_eb ~nu_eL";
   ParticleNeutralino3.Array_Decays[26][0] = PDGnumuon; ParticleNeutralino3.Array_Decays[26][1] = -PDGnusmuonL; ParticleNeutralino3.Array_Decays[26][2] = neutralino3amplitudesnumubarnumu; ParticleNeutralino3.Array_Decays[26][3] = 2; ParticleNeutralino3.Array_Comments[26] = "# ~chi_30 -> nu_mu ~nu_muL*";
   ParticleNeutralino3.Array_Decays[27][0] = -PDGnumuon; ParticleNeutralino3.Array_Decays[27][1] = PDGnusmuonL; ParticleNeutralino3.Array_Decays[27][2] = neutralino3amplitudesnumunumubar; ParticleNeutralino3.Array_Decays[27][3] = 2; ParticleNeutralino3.Array_Comments[27] = "# ~chi_30 -> nu_mub ~nu_muL";
   ParticleNeutralino3.Array_Decays[28][0] = PDGtop; ParticleNeutralino3.Array_Decays[28][1] = -PDGstop1; ParticleNeutralino3.Array_Decays[28][2] = neutralino3amplitudetopstop1bar; ParticleNeutralino3.Array_Decays[28][3] = 2; ParticleNeutralino3.Array_Comments[28] = "# ~chi_30 -> t ~t_1*";
   ParticleNeutralino3.Array_Decays[29][0] = PDGtop; ParticleNeutralino3.Array_Decays[29][1] = -PDGstop2; ParticleNeutralino3.Array_Decays[29][2] = neutralino3amplitudetopstop2bar; ParticleNeutralino3.Array_Decays[29][3] = 2; ParticleNeutralino3.Array_Comments[29] = "# ~chi_30 -> t ~t_2*";
   ParticleNeutralino3.Array_Decays[30][0] = -PDGtop; ParticleNeutralino3.Array_Decays[30][1] = PDGstop1; ParticleNeutralino3.Array_Decays[30][2] = neutralino3amplitudetopbarstop1; ParticleNeutralino3.Array_Decays[30][3] = 2; ParticleNeutralino3.Array_Comments[30] = "# ~chi_30 -> tb ~t_1";
   ParticleNeutralino3.Array_Decays[31][0] = -PDGtop; ParticleNeutralino3.Array_Decays[31][1] = PDGstop2; ParticleNeutralino3.Array_Decays[31][2] = neutralino3amplitudetopbarstop2; ParticleNeutralino3.Array_Decays[31][3] = 2; ParticleNeutralino3.Array_Comments[31] = "# ~chi_30 -> tb ~t_2";
   ParticleNeutralino3.Array_Decays[32][0] = PDGbottom; ParticleNeutralino3.Array_Decays[32][1] = -PDGsbottom1; ParticleNeutralino3.Array_Decays[32][2] = neutralino3amplitudebottomsbottom1bar; ParticleNeutralino3.Array_Decays[32][3] = 2; ParticleNeutralino3.Array_Comments[32] = "# ~chi_30 -> b ~b_1*";
   ParticleNeutralino3.Array_Decays[33][0] = PDGbottom; ParticleNeutralino3.Array_Decays[33][1] = -PDGsbottom2; ParticleNeutralino3.Array_Decays[33][2] = neutralino3amplitudebottomsbottom2bar; ParticleNeutralino3.Array_Decays[33][3] = 2; ParticleNeutralino3.Array_Comments[33] = "# ~chi_30 -> b ~b_2*";
   ParticleNeutralino3.Array_Decays[34][0] = -PDGbottom; ParticleNeutralino3.Array_Decays[34][1] = PDGsbottom1; ParticleNeutralino3.Array_Decays[34][2] = neutralino3amplitudebottombarsbottom1; ParticleNeutralino3.Array_Decays[34][3] = 2; ParticleNeutralino3.Array_Comments[34] = "# ~chi_30 -> bb ~b_1";
   ParticleNeutralino3.Array_Decays[35][0] = -PDGbottom; ParticleNeutralino3.Array_Decays[35][1] = PDGsbottom2; ParticleNeutralino3.Array_Decays[35][2] = neutralino3amplitudebottombarsbottom2; ParticleNeutralino3.Array_Decays[35][3] = 2; ParticleNeutralino3.Array_Comments[35] = "# ~chi_30 -> bb ~b_2";
   ParticleNeutralino3.Array_Decays[36][0] = -PDGstau1; ParticleNeutralino3.Array_Decays[36][1] = PDGtau; ParticleNeutralino3.Array_Decays[36][2] = neutralino3amplitudetaustau1bar; ParticleNeutralino3.Array_Decays[36][3] = 2; ParticleNeutralino3.Array_Comments[36] = "# ~chi_30 -> tau- ~tau_1+";
   ParticleNeutralino3.Array_Decays[37][0] = -PDGstau2; ParticleNeutralino3.Array_Decays[37][1] = PDGtau; ParticleNeutralino3.Array_Decays[37][2] = neutralino3amplitudetaustau2bar; ParticleNeutralino3.Array_Decays[37][3] = 2; ParticleNeutralino3.Array_Comments[37] = "# ~chi_30 -> tau- ~tau_2+";
   ParticleNeutralino3.Array_Decays[38][0] = PDGstau1; ParticleNeutralino3.Array_Decays[38][1] = -PDGtau; ParticleNeutralino3.Array_Decays[38][2] = neutralino3amplitudetaubarstau1; ParticleNeutralino3.Array_Decays[38][3] = 2; ParticleNeutralino3.Array_Comments[38] = "# ~chi_30 -> tau+ ~tau_1-";
   ParticleNeutralino3.Array_Decays[39][0] = PDGstau2; ParticleNeutralino3.Array_Decays[39][1] = -PDGtau; ParticleNeutralino3.Array_Decays[39][2] = neutralino3amplitudetaubarstau2; ParticleNeutralino3.Array_Decays[39][3] = 2; ParticleNeutralino3.Array_Comments[39] = "# ~chi_30 -> tau+ ~tau_2-"; 
   ParticleNeutralino3.Array_Decays[40][0] = PDGnutau; ParticleNeutralino3.Array_Decays[40][1] = -PDGnustauL; ParticleNeutralino3.Array_Decays[40][2] = neutralino3amplitudenutausnutaubar; ParticleNeutralino3.Array_Decays[40][3] = 2; ParticleNeutralino3.Array_Comments[40] = "# ~chi_30 -> nu_tau ~nu_tauL*";
   ParticleNeutralino3.Array_Decays[41][0] = -PDGnutau; ParticleNeutralino3.Array_Decays[41][1] = PDGnustauL; ParticleNeutralino3.Array_Decays[41][2] = neutralino3amplitudenutaubarsnutau; ParticleNeutralino3.Array_Decays[41][3] = 2; ParticleNeutralino3.Array_Comments[41] = "# ~chi_30 -> nu_taub ~nu_tauL";
   ParticleNeutralino3.Array_Decays[42][0] = PDGWplus; ParticleNeutralino3.Array_Decays[42][1] = -PDGchargino1; ParticleNeutralino3.Array_Decays[42][2] = neutralino3amplitudeWbosonpluscharginoW1; ParticleNeutralino3.Array_Decays[42][3] = 2; ParticleNeutralino3.Array_Comments[42] = "# ~chi_30 -> W+ ~chi_1-";
   ParticleNeutralino3.Array_Decays[43][0] = PDGWplus; ParticleNeutralino3.Array_Decays[43][1] = -PDGchargino2; ParticleNeutralino3.Array_Decays[43][2] = neutralino3amplitudeWbosonpluscharginoW2; ParticleNeutralino3.Array_Decays[43][3] = 2; ParticleNeutralino3.Array_Comments[43] = "# ~chi_30 -> W+ ~chi_2-";
   ParticleNeutralino3.Array_Decays[44][0] = -PDGWplus; ParticleNeutralino3.Array_Decays[44][1] = PDGchargino1; ParticleNeutralino3.Array_Decays[44][2] = neutralino3amplitudeWbosonminuscharginoW1; ParticleNeutralino3.Array_Decays[44][3] = 2; ParticleNeutralino3.Array_Comments[44] = "# ~chi_30 -> W- ~chi_1+";
   ParticleNeutralino3.Array_Decays[45][0] = -PDGWplus; ParticleNeutralino3.Array_Decays[45][1] = PDGchargino2; ParticleNeutralino3.Array_Decays[45][2] = neutralino3amplitudeWbosonminuscharginoW2; ParticleNeutralino3.Array_Decays[45][3] = 2; ParticleNeutralino3.Array_Comments[45] = "# ~chi_30 -> W- ~chi_2+";
   ParticleNeutralino3.Array_Decays[46][0] = PDGHplus; ParticleNeutralino3.Array_Decays[46][1] = -PDGchargino1; ParticleNeutralino3.Array_Decays[46][2] = neutralino3amplitudeHpluscharginoW1; ParticleNeutralino3.Array_Decays[46][3] = 2; ParticleNeutralino3.Array_Comments[46] = "# ~chi_30 -> H+ ~chi_1-";
   ParticleNeutralino3.Array_Decays[47][0] = PDGHplus; ParticleNeutralino3.Array_Decays[47][1] = -PDGchargino2; ParticleNeutralino3.Array_Decays[47][2] = neutralino3amplitudeHpluscharginoW2; ParticleNeutralino3.Array_Decays[47][3] = 2; ParticleNeutralino3.Array_Comments[47] = "# ~chi_30 -> H+ ~chi_2-";
   ParticleNeutralino3.Array_Decays[48][0] = -PDGHplus; ParticleNeutralino3.Array_Decays[48][1] = PDGchargino1; ParticleNeutralino3.Array_Decays[48][2] = neutralino3amplitudeHminuscharginoW1; ParticleNeutralino3.Array_Decays[48][3] = 2; ParticleNeutralino3.Array_Comments[48] = "# ~chi_30 -> H- ~chi_1+";
   ParticleNeutralino3.Array_Decays[49][0] = -PDGHplus; ParticleNeutralino3.Array_Decays[49][1] = PDGchargino2; ParticleNeutralino3.Array_Decays[49][2] = neutralino3amplitudeHminuscharginoW2; ParticleNeutralino3.Array_Decays[49][3] = 2; ParticleNeutralino3.Array_Comments[49] = "# ~chi_30 -> H- ~chi_2+";
   ParticleNeutralino3.Array_Decays[50][0] = PDGZboson; ParticleNeutralino3.Array_Decays[50][1] = PDGneutralino1; ParticleNeutralino3.Array_Decays[50][2] = neutralino3amplitudeZbosonneutralino1; ParticleNeutralino3.Array_Decays[50][3] = 2; ParticleNeutralino3.Array_Comments[50] = "# ~chi_30 -> Z ~chi_10";
   ParticleNeutralino3.Array_Decays[51][0] = PDGZboson; ParticleNeutralino3.Array_Decays[51][1] = PDGneutralino2; ParticleNeutralino3.Array_Decays[51][2] = neutralino3amplitudeZbosonneutralino2; ParticleNeutralino3.Array_Decays[51][3] = 2; ParticleNeutralino3.Array_Comments[51] = "# ~chi_30 -> Z ~chi_20";
   ParticleNeutralino3.Array_Decays[52][0] = PDGZboson; ParticleNeutralino3.Array_Decays[52][1] = PDGneutralino4; ParticleNeutralino3.Array_Decays[52][2] = neutralino3amplitudeZbosonneutralino4; ParticleNeutralino3.Array_Decays[52][3] = 2; ParticleNeutralino3.Array_Comments[52] = "# ~chi_30 -> Z ~chi_40";
   ParticleNeutralino3.Array_Decays[53][0] = PDGh0; ParticleNeutralino3.Array_Decays[53][1] = PDGneutralino1; ParticleNeutralino3.Array_Decays[53][2] = neutralino3amplitudehneutralino1; ParticleNeutralino3.Array_Decays[53][3] = 2; ParticleNeutralino3.Array_Comments[53] = "# ~chi_30 -> h ~chi_10";
   ParticleNeutralino3.Array_Decays[54][0] = PDGh0; ParticleNeutralino3.Array_Decays[54][1] = PDGneutralino2; ParticleNeutralino3.Array_Decays[54][2] = neutralino3amplitudehneutralino2; ParticleNeutralino3.Array_Decays[54][3] = 2; ParticleNeutralino3.Array_Comments[54] = "# ~chi_30 -> h ~chi_20";
   ParticleNeutralino3.Array_Decays[55][0] = PDGh0; ParticleNeutralino3.Array_Decays[55][1] = PDGneutralino4; ParticleNeutralino3.Array_Decays[55][2] = neutralino3amplitudehneutralino4; ParticleNeutralino3.Array_Decays[55][3] = 2; ParticleNeutralino3.Array_Comments[55] = "# ~chi_30 -> h ~chi_40";
   ParticleNeutralino3.Array_Decays[56][0] = PDGH0; ParticleNeutralino3.Array_Decays[56][1] = PDGneutralino1; ParticleNeutralino3.Array_Decays[56][2] = neutralino3amplitudeHneutralino1; ParticleNeutralino3.Array_Decays[56][3] = 2; ParticleNeutralino3.Array_Comments[56] = "# ~chi_30 -> H ~chi_10";
   ParticleNeutralino3.Array_Decays[57][0] = PDGH0; ParticleNeutralino3.Array_Decays[57][1] = PDGneutralino2; ParticleNeutralino3.Array_Decays[57][2] = neutralino3amplitudeHneutralino2; ParticleNeutralino3.Array_Decays[57][3] = 2; ParticleNeutralino3.Array_Comments[57] = "# ~chi_30 -> H ~chi_20";
   ParticleNeutralino3.Array_Decays[58][0] = PDGH0; ParticleNeutralino3.Array_Decays[58][1] = PDGneutralino4; ParticleNeutralino3.Array_Decays[58][2] = neutralino3amplitudeHneutralino4; ParticleNeutralino3.Array_Decays[58][3] = 2; ParticleNeutralino3.Array_Comments[58] = "# ~chi_30 -> H ~chi_40";
   ParticleNeutralino3.Array_Decays[59][0] = PDGA0; ParticleNeutralino3.Array_Decays[59][1] = PDGneutralino1; ParticleNeutralino3.Array_Decays[59][2] = neutralino3amplitudeAneutralino1; ParticleNeutralino3.Array_Decays[59][3] = 2; ParticleNeutralino3.Array_Comments[59] = "# ~chi_30 -> A ~chi_10";
   ParticleNeutralino3.Array_Decays[60][0] = PDGA0; ParticleNeutralino3.Array_Decays[60][1] = PDGneutralino2; ParticleNeutralino3.Array_Decays[60][2] = neutralino3amplitudeAneutralino2; ParticleNeutralino3.Array_Decays[60][3] = 2; ParticleNeutralino3.Array_Comments[60] = "# ~chi_30 -> A ~chi_20";
   ParticleNeutralino3.Array_Decays[61][0] = PDGA0; ParticleNeutralino3.Array_Decays[61][1] = PDGneutralino4; ParticleNeutralino3.Array_Decays[61][2] = neutralino3amplitudeAneutralino4; ParticleNeutralino3.Array_Decays[61][3] = 2; ParticleNeutralino3.Array_Comments[61] = "# ~chi_30 -> A ~chi_40";
   
   ParticleNeutralino3.Array_Decays[62][0] = PDGH3; ParticleNeutralino3.Array_Decays[62][1] = PDGneutralino1; ParticleNeutralino3.Array_Decays[62][2] = neutralino3amplitudeH3neutralino1; ParticleNeutralino3.Array_Decays[62][3] = 2; ParticleNeutralino3.Array_Comments[62] = "# ~chi_30 -> H3 ~chi_10";
   ParticleNeutralino3.Array_Decays[63][0] = PDGH3; ParticleNeutralino3.Array_Decays[63][1] = PDGneutralino2; ParticleNeutralino3.Array_Decays[63][2] = neutralino3amplitudeH3neutralino2; ParticleNeutralino3.Array_Decays[63][3] = 2; ParticleNeutralino3.Array_Comments[63] = "# ~chi_30 -> H3 ~chi_20";
   ParticleNeutralino3.Array_Decays[64][0] = PDGA2; ParticleNeutralino3.Array_Decays[64][1] = PDGneutralino1; ParticleNeutralino3.Array_Decays[64][2] = neutralino3amplitudeA2neutralino1; ParticleNeutralino3.Array_Decays[64][3] = 2; ParticleNeutralino3.Array_Comments[64] = "# ~chi_30 -> A2 ~chi_10";
   ParticleNeutralino3.Array_Decays[65][0] = PDGA2; ParticleNeutralino3.Array_Decays[65][1] = PDGneutralino2; ParticleNeutralino3.Array_Decays[65][2] = neutralino3amplitudeA2neutralino2; ParticleNeutralino3.Array_Decays[65][3] = 2; ParticleNeutralino3.Array_Comments[65] = "# ~chi_30 -> A2 ~chi_20";
   
   ParticleNeutralino3.Array_Decays[66][0] = PDGphoton; ParticleNeutralino3.Array_Decays[66][1] = PDGgravitino; ParticleNeutralino3.Array_Decays[66][2] = neutralino3amplitudephotongravitino; ParticleNeutralino3.Array_Decays[66][3] = 2; ParticleNeutralino3.Array_Comments[66] = "# ~chi_30 -> gamma ~G";
   ParticleNeutralino3.Array_Decays[67][0] = PDGZboson; ParticleNeutralino3.Array_Decays[67][1] = PDGgravitino; ParticleNeutralino3.Array_Decays[67][2] = neutralino3amplitudeZgravitino; ParticleNeutralino3.Array_Decays[67][3] = 2; ParticleNeutralino3.Array_Comments[67] = "# ~chi_30 -> Z ~G";
   ParticleNeutralino3.Array_Decays[68][0] = PDGh0; ParticleNeutralino3.Array_Decays[68][1] = PDGgravitino; ParticleNeutralino3.Array_Decays[68][2] = neutralino3amplitudehgravitino; ParticleNeutralino3.Array_Decays[68][3] = 2; ParticleNeutralino3.Array_Comments[68] = "# ~chi_30 -> h ~G";
   ParticleNeutralino3.Array_Decays[69][0] = PDGH0; ParticleNeutralino3.Array_Decays[69][1] = PDGgravitino; ParticleNeutralino3.Array_Decays[69][2] = neutralino3amplitudeHgravitino; ParticleNeutralino3.Array_Decays[69][3] = 2; ParticleNeutralino3.Array_Comments[69] = "# ~chi_30 -> H ~G";
   ParticleNeutralino3.Array_Decays[70][0] = PDGA0; ParticleNeutralino3.Array_Decays[70][1] = PDGgravitino; ParticleNeutralino3.Array_Decays[70][2] = neutralino3amplitudeAgravitino; ParticleNeutralino1.Array_Decays[70][3] = 2; ParticleNeutralino3.Array_Comments[70] = "# ~chi_30 -> A ~G";
   
   ParticleNeutralino3.Array_Decays[71][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[71][1] = PDGup; ParticleNeutralino3.Array_Decays[71][4] = -PDGup; ParticleNeutralino3.Array_Decays[71][2] = neutralino3amplitudeneut1uubar; ParticleNeutralino3.Array_Decays[71][3] = 3; ParticleNeutralino3.Array_Comments[71] = "# ~chi_30 -> ~chi_10 u ubar";
   ParticleNeutralino3.Array_Decays[72][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[72][1] = PDGdown; ParticleNeutralino3.Array_Decays[72][4] = -PDGdown; ParticleNeutralino3.Array_Decays[72][2] = neutralino3amplitudeneut1ddbar; ParticleNeutralino3.Array_Decays[72][3] = 3; ParticleNeutralino3.Array_Comments[72] = "# ~chi_30 -> ~chi_10 d dbar";
   ParticleNeutralino3.Array_Decays[73][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[73][1] = PDGcharm; ParticleNeutralino3.Array_Decays[73][4] = -PDGcharm; ParticleNeutralino3.Array_Decays[73][2] = neutralino3amplitudeneut1ccbar; ParticleNeutralino3.Array_Decays[73][3] = 3; ParticleNeutralino3.Array_Comments[73] = "# ~chi_30 -> ~chi_10 c cbar";
   ParticleNeutralino3.Array_Decays[74][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[74][1] = PDGstrange; ParticleNeutralino3.Array_Decays[74][4] = -PDGstrange; ParticleNeutralino3.Array_Decays[74][2] = neutralino3amplitudeneut1ssbar; ParticleNeutralino3.Array_Decays[74][3] = 3; ParticleNeutralino3.Array_Comments[74] = "# ~chi_30 -> ~chi_10 s sbar";
   ParticleNeutralino3.Array_Decays[75][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[75][1] = PDGtop; ParticleNeutralino3.Array_Decays[75][4] = -PDGtop; ParticleNeutralino3.Array_Decays[75][2] = neutralino3amplitudeneut1ttbar; ParticleNeutralino3.Array_Decays[75][3] = 3; ParticleNeutralino3.Array_Comments[75] = "# ~chi_30 -> ~chi_10 t tbar";
   ParticleNeutralino3.Array_Decays[76][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[76][1] = PDGbottom; ParticleNeutralino3.Array_Decays[76][4] = -PDGbottom; ParticleNeutralino3.Array_Decays[76][2] = neutralino3amplitudeneut1bbbar; ParticleNeutralino3.Array_Decays[76][3] = 3; ParticleNeutralino3.Array_Comments[76] = "# ~chi_30 -> ~chi_10 b bbar";
   ParticleNeutralino3.Array_Decays[77][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[77][1] = PDGelectron; ParticleNeutralino3.Array_Decays[77][4] = -PDGelectron; ParticleNeutralino3.Array_Decays[77][2] = neutralino3amplitudeneut1eebar; ParticleNeutralino3.Array_Decays[77][3] = 3; ParticleNeutralino3.Array_Comments[77] = "# ~chi_30 -> ~chi_10 e- e+";
   ParticleNeutralino3.Array_Decays[78][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[78][1] = PDGmuon; ParticleNeutralino3.Array_Decays[78][4] = -PDGmuon; ParticleNeutralino3.Array_Decays[78][2] = neutralino3amplitudeneut1mumubar; ParticleNeutralino3.Array_Decays[78][3] = 3; ParticleNeutralino3.Array_Comments[78] = "# ~chi_30 -> ~chi_10 mu- mu+";
   ParticleNeutralino3.Array_Decays[79][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[79][1] = PDGtau; ParticleNeutralino3.Array_Decays[79][4] = -PDGtau; ParticleNeutralino3.Array_Decays[79][2] = neutralino3amplitudeneut1tautaubar; ParticleNeutralino3.Array_Decays[79][3] = 3; ParticleNeutralino3.Array_Comments[79] = "# ~chi_30 -> ~chi_10 tau- tau+";
   ParticleNeutralino3.Array_Decays[80][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[80][1] = PDGnuelectron; ParticleNeutralino3.Array_Decays[80][4] = -PDGnuelectron; ParticleNeutralino3.Array_Decays[80][2] = neutralino3amplitudeneut1nuenuebar; ParticleNeutralino3.Array_Decays[80][3] = 3; ParticleNeutralino3.Array_Comments[80] = "# ~chi_30 -> ~chi_10 nue nuebar";
   ParticleNeutralino3.Array_Decays[81][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[81][1] = PDGnumuon; ParticleNeutralino3.Array_Decays[81][4] = -PDGnumuon; ParticleNeutralino3.Array_Decays[81][2] = neutralino3amplitudeneut1numunumubar; ParticleNeutralino3.Array_Decays[81][3] = 3; ParticleNeutralino3.Array_Comments[81] = "# ~chi_30 -> ~chi_10 numu numubar";
   ParticleNeutralino3.Array_Decays[82][0] = PDGneutralino1; ParticleNeutralino3.Array_Decays[82][1] = PDGnutau; ParticleNeutralino3.Array_Decays[82][4] = -PDGnutau; ParticleNeutralino3.Array_Decays[82][2] = neutralino3amplitudeneut1nutaunutaubar; ParticleNeutralino3.Array_Decays[82][3] = 3; ParticleNeutralino3.Array_Comments[82] = "# ~chi_30 -> ~chi_10 nutau nutaubar";
   ParticleNeutralino3.Array_Decays[83][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[83][1] = PDGup; ParticleNeutralino3.Array_Decays[83][4] = -PDGup; ParticleNeutralino3.Array_Decays[83][2] = neutralino3amplitudeneut2uubar; ParticleNeutralino3.Array_Decays[83][3] = 3; ParticleNeutralino3.Array_Comments[83] = "# ~chi_30 -> ~chi_20 u ubar";
   ParticleNeutralino3.Array_Decays[84][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[84][1] = PDGdown; ParticleNeutralino3.Array_Decays[84][4] = -PDGdown; ParticleNeutralino3.Array_Decays[84][2] = neutralino3amplitudeneut2ddbar; ParticleNeutralino3.Array_Decays[84][3] = 3; ParticleNeutralino3.Array_Comments[84] = "# ~chi_30 -> ~chi_20 d dbar";
   ParticleNeutralino3.Array_Decays[85][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[85][1] = PDGcharm; ParticleNeutralino3.Array_Decays[85][4] = -PDGcharm; ParticleNeutralino3.Array_Decays[85][2] = neutralino3amplitudeneut2ccbar; ParticleNeutralino3.Array_Decays[85][3] = 3; ParticleNeutralino3.Array_Comments[85] = "# ~chi_30 -> ~chi_20 c cbar";
   ParticleNeutralino3.Array_Decays[86][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[86][1] = PDGstrange; ParticleNeutralino3.Array_Decays[86][4] = -PDGstrange; ParticleNeutralino3.Array_Decays[86][2] = neutralino3amplitudeneut2ssbar; ParticleNeutralino3.Array_Decays[86][3] = 3; ParticleNeutralino3.Array_Comments[86] = "# ~chi_30 -> ~chi_20 s sbar";
   ParticleNeutralino3.Array_Decays[87][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[87][1] = PDGtop; ParticleNeutralino3.Array_Decays[87][4] = -PDGtop; ParticleNeutralino3.Array_Decays[87][2] = neutralino3amplitudeneut2ttbar; ParticleNeutralino3.Array_Decays[87][3] = 3; ParticleNeutralino3.Array_Comments[87] = "# ~chi_30 -> ~chi_20 t tbar";
   ParticleNeutralino3.Array_Decays[88][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[88][1] = PDGbottom; ParticleNeutralino3.Array_Decays[88][4] = -PDGbottom; ParticleNeutralino3.Array_Decays[88][2] = neutralino3amplitudeneut2bbbar; ParticleNeutralino3.Array_Decays[88][3] = 3; ParticleNeutralino3.Array_Comments[88] = "# ~chi_30 -> ~chi_20 b bbar";
   ParticleNeutralino3.Array_Decays[89][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[89][1] = PDGelectron; ParticleNeutralino3.Array_Decays[89][4] = -PDGelectron; ParticleNeutralino3.Array_Decays[89][2] = neutralino3amplitudeneut2eebar; ParticleNeutralino3.Array_Decays[89][3] = 3; ParticleNeutralino3.Array_Comments[89] = "# ~chi_30 -> ~chi_20 e- e+";
   ParticleNeutralino3.Array_Decays[90][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[90][1] = PDGmuon; ParticleNeutralino3.Array_Decays[90][4] = -PDGmuon; ParticleNeutralino3.Array_Decays[90][2] = neutralino3amplitudeneut2mumubar; ParticleNeutralino3.Array_Decays[90][3] = 3; ParticleNeutralino3.Array_Comments[90] = "# ~chi_30 -> ~chi_20 mu- mu+";
   ParticleNeutralino3.Array_Decays[91][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[91][1] = PDGtau; ParticleNeutralino3.Array_Decays[91][4] = -PDGtau; ParticleNeutralino3.Array_Decays[91][2] = neutralino3amplitudeneut2tautaubar; ParticleNeutralino3.Array_Decays[91][3] = 3; ParticleNeutralino3.Array_Comments[91] = "# ~chi_30 -> ~chi_20 tau- tau+";
   ParticleNeutralino3.Array_Decays[92][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[92][1] = PDGnuelectron; ParticleNeutralino3.Array_Decays[92][4] = -PDGnuelectron; ParticleNeutralino3.Array_Decays[92][2] = neutralino3amplitudeneut2nuenuebar; ParticleNeutralino3.Array_Decays[92][3] = 3; ParticleNeutralino3.Array_Comments[92] = "# ~chi_30 -> ~chi_20 nue nuebar";
   ParticleNeutralino3.Array_Decays[93][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[93][1] = PDGnumuon; ParticleNeutralino3.Array_Decays[93][4] = -PDGnumuon; ParticleNeutralino3.Array_Decays[93][2] = neutralino3amplitudeneut2numunumubar; ParticleNeutralino3.Array_Decays[93][3] = 3; ParticleNeutralino3.Array_Comments[93] = "# ~chi_30 -> ~chi_20 numu numubar";
   ParticleNeutralino3.Array_Decays[94][0] = PDGneutralino2; ParticleNeutralino3.Array_Decays[94][1] = PDGnutau; ParticleNeutralino3.Array_Decays[94][4] = -PDGnutau; ParticleNeutralino3.Array_Decays[94][2] = neutralino3amplitudeneut2nutaunutaubar; ParticleNeutralino3.Array_Decays[94][3] = 3; ParticleNeutralino3.Array_Comments[94] = "# ~chi_30 -> ~chi_20 nutau nutaubar";
   
   ParticleNeutralino3.Array_Decays[95][0] = PDGchargino1; ParticleNeutralino3.Array_Decays[95][1] = PDGup; ParticleNeutralino3.Array_Decays[95][4] = -PDGdown; ParticleNeutralino3.Array_Decays[95][2] = neutralino3amplitudechargino1udbar; ParticleNeutralino3.Array_Decays[95][3] = 3; ParticleNeutralino3.Array_Comments[95] = "# ~chi_30 -> chi_1- u db";
   ParticleNeutralino3.Array_Decays[96][0] = PDGchargino1; ParticleNeutralino3.Array_Decays[96][1] = PDGcharm; ParticleNeutralino3.Array_Decays[96][4] = -PDGstrange; ParticleNeutralino3.Array_Decays[96][2] = neutralino3amplitudechargino1csbar; ParticleNeutralino3.Array_Decays[96][3] = 3; ParticleNeutralino3.Array_Comments[96] = "# ~chi_30 -> chi_1- c sb";
   ParticleNeutralino3.Array_Decays[97][0] = PDGchargino1; ParticleNeutralino3.Array_Decays[97][1] = PDGnuelectron; ParticleNeutralino3.Array_Decays[97][4] = -PDGelectron; ParticleNeutralino3.Array_Decays[97][2] = neutralino3amplitudechargino1enuebar; ParticleNeutralino3.Array_Decays[97][3] = 3; ParticleNeutralino3.Array_Comments[97] = "# ~chi_30 -> chi_1- nu_e eb";
   ParticleNeutralino3.Array_Decays[98][0] = PDGchargino1; ParticleNeutralino3.Array_Decays[98][1] = PDGnumuon; ParticleNeutralino3.Array_Decays[98][4] = -PDGmuon; ParticleNeutralino3.Array_Decays[98][2] = neutralino3amplitudechargino1munumubar; ParticleNeutralino3.Array_Decays[98][3] = 3; ParticleNeutralino3.Array_Comments[98] = "# ~chi_30 -> chi_1- nu_mu mub";
   ParticleNeutralino3.Array_Decays[99][0] = PDGchargino1; ParticleNeutralino3.Array_Decays[99][1] = PDGnutau; ParticleNeutralino3.Array_Decays[99][4] = -PDGtau; ParticleNeutralino3.Array_Decays[99][2] = neutralino3amplitudechargino1taunutaubar; ParticleNeutralino3.Array_Decays[99][3] = 3; ParticleNeutralino3.Array_Comments[99] = "# ~chi_30 -> chi_1- nu_tau taub";
   ParticleNeutralino3.Array_Decays[100][0] = PDGchargino2; ParticleNeutralino3.Array_Decays[100][1] = PDGup; ParticleNeutralino3.Array_Decays[100][4] = -PDGdown; ParticleNeutralino3.Array_Decays[100][2] = neutralino3amplitudechargino2udbar; ParticleNeutralino3.Array_Decays[100][3] = 3; ParticleNeutralino3.Array_Comments[100] = "# ~chi_30 -> chi_2- u dbar";
   ParticleNeutralino3.Array_Decays[101][0] = PDGchargino2; ParticleNeutralino3.Array_Decays[101][1] = PDGcharm; ParticleNeutralino3.Array_Decays[101][4] = -PDGstrange; ParticleNeutralino3.Array_Decays[101][2] = neutralino3amplitudechargino2csbar; ParticleNeutralino3.Array_Decays[101][3] = 3; ParticleNeutralino3.Array_Comments[101] = "# ~chi_30 -> chi_2- c sbar";
   ParticleNeutralino3.Array_Decays[102][0] = PDGchargino2; ParticleNeutralino3.Array_Decays[102][1] = PDGnuelectron; ParticleNeutralino3.Array_Decays[102][4] = -PDGelectron; ParticleNeutralino3.Array_Decays[102][2] = neutralino3amplitudechargino2enuebar; ParticleNeutralino3.Array_Decays[102][3] = 3; ParticleNeutralino3.Array_Comments[102] = "# ~chi_30 -> chi_2- nu_e eb";
   ParticleNeutralino3.Array_Decays[103][0] = PDGchargino2; ParticleNeutralino3.Array_Decays[103][1] = PDGnumuon; ParticleNeutralino3.Array_Decays[103][4] = -PDGmuon; ParticleNeutralino3.Array_Decays[103][2] = neutralino3amplitudechargino2munumubar; ParticleNeutralino3.Array_Decays[103][3] = 3; ParticleNeutralino3.Array_Comments[103] = "# ~chi_30 -> chi_2- nu_mu mub";
   ParticleNeutralino3.Array_Decays[104][0] = PDGchargino2; ParticleNeutralino3.Array_Decays[104][1] = PDGnutau; ParticleNeutralino3.Array_Decays[104][4] = -PDGtau; ParticleNeutralino3.Array_Decays[104][2] = neutralino3amplitudechargino2taunutaubar; ParticleNeutralino3.Array_Decays[104][3] = 3; ParticleNeutralino3.Array_Comments[104] = "# ~chi_30 -> chi_2- nu_tau taubar";

   for(int i = 0; i<ParticleNeutralino3.No_of_Decays; i++) {
     if (ParticleNeutralino3.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleNeutralino3.Array_Comments[i] << " is negative = " << ParticleNeutralino3.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleNeutralino3.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   
   
   double Neut3_No_1to2_Decays = 0;
   
   Neut3_No_1to2_Decays = ParticleNeutralino3.No_1to2_Decays + ParticleNeutralino3.No_grav_Decays + ParticleNeutralino3.No_NMSSM_Decays;
   
   for (int j = 0; j<Neut3_No_1to2_Decays; j++) {
     ParticleNeutralino3.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<Neut3_No_1to2_Decays; j++) {
     ParticleNeutralino3.two_width = ParticleNeutralino3.two_width + ParticleNeutralino3.Array_Decays[j][2];
   }
   for (int j=Neut3_No_1to2_Decays; j<ParticleNeutralino3.No_of_Decays; j++) {
     ParticleNeutralino3.three_width = ParticleNeutralino3.three_width + ParticleNeutralino3.Array_Decays[j][2];
   }
   
   if ( ParticleNeutralino3.three_width != ParticleNeutralino3.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for neutralino 3 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleNeutralino3.No_of_Decays = Neut3_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleNeutralino3.total_width = ParticleNeutralino3.two_width;
     }
   else {
     ParticleNeutralino3.total_width = ParticleNeutralino3.two_width + ParticleNeutralino3.three_width;
   }
   
   if ( ParticleNeutralino3.total_width != ParticleNeutralino3.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleNeutralino3.No_of_Decays; i++) {
       //   fout << i << " " << ParticleNeutralino3.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in Neutralino3 total width \n");
     }
 }
 
 
 ///Neutralino4 Decays
 
 double neutralino4amplitudeuLubar=0, neutralino4amplitudeuRubar=0, neutralino4amplitudeuLbaru=0, neutralino4amplitudeuRbaru=0, neutralino4amplitudedLdbar=0, neutralino4amplitudedRdbar=0, neutralino4amplitudedLbard=0, neutralino4amplitudedRbard=0, neutralino4amplitudecLcbar=0, neutralino4amplitudecRcbar=0, neutralino4amplitudecLbarc=0, neutralino4amplitudecRbarc=0, neutralino4amplitudesLsbar=0, neutralino4amplitudesRsbar=0, neutralino4amplitudesLbars=0, neutralino4amplitudesRbars=0, neutralino4amplitudeeLebar=0, neutralino4amplitudeeRebar=0, neutralino4amplitudeeLbare=0, neutralino4amplitudeeRbare=0, neutralino4amplitudemuLmubar=0, neutralino4amplitudemuRmubar=0, neutralino4amplitudemuLbarmu=0, neutralino4amplitudemuRbarmu=0,neutralino4amplitudesnuenuebar=0, neutralino4amplitudesnuebarnue=0, neutralino4amplitudesnumunumubar=0, neutralino4amplitudesnumubarnumu=0, neutralino4amplitudetopstop1bar=0, neutralino4amplitudetopstop2bar=0, neutralino4amplitudetopbarstop1=0, neutralino4amplitudetopbarstop2=0, neutralino4amplitudebottomsbottom1bar=0, neutralino4amplitudebottomsbottom2bar=0, neutralino4amplitudebottombarsbottom1=0, neutralino4amplitudebottombarsbottom2=0, neutralino4amplitudetaustau1bar=0, neutralino4amplitudetaustau2bar=0, neutralino4amplitudetaubarstau1=0, neutralino4amplitudetaubarstau2=0, neutralino4amplitudenutausnutaubar=0, neutralino4amplitudenutaubarsnutau=0, neutralino4amplitudeWbosonpluscharginoW1=0, neutralino4amplitudeWbosonpluscharginoW2=0, neutralino4amplitudeWbosonminuscharginoW1=0, neutralino4amplitudeWbosonminuscharginoW2=0, neutralino4amplitudeHpluscharginoW1=0, neutralino4amplitudeHpluscharginoW2=0, neutralino4amplitudeHminuscharginoW1=0, neutralino4amplitudeHminuscharginoW2=0, neutralino4amplitudeZbosonneutralino1=0, neutralino4amplitudeZbosonneutralino2=0, neutralino4amplitudeZbosonneutralino3=0, neutralino4amplitudehneutralino1=0, neutralino4amplitudehneutralino2=0, neutralino4amplitudehneutralino3=0, neutralino4amplitudeHneutralino1=0, neutralino4amplitudeHneutralino2=0, neutralino4amplitudeHneutralino3=0, neutralino4amplitudeAneutralino1=0, neutralino4amplitudeAneutralino2=0, neutralino4amplitudeAneutralino3=0, neutralino4amplitudephotongravitino=0, neutralino4amplitudeZgravitino=0, neutralino4amplitudehgravitino=0, neutralino4amplitudeHgravitino=0, neutralino4amplitudeAgravitino=0;
 
 double neutralino4amplitudeneut1uubar=0, neutralino4amplitudeneut1ddbar=0, neutralino4amplitudeneut1ccbar=0, neutralino4amplitudeneut1ssbar=0, neutralino4amplitudeneut1ttbar=0, neutralino4amplitudeneut1bbbar=0, neutralino4amplitudeneut1eebar=0, neutralino4amplitudeneut1mumubar=0, neutralino4amplitudeneut1tautaubar=0, neutralino4amplitudeneut1nuenuebar=0, neutralino4amplitudeneut1numunumubar=0, neutralino4amplitudeneut1nutaunutaubar=0, neutralino4amplitudeneut2uubar=0, neutralino4amplitudeneut2ddbar=0, neutralino4amplitudeneut2ccbar=0, neutralino4amplitudeneut2ssbar=0, neutralino4amplitudeneut2ttbar=0, neutralino4amplitudeneut2bbbar=0, neutralino4amplitudeneut2eebar=0, neutralino4amplitudeneut2mumubar=0, neutralino4amplitudeneut2tautaubar=0, neutralino4amplitudeneut2nuenuebar=0, neutralino4amplitudeneut2numunumubar=0, neutralino4amplitudeneut2nutaunutaubar=0, neutralino4amplitudeneut3uubar=0, neutralino4amplitudeneut3ddbar=0, neutralino4amplitudeneut3ccbar=0, neutralino4amplitudeneut3ssbar=0, neutralino4amplitudeneut3ttbar=0, neutralino4amplitudeneut3bbbar=0, neutralino4amplitudeneut3eebar=0, neutralino4amplitudeneut3mumubar=0, neutralino4amplitudeneut3tautaubar=0, neutralino4amplitudeneut3nuenuebar=0, neutralino4amplitudeneut3numunumubar=0, neutralino4amplitudeneut3nutaunutaubar=0, neutralino4amplitudechargino1udbar=0, neutralino4amplitudechargino1csbar=0, neutralino4amplitudechargino1enuebar=0, neutralino4amplitudechargino1munumubar=0, neutralino4amplitudechargino1taunutaubar=0, neutralino4amplitudechargino2udbar=0, neutralino4amplitudechargino2csbar=0, neutralino4amplitudechargino2enuebar=0, neutralino4amplitudechargino2munumubar=0, neutralino4amplitudechargino2taunutaubar=0;
 
 double neutralino4amplitudeH3neutralino1 = 0, neutralino4amplitudeH3neutralino2 = 0, neutralino4amplitudeH3neutralino3 = 0, neutralino4amplitudeA2neutralino1 = 0, neutralino4amplitudeA2neutralino2 = 0, neutralino4amplitudeA2neutralino3 = 0;
 
 if (flagneut4 == 1) {
   if (nmssmIsIt == false) {
     neutralino4amplitudeuLubar = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mup, mu(1,1), g, gp, mixNeut, 1, 'L', 4);
     neutralino4amplitudeuRubar = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mup, mu(2,1), g, gp, mixNeut, 1, 'R', 4);
     neutralino4amplitudeuLbaru = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mup, mu(1,1), g, gp, mixNeut, 1, 'L', 4); 
     neutralino4amplitudeuRbaru = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mup, mu(2,1), g, gp, mixNeut, 1, 'R', 4);
     neutralino4amplitudedLdbar = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mdo, md(1,1), g, gp, mixNeut, 2, 'L', 4);
     neutralino4amplitudedRdbar = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mdo, md(2,1), g, gp, mixNeut, 2, 'R', 4);
     neutralino4amplitudedLbard = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mdo, md(1,1), g, gp, mixNeut, 2, 'L', 4);
     neutralino4amplitudedRbard = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mdo, md(2,1), g, gp, mixNeut, 2, 'R', 4);
     neutralino4amplitudecLcbar = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mc, mu(1,2), g, gp, mixNeut, 1, 'L', 4);
     neutralino4amplitudecRcbar = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mc, mu(2,2), g, gp, mixNeut, 1, 'R', 4);
     neutralino4amplitudecLbarc = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mc, mu(1,2), g, gp, mixNeut, 1, 'L', 4);
     neutralino4amplitudecRbarc = neutralinoamplitudedecayquarksquarkLorR (mneut(4), mc, mu(2,2), g, gp, mixNeut, 1, 'R', 4);
     neutralino4amplitudesLsbar = neutralinoamplitudedecayquarksquarkLorR (mneut(4), ms, md(1,2), g, gp, mixNeut, 2, 'L', 4);
     neutralino4amplitudesRsbar = neutralinoamplitudedecayquarksquarkLorR (mneut(4), ms, md(2,2), g, gp, mixNeut, 2, 'R', 4);
     neutralino4amplitudesLbars = neutralinoamplitudedecayquarksquarkLorR (mneut(4), ms, md(1,2), g, gp, mixNeut, 2, 'L', 4);
     neutralino4amplitudesRbars = neutralinoamplitudedecayquarksquarkLorR (mneut(4), ms, md(2,2), g, gp, mixNeut, 2, 'R', 4);
     neutralino4amplitudeeLebar = neutralinoamplitudedecayleptonsleptonLorR (mneut(4), mel, me(1,1), g, gp, mixNeut, 'L', 4);
     neutralino4amplitudeeRebar = neutralinoamplitudedecayleptonsleptonLorR (mneut(4), mel, me(2,1), g, gp, mixNeut, 'R', 4);
     neutralino4amplitudeeLbare = neutralinoamplitudedecayleptonsleptonLorR (mneut(4), mel, me(1,1), g, gp, mixNeut, 'L', 4);
     neutralino4amplitudeeRbare = neutralinoamplitudedecayleptonsleptonLorR (mneut(4), mel, me(2,1), g, gp, mixNeut, 'R', 4);
     neutralino4amplitudemuLmubar = neutralinoamplitudedecayleptonsleptonLorR (mneut(4), mmu, me(1,2), g, gp, mixNeut, 'L', 4);
     neutralino4amplitudemuRmubar = neutralinoamplitudedecayleptonsleptonLorR (mneut(4), mmu, me(2,2), g, gp, mixNeut, 'R', 4);
     neutralino4amplitudemuLbarmu = neutralinoamplitudedecayleptonsleptonLorR (mneut(4), mmu, me(1,2), g, gp, mixNeut, 'L', 4);
     neutralino4amplitudemuRbarmu = neutralinoamplitudedecayleptonsleptonLorR (mneut(4), mmu, me(2,2), g, gp, mixNeut, 'R', 4);
     neutralino4amplitudesnuenuebar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(4), 0, msnu(1), g, gp, mixNeut, 4);
     neutralino4amplitudesnuebarnue = neutralinoamplitudedecayneutrinosneutrinoL (mneut(4), 0, msnu(1), g, gp, mixNeut, 4);
     neutralino4amplitudesnumunumubar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(4), 0, msnu(2), g, gp, mixNeut, 4);
     neutralino4amplitudesnumubarnumu = neutralinoamplitudedecayneutrinosneutrinoL (mneut(4), 0, msnu(2), g, gp, mixNeut, 4);
     neutralino4amplitudetopstop1bar = neutralinoamplitudedecaysquark3quarkmix (mneut(4), mt, mu(1,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 4);
     neutralino4amplitudetopstop2bar = neutralinoamplitudedecaysquark3quarkmix (mneut(4), mt, mu(2,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 4);
     neutralino4amplitudetopbarstop1 = neutralinoamplitudedecaysquark3quarkmix (mneut(4), mt, mu(1,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 1, 4);
     neutralino4amplitudetopbarstop2 = neutralinoamplitudedecaysquark3quarkmix (mneut(4), mt, mu(2,3), runmw, thetat, beta, mixNeut, g, gp, runmt, 1, 2, 4);
     neutralino4amplitudebottomsbottom1bar = neutralinoamplitudedecaysquark3quarkmix (mneut(4), mb, md(1,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 4);
     neutralino4amplitudebottomsbottom2bar = neutralinoamplitudedecaysquark3quarkmix (mneut(4), mb, md(2,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 4);
     neutralino4amplitudebottombarsbottom1 = neutralinoamplitudedecaysquark3quarkmix (mneut(4), mb, md(1,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 1, 4);
     neutralino4amplitudebottombarsbottom2 = neutralinoamplitudedecaysquark3quarkmix (mneut(4), mb, md(2,3), runmw, thetab, beta, mixNeut, g, gp, runmb, 2, 2, 4);
     neutralino4amplitudetaustau1bar = neutralinoamplitudedecaystautau (mneut(4), mtau, me(1,3), runmw, thetatau, beta, mixNeut, g, gp, 1, 4);
     neutralino4amplitudetaustau2bar = neutralinoamplitudedecaystautau (mneut(4), mtau, me(2,3), runmw, thetatau, beta, mixNeut, g, gp, 2, 4);
     neutralino4amplitudetaubarstau1 = neutralinoamplitudedecaystautau (mneut(4), mtau, me(1,3), runmw, thetatau, beta, mixNeut, g, gp, 1, 4);
     neutralino4amplitudetaubarstau2 = neutralinoamplitudedecaystautau (mneut(4), mtau, me(2,3), runmw, thetatau, beta, mixNeut, g, gp, 2, 4);
     neutralino4amplitudenutausnutaubar = neutralinoamplitudedecayneutrinosneutrinoL (mneut(4), 0, msnu(3), g, gp, mixNeut, 4);
     neutralino4amplitudenutaubarsnutau = neutralinoamplitudedecayneutrinosneutrinoL (mneut(4), 0, msnu(3), g, gp, mixNeut, 4);
     neutralino4amplitudeWbosonpluscharginoW1 = neutralinoamplitudedecaycharginoWboson (mneut(4), polemw, MCH1, g, thetaL2, thetaR2, mixNeut, 4, 1);
     neutralino4amplitudeWbosonpluscharginoW2 = neutralinoamplitudedecaycharginoWboson (mneut(4), polemw, MCH2, g, thetaL2, thetaR2, mixNeut, 4, 2);
     neutralino4amplitudeWbosonminuscharginoW1 = neutralinoamplitudedecaycharginoWboson (mneut(4), polemw, MCH1, g, thetaL2, thetaR2, mixNeut, 4, 1);
     neutralino4amplitudeWbosonminuscharginoW2 = neutralinoamplitudedecaycharginoWboson (mneut(4), polemw, MCH2, g, thetaL2, thetaR2, mixNeut, 4, 2);
     neutralino4amplitudeHpluscharginoW1 = neutralinoamplitudedecaycharginoHplus (mneut(4), mHpm, MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 4, 1);
     neutralino4amplitudeHpluscharginoW2 = neutralinoamplitudedecaycharginoHplus (mneut(4), mHpm, MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 4, 2);
     neutralino4amplitudeHminuscharginoW1 = neutralinoamplitudedecaycharginoHplus (mneut(4), mHpm, MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 4, 1);
     neutralino4amplitudeHminuscharginoW2 = neutralinoamplitudedecaycharginoHplus (mneut(4), mHpm, MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 4, 2);
     neutralino4amplitudeZbosonneutralino1 = neutralinoamplitudedecayneutralinoZboson (mneut(4), polemz, mneut(1), g, gp, mixNeut, 4, 1);
     neutralino4amplitudeZbosonneutralino2 = neutralinoamplitudedecayneutralinoZboson (mneut(4), polemz, mneut(2), g, gp, mixNeut, 4, 2);
     neutralino4amplitudeZbosonneutralino3 = neutralinoamplitudedecayneutralinoZboson (mneut(4), polemz, mneut(3), g, gp, mixNeut, 4, 3);
     neutralino4amplitudehneutralino1 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(4), mh0(1), mneut(1), g, gp, mixNeut, alpha, 4, 1, 'h');
     neutralino4amplitudehneutralino2 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(4), mh0(1), mneut(2), g, gp, mixNeut, alpha, 4, 2, 'h');
     neutralino4amplitudehneutralino3 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(4), mh0(1), mneut(3), g, gp, mixNeut, alpha, 4, 3, 'h');
     neutralino4amplitudeHneutralino1 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(4), mh0(2), mneut(1), g, gp, mixNeut, alpha, 4, 1, 'H');
     neutralino4amplitudeHneutralino2 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(4), mh0(2), mneut(2), g, gp, mixNeut, alpha, 4, 2, 'H');
     neutralino4amplitudeHneutralino3 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(4), mh0(2), mneut(3), g, gp, mixNeut, alpha, 4, 3, 'H');
     neutralino4amplitudeAneutralino1 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(4), mA0(1), mneut(1), g, gp, mixNeut, beta, 4, 1, 'A');
     neutralino4amplitudeAneutralino2 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(4), mA0(1), mneut(2), g, gp, mixNeut, beta, 4, 2, 'A');
     neutralino4amplitudeAneutralino3 = neutralinoamplitudedecayneutralinoneutHiggs (mneut(4), mA0(1), mneut(3), g, gp, mixNeut, beta, 4, 3, 'A');
          
     neutralino4amplitudephotongravitino = neutralinoamplitudedecayphotongravitino (mneut(4), mgravitino, MPlreduced, mixNeut, g, gp, 4, gravonoff, neutNLSP);
     neutralino4amplitudeZgravitino = neutralinoamplitudedecayZgravitino(mneut(4), polemz, mgravitino, MPlreduced, mixNeut, g, gp, beta, 4, gravonoff, neutNLSP);
     neutralino4amplitudehgravitino = neutralinoamplitudedecayphigravitino(mneut(4), mh0(1), mgravitino, MPlreduced, mixNeut, alpha, beta, 4, gravonoff, 'h', neutNLSP);
     neutralino4amplitudeHgravitino = neutralinoamplitudedecayphigravitino(mneut(4), mh0(2), mgravitino, MPlreduced, mixNeut, alpha, beta, 4, gravonoff, 'H', neutNLSP);
     neutralino4amplitudeAgravitino = neutralinoamplitudedecayphigravitino(mneut(4), mA0(1), mgravitino, MPlreduced, mixNeut, alpha, beta, 4, gravonoff, 'A', neutNLSP);
     
     neutralino4amplitudeneut1uubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), mu(1,1), mu(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mup, alphas, 0, runmw, g, gp, alpha, beta, runmu, mixNeut, 4, 1, onetothree, 'u');
     neutralino4amplitudeneut1ddbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), md(1,1), md(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mdo, alphas, 0, runmw, g, gp, alpha, beta, runmd, mixNeut, 4, 1, onetothree, 'd');
     neutralino4amplitudeneut1ccbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), mu(1,2), mu(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mc, alphas, 0, runmw, g, gp, alpha, beta, runmc, mixNeut, 4, 1, onetothree, 'u');
     neutralino4amplitudeneut1ssbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), md(1,2), md(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(1), ms, alphas, 0, runmw, g, gp, alpha, beta, runms, mixNeut, 4, 1, onetothree, 'd');
     neutralino4amplitudeneut1ttbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), mu(1,3), mu(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mt, alphas, thetat, runmw, g, gp, alpha, beta, runmt, mixNeut, 4, 1, onetothree, 'u');
     neutralino4amplitudeneut1bbbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), md(1,3), md(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mb, alphas, thetab, runmw, g, gp, alpha, beta, runmb, mixNeut, 4, 1, onetothree, 'd');
     neutralino4amplitudeneut1eebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), me(1,1), me(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mel, alphas, 0, runmw, g, gp, alpha, beta, runmel, mixNeut, 4, 1, onetothree, 'l');
     neutralino4amplitudeneut1mumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), me(1,2), me(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mmu, alphas, 0, runmw, g, gp, alpha, beta, runmmu, mixNeut, 4, 1, onetothree, 'l');
     neutralino4amplitudeneut1tautaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), me(1,3), me(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(1), mtau, alphas, thetatau-PI/2, runmw, g, gp, alpha, beta, runmtau, mixNeut, 4, 1, onetothree, 'l');
     neutralino4amplitudeneut1nuenuebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), msnu(1), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(1), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 4, 1, onetothree, 'n');
     neutralino4amplitudeneut1numunumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), msnu(2), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(1), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 4, 1, onetothree, 'n');
     neutralino4amplitudeneut1nutaunutaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), msnu(3), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(1), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 4, 1, onetothree, 'n');
     neutralino4amplitudeneut2uubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), mu(1,1), mu(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mup, alphas, 0, runmw, g, gp, alpha, beta, runmu, mixNeut, 4, 2, onetothree, 'u');
     neutralino4amplitudeneut2ddbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), md(1,1), md(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mdo, alphas, 0, runmw, g, gp, alpha, beta, runmd, mixNeut, 4, 2, onetothree, 'd');
     neutralino4amplitudeneut2ccbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), mu(1,2), mu(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mc, alphas, 0, runmw, g, gp, alpha, beta, runmc, mixNeut, 4, 2, onetothree, 'u');
     neutralino4amplitudeneut2ssbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), md(1,2), md(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(2), ms, alphas, 0, runmw, g, gp, alpha, beta, runms, mixNeut, 4, 2, onetothree, 'd');
     neutralino4amplitudeneut2ttbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), mu(1,3), mu(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mt, alphas, thetat, runmw, g, gp, alpha, beta, runmt, mixNeut, 4, 2, onetothree, 'u');
     neutralino4amplitudeneut2bbbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), md(1,3), md(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mb, alphas, thetab, runmw, g, gp, alpha, beta, runmb, mixNeut, 4, 2, onetothree, 'd');
     neutralino4amplitudeneut2eebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), me(1,1), me(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mel, alphas, 0, runmw, g, gp, alpha, beta, runmel, mixNeut, 4, 2, onetothree, 'l');
     neutralino4amplitudeneut2mumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), me(1,2), me(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mmu, alphas, 0, runmw, g, gp, alpha, beta, runmmu, mixNeut, 4, 2, onetothree, 'l');
     neutralino4amplitudeneut2tautaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), me(1,3), me(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(2), mtau, alphas, thetatau-PI/2, runmw, g, gp, alpha, beta, runmtau, mixNeut, 4, 2, onetothree, 'l');
     neutralino4amplitudeneut2nuenuebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), msnu(1), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(2), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 4, 2, onetothree, 'n');
     neutralino4amplitudeneut2numunumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), msnu(2), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(2), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 4, 2, onetothree, 'n');
     neutralino4amplitudeneut2nutaunutaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), msnu(3), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(2), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 4, 2, onetothree, 'n');
     neutralino4amplitudeneut3uubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), mu(1,1), mu(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(3), mup, alphas, 0, runmw, g, gp, alpha, beta, runmu, mixNeut, 4, 3, onetothree, 'u');
     neutralino4amplitudeneut3ddbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), md(1,1), md(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(3), mdo, alphas, 0, runmw, g, gp, alpha, beta, runmd, mixNeut, 4, 3, onetothree, 'd');
     neutralino4amplitudeneut3ccbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), mu(1,2), mu(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(3), mc, alphas, 0, runmw, g, gp, alpha, beta, runmc, mixNeut, 4, 3, onetothree, 'u');
     neutralino4amplitudeneut3ssbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), md(1,2), md(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(3), ms, alphas, 0, runmw, g, gp, alpha, beta, runms, mixNeut, 4, 3, onetothree, 'd');
     neutralino4amplitudeneut3ttbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), mu(1,3), mu(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(3), mt, alphas, thetat, runmw, g, gp, alpha, beta, runmt, mixNeut, 4, 3, onetothree, 'u');
     neutralino4amplitudeneut3bbbar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), md(1,3), md(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(3), mb, alphas, thetab, runmw, g, gp, alpha, beta, runmb, mixNeut, 4, 3, onetothree, 'd');
     neutralino4amplitudeneut3eebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), me(1,1), me(2,1), runmz, mh0(1), mh0(2), mA0(1), mneut(3), mel, alphas, 0, runmw, g, gp, alpha, beta, runmel, mixNeut, 4, 3, onetothree, 'l');
     neutralino4amplitudeneut3mumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), me(1,2), me(2,2), runmz, mh0(1), mh0(2), mA0(1), mneut(3), mmu, alphas, 0, runmw, g, gp, alpha, beta, runmmu, mixNeut, 4, 3, onetothree, 'l');
     neutralino4amplitudeneut3tautaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), me(1,3), me(2,3), runmz, mh0(1), mh0(2), mA0(1), mneut(3), mtau, alphas, thetatau-PI/2, runmw, g, gp, alpha, beta, runmtau, mixNeut, 4, 3, onetothree, 'l');
     neutralino4amplitudeneut3nuenuebar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), msnu(1), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(3), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 4, 3, onetothree, 'n');
     neutralino4amplitudeneut3numunumubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), msnu(2), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(3), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 4, 3, onetothree, 'n');
     neutralino4amplitudeneut3nutaunutaubar = neutralinoamplitudedecaydgaussneutralinoffbar (mneut(4), msnu(3), 100000000000, runmz, mh0(1), mh0(2), mA0(1), mneut(3), 0, alphas, 0, runmw, g, gp, alpha, beta, 0, mixNeut, 4, 3, onetothree, 'n');
     
     neutralino4amplitudechargino1udbar = neutralinoamplitudedecaycharginoffprimebar (mneut(4), mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mch(1), mup, mdo, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 4, 1, onetothree, 'q', 'n');
     neutralino4amplitudechargino1csbar = neutralinoamplitudedecaycharginoffprimebar (mneut(4), mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mch(1), mc, ms, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 4, 1, onetothree, 'q', 'n');
     neutralino4amplitudechargino1enuebar = neutralinoamplitudedecaycharginoffprimebar (mneut(4), me(1,1), me(2,1), msnu(1), 100000000000, polemw, mHpm, mch(1), 0, mel, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 4, 1, onetothree, 'l', 'n');
     neutralino4amplitudechargino1munumubar = neutralinoamplitudedecaycharginoffprimebar (mneut(4), me(1,2), me(2,2), msnu(2), 100000000000, polemw, mHpm, mch(1), 0, mmu, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 4, 1, onetothree, 'l', 'n');
     neutralino4amplitudechargino1taunutaubar = neutralinoamplitudedecaycharginoffprimebar (mneut(4), me(1,3), me(2,3), msnu(3), 100000000000, polemw, mHpm, mch(1), 0, mtau, 0, thetatau-PI/2, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 4, 1, onetothree, 'l', 'n');
     neutralino4amplitudechargino2udbar = neutralinoamplitudedecaycharginoffprimebar (mneut(4), mu(1,1), mu(2,1), md(1,1), md(2,1), polemw, mHpm, mch(2), mup, mdo, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmu, runmd, mixNeut, 4, 2, onetothree, 'q', 'n');
     neutralino4amplitudechargino2csbar = neutralinoamplitudedecaycharginoffprimebar (mneut(4), mu(1,2), mu(2,2), md(1,2), md(2,2), polemw, mHpm, mch(2), mc, ms, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, runmc, runms, mixNeut, 4, 2, onetothree, 'q', 'n');
     neutralino4amplitudechargino2enuebar = neutralinoamplitudedecaycharginoffprimebar (mneut(4), me(1,1), me(2,1), msnu(1), 100000000000, polemw, mHpm, mch(2), 0, mel, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmel, mixNeut, 4, 2, onetothree, 'l', 'n');
     neutralino4amplitudechargino2munumubar = neutralinoamplitudedecaycharginoffprimebar (mneut(4), me(1,2), me(2,2), msnu(2), 100000000000, polemw, mHpm, mch(2), 0, mmu, 0, 0, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmmu, mixNeut, 4, 2, onetothree, 'l', 'n');
     neutralino4amplitudechargino2taunutaubar = neutralinoamplitudedecaycharginoffprimebar (mneut(4), me(1,3), me(2,3), msnu(3), 100000000000, polemw, mHpm, mch(2), 0, mtau, 0, thetatau-PI/2, g, gp, alphas, beta, thetaL2, thetaR2, 0, runmtau, mixNeut, 4, 2, onetothree, 'l', 'n');
   }
   else if (nmssmIsIt == true){
     neutralino4amplitudeWbosonpluscharginoW1 = neutralinoamplitudedecaycharginoWNMSSM (mneut(4), MCH1, polemw, g, thetaL2, thetaR2, mixNeut, 4, 1);
     neutralino4amplitudeWbosonpluscharginoW2 = neutralinoamplitudedecaycharginoWNMSSM (mneut(4), MCH2, polemw, g, thetaL2, thetaR2, mixNeut, 4, 2);
     neutralino4amplitudeWbosonminuscharginoW1 = neutralino4amplitudeWbosonpluscharginoW1;
     neutralino4amplitudeWbosonminuscharginoW2 = neutralino4amplitudeWbosonpluscharginoW2;
     neutralino4amplitudeZbosonneutralino1 = neutralinoamplitudedecayneutralinoZNMSSM (mneut(4), mneut(1), polemz, g, gp, mixNeut, 4, 1);
     neutralino4amplitudeZbosonneutralino2 = neutralinoamplitudedecayneutralinoZNMSSM (mneut(4), mneut(2), polemz, g, gp, mixNeut, 4, 2);
     neutralino4amplitudeZbosonneutralino3 = neutralinoamplitudedecayneutralinoZNMSSM (mneut(4), mneut(3), polemz, g, gp, mixNeut, 4, 3);
     
     neutralino4amplitudeHpluscharginoW1 = neutralinoamplitudecharginoHpmNMSSM (mneut(4), MCH1, mHpm, g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 4, 1);
     neutralino4amplitudeHpluscharginoW2 = neutralinoamplitudecharginoHpmNMSSM (mneut(4), MCH2, mHpm, g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 4, 2);
     neutralino4amplitudeHminuscharginoW1 = neutralino4amplitudeHpluscharginoW1;
     neutralino4amplitudeHminuscharginoW2 = neutralino4amplitudeHpluscharginoW2;
     
     neutralino4amplitudehneutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(4), mneut(1), mh0(1), g, gp, lam, kappa, mixNeut, CPEMix, 4, 1, 1);
     neutralino4amplitudeHneutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(4), mneut(1), mh0(2), g, gp, lam, kappa, mixNeut, CPEMix, 4, 1, 2);
     neutralino4amplitudeH3neutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(4), mneut(1), mh0(3), g, gp, lam, kappa, mixNeut, CPEMix, 4, 1, 3);
     neutralino4amplitudehneutralino2 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(4), mneut(2), mh0(1), g, gp, lam, kappa, mixNeut, CPEMix, 4, 2, 1);
     neutralino4amplitudeHneutralino2 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(4), mneut(2), mh0(2), g, gp, lam, kappa, mixNeut, CPEMix, 4, 2, 2);
     neutralino4amplitudeH3neutralino2 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(4), mneut(2), mh0(3), g, gp, lam, kappa, mixNeut, CPEMix, 4, 2, 3);
     neutralino4amplitudehneutralino3 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(4), mneut(3), mh0(1), g, gp, lam, kappa, mixNeut, CPEMix, 4, 3, 1);
     neutralino4amplitudeHneutralino3 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(4), mneut(3), mh0(2), g, gp, lam, kappa, mixNeut, CPEMix, 4, 3, 2);
     neutralino4amplitudeH3neutralino3 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(4), mneut(3), mh0(3), g, gp, lam, kappa, mixNeut, CPEMix, 4, 3, 3);
     
     neutralino4amplitudeAneutralino1 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(4), mneut(1), mA0(1), g, gp, lam, kappa, mixNeut, CPOMix, 4, 1, 1);
     neutralino4amplitudeA2neutralino1 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(4), mneut(1), mA0(2), g, gp, lam, kappa, mixNeut, CPOMix, 4, 1, 2);
     neutralino4amplitudeAneutralino2 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(4), mneut(2), mA0(1), g, gp, lam, kappa, mixNeut, CPOMix, 4, 2, 1);
     neutralino4amplitudeA2neutralino2 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(4), mneut(2), mA0(2), g, gp, lam, kappa, mixNeut, CPOMix, 4, 2, 2);
     neutralino4amplitudeAneutralino3 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(4), mneut(3), mA0(1), g, gp, lam, kappa, mixNeut, CPOMix, 4, 3, 1);
     neutralino4amplitudeA2neutralino3 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(4), mneut(3), mA0(2), g, gp, lam, kappa, mixNeut, CPOMix, 4, 3, 2);
     
     neutralino4amplitudeuLubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), mu(1,1), mup, g, gp, mixNeut, 4, 'u', 'L');
     neutralino4amplitudeuRubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), mu(2,1), mup, g, gp, mixNeut, 4, 'u', 'R');
     neutralino4amplitudeuLbaru = neutralino4amplitudeuLubar;
     neutralino4amplitudeuRbaru = neutralino4amplitudeuRubar;
     neutralino4amplitudedLdbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), md(1,1), mdo, g, gp, mixNeut, 4, 'd', 'L');
     neutralino4amplitudedRdbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), md(2,1), mdo, g, gp, mixNeut, 4, 'd', 'R');
     neutralino4amplitudedLbard = neutralino4amplitudedLdbar;
     neutralino4amplitudedRbard = neutralino4amplitudedRdbar;
     neutralino4amplitudecLcbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), mu(1,2), mc, g, gp, mixNeut, 4, 'u', 'L');
     neutralino4amplitudecRcbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), mu(2,2), mc, g, gp, mixNeut, 4, 'u', 'R');
     neutralino4amplitudecLbarc = neutralino4amplitudecLcbar;
     neutralino4amplitudecRbarc = neutralino4amplitudecRcbar;
     neutralino4amplitudesLsbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), md(1,2), ms, g, gp, mixNeut, 4, 'd', 'L');
     neutralino4amplitudesRsbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), md(2,2), ms, g, gp, mixNeut, 4, 'd', 'R');
     neutralino4amplitudesLbars = neutralino4amplitudesLsbar;
     neutralino4amplitudesRbars = neutralino4amplitudesRsbar;
     neutralino4amplitudeeLebar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), me(1,1), mel, g, gp, mixNeut, 4, 'l', 'L');
     neutralino4amplitudeeRebar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), me(2,1), mel, g, gp, mixNeut, 4, 'l', 'R');
     neutralino4amplitudeeLbare = neutralino4amplitudeeLebar;
     neutralino4amplitudeeRbare = neutralino4amplitudeeRebar;
     neutralino4amplitudemuLmubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), me(1,2), mmu, g, gp, mixNeut, 4, 'l', 'L');
     neutralino4amplitudemuRmubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(4), me(2,2), mmu, g, gp, mixNeut, 4, 'l', 'R');
     neutralino4amplitudemuLbarmu = neutralino4amplitudemuLmubar;
     neutralino4amplitudemuRbarmu = neutralino4amplitudemuRmubar;
     
     neutralino4amplitudetopstop1bar = neutralinoamplitudestoptopNMSSM (mneut(4), mu(1,3), mt, g, gp, thetat, beta, runmw, mixNeut, runmt, 4, 1);
     neutralino4amplitudetopstop2bar = neutralinoamplitudestoptopNMSSM (mneut(4), mu(2,3), mt, g, gp, thetat, beta, runmw, mixNeut, runmt, 4, 2);
     neutralino4amplitudetopbarstop1 = neutralino4amplitudetopstop1bar;
     neutralino4amplitudetopbarstop2 = neutralino4amplitudetopstop2bar;
     neutralino4amplitudebottomsbottom1bar = neutralinoamplitudesbottombottomNMSSM (mneut(4), md(1,3), mb, g, gp, thetab, beta, runmw, mixNeut, runmb, 4, 1);
     neutralino4amplitudebottomsbottom2bar = neutralinoamplitudesbottombottomNMSSM (mneut(4), md(2,3), mb, g, gp, thetab, beta, runmw, mixNeut, runmb, 4, 2);
     neutralino4amplitudebottombarsbottom1 = neutralino4amplitudebottomsbottom1bar;
     neutralino4amplitudebottombarsbottom2 = neutralino4amplitudebottomsbottom2bar;
     
     neutralino4amplitudetaustau1bar = neutralinoamplitudestautauNMSSM (mneut(4), me(1,3), mtau, g, gp, thetatau, beta, runmw, mixNeut, runmtau, 4, 1);
     neutralino4amplitudetaustau2bar = neutralinoamplitudestautauNMSSM (mneut(4), me(2,3), mtau, g, gp, thetatau, beta, runmw, mixNeut, runmtau, 4, 2);
     neutralino4amplitudetaubarstau1 = neutralino4amplitudetaustau1bar;
     neutralino4amplitudetaubarstau2 = neutralino4amplitudetaustau2bar;
     
     neutralino4amplitudesnuenuebar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(4), msnu(1), 0, g, gp, mixNeut, 4);
     neutralino4amplitudesnuebarnue = neutralino4amplitudesnuenuebar;
     neutralino4amplitudesnumunumubar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(4), msnu(2), 0, g, gp, mixNeut, 4);
     neutralino4amplitudesnumubarnumu = neutralino4amplitudesnumunumubar;
     neutralino4amplitudenutausnutaubar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(4), msnu(3), 0, g, gp, mixNeut, 4);
     neutralino4amplitudenutaubarsnutau = neutralino4amplitudenutausnutaubar;
   }
   
   ParticleNeutralino4.Array_Decays[0][0] = -PDGup; ParticleNeutralino4.Array_Decays[0][1] = PDGsupL; ParticleNeutralino4.Array_Decays[0][2] = neutralino4amplitudeuLubar; ParticleNeutralino4.Array_Decays[0][3] = 2; ParticleNeutralino4.Array_Comments[0] = "# ~chi_40 -> ub ~u_L";
   ParticleNeutralino4.Array_Decays[1][0] = -PDGup; ParticleNeutralino4.Array_Decays[1][1] = PDGsupR; ParticleNeutralino4.Array_Decays[1][2] = neutralino4amplitudeuRubar; ParticleNeutralino4.Array_Decays[1][3] = 2; ParticleNeutralino4.Array_Comments[1] = "# ~chi_40 -> ub ~u_R";
   ParticleNeutralino4.Array_Decays[2][0] = PDGup; ParticleNeutralino4.Array_Decays[2][1] = -PDGsupL; ParticleNeutralino4.Array_Decays[2][2] = neutralino4amplitudeuLbaru; ParticleNeutralino4.Array_Decays[2][3] = 2; ParticleNeutralino4.Array_Comments[2] = "# ~chi_40 -> u ~u_L*";
   ParticleNeutralino4.Array_Decays[3][0] = PDGup; ParticleNeutralino4.Array_Decays[3][1] = -PDGsupR; ParticleNeutralino4.Array_Decays[3][2] = neutralino4amplitudeuRbaru; ParticleNeutralino4.Array_Decays[3][3] = 2; ParticleNeutralino4.Array_Comments[3] = "# ~chi_40 -> u ~u_R*";
   ParticleNeutralino4.Array_Decays[4][0] = -PDGdown; ParticleNeutralino4.Array_Decays[4][1] = PDGsdownL; ParticleNeutralino4.Array_Decays[4][2] = neutralino4amplitudedLdbar;ParticleNeutralino4.Array_Decays[4][3] = 2;  ParticleNeutralino4.Array_Comments[4] = "# ~chi_40 -> db ~d_L";
   ParticleNeutralino4.Array_Decays[5][0] = -PDGdown; ParticleNeutralino4.Array_Decays[5][1] = PDGsdownR; ParticleNeutralino4.Array_Decays[5][2] = neutralino4amplitudedRdbar; ParticleNeutralino4.Array_Decays[5][3] = 2; ParticleNeutralino4.Array_Comments[5] = "# ~chi_40 -> db ~d_R";
   ParticleNeutralino4.Array_Decays[6][0] = PDGdown; ParticleNeutralino4.Array_Decays[6][1] = -PDGsdownL; ParticleNeutralino4.Array_Decays[6][2] = neutralino4amplitudedLbard; ParticleNeutralino4.Array_Decays[6][3] = 2; ParticleNeutralino4.Array_Comments[6] = "# ~chi_40 -> d ~d_L*";
   ParticleNeutralino4.Array_Decays[7][0] = PDGdown; ParticleNeutralino4.Array_Decays[7][1] = -PDGsdownR; ParticleNeutralino4.Array_Decays[7][2] = neutralino4amplitudedRbard; ParticleNeutralino4.Array_Decays[7][3] = 2; ParticleNeutralino4.Array_Comments[7] = "# ~chi_40 -> d ~d_R*";
   ParticleNeutralino4.Array_Decays[8][0] = -PDGcharm; ParticleNeutralino4.Array_Decays[8][1] = PDGscharmL; ParticleNeutralino4.Array_Decays[8][2] = neutralino4amplitudecLcbar; ParticleNeutralino4.Array_Decays[8][3] = 2; ParticleNeutralino4.Array_Comments[8] = "# ~chi_40 -> cb ~c_L";
   ParticleNeutralino4.Array_Decays[9][0] = -PDGcharm; ParticleNeutralino4.Array_Decays[9][1] = PDGscharmR; ParticleNeutralino4.Array_Decays[9][2] = neutralino4amplitudecRcbar; ParticleNeutralino4.Array_Decays[9][3] = 2; ParticleNeutralino4.Array_Comments[9] = "# ~chi_40 -> cb ~c_R";
   ParticleNeutralino4.Array_Decays[10][0] = PDGcharm; ParticleNeutralino4.Array_Decays[10][1] = -PDGscharmL; ParticleNeutralino4.Array_Decays[10][2] = neutralino4amplitudecLbarc; ParticleNeutralino4.Array_Decays[10][3] = 2; ParticleNeutralino4.Array_Comments[10] = "# ~chi_40 -> c ~c_L*";
   ParticleNeutralino4.Array_Decays[11][0] = PDGcharm; ParticleNeutralino4.Array_Decays[11][1] = -PDGscharmR; ParticleNeutralino4.Array_Decays[11][2] = neutralino4amplitudecRbarc; ParticleNeutralino4.Array_Decays[11][3] = 2; ParticleNeutralino4.Array_Comments[11] = "# ~chi_40 -> c ~c_R*";
   ParticleNeutralino4.Array_Decays[12][0] = -PDGstrange; ParticleNeutralino4.Array_Decays[12][1] = PDGsstrangeL; ParticleNeutralino4.Array_Decays[12][2] = neutralino4amplitudesLsbar; ParticleNeutralino4.Array_Decays[12][3] = 2; ParticleNeutralino4.Array_Comments[12] = "# ~chi_40 -> sb ~s_L"; 
   ParticleNeutralino4.Array_Decays[13][0] = -PDGstrange; ParticleNeutralino4.Array_Decays[13][1] = PDGsstrangeR; ParticleNeutralino4.Array_Decays[13][2] = neutralino4amplitudesRsbar; ParticleNeutralino4.Array_Decays[13][3] = 2; ParticleNeutralino4.Array_Comments[13] = "# ~chi_40 -> sb ~s_R";
   ParticleNeutralino4.Array_Decays[14][0] = PDGstrange; ParticleNeutralino4.Array_Decays[14][1] = -PDGsstrangeL; ParticleNeutralino4.Array_Decays[14][2] = neutralino4amplitudesLbars; ParticleNeutralino4.Array_Decays[14][3] = 2; ParticleNeutralino4.Array_Comments[14] = "# ~chi_40 -> s ~s_L*";
   ParticleNeutralino4.Array_Decays[15][0] = PDGstrange; ParticleNeutralino4.Array_Decays[15][1] = -PDGsstrangeR; ParticleNeutralino4.Array_Decays[15][2] = neutralino4amplitudesRbars; ParticleNeutralino4.Array_Decays[15][3] = 2; ParticleNeutralino4.Array_Comments[15] = "# ~chi_40 -> s ~s_R*";
   ParticleNeutralino4.Array_Decays[16][0] = -PDGelectron; ParticleNeutralino4.Array_Decays[16][1] = PDGselectronL; ParticleNeutralino4.Array_Decays[16][2] = neutralino4amplitudeeLebar; ParticleNeutralino4.Array_Decays[16][3] = 2; ParticleNeutralino4.Array_Comments[16] = "# ~chi_40 -> e+ ~e_L-";
   ParticleNeutralino4.Array_Decays[17][0] = -PDGelectron; ParticleNeutralino4.Array_Decays[17][1] = PDGselectronR; ParticleNeutralino4.Array_Decays[17][2] = neutralino4amplitudeeRebar; ParticleNeutralino4.Array_Decays[17][3] = 2; ParticleNeutralino4.Array_Comments[17] = "# ~chi_40 -> e+ ~e_R-";
   ParticleNeutralino4.Array_Decays[18][0] = PDGelectron; ParticleNeutralino4.Array_Decays[18][1] = -PDGselectronL; ParticleNeutralino4.Array_Decays[18][2] = neutralino4amplitudeeLbare; ParticleNeutralino4.Array_Decays[18][3] = 2; ParticleNeutralino4.Array_Comments[18] = "# ~chi_40 -> e ~e_L+";
   ParticleNeutralino4.Array_Decays[19][0] = PDGelectron; ParticleNeutralino4.Array_Decays[19][1] = -PDGselectronR; ParticleNeutralino4.Array_Decays[19][2] = neutralino4amplitudeeRbare; ParticleNeutralino4.Array_Decays[19][3] = 2; ParticleNeutralino4.Array_Comments[19] = "# ~chi_40 -> e ~e_R+";   
   ParticleNeutralino4.Array_Decays[20][0] = -PDGmuon; ParticleNeutralino4.Array_Decays[20][1] = PDGsmuonL; ParticleNeutralino4.Array_Decays[20][2] = neutralino4amplitudemuLmubar; ParticleNeutralino4.Array_Decays[20][3] = 2; ParticleNeutralino4.Array_Comments[20] = "# ~chi_40 -> mu+ ~mu_L-";
   ParticleNeutralino4.Array_Decays[21][0] = -PDGmuon; ParticleNeutralino4.Array_Decays[21][1] = PDGsmuonR; ParticleNeutralino4.Array_Decays[21][2] = neutralino4amplitudemuRmubar; ParticleNeutralino4.Array_Decays[21][3] = 2; ParticleNeutralino4.Array_Comments[21] = "# ~chi_40 -> mu+ ~mu_R-";
   ParticleNeutralino4.Array_Decays[22][0] = PDGmuon; ParticleNeutralino4.Array_Decays[22][1] = -PDGsmuonL; ParticleNeutralino4.Array_Decays[22][2] = neutralino4amplitudemuLbarmu; ParticleNeutralino4.Array_Decays[22][3] = 2; ParticleNeutralino4.Array_Comments[22] = "# ~chi_40 -> mu- ~mu_L+";
   ParticleNeutralino4.Array_Decays[23][0] = PDGmuon; ParticleNeutralino4.Array_Decays[23][1] = -PDGsmuonR; ParticleNeutralino4.Array_Decays[23][2] = neutralino4amplitudemuRbarmu; ParticleNeutralino4.Array_Decays[23][3] = 2; ParticleNeutralino4.Array_Comments[23] = "# ~chi_40 -> mu- ~mu_R+";
   ParticleNeutralino4.Array_Decays[24][0] = PDGnuelectron; ParticleNeutralino4.Array_Decays[24][1] = -PDGnuselectronL; ParticleNeutralino4.Array_Decays[24][2] = neutralino4amplitudesnuebarnue; ParticleNeutralino4.Array_Decays[24][3] = 2; ParticleNeutralino4.Array_Comments[24] = "# ~chi_40 -> nu_e ~nu_eL*";
   ParticleNeutralino4.Array_Decays[25][0] = -PDGnuelectron; ParticleNeutralino4.Array_Decays[25][1] = PDGnuselectronL; ParticleNeutralino4.Array_Decays[25][2] = neutralino4amplitudesnuenuebar; ParticleNeutralino4.Array_Decays[25][3] = 2; ParticleNeutralino4.Array_Comments[25] = "# ~chi_40 -> nu_eb ~nu_eL";
   ParticleNeutralino4.Array_Decays[26][0] = PDGnumuon; ParticleNeutralino4.Array_Decays[26][1] = -PDGnusmuonL; ParticleNeutralino4.Array_Decays[26][2] = neutralino4amplitudesnumubarnumu; ParticleNeutralino4.Array_Decays[26][3] = 2; ParticleNeutralino4.Array_Comments[26] = "# ~chi_40 -> numu ~nu_muL*";
   ParticleNeutralino4.Array_Decays[27][0] = -PDGnumuon; ParticleNeutralino4.Array_Decays[27][1] = PDGnusmuonL; ParticleNeutralino4.Array_Decays[27][2] = neutralino4amplitudesnumunumubar; ParticleNeutralino4.Array_Decays[27][3] = 2; ParticleNeutralino4.Array_Comments[27] = "# ~chi_40 -> nu_mub ~nu_muL";
   ParticleNeutralino4.Array_Decays[28][0] = PDGtop; ParticleNeutralino4.Array_Decays[28][1] = -PDGstop1; ParticleNeutralino4.Array_Decays[28][2] = neutralino4amplitudetopstop1bar; ParticleNeutralino4.Array_Decays[28][3] = 2; ParticleNeutralino4.Array_Comments[28] = "# ~chi_40 -> t ~t_1*";
   ParticleNeutralino4.Array_Decays[29][0] = PDGtop; ParticleNeutralino4.Array_Decays[29][1] = -PDGstop2; ParticleNeutralino4.Array_Decays[29][2] = neutralino4amplitudetopstop2bar; ParticleNeutralino4.Array_Decays[29][3] = 2; ParticleNeutralino4.Array_Comments[29] = "# ~chi_40 -> t ~t_2*";
   ParticleNeutralino4.Array_Decays[30][0] = -PDGtop; ParticleNeutralino4.Array_Decays[30][1] = PDGstop1; ParticleNeutralino4.Array_Decays[30][2] = neutralino4amplitudetopbarstop1; ParticleNeutralino4.Array_Decays[30][3] = 2; ParticleNeutralino4.Array_Comments[30] = "# ~chi_40 -> tb ~t_1";
   ParticleNeutralino4.Array_Decays[31][0] = -PDGtop; ParticleNeutralino4.Array_Decays[31][1] = PDGstop2; ParticleNeutralino4.Array_Decays[31][2] = neutralino4amplitudetopbarstop2; ParticleNeutralino4.Array_Decays[31][3] = 2; ParticleNeutralino4.Array_Comments[31] = "# ~chi_40 -> tb ~t_2";
   ParticleNeutralino4.Array_Decays[32][0] = PDGbottom; ParticleNeutralino4.Array_Decays[32][1] = -PDGsbottom1; ParticleNeutralino4.Array_Decays[32][2] = neutralino4amplitudebottomsbottom1bar; ParticleNeutralino4.Array_Decays[32][3] = 2; ParticleNeutralino4.Array_Comments[32] = "# ~chi_40 -> b ~b_1*";
   ParticleNeutralino4.Array_Decays[33][0] = PDGbottom; ParticleNeutralino4.Array_Decays[33][1] = -PDGsbottom2; ParticleNeutralino4.Array_Decays[33][2] = neutralino4amplitudebottomsbottom2bar; ParticleNeutralino4.Array_Decays[33][3] = 2; ParticleNeutralino4.Array_Comments[33] = "# ~chi_40 -> b ~b_2*";
   ParticleNeutralino4.Array_Decays[34][0] = -PDGbottom; ParticleNeutralino4.Array_Decays[34][1] = PDGsbottom1; ParticleNeutralino4.Array_Decays[34][2] = neutralino4amplitudebottombarsbottom1; ParticleNeutralino4.Array_Decays[34][3] = 2; ParticleNeutralino4.Array_Comments[34] = "# ~chi_40 -> bb ~b_1";
   ParticleNeutralino4.Array_Decays[35][0] = -PDGbottom; ParticleNeutralino4.Array_Decays[35][1] = PDGsbottom2; ParticleNeutralino4.Array_Decays[35][2] = neutralino4amplitudebottombarsbottom2; ParticleNeutralino4.Array_Decays[35][3] = 2; ParticleNeutralino4.Array_Comments[35] = "# ~chi_40 -> bb ~b_2";
   ParticleNeutralino4.Array_Decays[36][0] = -PDGstau1; ParticleNeutralino4.Array_Decays[36][1] = PDGtau; ParticleNeutralino4.Array_Decays[36][2] = neutralino4amplitudetaustau1bar; ParticleNeutralino4.Array_Decays[36][3] = 2; ParticleNeutralino4.Array_Comments[36] = "# ~chi_40 -> tau- ~tau_1+";
   ParticleNeutralino4.Array_Decays[37][0] = -PDGstau2; ParticleNeutralino4.Array_Decays[37][1] = PDGtau; ParticleNeutralino4.Array_Decays[37][2] = neutralino4amplitudetaustau2bar; ParticleNeutralino4.Array_Decays[37][3] = 2; ParticleNeutralino4.Array_Comments[37] = "# ~chi_40 -> tau- ~tau_2+";
   ParticleNeutralino4.Array_Decays[38][0] = PDGstau1; ParticleNeutralino4.Array_Decays[38][1] = -PDGtau; ParticleNeutralino4.Array_Decays[38][2] = neutralino4amplitudetaubarstau1; ParticleNeutralino4.Array_Decays[38][3] = 2; ParticleNeutralino4.Array_Comments[38] = "# ~chi_40 -> tau+ ~tau_1-";
   ParticleNeutralino4.Array_Decays[39][0] = PDGstau2; ParticleNeutralino4.Array_Decays[39][1] = -PDGtau; ParticleNeutralino4.Array_Decays[39][2] = neutralino4amplitudetaubarstau2; ParticleNeutralino4.Array_Decays[39][3] = 2; ParticleNeutralino4.Array_Comments[39] = "# ~chi_40 -> tau+ ~tau_2-";
   ParticleNeutralino4.Array_Decays[40][0] = PDGnutau; ParticleNeutralino4.Array_Decays[40][1] = -PDGnustauL; ParticleNeutralino4.Array_Decays[40][2] = neutralino4amplitudenutausnutaubar; ParticleNeutralino4.Array_Decays[40][3] = 2; ParticleNeutralino4.Array_Comments[40] = "# ~chi_40 -> nu_tau ~nu_tauL*";
   ParticleNeutralino4.Array_Decays[41][0] = -PDGnutau; ParticleNeutralino4.Array_Decays[41][1] = PDGnustauL; ParticleNeutralino4.Array_Decays[41][2] = neutralino4amplitudenutaubarsnutau; ParticleNeutralino4.Array_Decays[41][3] = 2; ParticleNeutralino4.Array_Comments[41] = "# ~chi_40 -> nu_taub ~nu_tauL";
   ParticleNeutralino4.Array_Decays[42][0] = PDGWplus; ParticleNeutralino4.Array_Decays[42][1] = -PDGchargino1; ParticleNeutralino4.Array_Decays[42][2] = neutralino4amplitudeWbosonpluscharginoW1; ParticleNeutralino4.Array_Decays[42][3] = 2; ParticleNeutralino4.Array_Comments[42] = "# ~chi_40 -> W+ ~chi_1-";
   ParticleNeutralino4.Array_Decays[43][0] = PDGWplus; ParticleNeutralino4.Array_Decays[43][1] = -PDGchargino2; ParticleNeutralino4.Array_Decays[43][2] = neutralino4amplitudeWbosonpluscharginoW2; ParticleNeutralino4.Array_Decays[43][3] = 2; ParticleNeutralino4.Array_Comments[43] = "# ~chi_40 -> W+ ~chi_2-";
   ParticleNeutralino4.Array_Decays[44][0] = -PDGWplus; ParticleNeutralino4.Array_Decays[44][1] = PDGchargino1; ParticleNeutralino4.Array_Decays[44][2] = neutralino4amplitudeWbosonminuscharginoW1; ParticleNeutralino4.Array_Decays[44][3] = 2; ParticleNeutralino4.Array_Comments[44] = "# ~chi_40 -> W- ~chi_1+";
   ParticleNeutralino4.Array_Decays[45][0] = -PDGWplus; ParticleNeutralino4.Array_Decays[45][1] = PDGchargino2; ParticleNeutralino4.Array_Decays[45][2] = neutralino4amplitudeWbosonminuscharginoW2; ParticleNeutralino4.Array_Decays[45][3] = 2; ParticleNeutralino4.Array_Comments[45] = "# ~chi_40 -> W- ~chi_2+";
   ParticleNeutralino4.Array_Decays[46][0] = PDGHplus; ParticleNeutralino4.Array_Decays[46][1] = -PDGchargino1; ParticleNeutralino4.Array_Decays[46][2] = neutralino4amplitudeHpluscharginoW1; ParticleNeutralino4.Array_Decays[46][3] = 2; ParticleNeutralino4.Array_Comments[46] = "# ~chi_40 -> H+ ~chi_1-";
   ParticleNeutralino4.Array_Decays[47][0] = PDGHplus; ParticleNeutralino4.Array_Decays[47][1] = -PDGchargino2; ParticleNeutralino4.Array_Decays[47][2] = neutralino4amplitudeHpluscharginoW2; ParticleNeutralino4.Array_Decays[47][3] = 2; ParticleNeutralino4.Array_Comments[47] = "# ~chi_40 -> H+ ~chi_2-";
   ParticleNeutralino4.Array_Decays[48][0] = -PDGHplus; ParticleNeutralino4.Array_Decays[48][1] = PDGchargino1; ParticleNeutralino4.Array_Decays[48][2] = neutralino4amplitudeHminuscharginoW1; ParticleNeutralino4.Array_Decays[48][3] = 2; ParticleNeutralino4.Array_Comments[48] = "# ~chi_40 -> H- ~chi_1+";
   ParticleNeutralino4.Array_Decays[49][0] = -PDGHplus; ParticleNeutralino4.Array_Decays[49][1] = PDGchargino2; ParticleNeutralino4.Array_Decays[49][2] = neutralino4amplitudeHminuscharginoW2; ParticleNeutralino4.Array_Decays[49][3] = 2; ParticleNeutralino4.Array_Comments[49] = "# ~chi_40 -> H- ~chi_2+";
   ParticleNeutralino4.Array_Decays[50][0] = PDGZboson; ParticleNeutralino4.Array_Decays[50][1] = PDGneutralino1; ParticleNeutralino4.Array_Decays[50][2] = neutralino4amplitudeZbosonneutralino1; ParticleNeutralino4.Array_Decays[50][3] = 2; ParticleNeutralino4.Array_Comments[50] = "# ~chi_40 -> Z ~chi_10";
   ParticleNeutralino4.Array_Decays[51][0] = PDGZboson; ParticleNeutralino4.Array_Decays[51][1] = PDGneutralino2; ParticleNeutralino4.Array_Decays[51][2] = neutralino4amplitudeZbosonneutralino2; ParticleNeutralino4.Array_Decays[51][3] = 2; ParticleNeutralino4.Array_Comments[51] = "# ~chi_40 -> Z ~chi_20";
   ParticleNeutralino4.Array_Decays[52][0] = PDGZboson; ParticleNeutralino4.Array_Decays[52][1] = PDGneutralino3; ParticleNeutralino4.Array_Decays[52][2] = neutralino4amplitudeZbosonneutralino3; ParticleNeutralino4.Array_Decays[52][3] = 2; ParticleNeutralino4.Array_Comments[52] = "# ~chi_40 -> Z ~chi_30";
   ParticleNeutralino4.Array_Decays[53][0] = PDGh0; ParticleNeutralino4.Array_Decays[53][1] = PDGneutralino1; ParticleNeutralino4.Array_Decays[53][2] = neutralino4amplitudehneutralino1; ParticleNeutralino4.Array_Decays[53][3] = 2; ParticleNeutralino4.Array_Comments[53] = "# ~chi_40 -> h ~chi_10";
   ParticleNeutralino4.Array_Decays[54][0] = PDGh0; ParticleNeutralino4.Array_Decays[54][1] = PDGneutralino2; ParticleNeutralino4.Array_Decays[54][2] = neutralino4amplitudehneutralino2; ParticleNeutralino4.Array_Decays[54][3] = 2; ParticleNeutralino4.Array_Comments[54] = "# ~chi_40 -> h ~chi_20"; 
   ParticleNeutralino4.Array_Decays[55][0] = PDGh0; ParticleNeutralino4.Array_Decays[55][1] = PDGneutralino3; ParticleNeutralino4.Array_Decays[55][2] = neutralino4amplitudehneutralino3; ParticleNeutralino4.Array_Decays[55][3] = 2; ParticleNeutralino4.Array_Comments[55] = "# ~chi_40 -> h ~chi_30";
   ParticleNeutralino4.Array_Decays[56][0] = PDGH0; ParticleNeutralino4.Array_Decays[56][1] = PDGneutralino1; ParticleNeutralino4.Array_Decays[56][2] = neutralino4amplitudeHneutralino1; ParticleNeutralino4.Array_Decays[56][3] = 2; ParticleNeutralino4.Array_Comments[56] = "# ~chi_40 -> H ~chi_10";
   ParticleNeutralino4.Array_Decays[57][0] = PDGH0; ParticleNeutralino4.Array_Decays[57][1] = PDGneutralino2; ParticleNeutralino4.Array_Decays[57][2] = neutralino4amplitudeHneutralino2; ParticleNeutralino4.Array_Decays[57][3] = 2; ParticleNeutralino4.Array_Comments[57] = "# ~chi_40 -> H ~chi_20";
   ParticleNeutralino4.Array_Decays[58][0] = PDGH0; ParticleNeutralino4.Array_Decays[58][1] = PDGneutralino3; ParticleNeutralino4.Array_Decays[58][2] = neutralino4amplitudeHneutralino3; ParticleNeutralino4.Array_Decays[58][3] = 2; ParticleNeutralino4.Array_Comments[58] = "# ~chi_40 -> H ~chi_30";
   ParticleNeutralino4.Array_Decays[59][0] = PDGA0; ParticleNeutralino4.Array_Decays[59][1] = PDGneutralino1; ParticleNeutralino4.Array_Decays[59][2] = neutralino4amplitudeAneutralino1; ParticleNeutralino4.Array_Decays[59][3] = 2; ParticleNeutralino4.Array_Comments[59] = "# ~chi_40 -> A ~chi_10";
   ParticleNeutralino4.Array_Decays[60][0] = PDGA0; ParticleNeutralino4.Array_Decays[60][1] = PDGneutralino2; ParticleNeutralino4.Array_Decays[60][2] = neutralino4amplitudeAneutralino2; ParticleNeutralino4.Array_Decays[60][3] = 2; ParticleNeutralino4.Array_Comments[60] = "# ~chi_40 -> A ~chi_20";
   ParticleNeutralino4.Array_Decays[61][0] = PDGA0; ParticleNeutralino4.Array_Decays[61][1] = PDGneutralino3; ParticleNeutralino4.Array_Decays[61][2] = neutralino4amplitudeAneutralino3; ParticleNeutralino4.Array_Decays[61][3] = 2; ParticleNeutralino4.Array_Comments[61] = "# ~chi_40 -> A ~chi_30";
   
   ParticleNeutralino4.Array_Decays[62][0] = PDGH3; ParticleNeutralino4.Array_Decays[62][1] = PDGneutralino1; ParticleNeutralino4.Array_Decays[62][2] = neutralino4amplitudeH3neutralino1; ParticleNeutralino4.Array_Decays[62][3] = 2; ParticleNeutralino4.Array_Comments[62] = "# ~chi_40 -> H3 ~chi_10";
   ParticleNeutralino4.Array_Decays[63][0] = PDGH3; ParticleNeutralino4.Array_Decays[63][1] = PDGneutralino2; ParticleNeutralino4.Array_Decays[63][2] = neutralino4amplitudeH3neutralino2; ParticleNeutralino4.Array_Decays[63][3] = 2; ParticleNeutralino4.Array_Comments[63] = "# ~chi_40 -> H3 ~chi_20";
   ParticleNeutralino4.Array_Decays[64][0] = PDGH3; ParticleNeutralino4.Array_Decays[64][1] = PDGneutralino3; ParticleNeutralino4.Array_Decays[64][2] = neutralino4amplitudeH3neutralino3; ParticleNeutralino4.Array_Decays[64][3] = 2; ParticleNeutralino4.Array_Comments[64] = "# ~chi_40 -> H3 ~chi_30";
   ParticleNeutralino4.Array_Decays[65][0] = PDGA2; ParticleNeutralino4.Array_Decays[65][1] = PDGneutralino1; ParticleNeutralino4.Array_Decays[65][2] = neutralino4amplitudeA2neutralino1; ParticleNeutralino4.Array_Decays[65][3] = 2; ParticleNeutralino4.Array_Comments[65] = "# ~chi_40 -> A2 ~chi_10";
   ParticleNeutralino4.Array_Decays[66][0] = PDGA2; ParticleNeutralino4.Array_Decays[66][1] = PDGneutralino2; ParticleNeutralino4.Array_Decays[66][2] = neutralino4amplitudeA2neutralino2; ParticleNeutralino4.Array_Decays[66][3] = 2; ParticleNeutralino4.Array_Comments[66] = "# ~chi_40 -> A2 ~chi_20";
   ParticleNeutralino4.Array_Decays[67][0] = PDGA0; ParticleNeutralino4.Array_Decays[67][1] = PDGneutralino3; ParticleNeutralino4.Array_Decays[67][2] = neutralino4amplitudeA2neutralino3; ParticleNeutralino4.Array_Decays[67][3] = 2; ParticleNeutralino4.Array_Comments[67] = "# ~chi_40 -> A2 ~chi_30";
   
   
   ParticleNeutralino4.Array_Decays[68][0] = PDGphoton; ParticleNeutralino4.Array_Decays[68][1] = PDGgravitino; ParticleNeutralino4.Array_Decays[68][2] = neutralino4amplitudephotongravitino; ParticleNeutralino4.Array_Decays[68][3] = 2; ParticleNeutralino4.Array_Comments[68] = "# ~chi_40 -> gamma ~G";
   ParticleNeutralino4.Array_Decays[69][0] = PDGZboson; ParticleNeutralino4.Array_Decays[69][1] = PDGgravitino; ParticleNeutralino4.Array_Decays[69][2] = neutralino4amplitudeZgravitino; ParticleNeutralino4.Array_Decays[69][3] = 2; ParticleNeutralino4.Array_Comments[69] = "# ~chi_40 -> Z ~G";
   ParticleNeutralino4.Array_Decays[70][0] = PDGh0; ParticleNeutralino4.Array_Decays[70][1] = PDGgravitino; ParticleNeutralino4.Array_Decays[70][2] = neutralino4amplitudehgravitino; ParticleNeutralino4.Array_Decays[70][3] = 2; ParticleNeutralino4.Array_Comments[70] = "# ~chi_40 -> h ~G";
   ParticleNeutralino4.Array_Decays[71][0] = PDGH0; ParticleNeutralino4.Array_Decays[71][1] = PDGgravitino; ParticleNeutralino4.Array_Decays[71][2] = neutralino4amplitudeHgravitino; ParticleNeutralino4.Array_Decays[71][3] = 2; ParticleNeutralino4.Array_Comments[71] = "# ~chi_40 -> H ~G";
   ParticleNeutralino4.Array_Decays[72][0] = PDGA0; ParticleNeutralino4.Array_Decays[72][1] = PDGgravitino; ParticleNeutralino4.Array_Decays[72][2] = neutralino4amplitudeAgravitino; ParticleNeutralino4.Array_Decays[72][3] = 2; ParticleNeutralino4.Array_Comments[72] = "# ~chi_40 -> A ~G";
   
   ParticleNeutralino4.Array_Decays[73][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[73][1] = PDGup; ParticleNeutralino4.Array_Decays[73][4] = -PDGup; ParticleNeutralino4.Array_Decays[73][2] = neutralino4amplitudeneut1uubar; ParticleNeutralino4.Array_Decays[73][3] = 3; ParticleNeutralino4.Array_Comments[73] = "# ~chi_40 -> ~chi_10 u ubar";
   ParticleNeutralino4.Array_Decays[74][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[74][1] = PDGdown; ParticleNeutralino4.Array_Decays[74][4] = -PDGdown; ParticleNeutralino4.Array_Decays[74][2] = neutralino4amplitudeneut1ddbar; ParticleNeutralino4.Array_Decays[74][3] = 3; ParticleNeutralino4.Array_Comments[74] = "# ~chi_40 -> ~chi_10 d dbar";
   ParticleNeutralino4.Array_Decays[75][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[75][1] = PDGcharm; ParticleNeutralino4.Array_Decays[75][4] = -PDGcharm; ParticleNeutralino4.Array_Decays[75][2] = neutralino4amplitudeneut1ccbar; ParticleNeutralino4.Array_Decays[75][3] = 3; ParticleNeutralino4.Array_Comments[75] = "# ~chi_40 -> ~chi_10 c cbar";
   ParticleNeutralino4.Array_Decays[76][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[76][1] = PDGstrange; ParticleNeutralino4.Array_Decays[76][4] = -PDGstrange; ParticleNeutralino4.Array_Decays[76][2] = neutralino4amplitudeneut1ssbar; ParticleNeutralino4.Array_Decays[76][3] = 3; ParticleNeutralino4.Array_Comments[76] = "# ~chi_40 -> ~chi_10 s sbar";
   ParticleNeutralino4.Array_Decays[77][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[77][1] = PDGtop; ParticleNeutralino4.Array_Decays[77][4] = -PDGtop; ParticleNeutralino4.Array_Decays[77][2] = neutralino4amplitudeneut1ttbar; ParticleNeutralino4.Array_Decays[77][3] = 3; ParticleNeutralino4.Array_Comments[77] = "# ~chi_40 -> ~chi_10 t tbar";
   ParticleNeutralino4.Array_Decays[78][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[78][1] = PDGbottom; ParticleNeutralino4.Array_Decays[78][4] = -PDGbottom; ParticleNeutralino4.Array_Decays[78][2] = neutralino4amplitudeneut1bbbar; ParticleNeutralino4.Array_Decays[78][3] = 3; ParticleNeutralino4.Array_Comments[78] = "# ~chi_40 -> ~chi_10 b bbar";
   ParticleNeutralino4.Array_Decays[79][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[79][1] = PDGelectron; ParticleNeutralino4.Array_Decays[79][4] = -PDGelectron; ParticleNeutralino4.Array_Decays[79][2] = neutralino4amplitudeneut1eebar; ParticleNeutralino4.Array_Decays[79][3] = 3; ParticleNeutralino4.Array_Comments[79] = "# ~chi_40 -> ~chi_10 e- e+";
   ParticleNeutralino4.Array_Decays[80][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[80][1] = PDGmuon; ParticleNeutralino4.Array_Decays[80][4] = -PDGmuon; ParticleNeutralino4.Array_Decays[80][2] = neutralino4amplitudeneut1mumubar; ParticleNeutralino4.Array_Decays[80][3] = 3; ParticleNeutralino4.Array_Comments[80] = "# ~chi_40 -> ~chi_10 mu- mu+";
   ParticleNeutralino4.Array_Decays[81][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[81][1] = PDGtau; ParticleNeutralino4.Array_Decays[81][4] = -PDGtau; ParticleNeutralino4.Array_Decays[81][2] = neutralino4amplitudeneut1tautaubar; ParticleNeutralino4.Array_Decays[81][3] = 3; ParticleNeutralino4.Array_Comments[81] = "# ~chi_40 -> ~chi_10 tau- tau+";
   ParticleNeutralino4.Array_Decays[82][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[82][1] = PDGnuelectron; ParticleNeutralino4.Array_Decays[82][4] = -PDGnuelectron; ParticleNeutralino4.Array_Decays[82][2] = neutralino4amplitudeneut1nuenuebar; ParticleNeutralino4.Array_Decays[82][3] = 3; ParticleNeutralino4.Array_Comments[82] = "# ~chi_40 -> ~chi_10 nue nuebar";
   ParticleNeutralino4.Array_Decays[83][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[83][1] = PDGnumuon; ParticleNeutralino4.Array_Decays[83][4] = -PDGnumuon; ParticleNeutralino4.Array_Decays[83][2] = neutralino4amplitudeneut1numunumubar; ParticleNeutralino4.Array_Decays[83][3] = 3; ParticleNeutralino4.Array_Comments[83] = "# ~chi_40 -> ~chi_10 numu numubar";
   ParticleNeutralino4.Array_Decays[84][0] = PDGneutralino1; ParticleNeutralino4.Array_Decays[84][1] = PDGnutau; ParticleNeutralino4.Array_Decays[84][4] = -PDGnutau; ParticleNeutralino4.Array_Decays[84][2] = neutralino4amplitudeneut1nutaunutaubar; ParticleNeutralino4.Array_Decays[84][3] = 3; ParticleNeutralino4.Array_Comments[84] = "# ~chi_40 -> ~chi_10 nutau nutaubar";
   ParticleNeutralino4.Array_Decays[85][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[85][1] = PDGup; ParticleNeutralino4.Array_Decays[85][4] = -PDGup; ParticleNeutralino4.Array_Decays[85][2] = neutralino4amplitudeneut2uubar; ParticleNeutralino4.Array_Decays[85][3] = 3; ParticleNeutralino4.Array_Comments[85] = "# ~chi_40 -> ~chi_20 u ubar";
   ParticleNeutralino4.Array_Decays[86][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[86][1] = PDGdown; ParticleNeutralino4.Array_Decays[86][4] = -PDGdown; ParticleNeutralino4.Array_Decays[86][2] = neutralino4amplitudeneut2ddbar; ParticleNeutralino4.Array_Decays[86][3] = 3; ParticleNeutralino4.Array_Comments[86] = "# ~chi_40 -> ~chi_20 d dbar";
   ParticleNeutralino4.Array_Decays[87][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[87][1] = PDGcharm; ParticleNeutralino4.Array_Decays[87][4] = -PDGcharm; ParticleNeutralino4.Array_Decays[87][2] = neutralino4amplitudeneut2ccbar; ParticleNeutralino4.Array_Decays[87][3] = 3; ParticleNeutralino4.Array_Comments[87] = "# ~chi_40 -> ~chi_20 c cbar";
   ParticleNeutralino4.Array_Decays[88][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[88][1] = PDGstrange; ParticleNeutralino4.Array_Decays[88][4] = -PDGstrange; ParticleNeutralino4.Array_Decays[88][2] = neutralino4amplitudeneut2ssbar; ParticleNeutralino4.Array_Decays[88][3] = 3; ParticleNeutralino4.Array_Comments[88] = "# ~chi_40 -> ~chi_20 s sbar";
   ParticleNeutralino4.Array_Decays[89][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[89][1] = PDGtop; ParticleNeutralino4.Array_Decays[89][4] = -PDGtop; ParticleNeutralino4.Array_Decays[89][2] = neutralino4amplitudeneut2ttbar; ParticleNeutralino4.Array_Decays[89][3] = 3; ParticleNeutralino4.Array_Comments[89] = "# ~chi_40 -> ~chi_20 t tbar";
   ParticleNeutralino4.Array_Decays[90][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[90][1] = PDGbottom; ParticleNeutralino4.Array_Decays[90][4] = -PDGbottom; ParticleNeutralino4.Array_Decays[90][2] = neutralino4amplitudeneut2bbbar; ParticleNeutralino4.Array_Decays[90][3] = 3; ParticleNeutralino4.Array_Comments[90] = "# ~chi_40 -> ~chi_20 b bbar";
   ParticleNeutralino4.Array_Decays[91][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[91][1] = PDGelectron; ParticleNeutralino4.Array_Decays[91][4] = -PDGelectron; ParticleNeutralino4.Array_Decays[91][2] = neutralino4amplitudeneut2eebar; ParticleNeutralino4.Array_Decays[91][3] = 3; ParticleNeutralino4.Array_Comments[91] = "# ~chi_40 -> ~chi_20 e- e+";
   ParticleNeutralino4.Array_Decays[92][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[92][1] = PDGmuon; ParticleNeutralino4.Array_Decays[92][4] = -PDGmuon; ParticleNeutralino4.Array_Decays[92][2] = neutralino4amplitudeneut2mumubar; ParticleNeutralino4.Array_Decays[92][3] = 3; ParticleNeutralino4.Array_Comments[92] = "# ~chi_40 -> ~chi_20 mu- mu+";
   ParticleNeutralino4.Array_Decays[93][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[93][1] = PDGtau; ParticleNeutralino4.Array_Decays[93][4] = -PDGtau; ParticleNeutralino4.Array_Decays[93][2] = neutralino4amplitudeneut2tautaubar; ParticleNeutralino4.Array_Decays[93][3] = 3; ParticleNeutralino4.Array_Comments[93] = "# ~chi_40 -> ~chi_20 tau- tau+";
   ParticleNeutralino4.Array_Decays[94][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[94][1] = PDGnuelectron; ParticleNeutralino4.Array_Decays[94][4] = -PDGnuelectron; ParticleNeutralino4.Array_Decays[94][2] = neutralino4amplitudeneut2nuenuebar; ParticleNeutralino4.Array_Decays[94][3] = 3; ParticleNeutralino4.Array_Comments[94] = "# ~chi_40 -> ~chi_20 nue nuebar";
   ParticleNeutralino4.Array_Decays[95][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[95][1] = PDGnumuon; ParticleNeutralino4.Array_Decays[95][4] = -PDGnumuon; ParticleNeutralino4.Array_Decays[95][2] = neutralino4amplitudeneut2numunumubar; ParticleNeutralino4.Array_Decays[95][3] = 3; ParticleNeutralino4.Array_Comments[95] = "# ~chi_40 -> ~chi_20 numu numubar";
   ParticleNeutralino4.Array_Decays[96][0] = PDGneutralino2; ParticleNeutralino4.Array_Decays[96][1] = PDGnutau; ParticleNeutralino4.Array_Decays[96][4] = -PDGnutau; ParticleNeutralino4.Array_Decays[96][2] = neutralino4amplitudeneut2nutaunutaubar; ParticleNeutralino4.Array_Decays[96][3] = 3; ParticleNeutralino4.Array_Comments[96] = "# ~chi_40 -> ~chi_20 nutau nutaubar";
   ParticleNeutralino4.Array_Decays[97][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[97][1] = PDGup; ParticleNeutralino4.Array_Decays[97][4] = -PDGup; ParticleNeutralino4.Array_Decays[97][2] = neutralino4amplitudeneut3uubar; ParticleNeutralino4.Array_Decays[97][3] = 3; ParticleNeutralino4.Array_Comments[97] = "# ~chi_40 -> ~chi_30 u ubar";
   ParticleNeutralino4.Array_Decays[98][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[98][1] = PDGdown; ParticleNeutralino4.Array_Decays[98][4] = -PDGdown; ParticleNeutralino4.Array_Decays[98][2] = neutralino4amplitudeneut3ddbar; ParticleNeutralino4.Array_Decays[98][3] = 3; ParticleNeutralino4.Array_Comments[98] = "# ~chi_40 -> ~chi_30 d dbar";
   ParticleNeutralino4.Array_Decays[99][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[99][1] = PDGcharm; ParticleNeutralino4.Array_Decays[99][4] = -PDGcharm; ParticleNeutralino4.Array_Decays[99][2] = neutralino4amplitudeneut3ccbar; ParticleNeutralino4.Array_Decays[99][3] = 3; ParticleNeutralino4.Array_Comments[99] = "# ~chi_40 -> ~chi_30 c cbar";
   ParticleNeutralino4.Array_Decays[100][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[100][1] = PDGstrange; ParticleNeutralino4.Array_Decays[100][4] = -PDGstrange; ParticleNeutralino4.Array_Decays[100][2] = neutralino4amplitudeneut3ssbar; ParticleNeutralino4.Array_Decays[100][3] = 3; ParticleNeutralino4.Array_Comments[100] = "# ~chi_40 -> ~chi_30 s sbar";
   ParticleNeutralino4.Array_Decays[101][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[101][1] = PDGtop; ParticleNeutralino4.Array_Decays[101][4] = -PDGtop; ParticleNeutralino4.Array_Decays[101][2] = neutralino4amplitudeneut3ttbar; ParticleNeutralino4.Array_Decays[101][3] = 3; ParticleNeutralino4.Array_Comments[101] = "# ~chi_40 -> ~chi_30 t tbar";
   ParticleNeutralino4.Array_Decays[102][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[102][1] = PDGbottom; ParticleNeutralino4.Array_Decays[102][4] = -PDGbottom; ParticleNeutralino4.Array_Decays[102][2] = neutralino4amplitudeneut3bbbar; ParticleNeutralino4.Array_Decays[102][3] = 3; ParticleNeutralino4.Array_Comments[102] = "# ~chi_40 -> ~chi_30 b bbar";
   ParticleNeutralino4.Array_Decays[103][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[103][1] = PDGelectron; ParticleNeutralino4.Array_Decays[103][4] = -PDGelectron; ParticleNeutralino4.Array_Decays[103][2] = neutralino4amplitudeneut3eebar; ParticleNeutralino4.Array_Decays[103][3] = 3; ParticleNeutralino4.Array_Comments[103] = "# ~chi_40 -> ~chi_30 e- e+";
   ParticleNeutralino4.Array_Decays[104][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[104][1] = PDGmuon; ParticleNeutralino4.Array_Decays[104][4] = -PDGmuon; ParticleNeutralino4.Array_Decays[104][2] = neutralino4amplitudeneut3mumubar; ParticleNeutralino4.Array_Decays[104][3] = 3; ParticleNeutralino4.Array_Comments[104] = "# ~chi_40 -> ~chi_30 mu- mu+";
   ParticleNeutralino4.Array_Decays[105][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[105][1] = PDGtau; ParticleNeutralino4.Array_Decays[105][4] = -PDGtau; ParticleNeutralino4.Array_Decays[105][2] = neutralino4amplitudeneut3tautaubar; ParticleNeutralino4.Array_Decays[105][3] = 3; ParticleNeutralino4.Array_Comments[105] = "# ~chi_40 -> ~chi_30 tau- tau+";
   ParticleNeutralino4.Array_Decays[106][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[106][1] = PDGnuelectron; ParticleNeutralino4.Array_Decays[106][4] = -PDGnuelectron; ParticleNeutralino4.Array_Decays[106][2] = neutralino4amplitudeneut3nuenuebar; ParticleNeutralino4.Array_Decays[106][3] = 3; ParticleNeutralino4.Array_Comments[106] = "# ~chi_40 -> ~chi_30 nue nuebar";
   ParticleNeutralino4.Array_Decays[107][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[107][1] = PDGnumuon; ParticleNeutralino4.Array_Decays[107][4] = -PDGnumuon; ParticleNeutralino4.Array_Decays[107][2] = neutralino4amplitudeneut3numunumubar; ParticleNeutralino4.Array_Decays[107][3] = 3; ParticleNeutralino4.Array_Comments[107] = "# ~chi_40 -> ~chi_30 numu numubar";
   ParticleNeutralino4.Array_Decays[108][0] = PDGneutralino3; ParticleNeutralino4.Array_Decays[108][1] = PDGnutau; ParticleNeutralino4.Array_Decays[108][4] = -PDGnutau; ParticleNeutralino4.Array_Decays[108][2] = neutralino4amplitudeneut3nutaunutaubar; ParticleNeutralino4.Array_Decays[108][3] = 3; ParticleNeutralino4.Array_Comments[108] = "# ~chi_40 -> ~chi_30 nutau nutaubar";
   
   ParticleNeutralino4.Array_Decays[109][0] = PDGchargino1; ParticleNeutralino4.Array_Decays[109][1] = PDGup; ParticleNeutralino4.Array_Decays[109][4] = -PDGdown; ParticleNeutralino4.Array_Decays[109][2] = neutralino4amplitudechargino1udbar; ParticleNeutralino4.Array_Decays[109][3] = 3; ParticleNeutralino4.Array_Comments[109] = "# ~chi_40 -> chi_1- u db";
   ParticleNeutralino4.Array_Decays[110][0] = PDGchargino1; ParticleNeutralino4.Array_Decays[110][1] = PDGcharm; ParticleNeutralino4.Array_Decays[110][4] = -PDGstrange; ParticleNeutralino4.Array_Decays[110][2] = neutralino4amplitudechargino1csbar; ParticleNeutralino4.Array_Decays[110][3] = 3; ParticleNeutralino4.Array_Comments[110] = "# ~chi_40 -> chi_1- c sb";
   ParticleNeutralino4.Array_Decays[111][0] = PDGchargino1; ParticleNeutralino4.Array_Decays[111][1] = PDGnuelectron; ParticleNeutralino4.Array_Decays[111][4] = -PDGelectron; ParticleNeutralino4.Array_Decays[111][2] = neutralino4amplitudechargino1enuebar; ParticleNeutralino4.Array_Decays[111][3] = 3; ParticleNeutralino4.Array_Comments[111] = "# ~chi_40 -> chi_1- nu_e eb";
   ParticleNeutralino4.Array_Decays[112][0] = PDGchargino1; ParticleNeutralino4.Array_Decays[112][1] = PDGnumuon; ParticleNeutralino4.Array_Decays[112][4] = -PDGmuon; ParticleNeutralino4.Array_Decays[112][2] = neutralino4amplitudechargino1munumubar; ParticleNeutralino4.Array_Decays[112][3] = 3; ParticleNeutralino4.Array_Comments[112] = "# ~chi_40 -> chi_1- nu_mu mub";
   ParticleNeutralino4.Array_Decays[113][0] = PDGchargino1; ParticleNeutralino4.Array_Decays[113][1] = PDGnutau; ParticleNeutralino4.Array_Decays[113][4] = -PDGtau; ParticleNeutralino4.Array_Decays[113][2] = neutralino4amplitudechargino1taunutaubar; ParticleNeutralino4.Array_Decays[113][3] = 3; ParticleNeutralino4.Array_Comments[113] = "# ~chi_40 -> chi_1- nu_tau taub";
   ParticleNeutralino4.Array_Decays[114][0] = PDGchargino2; ParticleNeutralino4.Array_Decays[114][1] = PDGup; ParticleNeutralino4.Array_Decays[114][4] = -PDGdown; ParticleNeutralino4.Array_Decays[114][2] = neutralino4amplitudechargino2udbar; ParticleNeutralino4.Array_Decays[114][3] = 3; ParticleNeutralino4.Array_Comments[114] = "# ~chi_40 -> chi_2- u dbar";
   ParticleNeutralino4.Array_Decays[115][0] = PDGchargino2; ParticleNeutralino4.Array_Decays[115][1] = PDGcharm; ParticleNeutralino4.Array_Decays[115][4] = -PDGstrange; ParticleNeutralino4.Array_Decays[115][2] = neutralino4amplitudechargino2csbar; ParticleNeutralino4.Array_Decays[115][3] = 3; ParticleNeutralino4.Array_Comments[115] = "# ~chi_40 -> chi_2- c sbar";
   ParticleNeutralino4.Array_Decays[116][0] = PDGchargino2; ParticleNeutralino4.Array_Decays[116][1] = PDGnuelectron; ParticleNeutralino4.Array_Decays[116][4] = -PDGelectron; ParticleNeutralino4.Array_Decays[116][2] = neutralino4amplitudechargino2enuebar; ParticleNeutralino4.Array_Decays[116][3] = 3; ParticleNeutralino4.Array_Comments[116] = "# ~chi_40 -> chi_2- nu_e eb";
   ParticleNeutralino4.Array_Decays[117][0] = PDGchargino2; ParticleNeutralino4.Array_Decays[117][1] = PDGnumuon; ParticleNeutralino4.Array_Decays[117][4] = -PDGmuon; ParticleNeutralino4.Array_Decays[117][2] = neutralino4amplitudechargino2munumubar; ParticleNeutralino4.Array_Decays[117][3] = 3; ParticleNeutralino4.Array_Comments[117] = "# ~chi_40 -> chi_2- nu_mu mub";
   ParticleNeutralino4.Array_Decays[118][0] = PDGchargino2; ParticleNeutralino4.Array_Decays[118][1] = PDGnutau; ParticleNeutralino4.Array_Decays[118][4] = -PDGtau; ParticleNeutralino4.Array_Decays[118][2] = neutralino4amplitudechargino2taunutaubar; ParticleNeutralino4.Array_Decays[118][3] = 3; ParticleNeutralino4.Array_Comments[118] = "# ~chi_40 -> chi_2- nu_tau taubar";

   for(int i = 0; i<ParticleNeutralino4.No_of_Decays; i++) {
     if (ParticleNeutralino4.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleNeutralino4.Array_Comments[i] << " is negative = " << ParticleNeutralino4.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleNeutralino4.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }
   
   double Neut4_No_1to2_Decays = 0;
   
   Neut4_No_1to2_Decays = ParticleNeutralino4.No_1to2_Decays + ParticleNeutralino4.No_grav_Decays  + ParticleNeutralino4.No_NMSSM_Decays;
   
   for (int j = 0; j<Neut4_No_1to2_Decays; j++) {
     ParticleNeutralino4.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<Neut4_No_1to2_Decays; j++) {
     ParticleNeutralino4.two_width = ParticleNeutralino4.two_width + ParticleNeutralino4.Array_Decays[j][2];
   }
   for (int j=Neut4_No_1to2_Decays; j<ParticleNeutralino4.No_of_Decays; j++) {
     ParticleNeutralino4.three_width = ParticleNeutralino4.three_width + ParticleNeutralino4.Array_Decays[j][2];
   }
   
   if ( ParticleNeutralino4.three_width != ParticleNeutralino4.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for neutralino 4 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleNeutralino4.No_of_Decays = Neut4_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleNeutralino4.total_width = ParticleNeutralino4.two_width;
     }
   else {
     ParticleNeutralino4.total_width = ParticleNeutralino4.two_width + ParticleNeutralino4.three_width;
   }
   
   if ( ParticleNeutralino4.total_width != ParticleNeutralino4.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleNeutralino4.No_of_Decays; i++) {
       //   fout << i << " " << ParticleNeutralino4.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in Neutralino4 total width \n");
     }
      
 }
 
 if (nmssmIsIt == true) {
   ///Neutralino5 Decays
   double neutralino5amplitudeuLubar=0, neutralino5amplitudeuRubar=0, neutralino5amplitudeuLbaru=0, neutralino5amplitudeuRbaru=0, neutralino5amplitudedLdbar=0, neutralino5amplitudedRdbar=0, neutralino5amplitudedLbard=0, neutralino5amplitudedRbard=0, neutralino5amplitudecLcbar=0, neutralino5amplitudecRcbar=0, neutralino5amplitudecLbarc=0, neutralino5amplitudecRbarc=0, neutralino5amplitudesLsbar=0, neutralino5amplitudesRsbar=0, neutralino5amplitudesLbars=0, neutralino5amplitudesRbars=0, neutralino5amplitudeeLebar=0, neutralino5amplitudeeRebar=0, neutralino5amplitudeeLbare=0, neutralino5amplitudeeRbare=0, neutralino5amplitudemuLmubar=0, neutralino5amplitudemuRmubar=0, neutralino5amplitudemuLbarmu=0, neutralino5amplitudemuRbarmu=0,neutralino5amplitudesnuenuebar=0, neutralino5amplitudesnuebarnue=0, neutralino5amplitudesnumunumubar=0, neutralino5amplitudesnumubarnumu=0, neutralino5amplitudetopstop1bar=0, neutralino5amplitudetopstop2bar=0, neutralino5amplitudetopbarstop1=0, neutralino5amplitudetopbarstop2=0, neutralino5amplitudebottomsbottom1bar=0, neutralino5amplitudebottomsbottom2bar=0, neutralino5amplitudebottombarsbottom1=0, neutralino5amplitudebottombarsbottom2=0, neutralino5amplitudetaustau1bar=0, neutralino5amplitudetaustau2bar=0, neutralino5amplitudetaubarstau1=0, neutralino5amplitudetaubarstau2=0, neutralino5amplitudenutausnutaubar=0, neutralino5amplitudenutaubarsnutau=0, neutralino5amplitudeWbosonpluscharginoW1=0, neutralino5amplitudeWbosonpluscharginoW2=0, neutralino5amplitudeWbosonminuscharginoW1=0, neutralino5amplitudeWbosonminuscharginoW2=0, neutralino5amplitudeHpluscharginoW1=0, neutralino5amplitudeHpluscharginoW2=0, neutralino5amplitudeHminuscharginoW1=0, neutralino5amplitudeHminuscharginoW2=0, neutralino5amplitudeZbosonneutralino1=0, neutralino5amplitudeZbosonneutralino2=0, neutralino5amplitudeZbosonneutralino3=0, neutralino5amplitudeZbosonneutralino4=0, neutralino5amplitudehneutralino1=0, neutralino5amplitudehneutralino2=0, neutralino5amplitudehneutralino3=0, neutralino5amplitudehneutralino4=0, neutralino5amplitudeHneutralino1=0, neutralino5amplitudeHneutralino2=0, neutralino5amplitudeHneutralino3=0, neutralino5amplitudeHneutralino4=0, neutralino5amplitudeAneutralino1=0, neutralino5amplitudeAneutralino2=0, neutralino5amplitudeAneutralino3=0, neutralino5amplitudeAneutralino4=0, neutralino5amplitudeH3neutralino1 = 0, neutralino5amplitudeH3neutralino2 = 0, neutralino5amplitudeH3neutralino3 = 0, neutralino5amplitudeH3neutralino4 = 0, neutralino5amplitudeA2neutralino1=0, neutralino5amplitudeA2neutralino2=0, neutralino5amplitudeA2neutralino3=0, neutralino5amplitudeA2neutralino4=0;
   
   if (flagneut5 == 1) {
     neutralino5amplitudeWbosonpluscharginoW1 = neutralinoamplitudedecaycharginoWNMSSM (mneut(5), MCH1, polemw, g, thetaL2, thetaR2, mixNeut, 5, 1);
     neutralino5amplitudeWbosonpluscharginoW2 = neutralinoamplitudedecaycharginoWNMSSM (mneut(5), MCH2, polemw, g, thetaL2, thetaR2, mixNeut, 5, 2);
     neutralino5amplitudeWbosonminuscharginoW1 = neutralino5amplitudeWbosonpluscharginoW1;
     neutralino5amplitudeWbosonminuscharginoW2 = neutralino5amplitudeWbosonpluscharginoW2;
     neutralino5amplitudeZbosonneutralino1 = neutralinoamplitudedecayneutralinoZNMSSM (mneut(5), mneut(1), polemz, g, gp, mixNeut, 5, 1);
     neutralino5amplitudeZbosonneutralino2 = neutralinoamplitudedecayneutralinoZNMSSM (mneut(5), mneut(2), polemz, g, gp, mixNeut, 5, 2);
     neutralino5amplitudeZbosonneutralino3 = neutralinoamplitudedecayneutralinoZNMSSM (mneut(5), mneut(3), polemz, g, gp, mixNeut, 5, 3);
     neutralino5amplitudeZbosonneutralino4 = neutralinoamplitudedecayneutralinoZNMSSM (mneut(5), mneut(4), polemz, g, gp, mixNeut, 5, 4);
     
     neutralino5amplitudeHpluscharginoW1 = neutralinoamplitudecharginoHpmNMSSM (mneut(5), MCH1, mHpm, g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 5, 1);
     neutralino5amplitudeHpluscharginoW2 = neutralinoamplitudecharginoHpmNMSSM (mneut(5), MCH2, mHpm, g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 5, 2);
     neutralino5amplitudeHminuscharginoW1 = neutralino5amplitudeHpluscharginoW1;
     neutralino5amplitudeHminuscharginoW2 = neutralino5amplitudeHpluscharginoW2;
     
     neutralino5amplitudehneutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(1), mh0(1), g, gp, lam, kappa, mixNeut, CPEMix, 5, 1, 1);
     neutralino5amplitudeHneutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(1), mh0(2), g, gp, lam, kappa, mixNeut, CPEMix, 5, 1, 2);
     neutralino5amplitudeH3neutralino1 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(1), mh0(3), g, gp, lam, kappa, mixNeut, CPEMix, 5, 1, 3);
     neutralino5amplitudehneutralino2 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(2), mh0(1), g, gp, lam, kappa, mixNeut, CPEMix, 5, 2, 1);
     neutralino5amplitudeHneutralino2 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(2), mh0(2), g, gp, lam, kappa, mixNeut, CPEMix, 5, 2, 2);
     neutralino5amplitudeH3neutralino2 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(2), mh0(3), g, gp, lam, kappa, mixNeut, CPEMix, 5, 2, 3);
     neutralino5amplitudehneutralino3 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(3), mh0(1), g, gp, lam, kappa, mixNeut, CPEMix, 5, 3, 1);
     neutralino5amplitudeHneutralino3 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(3), mh0(2), g, gp, lam, kappa, mixNeut, CPEMix, 5, 3, 2);
     neutralino5amplitudeH3neutralino3 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(3), mh0(3), g, gp, lam, kappa, mixNeut, CPEMix, 5, 3, 3);
     neutralino5amplitudehneutralino4 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(4), mh0(1), g, gp, lam, kappa, mixNeut, CPEMix, 5, 4, 1);
     neutralino5amplitudeHneutralino4 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(4), mh0(2), g, gp, lam, kappa, mixNeut, CPEMix, 5, 4, 2);
     neutralino5amplitudeH3neutralino4 = neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (mneut(5), mneut(4), mh0(3), g, gp, lam, kappa, mixNeut, CPEMix, 5, 4, 3);
     
     neutralino5amplitudeAneutralino1 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(5), mneut(1), mA0(1), g, gp, lam, kappa, mixNeut, CPOMix, 5, 1, 1);
     neutralino5amplitudeA2neutralino1 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(5), mneut(1), mA0(2), g, gp, lam, kappa, mixNeut, CPOMix, 5, 1, 2);
     neutralino5amplitudeAneutralino2 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(5), mneut(2), mA0(1), g, gp, lam, kappa, mixNeut, CPOMix, 5, 2, 1);
     neutralino5amplitudeA2neutralino2 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(5), mneut(2), mA0(2), g, gp, lam, kappa, mixNeut, CPOMix, 5, 2, 2);
     neutralino5amplitudeAneutralino3 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(5), mneut(3), mA0(1), g, gp, lam, kappa, mixNeut, CPOMix, 5, 3, 1);
     neutralino5amplitudeA2neutralino3 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(5), mneut(3), mA0(2), g, gp, lam, kappa, mixNeut, CPOMix, 5, 3, 2);
     neutralino5amplitudeAneutralino4 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(5), mneut(4), mA0(1), g, gp, lam, kappa, mixNeut, CPOMix, 5, 4, 1);
     neutralino5amplitudeA2neutralino4 = neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (mneut(5), mneut(4), mA0(2), g, gp, lam, kappa, mixNeut, CPOMix, 5, 4, 2);
     
     neutralino5amplitudeuLubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), mu(1,1), mup, g, gp, mixNeut, 5, 'u', 'L');
     neutralino5amplitudeuRubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), mu(2,1), mup, g, gp, mixNeut, 5, 'u', 'R');
     neutralino5amplitudeuLbaru = neutralino5amplitudeuLubar;
     neutralino5amplitudeuRbaru = neutralino5amplitudeuRubar;
     neutralino5amplitudedLdbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), md(1,1), mdo, g, gp, mixNeut, 5, 'd', 'L');
     neutralino5amplitudedRdbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), md(2,1), mdo, g, gp, mixNeut, 5, 'd', 'R');
     neutralino5amplitudedLbard = neutralino5amplitudedLdbar;
     neutralino5amplitudedRbard = neutralino5amplitudedRdbar;
     neutralino5amplitudecLcbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), mu(1,2), mc, g, gp, mixNeut, 5, 'u', 'L');
     neutralino5amplitudecRcbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), mu(2,2), mc, g, gp, mixNeut, 5, 'u', 'R');
     neutralino5amplitudecLbarc = neutralino5amplitudecLcbar;
     neutralino5amplitudecRbarc = neutralino5amplitudecRcbar;
     neutralino5amplitudesLsbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), md(1,2), ms, g, gp, mixNeut, 5, 'd', 'L');
     neutralino5amplitudesRsbar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), md(2,2), ms, g, gp, mixNeut, 5, 'd', 'R');
     neutralino5amplitudesLbars = neutralino5amplitudesLsbar;
     neutralino5amplitudesRbars = neutralino5amplitudesRsbar;
     neutralino5amplitudeeLebar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), me(1,1), mel, g, gp, mixNeut, 5, 'l', 'L');
     neutralino5amplitudeeRebar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), me(2,1), mel, g, gp, mixNeut, 5, 'l', 'R');
     neutralino5amplitudeeLbare = neutralino5amplitudeeLebar;
     neutralino5amplitudeeRbare = neutralino5amplitudeeRebar;
     neutralino5amplitudemuLmubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), me(1,2), mmu, g, gp, mixNeut, 5, 'l', 'L');
     neutralino5amplitudemuRmubar = neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (mneut(5), me(2,2), mmu, g, gp, mixNeut, 5, 'l', 'R');
     neutralino5amplitudemuLbarmu = neutralino5amplitudemuLmubar;
     neutralino5amplitudemuRbarmu = neutralino5amplitudemuRmubar;
     
     neutralino5amplitudetopstop1bar = neutralinoamplitudestoptopNMSSM (mneut(5), mu(1,3), mt, g, gp, thetat, beta, runmw, mixNeut, runmt, 5, 1);
     neutralino5amplitudetopstop2bar = neutralinoamplitudestoptopNMSSM (mneut(5), mu(2,3), mt, g, gp, thetat, beta, runmw, mixNeut, runmt, 5, 2);
     neutralino5amplitudetopbarstop1 = neutralino5amplitudetopstop1bar;
     neutralino5amplitudetopbarstop2 = neutralino5amplitudetopstop2bar;
     neutralino5amplitudebottomsbottom1bar = neutralinoamplitudesbottombottomNMSSM (mneut(5), md(1,3), mb, g, gp, thetab, beta, runmw, mixNeut, runmb, 5, 1);
     neutralino5amplitudebottomsbottom2bar = neutralinoamplitudesbottombottomNMSSM (mneut(5), md(2,3), mb, g, gp, thetab, beta, runmw, mixNeut, runmb, 5, 2);
     neutralino5amplitudebottombarsbottom1 = neutralino5amplitudebottomsbottom1bar;
     neutralino5amplitudebottombarsbottom2 = neutralino5amplitudebottomsbottom2bar;
     
     neutralino5amplitudetaustau1bar = neutralinoamplitudestautauNMSSM (mneut(5), me(1,3), mtau, g, gp, thetatau, beta, runmw, mixNeut, runmtau, 5, 1);
     neutralino5amplitudetaustau2bar = neutralinoamplitudestautauNMSSM (mneut(5), me(2,3), mtau, g, gp, thetatau, beta, runmw, mixNeut, runmtau, 5, 2);
     neutralino5amplitudetaubarstau1 = neutralino5amplitudetaustau1bar;
     neutralino5amplitudetaubarstau2 = neutralino5amplitudetaustau2bar;
     
     neutralino5amplitudesnuenuebar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(5), msnu(1), 0, g, gp, mixNeut, 5);
     neutralino5amplitudesnuebarnue = neutralino5amplitudesnuenuebar;
     neutralino5amplitudesnumunumubar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(5), msnu(2), 0, g, gp, mixNeut, 5);
     neutralino5amplitudesnumubarnumu = neutralino5amplitudesnumunumubar;
     neutralino5amplitudenutausnutaubar = neutralinoamplitudestauneutrinotauneutrinoNMSSM (mneut(5), msnu(3), 0, g, gp, mixNeut, 5);
     neutralino5amplitudenutaubarsnutau = neutralino5amplitudenutausnutaubar;
     
     ParticleNeutralino5.Array_Decays[0][0] = -PDGup; ParticleNeutralino5.Array_Decays[0][1] = PDGsupL; ParticleNeutralino5.Array_Decays[0][2] = neutralino5amplitudeuLubar; ParticleNeutralino5.Array_Decays[0][3] = 2; ParticleNeutralino5.Array_Comments[0] = "# ~chi_50 -> ub ~u_L";
     ParticleNeutralino5.Array_Decays[1][0] = -PDGup; ParticleNeutralino5.Array_Decays[1][1] = PDGsupR; ParticleNeutralino5.Array_Decays[1][2] = neutralino5amplitudeuRubar; ParticleNeutralino5.Array_Decays[1][3] = 2; ParticleNeutralino5.Array_Comments[1] = "# ~chi_50 -> ub ~u_R";
     ParticleNeutralino5.Array_Decays[2][0] = PDGup; ParticleNeutralino5.Array_Decays[2][1] = -PDGsupL; ParticleNeutralino5.Array_Decays[2][2] = neutralino5amplitudeuLbaru; ParticleNeutralino5.Array_Decays[2][3] = 2; ParticleNeutralino5.Array_Comments[2] = "# ~chi_50 -> u ~u_L*";
     ParticleNeutralino5.Array_Decays[3][0] = PDGup; ParticleNeutralino5.Array_Decays[3][1] = -PDGsupR; ParticleNeutralino5.Array_Decays[3][2] = neutralino5amplitudeuRbaru; ParticleNeutralino5.Array_Decays[3][3] = 2; ParticleNeutralino5.Array_Comments[3] = "# ~chi_50 -> u ~u_R*";
     ParticleNeutralino5.Array_Decays[4][0] = -PDGdown; ParticleNeutralino5.Array_Decays[4][1] = PDGsdownL; ParticleNeutralino5.Array_Decays[4][2] = neutralino5amplitudedLdbar;ParticleNeutralino5.Array_Decays[4][3] = 2;  ParticleNeutralino5.Array_Comments[4] = "# ~chi_50 -> db ~d_L";
     ParticleNeutralino5.Array_Decays[5][0] = -PDGdown; ParticleNeutralino5.Array_Decays[5][1] = PDGsdownR; ParticleNeutralino5.Array_Decays[5][2] = neutralino5amplitudedRdbar; ParticleNeutralino5.Array_Decays[5][3] = 2; ParticleNeutralino5.Array_Comments[5] = "# ~chi_50 -> db ~d_R";
     ParticleNeutralino5.Array_Decays[6][0] = PDGdown; ParticleNeutralino5.Array_Decays[6][1] = -PDGsdownL; ParticleNeutralino5.Array_Decays[6][2] = neutralino5amplitudedLbard; ParticleNeutralino5.Array_Decays[6][3] = 2; ParticleNeutralino5.Array_Comments[6] = "# ~chi_50 -> d ~d_L*";
     ParticleNeutralino5.Array_Decays[7][0] = PDGdown; ParticleNeutralino5.Array_Decays[7][1] = -PDGsdownR; ParticleNeutralino5.Array_Decays[7][2] = neutralino5amplitudedRbard; ParticleNeutralino5.Array_Decays[7][3] = 2; ParticleNeutralino5.Array_Comments[7] = "# ~chi_50 -> d ~d_R*";
     ParticleNeutralino5.Array_Decays[8][0] = -PDGcharm; ParticleNeutralino5.Array_Decays[8][1] = PDGscharmL; ParticleNeutralino5.Array_Decays[8][2] = neutralino5amplitudecLcbar; ParticleNeutralino5.Array_Decays[8][3] = 2; ParticleNeutralino5.Array_Comments[8] = "# ~chi_50 -> cb ~c_L";
     ParticleNeutralino5.Array_Decays[9][0] = -PDGcharm; ParticleNeutralino5.Array_Decays[9][1] = PDGscharmR; ParticleNeutralino5.Array_Decays[9][2] = neutralino5amplitudecRcbar; ParticleNeutralino5.Array_Decays[9][3] = 2; ParticleNeutralino5.Array_Comments[9] = "# ~chi_50 -> cb ~c_R";
     ParticleNeutralino5.Array_Decays[10][0] = PDGcharm; ParticleNeutralino5.Array_Decays[10][1] = -PDGscharmL; ParticleNeutralino5.Array_Decays[10][2] = neutralino5amplitudecLbarc; ParticleNeutralino5.Array_Decays[10][3] = 2; ParticleNeutralino5.Array_Comments[10] = "# ~chi_50 -> c ~c_L*";
     ParticleNeutralino5.Array_Decays[11][0] = PDGcharm; ParticleNeutralino5.Array_Decays[11][1] = -PDGscharmR; ParticleNeutralino5.Array_Decays[11][2] = neutralino5amplitudecRbarc; ParticleNeutralino5.Array_Decays[11][3] = 2; ParticleNeutralino5.Array_Comments[11] = "# ~chi_50 -> c ~c_R*";
     ParticleNeutralino5.Array_Decays[12][0] = -PDGstrange; ParticleNeutralino5.Array_Decays[12][1] = PDGsstrangeL; ParticleNeutralino5.Array_Decays[12][2] = neutralino5amplitudesLsbar; ParticleNeutralino5.Array_Decays[12][3] = 2; ParticleNeutralino5.Array_Comments[12] = "# ~chi_50 -> sb ~s_L"; 
     ParticleNeutralino5.Array_Decays[13][0] = -PDGstrange; ParticleNeutralino5.Array_Decays[13][1] = PDGsstrangeR; ParticleNeutralino5.Array_Decays[13][2] = neutralino5amplitudesRsbar; ParticleNeutralino5.Array_Decays[13][3] = 2; ParticleNeutralino5.Array_Comments[13] = "# ~chi_50 -> sb ~s_R";
     ParticleNeutralino5.Array_Decays[14][0] = PDGstrange; ParticleNeutralino5.Array_Decays[14][1] = -PDGsstrangeL; ParticleNeutralino5.Array_Decays[14][2] = neutralino5amplitudesLbars; ParticleNeutralino5.Array_Decays[14][3] = 2; ParticleNeutralino5.Array_Comments[14] = "# ~chi_50 -> s ~s_L*";
     ParticleNeutralino5.Array_Decays[15][0] = PDGstrange; ParticleNeutralino5.Array_Decays[15][1] = -PDGsstrangeR; ParticleNeutralino5.Array_Decays[15][2] = neutralino5amplitudesRbars; ParticleNeutralino5.Array_Decays[15][3] = 2; ParticleNeutralino5.Array_Comments[15] = "# ~chi_50 -> s ~s_R*";
     ParticleNeutralino5.Array_Decays[16][0] = -PDGelectron; ParticleNeutralino5.Array_Decays[16][1] = PDGselectronL; ParticleNeutralino5.Array_Decays[16][2] = neutralino5amplitudeeLebar; ParticleNeutralino5.Array_Decays[16][3] = 2; ParticleNeutralino5.Array_Comments[16] = "# ~chi_50 -> e+ ~e_L-";
     ParticleNeutralino5.Array_Decays[17][0] = -PDGelectron; ParticleNeutralino5.Array_Decays[17][1] = PDGselectronR; ParticleNeutralino5.Array_Decays[17][2] = neutralino5amplitudeeRebar; ParticleNeutralino5.Array_Decays[17][3] = 2; ParticleNeutralino5.Array_Comments[17] = "# ~chi_50 -> e+ ~e_R-";
     ParticleNeutralino5.Array_Decays[18][0] = PDGelectron; ParticleNeutralino5.Array_Decays[18][1] = -PDGselectronL; ParticleNeutralino5.Array_Decays[18][2] = neutralino5amplitudeeLbare; ParticleNeutralino5.Array_Decays[18][3] = 2; ParticleNeutralino5.Array_Comments[18] = "# ~chi_50 -> e ~e_L+";
     ParticleNeutralino5.Array_Decays[19][0] = PDGelectron; ParticleNeutralino5.Array_Decays[19][1] = -PDGselectronR; ParticleNeutralino5.Array_Decays[19][2] = neutralino5amplitudeeRbare; ParticleNeutralino5.Array_Decays[19][3] = 2; ParticleNeutralino5.Array_Comments[19] = "# ~chi_50 -> e ~e_R+";   
     ParticleNeutralino5.Array_Decays[20][0] = -PDGmuon; ParticleNeutralino5.Array_Decays[20][1] = PDGsmuonL; ParticleNeutralino5.Array_Decays[20][2] = neutralino5amplitudemuLmubar; ParticleNeutralino5.Array_Decays[20][3] = 2; ParticleNeutralino5.Array_Comments[20] = "# ~chi_50 -> mu+ ~mu_L-";
     ParticleNeutralino5.Array_Decays[21][0] = -PDGmuon; ParticleNeutralino5.Array_Decays[21][1] = PDGsmuonR; ParticleNeutralino5.Array_Decays[21][2] = neutralino5amplitudemuRmubar; ParticleNeutralino5.Array_Decays[21][3] = 2; ParticleNeutralino5.Array_Comments[21] = "# ~chi_50 -> mu+ ~mu_R-";
     ParticleNeutralino5.Array_Decays[22][0] = PDGmuon; ParticleNeutralino5.Array_Decays[22][1] = -PDGsmuonL; ParticleNeutralino5.Array_Decays[22][2] = neutralino5amplitudemuLbarmu; ParticleNeutralino5.Array_Decays[22][3] = 2; ParticleNeutralino5.Array_Comments[22] = "# ~chi_50 -> mu- ~mu_L+";
     ParticleNeutralino5.Array_Decays[23][0] = PDGmuon; ParticleNeutralino5.Array_Decays[23][1] = -PDGsmuonR; ParticleNeutralino5.Array_Decays[23][2] = neutralino5amplitudemuRbarmu; ParticleNeutralino5.Array_Decays[23][3] = 2; ParticleNeutralino5.Array_Comments[23] = "# ~chi_50 -> mu- ~mu_R+";
     ParticleNeutralino5.Array_Decays[24][0] = PDGnuelectron; ParticleNeutralino5.Array_Decays[24][1] = -PDGnuselectronL; ParticleNeutralino5.Array_Decays[24][2] = neutralino5amplitudesnuebarnue; ParticleNeutralino5.Array_Decays[24][3] = 2; ParticleNeutralino5.Array_Comments[24] = "# ~chi_50 -> nu_e ~nu_eL*";
     ParticleNeutralino5.Array_Decays[25][0] = -PDGnuelectron; ParticleNeutralino5.Array_Decays[25][1] = PDGnuselectronL; ParticleNeutralino5.Array_Decays[25][2] = neutralino5amplitudesnuenuebar; ParticleNeutralino5.Array_Decays[25][3] = 2; ParticleNeutralino5.Array_Comments[25] = "# ~chi_50 -> nu_eb ~nu_eL";
     ParticleNeutralino5.Array_Decays[26][0] = PDGnumuon; ParticleNeutralino5.Array_Decays[26][1] = -PDGnusmuonL; ParticleNeutralino5.Array_Decays[26][2] = neutralino5amplitudesnumubarnumu; ParticleNeutralino5.Array_Decays[26][3] = 2; ParticleNeutralino5.Array_Comments[26] = "# ~chi_50 -> numu ~nu_muL*";
     ParticleNeutralino5.Array_Decays[27][0] = -PDGnumuon; ParticleNeutralino5.Array_Decays[27][1] = PDGnusmuonL; ParticleNeutralino5.Array_Decays[27][2] = neutralino5amplitudesnumunumubar; ParticleNeutralino5.Array_Decays[27][3] = 2; ParticleNeutralino5.Array_Comments[27] = "# ~chi_50 -> nu_mub ~nu_muL";
     ParticleNeutralino5.Array_Decays[28][0] = PDGtop; ParticleNeutralino5.Array_Decays[28][1] = -PDGstop1; ParticleNeutralino5.Array_Decays[28][2] = neutralino5amplitudetopstop1bar; ParticleNeutralino5.Array_Decays[28][3] = 2; ParticleNeutralino5.Array_Comments[28] = "# ~chi_50 -> t ~t_1*";
     ParticleNeutralino5.Array_Decays[29][0] = PDGtop; ParticleNeutralino5.Array_Decays[29][1] = -PDGstop2; ParticleNeutralino5.Array_Decays[29][2] = neutralino5amplitudetopstop2bar; ParticleNeutralino5.Array_Decays[29][3] = 2; ParticleNeutralino5.Array_Comments[29] = "# ~chi_50 -> t ~t_2*";
     ParticleNeutralino5.Array_Decays[30][0] = -PDGtop; ParticleNeutralino5.Array_Decays[30][1] = PDGstop1; ParticleNeutralino5.Array_Decays[30][2] = neutralino5amplitudetopbarstop1; ParticleNeutralino5.Array_Decays[30][3] = 2; ParticleNeutralino5.Array_Comments[30] = "# ~chi_50 -> tb ~t_1";
     ParticleNeutralino5.Array_Decays[31][0] = -PDGtop; ParticleNeutralino5.Array_Decays[31][1] = PDGstop2; ParticleNeutralino5.Array_Decays[31][2] = neutralino5amplitudetopbarstop2; ParticleNeutralino5.Array_Decays[31][3] = 2; ParticleNeutralino5.Array_Comments[31] = "# ~chi_50 -> tb ~t_2";
     ParticleNeutralino5.Array_Decays[32][0] = PDGbottom; ParticleNeutralino5.Array_Decays[32][1] = -PDGsbottom1; ParticleNeutralino5.Array_Decays[32][2] = neutralino5amplitudebottomsbottom1bar; ParticleNeutralino5.Array_Decays[32][3] = 2; ParticleNeutralino5.Array_Comments[32] = "# ~chi_50 -> b ~b_1*";
     ParticleNeutralino5.Array_Decays[33][0] = PDGbottom; ParticleNeutralino5.Array_Decays[33][1] = -PDGsbottom2; ParticleNeutralino5.Array_Decays[33][2] = neutralino5amplitudebottomsbottom2bar; ParticleNeutralino5.Array_Decays[33][3] = 2; ParticleNeutralino5.Array_Comments[33] = "# ~chi_50 -> b ~b_2*";
     ParticleNeutralino5.Array_Decays[34][0] = -PDGbottom; ParticleNeutralino5.Array_Decays[34][1] = PDGsbottom1; ParticleNeutralino5.Array_Decays[34][2] = neutralino5amplitudebottombarsbottom1; ParticleNeutralino5.Array_Decays[34][3] = 2; ParticleNeutralino5.Array_Comments[34] = "# ~chi_50 -> bb ~b_1";
     ParticleNeutralino5.Array_Decays[35][0] = -PDGbottom; ParticleNeutralino5.Array_Decays[35][1] = PDGsbottom2; ParticleNeutralino5.Array_Decays[35][2] = neutralino5amplitudebottombarsbottom2; ParticleNeutralino5.Array_Decays[35][3] = 2; ParticleNeutralino5.Array_Comments[35] = "# ~chi_50 -> bb ~b_2";
     ParticleNeutralino5.Array_Decays[36][0] = -PDGstau1; ParticleNeutralino5.Array_Decays[36][1] = PDGtau; ParticleNeutralino5.Array_Decays[36][2] = neutralino5amplitudetaustau1bar; ParticleNeutralino5.Array_Decays[36][3] = 2; ParticleNeutralino5.Array_Comments[36] = "# ~chi_50 -> tau- ~tau_1+";
     ParticleNeutralino5.Array_Decays[37][0] = -PDGstau2; ParticleNeutralino5.Array_Decays[37][1] = PDGtau; ParticleNeutralino5.Array_Decays[37][2] = neutralino5amplitudetaustau2bar; ParticleNeutralino5.Array_Decays[37][3] = 2; ParticleNeutralino5.Array_Comments[37] = "# ~chi_50 -> tau- ~tau_2+";
     ParticleNeutralino5.Array_Decays[38][0] = PDGstau1; ParticleNeutralino5.Array_Decays[38][1] = -PDGtau; ParticleNeutralino5.Array_Decays[38][2] = neutralino5amplitudetaubarstau1; ParticleNeutralino5.Array_Decays[38][3] = 2; ParticleNeutralino5.Array_Comments[38] = "# ~chi_50 -> tau+ ~tau_1-";
     ParticleNeutralino5.Array_Decays[39][0] = PDGstau2; ParticleNeutralino5.Array_Decays[39][1] = -PDGtau; ParticleNeutralino5.Array_Decays[39][2] = neutralino5amplitudetaubarstau2; ParticleNeutralino5.Array_Decays[39][3] = 2; ParticleNeutralino5.Array_Comments[39] = "# ~chi_50 -> tau+ ~tau_2-";
     ParticleNeutralino5.Array_Decays[40][0] = PDGnutau; ParticleNeutralino5.Array_Decays[40][1] = -PDGnustauL; ParticleNeutralino5.Array_Decays[40][2] = neutralino5amplitudenutausnutaubar; ParticleNeutralino5.Array_Decays[40][3] = 2; ParticleNeutralino5.Array_Comments[40] = "# ~chi_50 -> nu_tau ~nu_tauL*";
     ParticleNeutralino5.Array_Decays[41][0] = -PDGnutau; ParticleNeutralino5.Array_Decays[41][1] = PDGnustauL; ParticleNeutralino5.Array_Decays[41][2] = neutralino5amplitudenutaubarsnutau; ParticleNeutralino5.Array_Decays[41][3] = 2; ParticleNeutralino5.Array_Comments[41] = "# ~chi_50 -> nu_taub ~nu_tauL";
     ParticleNeutralino5.Array_Decays[42][0] = PDGWplus; ParticleNeutralino5.Array_Decays[42][1] = -PDGchargino1; ParticleNeutralino5.Array_Decays[42][2] = neutralino5amplitudeWbosonpluscharginoW1; ParticleNeutralino5.Array_Decays[42][3] = 2; ParticleNeutralino5.Array_Comments[42] = "# ~chi_50 -> W+ ~chi_1-";
     ParticleNeutralino5.Array_Decays[43][0] = PDGWplus; ParticleNeutralino5.Array_Decays[43][1] = -PDGchargino2; ParticleNeutralino5.Array_Decays[43][2] = neutralino5amplitudeWbosonpluscharginoW2; ParticleNeutralino5.Array_Decays[43][3] = 2; ParticleNeutralino5.Array_Comments[43] = "# ~chi_50 -> W+ ~chi_2-";
     ParticleNeutralino5.Array_Decays[44][0] = -PDGWplus; ParticleNeutralino5.Array_Decays[44][1] = PDGchargino1; ParticleNeutralino5.Array_Decays[44][2] = neutralino5amplitudeWbosonminuscharginoW1; ParticleNeutralino5.Array_Decays[44][3] = 2; ParticleNeutralino5.Array_Comments[44] = "# ~chi_50 -> W- ~chi_1+";
     ParticleNeutralino5.Array_Decays[45][0] = -PDGWplus; ParticleNeutralino5.Array_Decays[45][1] = PDGchargino2; ParticleNeutralino5.Array_Decays[45][2] = neutralino5amplitudeWbosonminuscharginoW2; ParticleNeutralino5.Array_Decays[45][3] = 2; ParticleNeutralino5.Array_Comments[45] = "# ~chi_50 -> W- ~chi_2+";
     ParticleNeutralino5.Array_Decays[46][0] = PDGHplus; ParticleNeutralino5.Array_Decays[46][1] = -PDGchargino1; ParticleNeutralino5.Array_Decays[46][2] = neutralino5amplitudeHpluscharginoW1; ParticleNeutralino5.Array_Decays[46][3] = 2; ParticleNeutralino5.Array_Comments[46] = "# ~chi_50 -> H+ ~chi_1-";
     ParticleNeutralino5.Array_Decays[47][0] = PDGHplus; ParticleNeutralino5.Array_Decays[47][1] = -PDGchargino2; ParticleNeutralino5.Array_Decays[47][2] = neutralino5amplitudeHpluscharginoW2; ParticleNeutralino5.Array_Decays[47][3] = 2; ParticleNeutralino5.Array_Comments[47] = "# ~chi_50 -> H+ ~chi_2-";
     ParticleNeutralino5.Array_Decays[48][0] = -PDGHplus; ParticleNeutralino5.Array_Decays[48][1] = PDGchargino1; ParticleNeutralino5.Array_Decays[48][2] = neutralino5amplitudeHminuscharginoW1; ParticleNeutralino5.Array_Decays[48][3] = 2; ParticleNeutralino5.Array_Comments[48] = "# ~chi_50 -> H- ~chi_1+";
     ParticleNeutralino5.Array_Decays[49][0] = -PDGHplus; ParticleNeutralino5.Array_Decays[49][1] = PDGchargino2; ParticleNeutralino5.Array_Decays[49][2] = neutralino5amplitudeHminuscharginoW2; ParticleNeutralino5.Array_Decays[49][3] = 2; ParticleNeutralino5.Array_Comments[49] = "# ~chi_50 -> H- ~chi_2+";
     ParticleNeutralino5.Array_Decays[50][0] = PDGZboson; ParticleNeutralino5.Array_Decays[50][1] = PDGneutralino1; ParticleNeutralino5.Array_Decays[50][2] = neutralino5amplitudeZbosonneutralino1; ParticleNeutralino5.Array_Decays[50][3] = 2; ParticleNeutralino5.Array_Comments[50] = "# ~chi_50 -> Z ~chi_10";
     ParticleNeutralino5.Array_Decays[51][0] = PDGZboson; ParticleNeutralino5.Array_Decays[51][1] = PDGneutralino2; ParticleNeutralino5.Array_Decays[51][2] = neutralino5amplitudeZbosonneutralino2; ParticleNeutralino5.Array_Decays[51][3] = 2; ParticleNeutralino5.Array_Comments[51] = "# ~chi_50 -> Z ~chi_20";
     ParticleNeutralino5.Array_Decays[52][0] = PDGZboson; ParticleNeutralino5.Array_Decays[52][1] = PDGneutralino3; ParticleNeutralino5.Array_Decays[52][2] = neutralino5amplitudeZbosonneutralino3; ParticleNeutralino5.Array_Decays[52][3] = 2; ParticleNeutralino5.Array_Comments[52] = "# ~chi_50 -> Z ~chi_30";
     ParticleNeutralino5.Array_Decays[53][0] = PDGZboson; ParticleNeutralino5.Array_Decays[53][1] = PDGneutralino4; ParticleNeutralino5.Array_Decays[53][2] = neutralino5amplitudeZbosonneutralino4; ParticleNeutralino5.Array_Decays[53][3] = 2; ParticleNeutralino5.Array_Comments[53] = "# ~chi_50 -> Z ~chi_40";
     ParticleNeutralino5.Array_Decays[54][0] = PDGh0; ParticleNeutralino5.Array_Decays[54][1] = PDGneutralino1; ParticleNeutralino5.Array_Decays[54][2] = neutralino5amplitudehneutralino1; ParticleNeutralino5.Array_Decays[54][3] = 2; ParticleNeutralino5.Array_Comments[54] = "# ~chi_50 -> h ~chi_10";
     ParticleNeutralino5.Array_Decays[55][0] = PDGh0; ParticleNeutralino5.Array_Decays[55][1] = PDGneutralino2; ParticleNeutralino5.Array_Decays[55][2] = neutralino5amplitudehneutralino2; ParticleNeutralino5.Array_Decays[55][3] = 2; ParticleNeutralino5.Array_Comments[55] = "# ~chi_50 -> h ~chi_20"; 
     ParticleNeutralino5.Array_Decays[56][0] = PDGh0; ParticleNeutralino5.Array_Decays[56][1] = PDGneutralino3; ParticleNeutralino5.Array_Decays[56][2] = neutralino5amplitudehneutralino3; ParticleNeutralino5.Array_Decays[56][3] = 2; ParticleNeutralino5.Array_Comments[56] = "# ~chi_50 -> h ~chi_30";
     ParticleNeutralino5.Array_Decays[57][0] = PDGh0; ParticleNeutralino5.Array_Decays[57][1] = PDGneutralino4; ParticleNeutralino5.Array_Decays[57][2] = neutralino5amplitudehneutralino4; ParticleNeutralino5.Array_Decays[57][3] = 2; ParticleNeutralino5.Array_Comments[57] = "# ~chi_50 -> h ~chi_40";
     ParticleNeutralino5.Array_Decays[58][0] = PDGH0; ParticleNeutralino5.Array_Decays[58][1] = PDGneutralino1; ParticleNeutralino5.Array_Decays[58][2] = neutralino5amplitudeHneutralino1; ParticleNeutralino5.Array_Decays[58][3] = 2; ParticleNeutralino5.Array_Comments[58] = "# ~chi_50 -> H ~chi_10";
     ParticleNeutralino5.Array_Decays[59][0] = PDGH0; ParticleNeutralino5.Array_Decays[59][1] = PDGneutralino2; ParticleNeutralino5.Array_Decays[59][2] = neutralino5amplitudeHneutralino2; ParticleNeutralino5.Array_Decays[59][3] = 2; ParticleNeutralino5.Array_Comments[59] = "# ~chi_50 -> H ~chi_20";
     ParticleNeutralino5.Array_Decays[60][0] = PDGH0; ParticleNeutralino5.Array_Decays[60][1] = PDGneutralino3; ParticleNeutralino5.Array_Decays[60][2] = neutralino5amplitudeHneutralino3; ParticleNeutralino5.Array_Decays[60][3] = 2; ParticleNeutralino5.Array_Comments[60] = "# ~chi_50 -> H ~chi_30";
     ParticleNeutralino5.Array_Decays[61][0] = PDGH0; ParticleNeutralino5.Array_Decays[61][1] = PDGneutralino4; ParticleNeutralino5.Array_Decays[61][2] = neutralino5amplitudeHneutralino4; ParticleNeutralino5.Array_Decays[61][3] = 2; ParticleNeutralino5.Array_Comments[61] = "# ~chi_50 -> H ~chi_40";
     ParticleNeutralino5.Array_Decays[62][0] = PDGA0; ParticleNeutralino5.Array_Decays[62][1] = PDGneutralino1; ParticleNeutralino5.Array_Decays[62][2] = neutralino5amplitudeAneutralino1; ParticleNeutralino5.Array_Decays[62][3] = 2; ParticleNeutralino5.Array_Comments[62] = "# ~chi_50 -> A ~chi_10";
     ParticleNeutralino5.Array_Decays[63][0] = PDGA0; ParticleNeutralino5.Array_Decays[63][1] = PDGneutralino2; ParticleNeutralino5.Array_Decays[63][2] = neutralino5amplitudeAneutralino2; ParticleNeutralino5.Array_Decays[63][3] = 2; ParticleNeutralino5.Array_Comments[63] = "# ~chi_50 -> A ~chi_20";
     ParticleNeutralino5.Array_Decays[64][0] = PDGA0; ParticleNeutralino5.Array_Decays[64][1] = PDGneutralino3; ParticleNeutralino5.Array_Decays[64][2] = neutralino5amplitudeAneutralino3; ParticleNeutralino5.Array_Decays[64][3] = 2; ParticleNeutralino5.Array_Comments[64] = "# ~chi_50 -> A ~chi_30";
     ParticleNeutralino5.Array_Decays[65][0] = PDGA0; ParticleNeutralino5.Array_Decays[65][1] = PDGneutralino4; ParticleNeutralino5.Array_Decays[65][2] = neutralino5amplitudeAneutralino4; ParticleNeutralino5.Array_Decays[65][3] = 2; ParticleNeutralino5.Array_Comments[65] = "# ~chi_50 -> A ~chi_40";
     ParticleNeutralino5.Array_Decays[66][0] = PDGH3; ParticleNeutralino5.Array_Decays[66][1] = PDGneutralino1; ParticleNeutralino5.Array_Decays[66][2] = neutralino5amplitudeH3neutralino1; ParticleNeutralino5.Array_Decays[66][3] = 2; ParticleNeutralino5.Array_Comments[66] = "# ~chi_50 -> H3 ~chi_10";
     ParticleNeutralino5.Array_Decays[67][0] = PDGH3; ParticleNeutralino5.Array_Decays[67][1] = PDGneutralino2; ParticleNeutralino5.Array_Decays[67][2] = neutralino5amplitudeH3neutralino2; ParticleNeutralino5.Array_Decays[67][3] = 2; ParticleNeutralino5.Array_Comments[67] = "# ~chi_50 -> H3 ~chi_20";
     ParticleNeutralino5.Array_Decays[68][0] = PDGH3; ParticleNeutralino5.Array_Decays[68][1] = PDGneutralino3; ParticleNeutralino5.Array_Decays[68][2] = neutralino5amplitudeH3neutralino3; ParticleNeutralino5.Array_Decays[68][3] = 2; ParticleNeutralino5.Array_Comments[68] = "# ~chi_50 -> H3 ~chi_30";
     ParticleNeutralino5.Array_Decays[69][0] = PDGH3; ParticleNeutralino5.Array_Decays[69][1] = PDGneutralino4; ParticleNeutralino5.Array_Decays[69][2] = neutralino5amplitudeH3neutralino4; ParticleNeutralino5.Array_Decays[69][3] = 2; ParticleNeutralino5.Array_Comments[69] = "# ~chi_50 -> H3 ~chi_40";
     ParticleNeutralino5.Array_Decays[70][0] = PDGA2; ParticleNeutralino5.Array_Decays[70][1] = PDGneutralino1; ParticleNeutralino5.Array_Decays[70][2] = neutralino5amplitudeA2neutralino1; ParticleNeutralino5.Array_Decays[70][3] = 2; ParticleNeutralino5.Array_Comments[70] = "# ~chi_50 -> A2 ~chi_10";
     ParticleNeutralino5.Array_Decays[71][0] = PDGA2; ParticleNeutralino5.Array_Decays[71][1] = PDGneutralino2; ParticleNeutralino5.Array_Decays[71][2] = neutralino5amplitudeA2neutralino2; ParticleNeutralino5.Array_Decays[71][3] = 2; ParticleNeutralino5.Array_Comments[71] = "# ~chi_50 -> A2 ~chi_20";
     ParticleNeutralino5.Array_Decays[72][0] = PDGA2; ParticleNeutralino5.Array_Decays[72][1] = PDGneutralino3; ParticleNeutralino5.Array_Decays[72][2] = neutralino5amplitudeA2neutralino3; ParticleNeutralino5.Array_Decays[72][3] = 2; ParticleNeutralino5.Array_Comments[72] = "# ~chi_50 -> A2 ~chi_30";
     ParticleNeutralino5.Array_Decays[73][0] = PDGA2; ParticleNeutralino5.Array_Decays[73][1] = PDGneutralino4; ParticleNeutralino5.Array_Decays[73][2] = neutralino5amplitudeA2neutralino4; ParticleNeutralino5.Array_Decays[73][3] = 2; ParticleNeutralino5.Array_Comments[73] = "# ~chi_50 -> A2 ~chi_40";

     for(int i = 0; i<ParticleNeutralino5.No_of_Decays; i++) {
       if (ParticleNeutralino5.Array_Decays[i][2] < 0) {
	 fout << "#warning! Partial Width for " << ParticleNeutralino5.Array_Comments[i] << " is negative = " << ParticleNeutralino5.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
	 ParticleNeutralino5.Array_Decays[i][2] = 0;
	 errorflag = -1;
       }
     }     
     
     double Neut5_No_1to2_Decays = 0;
     
     Neut5_No_1to2_Decays = ParticleNeutralino5.No_1to2_Decays + ParticleNeutralino5.No_grav_Decays + ParticleNeutralino5.No_NMSSM_Decays;     
     for (int j = 0; j<Neut5_No_1to2_Decays; j++) {
       ParticleNeutralino5.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
     }
     
     for (int j=0; j<Neut5_No_1to2_Decays; j++) {
       ParticleNeutralino5.two_width = ParticleNeutralino5.two_width + ParticleNeutralino5.Array_Decays[j][2];
     }
     for (int j=Neut5_No_1to2_Decays; j<ParticleNeutralino5.No_of_Decays; j++) {
       ParticleNeutralino5.three_width = ParticleNeutralino5.three_width + ParticleNeutralino5.Array_Decays[j][2]; ///0 as no 1->3 decays included yet in NMSSM
     }
     
     if ( ParticleNeutralino5.three_width != ParticleNeutralino5.three_width) /// Tests for a nan as only nans aren't equal to themselves
       {
	 fout << "# Three body decays give nan for neutralino 5 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
	 errorflag = -1;
	 ParticleNeutralino5.No_of_Decays = Neut5_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
	 ParticleNeutralino5.total_width = ParticleNeutralino5.two_width;
       }
     else {
       ParticleNeutralino5.total_width = ParticleNeutralino5.two_width + ParticleNeutralino5.three_width;
     }
     
     if ( ParticleNeutralino5.total_width != ParticleNeutralino5.total_width) /// Tests for a nan as only nans aren't equal to themselves
       {
	 errorflag = -1;
	 // for (int i = 0; i<ParticleNeutralino5.No_of_Decays; i++) {
	 //   fout << i << " " << ParticleNeutralino5.Array_Decays[i][2] << endl;
	 // }	  
	 throw( "nan in Neutralino5 total width \n");
       }
     
   }

 }

///Higgs decays
///higgsl decays

 double h0amplitudeuantiu=0, h0amplitudedantid=0, h0amplitudesantis=0, h0amplitudecantic=0, h0amplitudebantib=0, h0amplitudeeantie=0, h0amplitudemuantimu=0, h0amplitudetauantitau=0, h0amplitudeneutZ1neutZ1=0, h0amplitudeneutZ1neutZ2=0, h0amplitudeneutZ1neutZ3=0, h0amplitudeneutZ1neutZ4=0, h0amplitudeneutZ2neutZ2=0, h0amplitudeneutZ2neutZ3=0, h0amplitudeneutZ2neutZ4=0, h0amplitudeneutZ3neutZ3=0, h0amplitudeneutZ3neutZ4=0, h0amplitudeneutZ4neutZ4=0, h0amplitudecharW1charW1=0, h0amplitudecharW1charW2=0, h0amplitudecharW2charW2=0, h0amplitudehiggsAhiggsA=0, h0amplitudehiggsAZboson=0, h0amplitudesupLantisupL=0, h0amplitudesupLantisupR=0, h0amplitudesupRantisupL=0, h0amplitudesupRantisupR=0, h0amplitudesdownLantisdownL=0, h0amplitudesdownLantisdownR=0, h0amplitudesdownRantisdownL=0, h0amplitudesdownRantisdownR=0, h0amplitudescharmLantischarmL=0, h0amplitudescharmLantischarmR=0, h0amplitudescharmRantischarmL=0, h0amplitudescharmRantischarmR=0, h0amplitudesstrangeLantisstrangeL=0, h0amplitudesstrangeLantisstrangeR=0, h0amplitudesstrangeRantisstrangeL=0, h0amplitudesstrangeRantisstrangeR=0, h0amplitudesnueLantisnueL=0, h0amplitudeselectronLantiselectronL=0, h0amplitudeselectronRantiselectronR=0, h0amplitudeselectronLantiselectronR=0, h0amplitudeselectronRantiselectronL=0, h0amplitudesnumuLantisnumuL=0, h0amplitudesnutauLantisnutauL=0, h0amplitudesmuonLantismuonL=0, h0amplitudesmuonRantismuonR=0, h0amplitudesmuonLantismuonR=0, h0amplitudesmuonRantismuonL=0, h0amplitudestau1antistau1=0, h0amplitudestau2antistau2=0, h0amplitudestau1antistau2=0, h0amplitudestau2antistau1=0, h0amplitudestop1antistop1=0, h0amplitudestop1antistop2=0, h0amplitudestop2antistop1=0, h0amplitudestop2antistop2=0, h0amplitudesbottom1antisbottom1=0, h0amplitudesbottom1antisbottom2=0, h0amplitudesbottom2antisbottom1=0, h0amplitudesbottom2antisbottom2=0, h0amplitudegluongluon=0, h0amplitudegammagamma=0, h0amplitudeWW=0, h0amplitudeZZ=0, h0amplitudeZgamma=0,h0amplitudehiggsAhiggsA2=0, h0amplitudehiggsA2higgsA2=0, h0amplitudehiggsA2Zboson=0, h0amplitudeHpmHpm = 0, h0amplitudeWHpm = 0;
 
 double h0amplitudeneutZ1neutZ5 = 0, h0amplitudeneutZ2neutZ5 = 0, h0amplitudeneutZ3neutZ5 = 0, h0amplitudeneutZ4neutZ5 = 0, h0amplitudeneutZ5neutZ5 = 0;

 if (flagh1 == 1) {
   ///If QCDcorr off use these below -> running masses used to approximate some of QCD corrections at tree-level
   if (QCDcorr == false) {
     ///No decays to u or d as PWs to u and d are tiny as proportional to yukawas squared
     ///Use running masses here to try to approximate some of the correction (which aren't included)
     h0amplitudecantic = higgslorHamplitudedecayquarkantiquark (mh0(1), runmc, g, alpha, beta, runmw, 1, 'l', CPEMix, nmssmIsIt, QCDcorr, alphasAtMh);
     h0amplitudesantis = higgslorHamplitudedecayquarkantiquark (mh0(1), runms, g, alpha, beta, runmw, 0, 'l', CPEMix, nmssmIsIt, QCDcorr, alphasAtMh);
     h0amplitudebantib = higgslorHamplitudedecayquarkantiquark (mh0(1), runmb, g, alpha, beta, runmw, 0, 'l', CPEMix, nmssmIsIt, QCDcorr, alphasAtMh); 
     ///No decay to two tops as kinematically forbidden
   }
   ///With QCDcorr on then use actual pole masses for quarks as corrections being accounted for
   else {
     h0amplitudecantic = higgslorHamplitudedecayquarkantiquark (mh0(1), mcpole, g, alpha, beta, runmw, 1, 'l', CPEMix, nmssmIsIt, QCDcorr, alphasAtMh);
     h0amplitudesantis = higgslorHamplitudedecayquarkantiquark (mh0(1), mspole, g, alpha, beta, runmw, 0, 'l', CPEMix, nmssmIsIt, QCDcorr, alphasAtMh);
     ///mcpole and mspole set in decays.h, this values used are those appropriate for the scheme used for the h -> qq QCD corrections, as in hdecay
     h0amplitudebantib = higgslorHamplitudedecayquarkantiquark (mh0(1), mbPole, g, alpha, beta, runmw, 0, 'l', CPEMix, nmssmIsIt, QCDcorr, alphasAtMh); 
     ///No decay to two tops as kinematically forbidden
   }

   h0amplitudeeantie = higgslorHamplitudedecayquarkantiquark (mh0(1), runmel, g, alpha, beta, runmw, 0, 'l', CPEMix, nmssmIsIt, false, alphasAtMh)/3; ///0 as leptons are like down-type quarks, divided by three as Nc is three for quarks but 1 for leptons
   h0amplitudemuantimu = higgslorHamplitudedecayquarkantiquark (mh0(1), runmmu, g, alpha, beta, runmw, 0, 'l', CPEMix, nmssmIsIt, false, alphasAtMh)/3;
   h0amplitudetauantitau = higgslorHamplitudedecayquarkantiquark (mh0(1), runmtau, g, alpha, beta, runmw, 0, 'l', CPEMix, nmssmIsIt, false, alphasAtMh)/3;
   
   if (nmssmIsIt == false) {
     h0amplitudesupLantisupL = higgshamplitudedecay2squarksamehand (mh0(1), mu(1,1), mu(1,1), g, gp, alpha, beta, runmw, runmu, runmd, 1);
     h0amplitudesupRantisupR = higgshamplitudedecay2squarksamehand (mh0(1), mu(2,1), mu(2,1), g, gp, alpha, beta, runmw, runmu, runmd, 3);
     h0amplitudesdownLantisdownL = higgshamplitudedecay2squarksamehand (mh0(1), md(1,1), md(1,1), g, gp, alpha, beta, runmw, runmu, runmd, 2);
     h0amplitudesdownRantisdownR = higgshamplitudedecay2squarksamehand (mh0(1), md(2,1), mu(2,1), g, gp, alpha, beta, runmw, runmu, runmd, 4);
     h0amplitudescharmLantischarmL = higgshamplitudedecay2squarksamehand (mh0(1), mu(1,2), mu(1,2), g, gp, alpha, beta, runmw, runmc, runms, 1);
     h0amplitudescharmRantischarmR = higgshamplitudedecay2squarksamehand (mh0(1), mu(2,2), mu(2,2), g, gp, alpha, beta, runmw, runmc, runms, 3);
     h0amplitudesstrangeLantisstrangeL = higgshamplitudedecay2squarksamehand (mh0(1), md(2,1), md(2,1), g, gp, alpha, beta, runmw, runmc, runms, 2);
     h0amplitudesstrangeRantisstrangeR = higgshamplitudedecay2squarksamehand (mh0(1), md(2,2), md(2,2), g, gp, alpha, beta, runmw, runmc, runms, 4);
     h0amplitudesnueLantisnueL = higgshamplitudedecay2sleptonsamehand (mh0(1), msnu(1), msnu(1), g, gp, alpha, beta, runmw, runmel, 1);
     h0amplitudeselectronLantiselectronL = higgshamplitudedecay2sleptonsamehand (mh0(1), me(1,1), me(1,1), g, gp, alpha, beta, runmw, runmel, 2);
     h0amplitudeselectronRantiselectronR = higgshamplitudedecay2sleptonsamehand (mh0(1), me(2,1), me(2,1), g, gp, alpha, beta, runmw, runmel, 3);
     h0amplitudesnumuLantisnumuL = higgshamplitudedecay2sleptonsamehand (mh0(1), msnu(2), msnu(2), g, gp, alpha, beta, runmw, runmmu, 1);
     h0amplitudesmuonLantismuonL = higgshamplitudedecay2sleptonsamehand (mh0(1), me(1,2), me(1,2), g, gp, alpha, beta, runmw, runmmu, 2);
     h0amplitudesmuonRantismuonR = higgshamplitudedecay2sleptonsamehand (mh0(1), me(2,2), me(2,2), g, gp, alpha, beta, runmw, runmmu, 3);
     h0amplitudesnutauLantisnutauL = higgshamplitudedecay2sleptonsamehand (mh0(1), msnu(3), msnu(3), g, gp, alpha, beta, runmw, runmtau, 1);
     h0amplitudesupLantisupR = higgshamplitudedecay2squarkdiffhand (mh0(1), mu(1,1), mu(2,1), g, alpha, beta, runmw, runmu, runmd, greekmu, Au, Ad, 1);
     h0amplitudesupRantisupL = higgshamplitudedecay2squarkdiffhand (mh0(1), mu(1,1), mu(2,1), g, alpha, beta, runmw, runmu, runmd, greekmu, Au, Ad, 1);
     h0amplitudesdownLantisdownR = higgshamplitudedecay2squarkdiffhand (mh0(1), md(1,1), md(2,1), g, alpha, beta, runmw, runmu, runmd, greekmu, Au, Ad, 2);
     h0amplitudesdownRantisdownL = higgshamplitudedecay2squarkdiffhand (mh0(1), md(1,1), md(2,1), g, alpha, beta, runmw, runmu, runmd, greekmu, Au, Ad, 2);
     h0amplitudescharmLantischarmR = higgshamplitudedecay2squarkdiffhand (mh0(1), mu(1,2), mu(2,2), g, alpha, beta, runmw, runmc, runms, greekmu, Ac, As, 1);
     h0amplitudescharmRantischarmL = higgshamplitudedecay2squarkdiffhand (mh0(1), mu(1,2), mu(2,2), g, alpha, beta, runmw, runmc, runms, greekmu, Ac, As, 1);
     h0amplitudesstrangeLantisstrangeR = higgshamplitudedecay2squarkdiffhand (mh0(1), md(2,1), md(2,2), g, alpha, beta, runmw, runmc, runms, greekmu, Ac, As, 2);
     h0amplitudesstrangeRantisstrangeL = higgshamplitudedecay2squarkdiffhand (mh0(1), md(2,1), md(2,2), g, alpha, beta, runmw, runmc, runms, greekmu, Ac, As, 2);
     h0amplitudeselectronLantiselectronR = higgshamplitudedecay2sleptondiffhand (mh0(1), me(1,1), me(2,1), g, alpha, beta, runmw, runmel, greekmu, Ae, 1);
     h0amplitudeselectronRantiselectronL = higgshamplitudedecay2sleptondiffhand (mh0(1), me(2,1), me(1,1), g, alpha, beta, runmw, runmel, greekmu, Ae, 1);
     h0amplitudesmuonLantismuonR = higgshamplitudedecay2sleptondiffhand (mh0(1), me(1,2), me(2,2), g, alpha, beta, runmw, runmmu, greekmu, Amu, 1);
     h0amplitudesmuonRantismuonL = higgshamplitudedecay2sleptondiffhand (mh0(1), me(2,2), me(2,1), g, alpha, beta, runmw, runmmu, greekmu, Amu, 1);
     h0amplitudestop1antistop1 = higgshamplitudedecaystop1stop1 (mh0(1), mu(1,3), mu(1,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetat); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     h0amplitudestop2antistop2 = higgshamplitudedecaystop2stop2 (mh0(1), mu(2,3), mu(2,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetat); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     h0amplitudestop1antistop2 = higgshamplitudedecaystop1stop2 (mh0(1), mu(1,3), mu(2,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetat); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     h0amplitudestop2antistop1 = higgshamplitudedecaystop1stop2 (mh0(1), mu(2,3), mu(1,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetat); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     h0amplitudesbottom1antisbottom1 = higgshamplitudedecaysbottom1sbottom1(mh0(1), md(1,3), md(1,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetab); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     h0amplitudesbottom2antisbottom2 = higgshamplitudedecaysbottom2sbottom2(mh0(1), md(2,3), md(2,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetab); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     h0amplitudesbottom1antisbottom2 = higgshamplitudedecaysbottom1sbottom2 (mh0(1), md(1,3), md(2,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetab); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     h0amplitudesbottom2antisbottom1 = higgshamplitudedecaysbottom1sbottom2 (mh0(1), md(2,3), md(1,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetab); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     h0amplitudestau1antistau1 = higgshamplitudedecaystau1stau1 (mh0(1), me(1,3), me(1,3), g, gp, alpha, beta, runmw, runmtau, greekmu, Atau, thetatau);
     h0amplitudestau2antistau2 = higgshamplitudedecaystau2stau2 (mh0(1), me(2,3), me(2,3), g, gp, alpha, beta, runmw, runmtau, greekmu, Atau, thetatau);
     h0amplitudestau1antistau2 = higgshamplitudedecaystau1stau2 (mh0(1), me(1,3), me(2,3), g, gp, alpha, beta, runmw, runmtau, greekmu, Atau, thetatau);
     h0amplitudestau2antistau1 = higgshamplitudedecaystau1stau2 (mh0(1), me(2,3), me(1,3), g, gp, alpha, beta, runmw, runmtau, greekmu, Atau, thetatau);
     
     h0amplitudecharW1charW1 = higgsphiamplitudedecaysamechargino (mh0(1), mch(1), g, thetaL2, thetaR2, alpha, beta, 1, 'h');
     h0amplitudecharW2charW2 = higgsphiamplitudedecaysamechargino (mh0(1), mch(2), g, thetaL2, thetaR2, alpha, beta, 2, 'h');
     h0amplitudecharW1charW2 = higgsphiamplitudedecaydifchargino (mh0(1), mch(1), mch(2), g, thetaL2, thetaR2, alpha, beta, 'h');

     h0amplitudegammagamma = higgsesamplitudedecaygammagammatotal(mh0(1), g, gp, alphaAtMh, runmw, polemw, alpha, beta, mtAtMh, mbAtMh, mcAtMh, runmtau, mHpm, mu(1,3), mu(2,3), md(1,3), md(2,3), me(1,3), me(2,3), mch(1), mch(2), thetaL, thetaR, thetat, thetab, thetatau, greekmu, At, Ab, Atau, 'h'); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     h0amplitudegluongluon = higgsesamplitudedecaygluongluontotal(mh0(1), g, g3atmh0, gp, runmw, alpha, beta, mtPole, mbPole, mcpole, mu(1,3), mu(2,3), md(1,3), md(2,3), thetat, thetab, greekmu, At, Ab, ms, mu(1,2), mu(2,2), md(1,2), md(2,2), Ac, As, mup, mdo, mu(1,1), mu(2,1), md(1,1), md(2,1), Au, Ad, 'h', QCDcorr); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     h0amplitudeZgamma = higgsesamplitudedecayZbosonphotontotal(mh0(1), polemz, g, gp, alphaAtMh, polemw, runmw, alpha, beta, mtAtMh, mbAtMh, mcAtMh, msAtMh, mu(1,3), mu(2,3), md(1,3), md(2,3), mHpm, thetat, thetab, greekmu, At, Ab, 'h'); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW

     
     h0amplitudeneutZ1neutZ1 = higgsphiamplitudedecayneutralinoneutralino (mh0(1), mneut(1), mneut(1), g, tanthetaW, alpha, mixNeut, 1, 1, 'h');
     h0amplitudeneutZ1neutZ2 = higgsphiamplitudedecayneutralinoneutralino (mh0(1), mneut(1), mneut(2), g, tanthetaW, alpha, mixNeut, 1, 2, 'h');
     h0amplitudeneutZ1neutZ3 = higgsphiamplitudedecayneutralinoneutralino (mh0(1), mneut(1), mneut(3), g, tanthetaW, alpha, mixNeut, 1, 3, 'h');
     h0amplitudeneutZ1neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mh0(1), mneut(1), mneut(4), g, tanthetaW, alpha, mixNeut, 1, 4, 'h');
     h0amplitudeneutZ2neutZ2 = higgsphiamplitudedecayneutralinoneutralino (mh0(1), mneut(2), mneut(2), g, tanthetaW, alpha, mixNeut, 2, 2, 'h');
     h0amplitudeneutZ2neutZ3 = higgsphiamplitudedecayneutralinoneutralino (mh0(1), mneut(2), mneut(3), g, tanthetaW, alpha, mixNeut, 2, 3, 'h');
     h0amplitudeneutZ2neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mh0(1), mneut(2), mneut(4), g, tanthetaW, alpha, mixNeut, 2, 4, 'h');
     h0amplitudeneutZ3neutZ3 = higgsphiamplitudedecayneutralinoneutralino (mh0(1), mneut(3), mneut(3), g, tanthetaW, alpha, mixNeut, 3, 3, 'h');
     h0amplitudeneutZ3neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mh0(1), mneut(3), mneut(4), g, tanthetaW, alpha, mixNeut, 3, 4, 'h');
     h0amplitudeneutZ4neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mh0(1), mneut(4), mneut(4), g, tanthetaW, alpha, mixNeut, 4, 4, 'h');
     
     h0amplitudehiggsAhiggsA = higgshamplitudedecayAA (mh0(1), mA0(1), g, gp, alpha, beta, runmw);
     h0amplitudehiggsAZboson = higgshamplitudedecayhiggsAZboson (mh0(1), polemz, mA0(1), g, gp, alpha, beta);
 

   
   }
   
   else if (nmssmIsIt == true) {
     h0amplitudesupLantisupL = higgshamplitudedecay2squarksamehandNMSSM (mh0(1), mu(1,1), mu(1,1), g, gp, alpha, beta, runmw, runmu, CPEMix, 1);
     h0amplitudesupRantisupR = higgshamplitudedecay2squarksamehandNMSSM (mh0(1), mu(2,1), mu(2,1), g, gp, alpha, beta, runmw, runmu, CPEMix, 3);
     h0amplitudesdownLantisdownL = higgshamplitudedecay2squarksamehandNMSSM (mh0(1), md(1,1), md(1,1), g, gp, alpha, beta, runmw, runmd, CPEMix, 2);
     h0amplitudesdownRantisdownR = higgshamplitudedecay2squarksamehandNMSSM (mh0(1), md(2,1), mu(2,1), g, gp, alpha, beta, runmw, runmd, CPEMix, 4);
     h0amplitudescharmLantischarmL = higgshamplitudedecay2squarksamehandNMSSM (mh0(1), mu(1,2), mu(1,2), g, gp, alpha, beta, runmw, runmc, CPEMix, 1);
     h0amplitudescharmRantischarmR = higgshamplitudedecay2squarksamehandNMSSM (mh0(1), mu(2,2), mu(2,2), g, gp, alpha, beta, runmw, runmc, CPEMix, 3);
     h0amplitudesstrangeLantisstrangeL = higgshamplitudedecay2squarksamehandNMSSM (mh0(1), md(2,1), md(2,1), g, gp, alpha, beta, runmw, runms, CPEMix, 2);
     h0amplitudesstrangeRantisstrangeR = higgshamplitudedecay2squarksamehandNMSSM (mh0(1), md(2,2), md(2,2), g, gp, alpha, beta, runmw, runms, CPEMix, 4);
     h0amplitudesnueLantisnueL = higgshamplitudedecay2sleptonsamehandNMSSM (mh0(1), msnu(1), msnu(1), g, gp, alpha, beta, runmw, 0, CPEMix, 1);
     h0amplitudeselectronLantiselectronL = higgshamplitudedecay2sleptonsamehandNMSSM (mh0(1), me(1,1), me(1,1), g, gp, alpha, beta, runmw, runmel, CPEMix, 2);
     h0amplitudeselectronRantiselectronR = higgshamplitudedecay2sleptonsamehandNMSSM (mh0(1), me(2,1), me(2,1), g, gp, alpha, beta, runmw, runmel, CPEMix, 3);
     h0amplitudesnumuLantisnumuL = higgshamplitudedecay2sleptonsamehandNMSSM (mh0(1), msnu(2), msnu(2), g, gp, alpha, beta, runmw, 0, CPEMix, 1);
     h0amplitudesmuonLantismuonL = higgshamplitudedecay2sleptonsamehandNMSSM (mh0(1), me(1,2), me(1,2), g, gp, alpha, beta, runmw, runmmu, CPEMix, 2);
     h0amplitudesmuonRantismuonR = higgshamplitudedecay2sleptonsamehandNMSSM (mh0(1), me(2,2), me(2,2), g, gp, alpha, beta, runmw, runmmu, CPEMix, 3);
     h0amplitudesnutauLantisnutauL = higgshamplitudedecay2sleptonsamehandNMSSM (mh0(1), msnu(3), msnu(3), g, gp, alpha, beta, runmw, 0, CPEMix, 1);
     
     h0amplitudesupLantisupR = higgshamplitudedecay2squarkdiffhandNMSSM(mh0(1), mu(1,1), mu(2,1), g, alpha, beta, runmw, runmu, Au, mueff, lam, CPEMix, 1, 1);
     h0amplitudesupRantisupL = higgshamplitudedecay2squarkdiffhandNMSSM(mh0(1), mu(1,1), mu(2,1), g, alpha, beta, runmw, runmu, Au, mueff, lam, CPEMix, 1, 1);
     h0amplitudesdownLantisdownR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(1), md(1,1), md(2,1), g, alpha, beta, runmw, runmd, Ad, mueff, lam, CPEMix, 2, 1);
     h0amplitudesdownRantisdownL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(1), md(1,1), md(2,1), g, alpha, beta, runmw, runmd, Ad, mueff, lam, CPEMix, 2, 1);
     h0amplitudescharmLantischarmR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(1), mu(1,2), mu(2,2), g, alpha, beta, runmw, runmc, Ac, mueff, lam, CPEMix, 1, 1);
     h0amplitudescharmRantischarmL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(1), mu(1,2), mu(2,2), g, alpha, beta, runmw, runmc, Ac, mueff, lam, CPEMix, 1, 1);
     h0amplitudesstrangeLantisstrangeR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(1), md(2,1), md(2,2), g, alpha, beta, runmw, runms, As, mueff, lam, CPEMix, 2, 1);
     h0amplitudesstrangeRantisstrangeL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(1), md(2,1), md(2,2), g, alpha, beta, runmw, runms, As, mueff, lam, CPEMix, 2, 1);
     h0amplitudeselectronLantiselectronR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(1), me(1,1), me(2,1), g, alpha, beta, runmw, runmel, Ae, mueff, lam, CPEMix, 2, 1)/3;
     h0amplitudeselectronRantiselectronL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(1), me(2,1), me(1,1), g, alpha, beta, runmw, runmel, Ae, mueff, lam, CPEMix, 2, 1)/3;
     h0amplitudesmuonLantismuonR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(1), me(1,2), me(2,2), g, alpha, beta, runmw, runmmu, Amu, mueff, lam, CPEMix, 2, 1)/3;
     h0amplitudesmuonRantismuonL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(1), me(2,2), me(2,1), g, alpha, beta, runmw, runmmu, Amu, mueff, lam, CPEMix, 2, 1)/3;
     
     h0amplitudestop1antistop1 = higgsCPevenamplitudedecaystopistopiNMSSM (mh0(1), mu(1,3), thetat, runmt, g, gp, runmw, beta, CPEMix, At, mueff, lam, 1, 1);
     h0amplitudestop2antistop2 = higgsCPevenamplitudedecaystopistopiNMSSM (mh0(1), mu(2,3), thetat, runmt, g, gp, runmw, beta, CPEMix, At, mueff, lam, 2, 1) ;
     h0amplitudestop1antistop2 = higgsCPevenamplitudedecaystopistopjNMSSM (mh0(1), mu(1,3), mu(2,3), thetat, runmt, g, gp, runmw, beta, CPEMix, At, mueff, lam, 1); 
     h0amplitudestop2antistop1 = h0amplitudestop1antistop2;
     h0amplitudesbottom1antisbottom1 = higgsCPevenamplitudedecaysbottomisbottomiNMSSM (mh0(1), md(1,3), thetab, runmb, g, gp, runmw, beta,  CPEMix, Ab, mueff, lam, 1, 1);
     h0amplitudesbottom2antisbottom2 = higgsCPevenamplitudedecaysbottomisbottomiNMSSM (mh0(1), md(2,3), thetab, runmb, g, gp, runmw, beta,  CPEMix, Ab, mueff, lam, 2, 1);
     h0amplitudesbottom1antisbottom2 = higgsCPevenamplitudedecaysbottomisbottomjNMSSM (mh0(1), md(1,3), md(2,3), thetab, runmb, g, gp, runmw, beta, CPEMix, Ab, mueff, lam, 1);
     h0amplitudesbottom2antisbottom1 = h0amplitudesbottom1antisbottom2;
     h0amplitudestau1antistau1 = higgsCPevenamplitudedecaystauistauiNMSSM (mh0(1), me(1,3), thetatau - PI/2, runmtau, g, gp, runmw, beta,  CPEMix, Atau, mueff, lam, 1, 1); 
     h0amplitudestau2antistau2 = higgsCPevenamplitudedecaystauistauiNMSSM (mh0(1), me(2,3), thetatau - PI/2, runmtau, g, gp, runmw, beta,  CPEMix, Atau, mueff, lam, 2, 1);
     h0amplitudestau1antistau2 = higgsCPevenamplitudedecaystauistaujNMSSM (mh0(1), me(1,3), me(2,3), thetatau - PI/2, runmtau, g, gp, runmw, beta, CPEMix, Atau, mueff, lam, 1);
     h0amplitudestau2antistau1 = h0amplitudestau1antistau2;
     
     h0amplitudecharW1charW1 = higgsphiamplitudedecaysamecharginoNMSSM (mh0(1), mch(1), g, thetaL2, thetaR2, lam, CPEMix, 1, 1);
     h0amplitudecharW2charW2 = higgsphiamplitudedecaysamecharginoNMSSM (mh0(1), mch(2), g, thetaL2, thetaR2, lam, CPEMix, 2, 1);
     h0amplitudecharW1charW2 = higgsphiamplitudedecaydiffcharginoNMSSM (mh0(1), mch(1), mch(2), g, thetaL2, thetaR2, lam, CPEMix, 1);
     h0amplitudegammagamma = higgsCPevenamplitudedecaygammagammaNMSSM(mh0(1), mtAtMh, mbAtMh, mcAtMh, runmtau, runmw, mHpm, mch(1), mch(2), mu(1,2), mu(2,2), mu(1,3), mu(2,3), md(1,2), md(2,2), md(1,3), md(2,3), me(1,2), me(2,2), me(1,3), me(2,3), CPEMix, beta, g, gp, alphaAtMh,thetat, thetab, thetatau-PI/2, thetaL2, thetaR2, At, Ab, Atau, greekmu, mueff, lam, kappa, Alambda, 1); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     h0amplitudegluongluon = higgsCPevenamplitudedecaygluongluonNMSSM(mh0(1), mtPole, mbPole, mcpole, runmw, mu(1,2), mu(2,2), mu(1,3), mu(2,3), md(1,2), md(2,2), md(1,3), md(2,3), mu(1,1), mu(2,1), md(1,1), md(2,1), mtAtMh, mbAtMh, CPEMix, beta, g, gp, gs, alphasAtMh, thetat, thetab, thetaL2, thetaR2, At, Ab, greekmu, mueff, lam, kappa, Alambda, 1, QCDcorr); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     h0amplitudeZgamma = higgshamplitudedecayZgammaNMSSM (mh0(1), g, gp, alphaAtMh, runmw, polemz, mHpm, CPEMix, beta, mtAtMh, mbAtMh, mcAtMh, mch(1), mch(2), thetaL2, thetaR2, lam, kappa, Alambda, greekmu, mueff, 1); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     
     h0amplitudeneutZ1neutZ1 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(1), mneut(1), g, gp, CPEMix, mixNeut, lam, kappa, 1, 1, 1);
     h0amplitudeneutZ1neutZ2 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(1), mneut(2), g, gp, CPEMix, mixNeut, lam, kappa, 1, 2, 1);
     h0amplitudeneutZ1neutZ3 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(1), mneut(3), g, gp, CPEMix, mixNeut, lam, kappa, 1, 3, 1);
     h0amplitudeneutZ1neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(1), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 1, 4, 1);
     h0amplitudeneutZ1neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(1), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 1, 5, 1);
     h0amplitudeneutZ2neutZ2 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(2), mneut(2), g, gp, CPEMix, mixNeut, lam, kappa, 2, 2, 1);
     h0amplitudeneutZ2neutZ3 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(2), mneut(3), g, gp, CPEMix, mixNeut, lam, kappa, 2, 3, 1);
     h0amplitudeneutZ2neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(2), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 2, 4, 1);
     h0amplitudeneutZ2neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(2), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 2, 5, 1);
     h0amplitudeneutZ3neutZ3 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(3), mneut(3), g, gp, CPEMix, mixNeut, lam, kappa, 3, 3, 1);
     h0amplitudeneutZ3neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(3), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 3, 4, 1);
     h0amplitudeneutZ3neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(3), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 3, 5, 1);
     h0amplitudeneutZ4neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(4), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 4, 4, 1);
     h0amplitudeneutZ4neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(4), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 4, 5, 1);
     h0amplitudeneutZ5neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(1), mneut(5), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 5, 5, 1);
     
     h0amplitudehiggsAhiggsA = higgsCPevenamplitudedecayAANMSSM(mh0(1), mA0(1), mA0(1), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 1, 1, 1);
     h0amplitudehiggsAhiggsA2 = 2*higgsCPevenamplitudedecayAANMSSM(mh0(1), mA0(1), mA0(2), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 1, 1, 2);
     h0amplitudehiggsA2higgsA2 = higgsCPevenamplitudedecayAANMSSM(mh0(1), mA0(2), mA0(2), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 1, 2, 2);
     h0amplitudehiggsAZboson = higgsCPevenamplitudedecaypseudoscalarZNMSSM(mh0(1), mA0(1), polemz, g, gp, beta, CPEMix, CPOMix, 1, 1);
     h0amplitudehiggsA2Zboson = higgsCPevenamplitudedecaypseudoscalarZNMSSM(mh0(1), mA0(2), polemz, g, gp, beta, CPEMix, CPOMix, 1, 2);
     h0amplitudeHpmHpm = higgsCPevenamplitudedecayHpHmNMSSM (mh0(1), mHpm, runmw, g, gp, runmt, runmb, beta, lam, mueff, kappa, Alambda, CPEMix, 1); ///Note not included for MSSM as h0 is always lighter than Hpm in that case
     h0amplitudeWHpm = higgsCPevenamplitudedecayWHpmNMSSM (mh0(1), polemw, mHpm, beta, g, CPEMix, 1)*2; ///*2 as W+H- or W-H+
   }
   
   h0amplitudeWW = higgshamplitudedecayVV(mh0(1), polemw, polemz, g, gp, alpha, beta, 'W', CPEMix, nmssmIsIt)(1);
   h0amplitudeZZ = higgshamplitudedecayVV(mh0(1), polemw, polemz, g, gp, alpha, beta, 'Z', CPEMix, nmssmIsIt)(1);
   
   int h0WWcommentcode, h0ZZcommentcode, h0WWNDA=0, h0ZZNDA=0;
   h0WWcommentcode = higgshamplitudedecayVV(mh0(1), polemw, polemz, g, gp, alpha, beta, 'W', CPEMix, nmssmIsIt)(2);
   string h0WWcomment, h0ZZcomment;
   if (h0WWcommentcode == 1) {
     h0WWcomment = "# h -> WW* -> W f f'bar";
     h0WWNDA = 2; ///So read into other programs, e.g. PYTHIA, correctly
   }
   else if(h0WWcommentcode == 2) {
     h0WWcomment = "# h -> W+ W-";
     h0WWNDA = 2;
   }
   h0ZZcommentcode = higgshamplitudedecayVV(mh0(1), polemw, polemz, g, gp, alpha, beta, 'Z', CPEMix, nmssmIsIt)(2);
   if (h0ZZcommentcode == 1) {
     h0ZZcomment = "# h -> ZZ* -> Z f f'bar";
     h0ZZNDA = 2; ///So read into other programs, e.g. PYTHIA, correctly
   }
   else if(h0ZZcommentcode == 2) {
     h0ZZcomment = "# h -> Z Z";
     h0ZZNDA = 2;
   }
   
   Particlehiggsl.Array_Decays[0][0] = PDGup; Particlehiggsl.Array_Decays[0][1] = -PDGup; Particlehiggsl.Array_Decays[0][2] = h0amplitudeuantiu; Particlehiggsl.Array_Decays[0][3] = 2; Particlehiggsl.Array_Comments[0] = "# h -> u ub";
   Particlehiggsl.Array_Decays[1][0] = PDGdown; Particlehiggsl.Array_Decays[1][1] = -PDGdown; Particlehiggsl.Array_Decays[1][2] = h0amplitudedantid; Particlehiggsl.Array_Decays[1][3] = 2; Particlehiggsl.Array_Comments[1] = "# h -> d db";
   Particlehiggsl.Array_Decays[2][0] = PDGcharm; Particlehiggsl.Array_Decays[2][1] = -PDGcharm; Particlehiggsl.Array_Decays[2][2] = h0amplitudecantic; Particlehiggsl.Array_Decays[2][3] = 2; Particlehiggsl.Array_Comments[2] = "# h -> c cb";
   Particlehiggsl.Array_Decays[3][0] = PDGstrange; Particlehiggsl.Array_Decays[3][1] = -PDGstrange; Particlehiggsl.Array_Decays[3][2] = h0amplitudesantis; Particlehiggsl.Array_Decays[3][3] = 2; Particlehiggsl.Array_Comments[3] = "# h -> s sb";
   Particlehiggsl.Array_Decays[4][0] = PDGbottom; Particlehiggsl.Array_Decays[4][1] = -PDGbottom; Particlehiggsl.Array_Decays[4][2] = h0amplitudebantib; Particlehiggsl.Array_Decays[4][3] = 2; Particlehiggsl.Array_Comments[4] = "# h -> b bb";
   Particlehiggsl.Array_Decays[5][0] = PDGelectron; Particlehiggsl.Array_Decays[5][1] = -PDGelectron; Particlehiggsl.Array_Decays[5][2] = h0amplitudeeantie; Particlehiggsl.Array_Decays[5][3] = 2; Particlehiggsl.Array_Comments[5] = "# h -> e- e+";
   Particlehiggsl.Array_Decays[6][0] = PDGmuon; Particlehiggsl.Array_Decays[6][1] = -PDGmuon; Particlehiggsl.Array_Decays[6][2] = h0amplitudemuantimu; Particlehiggsl.Array_Decays[6][3] = 2; Particlehiggsl.Array_Comments[6] = "# h -> mu- mu+";
   Particlehiggsl.Array_Decays[7][0] = PDGtau; Particlehiggsl.Array_Decays[7][1] = -PDGtau; Particlehiggsl.Array_Decays[7][2] = h0amplitudetauantitau; Particlehiggsl.Array_Decays[7][3] = 2; Particlehiggsl.Array_Comments[7] = "# h-> tau- tau+";
   Particlehiggsl.Array_Decays[8][0] = PDGneutralino1; Particlehiggsl.Array_Decays[8][1] = PDGneutralino1; Particlehiggsl.Array_Decays[8][2] = h0amplitudeneutZ1neutZ1; Particlehiggsl.Array_Decays[8][3] = 2; Particlehiggsl.Array_Comments[8] = "# h -> ~chi_10 ~chi_10";
   Particlehiggsl.Array_Decays[9][0] = PDGneutralino1; Particlehiggsl.Array_Decays[9][1] = PDGneutralino2; Particlehiggsl.Array_Decays[9][2] = h0amplitudeneutZ1neutZ2; Particlehiggsl.Array_Decays[9][3] = 2; Particlehiggsl.Array_Comments[9] = "# h -> ~chi_10 ~chi_20";
   Particlehiggsl.Array_Decays[10][0] = PDGneutralino1; Particlehiggsl.Array_Decays[10][1] = PDGneutralino3; Particlehiggsl.Array_Decays[10][2] = h0amplitudeneutZ1neutZ3; Particlehiggsl.Array_Decays[10][3] = 2; Particlehiggsl.Array_Comments[10] = "# h -> ~chi_10 ~chi_30";
   Particlehiggsl.Array_Decays[11][0] = PDGneutralino1; Particlehiggsl.Array_Decays[11][1] = PDGneutralino4; Particlehiggsl.Array_Decays[11][2] = h0amplitudeneutZ1neutZ4; Particlehiggsl.Array_Decays[11][3] = 2; Particlehiggsl.Array_Comments[11] = "# h -> ~chi_10 ~chi_40";
   Particlehiggsl.Array_Decays[12][0] = PDGneutralino2; Particlehiggsl.Array_Decays[12][1] = PDGneutralino2; Particlehiggsl.Array_Decays[12][2] = h0amplitudeneutZ2neutZ2; Particlehiggsl.Array_Decays[12][3] = 2; Particlehiggsl.Array_Comments[12] = "# h -> ~chi_20 ~chi_20";
   Particlehiggsl.Array_Decays[13][0] = PDGneutralino2; Particlehiggsl.Array_Decays[13][1] = PDGneutralino3; Particlehiggsl.Array_Decays[13][2] = h0amplitudeneutZ2neutZ3; Particlehiggsl.Array_Decays[13][3] = 2; Particlehiggsl.Array_Comments[13] = "# h -> ~chi_20 ~chi_30";
   Particlehiggsl.Array_Decays[14][0] = PDGneutralino2; Particlehiggsl.Array_Decays[14][1] = PDGneutralino4; Particlehiggsl.Array_Decays[14][2] = h0amplitudeneutZ2neutZ4; Particlehiggsl.Array_Decays[14][3] = 2; Particlehiggsl.Array_Comments[14] = "# h -> ~chi_20 ~chi_40";
   Particlehiggsl.Array_Decays[15][0] = PDGneutralino3; Particlehiggsl.Array_Decays[15][1] = PDGneutralino3; Particlehiggsl.Array_Decays[15][2] = h0amplitudeneutZ3neutZ3; Particlehiggsl.Array_Decays[15][3] = 2; Particlehiggsl.Array_Comments[15] = "# h -> ~chi_30 ~chi_30";
   Particlehiggsl.Array_Decays[16][0] = PDGneutralino3; Particlehiggsl.Array_Decays[16][1] = PDGneutralino4; Particlehiggsl.Array_Decays[16][2] = h0amplitudeneutZ3neutZ4; Particlehiggsl.Array_Decays[16][3] = 2; Particlehiggsl.Array_Comments[16] = "# h -> ~chi_30 ~chi_40";
   Particlehiggsl.Array_Decays[17][0] = PDGneutralino4; Particlehiggsl.Array_Decays[17][1] = PDGneutralino4; Particlehiggsl.Array_Decays[17][2] = h0amplitudeneutZ4neutZ4; Particlehiggsl.Array_Decays[17][3] = 2; Particlehiggsl.Array_Comments[17] = "# h -> ~chi_40 ~chi_40";
   Particlehiggsl.Array_Decays[18][0] = PDGchargino1; Particlehiggsl.Array_Decays[18][1] = -PDGchargino1; Particlehiggsl.Array_Decays[18][2] = h0amplitudecharW1charW1; Particlehiggsl.Array_Decays[18][3] = 2; Particlehiggsl.Array_Comments[18] = "# h -> ~chi_1+ ~chi_1-";
   Particlehiggsl.Array_Decays[19][0] = PDGchargino2; Particlehiggsl.Array_Decays[19][1] = -PDGchargino2; Particlehiggsl.Array_Decays[19][2] = h0amplitudecharW2charW2; Particlehiggsl.Array_Decays[19][3] = 2; Particlehiggsl.Array_Comments[19] = "# h -> ~chi_2+ ~chi_2-";
   Particlehiggsl.Array_Decays[20][0] = PDGchargino1; Particlehiggsl.Array_Decays[20][1] = -PDGchargino2; Particlehiggsl.Array_Decays[20][2] = h0amplitudecharW1charW2; Particlehiggsl.Array_Decays[20][3] = 2; Particlehiggsl.Array_Comments[20] = "# h -> ~chi_1+ ~chi_2-";
   Particlehiggsl.Array_Decays[21][0] = PDGchargino2; Particlehiggsl.Array_Decays[21][1] = -PDGchargino1; Particlehiggsl.Array_Decays[21][2] = h0amplitudecharW1charW2; Particlehiggsl.Array_Decays[21][3] = 2; Particlehiggsl.Array_Comments[21] = "# h -> ~chi_2+ ~chi_1-"; ///amplitude same as decay to ~chi1+ and ~chi2- by CP invariance
   Particlehiggsl.Array_Decays[22][0] = PDGA0; Particlehiggsl.Array_Decays[22][1] = PDGA0; Particlehiggsl.Array_Decays[22][2] = h0amplitudehiggsAhiggsA; Particlehiggsl.Array_Decays[22][3] = 2; Particlehiggsl.Array_Comments[22] = "# h -> A A";
   Particlehiggsl.Array_Decays[23][0] = PDGZboson; Particlehiggsl.Array_Decays[23][1] = PDGA0; Particlehiggsl.Array_Decays[23][2] = h0amplitudehiggsAZboson; Particlehiggsl.Array_Decays[23][3] = 2; Particlehiggsl.Array_Comments[23] = "# h -> A Z";
   Particlehiggsl.Array_Decays[24][0] = PDGsupL; Particlehiggsl.Array_Decays[24][1] = -PDGsupL; Particlehiggsl.Array_Decays[24][2] = h0amplitudesupLantisupL; Particlehiggsl.Array_Decays[24][3] = 2; Particlehiggsl.Array_Comments[24] = "# h -> ~u_L ~u_L*";
   Particlehiggsl.Array_Decays[25][0] = PDGsupR; Particlehiggsl.Array_Decays[25][1] = -PDGsupR; Particlehiggsl.Array_Decays[25][2] = h0amplitudesupRantisupR; Particlehiggsl.Array_Decays[25][3] = 2; Particlehiggsl.Array_Comments[25] = "# h -> ~u_R ~u_R*";
   Particlehiggsl.Array_Decays[26][0] = PDGsupL; Particlehiggsl.Array_Decays[26][1] = -PDGsupR; Particlehiggsl.Array_Decays[26][2] = h0amplitudesupLantisupR; Particlehiggsl.Array_Decays[26][3] = 2; Particlehiggsl.Array_Comments[26] = "# h -> ~u_L ~u_R*";
   Particlehiggsl.Array_Decays[27][0] = PDGsupR; Particlehiggsl.Array_Decays[27][1] = -PDGsupL; Particlehiggsl.Array_Decays[27][2] = h0amplitudesupRantisupL; Particlehiggsl.Array_Decays[27][3] = 2; Particlehiggsl.Array_Comments[27] = "# h -> ~u_R ~u_L*";
   Particlehiggsl.Array_Decays[28][0] = PDGsdownL; Particlehiggsl.Array_Decays[28][1] = -PDGsdownL; Particlehiggsl.Array_Decays[28][2] = h0amplitudesdownLantisdownL; Particlehiggsl.Array_Decays[28][3] = 2; Particlehiggsl.Array_Comments[28] = "# h -> ~d_L ~d_L*";
   Particlehiggsl.Array_Decays[29][0] = PDGsdownR; Particlehiggsl.Array_Decays[29][1] = -PDGsdownR; Particlehiggsl.Array_Decays[29][2] = h0amplitudesdownRantisdownR; Particlehiggsl.Array_Decays[29][3] = 2; Particlehiggsl.Array_Comments[29] = "# h -> ~d_R ~d_R*";
   Particlehiggsl.Array_Decays[30][0] = PDGsdownL; Particlehiggsl.Array_Decays[30][1] = -PDGsdownR; Particlehiggsl.Array_Decays[30][2] = h0amplitudesdownLantisdownR; Particlehiggsl.Array_Decays[30][3] = 2; Particlehiggsl.Array_Comments[30] = "# h -> ~d_L ~d_R*";
   Particlehiggsl.Array_Decays[31][0] = PDGsdownR; Particlehiggsl.Array_Decays[31][1] = -PDGsdownL; Particlehiggsl.Array_Decays[31][2] = h0amplitudesdownRantisdownL; Particlehiggsl.Array_Decays[31][3] = 2; Particlehiggsl.Array_Comments[31] = "# h -> ~d_R ~d_L*";
   Particlehiggsl.Array_Decays[32][0] = PDGscharmL; Particlehiggsl.Array_Decays[32][1] = -PDGscharmL; Particlehiggsl.Array_Decays[32][2] = h0amplitudescharmLantischarmL; Particlehiggsl.Array_Decays[32][3] = 2; Particlehiggsl.Array_Comments[32] = "# h -> ~c_L ~c_L*";
   Particlehiggsl.Array_Decays[33][0] = PDGscharmR; Particlehiggsl.Array_Decays[33][1] = -PDGscharmR; Particlehiggsl.Array_Decays[33][2] = h0amplitudescharmRantischarmR; Particlehiggsl.Array_Decays[33][3] = 2; Particlehiggsl.Array_Comments[33] = "# h -> ~c_R ~c_R*";
   Particlehiggsl.Array_Decays[34][0] = PDGscharmL; Particlehiggsl.Array_Decays[34][1] = -PDGscharmR; Particlehiggsl.Array_Decays[34][2] = h0amplitudescharmLantischarmR; Particlehiggsl.Array_Decays[34][3] = 2; Particlehiggsl.Array_Comments[34] = "# h -> ~c_L ~c_R*";
   Particlehiggsl.Array_Decays[35][0] = PDGscharmR; Particlehiggsl.Array_Decays[35][1] = -PDGscharmL; Particlehiggsl.Array_Decays[35][2] = h0amplitudescharmRantischarmL; Particlehiggsl.Array_Decays[35][3] = 2; Particlehiggsl.Array_Comments[35] = "# h -> ~c_R ~c_L*";
   Particlehiggsl.Array_Decays[36][0] = PDGsstrangeL; Particlehiggsl.Array_Decays[36][1] = -PDGsstrangeL; Particlehiggsl.Array_Decays[36][2] = h0amplitudesstrangeLantisstrangeL; Particlehiggsl.Array_Decays[36][3] = 2; Particlehiggsl.Array_Comments[36] = "# h -> ~s_L ~s_L*";
   Particlehiggsl.Array_Decays[37][0] = PDGsstrangeR; Particlehiggsl.Array_Decays[37][1] = -PDGsstrangeR; Particlehiggsl.Array_Decays[37][2] = h0amplitudesstrangeRantisstrangeR; Particlehiggsl.Array_Decays[37][3] = 2; Particlehiggsl.Array_Comments[37] = "# h -> ~s_R ~s_R*";
   Particlehiggsl.Array_Decays[38][0] = PDGsstrangeL; Particlehiggsl.Array_Decays[38][1] = -PDGsstrangeR; Particlehiggsl.Array_Decays[38][2] = h0amplitudesstrangeLantisstrangeR; Particlehiggsl.Array_Decays[38][3] = 2; Particlehiggsl.Array_Comments[38] = "# h -> ~s_L ~s_R*";
   Particlehiggsl.Array_Decays[39][0] = PDGsstrangeR; Particlehiggsl.Array_Decays[39][1] = -PDGsstrangeL; Particlehiggsl.Array_Decays[39][2] = h0amplitudesstrangeRantisstrangeL; Particlehiggsl.Array_Decays[39][3] = 2; Particlehiggsl.Array_Comments[39] = "# h -> ~s_R ~s_L*";
   Particlehiggsl.Array_Decays[40][0] = PDGnuselectronL; Particlehiggsl.Array_Decays[40][1] = -PDGnuselectronL; Particlehiggsl.Array_Decays[40][2] = h0amplitudesnueLantisnueL; Particlehiggsl.Array_Decays[40][3] = 2; Particlehiggsl.Array_Comments[40] = "# h -> ~nu_eL ~nu_eL*";
   Particlehiggsl.Array_Decays[41][0] = PDGselectronL; Particlehiggsl.Array_Decays[41][1] = -PDGselectronL; Particlehiggsl.Array_Decays[41][2] = h0amplitudeselectronLantiselectronL; Particlehiggsl.Array_Decays[41][3] = 2; Particlehiggsl.Array_Comments[41] = "# h -> ~e_L- ~e_L+";
   Particlehiggsl.Array_Decays[42][0] = PDGselectronR; Particlehiggsl.Array_Decays[42][1] = -PDGselectronR; Particlehiggsl.Array_Decays[42][2] = h0amplitudeselectronRantiselectronR; Particlehiggsl.Array_Decays[42][3] = 2; Particlehiggsl.Array_Comments[42] = "# h -> ~e_R- ~e_R+";
   Particlehiggsl.Array_Decays[43][0] = PDGselectronL; Particlehiggsl.Array_Decays[43][1] = -PDGselectronR; Particlehiggsl.Array_Decays[43][2] = h0amplitudeselectronLantiselectronR; Particlehiggsl.Array_Decays[43][3] = 2; Particlehiggsl.Array_Comments[43] = "# h -> ~e_L- ~e_R+";
   Particlehiggsl.Array_Decays[44][0] = PDGselectronR; Particlehiggsl.Array_Decays[44][1] = -PDGselectronL; Particlehiggsl.Array_Decays[44][2] = h0amplitudeselectronRantiselectronL; Particlehiggsl.Array_Decays[44][3] = 2; Particlehiggsl.Array_Comments[44] = "# h -> ~e_R- ~e_L+";
   Particlehiggsl.Array_Decays[45][0] = PDGnusmuonL; Particlehiggsl.Array_Decays[45][1] = -PDGnusmuonL; Particlehiggsl.Array_Decays[45][2] = h0amplitudesnumuLantisnumuL; Particlehiggsl.Array_Decays[45][3] = 2; Particlehiggsl.Array_Comments[45] = "# h -> ~nu_muL ~nu_muL*";
   Particlehiggsl.Array_Decays[46][0] = PDGsmuonL; Particlehiggsl.Array_Decays[46][1] = -PDGsmuonL; Particlehiggsl.Array_Decays[46][2] = h0amplitudesmuonLantismuonL; Particlehiggsl.Array_Decays[46][3] = 2; Particlehiggsl.Array_Comments[46] = "# h -> ~mu_L- ~mu_L+";
   Particlehiggsl.Array_Decays[47][0] = PDGsmuonR; Particlehiggsl.Array_Decays[47][1] = -PDGsmuonR; Particlehiggsl.Array_Decays[47][2] = h0amplitudesmuonRantismuonR; Particlehiggsl.Array_Decays[47][3] = 2; Particlehiggsl.Array_Comments[47] = "# h -> ~mu_R- ~mu_R+";
   Particlehiggsl.Array_Decays[48][0] = PDGsmuonL; Particlehiggsl.Array_Decays[48][1] = -PDGsmuonR; Particlehiggsl.Array_Decays[48][2] = h0amplitudesmuonLantismuonR; Particlehiggsl.Array_Decays[48][3] = 2; Particlehiggsl.Array_Comments[48] = "# h -> ~mu_L- ~mu_R+";
   Particlehiggsl.Array_Decays[49][0] = PDGsmuonR; Particlehiggsl.Array_Decays[49][1] = -PDGsmuonL; Particlehiggsl.Array_Decays[49][2] = h0amplitudesmuonRantismuonL; Particlehiggsl.Array_Decays[49][3] = 2; Particlehiggsl.Array_Comments[49] = "# h -> ~mu_R- ~mu_L+";	
   Particlehiggsl.Array_Decays[50][0] = PDGnustauL; Particlehiggsl.Array_Decays[50][1] = -PDGnustauL; Particlehiggsl.Array_Decays[50][2] = h0amplitudesnutauLantisnutauL; Particlehiggsl.Array_Decays[50][3] = 2; Particlehiggsl.Array_Comments[50] = "# h -> ~nu_tauL ~nu_tauL*";   Particlehiggsl.Array_Decays[51][0] = PDGstop1; Particlehiggsl.Array_Decays[51][1] = -PDGstop1; Particlehiggsl.Array_Decays[51][2] = h0amplitudestop1antistop1; Particlehiggsl.Array_Decays[51][3] = 2; Particlehiggsl.Array_Comments[51] = "# h -> ~t_1 ~t_1*";	  
   Particlehiggsl.Array_Decays[52][0] = PDGstop2; Particlehiggsl.Array_Decays[52][1] = -PDGstop2; Particlehiggsl.Array_Decays[52][2] = h0amplitudestop2antistop2; Particlehiggsl.Array_Decays[52][3] = 2; Particlehiggsl.Array_Comments[52] = "# h -> ~t_2 ~t_2*";
   Particlehiggsl.Array_Decays[53][0] = PDGstop1; Particlehiggsl.Array_Decays[53][1] = -PDGstop2; Particlehiggsl.Array_Decays[53][2] = h0amplitudestop1antistop2; Particlehiggsl.Array_Decays[53][3] = 2; Particlehiggsl.Array_Comments[53] = "# h -> ~t_1 ~t_2*";	  
   Particlehiggsl.Array_Decays[54][0] = PDGstop2; Particlehiggsl.Array_Decays[54][1] = -PDGstop1; Particlehiggsl.Array_Decays[54][2] = h0amplitudestop2antistop1; Particlehiggsl.Array_Decays[54][3] = 2; Particlehiggsl.Array_Comments[54] = "# h -> ~t_2 ~t_1*";
   Particlehiggsl.Array_Decays[55][0] = PDGsbottom1; Particlehiggsl.Array_Decays[55][1] = -PDGsbottom1; Particlehiggsl.Array_Decays[55][2] = h0amplitudesbottom1antisbottom1; Particlehiggsl.Array_Decays[55][3] = 2; Particlehiggsl.Array_Comments[55] = "# h -> ~b_1 ~b_1*";	  
   Particlehiggsl.Array_Decays[56][0] = PDGsbottom2; Particlehiggsl.Array_Decays[56][1] = -PDGsbottom2; Particlehiggsl.Array_Decays[56][2] = h0amplitudesbottom2antisbottom2; Particlehiggsl.Array_Decays[56][3] = 2; Particlehiggsl.Array_Comments[56] = "# h -> ~b_2 ~b_2*";
   Particlehiggsl.Array_Decays[57][0] = PDGsbottom1; Particlehiggsl.Array_Decays[57][1] = -PDGsbottom2; Particlehiggsl.Array_Decays[57][2] = h0amplitudesbottom1antisbottom2; Particlehiggsl.Array_Decays[57][3] = 2; Particlehiggsl.Array_Comments[57] = "# h -> ~b_1 ~b_2*";	  
   Particlehiggsl.Array_Decays[58][0] = PDGsbottom2; Particlehiggsl.Array_Decays[58][1] = -PDGsbottom1; Particlehiggsl.Array_Decays[58][2] = h0amplitudesbottom2antisbottom1; Particlehiggsl.Array_Decays[58][3] = 2; Particlehiggsl.Array_Comments[58] = "# h -> ~b_2 ~b_1*";
   Particlehiggsl.Array_Decays[59][0] = PDGstau1; Particlehiggsl.Array_Decays[59][1] = -PDGstau1; Particlehiggsl.Array_Decays[59][2] = h0amplitudestau1antistau1; Particlehiggsl.Array_Decays[59][3] = 2; Particlehiggsl.Array_Comments[59] = "# h -> ~tau_1- ~tau_1+";
   Particlehiggsl.Array_Decays[60][0] = PDGstau2; Particlehiggsl.Array_Decays[60][1] = -PDGstau2; Particlehiggsl.Array_Decays[60][2] = h0amplitudestau2antistau2; Particlehiggsl.Array_Decays[60][3] = 2; Particlehiggsl.Array_Comments[60] = "# h -> ~tau_2- ~tau_2+";
   Particlehiggsl.Array_Decays[61][0] = PDGstau1; Particlehiggsl.Array_Decays[61][1] = -PDGstau2; Particlehiggsl.Array_Decays[61][2] = h0amplitudestau1antistau2; Particlehiggsl.Array_Decays[61][3] = 2; Particlehiggsl.Array_Comments[61] = "# h -> ~tau_1- ~tau_2+";
   Particlehiggsl.Array_Decays[62][0] = PDGstau2; Particlehiggsl.Array_Decays[62][1] = -PDGstau1; Particlehiggsl.Array_Decays[62][2] = h0amplitudestau2antistau1; Particlehiggsl.Array_Decays[62][3] = 2; Particlehiggsl.Array_Comments[62] = "# h -> ~tau_2- ~tau_1+";
   Particlehiggsl.Array_Decays[63][0] = PDGphoton; Particlehiggsl.Array_Decays[63][1] = PDGphoton; Particlehiggsl.Array_Decays[63][2] = h0amplitudegammagamma; Particlehiggsl.Array_Decays[63][3] = 2; Particlehiggsl.Array_Comments[63] = "# h -> gamma gamma";
   Particlehiggsl.Array_Decays[64][0] = PDGgluon; Particlehiggsl.Array_Decays[64][1] = PDGgluon; Particlehiggsl.Array_Decays[64][2] = h0amplitudegluongluon; Particlehiggsl.Array_Decays[64][3] = 2; Particlehiggsl.Array_Comments[64] = "# h -> gluon gluon";
   Particlehiggsl.Array_Decays[65][0] = PDGZboson; Particlehiggsl.Array_Decays[65][1] = PDGphoton; Particlehiggsl.Array_Decays[65][2] = h0amplitudeZgamma; Particlehiggsl.Array_Decays[65][3] = 2; Particlehiggsl.Array_Comments[65] = "# h -> Z gamma";
   Particlehiggsl.Array_Decays[66][0] = PDGWplus; Particlehiggsl.Array_Decays[66][1] = -PDGWplus; Particlehiggsl.Array_Decays[66][2] = h0amplitudeWW; Particlehiggsl.Array_Decays[66][3] = h0WWNDA; Particlehiggsl.Array_Comments[66] = h0WWcomment;
   Particlehiggsl.Array_Decays[67][0] = PDGZboson; Particlehiggsl.Array_Decays[67][1] = PDGZboson; Particlehiggsl.Array_Decays[67][2] = h0amplitudeZZ; Particlehiggsl.Array_Decays[67][3] = h0ZZNDA; Particlehiggsl.Array_Comments[67] = h0ZZcomment;
   
   Particlehiggsl.Array_Decays[68][0] = PDGneutralino1; Particlehiggsl.Array_Decays[68][1] = PDGneutralino5; Particlehiggsl.Array_Decays[68][2] = h0amplitudeneutZ1neutZ5; Particlehiggsl.Array_Decays[68][3] = 2; Particlehiggsl.Array_Comments[68] = "# h -> ~chi_10 ~chi_50";
   Particlehiggsl.Array_Decays[69][0] = PDGneutralino1; Particlehiggsl.Array_Decays[69][1] = PDGneutralino5; Particlehiggsl.Array_Decays[69][2] = h0amplitudeneutZ2neutZ5; Particlehiggsl.Array_Decays[69][3] = 2; Particlehiggsl.Array_Comments[69] = "# h -> ~chi_20 ~chi_50";
   Particlehiggsl.Array_Decays[70][0] = PDGneutralino1; Particlehiggsl.Array_Decays[70][1] = PDGneutralino5; Particlehiggsl.Array_Decays[70][2] = h0amplitudeneutZ3neutZ5; Particlehiggsl.Array_Decays[70][3] = 2; Particlehiggsl.Array_Comments[70] = "# h -> ~chi_30 ~chi_50";
   Particlehiggsl.Array_Decays[71][0] = PDGneutralino1; Particlehiggsl.Array_Decays[71][1] = PDGneutralino5; Particlehiggsl.Array_Decays[71][2] = h0amplitudeneutZ4neutZ5; Particlehiggsl.Array_Decays[71][3] = 2; Particlehiggsl.Array_Comments[71] = "# h -> ~chi_40 ~chi_50";
   Particlehiggsl.Array_Decays[72][0] = PDGneutralino1; Particlehiggsl.Array_Decays[72][1] = PDGneutralino5; Particlehiggsl.Array_Decays[72][2] = h0amplitudeneutZ5neutZ5; Particlehiggsl.Array_Decays[72][3] = 2; Particlehiggsl.Array_Comments[72] = "# h -> ~chi_50 ~chi_50";
   Particlehiggsl.Array_Decays[73][0] = PDGA0; Particlehiggsl.Array_Decays[73][1] = PDGA2; Particlehiggsl.Array_Decays[73][2] = h0amplitudehiggsAhiggsA2; Particlehiggsl.Array_Decays[73][3] = 2; Particlehiggsl.Array_Comments[73] = "# h -> A A2";

   Particlehiggsl.Array_Decays[75][0] = PDGA2; Particlehiggsl.Array_Decays[75][1] = PDGA2; Particlehiggsl.Array_Decays[75][2] = h0amplitudehiggsA2higgsA2; Particlehiggsl.Array_Decays[75][3] = 2; Particlehiggsl.Array_Comments[75] = "# h -> A2 A2";
   Particlehiggsl.Array_Decays[76][0] = PDGZboson; Particlehiggsl.Array_Decays[76][1] = PDGA2; Particlehiggsl.Array_Decays[76][2] = h0amplitudehiggsA2Zboson; Particlehiggsl.Array_Decays[76][3] = 2; Particlehiggsl.Array_Comments[76] = "# h -> A2 Z";
   Particlehiggsl.Array_Decays[77][0] = PDGHplus; Particlehiggsl.Array_Decays[77][1] = -PDGHplus; Particlehiggsl.Array_Decays[77][2] = h0amplitudeHpmHpm; Particlehiggsl.Array_Decays[77][3] = 2; Particlehiggsl.Array_Comments[77] = "# h -> H+ H-";
   Particlehiggsl.Array_Decays[78][0] = PDGHplus; Particlehiggsl.Array_Decays[78][1] = -PDGWplus; Particlehiggsl.Array_Decays[78][2] = h0amplitudeWHpm; Particlehiggsl.Array_Decays[78][3] = 2; Particlehiggsl.Array_Comments[78] = "# h -> H+- W-+";

   for(int i = 0; i<Particlehiggsl.No_of_Decays; i++) {
     if (Particlehiggsl.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << Particlehiggsl.Array_Comments[i] << " is negative = " << Particlehiggsl.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       Particlehiggsl.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }   
   
   double higgsl_No_1to2_Decays = 0;
   
   higgsl_No_1to2_Decays = Particlehiggsl.No_1to2_Decays + Particlehiggsl.No_NMSSM_Decays; /// No higgs NLSP decays to gravitinos as they are unimportant as would be swamped by
   
   for (int j = 0; j<higgsl_No_1to2_Decays; j++) {
     Particlehiggsl.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<higgsl_No_1to2_Decays; j++) {
     Particlehiggsl.two_width = Particlehiggsl.two_width + Particlehiggsl.Array_Decays[j][2];
   }
   for (int j=higgsl_No_1to2_Decays; j<Particlehiggsl.No_of_Decays; j++) {
     Particlehiggsl.three_width = Particlehiggsl.three_width + Particlehiggsl.Array_Decays[j][2];
   }
   
   for(int j=0; j<Particlehiggsl.No_of_Decays; j++) {
     Particlehiggsl.Array_Decays[j][4] = 0;
   }
   
   ///Could argue no need for test for nans here as the higgs 1 -> 3 decay formulae are all purely analytic algebraic expressions, therefore no numerical integration is involved so we can't get nans. Will check anyway as possibility of -ve sqrts in kinematics or -ve logs, or infs etc
   if ( Particlehiggsl.three_width != Particlehiggsl.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for higgsl - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       Particlehiggsl.No_of_Decays = higgsl_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       Particlehiggsl.total_width = Particlehiggsl.two_width;
     }
   else {
     Particlehiggsl.total_width = Particlehiggsl.two_width + Particlehiggsl.three_width;
   }
   
   if ( Particlehiggsl.total_width != Particlehiggsl.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<Particlehiggsl.No_of_Decays; i++) {
       //   fout << Particlehiggsl.Array_Decays[i][2] << endl;
       // }
       throw( "nan in lightest higgs total width \n");
     }
   
 }
 ///higgsH decays
  double H0amplitudeuantiu=0, H0amplitudedantid=0, H0amplitudesantis=0, H0amplitudecantic=0, H0amplitudebantib=0, H0amplitudetantit=0, H0amplitudeeantie=0, H0amplitudemuantimu=0, H0amplitudetauantitau=0, H0amplitudeneutZ1neutZ1=0, H0amplitudeneutZ1neutZ2=0, H0amplitudeneutZ1neutZ3=0, H0amplitudeneutZ1neutZ4=0, H0amplitudeneutZ2neutZ2=0, H0amplitudeneutZ2neutZ3=0, H0amplitudeneutZ2neutZ4=0, H0amplitudeneutZ3neutZ3=0, H0amplitudeneutZ3neutZ4=0, H0amplitudeneutZ4neutZ4=0, H0amplitudecharW1charW1=0, H0amplitudecharW1charW2=0, H0amplitudecharW2charW2=0, H0amplitudeh0h0=0, H0amplitudehiggsAhiggsA=0, H0amplitudeHplusHminus=0, H0amplitudehiggsAZboson=0, H0amplitudesupLantisupL=0, H0amplitudesupLantisupR=0, H0amplitudesupRantisupL=0, H0amplitudesupRantisupR=0, H0amplitudesdownLantisdownL=0, H0amplitudesdownLantisdownR=0, H0amplitudesdownRantisdownL=0, H0amplitudesdownRantisdownR=0, H0amplitudescharmLantischarmL=0, H0amplitudescharmLantischarmR=0, H0amplitudescharmRantischarmL=0, H0amplitudescharmRantischarmR=0, H0amplitudesstrangeLantisstrangeL=0, H0amplitudesstrangeLantisstrangeR=0, H0amplitudesstrangeRantisstrangeL=0, H0amplitudesstrangeRantisstrangeR=0, H0amplitudesnueLantisnueL=0, H0amplitudeselectronLantiselectronL=0, H0amplitudeselectronRantiselectronR=0, H0amplitudeselectronLantiselectronR=0, H0amplitudeselectronRantiselectronL=0, H0amplitudesnumuLantisnumuL=0, H0amplitudesnutauLantisnutauL=0, H0amplitudesmuonLantismuonL=0, H0amplitudesmuonRantismuonR=0, H0amplitudesmuonLantismuonR=0, H0amplitudesmuonRantismuonL=0, H0amplitudestau1antistau1=0, H0amplitudestau2antistau2=0, H0amplitudestau1antistau2=0, H0amplitudestau2antistau1=0, H0amplitudestop1antistop1=0, H0amplitudestop1antistop2=0, H0amplitudestop2antistop1=0, H0amplitudestop2antistop2=0, H0amplitudesbottom1antisbottom1=0, H0amplitudesbottom1antisbottom2=0, H0amplitudesbottom2antisbottom1=0, H0amplitudesbottom2antisbottom2=0, H0amplitudegluongluon=0, H0amplitudegammagamma=0, H0amplitudeWbosonWboson=0, H0amplitudeZbosonZboson=0, H0amplitudeZgamma=0, H0amplitudehiggsAhiggsA2=0, H0amplitudehiggsA2higgsA2=0, H0amplitudehiggsA2Zboson=0, H0amplitudeWHpm=0;

 double H0amplitudeneutZ1neutZ5 = 0, H0amplitudeneutZ2neutZ5 = 0, H0amplitudeneutZ3neutZ5 = 0, H0amplitudeneutZ4neutZ5 = 0, H0amplitudeneutZ5neutZ5 = 0;

 if (flagH2 == 1) {
   if (QCDcorr == false) {
     ///No decays to u or d as negligible as PW proportional to yukawas squared
     ///Use running masses here to try to approximate some of the correction (which aren't included)
     H0amplitudecantic = higgslorHamplitudedecayquarkantiquark (mh0(2), runmc, g, alpha, beta, runmw, 1, 'H', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH);
     H0amplitudesantis = higgslorHamplitudedecayquarkantiquark (mh0(2), runms, g, alpha, beta, runmw, 0, 'H', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH);
     H0amplitudebantib = higgslorHamplitudedecayquarkantiquark (mh0(2), runmb, g, alpha, beta, runmw, 0, 'H', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH); ///use runmb here to reduce necessary corrections
     H0amplitudetantit = higgslorHamplitudedecayquarkantiquark (mh0(2), runmt, g, alpha, beta, runmw, 1, 'H', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH); ///may rather use mtPole here as closer to answer with corrections than runmt gives
   }
   else {
     H0amplitudecantic = higgslorHamplitudedecayquarkantiquark (mh0(2), mcpole, g, alpha, beta, runmw, 1, 'H', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH);
     H0amplitudesantis = higgslorHamplitudedecayquarkantiquark (mh0(2), mspole, g, alpha, beta, runmw, 0, 'H', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH);
     ///mcpole and mspole set in decays.h, this values used are those appropriate for the scheme used for the h -> qq QCD corrections, as in hdecay
     H0amplitudebantib = higgslorHamplitudedecayquarkantiquark (mh0(2), mbPole, g, alpha, beta, runmw, 0, 'H', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH); 
     H0amplitudetantit = higgslorHamplitudedecayquarkantiquark (mh0(2), mtPole, g, alpha, beta, runmw, 1, 'H', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH);
   }

   H0amplitudeeantie = higgslorHamplitudedecayquarkantiquark (mh0(2), runmel, g, alpha, beta, runmw, 0, 'H', CPEMix, nmssmIsIt, false, alphasAtMh)/3; ///0 as leptons are like down-type quarks, divide by 3 as No of colours is 1 for leptons cf 3 for quarks
   H0amplitudemuantimu = higgslorHamplitudedecayquarkantiquark (mh0(2), runmmu, g, alpha, beta, runmw, 0, 'H', CPEMix, nmssmIsIt, false, alphasAtMh)/3;
   H0amplitudetauantitau = higgslorHamplitudedecayquarkantiquark (mh0(2), runmtau, g, alpha, beta, runmw, 0, 'H', CPEMix, nmssmIsIt, false, alphasAtMh)/3;
   
   if (nmssmIsIt == false) {
     H0amplitudesupLantisupL = higgsHamplitudedecay2squarksamehand (mh0(2), mu(1,1), mu(1,1), g, gp, alpha, beta, runmw, runmu, runmd, 1);
     H0amplitudesupRantisupR = higgsHamplitudedecay2squarksamehand (mh0(2), mu(2,1), mu(2,1), g, gp, alpha, beta, runmw, runmu, runmd, 3);
     H0amplitudesdownLantisdownL = higgsHamplitudedecay2squarksamehand (mh0(2), md(1,1), md(1,1), g, gp, alpha, beta, runmw, runmu, runmd, 2);
     H0amplitudesdownRantisdownR = higgsHamplitudedecay2squarksamehand (mh0(2), md(2,1), mu(2,1), g, gp, alpha, beta, runmw, runmu, runmd, 4);
     H0amplitudescharmLantischarmL = higgsHamplitudedecay2squarksamehand (mh0(2), mu(1,2), mu(1,2), g, gp, alpha, beta, runmw, runmc, runms, 1);
     H0amplitudescharmRantischarmR = higgsHamplitudedecay2squarksamehand (mh0(2), mu(2,2), mu(2,2), g, gp, alpha, beta, runmw, runmc, runms, 3);
     H0amplitudesstrangeLantisstrangeL = higgsHamplitudedecay2squarksamehand (mh0(2), md(2,1), md(2,1), g, gp, alpha, beta, runmw, runmc, runms, 2);
     H0amplitudesstrangeRantisstrangeR = higgsHamplitudedecay2squarksamehand (mh0(2), md(2,2), md(2,2), g, gp, alpha, beta, runmw, runmc, runms, 4);
     H0amplitudesnueLantisnueL = higgsHamplitudedecay2sleptonsamehand (mh0(2), msnu(1), msnu(1), g, gp, alpha, beta, runmw, runmel, 1);
     H0amplitudeselectronLantiselectronL = higgsHamplitudedecay2sleptonsamehand (mh0(2), me(1,1), me(1,1), g, gp, alpha, beta, runmw, runmel, 2);
     H0amplitudeselectronRantiselectronR = higgsHamplitudedecay2sleptonsamehand (mh0(2), me(2,1), me(2,1), g, gp, alpha, beta, runmw, runmel, 3);
     H0amplitudesnumuLantisnumuL = higgsHamplitudedecay2sleptonsamehand (mh0(2), msnu(2), msnu(2), g, gp, alpha, beta, runmw, runmmu, 1);
     H0amplitudesmuonLantismuonL = higgsHamplitudedecay2sleptonsamehand (mh0(2), me(1,2), me(1,2), g, gp, alpha, beta, runmw, runmmu, 2);
     H0amplitudesmuonRantismuonR = higgsHamplitudedecay2sleptonsamehand (mh0(2), me(2,2), me(2,2), g, gp, alpha, beta, runmw, runmmu, 3);
     H0amplitudesnutauLantisnutauL = higgsHamplitudedecay2sleptonsamehand (mh0(2), msnu(3), msnu(3), g, gp, alpha, beta, runmw, runmtau, 1);
     H0amplitudesupLantisupR = higgsHamplitudedecay2squarkdiffhand (mh0(2), mu(1,1), mu(2,1), g, alpha, beta, runmw, runmu, runmd, greekmu, Au, Ad, 1);
     H0amplitudesupRantisupL = higgsHamplitudedecay2squarkdiffhand (mh0(2), mu(1,1), mu(2,1), g, alpha, beta, runmw, runmu, runmd, greekmu, Au, Ad, 1);
     H0amplitudesdownLantisdownR = higgsHamplitudedecay2squarkdiffhand (mh0(2), md(1,1), md(2,1), g, alpha, beta, runmw, runmu, runmd, greekmu, Au, Ad, 2);
     H0amplitudesdownRantisdownL = higgsHamplitudedecay2squarkdiffhand (mh0(2), md(1,1), md(2,1), g, alpha, beta, runmw, runmu, runmd, greekmu, Au, Ad, 2);
     H0amplitudescharmLantischarmR = higgsHamplitudedecay2squarkdiffhand (mh0(2), mu(1,2), mu(2,2), g, alpha, beta, runmw, runmc, runms, greekmu, Ac, As, 1);
     H0amplitudescharmRantischarmL = higgsHamplitudedecay2squarkdiffhand (mh0(2), mu(1,2), mu(2,2), g, alpha, beta, runmw, runmc, runms, greekmu, Ac, As, 1);
     H0amplitudesstrangeLantisstrangeR = higgsHamplitudedecay2squarkdiffhand (mh0(2), md(2,1), md(2,2), g, alpha, beta, runmw, runmc, runms, greekmu, Ac, As, 2);
     H0amplitudesstrangeRantisstrangeL = higgsHamplitudedecay2squarkdiffhand (mh0(2), md(2,1), md(2,2), g, alpha, beta, runmw, runmc, runms, greekmu, Ac, As, 2);
     H0amplitudeselectronLantiselectronR = higgsHamplitudedecay2sleptondiffhand (mh0(2), me(1,1), me(2,1), g, alpha, beta, runmw, runmel, greekmu, Ae, 1);
     H0amplitudeselectronRantiselectronL = higgsHamplitudedecay2sleptondiffhand (mh0(2), me(2,1), me(1,1), g, alpha, beta, runmw, runmel, greekmu, Ae, 1);
     H0amplitudesmuonLantismuonR = higgsHamplitudedecay2sleptondiffhand (mh0(2), me(1,2), me(2,2), g, alpha, beta, runmw, runmmu, greekmu, Amu, 1);
     H0amplitudesmuonRantismuonL = higgsHamplitudedecay2sleptondiffhand (mh0(2), me(2,2), me(2,1), g, alpha, beta, runmw, runmmu, greekmu, Amu, 1);
     H0amplitudestop1antistop1 = higgsHamplitudedecaystop1stop1 (mh0(2), mu(1,3), mu(1,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetat); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     H0amplitudestop2antistop2 = higgsHamplitudedecaystop2stop2 (mh0(2), mu(2,3), mu(2,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetat); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     H0amplitudestop1antistop2 = higgsHamplitudedecaystop1stop2 (mh0(2), mu(1,3), mu(2,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetat); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     H0amplitudestop2antistop1 = higgsHamplitudedecaystop1stop2 (mh0(2), mu(2,3), mu(1,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetat); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     H0amplitudesbottom1antisbottom1 = higgsHamplitudedecaysbottom1sbottom1(mh0(2), md(1,3), md(1,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetab); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     H0amplitudesbottom2antisbottom2 = higgsHamplitudedecaysbottom2sbottom2(mh0(2), md(2,3), md(2,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetab); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     H0amplitudesbottom1antisbottom2 = higgsHamplitudedecaysbottom1sbottom2 (mh0(2), md(1,3), md(2,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetab); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     H0amplitudesbottom2antisbottom1 = higgsHamplitudedecaysbottom1sbottom2 (mh0(2), md(2,3), md(1,3), g, gp, alpha, beta, runmw, runmt, runmb, greekmu, At, Ab, thetab); ///use runmt and runmb here as mass used to set yukawa coupling, note however pole masses give greater agreement with susyhit as susyhit uses non-running masses here, similar for other higgs decays to third generation sfermions
     H0amplitudestau1antistau1 = higgsHamplitudedecaystau1stau1 (mh0(2), me(1,3), me(1,3), g, gp, alpha, beta, runmw, runmtau, greekmu, Atau, thetatau);
     H0amplitudestau2antistau2 = higgsHamplitudedecaystau2stau2 (mh0(2), me(2,3), me(2,3), g, gp, alpha, beta, runmw, runmtau, greekmu, Atau, thetatau);
     H0amplitudestau1antistau2 = higgsHamplitudedecaystau1stau2 (mh0(2), me(1,3), me(2,3), g, gp, alpha, beta, runmw, runmtau, greekmu, Atau, thetatau);
     H0amplitudestau2antistau1 = higgsHamplitudedecaystau1stau2 (mh0(2), me(2,3), me(1,3), g, gp, alpha, beta, runmw, runmtau, greekmu, Atau, thetatau);
     
     H0amplitudecharW1charW1 = higgsphiamplitudedecaysamechargino (mh0(2), mch(1), g, thetaL2, thetaR2, alpha, beta, 1, 'H');
     H0amplitudecharW2charW2 = higgsphiamplitudedecaysamechargino (mh0(2), mch(2), g, thetaL2, thetaR2, alpha, beta, 2, 'H');
     H0amplitudecharW1charW2 = higgsphiamplitudedecaydifchargino (mh0(2), mch(1), mch(2), g, thetaL2, thetaR2, alpha, beta, 'H');
     
     H0amplitudegammagamma = higgsesamplitudedecaygammagammatotal(mh0(2), g, gp, alphaAtMH, runmw, polemw, alpha, beta, mtAtMH, mbAtMH, mcAtMH, runmtau, mHpm, mu(1,3), mu(2,3), md(1,3), md(2,3), me(1,3), me(2,3), mch(1), mch(2), thetaL, thetaR, thetat, thetab, thetatau, greekmu, At, Ab, Atau, 'H'); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     H0amplitudegluongluon = higgsesamplitudedecaygluongluontotal(mh0(2), g, g3atmH0, gp, runmw, alpha, beta, mtPole, mbPole, mcpole, mu(1,3), mu(2,3), md(1,3), md(2,3), thetat, thetab, greekmu, At, Ab, mspole, mu(1,2), mu(2,2), md(1,2), md(2,2), Ac, As, runmu, runmd, mu(1,1), mu(2,1), md(1,1), md(2,1), Au, Ad, 'H', QCDcorr); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     H0amplitudeZgamma = higgsesamplitudedecayZbosonphotontotal(mh0(2), polemz, g, gp, alphaAtMH, polemw, runmw, alpha, beta, mtAtMH, mbAtMH, mcAtMH, msAtMH, mu(1,3), mu(2,3), md(1,3), md(2,3), mHpm, thetat, thetab, greekmu, At, Ab, 'H'); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW

     H0amplitudeneutZ1neutZ1 = higgsphiamplitudedecayneutralinoneutralino (mh0(2), mneut(1), mneut(1), g, tanthetaW, alpha, mixNeut, 1, 1, 'H');
     H0amplitudeneutZ1neutZ2 = higgsphiamplitudedecayneutralinoneutralino (mh0(2), mneut(1), mneut(2), g, tanthetaW, alpha, mixNeut, 1, 2, 'H');
     H0amplitudeneutZ1neutZ3 = higgsphiamplitudedecayneutralinoneutralino (mh0(2), mneut(1), mneut(3), g, tanthetaW, alpha, mixNeut, 1, 3, 'H');
     H0amplitudeneutZ1neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mh0(2), mneut(1), mneut(4), g, tanthetaW, alpha, mixNeut, 1, 4, 'H');
     H0amplitudeneutZ2neutZ2 = higgsphiamplitudedecayneutralinoneutralino (mh0(2), mneut(2), mneut(2), g, tanthetaW, alpha, mixNeut, 2, 2, 'H');
     H0amplitudeneutZ2neutZ3 = higgsphiamplitudedecayneutralinoneutralino (mh0(2), mneut(2), mneut(3), g, tanthetaW, alpha, mixNeut, 2, 3, 'H');
     H0amplitudeneutZ2neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mh0(2), mneut(2), mneut(4), g, tanthetaW, alpha, mixNeut, 2, 4, 'H');
     H0amplitudeneutZ3neutZ3 = higgsphiamplitudedecayneutralinoneutralino (mh0(2), mneut(3), mneut(3), g, tanthetaW, alpha, mixNeut, 3, 3, 'H');
     H0amplitudeneutZ3neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mh0(2), mneut(3), mneut(4), g, tanthetaW, alpha, mixNeut, 3, 4, 'H');
     H0amplitudeneutZ4neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mh0(2), mneut(4), mneut(4), g, tanthetaW, alpha, mixNeut, 4, 4, 'H');
     
     H0amplitudehiggsAhiggsA = higgsHamplitudedecayAA (mh0(2), mA0(1), g, gp, alpha, beta, runmw);
     H0amplitudehiggsAZboson = higgsHamplitudedecayhiggsAZboson (mh0(2), polemz, mA0(1), g, gp, alpha, beta);
     H0amplitudeHplusHminus = higgsHamplitudedecayHplusHminus (mh0(2), mHpm, g, gp, alpha, beta, runmw);
     H0amplitudeh0h0 = higgsHamplitudedecayhh (mh0(2), mh0(1), g, gp, alpha, beta, runmw);
   }
   else if (nmssmIsIt == true) {
     H0amplitudesupLantisupL = higgsHamplitudedecay2squarksamehandNMSSM (mh0(2), mu(1,1), mu(1,1), g, gp, alpha, beta, runmw, runmu, CPEMix, 1);
     H0amplitudesupRantisupR = higgsHamplitudedecay2squarksamehandNMSSM (mh0(2), mu(2,1), mu(2,1), g, gp, alpha, beta, runmw, runmu, CPEMix, 3);
     H0amplitudesdownLantisdownL = higgsHamplitudedecay2squarksamehandNMSSM (mh0(2), md(1,1), md(1,1), g, gp, alpha, beta, runmw, runmd, CPEMix, 2);
     H0amplitudesdownRantisdownR = higgsHamplitudedecay2squarksamehandNMSSM (mh0(2), md(2,1), mu(2,1), g, gp, alpha, beta, runmw, runmd, CPEMix, 4);
     H0amplitudescharmLantischarmL = higgsHamplitudedecay2squarksamehandNMSSM (mh0(2), mu(1,2), mu(1,2), g, gp, alpha, beta, runmw, runmc, CPEMix, 1);
     H0amplitudescharmRantischarmR = higgsHamplitudedecay2squarksamehandNMSSM (mh0(2), mu(2,2), mu(2,2), g, gp, alpha, beta, runmw, runmc, CPEMix, 3);
     H0amplitudesstrangeLantisstrangeL = higgsHamplitudedecay2squarksamehandNMSSM (mh0(2), md(2,1), md(2,1), g, gp, alpha, beta, runmw, runms, CPEMix, 2);
     H0amplitudesstrangeRantisstrangeR = higgsHamplitudedecay2squarksamehandNMSSM (mh0(2), md(2,2), md(2,2), g, gp, alpha, beta, runmw, runms, CPEMix, 4);
     H0amplitudesnueLantisnueL = higgsHamplitudedecay2sleptonsamehandNMSSM (mh0(2), msnu(1), msnu(1), g, gp, alpha, beta, runmw, 0, CPEMix, 1);
     H0amplitudeselectronLantiselectronL = higgsHamplitudedecay2sleptonsamehandNMSSM (mh0(2), me(1,1), me(1,1), g, gp, alpha, beta, runmw, runmel, CPEMix, 2);
     H0amplitudeselectronRantiselectronR = higgsHamplitudedecay2sleptonsamehandNMSSM (mh0(2), me(2,1), me(2,1), g, gp, alpha, beta, runmw, runmel, CPEMix, 3);
     H0amplitudesnumuLantisnumuL = higgsHamplitudedecay2sleptonsamehandNMSSM (mh0(2), msnu(2), msnu(2), g, gp, alpha, beta, runmw, 0, CPEMix, 1);
     H0amplitudesmuonLantismuonL = higgsHamplitudedecay2sleptonsamehandNMSSM (mh0(2), me(1,2), me(1,2), g, gp, alpha, beta, runmw, runmmu, CPEMix, 2);
     H0amplitudesmuonRantismuonR = higgsHamplitudedecay2sleptonsamehandNMSSM (mh0(2), me(2,2), me(2,2), g, gp, alpha, beta, runmw, runmmu, CPEMix, 3);
     H0amplitudesnutauLantisnutauL = higgsHamplitudedecay2sleptonsamehandNMSSM (mh0(2), msnu(3), msnu(3), g, gp, alpha, beta, runmw, 0, CPEMix, 1);
     
     H0amplitudesupLantisupR = higgshamplitudedecay2squarkdiffhandNMSSM(mh0(2), mu(1,1), mu(2,1), g, alpha, beta, runmw, runmu, Au, mueff, lam, CPEMix, 1, 2);
     H0amplitudesupRantisupL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), mu(1,1), mu(2,1), g, alpha, beta, runmw, runmu, Au, mueff, lam, CPEMix, 1, 2);
     H0amplitudesdownLantisdownR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), md(1,1), md(2,1), g, alpha, beta, runmw, runmd, Ad, mueff, lam, CPEMix, 2, 2);
     H0amplitudesdownRantisdownL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), md(1,1), md(2,1), g, alpha, beta, runmw, runmd, Ad, mueff, lam, CPEMix, 2, 2);
     H0amplitudescharmLantischarmR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), mu(1,2), mu(2,2), g, alpha, beta, runmw, runmc, Ac, mueff, lam, CPEMix, 1, 2);
     H0amplitudescharmRantischarmL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), mu(1,2), mu(2,2), g, alpha, beta, runmw, runmc, Ac, mueff, lam, CPEMix, 1, 2);
     H0amplitudesstrangeLantisstrangeR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), md(2,1), md(2,2), g, alpha, beta, runmw, runms, As, mueff, lam, CPEMix, 2, 2);
     H0amplitudesstrangeRantisstrangeL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), md(2,1), md(2,2), g, alpha, beta, runmw, runms, As, mueff, lam, CPEMix, 2, 2);
     H0amplitudeselectronLantiselectronR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), me(1,1), me(2,1), g, alpha, beta, runmw, runmel, Ae, mueff, lam, CPEMix, 2, 2)/3;
     H0amplitudeselectronRantiselectronL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), me(2,1), me(1,1), g, alpha, beta, runmw, runmel, Ae, mueff, lam, CPEMix, 2, 2)/3;
     H0amplitudesmuonLantismuonR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), me(1,2), me(2,2), g, alpha, beta, runmw, runmmu, Amu, mueff, lam, CPEMix, 2, 2)/3;
     H0amplitudesmuonRantismuonL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(2), me(2,2), me(2,1), g, alpha, beta, runmw, runmmu, Amu, mueff, lam, CPEMix, 2, 2)/3;
     
     H0amplitudestop1antistop1 = higgsCPevenamplitudedecaystopistopiNMSSM (mh0(2), mu(1,3), thetat, runmt, g, gp, runmw, beta, CPEMix, At, mueff, lam, 1, 2);
     H0amplitudestop2antistop2 = higgsCPevenamplitudedecaystopistopiNMSSM (mh0(2), mu(2,3), thetat, runmt, g, gp, runmw, beta, CPEMix, At, mueff, lam, 2, 2);
     H0amplitudestop1antistop2 = higgsCPevenamplitudedecaystopistopjNMSSM (mh0(2), mu(1,3), mu(2,3), thetat, runmt, g, gp, runmw, beta, CPEMix, At, mueff, lam, 2); 
     H0amplitudestop2antistop1 = h0amplitudestop1antistop2;
     H0amplitudesbottom1antisbottom1 = higgsCPevenamplitudedecaysbottomisbottomiNMSSM (mh0(2), md(1,3), thetab, runmb, g, gp, runmw, beta,  CPEMix, Ab, mueff, lam, 1, 2);
     H0amplitudesbottom2antisbottom2 = higgsCPevenamplitudedecaysbottomisbottomiNMSSM (mh0(2), md(2,3), thetab, runmb, g, gp, runmw, beta,  CPEMix, Ab, mueff, lam, 2, 2);
     H0amplitudesbottom1antisbottom2 = higgsCPevenamplitudedecaysbottomisbottomjNMSSM (mh0(2), md(1,3), md(2,3), thetab, runmb, g, gp, runmw, beta, CPEMix, Ab, mueff, lam, 2);
     H0amplitudesbottom2antisbottom1 = h0amplitudesbottom1antisbottom2;
     H0amplitudestau1antistau1 = higgsCPevenamplitudedecaystauistauiNMSSM (mh0(2), me(1,3), thetatau - PI/2, runmtau, g, gp, runmw, beta,  CPEMix, Atau, mueff, lam, 1, 2); 
     H0amplitudestau2antistau2 = higgsCPevenamplitudedecaystauistauiNMSSM (mh0(2), me(2,3), thetatau - PI/2, runmtau, g, gp, runmw, beta,  CPEMix, Atau, mueff, lam, 2, 2);
     H0amplitudestau1antistau2 = higgsCPevenamplitudedecaystauistaujNMSSM (mh0(2), me(1,3), me(2,3), thetatau - PI/2, runmtau, g, gp, runmw, beta, CPEMix, Atau, mueff, lam, 2);
     H0amplitudestau2antistau1 = h0amplitudestau1antistau2;
     
     H0amplitudecharW1charW1 = higgsphiamplitudedecaysamecharginoNMSSM (mh0(2), mch(1), g, thetaL2, thetaR2, lam, CPEMix, 1, 2);
     H0amplitudecharW2charW2 = higgsphiamplitudedecaysamecharginoNMSSM (mh0(2), mch(2), g, thetaL2, thetaR2, lam, CPEMix, 2, 2);
     H0amplitudecharW1charW2 = higgsphiamplitudedecaydiffcharginoNMSSM (mh0(2), mch(1), mch(2), g, thetaL2, thetaR2, lam, CPEMix, 2);
     H0amplitudegammagamma = higgsCPevenamplitudedecaygammagammaNMSSM(mh0(2), mtAtMH, mbAtMH, mcAtMH, runmtau, runmw, mHpm, mch(1), mch(2), mu(1,2), mu(2,2), mu(1,3), mu(2,3), md(1,2), md(2,2), md(1,3), md(2,3), me(1,2), me(2,2), me(1,3), me(2,3), CPEMix, beta, g, gp, alphaAtMH, thetat, thetab, thetatau-PI/2, thetaL2, thetaR2, At, Ab, Atau, greekmu, mueff, lam, kappa, Alambda, 2); //Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     H0amplitudegluongluon = higgsCPevenamplitudedecaygluongluonNMSSM(mh0(2), mtPole, mbPole, mcpole, runmw, mu(1,2), mu(2,2), mu(1,3), mu(2,3), md(1,2), md(2,2), md(1,3), md(2,3), mu(1,1), mu(2,1), md(1,1), md(2,1), mtAtMH, mbAtMH, CPEMix, beta, g, gp, gs, alphasAtMH, thetat, thetab, thetaL2, thetaR2, At, Ab, greekmu, mueff, lam, kappa, Alambda, 2, QCDcorr);///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     H0amplitudeZgamma = higgshamplitudedecayZgammaNMSSM (mh0(2), g, gp, alphaAtMH, runmw, polemz, mHpm, CPEMix, beta, mtAtMH, mbAtMH, mcAtMH, mch(1), mch(2), thetaL2, thetaR2, lam, kappa, Alambda, greekmu, mueff, 2);///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     
     H0amplitudeneutZ1neutZ1 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(1), mneut(1), g, gp, CPEMix, mixNeut, lam, kappa, 1, 1, 2);
     H0amplitudeneutZ1neutZ2 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(1), mneut(2), g, gp, CPEMix, mixNeut, lam, kappa, 1, 2, 2);
     H0amplitudeneutZ1neutZ3 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(1), mneut(3), g, gp, CPEMix, mixNeut, lam, kappa, 1, 3, 2);
     H0amplitudeneutZ1neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(1), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 1, 4, 2);
     H0amplitudeneutZ1neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(1), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 1, 5, 2);
     H0amplitudeneutZ2neutZ2 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(2), mneut(2), g, gp, CPEMix, mixNeut, lam, kappa, 2, 2, 2);
     H0amplitudeneutZ2neutZ3 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(2), mneut(3), g, gp, CPEMix, mixNeut, lam, kappa, 2, 3, 2);
     H0amplitudeneutZ2neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(2), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 2, 4, 2);
     H0amplitudeneutZ2neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(2), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 2, 5, 2);
     H0amplitudeneutZ3neutZ3 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(3), mneut(3), g, gp, CPEMix, mixNeut, lam, kappa, 3, 3, 2);
     H0amplitudeneutZ3neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(3), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 3, 4, 2);
     H0amplitudeneutZ3neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(3), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 3, 5, 2);
     H0amplitudeneutZ4neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(4), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 4, 4, 2);
     H0amplitudeneutZ4neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(4), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 4, 5, 2);
     H0amplitudeneutZ5neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(2), mneut(5), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 5, 5, 2);
     H0amplitudehiggsAhiggsA = higgsCPevenamplitudedecayAANMSSM(mh0(2), mA0(1), mA0(1), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 2, 1, 1);
     H0amplitudehiggsAhiggsA2 = 2*higgsCPevenamplitudedecayAANMSSM(mh0(2), mA0(1), mA0(2), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 2, 1, 2);
     H0amplitudehiggsA2higgsA2 = higgsCPevenamplitudedecayAANMSSM(mh0(2), mA0(2), mA0(2), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 2, 2, 2);
     H0amplitudehiggsAZboson = higgsCPevenamplitudedecaypseudoscalarZNMSSM(mh0(2), mA0(1), polemz, g, gp, beta, CPEMix, CPOMix, 2, 1);
     H0amplitudehiggsA2Zboson = higgsCPevenamplitudedecaypseudoscalarZNMSSM(mh0(2), mA0(2), polemz, g, gp, beta, CPEMix, CPOMix, 2, 2);
     H0amplitudeHplusHminus = higgsCPevenamplitudedecayHpHmNMSSM (mh0(2), mHpm, runmw, g, gp, runmt, runmb, beta, lam, mueff, kappa, Alambda, CPEMix, 2);
     H0amplitudeh0h0 = higgsCPevenamplitudedecayhhorhHorHHNMSSM(mh0(2), mh0(1), mh0(1), g, gp, runmw, runmt, runmb, beta, lam,  Alambda, kappa, Akappa, mueff, CPEMix, CPOMix, 1, 1, 2);
     H0amplitudeWHpm = higgsCPevenamplitudedecayWHpmNMSSM (mh0(2), polemw, mHpm, beta, g, CPEMix, 2)*2; ///*2 as W+H- or W-H+
   }
  
   H0amplitudeWbosonWboson = higgsHamplitudedecayVV(mh0(2), polemw, polemz, g, gp, alpha, beta, 'W', CPEMix, nmssmIsIt)(1);
   H0amplitudeZbosonZboson = higgsHamplitudedecayVV(mh0(2), polemw, polemz, g, gp, alpha, beta, 'Z', CPEMix, nmssmIsIt)(1);
   
   int H0WWcommentcode, H0ZZcommentcode, H0WWNDA=0, H0ZZNDA=0;
   H0WWcommentcode = higgsHamplitudedecayVV(mh0(2), polemw, polemz, g, gp, alpha, beta, 'W', CPEMix, nmssmIsIt)(2);
   string H0WWcomment, H0ZZcomment;
   if (H0WWcommentcode == 1) {
     H0WWcomment = "# H -> WW* -> W f f'bar";
     H0WWNDA = 2; ///So read into other programs, e.g. PYTHIA, correctly
   }
   else if(H0WWcommentcode == 2) {
     H0WWcomment = "# H -> W+ W-";
     H0WWNDA = 2;
   }
   
   H0ZZcommentcode = higgsHamplitudedecayVV(mh0(2), polemw, polemz, g, gp, alpha, beta, 'Z', CPEMix, nmssmIsIt)(2);
   if (H0ZZcommentcode == 1) {
     H0ZZcomment = "# H -> ZZ* -> Z f f'bar";
     H0ZZNDA = 2; ///So read into other programs, e.g. PYTHIA, correctly
   }
   else if(H0ZZcommentcode == 2) {
     H0ZZcomment = "# H -> Z Z";
     H0ZZNDA = 2;
   }
   
   ParticleHiggsH.Array_Decays[0][0] = PDGup; ParticleHiggsH.Array_Decays[0][1] = -PDGup; ParticleHiggsH.Array_Decays[0][2] = H0amplitudeuantiu; ParticleHiggsH.Array_Decays[0][3] = 2; ParticleHiggsH.Array_Comments[0] = "# H -> u ub";
   ParticleHiggsH.Array_Decays[1][0] = PDGdown; ParticleHiggsH.Array_Decays[1][1] = -PDGdown; ParticleHiggsH.Array_Decays[1][2] = H0amplitudedantid; ParticleHiggsH.Array_Decays[1][3] = 2; ParticleHiggsH.Array_Comments[1] = "# H -> d db";
   ParticleHiggsH.Array_Decays[2][0] = PDGcharm; ParticleHiggsH.Array_Decays[2][1] = -PDGcharm; ParticleHiggsH.Array_Decays[2][2] = H0amplitudecantic; ParticleHiggsH.Array_Decays[2][3] = 2; ParticleHiggsH.Array_Comments[2] = "# H -> c cb";
   ParticleHiggsH.Array_Decays[3][0] = PDGstrange; ParticleHiggsH.Array_Decays[3][1] = -PDGstrange; ParticleHiggsH.Array_Decays[3][2] = H0amplitudesantis; ParticleHiggsH.Array_Decays[3][3] = 2; ParticleHiggsH.Array_Comments[3] = "# H -> s sb";
   ParticleHiggsH.Array_Decays[4][0] = PDGbottom; ParticleHiggsH.Array_Decays[4][1] = -PDGbottom; ParticleHiggsH.Array_Decays[4][2] = H0amplitudebantib; ParticleHiggsH.Array_Decays[4][3] = 2; ParticleHiggsH.Array_Comments[4] = "# H -> b bb";
   ParticleHiggsH.Array_Decays[5][0] = PDGtop; ParticleHiggsH.Array_Decays[5][1] = -PDGtop; ParticleHiggsH.Array_Decays[5][2] = H0amplitudetantit; ParticleHiggsH.Array_Decays[5][3] = 2; ParticleHiggsH.Array_Comments[5] = "# H -> t tb";
   ParticleHiggsH.Array_Decays[6][0] = PDGelectron; ParticleHiggsH.Array_Decays[6][1] = -PDGelectron; ParticleHiggsH.Array_Decays[6][2] = H0amplitudeeantie; ParticleHiggsH.Array_Decays[6][3] = 2; ParticleHiggsH.Array_Comments[6] = "# H -> e- e+";
   ParticleHiggsH.Array_Decays[7][0] = PDGmuon; ParticleHiggsH.Array_Decays[7][1] = -PDGmuon; ParticleHiggsH.Array_Decays[7][2] = H0amplitudemuantimu; ParticleHiggsH.Array_Decays[7][3] = 2; ParticleHiggsH.Array_Comments[7] = "# H -> mu- mu+";
   ParticleHiggsH.Array_Decays[8][0] = PDGtau; ParticleHiggsH.Array_Decays[8][1] = -PDGtau; ParticleHiggsH.Array_Decays[8][2] = H0amplitudetauantitau; ParticleHiggsH.Array_Decays[8][3] = 2; ParticleHiggsH.Array_Comments[8] = "# H -> tau- tau+";
   ParticleHiggsH.Array_Decays[9][0] = PDGneutralino1; ParticleHiggsH.Array_Decays[9][1] = PDGneutralino1; ParticleHiggsH.Array_Decays[9][2] = H0amplitudeneutZ1neutZ1; ParticleHiggsH.Array_Decays[9][3] = 2; ParticleHiggsH.Array_Comments[9] = "# H -> ~chi_10 ~chi_10";
   ParticleHiggsH.Array_Decays[10][0] = PDGneutralino1; ParticleHiggsH.Array_Decays[10][1] = PDGneutralino2; ParticleHiggsH.Array_Decays[10][2] = H0amplitudeneutZ1neutZ2; ParticleHiggsH.Array_Decays[10][3] = 2; ParticleHiggsH.Array_Comments[10] = "# H -> ~chi_10 ~chi_20";
   ParticleHiggsH.Array_Decays[11][0] = PDGneutralino1; ParticleHiggsH.Array_Decays[11][1] = PDGneutralino3; ParticleHiggsH.Array_Decays[11][2] = H0amplitudeneutZ1neutZ3; ParticleHiggsH.Array_Decays[11][3] = 2; ParticleHiggsH.Array_Comments[11] = "# H -> ~chi_10 ~chi_30";
   ParticleHiggsH.Array_Decays[12][0] = PDGneutralino1; ParticleHiggsH.Array_Decays[12][1] = PDGneutralino4; ParticleHiggsH.Array_Decays[12][2] = H0amplitudeneutZ1neutZ4; ParticleHiggsH.Array_Decays[12][3] = 2; ParticleHiggsH.Array_Comments[12] = "# H -> ~chi_10 ~chi_40";
   ParticleHiggsH.Array_Decays[13][0] = PDGneutralino2; ParticleHiggsH.Array_Decays[13][1] = PDGneutralino2; ParticleHiggsH.Array_Decays[13][2] = H0amplitudeneutZ2neutZ2; ParticleHiggsH.Array_Decays[13][3] = 2; ParticleHiggsH.Array_Comments[13] = "# H -> ~chi_20 ~chi_20";
   ParticleHiggsH.Array_Decays[14][0] = PDGneutralino2; ParticleHiggsH.Array_Decays[14][1] = PDGneutralino3; ParticleHiggsH.Array_Decays[14][2] = H0amplitudeneutZ2neutZ3; ParticleHiggsH.Array_Decays[14][3] = 2; ParticleHiggsH.Array_Comments[14] = "# H -> ~chi_20 ~chi_30";
   ParticleHiggsH.Array_Decays[15][0] = PDGneutralino2; ParticleHiggsH.Array_Decays[15][1] = PDGneutralino4; ParticleHiggsH.Array_Decays[15][2] = H0amplitudeneutZ2neutZ4; ParticleHiggsH.Array_Decays[15][3] = 2; ParticleHiggsH.Array_Comments[15] = "# H -> ~chi_20 ~chi_40";
   ParticleHiggsH.Array_Decays[16][0] = PDGneutralino3; ParticleHiggsH.Array_Decays[16][1] = PDGneutralino3; ParticleHiggsH.Array_Decays[16][2] = H0amplitudeneutZ3neutZ3; ParticleHiggsH.Array_Decays[16][3] = 2; ParticleHiggsH.Array_Comments[16] = "# H -> ~chi_30 ~chi_30";
   ParticleHiggsH.Array_Decays[17][0] = PDGneutralino3; ParticleHiggsH.Array_Decays[17][1] = PDGneutralino4; ParticleHiggsH.Array_Decays[17][2] = H0amplitudeneutZ3neutZ4; ParticleHiggsH.Array_Decays[17][3] = 2; ParticleHiggsH.Array_Comments[17] = "# H -> ~chi_30 ~chi_40";
   ParticleHiggsH.Array_Decays[18][0] = PDGneutralino4; ParticleHiggsH.Array_Decays[18][1] = PDGneutralino4; ParticleHiggsH.Array_Decays[18][2] = H0amplitudeneutZ4neutZ4; ParticleHiggsH.Array_Decays[18][3] = 2; ParticleHiggsH.Array_Comments[18] = "# H -> ~chi_40 ~chi_40";
   ParticleHiggsH.Array_Decays[19][0] = PDGchargino1; ParticleHiggsH.Array_Decays[19][1] = -PDGchargino1; ParticleHiggsH.Array_Decays[19][2] = H0amplitudecharW1charW1; ParticleHiggsH.Array_Decays[19][3] = 2; ParticleHiggsH.Array_Comments[19] = "# H -> ~chi_1+ ~chi_1-";
   ParticleHiggsH.Array_Decays[20][0] = PDGchargino2; ParticleHiggsH.Array_Decays[20][1] = -PDGchargino2; ParticleHiggsH.Array_Decays[20][2] = H0amplitudecharW2charW2; ParticleHiggsH.Array_Decays[20][3] = 2; ParticleHiggsH.Array_Comments[20] = "# H -> ~chi_2+ ~chi_2-";
   ParticleHiggsH.Array_Decays[21][0] = PDGchargino1; ParticleHiggsH.Array_Decays[21][1] = -PDGchargino2; ParticleHiggsH.Array_Decays[21][2] = H0amplitudecharW1charW2; ParticleHiggsH.Array_Decays[21][3] = 2; ParticleHiggsH.Array_Comments[21] = "# H -> ~chi_1+ ~chi_2-";
   ParticleHiggsH.Array_Decays[22][0] = PDGchargino2; ParticleHiggsH.Array_Decays[22][1] = -PDGchargino1; ParticleHiggsH.Array_Decays[22][2] = H0amplitudecharW1charW2; ParticleHiggsH.Array_Decays[22][3] = 2; ParticleHiggsH.Array_Comments[22] = "# H -> ~chi_2+ ~chi_1-"; ///amplitude same as decay to ~chi1+ and ~chi2- by CP invariance
   ParticleHiggsH.Array_Decays[23][0] = PDGh0; ParticleHiggsH.Array_Decays[23][1] = PDGh0; ParticleHiggsH.Array_Decays[23][2] = H0amplitudeh0h0; ParticleHiggsH.Array_Decays[23][3] = 2; ParticleHiggsH.Array_Comments[23] = "# H -> h h";
   ParticleHiggsH.Array_Decays[24][0] = PDGA0; ParticleHiggsH.Array_Decays[24][1] = PDGA0; ParticleHiggsH.Array_Decays[24][2] = H0amplitudehiggsAhiggsA; ParticleHiggsH.Array_Decays[24][3] = 2; ParticleHiggsH.Array_Comments[24] = "# H -> A A";
   ParticleHiggsH.Array_Decays[25][0] = PDGHplus; ParticleHiggsH.Array_Decays[25][1] = -PDGHplus; ParticleHiggsH.Array_Decays[25][2] = H0amplitudeHplusHminus; ParticleHiggsH.Array_Decays[25][3] = 2; ParticleHiggsH.Array_Comments[25] = "# H -> H+ H-";
   ParticleHiggsH.Array_Decays[26][0] = PDGA0; ParticleHiggsH.Array_Decays[26][1] = PDGZboson; ParticleHiggsH.Array_Decays[26][2] = H0amplitudehiggsAZboson; ParticleHiggsH.Array_Decays[26][3] = 2; ParticleHiggsH.Array_Comments[26] = "# H -> A Z";
   ParticleHiggsH.Array_Decays[27][0] = PDGsupL; ParticleHiggsH.Array_Decays[27][1] = -PDGsupL; ParticleHiggsH.Array_Decays[27][2] = H0amplitudesupLantisupL; ParticleHiggsH.Array_Decays[27][3] = 2; ParticleHiggsH.Array_Comments[27] = "# H -> ~u_L ~u_L*";
   ParticleHiggsH.Array_Decays[28][0] = PDGsupR; ParticleHiggsH.Array_Decays[28][1] = -PDGsupR; ParticleHiggsH.Array_Decays[28][2] = H0amplitudesupRantisupR; ParticleHiggsH.Array_Decays[28][3] = 2; ParticleHiggsH.Array_Comments[28] = "# H -> ~u_R ~u_R*";
   ParticleHiggsH.Array_Decays[29][0] = PDGsupL; ParticleHiggsH.Array_Decays[29][1] = -PDGsupR; ParticleHiggsH.Array_Decays[29][2] = H0amplitudesupLantisupR; ParticleHiggsH.Array_Decays[29][3] = 2; ParticleHiggsH.Array_Comments[29] = "# H -> ~u_L ~u_R*";
   ParticleHiggsH.Array_Decays[30][0] = PDGsupR; ParticleHiggsH.Array_Decays[30][1] = -PDGsupL; ParticleHiggsH.Array_Decays[30][2] = H0amplitudesupRantisupL; ParticleHiggsH.Array_Decays[30][3] = 2; ParticleHiggsH.Array_Comments[30] = "# H -> ~u_R ~u_L*";
   ParticleHiggsH.Array_Decays[31][0] = PDGsdownL; ParticleHiggsH.Array_Decays[31][1] = -PDGsdownL; ParticleHiggsH.Array_Decays[31][2] = H0amplitudesdownLantisdownL; ParticleHiggsH.Array_Decays[31][3] = 2; ParticleHiggsH.Array_Comments[31] = "# H -> ~d_L ~d_L*";
   ParticleHiggsH.Array_Decays[32][0] = PDGsdownR; ParticleHiggsH.Array_Decays[32][1] = -PDGsdownR; ParticleHiggsH.Array_Decays[32][2] = H0amplitudesdownRantisdownR; ParticleHiggsH.Array_Decays[32][3] = 2; ParticleHiggsH.Array_Comments[32] = "# H -> ~d_R ~d_R*";
   ParticleHiggsH.Array_Decays[33][0] = PDGsdownL; ParticleHiggsH.Array_Decays[33][1] = -PDGsdownR; ParticleHiggsH.Array_Decays[33][2] = H0amplitudesdownLantisdownR; ParticleHiggsH.Array_Decays[33][3] = 2; ParticleHiggsH.Array_Comments[33] = "# H -> ~d_L ~d_R*";
   ParticleHiggsH.Array_Decays[34][0] = PDGsdownR; ParticleHiggsH.Array_Decays[34][1] = -PDGsdownL; ParticleHiggsH.Array_Decays[34][2] = H0amplitudesdownRantisdownL; ParticleHiggsH.Array_Decays[34][3] = 2; ParticleHiggsH.Array_Comments[34] = "# H -> ~d_R ~d_L*";
   ParticleHiggsH.Array_Decays[35][0] = PDGscharmL; ParticleHiggsH.Array_Decays[35][1] = -PDGscharmL; ParticleHiggsH.Array_Decays[35][2] = H0amplitudescharmLantischarmL; ParticleHiggsH.Array_Decays[35][3] = 2; ParticleHiggsH.Array_Comments[35] = "# H -> ~c_L ~c_L*";
   ParticleHiggsH.Array_Decays[36][0] = PDGscharmR; ParticleHiggsH.Array_Decays[36][1] = -PDGscharmR; ParticleHiggsH.Array_Decays[36][2] = H0amplitudescharmRantischarmR; ParticleHiggsH.Array_Decays[36][3] = 2; ParticleHiggsH.Array_Comments[36] = "# H -> ~c_R ~c_R*";
   ParticleHiggsH.Array_Decays[37][0] = PDGscharmL; ParticleHiggsH.Array_Decays[37][1] = -PDGscharmR; ParticleHiggsH.Array_Decays[37][2] = H0amplitudescharmLantischarmR; ParticleHiggsH.Array_Decays[37][3] = 2; ParticleHiggsH.Array_Comments[37] = "# H -> ~c_L ~c_R*";
   ParticleHiggsH.Array_Decays[38][0] = PDGscharmR; ParticleHiggsH.Array_Decays[38][1] = -PDGscharmL; ParticleHiggsH.Array_Decays[38][2] = H0amplitudescharmRantischarmL; ParticleHiggsH.Array_Decays[38][3] = 2; ParticleHiggsH.Array_Comments[38] = "# H -> ~c_R ~c_L*";
   ParticleHiggsH.Array_Decays[39][0] = PDGsstrangeL; ParticleHiggsH.Array_Decays[39][1] = -PDGsstrangeL; ParticleHiggsH.Array_Decays[39][2] = H0amplitudesstrangeLantisstrangeL; ParticleHiggsH.Array_Decays[39][3] = 2; ParticleHiggsH.Array_Comments[39] = "# H -> ~s_L ~s_L*";
   ParticleHiggsH.Array_Decays[40][0] = PDGsstrangeR; ParticleHiggsH.Array_Decays[40][1] = -PDGsstrangeR; ParticleHiggsH.Array_Decays[40][2] = H0amplitudesstrangeRantisstrangeR; ParticleHiggsH.Array_Decays[40][3] = 2; ParticleHiggsH.Array_Comments[40] = "# H -> ~s_R ~s_R*";
   ParticleHiggsH.Array_Decays[41][0] = PDGsstrangeL; ParticleHiggsH.Array_Decays[41][1] = -PDGsstrangeR; ParticleHiggsH.Array_Decays[41][2] = H0amplitudesstrangeLantisstrangeR; ParticleHiggsH.Array_Decays[41][3] = 2; ParticleHiggsH.Array_Comments[41] = "# H -> ~s_L ~s_R*";
   ParticleHiggsH.Array_Decays[42][0] = PDGsstrangeR; ParticleHiggsH.Array_Decays[42][1] = -PDGsstrangeL; ParticleHiggsH.Array_Decays[42][2] = H0amplitudesstrangeRantisstrangeL; ParticleHiggsH.Array_Decays[42][3] = 2; ParticleHiggsH.Array_Comments[42] = "# H -> ~s_R ~s_L*";
   ParticleHiggsH.Array_Decays[43][0] = PDGnuselectronL; ParticleHiggsH.Array_Decays[43][1] = -PDGnuselectronL; ParticleHiggsH.Array_Decays[43][2] = H0amplitudesnueLantisnueL; ParticleHiggsH.Array_Decays[43][3] = 2; ParticleHiggsH.Array_Comments[43] = "# H -> ~nu_eL ~nu_eL*";
   ParticleHiggsH.Array_Decays[44][0] = PDGselectronL; ParticleHiggsH.Array_Decays[44][1] = -PDGselectronL; ParticleHiggsH.Array_Decays[44][2] = H0amplitudeselectronLantiselectronL; ParticleHiggsH.Array_Decays[44][3] = 2; ParticleHiggsH.Array_Comments[44] = "# H -> ~e_L- ~e_L+";
   ParticleHiggsH.Array_Decays[45][0] = PDGselectronR; ParticleHiggsH.Array_Decays[45][1] = -PDGselectronR; ParticleHiggsH.Array_Decays[45][2] = H0amplitudeselectronRantiselectronR; ParticleHiggsH.Array_Decays[45][3] = 2; ParticleHiggsH.Array_Comments[45] = "# H -> ~e_R- ~e_R+";
   ParticleHiggsH.Array_Decays[46][0] = PDGselectronL; ParticleHiggsH.Array_Decays[46][1] = -PDGselectronR; ParticleHiggsH.Array_Decays[46][2] = H0amplitudeselectronLantiselectronR; ParticleHiggsH.Array_Decays[46][3] = 2; ParticleHiggsH.Array_Comments[46] = "# H -> ~e_L- ~e_R+";
   ParticleHiggsH.Array_Decays[47][0] = PDGselectronR; ParticleHiggsH.Array_Decays[47][1] = -PDGselectronL; ParticleHiggsH.Array_Decays[47][2] = H0amplitudeselectronRantiselectronL; ParticleHiggsH.Array_Decays[47][3] = 2; ParticleHiggsH.Array_Comments[47] = "# H -> ~e_R- ~e_L+";
   ParticleHiggsH.Array_Decays[48][0] = PDGnusmuonL; ParticleHiggsH.Array_Decays[48][1] = -PDGnusmuonL; ParticleHiggsH.Array_Decays[48][2] = H0amplitudesnumuLantisnumuL; ParticleHiggsH.Array_Decays[48][3] = 2; ParticleHiggsH.Array_Comments[48] = "# H -> ~nu_muL ~nu_muL*";
   ParticleHiggsH.Array_Decays[49][0] = PDGsmuonL; ParticleHiggsH.Array_Decays[49][1] = -PDGsmuonL; ParticleHiggsH.Array_Decays[49][2] = H0amplitudesmuonLantismuonL; ParticleHiggsH.Array_Decays[49][3] = 2; ParticleHiggsH.Array_Comments[49] = "# H -> ~mu_L- ~mu_L+";
   ParticleHiggsH.Array_Decays[50][0] = PDGsmuonR; ParticleHiggsH.Array_Decays[50][1] = -PDGsmuonR; ParticleHiggsH.Array_Decays[50][2] = H0amplitudesmuonRantismuonR; ParticleHiggsH.Array_Decays[50][3] = 2; ParticleHiggsH.Array_Comments[50] = "# H -> ~mu_R- ~mu_R+";
   ParticleHiggsH.Array_Decays[51][0] = PDGsmuonL; ParticleHiggsH.Array_Decays[51][1] = -PDGsmuonR; ParticleHiggsH.Array_Decays[51][2] = H0amplitudesmuonLantismuonR; ParticleHiggsH.Array_Decays[51][3] = 2; ParticleHiggsH.Array_Comments[51] = "# H -> ~mu_L- ~mu_R+";
   ParticleHiggsH.Array_Decays[52][0] = PDGsmuonR; ParticleHiggsH.Array_Decays[52][1] = -PDGsmuonL; ParticleHiggsH.Array_Decays[52][2] = H0amplitudesmuonRantismuonL; ParticleHiggsH.Array_Decays[52][3] = 2; ParticleHiggsH.Array_Comments[52] = "# H -> ~mu_R- ~mu_L+";	
   ParticleHiggsH.Array_Decays[53][0] = PDGnustauL; ParticleHiggsH.Array_Decays[53][1] = -PDGnustauL; ParticleHiggsH.Array_Decays[53][2] = H0amplitudesnutauLantisnutauL; ParticleHiggsH.Array_Decays[53][3] = 2; ParticleHiggsH.Array_Comments[53] = "# H -> ~nu_tauL ~nu_tauL*";
   ParticleHiggsH.Array_Decays[54][0] = PDGstop1; ParticleHiggsH.Array_Decays[54][1] = -PDGstop1; ParticleHiggsH.Array_Decays[54][2] = H0amplitudestop1antistop1; ParticleHiggsH.Array_Decays[54][3] = 2; ParticleHiggsH.Array_Comments[54] = "# H -> ~t_1 ~t_1*";	  
   ParticleHiggsH.Array_Decays[55][0] = PDGstop2; ParticleHiggsH.Array_Decays[55][1] = -PDGstop2; ParticleHiggsH.Array_Decays[55][2] = H0amplitudestop2antistop2; ParticleHiggsH.Array_Decays[55][3] = 2; ParticleHiggsH.Array_Comments[55] = "# H -> ~t_2 ~t_2*";
   ParticleHiggsH.Array_Decays[56][0] = PDGstop1; ParticleHiggsH.Array_Decays[56][1] = -PDGstop2; ParticleHiggsH.Array_Decays[56][2] = H0amplitudestop1antistop2; ParticleHiggsH.Array_Decays[56][3] = 2; ParticleHiggsH.Array_Comments[56] = "# H -> ~t_1 ~t_2*";	  
   ParticleHiggsH.Array_Decays[57][0] = PDGstop2; ParticleHiggsH.Array_Decays[57][1] = -PDGstop1; ParticleHiggsH.Array_Decays[57][2] = H0amplitudestop2antistop1; ParticleHiggsH.Array_Decays[57][3] = 2; ParticleHiggsH.Array_Comments[57] = "# H -> ~t_2 ~t_1*";
   ParticleHiggsH.Array_Decays[58][0] = PDGsbottom1; ParticleHiggsH.Array_Decays[58][1] = -PDGsbottom1; ParticleHiggsH.Array_Decays[58][2] = H0amplitudesbottom1antisbottom1; ParticleHiggsH.Array_Decays[58][3] = 2; ParticleHiggsH.Array_Comments[58] = "# H -> ~b_1 ~b_1*";	  
   ParticleHiggsH.Array_Decays[59][0] = PDGsbottom2; ParticleHiggsH.Array_Decays[59][1] = -PDGsbottom2; ParticleHiggsH.Array_Decays[59][2] = H0amplitudesbottom2antisbottom2; ParticleHiggsH.Array_Decays[59][3] = 2; ParticleHiggsH.Array_Comments[59] = "# H -> ~b_2 ~b_2*";
   ParticleHiggsH.Array_Decays[60][0] = PDGsbottom1; ParticleHiggsH.Array_Decays[60][1] = -PDGsbottom2; ParticleHiggsH.Array_Decays[60][2] = H0amplitudesbottom1antisbottom2; ParticleHiggsH.Array_Decays[60][3] = 2; ParticleHiggsH.Array_Comments[60] = "# H -> ~b_1 ~b_2*";	  
   ParticleHiggsH.Array_Decays[61][0] = PDGsbottom2; ParticleHiggsH.Array_Decays[61][1] = -PDGsbottom1; ParticleHiggsH.Array_Decays[61][2] = H0amplitudesbottom2antisbottom1; ParticleHiggsH.Array_Decays[61][3] = 2; ParticleHiggsH.Array_Comments[61] = "# H -> ~b_2 ~b_1*";
   ParticleHiggsH.Array_Decays[62][0] = PDGstau1; ParticleHiggsH.Array_Decays[62][1] = -PDGstau1; ParticleHiggsH.Array_Decays[62][2] = H0amplitudestau1antistau1; ParticleHiggsH.Array_Decays[62][3] = 2; ParticleHiggsH.Array_Comments[62] = "# H -> ~tau_1- ~tau_1+";	  
   ParticleHiggsH.Array_Decays[63][0] = PDGstau2; ParticleHiggsH.Array_Decays[63][1] = -PDGstau2; ParticleHiggsH.Array_Decays[63][2] = H0amplitudestau2antistau2; ParticleHiggsH.Array_Decays[63][3] = 2; ParticleHiggsH.Array_Comments[63] = "# H -> ~tau_2- ~tau_2+";
   ParticleHiggsH.Array_Decays[64][0] = PDGstau1; ParticleHiggsH.Array_Decays[64][1] = -PDGstau2; ParticleHiggsH.Array_Decays[64][2] = H0amplitudestau1antistau2; ParticleHiggsH.Array_Decays[64][3] = 2; ParticleHiggsH.Array_Comments[64] = "# H -> ~tau_1- ~tau_2+";	  
   ParticleHiggsH.Array_Decays[65][0] = PDGstau2; ParticleHiggsH.Array_Decays[65][1] = -PDGstau1; ParticleHiggsH.Array_Decays[65][2] = H0amplitudestau2antistau1; ParticleHiggsH.Array_Decays[65][3] = 2; ParticleHiggsH.Array_Comments[65] = "# H -> ~tau_2- ~tau_1+";
   ParticleHiggsH.Array_Decays[66][0] = PDGphoton; ParticleHiggsH.Array_Decays[66][1] = PDGphoton; ParticleHiggsH.Array_Decays[66][2] = H0amplitudegammagamma; ParticleHiggsH.Array_Decays[66][3] = 2; ParticleHiggsH.Array_Comments[66] = "# H -> gamma gamma";
   ParticleHiggsH.Array_Decays[67][0] = PDGgluon; ParticleHiggsH.Array_Decays[67][1] = PDGgluon; ParticleHiggsH.Array_Decays[67][2] = H0amplitudegluongluon; ParticleHiggsH.Array_Decays[67][3] = 2; ParticleHiggsH.Array_Comments[67] = "# H -> gluon gluon";
   ParticleHiggsH.Array_Decays[68][0] = PDGZboson; ParticleHiggsH.Array_Decays[68][1] = PDGphoton; ParticleHiggsH.Array_Decays[68][2] = H0amplitudeZgamma; ParticleHiggsH.Array_Decays[68][3] = 2; ParticleHiggsH.Array_Comments[68] = "# H -> Z gamma";
   ParticleHiggsH.Array_Decays[69][0] = PDGWplus; ParticleHiggsH.Array_Decays[69][1] = -PDGWplus; ParticleHiggsH.Array_Decays[69][2] = H0amplitudeWbosonWboson; ParticleHiggsH.Array_Decays[69][3] = H0WWNDA; ParticleHiggsH.Array_Comments[69] = H0WWcomment;
   ParticleHiggsH.Array_Decays[70][0] = PDGZboson; ParticleHiggsH.Array_Decays[70][1] = -PDGZboson; ParticleHiggsH.Array_Decays[70][2] = H0amplitudeZbosonZboson; ParticleHiggsH.Array_Decays[70][3] = H0ZZNDA; ParticleHiggsH.Array_Comments[70] = H0ZZcomment;
   
   ParticleHiggsH.Array_Decays[71][0] = PDGneutralino1; ParticleHiggsH.Array_Decays[71][1] = PDGneutralino5; ParticleHiggsH.Array_Decays[71][2] = H0amplitudeneutZ1neutZ5; ParticleHiggsH.Array_Decays[71][3] = 2; ParticleHiggsH.Array_Comments[71] = "# H -> ~chi_10 ~chi_50";
   ParticleHiggsH.Array_Decays[72][0] = PDGneutralino1; ParticleHiggsH.Array_Decays[72][1] = PDGneutralino5; ParticleHiggsH.Array_Decays[72][2] = H0amplitudeneutZ2neutZ5; ParticleHiggsH.Array_Decays[72][3] = 2; ParticleHiggsH.Array_Comments[72] = "# H -> ~chi_20 ~chi_50";
   ParticleHiggsH.Array_Decays[73][0] = PDGneutralino1; ParticleHiggsH.Array_Decays[73][1] = PDGneutralino5; ParticleHiggsH.Array_Decays[73][2] = H0amplitudeneutZ3neutZ5; ParticleHiggsH.Array_Decays[73][3] = 2; ParticleHiggsH.Array_Comments[73] = "# H -> ~chi_30 ~chi_50";
   ParticleHiggsH.Array_Decays[74][0] = PDGneutralino1; ParticleHiggsH.Array_Decays[74][1] = PDGneutralino5; ParticleHiggsH.Array_Decays[74][2] = H0amplitudeneutZ4neutZ5; ParticleHiggsH.Array_Decays[74][3] = 2; ParticleHiggsH.Array_Comments[74] = "# H -> ~chi_40 ~chi_50";
   ParticleHiggsH.Array_Decays[75][0] = PDGneutralino1; ParticleHiggsH.Array_Decays[75][1] = PDGneutralino5; ParticleHiggsH.Array_Decays[75][2] = H0amplitudeneutZ5neutZ5; ParticleHiggsH.Array_Decays[75][3] = 2; ParticleHiggsH.Array_Comments[75] = "# H -> ~chi_50 ~chi_50";
   ParticleHiggsH.Array_Decays[76][0] = PDGA0; ParticleHiggsH.Array_Decays[76][1] = PDGA2; ParticleHiggsH.Array_Decays[76][2] = H0amplitudehiggsAhiggsA2; ParticleHiggsH.Array_Decays[76][3] = 2; ParticleHiggsH.Array_Comments[76] = "# H -> A A2";

   ParticleHiggsH.Array_Decays[78][0] = PDGA2; ParticleHiggsH.Array_Decays[78][1] = PDGA2; ParticleHiggsH.Array_Decays[78][2] = H0amplitudehiggsA2higgsA2; ParticleHiggsH.Array_Decays[78][3] = 2; ParticleHiggsH.Array_Comments[78] = "# H -> A2 A2";
   ParticleHiggsH.Array_Decays[79][0] = PDGZboson; ParticleHiggsH.Array_Decays[79][1] = PDGA2; ParticleHiggsH.Array_Decays[79][2] = H0amplitudehiggsA2Zboson; ParticleHiggsH.Array_Decays[79][3] = 2; ParticleHiggsH.Array_Comments[79] = "# H -> A2 Z";
   ParticleHiggsH.Array_Decays[80][0] = PDGHplus; ParticleHiggsH.Array_Decays[80][1] = -PDGWplus; ParticleHiggsH.Array_Decays[80][2] = H0amplitudeWHpm; ParticleHiggsH.Array_Decays[80][3] = 2; ParticleHiggsH.Array_Comments[80] = "# H -> H+- W-+";
   
   for(int i = 0; i<ParticleHiggsH.No_of_Decays; i++) {
     if (ParticleHiggsH.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleHiggsH.Array_Comments[i] << " is negative = " << ParticleHiggsH.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleHiggsH.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }    
   
   double HiggsH_No_1to2_Decays = 0;
   
   HiggsH_No_1to2_Decays = ParticleHiggsH.No_1to2_Decays + ParticleHiggsH.No_NMSSM_Decays; /// As higgsH can't be NLSP as heavier than higgsl
   
   for (int j = 0; j<HiggsH_No_1to2_Decays; j++) {
     ParticleHiggsH.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<HiggsH_No_1to2_Decays; j++) {
     ParticleHiggsH.two_width = ParticleHiggsH.two_width + ParticleHiggsH.Array_Decays[j][2];
   }
   for (int j=HiggsH_No_1to2_Decays; j<ParticleHiggsH.No_of_Decays; j++) {
     ParticleHiggsH.three_width = ParticleHiggsH.three_width + ParticleHiggsH.Array_Decays[j][2];
   }
   
   for(int j=0; j<ParticleHiggsH.No_of_Decays; j++) {
     ParticleHiggsH.Array_Decays[j][4] = 0;
   }
   
   ///Could argue no need for test for nans here as the higgs 1 -> 3 decay formulae are all purely analytic algebraic expressions, therefore no numerical integration is involved so we can't get nans. Will check anyway as possibility of -ve sqrts in kinematics or -ve logs, or infs etc
   if ( ParticleHiggsH.three_width != ParticleHiggsH.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for HiggsH - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleHiggsH.No_of_Decays = HiggsH_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleHiggsH.total_width = ParticleHiggsH.two_width;
     }
   else {
     ParticleHiggsH.total_width = ParticleHiggsH.two_width + ParticleHiggsH.three_width;
   }
   
   if ( ParticleHiggsH.total_width != ParticleHiggsH.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleHiggsH.No_of_Decays; i++) {
       //   fout << i << " " << ParticleHiggsH.Array_Decays[i][2] << endl;
       // }	  
       throw( "nan in H0 (second heaviest higgs) total width \n");
     }
 } 
 
 if (nmssmIsIt == true) {
   
   ///higgsH3 decays
   
   double H03amplitudeuantiu=0, H03amplitudedantid=0, H03amplitudesantis=0, H03amplitudecantic=0, H03amplitudebantib=0, H03amplitudetantit=0, H03amplitudeeantie=0, H03amplitudemuantimu=0, H03amplitudetauantitau=0, H03amplitudeneutZ1neutZ1=0, H03amplitudeneutZ1neutZ2=0, H03amplitudeneutZ1neutZ3=0, H03amplitudeneutZ1neutZ4=0, H03amplitudeneutZ2neutZ2=0, H03amplitudeneutZ2neutZ3=0, H03amplitudeneutZ2neutZ4=0, H03amplitudeneutZ3neutZ3=0, H03amplitudeneutZ3neutZ4=0, H03amplitudeneutZ4neutZ4=0, H03amplitudecharW1charW1=0, H03amplitudecharW1charW2=0, H03amplitudecharW2charW2=0, H03amplitudeh0h0=0, H03amplitudehiggsAhiggsA=0, H03amplitudeHplusHminus=0, H03amplitudehiggsAZboson=0, H03amplitudesupLantisupL=0, H03amplitudesupLantisupR=0, H03amplitudesupRantisupL=0, H03amplitudesupRantisupR=0, H03amplitudesdownLantisdownL=0, H03amplitudesdownLantisdownR=0, H03amplitudesdownRantisdownL=0, H03amplitudesdownRantisdownR=0, H03amplitudescharmLantischarmL=0, H03amplitudescharmLantischarmR=0, H03amplitudescharmRantischarmL=0, H03amplitudescharmRantischarmR=0, H03amplitudesstrangeLantisstrangeL=0, H03amplitudesstrangeLantisstrangeR=0, H03amplitudesstrangeRantisstrangeL=0, H03amplitudesstrangeRantisstrangeR=0, H03amplitudesnueLantisnueL=0, H03amplitudeselectronLantiselectronL=0, H03amplitudeselectronRantiselectronR=0, H03amplitudeselectronLantiselectronR=0, H03amplitudeselectronRantiselectronL=0, H03amplitudesnumuLantisnumuL=0, H03amplitudesnutauLantisnutauL=0, H03amplitudesmuonLantismuonL=0, H03amplitudesmuonRantismuonR=0, H03amplitudesmuonLantismuonR=0, H03amplitudesmuonRantismuonL=0, H03amplitudestau1antistau1=0, H03amplitudestau2antistau2=0, H03amplitudestau1antistau2=0, H03amplitudestau2antistau1=0, H03amplitudestop1antistop1=0, H03amplitudestop1antistop2=0, H03amplitudestop2antistop1=0, H03amplitudestop2antistop2=0, H03amplitudesbottom1antisbottom1=0, H03amplitudesbottom1antisbottom2=0, H03amplitudesbottom2antisbottom1=0, H03amplitudesbottom2antisbottom2=0, H03amplitudegluongluon=0, H03amplitudegammagamma=0, H03amplitudeZgamma=0, H03amplitudeneutZ1neutZ5 = 0, H03amplitudeneutZ2neutZ5 = 0, H03amplitudeneutZ3neutZ5 = 0, H03amplitudeneutZ4neutZ5 = 0, H03amplitudeneutZ5neutZ5 = 0, H03amplitudeWW = 0, H03amplitudeZZ = 0, H03amplitudehiggsAhiggsA2 = 0, H03amplitudehiggsA2higgsA2 = 0, H03amplitudehiggsA2Zboson = 0, H03amplitudeh0H0 = 0, H03amplitudeH0H0 = 0, H03amplitudeWHpm = 0;
   
   if (flagH3 == 1) {
     
     if (QCDcorr == false) {
       ///No decays to u or d as PWs to u and d are tiny as proportional to yukawas squared
       ///Use running masses here to try to approximate some of the correction (which aren't included)
       H03amplitudecantic = higgslorHamplitudedecayquarkantiquark (mh0(3), runmc, g, alpha, beta, runmw, 1, 'M', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH3);
       H03amplitudesantis = higgslorHamplitudedecayquarkantiquark (mh0(3), runms, g, alpha, beta, runmw, 0, 'M', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH3);
       H03amplitudebantib = higgslorHamplitudedecayquarkantiquark (mh0(3), runmb, g, alpha, beta, runmw, 0, 'M', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH3); ///use runmb here to reduce necessary corrections
       H03amplitudetantit = higgslorHamplitudedecayquarkantiquark (mh0(3), runmt, g, alpha, beta, runmw, 1, 'M', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH3); ///may wish to use mtPole here rather than runmt as gives closer to answer with corrections
     }
     else {
       H03amplitudecantic = higgslorHamplitudedecayquarkantiquark (mh0(3), mcpole, g, alpha, beta, runmw, 1, 'M', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH3);
       H03amplitudesantis = higgslorHamplitudedecayquarkantiquark (mh0(3), mspole, g, alpha, beta, runmw, 0, 'M', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH3);
       ///mcpole and mspole set in decays.h, this values used are those appropriate for the scheme used for the h -> qq QCD corrections, as in hdecay
       H03amplitudebantib = higgslorHamplitudedecayquarkantiquark (mh0(3), mbPole, g, alpha, beta, runmw, 0, 'M', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH3); 
       H03amplitudetantit = higgslorHamplitudedecayquarkantiquark (mh0(3), mtPole, g, alpha, beta, runmw, 1, 'M', CPEMix, nmssmIsIt, QCDcorr, alphasAtMH3); 
     }
     
     H03amplitudeeantie = higgslorHamplitudedecayquarkantiquark (mh0(3), runmel, g, alpha, beta, runmw, 0, 'M', CPEMix, nmssmIsIt, false, alphasAtMH3)/3; ///0 as leptons are like down-type quarks, divide by 3 as No of colours is 1 for leptons cf 3 for quarks
     H03amplitudemuantimu = higgslorHamplitudedecayquarkantiquark (mh0(3), runmmu, g, alpha, beta, runmw, 0, 'M', CPEMix, nmssmIsIt, false, alphasAtMH3)/3;
     H03amplitudetauantitau = higgslorHamplitudedecayquarkantiquark (mh0(3), runmtau, g, alpha, beta, runmw, 0, 'M', CPEMix, nmssmIsIt, false, alphasAtMH3)/3;
  
     H03amplitudesupLantisupL = higgsH3amplitudedecay2squarksamehandNMSSM (mh0(3), mu(1,1), mu(1,1), g, gp, alpha, beta, runmw, runmu, CPEMix, 1);
     H03amplitudesupRantisupR = higgsH3amplitudedecay2squarksamehandNMSSM (mh0(3), mu(2,1), mu(2,1), g, gp, alpha, beta, runmw, runmu, CPEMix, 3);
     H03amplitudesdownLantisdownL = higgsH3amplitudedecay2squarksamehandNMSSM (mh0(3), md(1,1), md(1,1), g, gp, alpha, beta, runmw, runmd, CPEMix, 2);
     H03amplitudesdownRantisdownR = higgsH3amplitudedecay2squarksamehandNMSSM (mh0(3), md(2,1), mu(2,1), g, gp, alpha, beta, runmw, runmd, CPEMix, 4);
     H03amplitudescharmLantischarmL = higgsH3amplitudedecay2squarksamehandNMSSM (mh0(3), mu(1,2), mu(1,2), g, gp, alpha, beta, runmw, runmc, CPEMix, 1);
     H03amplitudescharmRantischarmR = higgsH3amplitudedecay2squarksamehandNMSSM (mh0(3), mu(2,2), mu(2,2), g, gp, alpha, beta, runmw, runmc, CPEMix, 3);
     H03amplitudesstrangeLantisstrangeL = higgsH3amplitudedecay2squarksamehandNMSSM (mh0(3), md(2,1), md(2,1), g, gp, alpha, beta, runmw, runms, CPEMix, 2);
     H03amplitudesstrangeRantisstrangeR = higgsH3amplitudedecay2squarksamehandNMSSM (mh0(3), md(2,2), md(2,2), g, gp, alpha, beta, runmw, runms, CPEMix, 4);
     
     H03amplitudesnueLantisnueL = higgsH3amplitudedecay2sleptonsamehandNMSSM (mh0(3), msnu(1), msnu(1), g, gp, alpha, beta, runmw, 0, CPEMix, 1);
     H03amplitudeselectronLantiselectronL = higgsH3amplitudedecay2sleptonsamehandNMSSM (mh0(3), me(1,1), me(1,1), g, gp, alpha, beta, runmw, runmel, CPEMix, 2);
     H03amplitudeselectronRantiselectronR = higgsH3amplitudedecay2sleptonsamehandNMSSM (mh0(3), me(2,1), me(2,1), g, gp, alpha, beta, runmw, runmel, CPEMix, 3);
     H03amplitudesnumuLantisnumuL = higgsH3amplitudedecay2sleptonsamehandNMSSM (mh0(3), msnu(2), msnu(2), g, gp, alpha, beta, runmw, 0, CPEMix, 1);
     H03amplitudesmuonLantismuonL = higgsH3amplitudedecay2sleptonsamehandNMSSM (mh0(3), me(1,2), me(1,2), g, gp, alpha, beta, runmw, runmmu, CPEMix, 2);
     H03amplitudesmuonRantismuonR = higgsH3amplitudedecay2sleptonsamehandNMSSM (mh0(3), me(2,2), me(2,2), g, gp, alpha, beta, runmw, runmmu, CPEMix, 3);
     H03amplitudesnutauLantisnutauL = higgsH3amplitudedecay2sleptonsamehandNMSSM (mh0(3), msnu(3), msnu(3), g, gp, alpha, beta, runmw, 0, CPEMix, 1);
     
     H03amplitudesupLantisupR = higgshamplitudedecay2squarkdiffhandNMSSM(mh0(3), mu(1,1), mu(2,1), g, alpha, beta, runmw, runmu, Au, mueff, lam, CPEMix, 1, 3);
     H03amplitudesupRantisupL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), mu(1,1), mu(2,1), g, alpha, beta, runmw, runmu, Au, mueff, lam, CPEMix, 1, 3);
     H03amplitudesdownLantisdownR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), md(1,1), md(2,1), g, alpha, beta, runmw, runmd, Ad, mueff, lam, CPEMix, 2, 3);
     H03amplitudesdownRantisdownL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), md(1,1), md(2,1), g, alpha, beta, runmw, runmd, Ad, mueff, lam, CPEMix, 2, 3);
     H03amplitudescharmLantischarmR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), mu(1,2), mu(2,2), g, alpha, beta, runmw, runmc, Ac, mueff, lam, CPEMix, 1, 3);
     H03amplitudescharmRantischarmL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), mu(1,2), mu(2,2), g, alpha, beta, runmw, runmc, Ac, mueff, lam, CPEMix, 1, 3);
     H03amplitudesstrangeLantisstrangeR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), md(2,1), md(2,2), g, alpha, beta, runmw, runms, As, mueff, lam, CPEMix, 2, 3);
     H03amplitudesstrangeRantisstrangeL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), md(2,1), md(2,2), g, alpha, beta, runmw, runms, As, mueff, lam, CPEMix, 2, 3);
     H03amplitudeselectronLantiselectronR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), me(1,1), me(2,1), g, alpha, beta, runmw, runmel, Ae, mueff, lam, CPEMix, 2, 3)/3;
     H03amplitudeselectronRantiselectronL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), me(2,1), me(1,1), g, alpha, beta, runmw, runmel, Ae, mueff, lam, CPEMix, 2, 3)/3;
     H03amplitudesmuonLantismuonR = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), me(1,2), me(2,2), g, alpha, beta, runmw, runmmu, Amu, mueff, lam, CPEMix, 2, 3)/3;
     H03amplitudesmuonRantismuonL = higgshamplitudedecay2squarkdiffhandNMSSM (mh0(3), me(2,2), me(2,1), g, alpha, beta, runmw, runmmu, Amu, mueff, lam, CPEMix, 2, 3)/3;
     H03amplitudestop1antistop1 = higgsCPevenamplitudedecaystopistopiNMSSM (mh0(2), mu(1,3), thetat, runmt, g, gp, runmw, beta, CPEMix, At, mueff, lam, 1, 3);
     H03amplitudestop2antistop2 = higgsCPevenamplitudedecaystopistopiNMSSM (mh0(2), mu(2,3), thetat, runmt, g, gp, runmw, beta, CPEMix, At, mueff, lam, 2, 3) ;
     H03amplitudestop1antistop2 = higgsCPevenamplitudedecaystopistopjNMSSM (mh0(3), mu(1,3), mu(2,3), thetat, runmt, g, gp, runmw, beta, CPEMix, At, mueff, lam, 3); 
     H03amplitudestop2antistop1 = h0amplitudestop1antistop2;
     H03amplitudesbottom1antisbottom1 = higgsCPevenamplitudedecaysbottomisbottomiNMSSM (mh0(3), md(1,3), thetab, runmb, g, gp, runmw, beta,  CPEMix, Ab, mueff, lam, 1, 3);
     H03amplitudesbottom2antisbottom2 = higgsCPevenamplitudedecaysbottomisbottomiNMSSM (mh0(3), md(2,3), thetab, runmb, g, gp, runmw, beta,  CPEMix, Ab, mueff, lam, 2, 3);
     H03amplitudesbottom1antisbottom2 = higgsCPevenamplitudedecaysbottomisbottomjNMSSM (mh0(3), md(1,3), md(2,3), thetab, runmb, g, gp, runmw, beta, CPEMix, Ab, mueff, lam, 3);
     H03amplitudesbottom2antisbottom1 = h0amplitudesbottom1antisbottom2;
     H03amplitudestau1antistau1 = higgsCPevenamplitudedecaystauistauiNMSSM (mh0(3), me(1,3), thetatau - PI/2, runmtau, g, gp, runmw, beta,  CPEMix, Atau, mueff, lam, 1, 3); 
     H03amplitudestau2antistau2 = higgsCPevenamplitudedecaystauistauiNMSSM (mh0(3), me(2,3), thetatau - PI/2, runmtau, g, gp, runmw, beta,  CPEMix, Atau, mueff, lam, 2, 3);
     H03amplitudestau1antistau2 = higgsCPevenamplitudedecaystauistaujNMSSM (mh0(3), me(1,3), me(2,3), thetatau - PI/2, runmtau, g, gp, runmw, beta, CPEMix, Atau, mueff, lam, 3);
     H03amplitudestau2antistau1 = h0amplitudestau1antistau2;
     
     H03amplitudecharW1charW1 = higgsphiamplitudedecaysamecharginoNMSSM (mh0(3), mch(1), g, thetaL2, thetaR2, lam, CPEMix, 1, 3);
     H03amplitudecharW2charW2 = higgsphiamplitudedecaysamecharginoNMSSM (mh0(3), mch(2), g, thetaL2, thetaR2, lam, CPEMix, 2, 3);
     H03amplitudecharW1charW2 = higgsphiamplitudedecaydiffcharginoNMSSM (mh0(3), -mch(1), mch(2), g, thetaL2, thetaR2, lam, CPEMix, 3);
     
     H03amplitudegammagamma = higgsCPevenamplitudedecaygammagammaNMSSM(mh0(3), mtAtMH3, mbAtMH3, mcAtMH3, runmtau, runmw, mHpm, mch(1), mch(2), mu(1,2), mu(2,2), mu(1,3), mu(2,3), md(1,2), md(2,2), md(1,3), md(2,3), me(1,2), me(2,2), me(1,3), me(2,3), CPEMix, beta, g, gp, alphaAtMH3, thetat, thetab, thetatau-PI/2, thetaL2, thetaR2, At, Ab, Atau, greekmu, mueff, lam, kappa, Alambda, 3); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     H03amplitudegluongluon = higgsCPevenamplitudedecaygluongluonNMSSM(mh0(3), mtPole, mbPole, mcpole, runmw, mu(1,2), mu(2,2), mu(1,3), mu(2,3), md(1,2), md(2,2), md(1,3), md(2,3), mu(1,1), mu(2,1), md(1,1), md(2,1), mtAtMH3, mbAtMH3, CPEMix, beta, g, gp, gs, alphasAtMH3, thetat, thetab, thetaL2, thetaR2, At, Ab, greekmu, mueff, lam, kappa, Alambda, 3, QCDcorr); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     H03amplitudeZgamma = higgshamplitudedecayZgammaNMSSM (mh0(3), g, gp, alphaAtMH3, runmw, polemz, mHpm, CPEMix, beta, mtAtMH3, mbAtMH3, mcAtMH3, mch(1), mch(2), thetaL2, thetaR2, lam, kappa, Alambda, greekmu, mueff, 3); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     
     H03amplitudeneutZ1neutZ1 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(1), mneut(1), g, gp, CPEMix, mixNeut, lam, kappa, 1, 1, 3);
     H03amplitudeneutZ1neutZ2 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(1), mneut(2), g, gp, CPEMix, mixNeut, lam, kappa, 1, 2, 3);
     H03amplitudeneutZ1neutZ3 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(1), mneut(3), g, gp, CPEMix, mixNeut, lam, kappa, 1, 3, 3);
     H03amplitudeneutZ1neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(1), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 1, 4, 3);
     H03amplitudeneutZ1neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(1), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 1, 5, 3);
     H03amplitudeneutZ2neutZ2 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(2), mneut(2), g, gp, CPEMix, mixNeut, lam, kappa, 2, 2, 3);
     H03amplitudeneutZ2neutZ3 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(2), mneut(3), g, gp, CPEMix, mixNeut, lam, kappa, 2, 3, 3);
     H03amplitudeneutZ2neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(2), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 2, 4, 3);
     H03amplitudeneutZ2neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(2), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 2, 5, 3);
     H03amplitudeneutZ3neutZ3 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(3), mneut(3), g, gp, CPEMix, mixNeut, lam, kappa, 3, 3, 3);
     H03amplitudeneutZ3neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(3), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 3, 4, 3);
     H03amplitudeneutZ3neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(3), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 3, 5, 3);
     H03amplitudeneutZ4neutZ4 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(4), mneut(4), g, gp, CPEMix, mixNeut, lam, kappa, 4, 4, 3);
     H03amplitudeneutZ4neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(4), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 4, 5, 3);
     H03amplitudeneutZ5neutZ5 = higgshamplitudedecayneutineutjNMSSM (mh0(3), mneut(5), mneut(5), g, gp, CPEMix, mixNeut, lam, kappa, 5, 5, 3);
     
     H03amplitudehiggsAhiggsA = higgsCPevenamplitudedecayAANMSSM(mh0(3), mA0(1), mA0(1), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 3, 1, 1);
     H03amplitudehiggsAhiggsA2 = 2*higgsCPevenamplitudedecayAANMSSM(mh0(3), mA0(1), mA0(2), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 3, 1, 2);
     H03amplitudehiggsA2higgsA2 = higgsCPevenamplitudedecayAANMSSM(mh0(3), mA0(2), mA0(2), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 3, 2, 2);
      H03amplitudehiggsAZboson = higgsCPevenamplitudedecaypseudoscalarZNMSSM(mh0(3), mA0(1), polemz, g, gp, beta, CPEMix, CPOMix, 3, 1);
      H03amplitudehiggsA2Zboson = higgsCPevenamplitudedecaypseudoscalarZNMSSM(mh0(3), mA0(2), polemz, g, gp, beta, CPEMix, CPOMix, 3, 2);
      H03amplitudeHplusHminus = higgsCPevenamplitudedecayHpHmNMSSM (mh0(3), mHpm, runmw, g, gp, runmt, runmb, beta, lam, mueff, kappa, Alambda, CPEMix, 3);
      H03amplitudeh0h0 = higgsCPevenamplitudedecayhhorhHorHHNMSSM(mh0(3), mh0(1), mh0(1), g, gp, runmw, runmt, runmb, beta, lam,  Alambda, kappa, Akappa, mueff, CPEMix, CPOMix, 1, 1, 3);
      H03amplitudeh0H0 = higgsCPevenamplitudedecayhhorhHorHHNMSSM(mh0(3), mh0(1), mh0(2), g, gp, runmw, runmt, runmb, beta, lam,  Alambda, kappa, Akappa, mueff, CPEMix, CPOMix, 1, 2, 3);
      H03amplitudeH0H0 = higgsCPevenamplitudedecayhhorhHorHHNMSSM(mh0(3), mh0(2), mh0(2), g, gp, runmw, runmt, runmb, beta, lam,  Alambda, kappa, Akappa, mueff, CPEMix, CPOMix, 2, 2, 3);
      H03amplitudeWHpm = higgsCPevenamplitudedecayWHpmNMSSM (mh0(3), polemw, mHpm, beta, g, CPEMix, 3)*2; ///*2 as W+H- or W-H+
      H03amplitudeWW = higgsH3amplitudedecayVVNMSSM(mh0(3), polemw, polemz, g, gp, alpha, beta, 'W', CPEMix, nmssmIsIt)(1);
      H03amplitudeZZ = higgsH3amplitudedecayVVNMSSM(mh0(3), polemw, polemz, g, gp, alpha, beta, 'Z', CPEMix, nmssmIsIt)(1);
      
      int H03WWcommentcode, H03ZZcommentcode, H03WWNDA=0, H03ZZNDA=0;
      H03WWcommentcode = higgsH3amplitudedecayVVNMSSM(mh0(3), polemw, polemz, g, gp, alpha, beta, 'W', CPEMix, nmssmIsIt)(2);
      string H03WWcomment, H03ZZcomment;
      if (H03WWcommentcode == 1) {
	H03WWcomment = "# H3 -> WW* -> W f f'bar";
	H03WWNDA = 2; ///So read into other programs, e.g. PYTHIA, correctly
      }
      else if(H03WWcommentcode == 2) {
	H03WWcomment = "# H3 -> W+ W-";
	H03WWNDA = 2;
      }
      H03ZZcommentcode = higgsH3amplitudedecayVVNMSSM(mh0(3), polemw, polemz, g, gp, alpha, beta, 'Z', CPEMix, nmssmIsIt)(2);
      if (H03ZZcommentcode == 1) {
	H03ZZcomment = "# H3 -> ZZ* -> Z f f'bar";
	H03ZZNDA = 2; ///So read into other programs, e.g. PYTHIA, correctly
      }
      else if(H03ZZcommentcode == 2) {
	H03ZZcomment = "# H3 -> Z Z";
	H03ZZNDA = 2; 
      }
      
      
      ParticleHiggsH3.Array_Decays[0][0] = PDGup; ParticleHiggsH3.Array_Decays[0][1] = -PDGup; ParticleHiggsH3.Array_Decays[0][2] = H03amplitudeuantiu; ParticleHiggsH3.Array_Decays[0][3] = 2; ParticleHiggsH3.Array_Comments[0] = "# H3 -> u ub";
      ParticleHiggsH3.Array_Decays[1][0] = PDGdown; ParticleHiggsH3.Array_Decays[1][1] = -PDGdown; ParticleHiggsH3.Array_Decays[1][2] = H03amplitudedantid; ParticleHiggsH3.Array_Decays[1][3] = 2; ParticleHiggsH3.Array_Comments[1] = "# H3 -> d db";
      ParticleHiggsH3.Array_Decays[2][0] = PDGcharm; ParticleHiggsH3.Array_Decays[2][1] = -PDGcharm; ParticleHiggsH3.Array_Decays[2][2] = H03amplitudecantic; ParticleHiggsH3.Array_Decays[2][3] = 2; ParticleHiggsH3.Array_Comments[2] = "# H3 -> c cb";
      ParticleHiggsH3.Array_Decays[3][0] = PDGstrange; ParticleHiggsH3.Array_Decays[3][1] = -PDGstrange; ParticleHiggsH3.Array_Decays[3][2] = H03amplitudesantis; ParticleHiggsH3.Array_Decays[3][3] = 2; ParticleHiggsH3.Array_Comments[3] = "# H3 -> s sb";
      ParticleHiggsH3.Array_Decays[4][0] = PDGbottom; ParticleHiggsH3.Array_Decays[4][1] = -PDGbottom; ParticleHiggsH3.Array_Decays[4][2] = H03amplitudebantib; ParticleHiggsH3.Array_Decays[4][3] = 2; ParticleHiggsH3.Array_Comments[4] = "# H3 -> b bb";
      ParticleHiggsH3.Array_Decays[5][0] = PDGtop; ParticleHiggsH3.Array_Decays[5][1] = -PDGtop; ParticleHiggsH3.Array_Decays[5][2] = H03amplitudetantit; ParticleHiggsH3.Array_Decays[5][3] = 2; ParticleHiggsH3.Array_Comments[5] = "# H3 -> t tb";
      ParticleHiggsH3.Array_Decays[6][0] = PDGelectron; ParticleHiggsH3.Array_Decays[6][1] = -PDGelectron; ParticleHiggsH3.Array_Decays[6][2] = H03amplitudeeantie; ParticleHiggsH3.Array_Decays[6][3] = 2; ParticleHiggsH3.Array_Comments[6] = "# H3 -> e- e+";
      ParticleHiggsH3.Array_Decays[7][0] = PDGmuon; ParticleHiggsH3.Array_Decays[7][1] = -PDGmuon; ParticleHiggsH3.Array_Decays[7][2] = H03amplitudemuantimu; ParticleHiggsH3.Array_Decays[7][3] = 2; ParticleHiggsH3.Array_Comments[7] = "# H3 -> mu- mu+";
      ParticleHiggsH3.Array_Decays[8][0] = PDGtau; ParticleHiggsH3.Array_Decays[8][1] = -PDGtau; ParticleHiggsH3.Array_Decays[8][2] = H03amplitudetauantitau; ParticleHiggsH3.Array_Decays[8][3] = 2; ParticleHiggsH3.Array_Comments[8] = "# H3 -> tau- tau+";
      
      ParticleHiggsH3.Array_Decays[9][0] = PDGneutralino1; ParticleHiggsH3.Array_Decays[9][1] = PDGneutralino1; ParticleHiggsH3.Array_Decays[9][2] = H03amplitudeneutZ1neutZ1; ParticleHiggsH3.Array_Decays[9][3] = 2; ParticleHiggsH3.Array_Comments[9] = "# H3 -> ~chi_10 ~chi_10";
      ParticleHiggsH3.Array_Decays[10][0] = PDGneutralino1; ParticleHiggsH3.Array_Decays[10][1] = PDGneutralino2; ParticleHiggsH3.Array_Decays[10][2] = H03amplitudeneutZ1neutZ2; ParticleHiggsH3.Array_Decays[10][3] = 2; ParticleHiggsH3.Array_Comments[10] = "# H3 -> ~chi_10 ~chi_20";
      ParticleHiggsH3.Array_Decays[11][0] = PDGneutralino1; ParticleHiggsH3.Array_Decays[11][1] = PDGneutralino3; ParticleHiggsH3.Array_Decays[11][2] = H03amplitudeneutZ1neutZ3; ParticleHiggsH3.Array_Decays[11][3] = 2; ParticleHiggsH3.Array_Comments[11] = "# H3 -> ~chi_10 ~chi_30";
      ParticleHiggsH3.Array_Decays[12][0] = PDGneutralino1; ParticleHiggsH3.Array_Decays[12][1] = PDGneutralino4; ParticleHiggsH3.Array_Decays[12][2] = H03amplitudeneutZ1neutZ4; ParticleHiggsH3.Array_Decays[12][3] = 2; ParticleHiggsH3.Array_Comments[12] = "# H3 -> ~chi_10 ~chi_40";
      ParticleHiggsH3.Array_Decays[13][0] = PDGneutralino2; ParticleHiggsH3.Array_Decays[13][1] = PDGneutralino2; ParticleHiggsH3.Array_Decays[13][2] = H03amplitudeneutZ2neutZ2; ParticleHiggsH3.Array_Decays[13][3] = 2; ParticleHiggsH3.Array_Comments[13] = "# H3 -> ~chi_20 ~chi_20";
      ParticleHiggsH3.Array_Decays[14][0] = PDGneutralino2; ParticleHiggsH3.Array_Decays[14][1] = PDGneutralino3; ParticleHiggsH3.Array_Decays[14][2] = H03amplitudeneutZ2neutZ3; ParticleHiggsH3.Array_Decays[14][3] = 2; ParticleHiggsH3.Array_Comments[14] = "# H3 -> ~chi_20 ~chi_30";
      ParticleHiggsH3.Array_Decays[15][0] = PDGneutralino2; ParticleHiggsH3.Array_Decays[15][1] = PDGneutralino4; ParticleHiggsH3.Array_Decays[15][2] = H03amplitudeneutZ2neutZ4; ParticleHiggsH3.Array_Decays[15][3] = 2; ParticleHiggsH3.Array_Comments[15] = "# H3 -> ~chi_20 ~chi_40";
      ParticleHiggsH3.Array_Decays[16][0] = PDGneutralino3; ParticleHiggsH3.Array_Decays[16][1] = PDGneutralino3; ParticleHiggsH3.Array_Decays[16][2] = H03amplitudeneutZ3neutZ3; ParticleHiggsH3.Array_Decays[16][3] = 2; ParticleHiggsH3.Array_Comments[16] = "# H3 -> ~chi_30 ~chi_30";
      ParticleHiggsH3.Array_Decays[17][0] = PDGneutralino3; ParticleHiggsH3.Array_Decays[17][1] = PDGneutralino4; ParticleHiggsH3.Array_Decays[17][2] = H03amplitudeneutZ3neutZ4; ParticleHiggsH3.Array_Decays[17][3] = 2; ParticleHiggsH3.Array_Comments[17] = "# H3 -> ~chi_30 ~chi_40";
      ParticleHiggsH3.Array_Decays[18][0] = PDGneutralino4; ParticleHiggsH3.Array_Decays[18][1] = PDGneutralino4; ParticleHiggsH3.Array_Decays[18][2] = H03amplitudeneutZ4neutZ4; ParticleHiggsH3.Array_Decays[18][3] = 2; ParticleHiggsH3.Array_Comments[18] = "# H3 -> ~chi_40 ~chi_40";
      
      ParticleHiggsH3.Array_Decays[19][0] = PDGchargino1; ParticleHiggsH3.Array_Decays[19][1] = -PDGchargino1; ParticleHiggsH3.Array_Decays[19][2] = H03amplitudecharW1charW1; ParticleHiggsH3.Array_Decays[19][3] = 2; ParticleHiggsH3.Array_Comments[19] = "# H3 -> ~chi_1+ ~chi_1-";
      ParticleHiggsH3.Array_Decays[20][0] = PDGchargino2; ParticleHiggsH3.Array_Decays[20][1] = -PDGchargino2; ParticleHiggsH3.Array_Decays[20][2] = H03amplitudecharW2charW2; ParticleHiggsH3.Array_Decays[20][3] = 2; ParticleHiggsH3.Array_Comments[20] = "# H3 -> ~chi_2+ ~chi_2-";
      ParticleHiggsH3.Array_Decays[21][0] = PDGchargino1; ParticleHiggsH3.Array_Decays[21][1] = -PDGchargino2; ParticleHiggsH3.Array_Decays[21][2] = H03amplitudecharW1charW2; ParticleHiggsH3.Array_Decays[21][3] = 2; ParticleHiggsH3.Array_Comments[21] = "# H3 -> ~chi_1+ ~chi_2-";
      ParticleHiggsH3.Array_Decays[22][0] = PDGchargino2; ParticleHiggsH3.Array_Decays[22][1] = -PDGchargino1; ParticleHiggsH3.Array_Decays[22][2] = H03amplitudecharW1charW2; ParticleHiggsH3.Array_Decays[22][3] = 2; ParticleHiggsH3.Array_Comments[22] = "# H3 -> ~chi_2+ ~chi_1-"; ///amplitude same as decay to W1+ and W2- by CP invariance
      
      ParticleHiggsH3.Array_Decays[23][0] = PDGh0; ParticleHiggsH3.Array_Decays[23][1] = PDGh0; ParticleHiggsH3.Array_Decays[23][2] = H03amplitudeh0h0; ParticleHiggsH3.Array_Decays[23][3] = 2; ParticleHiggsH3.Array_Comments[23] = "# H3 -> h h";
      ParticleHiggsH3.Array_Decays[24][0] = PDGA0; ParticleHiggsH3.Array_Decays[24][1] = PDGA0; ParticleHiggsH3.Array_Decays[24][2] = H03amplitudehiggsAhiggsA; ParticleHiggsH3.Array_Decays[24][3] = 2; ParticleHiggsH3.Array_Comments[24] = "# H3 -> A A";
      ParticleHiggsH3.Array_Decays[25][0] = PDGHplus; ParticleHiggsH3.Array_Decays[25][1] = -PDGHplus; ParticleHiggsH3.Array_Decays[25][2] = H03amplitudeHplusHminus; ParticleHiggsH3.Array_Decays[25][3] = 2; ParticleHiggsH3.Array_Comments[25] = "# H3 -> H+ H-";
      ParticleHiggsH3.Array_Decays[26][0] = PDGA0; ParticleHiggsH3.Array_Decays[26][1] = PDGZboson; ParticleHiggsH3.Array_Decays[26][2] = H03amplitudehiggsAZboson; ParticleHiggsH3.Array_Decays[26][3] = 2; ParticleHiggsH3.Array_Comments[26] = "# H3 -> A Z";
      
      ParticleHiggsH3.Array_Decays[27][0] = PDGsupL; ParticleHiggsH3.Array_Decays[27][1] = -PDGsupL; ParticleHiggsH3.Array_Decays[27][2] = H03amplitudesupLantisupL; ParticleHiggsH3.Array_Decays[27][3] = 2; ParticleHiggsH3.Array_Comments[27] = "# H3 -> ~u_L ~u_L*";
      ParticleHiggsH3.Array_Decays[28][0] = PDGsupR; ParticleHiggsH3.Array_Decays[28][1] = -PDGsupR; ParticleHiggsH3.Array_Decays[28][2] = H03amplitudesupRantisupR; ParticleHiggsH3.Array_Decays[28][3] = 2; ParticleHiggsH3.Array_Comments[28] = "# H3 -> ~u_R ~u_R*";
      ParticleHiggsH3.Array_Decays[29][0] = PDGsupL; ParticleHiggsH3.Array_Decays[29][1] = -PDGsupR; ParticleHiggsH3.Array_Decays[29][2] = H03amplitudesupLantisupR; ParticleHiggsH3.Array_Decays[29][3] = 2; ParticleHiggsH3.Array_Comments[29] = "# H3 -> ~u_L ~u_R*";
      ParticleHiggsH3.Array_Decays[30][0] = PDGsupR; ParticleHiggsH3.Array_Decays[30][1] = -PDGsupL; ParticleHiggsH3.Array_Decays[30][2] = H03amplitudesupRantisupL; ParticleHiggsH3.Array_Decays[30][3] = 2; ParticleHiggsH3.Array_Comments[30] = "# H3 -> ~u_R ~u_L*";
      ParticleHiggsH3.Array_Decays[31][0] = PDGsdownL; ParticleHiggsH3.Array_Decays[31][1] = -PDGsdownL; ParticleHiggsH3.Array_Decays[31][2] = H03amplitudesdownLantisdownL; ParticleHiggsH3.Array_Decays[31][3] = 2; ParticleHiggsH3.Array_Comments[31] = "# H3 -> ~d_L ~d_L*";
      ParticleHiggsH3.Array_Decays[32][0] = PDGsdownR; ParticleHiggsH3.Array_Decays[32][1] = -PDGsdownR; ParticleHiggsH3.Array_Decays[32][2] = H03amplitudesdownRantisdownR; ParticleHiggsH3.Array_Decays[32][3] = 2; ParticleHiggsH3.Array_Comments[32] = "# H3 -> ~d_R ~d_R*";
      ParticleHiggsH3.Array_Decays[33][0] = PDGsdownL; ParticleHiggsH3.Array_Decays[33][1] = -PDGsdownR; ParticleHiggsH3.Array_Decays[33][2] = H03amplitudesdownLantisdownR; ParticleHiggsH3.Array_Decays[33][3] = 2; ParticleHiggsH3.Array_Comments[33] = "# H3 -> ~d_L ~d_R*";
      ParticleHiggsH3.Array_Decays[34][0] = PDGsdownR; ParticleHiggsH3.Array_Decays[34][1] = -PDGsdownL; ParticleHiggsH3.Array_Decays[34][2] = H03amplitudesdownRantisdownL; ParticleHiggsH3.Array_Decays[34][3] = 2; ParticleHiggsH3.Array_Comments[34] = "# H3 -> ~d_R ~d_L*";
      ParticleHiggsH3.Array_Decays[35][0] = PDGscharmL; ParticleHiggsH3.Array_Decays[35][1] = -PDGscharmL; ParticleHiggsH3.Array_Decays[35][2] = H03amplitudescharmLantischarmL; ParticleHiggsH3.Array_Decays[35][3] = 2; ParticleHiggsH3.Array_Comments[35] = "# H3 -> ~c_L ~c_L*";
      ParticleHiggsH3.Array_Decays[36][0] = PDGscharmR; ParticleHiggsH3.Array_Decays[36][1] = -PDGscharmR; ParticleHiggsH3.Array_Decays[36][2] = H03amplitudescharmRantischarmR; ParticleHiggsH3.Array_Decays[36][3] = 2; ParticleHiggsH3.Array_Comments[36] = "# H3 -> ~c_R ~c_R*";
      ParticleHiggsH3.Array_Decays[37][0] = PDGscharmL; ParticleHiggsH.Array_Decays[37][1] = -PDGscharmR; ParticleHiggsH3.Array_Decays[37][2] = H03amplitudescharmLantischarmR; ParticleHiggsH3.Array_Decays[37][3] = 2; ParticleHiggsH3.Array_Comments[37] = "# H3 -> ~c_L ~c_R*";
      ParticleHiggsH3.Array_Decays[38][0] = PDGscharmR; ParticleHiggsH3.Array_Decays[38][1] = -PDGscharmL; ParticleHiggsH3.Array_Decays[38][2] = H03amplitudescharmRantischarmL; ParticleHiggsH3.Array_Decays[38][3] = 2; ParticleHiggsH3.Array_Comments[38] = "# H3 -> ~c_R ~c_L*";
      ParticleHiggsH3.Array_Decays[39][0] = PDGsstrangeL; ParticleHiggsH3.Array_Decays[39][1] = -PDGsstrangeL; ParticleHiggsH3.Array_Decays[39][2] = H03amplitudesstrangeLantisstrangeL; ParticleHiggsH3.Array_Decays[39][3] = 2; ParticleHiggsH3.Array_Comments[39] = "# H3 -> ~s_L ~s_L*";
      ParticleHiggsH3.Array_Decays[40][0] = PDGsstrangeR; ParticleHiggsH3.Array_Decays[40][1] = -PDGsstrangeR; ParticleHiggsH3.Array_Decays[40][2] = H03amplitudesstrangeRantisstrangeR; ParticleHiggsH3.Array_Decays[40][3] = 2; ParticleHiggsH3.Array_Comments[40] = "# H3 -> ~s_R ~s_R*";
      ParticleHiggsH3.Array_Decays[41][0] = PDGsstrangeL; ParticleHiggsH3.Array_Decays[41][1] = -PDGsstrangeR; ParticleHiggsH3.Array_Decays[41][2] = H03amplitudesstrangeLantisstrangeR; ParticleHiggsH3.Array_Decays[41][3] = 2; ParticleHiggsH3.Array_Comments[41] = "# H3 -> ~s_L ~s_R*";
      ParticleHiggsH3.Array_Decays[42][0] = PDGsstrangeR; ParticleHiggsH3.Array_Decays[42][1] = -PDGsstrangeL; ParticleHiggsH3.Array_Decays[42][2] = H03amplitudesstrangeRantisstrangeL; ParticleHiggsH3.Array_Decays[42][3] = 2; ParticleHiggsH3.Array_Comments[42] = "# H3 -> ~s_R ~s_L*";
      
      ParticleHiggsH3.Array_Decays[43][0] = PDGnuselectronL; ParticleHiggsH3.Array_Decays[43][1] = -PDGnuselectronL; ParticleHiggsH3.Array_Decays[43][2] = H03amplitudesnueLantisnueL; ParticleHiggsH3.Array_Decays[43][3] = 2; ParticleHiggsH3.Array_Comments[43] = "# H3 -> ~nu_eL ~nu_eL*";
      ParticleHiggsH3.Array_Decays[44][0] = PDGselectronL; ParticleHiggsH3.Array_Decays[44][1] = -PDGselectronL; ParticleHiggsH3.Array_Decays[44][2] = H03amplitudeselectronLantiselectronL; ParticleHiggsH3.Array_Decays[44][3] = 2; ParticleHiggsH3.Array_Comments[44] = "# H3 -> ~e_L- ~e_L+";
      ParticleHiggsH3.Array_Decays[45][0] = PDGselectronR; ParticleHiggsH3.Array_Decays[45][1] = -PDGselectronR; ParticleHiggsH3.Array_Decays[45][2] = H03amplitudeselectronRantiselectronR; ParticleHiggsH3.Array_Decays[45][3] = 2; ParticleHiggsH3.Array_Comments[45] = "# H3 -> ~e_R- ~e_R+";
      ParticleHiggsH3.Array_Decays[46][0] = PDGselectronL; ParticleHiggsH3.Array_Decays[46][1] = -PDGselectronR; ParticleHiggsH3.Array_Decays[46][2] = H03amplitudeselectronLantiselectronR; ParticleHiggsH3.Array_Decays[46][3] = 2; ParticleHiggsH3.Array_Comments[46] = "# H3 -> ~e_L- ~e_R+";
      ParticleHiggsH3.Array_Decays[47][0] = PDGselectronR; ParticleHiggsH3.Array_Decays[47][1] = -PDGselectronL; ParticleHiggsH3.Array_Decays[47][2] = H03amplitudeselectronRantiselectronL; ParticleHiggsH3.Array_Decays[47][3] = 2; ParticleHiggsH3.Array_Comments[47] = "# H3 -> ~e_R- ~e_L+";
      ParticleHiggsH3.Array_Decays[48][0] = PDGnusmuonL; ParticleHiggsH3.Array_Decays[48][1] = -PDGnusmuonL; ParticleHiggsH3.Array_Decays[48][2] = H03amplitudesnumuLantisnumuL; ParticleHiggsH3.Array_Decays[48][3] = 2; ParticleHiggsH3.Array_Comments[48] = "# H3 -> ~nu_muL ~nu_muL*";
      ParticleHiggsH3.Array_Decays[49][0] = PDGsmuonL; ParticleHiggsH3.Array_Decays[49][1] = -PDGsmuonL; ParticleHiggsH3.Array_Decays[49][2] = H03amplitudesmuonLantismuonL; ParticleHiggsH3.Array_Decays[49][3] = 2; ParticleHiggsH3.Array_Comments[49] = "# H3 -> ~mu_L- ~mu_L+";
      ParticleHiggsH3.Array_Decays[50][0] = PDGsmuonR; ParticleHiggsH3.Array_Decays[50][1] = -PDGsmuonR; ParticleHiggsH3.Array_Decays[50][2] = H03amplitudesmuonRantismuonR; ParticleHiggsH3.Array_Decays[50][3] = 2; ParticleHiggsH3.Array_Comments[50] = "# H3 -> ~mu_R- ~mu_R+";
      ParticleHiggsH3.Array_Decays[51][0] = PDGsmuonL; ParticleHiggsH3.Array_Decays[51][1] = -PDGsmuonR; ParticleHiggsH3.Array_Decays[51][2] = H03amplitudesmuonLantismuonR; ParticleHiggsH3.Array_Decays[51][3] = 2; ParticleHiggsH3.Array_Comments[51] = "# H3 -> ~mu_L- ~mu_R+";
      ParticleHiggsH3.Array_Decays[52][0] = PDGsmuonR; ParticleHiggsH3.Array_Decays[52][1] = -PDGsmuonL; ParticleHiggsH3.Array_Decays[52][2] = H03amplitudesmuonRantismuonL; ParticleHiggsH3.Array_Decays[52][3] = 2; ParticleHiggsH3.Array_Comments[52] = "# H3 -> ~mu_R- ~mu_L+";	
      ParticleHiggsH3.Array_Decays[53][0] = PDGnustauL; ParticleHiggsH3.Array_Decays[53][1] = -PDGnustauL; ParticleHiggsH3.Array_Decays[53][2] = H03amplitudesnutauLantisnutauL; ParticleHiggsH3.Array_Decays[53][3] = 2; ParticleHiggsH3.Array_Comments[53] = "# H3 -> ~nu_tauL ~nu_tauL*";
      ParticleHiggsH3.Array_Decays[54][0] = PDGstop1; ParticleHiggsH3.Array_Decays[54][1] = -PDGstop1; ParticleHiggsH3.Array_Decays[54][2] = H03amplitudestop1antistop1; ParticleHiggsH3.Array_Decays[54][3] = 2; ParticleHiggsH3.Array_Comments[54] = "# H3 -> ~t_1 ~t_1*";	  
      ParticleHiggsH3.Array_Decays[55][0] = PDGstop2; ParticleHiggsH3.Array_Decays[55][1] = -PDGstop2; ParticleHiggsH3.Array_Decays[55][2] = H03amplitudestop2antistop2; ParticleHiggsH3.Array_Decays[55][3] = 2; ParticleHiggsH3.Array_Comments[55] = "# H3 -> ~t_2 ~t_2*";
      ParticleHiggsH3.Array_Decays[56][0] = PDGstop1; ParticleHiggsH3.Array_Decays[56][1] = -PDGstop2; ParticleHiggsH3.Array_Decays[56][2] = H03amplitudestop1antistop2; ParticleHiggsH3.Array_Decays[56][3] = 2; ParticleHiggsH3.Array_Comments[56] = "# H3 -> ~t_1 ~t_2*";	  
      ParticleHiggsH3.Array_Decays[57][0] = PDGstop2; ParticleHiggsH3.Array_Decays[57][1] = -PDGstop1; ParticleHiggsH3.Array_Decays[57][2] = H03amplitudestop2antistop1; ParticleHiggsH3.Array_Decays[57][3] = 2; ParticleHiggsH3.Array_Comments[57] = "# H3 -> ~t_2 ~t_1*";
      ParticleHiggsH3.Array_Decays[58][0] = PDGsbottom1; ParticleHiggsH3.Array_Decays[58][1] = -PDGsbottom1; ParticleHiggsH3.Array_Decays[58][2] = H03amplitudesbottom1antisbottom1; ParticleHiggsH.Array_Decays[58][3] = 2; ParticleHiggsH3.Array_Comments[58] = "# H3 -> ~b_1 ~b_1*";	  
      ParticleHiggsH3.Array_Decays[59][0] = PDGsbottom2; ParticleHiggsH3.Array_Decays[59][1] = -PDGsbottom2; ParticleHiggsH3.Array_Decays[59][2] = H03amplitudesbottom2antisbottom2; ParticleHiggsH3.Array_Decays[59][3] = 2; ParticleHiggsH3.Array_Comments[59] = "# H3 -> ~b_2 ~b_2*";
      ParticleHiggsH3.Array_Decays[60][0] = PDGsbottom1; ParticleHiggsH3.Array_Decays[60][1] = -PDGsbottom2; ParticleHiggsH3.Array_Decays[60][2] = H03amplitudesbottom1antisbottom2; ParticleHiggsH3.Array_Decays[60][3] = 2; ParticleHiggsH3.Array_Comments[60] = "# H3 -> ~b_1 ~b_2*";	  
      ParticleHiggsH3.Array_Decays[61][0] = PDGsbottom2; ParticleHiggsH3.Array_Decays[61][1] = -PDGsbottom1; ParticleHiggsH3.Array_Decays[61][2] = H03amplitudesbottom2antisbottom1; ParticleHiggsH3.Array_Decays[61][3] = 2; ParticleHiggsH3.Array_Comments[61] = "# H3 -> ~b_2 ~b_1*";
      ParticleHiggsH3.Array_Decays[62][0] = PDGstau1; ParticleHiggsH3.Array_Decays[62][1] = -PDGstau1; ParticleHiggsH3.Array_Decays[62][2] = H03amplitudestau1antistau1; ParticleHiggsH3.Array_Decays[62][3] = 2; ParticleHiggsH3.Array_Comments[62] = "# H3 -> ~tau_1- ~tau_1+";	  
      ParticleHiggsH3.Array_Decays[63][0] = PDGstau2; ParticleHiggsH3.Array_Decays[63][1] = -PDGstau2; ParticleHiggsH3.Array_Decays[63][2] = H03amplitudestau2antistau2; ParticleHiggsH3.Array_Decays[63][3] = 2; ParticleHiggsH3.Array_Comments[63] = "# H3 -> ~tau_2- ~tau_2+";
      ParticleHiggsH3.Array_Decays[64][0] = PDGstau1; ParticleHiggsH3.Array_Decays[64][1] = -PDGstau2; ParticleHiggsH3.Array_Decays[64][2] = H03amplitudestau1antistau2; ParticleHiggsH3.Array_Decays[64][3] = 2; ParticleHiggsH3.Array_Comments[64] = "# H3 -> ~tau_1- ~tau_2+";	  
      ParticleHiggsH3.Array_Decays[65][0] = PDGstau2; ParticleHiggsH3.Array_Decays[65][1] = -PDGstau1; ParticleHiggsH3.Array_Decays[65][2] = H03amplitudestau2antistau1; ParticleHiggsH3.Array_Decays[65][3] = 2; ParticleHiggsH3.Array_Comments[65] = "# H3 -> ~tau_2- ~tau_1+";
      
      ParticleHiggsH3.Array_Decays[66][0] = PDGphoton; ParticleHiggsH3.Array_Decays[66][1] = PDGphoton; ParticleHiggsH3.Array_Decays[66][2] = H03amplitudegammagamma; ParticleHiggsH3.Array_Decays[66][3] = 2; ParticleHiggsH3.Array_Comments[66] = "# H3 -> gamma gamma";
      ParticleHiggsH3.Array_Decays[67][0] = PDGgluon; ParticleHiggsH3.Array_Decays[67][1] = PDGgluon; ParticleHiggsH3.Array_Decays[67][2] = H03amplitudegluongluon; ParticleHiggsH3.Array_Decays[67][3] = 2; ParticleHiggsH3.Array_Comments[67] = "# H3 -> gluon gluon";
      ParticleHiggsH3.Array_Decays[68][0] = PDGZboson; ParticleHiggsH3.Array_Decays[68][1] = PDGphoton; ParticleHiggsH3.Array_Decays[68][2] = H03amplitudeZgamma; ParticleHiggsH3.Array_Decays[68][3] = 2; ParticleHiggsH3.Array_Comments[68] = "# H3 -> Z gamma";
      ParticleHiggsH3.Array_Decays[69][0] = PDGWplus; ParticleHiggsH3.Array_Decays[69][1] = -PDGWplus; ParticleHiggsH3.Array_Decays[69][2] = H03amplitudeWW; ParticleHiggsH3.Array_Decays[69][3] = H03WWNDA; ParticleHiggsH3.Array_Comments[69] = H03WWcomment;
      ParticleHiggsH3.Array_Decays[70][0] = PDGZboson; ParticleHiggsH3.Array_Decays[70][1] = -PDGZboson; ParticleHiggsH3.Array_Decays[70][2] = H03amplitudeZZ; ParticleHiggsH3.Array_Decays[70][3] = H03ZZNDA; ParticleHiggsH3.Array_Comments[70] = H03ZZcomment;
      
      ParticleHiggsH3.Array_Decays[71][0] = PDGneutralino1; ParticleHiggsH3.Array_Decays[71][1] = PDGneutralino5; ParticleHiggsH3.Array_Decays[71][2] = H03amplitudeneutZ1neutZ5; ParticleHiggsH3.Array_Decays[71][3] = 2; ParticleHiggsH3.Array_Comments[71] = "# H3 -> ~chi_10 ~chi_50";
      ParticleHiggsH3.Array_Decays[72][0] = PDGneutralino2; ParticleHiggsH3.Array_Decays[72][1] = PDGneutralino5; ParticleHiggsH3.Array_Decays[72][2] = H03amplitudeneutZ2neutZ5; ParticleHiggsH3.Array_Decays[72][3] = 2; ParticleHiggsH3.Array_Comments[72] = "# H3 -> ~chi_20 ~chi_50";
      ParticleHiggsH3.Array_Decays[73][0] = PDGneutralino3; ParticleHiggsH3.Array_Decays[73][1] = PDGneutralino5; ParticleHiggsH3.Array_Decays[73][2] = H03amplitudeneutZ3neutZ5; ParticleHiggsH3.Array_Decays[73][3] = 2; ParticleHiggsH3.Array_Comments[73] = "# H3 -> ~chi_30 ~chi_50";
      ParticleHiggsH3.Array_Decays[74][0] = PDGneutralino4; ParticleHiggsH3.Array_Decays[74][1] = PDGneutralino5; ParticleHiggsH3.Array_Decays[74][2] = H03amplitudeneutZ4neutZ5; ParticleHiggsH3.Array_Decays[74][3] = 2; ParticleHiggsH3.Array_Comments[74] = "# H3 -> ~chi_40 ~chi_50";
      ParticleHiggsH3.Array_Decays[75][0] = PDGneutralino5; ParticleHiggsH3.Array_Decays[75][1] = PDGneutralino5; ParticleHiggsH3.Array_Decays[75][2] = H03amplitudeneutZ5neutZ5; ParticleHiggsH3.Array_Decays[75][3] = 2; ParticleHiggsH3.Array_Comments[75] = "# H3 -> ~chi_50 ~chi_50";
      ParticleHiggsH3.Array_Decays[76][0] = PDGA0; ParticleHiggsH3.Array_Decays[76][1] = PDGA2; ParticleHiggsH3.Array_Decays[76][2] = H03amplitudehiggsAhiggsA2; ParticleHiggsH3.Array_Decays[76][3] = 2; ParticleHiggsH3.Array_Comments[76] = "# H3 -> A A2";

      ParticleHiggsH3.Array_Decays[78][0] = PDGA2; ParticleHiggsH3.Array_Decays[78][1] = PDGA2; ParticleHiggsH3.Array_Decays[78][2] = H03amplitudehiggsA2higgsA2; ParticleHiggsH3.Array_Decays[78][3] = 2; ParticleHiggsH3.Array_Comments[78] = "# H3 -> A2 A2";
      ParticleHiggsH3.Array_Decays[79][0] = PDGZboson; ParticleHiggsH3.Array_Decays[79][1] = PDGA2; ParticleHiggsH3.Array_Decays[79][2] = H03amplitudehiggsA2Zboson; ParticleHiggsH3.Array_Decays[79][3] = 2; ParticleHiggsH3.Array_Comments[79] = "# H3 -> A2 Z";
      ParticleHiggsH3.Array_Decays[80][0] = PDGh0; ParticleHiggsH3.Array_Decays[80][1] = PDGH0; ParticleHiggsH3.Array_Decays[80][2] = H03amplitudeh0H0; ParticleHiggsH3.Array_Decays[80][3] = 2; ParticleHiggsH3.Array_Comments[80] = "# H3 -> h H";
      ParticleHiggsH3.Array_Decays[81][0] = PDGH0; ParticleHiggsH3.Array_Decays[81][1] = PDGH0; ParticleHiggsH3.Array_Decays[81][2] = H03amplitudeH0H0; ParticleHiggsH3.Array_Decays[81][3] = 2; ParticleHiggsH3.Array_Comments[81] = "# H3 -> H H";
      ParticleHiggsH3.Array_Decays[82][0] = PDGHplus; ParticleHiggsH3.Array_Decays[82][1] = -PDGWplus; ParticleHiggsH3.Array_Decays[82][2] = H03amplitudeWHpm; ParticleHiggsH3.Array_Decays[82][3] = 2; ParticleHiggsH3.Array_Comments[82] = "# H3 -> H+- W-+";

      for(int i = 0; i<ParticleHiggsH3.No_of_Decays; i++) {
	if (ParticleHiggsH3.Array_Decays[i][2] < 0) {
	  fout << "#warning! Partial Width for " << ParticleHiggsH3.Array_Comments[i] << " is negative = " << ParticleHiggsH3.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
	  ParticleHiggsH3.Array_Decays[i][2] = 0;
	  errorflag = -1;
	}
   }          
      
      double HiggsH3_No_1to2_Decays = 0;
      
      HiggsH3_No_1to2_Decays = ParticleHiggsH3.No_1to2_Decays + ParticleHiggsH3.No_NMSSM_Decays; 
      
      for (int j = 0; j<HiggsH3_No_1to2_Decays; j++) {
	ParticleHiggsH3.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
      }
      
      for (int j=0; j<HiggsH3_No_1to2_Decays; j++) {
	ParticleHiggsH3.two_width = ParticleHiggsH3.two_width + ParticleHiggsH3.Array_Decays[j][2];
      }
      for (int j=HiggsH3_No_1to2_Decays; j<ParticleHiggsH3.No_of_Decays; j++) {
	ParticleHiggsH3.three_width = ParticleHiggsH3.three_width + ParticleHiggsH3.Array_Decays[j][2];
      }
      
      for(int j=0; j<ParticleHiggsH3.No_of_Decays; j++) {
	ParticleHiggsH3.Array_Decays[j][4] = 0;
      }
      
      ///Could argue no need for test for nans here as the higgs 1 -> 3 decay formulae are all purely analytic algebraic expressions, therefore no numerical integration is involved so we can't get nans. Will check anyway as possibility of -ve sqrts in kinematics or -ve logs, or infs etc
   if ( ParticleHiggsH3.three_width != ParticleHiggsH3.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for HiggsH3 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleHiggsH3.No_of_Decays = HiggsH3_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleHiggsH3.total_width = ParticleHiggsH3.two_width;
     }
   else {
     ParticleHiggsH3.total_width = ParticleHiggsH3.two_width + ParticleHiggsH3.three_width;
   }

   if ( ParticleHiggsH3.total_width != ParticleHiggsH3.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleHiggsH3.No_of_Decays; i++) {
       //   fout << ParticleHiggsH3.Array_Decays[i][2] << endl;
       // }	    
       throw( "nan in H03 heaviest higgs total width \n");
     }
   
   }
 }

 ///higgsA decays

 double A0amplitudeuantiu=0, A0amplitudedantid=0, A0amplitudesantis=0, A0amplitudecantic=0, A0amplitudebantib=0, A0amplitudetantit=0, A0amplitudeeantie=0, A0amplitudemuantimu=0, A0amplitudetauantitau=0, A0amplitudeneutZ1neutZ1=0, A0amplitudeneutZ1neutZ2=0, A0amplitudeneutZ1neutZ3=0, A0amplitudeneutZ1neutZ4=0, A0amplitudeneutZ2neutZ2=0, A0amplitudeneutZ2neutZ3=0, A0amplitudeneutZ2neutZ4=0, A0amplitudeneutZ3neutZ3=0, A0amplitudeneutZ3neutZ4=0, A0amplitudeneutZ4neutZ4=0, A0amplitudecharW1charW1=0, A0amplitudecharW1charW2=0, A0amplitudecharW2charW2=0, A0amplitudehiggshZboson=0, A0amplitudehiggsHZboson=0, A0amplitudesupLantisupR=0, A0amplitudesupRantisupL=0, A0amplitudesdownLantisdownR=0, A0amplitudesdownRantisdownL=0, A0amplitudescharmLantischarmR=0, A0amplitudescharmRantischarmL=0, A0amplitudesstrangeLantisstrangeR=0, A0amplitudesstrangeRantisstrangeL=0, A0amplitudestop1antistop2=0, A0amplitudestop2antistop1=0, A0amplitudesbottom1antisbottom2=0, A0amplitudesbottom2antisbottom1=0, A0amplitudeselectronLselectronR=0, A0amplitudeselectronRselectronL=0, A0amplitudesmuonLsmuonR=0, A0amplitudesmuonRsmuonL=0, A0amplitudestau1stau2=0, A0amplitudestau2stau1=0, A0amplitudegluongluon=0, A0amplitudegammagamma=0, A0amplitudeZgamma=0, A0amplitudeWHpm=0;

 double A0amplitudeneutZ1neutZ5=0, A0amplitudeneutZ2neutZ5=0, A0amplitudeneutZ3neutZ5=0, A0amplitudeneutZ4neutZ5=0, A0amplitudeneutZ5neutZ5=0, A0amplitudehiggsH3Zboson=0;

 if (flagA1 == 1) {
   if (nmssmIsIt == false) {
     if (QCDcorr == false) {
       ///No decays to u or d as PWs to u and d are tiny as proportional to yukawas squared
       ///Use running masses here to try to approximate some of the correction (which aren't included)
       A0amplitudecantic = higgsAamplitudedecayquarkantiquark (mA0(1), runmc, g, beta, runmw, 1, QCDcorr, alphasAtMA);
       A0amplitudesantis = higgsAamplitudedecayquarkantiquark (mA0(1), runms, g, beta, runmw, 0, QCDcorr, alphasAtMA);
       A0amplitudebantib = higgsAamplitudedecayquarkantiquark (mA0(1), runmb, g, beta, runmw, 0, QCDcorr, alphasAtMA);
       A0amplitudetantit = higgsAamplitudedecayquarkantiquark (mA0(1), runmt, g, beta, runmw, 1, QCDcorr, alphasAtMA); 
     }
     else {
       A0amplitudecantic = higgsAamplitudedecayquarkantiquark (mA0(1), mcpole, g, beta, runmw, 1, QCDcorr, alphasAtMA);
       A0amplitudesantis = higgsAamplitudedecayquarkantiquark (mA0(1), mspole, g, beta, runmw, 0, QCDcorr, alphasAtMA);
       ///mcpole and mspole set in decays.h, this values used are those appropriate for the scheme used for the h -> qq QCD corrections, as in hdecay
       A0amplitudebantib = higgsAamplitudedecayquarkantiquark (mA0(1), mbPole, g, beta, runmw, 0, QCDcorr, alphasAtMA);
       A0amplitudetantit = higgsAamplitudedecayquarkantiquark (mA0(1), mtPole, g, beta, runmw, 1, QCDcorr, alphasAtMA);
     }

     A0amplitudeeantie = higgsAamplitudedecayquarkantiquark (mA0(1), runmel, g, beta, runmw, 0, false, alphasAtMA)/3; ///0 as leptons are like down-type quarks, divide by 3 as No of colours is 1 for leptons cf 3 for quarks
     A0amplitudemuantimu = higgsAamplitudedecayquarkantiquark (mA0(1), runmmu, g, beta, runmw, 0, false, alphasAtMA)/3;
     A0amplitudetauantitau = higgsAamplitudedecayquarkantiquark (mA0(1), runmtau, g, beta, runmw, 0, false, alphasAtMA)/3; 
     A0amplitudeneutZ1neutZ1 = higgsphiamplitudedecayneutralinoneutralino (mA0(1), mneut(1), mneut(1), g, tanthetaW, beta, mixNeut, 1, 1, 'A');
     A0amplitudeneutZ1neutZ2 = higgsphiamplitudedecayneutralinoneutralino (mA0(1), mneut(1), mneut(2), g, tanthetaW, beta, mixNeut, 1, 2, 'A');
     A0amplitudeneutZ1neutZ3 = higgsphiamplitudedecayneutralinoneutralino (mA0(1), mneut(1), mneut(3), g, tanthetaW, beta, mixNeut, 1, 3, 'A');
     A0amplitudeneutZ1neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mA0(1), mneut(1), mneut(4), g, tanthetaW, beta, mixNeut, 1, 4, 'A');
     A0amplitudeneutZ2neutZ2 = higgsphiamplitudedecayneutralinoneutralino (mA0(1), mneut(2), mneut(2), g, tanthetaW, beta, mixNeut, 2, 2, 'A');
     A0amplitudeneutZ2neutZ3 = higgsphiamplitudedecayneutralinoneutralino (mA0(1), mneut(2), mneut(3), g, tanthetaW, beta, mixNeut, 2, 3, 'A');
     A0amplitudeneutZ2neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mA0(1), mneut(2), mneut(4), g, tanthetaW, beta, mixNeut, 2, 4, 'A');
     A0amplitudeneutZ3neutZ3 = higgsphiamplitudedecayneutralinoneutralino (mA0(1), mneut(3), mneut(3), g, tanthetaW, beta, mixNeut, 3, 3, 'A');
     A0amplitudeneutZ3neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mA0(1), mneut(3), mneut(4), g, tanthetaW, beta, mixNeut, 3, 4, 'A');
     A0amplitudeneutZ4neutZ4 = higgsphiamplitudedecayneutralinoneutralino (mA0(1), mneut(4), mneut(4), g, tanthetaW, beta, mixNeut, 4, 4, 'A');
     A0amplitudecharW1charW1 = higgsphiamplitudedecaysamechargino (mA0(1), mch(1), g, thetaL2, thetaR2, alpha, beta, 1, 'A');
     A0amplitudecharW2charW2 = higgsphiamplitudedecaysamechargino (mA0(1), mch(2), g, thetaL2, thetaR2, alpha, beta, 2, 'A');
     A0amplitudecharW1charW2 = higgsphiamplitudedecaydifchargino (mA0(1), mch(1), mch(2), g, thetaL2, thetaR2, alpha, beta, 'A');
     A0amplitudehiggshZboson = higgsAamplitudedecayhiggshZboson (mA0(1), polemz, mh0(1), g, gp, alpha, beta);
     
     ///In general may wish to not allow the A -> HZ decay for heavy higgs as it's ruled out by SUSY constraints on the mass spectrum (Djouadi Tome II)?
     A0amplitudehiggsHZboson = higgsAamplitudedecayhiggsHZboson (mA0(1), polemz, mh0(2), g, gp, alpha, beta);
     A0amplitudesupLantisupR = 3*higgsAamplitudedecaysfermions (mA0(1), mu(1,1), mu(2,1), g, runmw, runmu, greekmu, Au, beta, 'u');
     A0amplitudesupRantisupL = 3*higgsAamplitudedecaysfermions (mA0(1), mu(2,1), mu(1,1), g, runmw, runmu, greekmu, Au, beta, 'u');
     A0amplitudesdownLantisdownR = 3*higgsAamplitudedecaysfermions (mA0(1), md(1,1), md(2,1), g, runmw, runmd, greekmu, Ad, beta, 'd');
     A0amplitudesdownRantisdownL = 3*higgsAamplitudedecaysfermions (mA0(1), md(2,1), md(1,1), g, runmw, runmd, greekmu, Ad, beta, 'd');
     A0amplitudescharmLantischarmR = 3*higgsAamplitudedecaysfermions(mA0(1), mu(1,2), mu(2,2), g, runmw, runmc, greekmu, Ac, beta, 'u');
     A0amplitudescharmRantischarmL = 3*higgsAamplitudedecaysfermions(mA0(1), mu(2,2), mu(1,2), g, runmw, runmc, greekmu, Ac, beta, 'u');
     A0amplitudesstrangeLantisstrangeR = 3*higgsAamplitudedecaysfermions(mA0(1), md(1,2), md(2,2), g, runmw, runms, greekmu, As, beta, 'd');
     A0amplitudesstrangeRantisstrangeL = 3*higgsAamplitudedecaysfermions(mA0(1), md(2,2), md(1,2), g, runmw, runms, greekmu, As, beta, 'd');
     A0amplitudestop1antistop2 = 3*higgsAamplitudedecaysfermions(mA0(1), mu(1,3), mu(2,3), g, runmw, runmt, greekmu, At, beta, 'u');
     A0amplitudestop2antistop1 = 3*higgsAamplitudedecaysfermions(mA0(1), mu(2,3), mu(1,3), g, runmw, runmt, greekmu, At, beta, 'u');
     A0amplitudesbottom1antisbottom2 = 3*higgsAamplitudedecaysfermions(mA0(1), md(1,3), md(2,3), g, runmw, runmb, greekmu, Ab, beta, 'd');
     A0amplitudesbottom2antisbottom1 = 3*higgsAamplitudedecaysfermions(mA0(1), md(2,3), md(1,3), g, runmw, runmb, greekmu, Ab, beta, 'd');
     A0amplitudeselectronLselectronR = higgsAamplitudedecaysfermions(mA0(1), me(1,1), me(2,1), g, runmw, runmel, greekmu, Ae, beta, 'd');
     A0amplitudeselectronRselectronL = higgsAamplitudedecaysfermions(mA0(1), me(2,1), me(1,1), g, runmw, runmel, greekmu, Ae, beta, 'd');
     A0amplitudesmuonLsmuonR = higgsAamplitudedecaysfermions(mA0(1), me(1,2), me(2,2), g, runmw, runmmu, greekmu, Amu, beta, 'd');
     A0amplitudesmuonRsmuonL = higgsAamplitudedecaysfermions(mA0(1), me(2,2), me(1,2), g, runmw, runmmu, greekmu, Amu, beta, 'd');
     A0amplitudestau1stau2 = higgsAamplitudedecaysfermions(mA0(1), me(1,3), me(2,3), g, runmw, runmtau, greekmu, Atau, beta, 'd');
     A0amplitudestau2stau1 = higgsAamplitudedecaysfermions(mA0(1), me(2,3), me(1,3), g, runmw, runmtau, greekmu, Atau, beta, 'd');
     
     A0amplitudegammagamma = higgsesamplitudedecaygammagammatotal(mA0(1), g, gp, alphaAtMA, runmw, polemw, alpha, beta, mtAtMA, mbAtMA, mcAtMA, runmtau, mHpm, mu(1,3), mu(2,3), md(1,3), md(2,3), me(1,3), me(2,3), mch(1), mch(2), thetaL, thetaR, thetat, thetab, thetatau, greekmu, At, Ab, Atau, 'A'); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     A0amplitudegluongluon = higgsesamplitudedecaygluongluontotal(mA0(1), g, g3atmA0, gp, runmw, alpha, beta, mtPole, mbPole, mcpole, mu(1,3), mu(2,3), md(1,3), md(2,3), thetat, thetab, greekmu, At, Ab, runms, mu(1,2), mu(2,2), md(1,2), md(2,2), Ac, As, runmu, runmd, mu(1,1), mu(2,1), md(1,1), md(2,1), Au, Ad, 'A', QCDcorr);///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     A0amplitudeZgamma = higgsesamplitudedecayZbosonphotontotal(mA0(1), polemz, g, gp, alphaAtMA, polemw, runmw, alpha, beta, mtAtMA, mbAtMA, mcAtMA, msAtMA, mu(1,3), mu(2,3), md(1,3), md(2,3), mHpm, thetat, thetab, greekmu, At, Ab, 'A');///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
   
     A0amplitudeWHpm = higgsAamplitudedecayHpmWboson(mA0(1), polemw, mHpm, g, thetaA, 1, nmssmIsIt);
   }
   
   else if (nmssmIsIt == true){ ///NMSSM so need to modify by elements of pseudoscalar mixing matrix
     if (QCDcorr == false) {
       ///No decays to u or d as PWs to u and d are tiny as proportional to yukawas squared
       ///Use running masses here to try to approximate some of the correction (which aren't included)
       A0amplitudecantic = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(1), runmc, beta, CPOMix, 1, 1, QCDcorr, alphasAtMA);
       A0amplitudesantis = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(1), runms, beta, CPOMix, 0, 1, QCDcorr, alphasAtMA);
       A0amplitudebantib = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(1), runmb, beta, CPOMix, 0, 1, QCDcorr, alphasAtMA);
       A0amplitudetantit = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(1), runmt, beta, CPOMix, 1, 1, QCDcorr, alphasAtMA);
     }
     else {
       A0amplitudecantic = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(1), mcpole, beta, CPOMix, 1, 1, QCDcorr, alphasAtMA);
       A0amplitudesantis = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(1), mspole, beta, CPOMix, 0, 1, QCDcorr, alphasAtMA);
       ///mcpole and mspole set in decays.h, this values used are those appropriate for the scheme used for the h -> qq QCD corrections, as in hdecay
       A0amplitudebantib = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(1), mbPole, beta, CPOMix, 0, 1, QCDcorr, alphasAtMA);
       A0amplitudetantit = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(1), mtPole, beta, CPOMix, 1, 1, QCDcorr, alphasAtMA);
     }
     
     A0amplitudeeantie = higgsAamplitudedecayquarkantiquarkNMSSM(mA0(1), runmel, beta, CPOMix, 0, 1, false, alphasAtMA)/3; ///0 as leptons are like down-type quarks, divide by 3 as No of colours is 1 for leptons cf 3 for quarks
     A0amplitudemuantimu = higgsAamplitudedecayquarkantiquarkNMSSM(mA0(1), runmmu, beta, CPOMix, 0, 1, false, alphasAtMA)/3;
     A0amplitudetauantitau = higgsAamplitudedecayquarkantiquarkNMSSM(mA0(1), runmtau, beta, CPOMix, 0, 1, false, alphasAtMA)/3;
     
     A0amplitudeneutZ1neutZ1 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(1), mneut(1), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 1, 1, 1);
     A0amplitudeneutZ1neutZ2 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(1), mneut(2), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 1, 2, 1);
     A0amplitudeneutZ1neutZ3 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(1), mneut(3), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 1, 3, 1);
     A0amplitudeneutZ1neutZ4 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(1), mneut(4), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 1, 4, 1);
     A0amplitudeneutZ1neutZ5 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(1), mneut(5), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 1, 5, 1);
     A0amplitudeneutZ2neutZ2 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(2), mneut(2), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 2, 2, 1);
     A0amplitudeneutZ2neutZ3 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(2), mneut(3), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 2, 3, 1);
     A0amplitudeneutZ2neutZ4 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(2), mneut(4), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 2, 4, 1);
     A0amplitudeneutZ2neutZ5 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(2), mneut(5), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 2, 5, 1);
     A0amplitudeneutZ3neutZ3 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(3), mneut(3), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 3, 3, 1);
     A0amplitudeneutZ3neutZ4 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(3), mneut(4), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 3, 4, 1);
     A0amplitudeneutZ3neutZ5 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(3), mneut(5), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 3, 5, 1);
     A0amplitudeneutZ4neutZ4 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(4), mneut(4), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 4, 4, 1);
     A0amplitudeneutZ4neutZ5 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(4), mneut(5), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 4, 5, 1);
     A0amplitudeneutZ5neutZ5 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(1), mneut(5), mneut(5), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 5, 5, 1);
     A0amplitudecharW1charW1 = higgsAamplitudedecaysamecharginoNMSSM (mA0(1), mch(1), g, thetaL2, thetaR2, alpha, lam, CPOMix, 1, 1);
     A0amplitudecharW2charW2 = higgsAamplitudedecaysamecharginoNMSSM (mA0(1), mch(2), g, thetaL2, thetaR2, alpha, lam, CPOMix, 2, 1);
     A0amplitudecharW1charW2 = higgsAamplitudedecaydifcharginoNMSSM (mA0(1), mch(1), mch(2), g, thetaL2, thetaR2, alpha, lam, CPOMix, 1);
     A0amplitudehiggshZboson = higgsAamplitudedecayhiggshorHZbosonNMSSM (mA0(1), polemz, mh0(1), g, gp, alpha, beta, thetaA, CPEMix, 1, 1);
     A0amplitudehiggsHZboson = higgsAamplitudedecayhiggshorHZbosonNMSSM (mA0(1), polemz, mh0(2), g, gp, alpha, beta, thetaA, CPEMix, 1, 2);
     A0amplitudehiggsH3Zboson = higgsAamplitudedecayhiggshorHZbosonNMSSM (mA0(1), polemz, mh0(3), g, gp, alpha, beta, thetaA, CPEMix, 1, 3);   
     A0amplitudesupLantisupR = 3*higgsAamplitudedecaysfermionsNMSSM (mA0(1), mu(1,1), mu(2,1), g, runmw, runmu, Au, beta, lam, mueff, CPOMix, 'u', 1);
     A0amplitudesupRantisupL = 3*higgsAamplitudedecaysfermionsNMSSM (mA0(1), mu(2,1), mu(1,1), g, runmw, runmu, Au, beta, lam, mueff, CPOMix, 'u', 1);
     A0amplitudesdownLantisdownR = 3*higgsAamplitudedecaysfermionsNMSSM (mA0(1), md(1,1), md(2,1), g, runmw, runmd, Ad, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudesdownRantisdownL = 3*higgsAamplitudedecaysfermionsNMSSM (mA0(1), md(2,1), md(1,1), g, runmw, runmd, Ad, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudescharmLantischarmR = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(1), mu(1,2), mu(2,2), g, runmw, runmc, Ac, beta, lam, mueff, CPOMix, 'u', 1);
     A0amplitudescharmRantischarmL = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(1), mu(2,2), mu(1,2), g, runmw, runmc, Ac, beta, lam, mueff, CPOMix, 'u', 1);
     A0amplitudesstrangeLantisstrangeR = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(1), md(1,2), md(2,2), g, runmw, runms, As, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudesstrangeRantisstrangeL = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(1), md(2,2), md(1,2), g, runmw, runms, As, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudestop1antistop2 = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(1), mu(1,3), mu(2,3), g, runmw, runmt, At, beta, lam, mueff, CPOMix, 'u', 1);
     A0amplitudestop2antistop1 = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(1), mu(2,3), mu(1,3), g, runmw, runmt, At, beta, lam, mueff, CPOMix, 'u', 1);
     A0amplitudesbottom1antisbottom2 = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(1), md(1,3), md(2,3), g, runmw, runmb, Ab, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudesbottom2antisbottom1 = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(1), md(2,3), md(1,3), g, runmw, runmb, Ab, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudeselectronLselectronR = higgsAamplitudedecaysfermionsNMSSM(mA0(1), me(1,1), me(2,1), g, runmw, runmel, Ae, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudeselectronRselectronL = higgsAamplitudedecaysfermionsNMSSM(mA0(1), me(2,1), me(1,1), g, runmw, runmel, Ae, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudesmuonLsmuonR = higgsAamplitudedecaysfermionsNMSSM(mA0(1), me(1,2), me(2,2), g, runmw, runmmu, Amu, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudesmuonRsmuonL = higgsAamplitudedecaysfermionsNMSSM(mA0(1), me(2,2), me(1,2), g, runmw, runmmu, Amu, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudestau1stau2 = higgsAamplitudedecaysfermionsNMSSM(mA0(1), me(1,3), me(2,3), g, runmw, runmtau, Atau, beta, lam, mueff, CPOMix, 'd', 1);
     A0amplitudestau2stau1 = higgsAamplitudedecaysfermionsNMSSM(mA0(1), me(2,3), me(1,3), g, runmw, runmtau, Atau, beta, lam, mueff, CPOMix, 'd', 1); 

   A0amplitudegammagamma = higgsAamplitudedecaygammagammaNMSSM (mA0(1), g, gp, alphaAtMA, runmw, CPOMix, beta, mtAtMA, mbAtMA, mcAtMA, runmtau, mch(1), mch(2), thetaL2, thetaR2, lam, 1); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
   A0amplitudegluongluon = higgsAamplitudedecaygluongluonNMSSM (mA0(1), g, gs, alphasAtMA, runmw, CPOMix, beta, mtPole, mbPole, mcpole, lam, 1, QCDcorr); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
   A0amplitudeZgamma = higgsAamplitudedecayZgammaNMSSM (mA0(1), g, gp, alphaAtMA, runmw, polemz, CPOMix, beta,mtAtMA, mbAtMA, mcAtMA, mch(1), mch(2), thetaL2, thetaR2, lam, 1); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
   
   A0amplitudeWHpm = higgsAamplitudedecayHpmWboson(mA0(1), polemw, mHpm, g, thetaA, 1, nmssmIsIt);
   }
   
   ParticleHiggsA.Array_Decays[0][0] = PDGup; ParticleHiggsA.Array_Decays[0][1] = -PDGup; ParticleHiggsA.Array_Decays[0][2] = A0amplitudeuantiu; ParticleHiggsA.Array_Decays[0][3] = 2; ParticleHiggsA.Array_Comments[0] = "# A -> u ub";
   ParticleHiggsA.Array_Decays[1][0] = PDGdown; ParticleHiggsA.Array_Decays[1][1] = -PDGdown; ParticleHiggsA.Array_Decays[1][2] = A0amplitudedantid; ParticleHiggsA.Array_Decays[1][3] = 2; ParticleHiggsA.Array_Comments[1] = "# A -> d db";
   ParticleHiggsA.Array_Decays[2][0] = PDGcharm; ParticleHiggsA.Array_Decays[2][1] = -PDGcharm; ParticleHiggsA.Array_Decays[2][2] = A0amplitudecantic; ParticleHiggsA.Array_Decays[2][3] = 2; ParticleHiggsA.Array_Comments[2] = "# A -> c cb";
   ParticleHiggsA.Array_Decays[3][0] = PDGstrange; ParticleHiggsA.Array_Decays[3][1] = -PDGstrange; ParticleHiggsA.Array_Decays[3][2] = A0amplitudesantis; ParticleHiggsA.Array_Decays[3][3] = 2; ParticleHiggsA.Array_Comments[3] = "# A -> s sb";
   ParticleHiggsA.Array_Decays[4][0] = PDGbottom; ParticleHiggsA.Array_Decays[4][1] = -PDGbottom; ParticleHiggsA.Array_Decays[4][2] = A0amplitudebantib; ParticleHiggsA.Array_Decays[4][3] = 2; ParticleHiggsA.Array_Comments[4] = "# A -> b bb";
   ParticleHiggsA.Array_Decays[5][0] = PDGtop; ParticleHiggsA.Array_Decays[5][1] = -PDGtop; ParticleHiggsA.Array_Decays[5][2] = A0amplitudetantit; ParticleHiggsA.Array_Decays[5][3] = 2; ParticleHiggsA.Array_Comments[5] = "# A -> t tb";
   ParticleHiggsA.Array_Decays[6][0] = PDGelectron; ParticleHiggsA.Array_Decays[6][1] = -PDGelectron; ParticleHiggsA.Array_Decays[6][2] = A0amplitudeeantie; ParticleHiggsA.Array_Decays[6][3] = 2; ParticleHiggsA.Array_Comments[6] = "# A -> e- e+";
   ParticleHiggsA.Array_Decays[7][0] = PDGmuon; ParticleHiggsA.Array_Decays[7][1] = -PDGmuon; ParticleHiggsA.Array_Decays[7][2] = A0amplitudemuantimu; ParticleHiggsA.Array_Decays[7][3] = 2; ParticleHiggsA.Array_Comments[7] = "# A -> mu- mu+";
   ParticleHiggsA.Array_Decays[8][0] = PDGtau; ParticleHiggsA.Array_Decays[8][1] = -PDGtau; ParticleHiggsA.Array_Decays[8][2] = A0amplitudetauantitau; ParticleHiggsA.Array_Decays[8][3] = 2; ParticleHiggsA.Array_Comments[8] = "# A -> tau- tau+";
   ParticleHiggsA.Array_Decays[9][0] = PDGneutralino1; ParticleHiggsA.Array_Decays[9][1] = PDGneutralino1; ParticleHiggsA.Array_Decays[9][2] = A0amplitudeneutZ1neutZ1; ParticleHiggsA.Array_Decays[9][3] = 2; ParticleHiggsA.Array_Comments[9] = "# A -> ~chi_10 ~chi_10";
   ParticleHiggsA.Array_Decays[10][0] = PDGneutralino1; ParticleHiggsA.Array_Decays[10][1] = PDGneutralino2; ParticleHiggsA.Array_Decays[10][2] = A0amplitudeneutZ1neutZ2; ParticleHiggsA.Array_Decays[10][3] = 2; ParticleHiggsA.Array_Comments[10] = "# A -> ~chi_10 ~chi_20";
   ParticleHiggsA.Array_Decays[11][0] = PDGneutralino1; ParticleHiggsA.Array_Decays[11][1] = PDGneutralino3; ParticleHiggsA.Array_Decays[11][2] = A0amplitudeneutZ1neutZ3; ParticleHiggsA.Array_Decays[11][3] = 2; ParticleHiggsA.Array_Comments[11] = "# A -> ~chi_10 ~chi_30";
   ParticleHiggsA.Array_Decays[12][0] = PDGneutralino1; ParticleHiggsA.Array_Decays[12][1] = PDGneutralino4; ParticleHiggsA.Array_Decays[12][2] = A0amplitudeneutZ1neutZ4; ParticleHiggsA.Array_Decays[12][3] = 2; ParticleHiggsA.Array_Comments[12] = "# A -> ~chi_10 ~chi_40";
   ParticleHiggsA.Array_Decays[13][0] = PDGneutralino2; ParticleHiggsA.Array_Decays[13][1] = PDGneutralino2; ParticleHiggsA.Array_Decays[13][2] = A0amplitudeneutZ2neutZ2; ParticleHiggsA.Array_Decays[13][3] = 2; ParticleHiggsA.Array_Comments[13] = "# A -> ~chi_20 ~chi_20";
   ParticleHiggsA.Array_Decays[14][0] = PDGneutralino2; ParticleHiggsA.Array_Decays[14][1] = PDGneutralino3; ParticleHiggsA.Array_Decays[14][2] = A0amplitudeneutZ2neutZ3; ParticleHiggsA.Array_Decays[14][3] = 2; ParticleHiggsA.Array_Comments[14] = "# A -> ~chi_20 ~chi_30";
   ParticleHiggsA.Array_Decays[15][0] = PDGneutralino2; ParticleHiggsA.Array_Decays[15][1] = PDGneutralino4; ParticleHiggsA.Array_Decays[15][2] = A0amplitudeneutZ2neutZ4; ParticleHiggsA.Array_Decays[15][3] = 2; ParticleHiggsA.Array_Comments[15] = "# A -> ~chi_20 ~chi_40";
   ParticleHiggsA.Array_Decays[16][0] = PDGneutralino3; ParticleHiggsA.Array_Decays[16][1] = PDGneutralino3; ParticleHiggsA.Array_Decays[16][2] = A0amplitudeneutZ3neutZ3; ParticleHiggsA.Array_Decays[16][3] = 2; ParticleHiggsA.Array_Comments[16] = "# A -> ~chi_30 ~chi_30";
   ParticleHiggsA.Array_Decays[17][0] = PDGneutralino3; ParticleHiggsA.Array_Decays[17][1] = PDGneutralino4; ParticleHiggsA.Array_Decays[17][2] = A0amplitudeneutZ3neutZ4; ParticleHiggsA.Array_Decays[17][3] = 2; ParticleHiggsA.Array_Comments[17] = "# A -> ~chi_30 ~chi_40";
   ParticleHiggsA.Array_Decays[18][0] = PDGneutralino4; ParticleHiggsA.Array_Decays[18][1] = PDGneutralino4; ParticleHiggsA.Array_Decays[18][2] = A0amplitudeneutZ4neutZ4; ParticleHiggsA.Array_Decays[18][3] = 2; ParticleHiggsA.Array_Comments[18] = "# A -> ~chi_40 ~chi_40";
   ParticleHiggsA.Array_Decays[19][0] = PDGchargino1; ParticleHiggsA.Array_Decays[19][1] = -PDGchargino1; ParticleHiggsA.Array_Decays[19][2] = A0amplitudecharW1charW1; ParticleHiggsA.Array_Decays[19][3] = 2; ParticleHiggsA.Array_Comments[19] = "# A -> ~chi_1+ ~chi_1-";
   ParticleHiggsA.Array_Decays[20][0] = PDGchargino2; ParticleHiggsA.Array_Decays[20][1] = -PDGchargino2; ParticleHiggsA.Array_Decays[20][2] = A0amplitudecharW2charW2; ParticleHiggsA.Array_Decays[20][3] = 2; ParticleHiggsA.Array_Comments[20] = "# A -> ~chi_2+ ~chi_2-";
   ParticleHiggsA.Array_Decays[21][0] = PDGchargino1; ParticleHiggsA.Array_Decays[21][1] = -PDGchargino2; ParticleHiggsA.Array_Decays[21][2] = A0amplitudecharW1charW2; ParticleHiggsA.Array_Decays[21][3] = 2; ParticleHiggsA.Array_Comments[21] = "# A -> ~chi_1+ ~chi_2-";
   ParticleHiggsA.Array_Decays[22][0] = PDGchargino2; ParticleHiggsA.Array_Decays[22][1] = -PDGchargino1; ParticleHiggsA.Array_Decays[22][2] = A0amplitudecharW1charW2; ParticleHiggsA.Array_Decays[22][3] = 2; ParticleHiggsA.Array_Comments[22] = "# A -> ~chi_2+ ~chi_1-"; ///amplitude same as decay to W1+ and W2- by CP invariance
   ParticleHiggsA.Array_Decays[23][0] = PDGZboson; ParticleHiggsA.Array_Decays[23][1] = PDGh0; ParticleHiggsA.Array_Decays[23][2] = A0amplitudehiggshZboson; ParticleHiggsA.Array_Decays[23][3] = 2; ParticleHiggsA.Array_Comments[23] = "# A -> h Z";
   ParticleHiggsA.Array_Decays[24][0] = PDGZboson; ParticleHiggsA.Array_Decays[24][1] = PDGH0; ParticleHiggsA.Array_Decays[24][2] = A0amplitudehiggsHZboson; ParticleHiggsA.Array_Decays[24][3] = 2; ParticleHiggsA.Array_Comments[24] = "# A -> H Z";
   ParticleHiggsA.Array_Decays[25][0] = PDGsupL; ParticleHiggsA.Array_Decays[25][1] = PDGsupR; ParticleHiggsA.Array_Decays[25][2] = A0amplitudesupLantisupR; ParticleHiggsA.Array_Decays[25][3] = 2; ParticleHiggsA.Array_Comments[25] = "# A-> ~u_L ~u_R*"; 
   ParticleHiggsA.Array_Decays[26][0] = PDGsupR; ParticleHiggsA.Array_Decays[26][1] = PDGsupL; ParticleHiggsA.Array_Decays[26][2] = A0amplitudesupRantisupL; ParticleHiggsA.Array_Decays[26][3] = 2; ParticleHiggsA.Array_Comments[26] = "# A-> ~u_R ~u_L*";
   ParticleHiggsA.Array_Decays[27][0] = PDGsdownL; ParticleHiggsA.Array_Decays[27][1] = PDGsdownR; ParticleHiggsA.Array_Decays[27][2] = A0amplitudesdownLantisdownR; ParticleHiggsA.Array_Decays[27][3] = 2; ParticleHiggsA.Array_Comments[27] = "# A-> ~d_L ~d_R*";
   ParticleHiggsA.Array_Decays[28][0] = PDGsdownR; ParticleHiggsA.Array_Decays[28][1] = PDGsdownL; ParticleHiggsA.Array_Decays[28][2] = A0amplitudesdownRantisdownL; ParticleHiggsA.Array_Decays[28][3] = 2; ParticleHiggsA.Array_Comments[28] = "# A-> ~d_R ~d_L*";
   ParticleHiggsA.Array_Decays[29][0] = PDGscharmL; ParticleHiggsA.Array_Decays[29][1] = PDGscharmR; ParticleHiggsA.Array_Decays[29][2] = A0amplitudescharmLantischarmR; ParticleHiggsA.Array_Decays[29][3] = 2; ParticleHiggsA.Array_Comments[29] = "# A-> ~c_L ~c_R*"; 
   ParticleHiggsA.Array_Decays[30][0] = PDGscharmR; ParticleHiggsA.Array_Decays[30][1] = PDGscharmL; ParticleHiggsA.Array_Decays[30][2] = A0amplitudescharmRantischarmL; ParticleHiggsA.Array_Decays[30][3] = 2; ParticleHiggsA.Array_Comments[30] = "# A-> ~c_R ~c_L*";
   ParticleHiggsA.Array_Decays[31][0] = PDGsstrangeL; ParticleHiggsA.Array_Decays[31][1] = PDGsstrangeR; ParticleHiggsA.Array_Decays[31][2] = A0amplitudesstrangeLantisstrangeR; ParticleHiggsA.Array_Decays[31][3] = 2; ParticleHiggsA.Array_Comments[31] = "# A-> ~s_L ~s_R*";
   ParticleHiggsA.Array_Decays[32][0] = PDGsstrangeR; ParticleHiggsA.Array_Decays[32][1] = PDGsstrangeL; ParticleHiggsA.Array_Decays[32][2] = A0amplitudesstrangeRantisstrangeL; ParticleHiggsA.Array_Decays[32][3] = 2; ParticleHiggsA.Array_Comments[32] = "# A-> ~s_R ~s_L*";
   ParticleHiggsA.Array_Decays[33][0] = PDGstop1; ParticleHiggsA.Array_Decays[33][1] = PDGstop2; ParticleHiggsA.Array_Decays[33][2] = A0amplitudestop1antistop2; ParticleHiggsA.Array_Decays[33][3] = 2; ParticleHiggsA.Array_Comments[33] = "# A-> ~t_1 ~t_2*"; 
   ParticleHiggsA.Array_Decays[34][0] = PDGstop2; ParticleHiggsA.Array_Decays[34][1] = PDGstop1; ParticleHiggsA.Array_Decays[34][2] = A0amplitudestop2antistop1; ParticleHiggsA.Array_Decays[34][3] = 2; ParticleHiggsA.Array_Comments[34] = "# A-> ~t_2 ~t_1*";
   ParticleHiggsA.Array_Decays[35][0] = PDGsbottom1; ParticleHiggsA.Array_Decays[35][1] = PDGsbottom2; ParticleHiggsA.Array_Decays[35][2] = A0amplitudesbottom1antisbottom2; ParticleHiggsA.Array_Decays[35][3] = 2; ParticleHiggsA.Array_Comments[35] = "# A-> ~b_1 ~b_2*";
   ParticleHiggsA.Array_Decays[36][0] = PDGsbottom2; ParticleHiggsA.Array_Decays[36][1] = PDGsbottom1; ParticleHiggsA.Array_Decays[36][2] = A0amplitudesbottom2antisbottom1; ParticleHiggsA.Array_Decays[36][3] = 2; ParticleHiggsA.Array_Comments[36] = "# A-> ~b_2 ~b_1*";
   ParticleHiggsA.Array_Decays[37][0] = PDGselectronL; ParticleHiggsA.Array_Decays[37][1] = PDGselectronR; ParticleHiggsA.Array_Decays[37][2] = A0amplitudeselectronLselectronR; ParticleHiggsA.Array_Decays[37][3] = 2; ParticleHiggsA.Array_Comments[37] = "# A-> ~e_L- ~e_R+"; 
   ParticleHiggsA.Array_Decays[38][0] = PDGselectronR; ParticleHiggsA.Array_Decays[38][1] = PDGselectronL; ParticleHiggsA.Array_Decays[38][2] = A0amplitudeselectronRselectronL; ParticleHiggsA.Array_Decays[38][3] = 2; ParticleHiggsA.Array_Comments[38] = "# A-> ~e_R- ~e_L+";
   ParticleHiggsA.Array_Decays[39][0] = PDGsmuonL; ParticleHiggsA.Array_Decays[39][1] = PDGsmuonR; ParticleHiggsA.Array_Decays[39][2] = A0amplitudesmuonLsmuonR; ParticleHiggsA.Array_Decays[39][3] = 2; ParticleHiggsA.Array_Comments[39] = "# A-> ~mu_L- ~mu_R+";
   ParticleHiggsA.Array_Decays[40][0] = PDGsmuonR; ParticleHiggsA.Array_Decays[40][1] = PDGsmuonL; ParticleHiggsA.Array_Decays[40][2] = A0amplitudesmuonRsmuonL; ParticleHiggsA.Array_Decays[40][3] = 2; ParticleHiggsA.Array_Comments[40] = "# A-> ~mu_R- ~mu_L+";
   ParticleHiggsA.Array_Decays[41][0] = PDGstau1; ParticleHiggsA.Array_Decays[41][1] = PDGstau2; ParticleHiggsA.Array_Decays[41][2] = A0amplitudestau1stau2; ParticleHiggsA.Array_Decays[41][3] = 2; ParticleHiggsA.Array_Comments[41] = "# A-> ~tau_1- ~tau_2+";
   ParticleHiggsA.Array_Decays[42][0] = PDGstau2; ParticleHiggsA.Array_Decays[42][1] = PDGstau1; ParticleHiggsA.Array_Decays[42][2] = A0amplitudestau2stau1; ParticleHiggsA.Array_Decays[42][3] = 2; ParticleHiggsA.Array_Comments[42] = "# A-> ~tau_2- ~tau_1+";
   ParticleHiggsA.Array_Decays[43][0] = PDGphoton; ParticleHiggsA.Array_Decays[43][1] = PDGphoton; ParticleHiggsA.Array_Decays[43][2] = A0amplitudegammagamma; ParticleHiggsA.Array_Decays[43][3] = 2; ParticleHiggsA.Array_Comments[43] = "# A-> gamma gamma";
   ParticleHiggsA.Array_Decays[44][0] = PDGgluon; ParticleHiggsA.Array_Decays[44][1] = PDGgluon; ParticleHiggsA.Array_Decays[44][2] = A0amplitudegluongluon; ParticleHiggsA.Array_Decays[44][3] = 2; ParticleHiggsA.Array_Comments[44] = "# A -> gluon gluon";
   ParticleHiggsA.Array_Decays[45][0] = PDGZboson; ParticleHiggsA.Array_Decays[45][1] = PDGphoton; ParticleHiggsA.Array_Decays[45][2] = A0amplitudeZgamma; ParticleHiggsA.Array_Decays[45][3] = 2; ParticleHiggsA.Array_Comments[45] = "# A -> Z gamma"; 
   ParticleHiggsA.Array_Decays[46][0] = PDGWplus; ParticleHiggsA.Array_Decays[46][1] = -PDGHplus; ParticleHiggsA.Array_Decays[46][2] = A0amplitudeWHpm*2; ParticleHiggsA.Array_Decays[46][3] = 2; ParticleHiggsA.Array_Comments[46] = "# A -> W+- H-+"; ///*2 so includes A -> W- H+
   
   ParticleHiggsA.Array_Decays[47][0] = PDGneutralino1; ParticleHiggsA.Array_Decays[47][1] = PDGneutralino5; ParticleHiggsA.Array_Decays[47][2] = A0amplitudeneutZ1neutZ5; ParticleHiggsA.Array_Decays[47][3] = 2; ParticleHiggsA.Array_Comments[47] = "# A -> ~chi_10 ~chi_50";
   ParticleHiggsA.Array_Decays[48][0] = PDGneutralino1; ParticleHiggsA.Array_Decays[48][1] = PDGneutralino5; ParticleHiggsA.Array_Decays[48][2] = A0amplitudeneutZ2neutZ5; ParticleHiggsA.Array_Decays[48][3] = 2; ParticleHiggsA.Array_Comments[48] = "# A -> ~chi_20 ~chi_50";
   ParticleHiggsA.Array_Decays[49][0] = PDGneutralino1; ParticleHiggsA.Array_Decays[49][1] = PDGneutralino5; ParticleHiggsA.Array_Decays[49][2] = A0amplitudeneutZ3neutZ5; ParticleHiggsA.Array_Decays[49][3] = 2; ParticleHiggsA.Array_Comments[49] = "# A -> ~chi_30 ~chi_50";
   ParticleHiggsA.Array_Decays[50][0] = PDGneutralino1; ParticleHiggsA.Array_Decays[50][1] = PDGneutralino5; ParticleHiggsA.Array_Decays[50][2] = A0amplitudeneutZ4neutZ5; ParticleHiggsA.Array_Decays[50][3] = 2; ParticleHiggsA.Array_Comments[50] = "# A -> ~chi_40 ~chi_50";
   ParticleHiggsA.Array_Decays[51][0] = PDGneutralino1; ParticleHiggsA.Array_Decays[51][1] = PDGneutralino5; ParticleHiggsA.Array_Decays[51][2] = A0amplitudeneutZ5neutZ5; ParticleHiggsA.Array_Decays[51][3] = 2; ParticleHiggsA.Array_Comments[51] = "# A -> ~chi_50 ~chi_50";
   ParticleHiggsA.Array_Decays[52][0] = PDGZboson; ParticleHiggsA.Array_Decays[52][1] = PDGH3; ParticleHiggsA.Array_Decays[52][2] = A0amplitudehiggsH3Zboson; ParticleHiggsA.Array_Decays[52][3] = 2; ParticleHiggsA.Array_Comments[52] = "# A -> H3 Z";

   for(int i = 0; i<ParticleHiggsA.No_of_Decays; i++) {
     if (ParticleHiggsA.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleHiggsA.Array_Comments[i] << " is negative = " << ParticleHiggsA.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleHiggsA.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }       
   
   double HiggsA_No_1to2_Decays = 0;
   
   HiggsA_No_1to2_Decays = ParticleHiggsA.No_1to2_Decays + ParticleHiggsA.No_NMSSM_Decays; /// As higgsA can't be NLSP as heavier than higgsl
   
   for (int j = 0; j<HiggsA_No_1to2_Decays; j++) {
     ParticleHiggsA.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<HiggsA_No_1to2_Decays; j++) {
     ParticleHiggsA.two_width = ParticleHiggsA.two_width + ParticleHiggsA.Array_Decays[j][2];
   }
   for (int j=HiggsA_No_1to2_Decays; j<ParticleHiggsA.No_of_Decays; j++) {
     ParticleHiggsA.three_width = ParticleHiggsA.three_width + ParticleHiggsA.Array_Decays[j][2];
   }
   
   for(int j=0; j<ParticleHiggsA.No_of_Decays; j++) {
     ParticleHiggsA.Array_Decays[j][4] = 0;
   }

   ///Note no 3 body decays for HiggsA
   if ( ParticleHiggsA.three_width != ParticleHiggsA.three_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       fout << "# Three body decays give nan for HiggsA - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
       errorflag = -1;
       ParticleHiggsA.No_of_Decays = HiggsA_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
       ParticleHiggsA.total_width = ParticleHiggsA.two_width;
     }
   else {
     ParticleHiggsA.total_width = ParticleHiggsA.two_width + ParticleHiggsA.three_width;
   }

   if ( ParticleHiggsA.total_width != ParticleHiggsA.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleHiggsA.No_of_Decays; i++) {
       //   fout << ParticleHiggsA.Array_Decays[i][2] << endl;
       // }	    
       throw( "nan in A0 lightest pseudoscalar higgs total width \n");
     }
   
 }

 ///NMSSM Pseudoscalar2 (A2) decays:
 if(nmssmIsIt == true) { ///NMSSM
   double A02amplitudeuantiu=0, A02amplitudedantid=0, A02amplitudesantis=0, A02amplitudecantic=0, A02amplitudebantib=0, A02amplitudetantit=0, A02amplitudeeantie=0, A02amplitudemuantimu=0, A02amplitudetauantitau=0, A02amplitudeneutZ1neutZ1=0, A02amplitudeneutZ1neutZ2=0, A02amplitudeneutZ1neutZ3=0, A02amplitudeneutZ1neutZ4=0, A02amplitudeneutZ2neutZ2=0, A02amplitudeneutZ2neutZ3=0, A02amplitudeneutZ2neutZ4=0, A02amplitudeneutZ3neutZ3=0, A02amplitudeneutZ3neutZ4=0, A02amplitudeneutZ4neutZ4=0, A02amplitudecharW1charW1=0, A02amplitudecharW1charW2=0, A02amplitudecharW2charW2=0, A02amplitudehiggshZboson=0, A02amplitudehiggsHZboson=0, A02amplitudesupLantisupR=0, A02amplitudesupRantisupL=0, A02amplitudesdownLantisdownR=0, A02amplitudesdownRantisdownL=0, A02amplitudescharmLantischarmR=0, A02amplitudescharmRantischarmL=0, A02amplitudesstrangeLantisstrangeR=0, A02amplitudesstrangeRantisstrangeL=0, A02amplitudestop1antistop2=0, A02amplitudestop2antistop1=0, A02amplitudesbottom1antisbottom2=0, A02amplitudesbottom2antisbottom1=0, A02amplitudeselectronLselectronR=0, A02amplitudeselectronRselectronL=0, A02amplitudesmuonLsmuonR=0, A02amplitudesmuonRsmuonL=0, A02amplitudestau1stau2=0, A02amplitudestau2stau1=0, A02amplitudegluongluon=0, A02amplitudegammagamma=0, A02amplitudeZgamma=0, A02amplitudeWHpm=0;
   
   double A02amplitudeneutZ1neutZ5=0, A02amplitudeneutZ2neutZ5=0, A02amplitudeneutZ3neutZ5=0, A02amplitudeneutZ4neutZ5=0, A02amplitudeneutZ5neutZ5=0, A02amplitudehiggsH3Zboson=0, A02amplitudehiggshA0=0, A02amplitudehiggsHA0=0, A02amplitudehiggsH3A0=0;

   if (flagA2 == 1) {
     if (QCDcorr == false) {
       ///No decays to u or d as PWs to u and d are tiny as proportional to yukawas squared
       ///Use running masses here to try to approximate some of the correction (which aren't included)
       A02amplitudecantic = higgsAamplitudedecayquarkantiquarkNMSSM(mA0(2), runmc, beta, CPOMix, 1, 2, QCDcorr, alphasAtMA2);
       A02amplitudesantis = higgsAamplitudedecayquarkantiquarkNMSSM(mA0(2), runms, beta, CPOMix, 0, 2, QCDcorr, alphasAtMA2);
       A02amplitudebantib = higgsAamplitudedecayquarkantiquarkNMSSM(mA0(2), runmb, beta, CPOMix, 0, 2, QCDcorr, alphasAtMA2);
       A02amplitudetantit = higgsAamplitudedecayquarkantiquarkNMSSM(mA0(2), runmt, beta, CPOMix, 1, 2, QCDcorr, alphasAtMA2);
     }
     else {
       A02amplitudecantic = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(2), mcpole, beta, CPOMix, 1, 2, QCDcorr, alphasAtMA2);
       A02amplitudesantis = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(2), mspole, beta, CPOMix, 0, 2, QCDcorr, alphasAtMA2);
       A02amplitudebantib = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(2), mbPole, beta, CPOMix, 0, 2, QCDcorr, alphasAtMA2);
       ///mcpole and mspole set in decays.h, this values used are those appropriate for the scheme used for the h -> qq QCD corrections, as in hdecay
       A02amplitudetantit = higgsAamplitudedecayquarkantiquarkNMSSM (mA0(2), mtPole, beta, CPOMix, 1, 2, QCDcorr, alphasAtMA2);
     }
     
     A02amplitudeeantie = higgsAamplitudedecayquarkantiquarkNMSSM(mA0(2), runmel, beta, CPOMix, 0, 2, false, alphasAtMA2)/3; ///0 as leptons are like down-type quarks, divide by 3 as No of colours is 1 for leptons cf 3 for quarks
     A02amplitudemuantimu = higgsAamplitudedecayquarkantiquarkNMSSM(mA0(2), runmmu, beta, CPOMix, 0, 2, false, alphasAtMA2)/3;
     A02amplitudetauantitau = higgsAamplitudedecayquarkantiquarkNMSSM(mA0(2), runmtau, beta, CPOMix, 0, 2, false, alphasAtMA2)/3;
     
     A02amplitudeneutZ1neutZ1 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(1), mneut(1), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 1, 1, 2);
     A02amplitudeneutZ1neutZ2 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(1), mneut(2), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 1, 2, 2);
     A02amplitudeneutZ1neutZ3 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(1), mneut(3), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 1, 3, 2);
     A02amplitudeneutZ1neutZ4 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(1), mneut(4), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 1, 4, 2);
     A02amplitudeneutZ1neutZ5 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(1), mneut(5), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 1, 5, 2);
     A02amplitudeneutZ2neutZ2 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(2), mneut(2), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 2, 2, 2);
     A02amplitudeneutZ2neutZ3 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(2), mneut(3), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 2, 3, 2);
     A02amplitudeneutZ2neutZ4 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(2), mneut(4), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 2, 4, 2);
     A02amplitudeneutZ2neutZ5 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(2), mneut(5), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 2, 5, 2);
     A02amplitudeneutZ3neutZ3 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(3), mneut(3), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 3, 3, 2);
     A02amplitudeneutZ3neutZ4 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(3), mneut(4), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 3, 4, 2);
     A02amplitudeneutZ3neutZ5 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(3), mneut(5), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 3, 5, 2);
     A02amplitudeneutZ4neutZ4 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(4), mneut(4), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 4, 4, 2);
     A02amplitudeneutZ4neutZ5 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(4), mneut(5), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 4, 5, 2);
     A02amplitudeneutZ5neutZ5 = higgsAamplitudedecayneutralinoneutralinoNMSSM (mA0(2), mneut(5), mneut(5), g, tanthetaW, lam, kappa, CPOMix, mixNeut, 5, 5, 2);
     A02amplitudecharW1charW1 = higgsAamplitudedecaysamecharginoNMSSM (mA0(2), mch(1), g, thetaL2, thetaR2, alpha, lam, CPOMix, 1, 2);
     A02amplitudecharW2charW2 = higgsAamplitudedecaysamecharginoNMSSM (mA0(2), mch(2), g, thetaL2, thetaR2, alpha, lam, CPOMix, 2, 2);
     A02amplitudecharW1charW2 = higgsAamplitudedecaydifcharginoNMSSM (mA0(2), mch(1), mch(2), g, thetaL2, thetaR2, alpha, lam, CPOMix, 2);
     A02amplitudehiggshZboson = higgsAamplitudedecayhiggshorHZbosonNMSSM (mA0(2), polemz, mh0(1), g, gp, alpha, beta, thetaA, CPEMix, 2, 1);
     A02amplitudehiggsHZboson = higgsAamplitudedecayhiggshorHZbosonNMSSM (mA0(2), polemz, mh0(2), g, gp, alpha, beta, thetaA, CPEMix, 2, 2);
     A02amplitudehiggsH3Zboson = higgsAamplitudedecayhiggshorHZbosonNMSSM (mA0(2), polemz, mh0(3), g, gp, alpha, beta, thetaA, CPEMix, 2, 3);     
     A02amplitudesupLantisupR = 3*higgsAamplitudedecaysfermionsNMSSM (mA0(2), mu(1,1), mu(2,1), g, runmw, runmu, Au, beta, lam, mueff, CPOMix, 'u', 2);
     A02amplitudesupRantisupL = 3*higgsAamplitudedecaysfermionsNMSSM (mA0(2), mu(2,1), mu(1,1), g, runmw, runmu, Au, beta, lam, mueff, CPOMix, 'u', 2);
     A02amplitudesdownLantisdownR = 3*higgsAamplitudedecaysfermionsNMSSM (mA0(2), md(1,1), md(2,1), g, runmw, runmd, Ad, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudesdownRantisdownL = 3*higgsAamplitudedecaysfermionsNMSSM (mA0(2), md(2,1), md(1,1), g, runmw, runmd, Ad, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudescharmLantischarmR = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(2), mu(1,2), mu(2,2), g, runmw, runmc, Ac, beta, lam, mueff, CPOMix, 'u', 2);
     A02amplitudescharmRantischarmL = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(2), mu(2,2), mu(1,2), g, runmw, runmc, Ac, beta, lam, mueff, CPOMix, 'u', 2);
     A02amplitudesstrangeLantisstrangeR = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(2), md(1,2), md(2,2), g, runmw, runms, As, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudesstrangeRantisstrangeL = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(2), md(2,2), md(1,2), g, runmw, runms, As, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudestop1antistop2 = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(2), mu(1,3), mu(2,3), g, runmw, runmt, At, beta, lam, mueff, CPOMix, 'u', 2);
     A02amplitudestop2antistop1 = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(2), mu(2,3), mu(1,3), g, runmw, runmt, At, beta, lam, mueff, CPOMix, 'u', 2);
     A02amplitudesbottom1antisbottom2 = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(2), md(1,3), md(2,3), g, runmw, runmb, Ab, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudesbottom2antisbottom1 = 3*higgsAamplitudedecaysfermionsNMSSM(mA0(2), md(2,3), md(1,3), g, runmw, runmb, Ab, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudeselectronLselectronR = higgsAamplitudedecaysfermionsNMSSM(mA0(2), me(1,1), me(2,1), g, runmw, runmel, Ae, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudeselectronRselectronL = higgsAamplitudedecaysfermionsNMSSM(mA0(2), me(2,1), me(1,1), g, runmw, runmel, Ae, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudesmuonLsmuonR = higgsAamplitudedecaysfermionsNMSSM(mA0(2), me(1,2), me(2,2), g, runmw, runmmu, Amu, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudesmuonRsmuonL = higgsAamplitudedecaysfermionsNMSSM(mA0(2), me(2,2), me(1,2), g, runmw, runmmu, Amu, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudestau1stau2 = higgsAamplitudedecaysfermionsNMSSM(mA0(2), me(1,3), me(2,3), g, runmw, runmtau, Atau, beta, lam, mueff, CPOMix, 'd', 2);
     A02amplitudestau2stau1 = higgsAamplitudedecaysfermionsNMSSM(mA0(2), me(2,3), me(1,3), g, runmw, runmtau, Atau, beta, lam, mueff, CPOMix, 'd', 2); 
     
     A02amplitudegluongluon = higgsAamplitudedecaygluongluonNMSSM (mA0(2), g, gs, alphasAtMA2, runmw, CPOMix, beta, mtPole, mbPole, mcpole, lam, 2, QCDcorr); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     A02amplitudegammagamma = higgsAamplitudedecaygammagammaNMSSM (mA0(2), g, gp, alphaAtMA2, runmw, CPOMix, beta, mtAtMA2, mbAtMA2, mcAtMA2, runmtau, mch(1), mch(2), thetaL2, thetaR2, lam, 2); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW
     A02amplitudeZgamma = higgsAamplitudedecayZgammaNMSSM (mA0(2), g, gp, alphaAtMA2, runmw, polemz, CPOMix, beta, mtAtMA2, mbAtMA2, mcAtMA2, mch(1), mch(2), thetaL2, thetaR2, lam, 2); ///Use quark masses and gauge couplings run to the mass of the decaying higgs, exact scale these were evaluated may significantly alter the PW

     A02amplitudeWHpm = higgsAamplitudedecayHpmWboson(mA0(2), polemw, mHpm, g, thetaA, 2, nmssmIsIt);
     A02amplitudehiggshA0 = higgsA2amplitudedecayA1CPevenNMSSM(mA0(2), mA0(1), mh0(1), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 1);
     A02amplitudehiggsHA0 = higgsA2amplitudedecayA1CPevenNMSSM(mA0(2), mA0(1), mh0(2), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 2);
     A02amplitudehiggsH3A0 = higgsA2amplitudedecayA1CPevenNMSSM(mA0(2), mA0(1), mh0(3), runmw, runmt, runmb, g, gp, beta, CPEMix, CPOMix, lam, kappa, Alambda, Akappa, mueff, 3);
     
     ParticleHiggsA2.Array_Decays[0][0] = PDGup; ParticleHiggsA2.Array_Decays[0][1] = -PDGup; ParticleHiggsA2.Array_Decays[0][2] = A02amplitudeuantiu; ParticleHiggsA2.Array_Decays[0][3] = 2; ParticleHiggsA2.Array_Comments[0] = "# A2 -> u ub";
     ParticleHiggsA2.Array_Decays[1][0] = PDGdown; ParticleHiggsA2.Array_Decays[1][1] = -PDGdown; ParticleHiggsA2.Array_Decays[1][2] = A02amplitudedantid; ParticleHiggsA2.Array_Decays[1][3] = 2; ParticleHiggsA2.Array_Comments[1] = "# A2 -> d db";
     ParticleHiggsA2.Array_Decays[2][0] = PDGcharm; ParticleHiggsA2.Array_Decays[2][1] = -PDGcharm; ParticleHiggsA2.Array_Decays[2][2] = A02amplitudecantic; ParticleHiggsA2.Array_Decays[2][3] = 2; ParticleHiggsA2.Array_Comments[2] = "# A2 -> c cb";
     ParticleHiggsA2.Array_Decays[3][0] = PDGstrange; ParticleHiggsA2.Array_Decays[3][1] = -PDGstrange; ParticleHiggsA2.Array_Decays[3][2] = A02amplitudesantis; ParticleHiggsA2.Array_Decays[3][3] = 2; ParticleHiggsA2.Array_Comments[3] = "# A2 -> s sb";
     ParticleHiggsA2.Array_Decays[4][0] = PDGbottom; ParticleHiggsA2.Array_Decays[4][1] = -PDGbottom; ParticleHiggsA2.Array_Decays[4][2] = A02amplitudebantib; ParticleHiggsA2.Array_Decays[4][3] = 2; ParticleHiggsA2.Array_Comments[4] = "# A2 -> b bb";
     ParticleHiggsA2.Array_Decays[5][0] = PDGtop; ParticleHiggsA2.Array_Decays[5][1] = -PDGtop; ParticleHiggsA2.Array_Decays[5][2] = A02amplitudetantit; ParticleHiggsA2.Array_Decays[5][3] = 2; ParticleHiggsA2.Array_Comments[5] = "# A2 -> t tb";
     ParticleHiggsA2.Array_Decays[6][0] = PDGelectron; ParticleHiggsA2.Array_Decays[6][1] = -PDGelectron; ParticleHiggsA2.Array_Decays[6][2] = A02amplitudeeantie; ParticleHiggsA2.Array_Decays[6][3] = 2; ParticleHiggsA2.Array_Comments[6] = "# A2 -> e- e+";
     ParticleHiggsA2.Array_Decays[7][0] = PDGmuon; ParticleHiggsA2.Array_Decays[7][1] = -PDGmuon; ParticleHiggsA2.Array_Decays[7][2] = A02amplitudemuantimu; ParticleHiggsA2.Array_Decays[7][3] = 2; ParticleHiggsA2.Array_Comments[7] = "# A2 -> mu- mu+";
     ParticleHiggsA2.Array_Decays[8][0] = PDGtau; ParticleHiggsA2.Array_Decays[8][1] = -PDGtau; ParticleHiggsA2.Array_Decays[8][2] = A02amplitudetauantitau; ParticleHiggsA2.Array_Decays[8][3] = 2; ParticleHiggsA2.Array_Comments[8] = "# A2 -> tau- tau+";
     
     ParticleHiggsA2.Array_Decays[9][0] = PDGneutralino1; ParticleHiggsA2.Array_Decays[9][1] = PDGneutralino1; ParticleHiggsA2.Array_Decays[9][2] = A02amplitudeneutZ1neutZ1; ParticleHiggsA2.Array_Decays[9][3] = 2; ParticleHiggsA2.Array_Comments[9] = "# A2 -> ~chi_10 ~chi_10";
     ParticleHiggsA2.Array_Decays[10][0] = PDGneutralino1; ParticleHiggsA2.Array_Decays[10][1] = PDGneutralino2; ParticleHiggsA2.Array_Decays[10][2] = A02amplitudeneutZ1neutZ2; ParticleHiggsA2.Array_Decays[10][3] = 2; ParticleHiggsA2.Array_Comments[10] = "# A2 -> ~chi_10 ~chi_20";
     ParticleHiggsA2.Array_Decays[11][0] = PDGneutralino1; ParticleHiggsA2.Array_Decays[11][1] = PDGneutralino3; ParticleHiggsA2.Array_Decays[11][2] = A02amplitudeneutZ1neutZ3; ParticleHiggsA2.Array_Decays[11][3] = 2; ParticleHiggsA2.Array_Comments[11] = "# A2 -> ~chi_10 ~chi_30";
     ParticleHiggsA2.Array_Decays[12][0] = PDGneutralino1; ParticleHiggsA2.Array_Decays[12][1] = PDGneutralino4; ParticleHiggsA2.Array_Decays[12][2] = A02amplitudeneutZ1neutZ4; ParticleHiggsA2.Array_Decays[12][3] = 2; ParticleHiggsA2.Array_Comments[12] = "# A2 -> ~chi_10 ~chi_40";
     ParticleHiggsA2.Array_Decays[13][0] = PDGneutralino2; ParticleHiggsA2.Array_Decays[13][1] = PDGneutralino2; ParticleHiggsA2.Array_Decays[13][2] = A02amplitudeneutZ2neutZ2; ParticleHiggsA2.Array_Decays[13][3] = 2; ParticleHiggsA2.Array_Comments[13] = "# A2 -> ~chi_20 ~chi_20";
     ParticleHiggsA2.Array_Decays[14][0] = PDGneutralino2; ParticleHiggsA2.Array_Decays[14][1] = PDGneutralino3; ParticleHiggsA2.Array_Decays[14][2] = A02amplitudeneutZ2neutZ3; ParticleHiggsA2.Array_Decays[14][3] = 2; ParticleHiggsA2.Array_Comments[14] = "# A2 -> ~chi_20 ~chi_30";
     ParticleHiggsA2.Array_Decays[15][0] = PDGneutralino2; ParticleHiggsA2.Array_Decays[15][1] = PDGneutralino4; ParticleHiggsA2.Array_Decays[15][2] = A02amplitudeneutZ2neutZ4; ParticleHiggsA2.Array_Decays[15][3] = 2; ParticleHiggsA2.Array_Comments[15] = "# A2 -> ~chi_20 ~chi_40";
     ParticleHiggsA2.Array_Decays[16][0] = PDGneutralino3; ParticleHiggsA2.Array_Decays[16][1] = PDGneutralino3; ParticleHiggsA2.Array_Decays[16][2] = A02amplitudeneutZ3neutZ3; ParticleHiggsA2.Array_Decays[16][3] = 2; ParticleHiggsA2.Array_Comments[16] = "# A2 -> ~chi_30 ~chi_30";
     ParticleHiggsA2.Array_Decays[17][0] = PDGneutralino3; ParticleHiggsA2.Array_Decays[17][1] = PDGneutralino4; ParticleHiggsA2.Array_Decays[17][2] = A02amplitudeneutZ3neutZ4; ParticleHiggsA2.Array_Decays[17][3] = 2; ParticleHiggsA2.Array_Comments[17] = "# A2 -> ~chi_30 ~chi_40";
     ParticleHiggsA2.Array_Decays[18][0] = PDGneutralino4; ParticleHiggsA2.Array_Decays[18][1] = PDGneutralino4; ParticleHiggsA2.Array_Decays[18][2] = A02amplitudeneutZ4neutZ4; ParticleHiggsA2.Array_Decays[18][3] = 2; ParticleHiggsA2.Array_Comments[18] = "# A2 -> ~chi_40 ~chi_40";
     
     ParticleHiggsA2.Array_Decays[19][0] = PDGchargino1; ParticleHiggsA2.Array_Decays[19][1] = -PDGchargino1; ParticleHiggsA2.Array_Decays[19][2] = A02amplitudecharW1charW1; ParticleHiggsA2.Array_Decays[19][3] = 2; ParticleHiggsA2.Array_Comments[19] = "# A2 -> ~chi_1+ ~chi_1-";
     ParticleHiggsA2.Array_Decays[20][0] = PDGchargino2; ParticleHiggsA2.Array_Decays[20][1] = -PDGchargino2; ParticleHiggsA2.Array_Decays[20][2] = A02amplitudecharW2charW2; ParticleHiggsA2.Array_Decays[20][3] = 2; ParticleHiggsA2.Array_Comments[20] = "# A2 -> ~chi_2+ ~chi_2-";
     ParticleHiggsA2.Array_Decays[21][0] = PDGchargino1; ParticleHiggsA2.Array_Decays[21][1] = -PDGchargino2; ParticleHiggsA2.Array_Decays[21][2] = A02amplitudecharW1charW2; ParticleHiggsA2.Array_Decays[21][3] = 2; ParticleHiggsA2.Array_Comments[21] = "# A2 -> ~chi_1+ ~chi_2-";
     ParticleHiggsA2.Array_Decays[22][0] = PDGchargino2; ParticleHiggsA2.Array_Decays[22][1] = -PDGchargino1; ParticleHiggsA2.Array_Decays[22][2] = A02amplitudecharW1charW2; ParticleHiggsA2.Array_Decays[22][3] = 2; ParticleHiggsA2.Array_Comments[22] = "# A2 -> ~chi_2+ ~chi_1-"; ///amplitude same as decay to W1+ and W2- by CP invariance

     
     ParticleHiggsA2.Array_Decays[23][0] = PDGZboson; ParticleHiggsA2.Array_Decays[23][1] = PDGh0; ParticleHiggsA2.Array_Decays[23][2] = A02amplitudehiggshZboson; ParticleHiggsA2.Array_Decays[23][3] = 2; ParticleHiggsA2.Array_Comments[23] = "# A2 -> h Z";
     ParticleHiggsA2.Array_Decays[24][0] = PDGZboson; ParticleHiggsA2.Array_Decays[24][1] = PDGH0; ParticleHiggsA2.Array_Decays[24][2] = A02amplitudehiggsHZboson; ParticleHiggsA2.Array_Decays[24][3] = 2; ParticleHiggsA2.Array_Comments[24] = "# A2 -> H Z";
     
     ParticleHiggsA2.Array_Decays[25][0] = PDGsupL; ParticleHiggsA2.Array_Decays[25][1] = PDGsupR; ParticleHiggsA2.Array_Decays[25][2] = A02amplitudesupLantisupR; ParticleHiggsA2.Array_Decays[25][3] = 2; ParticleHiggsA2.Array_Comments[25] = "# A2-> ~u_L ~u_R*"; 
     ParticleHiggsA.Array_Decays[26][0] = PDGsupR; ParticleHiggsA2.Array_Decays[26][1] = PDGsupL; ParticleHiggsA2.Array_Decays[26][2] = A02amplitudesupRantisupL; ParticleHiggsA2.Array_Decays[26][3] = 2; ParticleHiggsA2.Array_Comments[26] = "# A2-> ~u_R ~u_L*";
     ParticleHiggsA2.Array_Decays[27][0] = PDGsdownL; ParticleHiggsA2.Array_Decays[27][1] = PDGsdownR; ParticleHiggsA2.Array_Decays[27][2] = A02amplitudesdownLantisdownR; ParticleHiggsA2.Array_Decays[27][3] = 2; ParticleHiggsA2.Array_Comments[27] = "# A2-> ~d_L ~d_R*";
     ParticleHiggsA2.Array_Decays[28][0] = PDGsdownR; ParticleHiggsA2.Array_Decays[28][1] = PDGsdownL; ParticleHiggsA2.Array_Decays[28][2] = A02amplitudesdownRantisdownL; ParticleHiggsA2.Array_Decays[28][3] = 2; ParticleHiggsA2.Array_Comments[28] = "# A2-> ~d_R ~d_L*";
     ParticleHiggsA2.Array_Decays[29][0] = PDGscharmL; ParticleHiggsA2.Array_Decays[29][1] = PDGscharmR; ParticleHiggsA2.Array_Decays[29][2] = A02amplitudescharmLantischarmR; ParticleHiggsA2.Array_Decays[29][3] = 2; ParticleHiggsA2.Array_Comments[29] = "# A2-> ~c_L ~c_R*"; 
     ParticleHiggsA2.Array_Decays[30][0] = PDGscharmR; ParticleHiggsA2.Array_Decays[30][1] = PDGscharmL; ParticleHiggsA2.Array_Decays[30][2] = A02amplitudescharmRantischarmL; ParticleHiggsA2.Array_Decays[30][3] = 2; ParticleHiggsA2.Array_Comments[30] = "# A2-> ~c_R ~c_L*";
     ParticleHiggsA2.Array_Decays[31][0] = PDGsstrangeL; ParticleHiggsA2.Array_Decays[31][1] = PDGsstrangeR; ParticleHiggsA2.Array_Decays[31][2] = A02amplitudesstrangeLantisstrangeR; ParticleHiggsA2.Array_Decays[31][3] = 2; ParticleHiggsA2.Array_Comments[31] = "# A2-> ~s_L ~s_R*";
     ParticleHiggsA2.Array_Decays[32][0] = PDGsstrangeR; ParticleHiggsA2.Array_Decays[32][1] = PDGsstrangeL; ParticleHiggsA2.Array_Decays[32][2] = A02amplitudesstrangeRantisstrangeL; ParticleHiggsA2.Array_Decays[32][3] = 2; ParticleHiggsA2.Array_Comments[32] = "# A2-> ~s_R ~s_L*";
     ParticleHiggsA2.Array_Decays[33][0] = PDGstop1; ParticleHiggsA2.Array_Decays[33][1] = PDGstop2; ParticleHiggsA2.Array_Decays[33][2] = A02amplitudestop1antistop2; ParticleHiggsA2.Array_Decays[33][3] = 2; ParticleHiggsA2.Array_Comments[33] = "# A2-> ~t_1 ~t_2*"; 
     ParticleHiggsA2.Array_Decays[34][0] = PDGstop2; ParticleHiggsA2.Array_Decays[34][1] = PDGstop1; ParticleHiggsA2.Array_Decays[34][2] = A02amplitudestop2antistop1; ParticleHiggsA2.Array_Decays[34][3] = 2; ParticleHiggsA2.Array_Comments[34] = "# A2-> ~t_2 ~t_1*";
     ParticleHiggsA2.Array_Decays[35][0] = PDGsbottom1; ParticleHiggsA2.Array_Decays[35][1] = PDGsbottom2; ParticleHiggsA2.Array_Decays[35][2] = A02amplitudesbottom1antisbottom2; ParticleHiggsA2.Array_Decays[35][3] = 2; ParticleHiggsA2.Array_Comments[35] = "# A2-> ~b_1 ~b_2*";
     ParticleHiggsA2.Array_Decays[36][0] = PDGsbottom2; ParticleHiggsA2.Array_Decays[36][1] = PDGsbottom1; ParticleHiggsA2.Array_Decays[36][2] = A02amplitudesbottom2antisbottom1; ParticleHiggsA2.Array_Decays[36][3] = 2; ParticleHiggsA2.Array_Comments[36] = "# A2-> ~b_2 ~b_1*";
     ParticleHiggsA2.Array_Decays[37][0] = PDGselectronL; ParticleHiggsA2.Array_Decays[37][1] = PDGselectronR; ParticleHiggsA2.Array_Decays[37][2] = A02amplitudeselectronLselectronR; ParticleHiggsA2.Array_Decays[37][3] = 2; ParticleHiggsA2.Array_Comments[37] = "# A2-> ~e_L- ~e_R+"; 
     ParticleHiggsA.Array_Decays[38][0] = PDGselectronR; ParticleHiggsA2.Array_Decays[38][1] = PDGselectronL; ParticleHiggsA2.Array_Decays[38][2] = A02amplitudeselectronRselectronL; ParticleHiggsA2.Array_Decays[38][3] = 2; ParticleHiggsA2.Array_Comments[38] = "# A2-> ~e_R- ~e_L+";
     ParticleHiggsA2.Array_Decays[39][0] = PDGsmuonL; ParticleHiggsA2.Array_Decays[39][1] = PDGsmuonR; ParticleHiggsA2.Array_Decays[39][2] = A02amplitudesmuonLsmuonR; ParticleHiggsA2.Array_Decays[39][3] = 2; ParticleHiggsA2.Array_Comments[39] = "# A2-> ~mu_L- ~mu_R+";
     ParticleHiggsA2.Array_Decays[40][0] = PDGsmuonR; ParticleHiggsA2.Array_Decays[40][1] = PDGsmuonL; ParticleHiggsA2.Array_Decays[40][2] = A02amplitudesmuonRsmuonL; ParticleHiggsA2.Array_Decays[40][3] = 2; ParticleHiggsA2.Array_Comments[40] = "# A2-> ~mu_R- ~mu_L+";
     ParticleHiggsA2.Array_Decays[41][0] = PDGstau1; ParticleHiggsA2.Array_Decays[41][1] = PDGstau2; ParticleHiggsA2.Array_Decays[41][2] = A02amplitudestau1stau2; ParticleHiggsA2.Array_Decays[41][3] = 2; ParticleHiggsA2.Array_Comments[41] = "# A2-> ~tau_1- ~tau_2+";
     ParticleHiggsA2.Array_Decays[42][0] = PDGstau2; ParticleHiggsA2.Array_Decays[42][1] = PDGstau1; ParticleHiggsA2.Array_Decays[42][2] = A02amplitudestau2stau1; ParticleHiggsA2.Array_Decays[42][3] = 2; ParticleHiggsA2.Array_Comments[42] = "# A2-> ~tau_2- ~tau_1+";
     
     ParticleHiggsA2.Array_Decays[43][0] = PDGphoton; ParticleHiggsA2.Array_Decays[43][1] = PDGphoton; ParticleHiggsA2.Array_Decays[43][2] = A02amplitudegammagamma; ParticleHiggsA2.Array_Decays[43][3] = 2; ParticleHiggsA2.Array_Comments[43] = "# A2-> gamma gamma";
     ParticleHiggsA2.Array_Decays[44][0] = PDGgluon; ParticleHiggsA2.Array_Decays[44][1] = PDGgluon; ParticleHiggsA2.Array_Decays[44][2] = A02amplitudegluongluon; ParticleHiggsA2.Array_Decays[44][3] = 2; ParticleHiggsA2.Array_Comments[44] = "# A2-> gluon gluon";
     ParticleHiggsA2.Array_Decays[45][0] = PDGZboson; ParticleHiggsA2.Array_Decays[45][1] = PDGphoton; ParticleHiggsA2.Array_Decays[45][2] = A02amplitudeZgamma; ParticleHiggsA2.Array_Decays[45][3] = 2; ParticleHiggsA2.Array_Comments[45] = "# A2 -> Z gamma"; 
     ParticleHiggsA2.Array_Decays[46][0] = PDGWplus; ParticleHiggsA2.Array_Decays[46][1] = -PDGHplus; ParticleHiggsA2.Array_Decays[46][2] = A02amplitudeWHpm*2; ParticleHiggsA2.Array_Decays[46][3] = 2; ParticleHiggsA2.Array_Comments[46] = "# A2 -> W+- H-+"; ///*2 so includes A2 -> W- H+
     
     ParticleHiggsA2.Array_Decays[47][0] = PDGneutralino1; ParticleHiggsA2.Array_Decays[47][1] = PDGneutralino5; ParticleHiggsA2.Array_Decays[47][2] = A02amplitudeneutZ1neutZ5; ParticleHiggsA2.Array_Decays[47][3] = 2; ParticleHiggsA2.Array_Comments[47] = "# A2 -> ~chi_10 ~chi_50";
     ParticleHiggsA2.Array_Decays[48][0] = PDGneutralino1; ParticleHiggsA2.Array_Decays[48][1] = PDGneutralino5; ParticleHiggsA2.Array_Decays[48][2] = A02amplitudeneutZ2neutZ5; ParticleHiggsA2.Array_Decays[48][3] = 2; ParticleHiggsA2.Array_Comments[48] = "# A2 -> ~chi_20 ~chi_50";
     ParticleHiggsA2.Array_Decays[49][0] = PDGneutralino1; ParticleHiggsA2.Array_Decays[49][1] = PDGneutralino5; ParticleHiggsA2.Array_Decays[49][2] = A02amplitudeneutZ3neutZ5; ParticleHiggsA2.Array_Decays[49][3] = 2; ParticleHiggsA2.Array_Comments[49] = "# A2 -> ~chi_30 ~chi_50";
     ParticleHiggsA2.Array_Decays[50][0] = PDGneutralino1; ParticleHiggsA2.Array_Decays[50][1] = PDGneutralino5; ParticleHiggsA2.Array_Decays[50][2] = A02amplitudeneutZ4neutZ5; ParticleHiggsA2.Array_Decays[50][3] = 2; ParticleHiggsA2.Array_Comments[50] = "# A2 -> ~chi_40 ~chi_50";
     ParticleHiggsA2.Array_Decays[51][0] = PDGneutralino1; ParticleHiggsA2.Array_Decays[51][1] = PDGneutralino5; ParticleHiggsA2.Array_Decays[51][2] = A02amplitudeneutZ5neutZ5; ParticleHiggsA2.Array_Decays[51][3] = 2; ParticleHiggsA2.Array_Comments[51] = "# A2 -> ~chi_50 ~chi_50";
     ParticleHiggsA2.Array_Decays[52][0] = PDGZboson; ParticleHiggsA2.Array_Decays[52][1] = PDGH3; ParticleHiggsA2.Array_Decays[52][2] = A02amplitudehiggsH3Zboson; ParticleHiggsA2.Array_Decays[52][3] = 2; ParticleHiggsA2.Array_Comments[52] = "# A2 -> H3 Z";
     ParticleHiggsA2.Array_Decays[53][0] = PDGh0; ParticleHiggsA2.Array_Decays[53][1] = PDGA0; ParticleHiggsA2.Array_Decays[53][2] = A02amplitudehiggshA0; ParticleHiggsA2.Array_Decays[53][3] = 2; ParticleHiggsA2.Array_Comments[53] = "# A2 -> h A";
     ParticleHiggsA2.Array_Decays[54][0] = PDGH0; ParticleHiggsA2.Array_Decays[54][1] = PDGA0; ParticleHiggsA2.Array_Decays[54][2] = A02amplitudehiggsHA0; ParticleHiggsA2.Array_Decays[54][3] = 2; ParticleHiggsA2.Array_Comments[54] = "# A2 -> H A";
     ParticleHiggsA2.Array_Decays[55][0] = PDGH3; ParticleHiggsA2.Array_Decays[55][1] = PDGA0; ParticleHiggsA2.Array_Decays[55][2] = A02amplitudehiggsH3A0; ParticleHiggsA2.Array_Decays[55][3] = 2; ParticleHiggsA2.Array_Comments[55] = "# A2 -> H3 A";

     for(int i = 0; i<ParticleHiggsA2.No_of_Decays; i++) {
       if (ParticleHiggsA2.Array_Decays[i][2] < 0) {
	 fout << "#warning! Partial Width for " << ParticleHiggsA2.Array_Comments[i] << " is negative = " << ParticleHiggsA2.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
	 ParticleHiggsA2.Array_Decays[i][2] = 0;
	 errorflag = -1;
       }
     }       
     
     double HiggsA2_No_1to2_Decays = 0;
     ParticleHiggsA2.two_width = 0;
     ParticleHiggsA2.three_width = 0;
     ParticleHiggsA2.total_width = 0;
     
     HiggsA2_No_1to2_Decays = ParticleHiggsA2.No_1to2_Decays + ParticleHiggsA2.No_NMSSM_Decays; 
     
     for (int j = 0; j<HiggsA2_No_1to2_Decays; j++) {
       ParticleHiggsA2.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
     }
     
     for (int j=0; j<HiggsA2_No_1to2_Decays; j++) {
       ParticleHiggsA2.two_width = ParticleHiggsA2.two_width + ParticleHiggsA2.Array_Decays[j][2];
     }
     for (int j=HiggsA2_No_1to2_Decays; j<ParticleHiggsA2.No_of_Decays; j++) {
       ParticleHiggsA2.three_width = ParticleHiggsA2.three_width + ParticleHiggsA2.Array_Decays[j][2];
     }
     
     for(int j=0; j<ParticleHiggsA2.No_of_Decays; j++) {
       ParticleHiggsA2.Array_Decays[j][4] = 0;
     }
     
     ///Note no 3 body decays for HiggsA2
     if ( ParticleHiggsA2.three_width != ParticleHiggsA2.three_width) /// Tests for a nan as only nans aren't equal to themselves
       {
	 fout << "# Three body decays give nan for HiggsA2 - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
	 errorflag = -1;
	 ParticleHiggsA2.No_of_Decays = HiggsA2_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
	 ParticleHiggsA2.total_width = ParticleHiggsA2.two_width;
       }
     else {
       ParticleHiggsA2.total_width = ParticleHiggsA2.two_width + ParticleHiggsA2.three_width;
     }
     
     if ( ParticleHiggsA2.total_width != ParticleHiggsA2.total_width) /// Tests for a nan as only nans aren't equal to themselves
       {
	 errorflag = -1;
	 // for (int i = 0; i<ParticleHiggsA2.No_of_Decays; i++) {
	 //   fout << ParticleHiggsA2.Array_Decays[i][2] << endl;
	 // }	    
	 throw( "nan in A02 heaviest pseudoscalar higgs total width \n");
       }
          
   }
 }
 
 ///Charged higgs+ decays
 
 double Hplusamplitudeupdown=0, Hplusamplitudecharmstrange=0, Hplusamplitudetopbottom=0, Hplusamplitudeelectronelectronneutrino=0, Hplusamplitudemuonmuonneutrino=0, Hplusamplitudetautauneutrino=0, Hplusamplitudeupstrange=0, Hplusamplitudeupbottom=0, Hplusamplitudecharmdown=0, Hplusamplitudecharmbottom=0, Hplusamplitudetopdown=0, Hplusamplitudetopstrange=0, HplusamplitudeneutZ1charW1=0, HplusamplitudeneutZ1charW2=0, HplusamplitudeneutZ2charW1=0, HplusamplitudeneutZ2charW2=0, HplusamplitudeneutZ3charW1=0, HplusamplitudeneutZ3charW2=0, HplusamplitudeneutZ4charW1=0, HplusamplitudeneutZ4charW2=0, HplusamplitudehW=0, HplusamplitudesupLsdownL=0, HplusamplitudesupRsdownL=0, HplusamplitudesupLsdownR=0, HplusamplitudesupRsdownR=0, HplusamplitudescharmLsstrangeL=0, HplusamplitudescharmLsstrangeR=0, HplusamplitudescharmRsstrangeL=0, HplusamplitudescharmRsstrangeR=0, Hplusamplitudestop1sbottom1=0, Hplusamplitudestop1sbottom2=0, Hplusamplitudestop2sbottom1=0, Hplusamplitudestop2sbottom2=0, HplusamplitudeselectronLsnue=0, HplusamplitudeselectronRsnue=0, HplusamplitudesmuonLsnumu=0, HplusamplitudesmuonRsnumu=0, Hplusamplitudestau1snutau=0, Hplusamplitudestau2snutau=0; ///note Hplus -> H W is kinematically forbidden if we take tree-level mass formulae in the MSSM

 double HplusamplitudeWH = 0, HplusamplitudeWH3 = 0, HplusamplitudeWA = 0, HplusamplitudeWA2 = 0, HplusamplitudeneutZ5charW1 = 0, HplusamplitudeneutZ5charW2 = 0;///for NMSSM

 if (flagHpm == 1) {
   Hplusamplitudeupdown = higgsHplusamplitudedecayquarkantiquark (mHpm, mup, mdo, g, runmw, beta, VCKM, 1, 1);
   Hplusamplitudecharmstrange = higgsHplusamplitudedecayquarkantiquark (mHpm, mc, ms, g, runmw, beta, VCKM, 2, 2);
   Hplusamplitudetopbottom = higgsHplusamplitudedecayquarkantiquark (mHpm, mt, mb, g, runmw, beta, VCKM, 3, 3);
   Hplusamplitudeupstrange = higgsHplusamplitudedecayquarkantiquark (mHpm, mup, ms, g, runmw, beta, VCKM, 1, 2);
   Hplusamplitudeupbottom = higgsHplusamplitudedecayquarkantiquark (mHpm, mup, mb, g, runmw, beta, VCKM, 1, 3);
   Hplusamplitudecharmdown = higgsHplusamplitudedecayquarkantiquark (mHpm, mc, mdo, g, runmw, beta, VCKM, 2, 1);
   Hplusamplitudecharmbottom = higgsHplusamplitudedecayquarkantiquark (mHpm, mc, mb, g, runmw, beta, VCKM, 2, 3);
   Hplusamplitudetopdown = higgsHplusamplitudedecayquarkantiquark (mHpm, mt, mdo, g, runmw, beta, VCKM, 3, 1);
   Hplusamplitudetopstrange = higgsHplusamplitudedecayquarkantiquark (mHpm, mt, ms, g, runmw, beta, VCKM, 3, 2);
   Hplusamplitudeelectronelectronneutrino = higgsHplusamplitudedecayquarkantiquark(mHpm, 0, mel, g, runmw, beta, I3, 1, 1)/3;
   Hplusamplitudemuonmuonneutrino = higgsHplusamplitudedecayquarkantiquark(mHpm, 0, mmu, g, runmw, beta, I3, 2, 2)/3;
   Hplusamplitudetautauneutrino = higgsHplusamplitudedecayquarkantiquark(mHpm, 0, mtau, g, runmw, beta, I3, 3, 3)/3;
   
   HplusamplitudesupLsdownL = higgsHplusamplitudedecaysquarksquark (mHpm, mu(1,1), md(1,1), g, beta, runmw, runmu, runmd, greekmu, Au, Ad) (1);
   HplusamplitudesupRsdownR = higgsHplusamplitudedecaysquarksquark (mHpm, mu(2,1), md(2,1), g, beta, runmw, runmu, runmd, greekmu, Au, Ad) (2);
   HplusamplitudesupLsdownR = higgsHplusamplitudedecaysquarksquark (mHpm, mu(1,1), md(2,1), g, beta, runmw, runmu, runmd, greekmu, Au, Ad) (3);
   HplusamplitudesupRsdownL = higgsHplusamplitudedecaysquarksquark (mHpm, mu(2,1), md(2,1), g, beta, runmw, runmu, runmd, greekmu, Au, Ad) (4);
   HplusamplitudescharmLsstrangeL = higgsHplusamplitudedecaysquarksquark (mHpm, mu(1,2), md(1,2), g, beta, runmw, runmc, runms, greekmu, Ac, As) (1);
   HplusamplitudescharmRsstrangeR = higgsHplusamplitudedecaysquarksquark (mHpm, mu(2,2), md(2,2), g, beta, runmw, runmc, runms, greekmu, Ac, As) (2);
   HplusamplitudescharmLsstrangeR = higgsHplusamplitudedecaysquarksquark (mHpm, mu(1,2), md(2,2), g, beta, runmw, runmc, runms, greekmu, Ac, As) (3);
   HplusamplitudescharmRsstrangeL = higgsHplusamplitudedecaysquarksquark (mHpm, mu(2,2), md(1,2), g, beta, runmw, runmc, runms, greekmu, Ac, As) (4);
   HplusamplitudeselectronLsnue = higgsHplusamplitudedecaysquarksquark (mHpm, msnu(1), me(1,1), g, beta, runmw, 0, runmel, greekmu, 0, Ae) (1)/3;
   HplusamplitudeselectronRsnue = higgsHplusamplitudedecaysquarksquark (mHpm, msnu(1), me(2,1), g, beta, runmw, 0, runmel, greekmu, 0, Ae) (3)/3;
   HplusamplitudesmuonLsnumu = higgsHplusamplitudedecaysquarksquark (mHpm, msnu(2), me(1,2), g, beta, runmw, 0, runmmu, greekmu, 0, Amu) (1)/3;
   HplusamplitudesmuonRsnumu = higgsHplusamplitudedecaysquarksquark (mHpm, msnu(2), me(2,2), g, beta, runmw, 0, runmmu, greekmu, 0, Amu) (3)/3;
   Hplusamplitudestop1sbottom1 = higgsHplusamplitudedecaysquarksquarkmix (mHpm, mu(1,3), md(1,3), g, beta, runmw, runmt, runmb, greekmu,  At, Ab, thetat, thetab) (1); 
   Hplusamplitudestop2sbottom2 = higgsHplusamplitudedecaysquarksquarkmix (mHpm, mu(2,3), md(2,3), g, beta, runmw, runmt, runmb, greekmu,  At, Ab, thetat, thetab) (2);  
   Hplusamplitudestop1sbottom2 = higgsHplusamplitudedecaysquarksquarkmix (mHpm, mu(1,3), md(2,3), g, beta, runmw, runmt, runmb, greekmu,  At, Ab, thetat, thetab) (3);
   Hplusamplitudestop2sbottom1 = higgsHplusamplitudedecaysquarksquarkmix (mHpm, mu(2,3), md(1,3), g, beta, runmw, runmt, runmb, greekmu,  At, Ab, thetat, thetab) (4);
   Hplusamplitudestau1snutau = higgsHplusamplitudedecaysquarksquarkmix (mHpm, msnu(3), me(1,3), g, beta, runmw, 0, runmtau, greekmu, 0, Atau, 0, thetatau-PI/2) (1)/3;
   Hplusamplitudestau2snutau = higgsHplusamplitudedecaysquarksquarkmix (mHpm, msnu(3), me(2,3), g, beta, runmw, 0, runmtau, greekmu, 0, Atau, 0, thetatau-PI/2) (3)/3;
   
   if (nmssmIsIt == false) { 
     HplusamplitudehW = higgsHplusamplitudedecayWbosonhiggsh (mHpm, polemw, mh0(1), g, alpha, beta);
     HplusamplitudeneutZ1charW1 = higgsHplusamplitudedecayneutralinochargino (mHpm, mneut(1), MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 1, 1);
     HplusamplitudeneutZ1charW2 = higgsHplusamplitudedecayneutralinochargino (mHpm, mneut(1), MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 1, 2);
     HplusamplitudeneutZ2charW1 = higgsHplusamplitudedecayneutralinochargino (mHpm, mneut(2), MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 2, 1);
     HplusamplitudeneutZ2charW2 = higgsHplusamplitudedecayneutralinochargino (mHpm, mneut(2), MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 2, 2);
     HplusamplitudeneutZ3charW1 = higgsHplusamplitudedecayneutralinochargino (mHpm, mneut(3), MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 3, 1);
     HplusamplitudeneutZ3charW2 = higgsHplusamplitudedecayneutralinochargino (mHpm, mneut(3), MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 3, 2);
     HplusamplitudeneutZ4charW1 = higgsHplusamplitudedecayneutralinochargino (mHpm, mneut(4), MCH1, g, gp, beta, thetaL2, thetaR2, mixNeut, 4, 1);
     HplusamplitudeneutZ4charW2 = higgsHplusamplitudedecayneutralinochargino (mHpm, mneut(4), MCH2, g, gp, beta, thetaL2, thetaR2, mixNeut, 4, 2);
   }
   
   else if (nmssmIsIt == true) {
     HplusamplitudehW = higgsCPevenamplitudedecayWHpmNMSSM (mHpm, polemw, mh0(1), beta, g, CPEMix, 1);
     HplusamplitudeWH = higgsCPevenamplitudedecayWHpmNMSSM (mHpm, polemw, mh0(2), beta, g, CPEMix, 2);
     HplusamplitudeWH3 = higgsCPevenamplitudedecayWHpmNMSSM (mHpm, polemw, mh0(3), beta, g, CPEMix, 3);
     HplusamplitudeWA = higgsAamplitudedecayHpmWboson(mHpm, polemw, mA0(1), g, thetaA, 1, nmssmIsIt);
     HplusamplitudeWA2 = higgsAamplitudedecayHpmWboson(mHpm, polemw, mA0(2), g, thetaA, 2, nmssmIsIt);
     HplusamplitudeneutZ1charW1 = HpmamplitudecharginojneutralinoiNMSSM (mHpm, MCH1, mneut(1), g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 1, 1);
     HplusamplitudeneutZ1charW2 = HpmamplitudecharginojneutralinoiNMSSM (mHpm, MCH2, mneut(1), g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 1, 2);
     HplusamplitudeneutZ2charW1 = HpmamplitudecharginojneutralinoiNMSSM (mHpm, MCH1, mneut(2), g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 2, 1);
     HplusamplitudeneutZ2charW2 = HpmamplitudecharginojneutralinoiNMSSM (mHpm, MCH2, mneut(2), g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 2, 2);
     HplusamplitudeneutZ3charW1 = HpmamplitudecharginojneutralinoiNMSSM (mHpm, MCH1, mneut(3), g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 3, 1);
     HplusamplitudeneutZ3charW2 = HpmamplitudecharginojneutralinoiNMSSM (mHpm, MCH2, mneut(3), g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 3, 2);
     HplusamplitudeneutZ4charW1 = HpmamplitudecharginojneutralinoiNMSSM (mHpm, MCH1, mneut(4), g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 4, 1);
     HplusamplitudeneutZ4charW2 = HpmamplitudecharginojneutralinoiNMSSM (mHpm, MCH2, mneut(4), g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 4, 2);
     HplusamplitudeneutZ5charW1 = HpmamplitudecharginojneutralinoiNMSSM (mHpm, MCH1, mneut(5), g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 5, 1);
     HplusamplitudeneutZ5charW2 = HpmamplitudecharginojneutralinoiNMSSM (mHpm, MCH2, mneut(5), g, gp, beta, thetaL2, thetaR2, lam, mixNeut, 5, 2);
   }
   
   ParticleHiggsplus.Array_Decays[0][0] = PDGup; ParticleHiggsplus.Array_Decays[0][1] = -PDGdown; ParticleHiggsplus.Array_Decays[0][2] = Hplusamplitudeupdown; ParticleHiggsplus.Array_Decays[0][3] = 2; ParticleHiggsplus.Array_Comments[0] = "# H+ -> u db";
   ParticleHiggsplus.Array_Decays[1][0] = PDGcharm; ParticleHiggsplus.Array_Decays[1][1] = -PDGstrange; ParticleHiggsplus.Array_Decays[1][2] = Hplusamplitudecharmstrange; ParticleHiggsplus.Array_Decays[1][3] = 2; ParticleHiggsplus.Array_Comments[1] = "# H+ -> c sb";
   ParticleHiggsplus.Array_Decays[2][0] = PDGtop; ParticleHiggsplus.Array_Decays[2][1] = -PDGbottom; ParticleHiggsplus.Array_Decays[2][2] = Hplusamplitudetopbottom; ParticleHiggsplus.Array_Decays[2][3] = 2; ParticleHiggsplus.Array_Comments[2] = "# H+ -> t bb";
   ParticleHiggsplus.Array_Decays[3][0] = PDGup; ParticleHiggsplus.Array_Decays[3][1] = -PDGstrange; ParticleHiggsplus.Array_Decays[3][2] = Hplusamplitudeupstrange; ParticleHiggsplus.Array_Decays[3][3] = 2; ParticleHiggsplus.Array_Comments[3] = "# H+ -> u sb";
   ParticleHiggsplus.Array_Decays[4][0] = PDGup; ParticleHiggsplus.Array_Decays[4][1] = -PDGbottom; ParticleHiggsplus.Array_Decays[4][2] = Hplusamplitudeupbottom; ParticleHiggsplus.Array_Decays[4][3] = 2; ParticleHiggsplus.Array_Comments[4] = "# H+ -> u bb";
   ParticleHiggsplus.Array_Decays[5][0] = PDGcharm; ParticleHiggsplus.Array_Decays[5][1] = -PDGdown; ParticleHiggsplus.Array_Decays[5][2] = Hplusamplitudecharmdown; ParticleHiggsplus.Array_Decays[5][3] = 2; ParticleHiggsplus.Array_Comments[5] = "# H+ -> c db";
   ParticleHiggsplus.Array_Decays[6][0] = PDGcharm; ParticleHiggsplus.Array_Decays[6][1] = -PDGbottom; ParticleHiggsplus.Array_Decays[6][2] = Hplusamplitudecharmbottom; ParticleHiggsplus.Array_Decays[6][3] = 2; ParticleHiggsplus.Array_Comments[6] = "# H+ -> c bb";
   ParticleHiggsplus.Array_Decays[7][0] = PDGtop; ParticleHiggsplus.Array_Decays[7][1] = -PDGdown; ParticleHiggsplus.Array_Decays[7][2] = Hplusamplitudetopdown; ParticleHiggsplus.Array_Decays[7][3] = 2; ParticleHiggsplus.Array_Comments[7] = "# H+ -> t db";
   ParticleHiggsplus.Array_Decays[8][0] = PDGtop; ParticleHiggsplus.Array_Decays[8][1] = -PDGstrange; ParticleHiggsplus.Array_Decays[8][2] = Hplusamplitudetopstrange; ParticleHiggsplus.Array_Decays[8][3] = 2; ParticleHiggsplus.Array_Comments[8] = "# H+ -> t sb";
   ParticleHiggsplus.Array_Decays[9][0] = PDGnuelectron; ParticleHiggsplus.Array_Decays[9][1] = -PDGelectron; ParticleHiggsplus.Array_Decays[9][2] = Hplusamplitudeelectronelectronneutrino; ParticleHiggsplus.Array_Decays[9][3] = 2; ParticleHiggsplus.Array_Comments[9] = "# H+ -> e+ nu_e";
   ParticleHiggsplus.Array_Decays[10][0] = PDGnumuon; ParticleHiggsplus.Array_Decays[10][1] = -PDGmuon; ParticleHiggsplus.Array_Decays[10][2] = Hplusamplitudemuonmuonneutrino; ParticleHiggsplus.Array_Decays[10][3] = 2; ParticleHiggsplus.Array_Comments[10] = "# H+ -> mu+ nu_mu";
   ParticleHiggsplus.Array_Decays[11][0] = PDGnutau; ParticleHiggsplus.Array_Decays[11][1] = -PDGtau; ParticleHiggsplus.Array_Decays[11][2] = Hplusamplitudetautauneutrino; ParticleHiggsplus.Array_Decays[11][3] = 2; ParticleHiggsplus.Array_Comments[11] = "# H+ -> tau+ nu_tau";
   ParticleHiggsplus.Array_Decays[12][0] = PDGneutralino1; ParticleHiggsplus.Array_Decays[12][1] = PDGchargino1; ParticleHiggsplus.Array_Decays[12][2] = HplusamplitudeneutZ1charW1; ParticleHiggsplus.Array_Decays[12][3] = 2; ParticleHiggsplus.Array_Comments[12] = "# H+ -> ~chi_10 ~chi_1+";
   ParticleHiggsplus.Array_Decays[13][0] = PDGneutralino1; ParticleHiggsplus.Array_Decays[13][1] = PDGchargino2; ParticleHiggsplus.Array_Decays[13][2] = HplusamplitudeneutZ1charW2; ParticleHiggsplus.Array_Decays[13][3] = 2; ParticleHiggsplus.Array_Comments[13] = "# H+ -> ~chi_10 ~chi_2+";
   ParticleHiggsplus.Array_Decays[14][0] = PDGneutralino2; ParticleHiggsplus.Array_Decays[14][1] = PDGchargino1; ParticleHiggsplus.Array_Decays[14][2] = HplusamplitudeneutZ2charW1; ParticleHiggsplus.Array_Decays[14][3] = 2; ParticleHiggsplus.Array_Comments[14] = "# H+ -> ~chi_20 ~chi_1+";
   ParticleHiggsplus.Array_Decays[15][0] = PDGneutralino2; ParticleHiggsplus.Array_Decays[15][1] = PDGchargino2; ParticleHiggsplus.Array_Decays[15][2] = HplusamplitudeneutZ2charW2; ParticleHiggsplus.Array_Decays[15][3] = 2; ParticleHiggsplus.Array_Comments[15] = "# H+ -> ~chi_20 ~chi_2+";
   ParticleHiggsplus.Array_Decays[16][0] = PDGneutralino3; ParticleHiggsplus.Array_Decays[16][1] = PDGchargino1; ParticleHiggsplus.Array_Decays[16][2] = HplusamplitudeneutZ3charW1; ParticleHiggsplus.Array_Decays[16][3] = 2; ParticleHiggsplus.Array_Comments[16] = "# H+ -> ~chi_30 ~chi_1+";
   ParticleHiggsplus.Array_Decays[17][0] = PDGneutralino3; ParticleHiggsplus.Array_Decays[17][1] = PDGchargino2; ParticleHiggsplus.Array_Decays[17][2] = HplusamplitudeneutZ3charW2; ParticleHiggsplus.Array_Decays[17][3] = 2; ParticleHiggsplus.Array_Comments[17] = "# H+ -> ~chi_30 ~chi_2+";
   ParticleHiggsplus.Array_Decays[18][0] = PDGneutralino4; ParticleHiggsplus.Array_Decays[18][1] = PDGchargino1; ParticleHiggsplus.Array_Decays[18][2] = HplusamplitudeneutZ4charW1; ParticleHiggsplus.Array_Decays[18][3] = 2; ParticleHiggsplus.Array_Comments[18] = "# H+ -> ~chi_40 ~chi_1+";
   ParticleHiggsplus.Array_Decays[19][0] = PDGneutralino4; ParticleHiggsplus.Array_Decays[19][1] = PDGchargino2; ParticleHiggsplus.Array_Decays[19][2] = HplusamplitudeneutZ4charW2; ParticleHiggsplus.Array_Decays[19][3] = 2; ParticleHiggsplus.Array_Comments[19] = "# H+ -> ~chi_40 ~chi_2+";
   ParticleHiggsplus.Array_Decays[20][0] = PDGWplus; ParticleHiggsplus.Array_Decays[20][1] = PDGh0; ParticleHiggsplus.Array_Decays[20][2] = HplusamplitudehW; ParticleHiggsplus.Array_Decays[20][3] = 2; ParticleHiggsplus.Array_Comments[20] = "# H+ -> W+ h";  
   ParticleHiggsplus.Array_Decays[21][0] = PDGsupL; ParticleHiggsplus.Array_Decays[21][1] = -PDGsdownL; ParticleHiggsplus.Array_Decays[21][2] = HplusamplitudesupLsdownL; ParticleHiggsplus.Array_Decays[21][3] = 2; ParticleHiggsplus.Array_Comments[21] = "# H+ -> ~u_L d_L*"; 
   ParticleHiggsplus.Array_Decays[22][0] = PDGsupR; ParticleHiggsplus.Array_Decays[22][1] = -PDGsdownR; ParticleHiggsplus.Array_Decays[22][2] = HplusamplitudesupRsdownR; ParticleHiggsplus.Array_Decays[22][3] = 2; ParticleHiggsplus.Array_Comments[22] = "# H+ -> ~u_R d_R*";
   ParticleHiggsplus.Array_Decays[23][0] = PDGsupL; ParticleHiggsplus.Array_Decays[23][1] = -PDGsdownR; ParticleHiggsplus.Array_Decays[23][2] = HplusamplitudesupLsdownR; ParticleHiggsplus.Array_Decays[23][3] = 2; ParticleHiggsplus.Array_Comments[23] = "# H+ -> ~u_L d_R*";
   ParticleHiggsplus.Array_Decays[24][0] = PDGsupR; ParticleHiggsplus.Array_Decays[24][1] = -PDGsdownL; ParticleHiggsplus.Array_Decays[24][2] = HplusamplitudesupRsdownL; ParticleHiggsplus.Array_Decays[24][3] = 2; ParticleHiggsplus.Array_Comments[24] = "# H+ -> ~u_R d_L*";
   ParticleHiggsplus.Array_Decays[25][0] = PDGscharmL; ParticleHiggsplus.Array_Decays[25][1] = -PDGsstrangeL; ParticleHiggsplus.Array_Decays[25][2] = HplusamplitudescharmLsstrangeL; ParticleHiggsplus.Array_Decays[25][3] = 2; ParticleHiggsplus.Array_Comments[25] = "# H+ -> ~c_L s_L*"; 
   ParticleHiggsplus.Array_Decays[26][0] = PDGscharmR; ParticleHiggsplus.Array_Decays[26][1] = -PDGsstrangeR; ParticleHiggsplus.Array_Decays[26][2] = HplusamplitudescharmRsstrangeR; ParticleHiggsplus.Array_Decays[26][3] = 2; ParticleHiggsplus.Array_Comments[26] = "# H+ -> ~c_R s_R*";
   ParticleHiggsplus.Array_Decays[27][0] = PDGscharmL; ParticleHiggsplus.Array_Decays[27][1] = -PDGsstrangeR; ParticleHiggsplus.Array_Decays[27][2] = HplusamplitudescharmLsstrangeR; ParticleHiggsplus.Array_Decays[27][3] = 2; ParticleHiggsplus.Array_Comments[27] = "# H+ -> ~c_L s_R*";
   ParticleHiggsplus.Array_Decays[28][0] = PDGscharmR; ParticleHiggsplus.Array_Decays[28][1] = -PDGsstrangeL; ParticleHiggsplus.Array_Decays[28][2] = HplusamplitudescharmRsstrangeL; ParticleHiggsplus.Array_Decays[28][3] = 2; ParticleHiggsplus.Array_Comments[28] = "# H+ -> ~c_R s_L*";
   ParticleHiggsplus.Array_Decays[29][0] = PDGnuselectronL; ParticleHiggsplus.Array_Decays[29][1] = -PDGselectronL; ParticleHiggsplus.Array_Decays[29][2] = HplusamplitudeselectronLsnue; ParticleHiggsplus.Array_Decays[29][3] = 2; ParticleHiggsplus.Array_Comments[29] = "# H+ -> ~e_L+ nu_eL";
   ParticleHiggsplus.Array_Decays[30][0] = PDGnuselectronL; ParticleHiggsplus.Array_Decays[30][1] = -PDGselectronR; ParticleHiggsplus.Array_Decays[30][2] = HplusamplitudeselectronRsnue; ParticleHiggsplus.Array_Decays[30][3] = 2; ParticleHiggsplus.Array_Comments[30] = "# H+ -> ~e_R+ nu_eL";
   ParticleHiggsplus.Array_Decays[31][0] = PDGnusmuonL; ParticleHiggsplus.Array_Decays[31][1] = -PDGsmuonL; ParticleHiggsplus.Array_Decays[31][2] = HplusamplitudesmuonLsnumu; ParticleHiggsplus.Array_Decays[31][3] = 2; ParticleHiggsplus.Array_Comments[31] = "# H+ -> ~mu_L+ nu_muL";
   ParticleHiggsplus.Array_Decays[32][0] = PDGnusmuonL; ParticleHiggsplus.Array_Decays[32][1] = -PDGsmuonR; ParticleHiggsplus.Array_Decays[32][2] = HplusamplitudesmuonRsnumu; ParticleHiggsplus.Array_Decays[32][3] = 2; ParticleHiggsplus.Array_Comments[32] = "# H+ -> ~mu_R+ nu_muL";
   ParticleHiggsplus.Array_Decays[33][0] = PDGstop1; ParticleHiggsplus.Array_Decays[33][1] = -PDGsbottom1; ParticleHiggsplus.Array_Decays[33][2] = Hplusamplitudestop1sbottom1; ParticleHiggsplus.Array_Decays[33][3] = 2; ParticleHiggsplus.Array_Comments[33] = "# H+ -> ~t_1 b_1*";
   ParticleHiggsplus.Array_Decays[34][0] = PDGstop2; ParticleHiggsplus.Array_Decays[34][1] = -PDGsbottom2; ParticleHiggsplus.Array_Decays[34][2] = Hplusamplitudestop2sbottom2; ParticleHiggsplus.Array_Decays[34][3] = 2; ParticleHiggsplus.Array_Comments[34] = "# H+ -> ~t_2 b_2*";
   ParticleHiggsplus.Array_Decays[35][0] = PDGstop1; ParticleHiggsplus.Array_Decays[35][1] = -PDGsbottom2; ParticleHiggsplus.Array_Decays[35][2] = Hplusamplitudestop1sbottom2; ParticleHiggsplus.Array_Decays[35][3] = 2; ParticleHiggsplus.Array_Comments[35] = "# H+ -> ~t_1 b_2*";
   ParticleHiggsplus.Array_Decays[36][0] = PDGstop2; ParticleHiggsplus.Array_Decays[36][1] = -PDGsbottom1; ParticleHiggsplus.Array_Decays[36][2] = Hplusamplitudestop2sbottom1; ParticleHiggsplus.Array_Decays[36][3] = 2; ParticleHiggsplus.Array_Comments[36] = "# H+ -> ~t_2 b_1*";
   ParticleHiggsplus.Array_Decays[37][0] = PDGnustauL; ParticleHiggsplus.Array_Decays[37][1] = -PDGstau1; ParticleHiggsplus.Array_Decays[37][2] = Hplusamplitudestau1snutau; ParticleHiggsplus.Array_Decays[37][3] = 2; ParticleHiggsplus.Array_Comments[37] = "# H+ -> ~tau_1+ nu_tauL";
   ParticleHiggsplus.Array_Decays[38][0] = PDGnustauL; ParticleHiggsplus.Array_Decays[38][1] = -PDGstau2; ParticleHiggsplus.Array_Decays[38][2] = Hplusamplitudestau2snutau; ParticleHiggsplus.Array_Decays[38][3] = 2; ParticleHiggsplus.Array_Comments[38] = "# H+ -> ~tau_2+ nu_tauL";
   
   ///NMSSM
   ParticleHiggsplus.Array_Decays[39][0] = PDGWplus; ParticleHiggsplus.Array_Decays[39][1] = PDGH0; ParticleHiggsplus.Array_Decays[39][2] = HplusamplitudeWH; ParticleHiggsplus.Array_Decays[39][3] = 2; ParticleHiggsplus.Array_Comments[39] = "# H+ -> W+ H"; 
   ParticleHiggsplus.Array_Decays[40][0] = PDGWplus; ParticleHiggsplus.Array_Decays[40][1] = PDGH3; ParticleHiggsplus.Array_Decays[40][2] = HplusamplitudeWH3; ParticleHiggsplus.Array_Decays[40][3] = 2; ParticleHiggsplus.Array_Comments[40] = "# H+ -> W+ H3";  
   ParticleHiggsplus.Array_Decays[41][0] = PDGWplus; ParticleHiggsplus.Array_Decays[41][1] = PDGA0; ParticleHiggsplus.Array_Decays[41][2] = HplusamplitudeWA; ParticleHiggsplus.Array_Decays[41][3] = 2; ParticleHiggsplus.Array_Comments[41] = "# H+ -> W+ A"; 
   ParticleHiggsplus.Array_Decays[42][0] = PDGWplus; ParticleHiggsplus.Array_Decays[42][1] = PDGA2; ParticleHiggsplus.Array_Decays[42][2] = HplusamplitudeWA2; ParticleHiggsplus.Array_Decays[42][3] = 2; ParticleHiggsplus.Array_Comments[42] = "# H+ -> W+ A2";
   ParticleHiggsplus.Array_Decays[43][0] = PDGneutralino5; ParticleHiggsplus.Array_Decays[43][1] = PDGchargino1; ParticleHiggsplus.Array_Decays[43][2] = HplusamplitudeneutZ5charW1; ParticleHiggsplus.Array_Decays[43][3] = 2; ParticleHiggsplus.Array_Comments[43] = "# H+ -> ~chi_50 ~chi_1+";
   ParticleHiggsplus.Array_Decays[44][0] = PDGneutralino5; ParticleHiggsplus.Array_Decays[44][1] = PDGchargino2; ParticleHiggsplus.Array_Decays[44][2] = HplusamplitudeneutZ5charW2; ParticleHiggsplus.Array_Decays[44][3] = 2; ParticleHiggsplus.Array_Comments[44] = "# H+ -> ~chi_50 ~chi_2+";

   for(int i = 0; i<ParticleHiggsplus.No_of_Decays; i++) {
     if (ParticleHiggsplus.Array_Decays[i][2] < 0) {
       fout << "#warning! Partial Width for " << ParticleHiggsplus.Array_Comments[i] << " is negative = " << ParticleHiggsplus.Array_Decays[i][2] << " and is set to zero here" << endl; //most likely due to numerical precision issues and so for small partial widths which are negligible if they are not the only decay modes available, contact the authors for more information or if you are concerned about your specific point producing a negative partial width in error
       ParticleHiggsplus.Array_Decays[i][2] = 0;
       errorflag = -1;
     }
   }       
   
   double Higgsplus_No_1to2_Decays = 0;
   
   Higgsplus_No_1to2_Decays = ParticleHiggsplus.No_1to2_Decays + ParticleHiggsplus.No_NMSSM_Decays; /// As higgsplus can't be NLSP as heavier than higgsl
   
   for (int j = 0; j<Higgsplus_No_1to2_Decays; j++) {
     ParticleHiggsplus.Array_Decays[j][4] = 0; ///0 indicates no 3rd daughter so 1->2 decay.
   }
   
   for (int j=0; j<Higgsplus_No_1to2_Decays; j++) {
     ParticleHiggsplus.two_width = ParticleHiggsplus.two_width + ParticleHiggsplus.Array_Decays[j][2];
   }
   
   for (int j=Higgsplus_No_1to2_Decays; j<ParticleHiggsplus.No_of_Decays; j++) {
     ParticleHiggsplus.three_width = ParticleHiggsplus.three_width + ParticleHiggsplus.Array_Decays[j][2];
   }
   
   for(int j=0; j<ParticleHiggsplus.No_of_Decays; j++) {
     ParticleHiggsplus.Array_Decays[j][4] = 0;
   }

   ///Note no 3 body decays for Higgsplus
   if ( ParticleHiggsplus.three_width != ParticleHiggsplus.three_width) /// Tests for a nan as only nans aren't equal to themselves
       {
	 fout << "# Three body decays give nan for Higgsplus - problem! Therefore total and partial widths and branching ratios output only includes 1->2 decays" << endl;
	 errorflag = -1;
	 ParticleHiggsplus.No_of_Decays = Higgsplus_No_1to2_Decays; ///So only 1 to 2 decays are output if a 1 to 3 decay gives a nan
	 ParticleHiggsplus.total_width = ParticleHiggsplus.two_width;
       }
   else {
     ParticleHiggsplus.total_width = ParticleHiggsplus.two_width + ParticleHiggsplus.three_width;
   }
   
   if ( ParticleHiggsplus.total_width != ParticleHiggsplus.total_width) /// Tests for a nan as only nans aren't equal to themselves
     {
       errorflag = -1;
       // for (int i = 0; i<ParticleHiggsplus.No_of_Decays; i++) {
       //   fout << ParticleHiggsplus.Array_Decays[i][2] << endl;
       // }	    
       throw( "nan in H+ total width \n");
     }
 }

 /// Construct decay table
 decayTable.push_back(ParticleGluino);
 decayTable.push_back(ParticleSdownL);
 decayTable.push_back(ParticleSdownR);
 decayTable.push_back(ParticleSupL);
 decayTable.push_back(ParticleSupR);
 decayTable.push_back(ParticleSstrangeL);
 decayTable.push_back(ParticleSstrangeR);
 decayTable.push_back(ParticleScharmL);
 decayTable.push_back(ParticleScharmR);
 decayTable.push_back(ParticleSbottom1);
 decayTable.push_back(ParticleSbottom2);
 decayTable.push_back(ParticleStop1);
 decayTable.push_back(ParticleStop2);
 decayTable.push_back(ParticleSelectronL);
 decayTable.push_back(ParticleSelectronR);
 decayTable.push_back(ParticleSmuonL);
 decayTable.push_back(ParticleSmuonR);
 decayTable.push_back(ParticleSnue);
 decayTable.push_back(ParticleSnumu);
 decayTable.push_back(ParticleStau1);
 decayTable.push_back(ParticleStau2);
 decayTable.push_back(ParticleSnutau);
 decayTable.push_back(ParticleChargino1);
 decayTable.push_back(ParticleChargino2);
 decayTable.push_back(ParticleNeutralino1);
 decayTable.push_back(ParticleNeutralino2);
 decayTable.push_back(ParticleNeutralino3);
 decayTable.push_back(ParticleNeutralino4);
 decayTable.push_back(Particlehiggsl);
 decayTable.push_back(ParticleHiggsH);
 decayTable.push_back(ParticleHiggsA);
 decayTable.push_back(ParticleHiggsplus);
 if (nmssmIsIt) {
   decayTable.push_back(ParticleHiggsA2);
   decayTable.push_back(ParticleHiggsH3);
   decayTable.push_back(ParticleNeutralino5);
 }

 /// calculate branching ratios
  vector<Particle>::iterator ii;
  for (ii=decayTable.begin(); ii<decayTable.end(); ii++)
    for (int i=0; i<ii->No_of_Decays; i++)
      ii->Array_Decays[i][5] = ii->Array_Decays[i][2] / ii->total_width;

  /// re-set problem flags to their values before decays were calculated
  r->setProblem(rProb);
  
 return errorflag;
}
