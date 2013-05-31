/** \file nmssmsoftsusy.cpp
    Project: NMSSMSOFTSUSY 
    Author: Ben Allanach, Peter Athron, Lewis Tunstall, Alexander Voigt 
    Manual: TBW
    Webpage:  https://github.com/Expander/softsusy.git 
*/

#include "nmssmsoftsusy.h"

#ifdef NMSSMSOFTSUSY_H

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
void NmssmSoftsusy::treeUpSquark(DoubleMatrix & mass, double mtrun, 
				double pizztMS, double sinthDRbarMS, 
				int family) { 
   //PA: only modification is to add lambda * s / root to mu
   double lam = displayLambda(), svev = displaySvev(), tanb = displayTanb();
   Softsusy::treeUpSquark(mass, mtrun, pizztMS, sinthDRbarMS, family);
   if (family == 3){ 
   mass(1, 2) = mass(1, 2) -  mtrun * lam * svev / (root2 * tanb);
   mass(2, 1) = mass(1, 2);
   }
   cout << "Dev: mass = "  << mass << endl;
   cout << "mtrun = "  << mtrun << endl;
  cout << "lam = " << lam << endl;
  cout << "mueff = " <<displaySusyMu() + lam * svev / root2  << endl;
  cout<< "displaySusyMu() = " << displaySusyMu() << endl;
  cout << "EPSTOL = "  << EPSTOL << endl;
  
  cout << "displaySoftA(UA, 3, 3) = "  << displaySoftA(UA, 3, 3) << endl;
  cout << "tanb = "  << tanb << endl;
  cout << "mt2 = "  << sqr(mtrun) << endl;
}



void NmssmSoftsusy::treeDownSquark(DoubleMatrix & mass, double mbrun, 
				double pizztMS, double sinthDRbarMS, 
				int family) { 
   //PA: only modification is to add lambda * s / root to mu
   double lam = displayLambda(), svev = displaySvev(), tanb = displayTanb();
   Softsusy::treeDownSquark(mass, mbrun, pizztMS, sinthDRbarMS, family);
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
   Softsusy::treeChargedSlepton(mass, mtaurun, pizztMS, sinthDRbarMS, family);
   if (family == 3) { 
   mass(1, 2) = mass(1, 2) -  mtaurun * lam * svev * tanb / (root2);
   mass(2, 1) = mass(1, 2);
   }
}

/// LCT: new routine to allocate NMSSM chargino masses
void NmssmSoftsusy::calcDrBarCharginos(DoubleMatrix & mass, double beta, double mw) {	
   double lam = displayLambda(), svev = displaySvev();
   
   Softsusy::calcDrBarCharginos(mass, beta, mw);
   mass(2, 2) = mass(2, 2) + lam * svev / root2;
}

/// LCT: new routine for NMSSM neutralino masses
void NmssmSoftsusy::calcDrBarNeutralinos(DoubleMatrix & mass, double beta, double mz, double mw, double sinthDRbar) {
	double lam = displayLambda(), kap = displayKappa();
	double smu = displaySusyMu(), mupr = displayMupr();
	double cosb = cos(beta), sinb = sin(beta);
	double vev = displayHvev(), svev = displaySvev();
	
	/// Call MSSM 4 x4 neutralino mass matrix
        Softsusy::calcDrBarNeutralinos(mass, beta, mz, mw, sinthDRbar);
        
	/// Fill remaining values
	mass(3, 4) = mass(3, 4) - lam * svev / root2;
	mass(3, 5) = - lam * vev * sinb / root2;
	mass(4, 5) = - lam * vev * cosb / root2;
	mass(5, 5) = root2 * kap * svev + mupr;
	
	/// symmetrise tree-level
	mass.symmetrise();

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
  //  cout << "mSq: " << mSq << endl;
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
  double sw2, guL, gdL, geL, guR, gdR, geR;
  double vev = displayHvev();
  
  Softsusy::setNeutCurrCouplings(sinthDRbar, sw2, guL, gdL, geL, guR, gdR, geR);
  Softsusy::calcDRTrilinears(eg, vev, beta);
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
  DoubleMatrix mCh(2, 2);   
  calcDrBarCharginos(mCh, beta, mw);
  eg.mch = mCh.asy2by2(eg.thetaL, eg.thetaR);
  eg.mpzCharginos();
  
  DoubleMatrix mNeut(5, 5);
  calcDrBarNeutralinos(mNeut, beta, mz, mw, sinthDRbar);
  if (mNeut.diagonaliseSym(eg.mixNeut, eg.mneut) > TOLERANCE *
       1.0e-3) {  
      ostringstream ii;
      ii << "accuracy bad in neutralino diagonalisation"<< flush;
      throw ii.str(); 
   }
 
  eg.mpzNeutralinos();
  eg.mw = mw;
  eg.mz = mz;
 
  calcDrBarHiggs(beta, mz2, mw2, sinthDRbar, eg);  
  
  setDrBarPars(eg);
  
  return;

}

#endif
